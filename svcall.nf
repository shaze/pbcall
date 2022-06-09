
// (c) University of the Witwatersrand, Johannesburg, 2022
// Author Scott Hazelhurst
// Released under MIT Licence




include { cuteSVAsIndividuals as cute  } from './callers/cuteSV.nf'
include { pbsvIndividuals  } from './callers/pbsv.nf'

nextflow.enable.dsl=2


params.output = "sample"
params.has_bam = false


trf = file(params.tandem_example)

ref_seq=file(params.ref_seq)
ref_fai=file(params.ref_fai)
ref_mmi=file(params.ref_mmi)

params.chrom_prefix="chr"

fnames = file(params.input,type:'dir')


def getBAMs(fnames) {
    got_bams=true
    fnames.each { 
	 base=it.baseName; 
	 if (! file("${params.bam}/${base}.bam").exists()) {
	   println "No file ${params.bam}/${base}.bam";
	   got_bams=false
	 } 
	 if (! file("${params.bam}/${base}.bam.bai").exists())  {
	   println "No file ${params.bam}/${base}.bam.bai";
	   got_bams=false
	 }
    }
    if (! got_bams) {
       println "Promised bams but at least some of them or their index files missing"
       System.exit(18)
    }
    globs = fname .collect { "${it}*" }.join(",")
    bf = Channel.fromFilePairs("${params.bam}/{${globs}}",size:3) \
	    { file -> file.simpleName }.\
	    map { b, f -> [b,f[0],f[1]] }
    return bf;
}


   // Create a BAM file from the reads aligning to the reference genome
   // This is done per sample
process pbio_bamify {
     cpus params.bamify_cpus
     memory params.bamify_mem
     maxForks 8
     errorStrategy 'finish'
     input:
       path(fq)
     output:
      tuple val(the_id), file("${the_id}.bam"), file("${the_id}.bam.bai"), emit: 'bam'
     storeDir params.bam
     script:
       the_id = fq.baseName
       """
       ls $fq/*.f*q.gz   > files.fofn
       pbmm2 align $ref_mmi  files.fofn ${the_id}.bam --preset HIFI --rg '@RG\tID:$the_id\tSM:$the_id' --unmapped --sort 
       """
}


autosomes = 1..22
chroms    = autosomes +  ['X','Y','M'] 




workflow {
    if (params.has_bam) {
       bf = getBAMs(fnames)
    }  else {
	pbio_bamify(Channel.fromPath(fnames))
        bf=pbio_bamify.out.bam
    }
    pbsv(bams,trf)
    cute(bams)
    svim(bams)
}
