
// (c) University of the Witwatersrand, Johannesburg
// Author Scott Hazelhurst
// Released under https://creativecommons.org/licenses/by-sa/4.0/


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
    globs = fnames .collect { "${it}*" }.join(",")
    println(globs)
    bf = Channel.fromFilePairs("${params.bam}/*bam*",size:2) \
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

// Do an initial calling of the variants
// Note that this is done per chromosome
// -- so for each sample of input we have multiple processes running
process deepcall {
  maxForks 40  
  label 'deepvariant'
  cpus params.call_cpus
  memory params.call_mem
  errorStrategy 'finish'
  if (params.constraint)
     clusterOptions="--constraint=${params.constraint}"
  input:
     tuple val(the_id), file(bam), file(bai) 
     file ref_seq
     file ref_fai
     each chroms
  output:
    tuple val(the_id), val(chroms), file("${vcf}.vcf.gz"), file(tbi), \
	file("${vcf}.gvcf.gz"), file("${vcf}.gvcf.gz.tbi"),\
	file(bam), file(bai)
    publishDir "$params.out/gvcf/$c",  pattern: "*.gvcf.gz"
    publishDir "$params.out/vcf/$c",  pattern: "*.vcf.gz"     
  script:
     name = bam.simpleName
     if (['X','Y','M'].contains(chroms))
        c=chroms
     else
         c = "${chroms}".padLeft(2,"0")
         vcf = "${name}-${c}"
         tbi = "${vcf}.vcf.gz.tbi"
     """
      hostname
      /opt/deepvariant/bin/run_deepvariant  \
                --model_type PACBIO \
                --ref $ref_seq  \
                --reads  $bam  \
                --output_vcf ${vcf}.vcf.gz \
                --output_gvcf=${vcf}.gvcf.gz \
                --num_shards ${params.call_cpus} \
                --regions chr$chroms
      """
}


    

workflow {
    if (params.has_bam) {
       bf = getBAMs(fnames)
    }  else {
	pbio_bamify(Channel.fromPath(fnames))
        bf=pbio_bamify.out.bam
    }
    deepcall(bf,ref_seq,ref_fai,chroms)
}
