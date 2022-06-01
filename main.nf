
// (c) University of the Witwatersrand, Johannesburg
// Author Scott Hazelhurst
// Released under https://creativecommons.org/licenses/by-sa/4.0/


nextflow.enable.dsl=2

params.output = "sample"
params.has_bam = false
trf = file(params.tandem_example)

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
     maxForks 5
     errorStrategy 'finish'
     input:
       path (fq)
     output:
       tuple val(the_id), file("${the_id}.bam"), file("${the_id}.bam.bai"), emit 'bam'
     publishDir params.bam
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
  label 'deepvariant'
  cpus params.call_cpus
  memory params.call_mem
  errorStrategy 'finish'
  if (params.constraint)
     clusterOptions="--constraint=${params.constraint}"
  input:
     tuple val(the_id), file(bam), file(bai) from bam_file1
     file ref_seq
     file ref_fai
     each chroms
  output:
    tuple val(the_id), val(chrom), file(vcf), file(tbi), \
	  file(bam), file(bai) 
  script:
     name = bam.simpleName
     if (['X','Y','M'].contains(chrom))
        c=chrom
     else
        c = "${chrom}".padLeft(2,"0")
     vcf = "${name}-${c}-unphased.vcf.gz"
     tbi = "${vcf}.tbi"
     """
      hostname
      /opt/deepvariant/bin/run_deepvariant  \
                --model_type PACBIO \
                --ref $ref_seq  \
                --reads  $bam  \
                --output_vcf $vcf \
                --num_shards 16 \
                --regions chr$chrom

      """
}

// Now we phase the VCFs using the BAM file as input 
// Note that we have to match the BAM samples and VCF samples

process phaseVCFs {
  input:
    tuple val(name), val(chrom), file(unphased), file(unphased_tbi), \
	  file(bam), file(bai) 
     file(ref_seq)
     file(ref_fai)
  output:
    tuple val(name), val(chrom), file(phased_vcf), file("${phased_vcf}.tbi"), \
	file(bam), file(bai), emit : 'phased_vcf'
  script:
   phased_vcf = "${name}_${chrom}_phased.vcf.gz"
   """
      hostname
      whatshap phase \
        --output $phased_vcf \
        --reference  $ref_seq \
        --chromosome $chrom \
        $unphased \
        $bam
      tabix -p vcf $phased_vcf
   """
}


// Using the phased vcf run haplotag in the bam file to get the haplotagged bams
process haplotag {
  maxForks 12
  input:
    tuple val(name), val(chrom), file(phased_vcf), file(tbi), \
	  file(bam), file(bai) 
     file(ref_seq)
     file(ref_fai)
  output:
    tuple val(name), val(chrom), file(phased_bam), file("${phased_bam}.bai"),
	emit: 'haplotag'
  script:
   phased_bam = "${name}-${chrom}-haplotagged.bam"
   """
    whatshap haplotag \
        --output $phased_bam \
        --reference $ref_seq \
        $phased_vcf \
        $bam
    samtools index $phased_bam
     """
}


// We now re-run deepvariant on the phased bam file
process deepcall2 {
  label 'deepvariant'
  cpus params.call_cpus
  memory params.call_mem
  errorStrategy 'finish'
  if (params.constraint)
     clusterOptions="--constraint=${params.constraint}"
  input:
  tuple val(the_id), val(chrom), file(phased_bam), file(bai) 
     file(ref_seq)
     file(ref_fai)
  output:
    tuple val(the_id), file(phased_vcf), file(tbi), emit: 'called'
  publishDir params.vcf
  script:
     name = phased_bam.simpleName
     if (['X','Y','M'].contains(chrom))
        c=chrom
     else
        c = "${chrom}".padLeft(2,"0")
     phased_vcf = "${name}-${c}-phased.vcf.gz"
     tbi = "${phased_vcf}.tbi"
     """
      hostname
      run_deepvariant \
       --model_type PACBIO \
       --ref $ref_seq \
       --reads $phased_bam \
       --use_hp_information \
       --output_vcf $phased_vcf \
       --num_shards 16 \
       --regions ${params.chrom_prefix}$chrom
      """
}

process discover {

  input:
     tuple val(base), file(bam), file(bai)
     file(trf)
  output:
     file("${base}.svsig.gz") 
  memory 4.GB
  cpus 1
  errorStrategy 'finish'
  script:
     base = bam.simpleName
     """
       pbsv discover --tandem-repeats  $trf $bam ${base}.svsig.gz
     """
}

process pbsvcall {
  input:
     file(sigs) 
     file(ref_seq)
  output:
     tuple file("${vcf}.gz"), file("${vcf}.gz.tbi") into pbs_call_ch
  memory "12GB"
  errorStrategy 'finish'
  publishDir "${params.vcf}/"
  script:
  vcf = "${params.out}.vcf"
     """
       hostname
       /usr/bin/time -f "%e %M" pbsv call -j 8 --ccs  $ref_seq $sigs $vcf
       bgzip $vcf
       tabix ${vcf}.gz 
     """
}


workflow {
    if (params.has_bam) {
       bf = getBAMs(fnames)
    }  else {
	pbio_bamify(fnames)
        bf=pbio_bamify.out.bam
    }
    deepcall(bf,ref_seq,ref_fai,chroms)
    phaseVCFs(deepcall.out.bam,ref_seq,ref_fai)
    haplotag(phaseVCFs.out.phased_vcf, ref_seq, ref_fai)
    deepcalll(haplotag.out.haplotag, ref_seq, ref_fai)
    discover(bf,trf)
    pbsvcall(bf.out.collect(), ref_seq)
}
