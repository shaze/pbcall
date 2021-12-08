
// (c) University of the Witwatersrand, Johannesburg
// Author Scott Hazelhurst
// Released under https://creativecommons.org/licenses/by-sa/4.0/




params.output = "sample"
params.has_bam = false
trf = file(params.tandem_example)

params.chrom_prefix="chr"

fnames = params.input.tokenize(",")


ref_seq = file(params.ref_seq)
ref_fai = file(params.ref_fai)
ref_mmi = file(params.ref_mmi)



if (params.has_bam) {
  got_bams = true
  fnames.each { 
     base=file(it).baseName; 
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
   println "You promised me bams but at least some of them or their index files missing"
   System.exit(18)
  }
  globs = fname .collect { "${it}*" }.join(",")
  bf = Channel.fromFilePairs("${params.bam}/{${globs}}",size:3) { file -> file.simpleName }.map { b, f -> [b,f[0],f[1]] }
  (bam_file1, bam_file3) = bf.into (2)
} else {

  input_ch = Channel.fromPath(fnames)



   // Create a BAM file from the reads aligning to the reference genome
   // This is done per sample
   process pbio_bamify {
     cpus params.bamify_cpus
     memory params.bamify_mem
     maxForks 5
     errorStrategy 'finish'
     input:
       path(fq) from input_ch
     output:
       set val(the_id), file("${the_id}.bam"), file("${the_id}.bam.bai") into (bam_file1, bam_file2, bam_file3)
     publishDir params.bam
     script:
       the_id = fq.baseName
       """
       ls $fq/*.f*q.gz   > files.fofn
       pbmm2 align $ref_mmi  files.fofn ${the_id}.bam --preset HIFI --rg '@RG\tID:$the_id\tSM:$the_id' --unmapped --sort 
       """
   }
}

autosomes = 1..22
chroms    = autosomes +  ['X','Y','M'] 

// Do an initial calling of the variants
// Note that this is done per chromosome -- so for each sample of input we have multiple processes running
process deepcall {
  label 'deepvariant'
  cpus params.call_cpus
  memory params.call_mem
  errorStrategy 'finish'
  if (params.constraint)
     clusterOptions="--constraint=${params.constraint}"
  input:
     set val(the_id), file(bam), file(bai) from bam_file1
     file ref_seq
     file ref_fai
     each chrom from chroms
  output:
     set val(the_id), val(chrom), file(vcf), file(tbi), file(bam), file(bai) into unphased_ch
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
                --regions ${params.chrom_prefix}$chrom

      """
}

// Now we phase the VCFs using the BAM file as input 
// Note that we have to match the BAM samples and VCF samples

process phaseVCFs {
  input:
     set val(name), val(chrom), file(unphased), file(unphased_tbi), file(bam), file(bai) from unphased_ch
     file(ref_seq)
     file(ref_fai)
  output:
     set val(name), val(chrom), file(phased_vcf), file("${phased_vcf}.tbi"), file(bam), file(bai) into init_phased_vcf_ch
  script:
   phased_vcf = "${name}_${chrom}_phased.vcf.gz"
   """
      hostname
      /opt/exp_soft/python37/bin/whatshap phase \
        --output $phased_vcf \
        --reference  $ref_seq \
        --chromosome ${params.chrom_prefix}$chrom \
        $unphased \
        $bam
      tabix -p vcf $phased_vcf
   """
}


// Using the phased vcf run haplotag in the bam file to get the haplotagged bams
process haplotag {
  maxForks 12
  input:
  set val(name), val(chrom), file(phased_vcf), file(tbi), file(bam), file(bai) from init_phased_vcf_ch
     file(ref_seq)
     file(ref_fai)
  output:
     set val(name), val(chrom), file(phased_bam), file("${phased_bam}.bai") into tagged_bam_ch
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
  set val(the_id), val(chrom), file(phased_bam), file(bai)  from tagged_bam_ch
     file(ref_seq)
     file(ref_fai)
  output:
     set file(phased_vcf), file(tbi) into phased_vcf_ch
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
     set val(base), file(bam), file(bai) from bam_file3
     file(trf)
  output:
     file("${base}.svsig.gz") into svsig_ch
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
     file(sigs) from svsig_ch.collect()
     file(ref_seq)
  output:
     set file("${vcf}.gz"), file("${vcf}.gz.tbi") into pbs_call_ch
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


