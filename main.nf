
// (c) University of the Witwatersrand, Johannesburg
// Author Scott Hazelhurst
// Released under https://creativecommons.org/licenses/by-sa/4.0/




params.output = "sample"
params.has_bam = false

include { GLnexus } from "./modules/glnexus.nf"
include { index; norm; bcf_filter; bcf_exclude_regions  } from "./modules/qc.nf"


params.chrom_prefix="chr"


stats_report =  params.vcf_stats_report ? " --vcf_stats_report=true " : ""

def get_bams (fnames) {
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
  bf = Channel.fromFilePairs("${params.bam}/${globs}",size:2) { file -> file.simpleName }.map { b, f -> [b,f[0],f[1]] }
  return bf
}



// Create a BAM file from the reads aligning to the reference genome
// This is done per sample
process fastp {
   cpus 12
   input:
      path(fq)
   output:
      path("output/${qc_fq}.fq.gz"), emit: fq
      path("output/*.{json,html}"), emit:reports
   if (params.fast_qc) {
     storeDir params.fast_qc
   }
   script:
      qc_fq=fq.simpleName
      """
      mkdir -p output
      fastplong -w ${task.cpus} -e ${params.min_read_fq_qc} -i $fq/*f*q.gz \
                     -j output/${qc_fq}.json -h output/${qc_fq}.html -o output/${qc_fq}.fq.gz
      """
}

process pbio_bamify {
   cpus params.bamify_cpus
   memory params.bamify_mem
   maxForks 5
   errorStrategy 'finish'
   input:
      tuple path(ref), path(ref_mmi), path(fq)
   output:
      tuple val(the_id), path("${the_id}.bam"), path("${the_id}.bam.bai") 
   publishDir params.bam, mode: 'copy'
   script:
   the_id = fq.simpleName
   """
   hostname
   pbmm2 align ${ref_mmi}  $fq  ${the_id}.bam -j ${task.cpus} \
           --preset HIFI --rg '@RG\\tID:$the_id\\tSM:$the_id' --unmapped --sort 
   """
}

process create_index {
   cpus 2
   memory '24.GB'
   input:
      path(ref)
   output:
      tuple path(ref), path(refi)
   storeDir params.ref_dir
   script:
     base=ref.baseName
     refi=base+".mmi"
     """
        hostname
        pbmm2 index $ref $refi --preset HIFI
     """
}

workflow make_bams {
   take:
      fastq
   main:
      fastp(fastq)
      fq_qc = fastp.out.fq
      ref=Channel.fromPath("${params.ref_dir}/${params.ref}")
      create_index(ref) \
      | combine(fq_qc)    \
      | pbio_bamify
  emit:
     pbio_bamify.out
}



autosomes = 1..22
chroms    = autosomes +  ['X','Y','M'] 

// Do an initial calling of the variants
// Note that this is done per chromosome --
// so for each sample of input we have multiple processes running
process deepcall {
  label 'deepvariant'
  cpus params.call_cpus
  memory params.call_mem
  errorStrategy 'finish'
  if (params.constraint)
     clusterOptions="--constraint=${params.constraint}"
  input:
     tuple  path(ref_dir), val(the_id), path(bam), path(bai)
     each chrom
  output:
     tuple val(c), path(vcf), path(tbi)
  publishDir params.output_vcf
  script:
     par_regions  =  params.par_regions_bed ? \
                  " --par_regions_bed=${ref_dir}/${params.par_regions_bed} " : ""
     name = bam.simpleName
     if (['X','Y','M'].contains(chrom))
        c=chrom
     else
        c = "${chrom}".padLeft(2,"0")
     vcf = "${name}-${c}-unphased.vcf.gz"
     tbi = "${vcf}.tbi"
     if (chrom=="M") {
        shards=1;
	extra= ' --call_variants_extra_args="allow_empty_examples=true" '
     } else {
        shards=16;
	extra=""
     }
     gvcf = params.output_gvcf ? " --output_gvcf ${name}-${c}.gvcf.gz " : ""
     
     """
      hostname
      /opt/deepvariant/bin/run_deepvariant  \
                --model_type PACBIO \
                $extra \
                --ref $ref_dir/${params.ref}  \
                --reads  $bam  \
                $gvcf \
                $par_regions \
                $stats_report \
                --output_vcf $vcf \
                --num_shards $shards \
                --regions chr$chrom
      """
}


process discover {

  input:
     tuple path(trf), val(base), path(bam), path(bai) 
  output:
     path("${base}.svsig.gz") 
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
  cpus 12
  input:
     tuple path(ref_dir), path (data)
  output:
     tuple  file("${vcf}.gz"), file("${vcf}.gz.tbi")
  memory "12GB"
  errorStrategy 'finish'
  publishDir "${params.output_vcf}/", pattern: ".vcf.gz", mode: 'copy'
  publishDir "${params.output_gvcf}/", pattern: ".gvcf.gz", enabled: params.output_gvcf 
  script:
    base = data.simpleName
    vcf = "${base}.vcf"
     """
       hostname
       /usr/bin/time -f "%e %M" pbsv call -j 12 --ccs ${ref_dir}/${params.ref}  *svsig.gz  $vcf
       sed -i  '/\\x00/d'  $vcf
       bgzip $vcf
       tabix ${vcf}.gz 
     """
}

workflow {
  fqs = Channel.fromPath(params.fqs,type:'dir')
  trf        = Channel.fromPath(params.tandem_example)
  ref_ch     = Channel.fromPath(params.ref_dir,type:"dir")
  exclude_ch = Channel.fromPath(params.exclude_regions)
  bf = Channel.fromFilePairs("${params.bam}/*.bam*",size:2) { file -> file.simpleName }.map { b, f -> [b,f[0],f[1]] }
  main:
    bam_ch = params.no_align ? bf : make_bams(fqs)
    samples = ref_ch.combine(bam_ch)
    if (!params.skip_sv) {
      trf.combine(bam_ch) | discover
      ref_ch.combine(discover.out) | pbsvcall
    }
    deepcall(samples,chroms)
       | groupTuple
       | GLnexus
       | bcf_filter
    ref_ch.combine(bcf_filter.out) | norm        | combine(exclude_ch)
       | bcf_exclude_regions

}
