
// (c) University of the Witwatersrand, Johannesburg
// Author Scott Hazelhurst
// Released under https://creativecommons.org/licenses/by-sa/4.0/

ref = file(params.ref)
params.output = "sample"
container = file(params.container)
input_ch = Channel.fromPath("${params.input}")


process create_bam {
  cpus 16
  memory "32G"
  input:
    path(fq) from input_ch
  output:
    set file("${the_id}.bam"), file("${the_id}.bam.bai") into bam_file
  publishDir params.output
  script:
    the_id = fq.baseName
    """
    ls $fq/*.fa*gz   > files.fofn
    pbmm2 align $ref  files.fofn ${the_id}.bam --sort 
    """
}

process call {
  cpus 16
  memory "32G"
  input:
     """
     singularity exec -B $ref:/mnt/ --bind /usr/lib/locale/  \
                $container   \
                /opt/deepvariant/bin/run_deepvariant  \
                --model_type PACBIO \
                --ref $ref/GRCh38_no_alt_analysis_set.fasta  \\
                --reads  $bam  \
                --output_vcf ${out}.vcf.gz \
                --num_shards 16 \
                --regions chr7
      """
}
