// (c) University of the Witwatersrand, Johannesburg, 2022
// Scott Hazelhurst
// MIT Licence


params.cute_args=" --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9  --max_cluster_bias_DEL 1000  --diff_ratio_merging_DEL  0.5 "

include { survivor as survivorMultSamples } from \
    '../misc/survivor.nf'  

include { compressTuple as compress } from \
    '../misc/vcf.nf'  


ref=file(params.ref_seq)


process cuteSVs {
    cpus 16
    input:
      tuple val(b),  path(bam), path(bai)
      path(ref)
    output:
    path output
    script:
    base   = bam.simpleName
    output = "${base}.vcf"
    """
    mkdir tempwork
    cuteSV ${params.cute_args} $bam $ref $output tempwork
    """
}





workflow cuteSVAsIndividuals {
    take: bams
    main:
      cuteSVs(bams) | compress
    emit:
      compress.out 
}

workflow cuteSVCombine {
    take: bams
    main:
      cuteSVs(bams,ref) | collect | survivorMultSamples
    emit:
       survivorMultSamples.out
}

