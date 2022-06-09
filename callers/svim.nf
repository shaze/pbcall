// (c) University of the Witwatersrand, Johannesburg, 2022
// Scott Hazelhurst
// MIT Licence


params.svim_args=" --min_sv_size 40 "

include { compressTuple as compress } from \
    '../misc/vcf.nf'  


ref=file(params.ref_seq)




process svimAlign {
    cpus 1
    input:
      tuple val(b),  path(bam), path(bai)
      path(ref)
    output:
    path output
    script:
    base   = bam.simpleName
    output = "work/${base}.vcf"
    """
    svim alignment  ${params.svim_args} work $bam $ref 
    """
}



workflow svim {
    take: bams
    main:
      svimAlign(bams) | compress
    emit:
      compress.out 
}
