//  (C) University of the Witwatersrand, Johannesburg, 2022
//  Scott   Hazelhurst



include { tabixTuple as tabix; splitVCF } from \
    '../misc/vcf.nf'  


ref=file(params.ref_seq)


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
     tuple file("${vcf}.gz"), file("${vcf}.gz.tbi") 
  memory "12GB"
  errorStrategy 'finish'
  script:
  vcf = "${params.out}.vcf"
     """
       hostname
       /usr/bin/time -f "%e %M" pbsv call -j 8 --ccs  $ref_seq $sigs $vcf
       bgzip $vcf
       tabix ${vcf}.gz 
     """
}


workflow pbsv {
    take: bams
    take: trf
    main:
       discover(bams,trf) | pbsvcall 
    emit:
       discover.out
}


workflow pbsvIndividuals {
    take: bams
    take: trf
    main:
       discover(bams,trf) | pbsvcall | splitVCF | tabix
    emit:
       tabix.out
}

