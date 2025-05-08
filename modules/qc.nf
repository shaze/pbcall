

process bcf_exclude_regions {
   input:
       tuple path(bcf), path(exclude_regions)
   output:
       tuple path(out), path(outcsi)
   publishDir params.joint_dir, mode: 'copy'
   script:
       out=bcf.simpleName.replaceAll("-norm","-clean")+".bcf"
       outcsi=out+".csi"
       """
       bcftools view $bcf -T^${exclude_regions} -Ob -o $out
       bcftools index $out
       """
}


process bcf_filter {
   input:
      file(bcf)
   output:
      file(filterbcf)
   script:
      filterbcf="${bcf.simpleName}-filter.bcf"
   """
     hostname
     bcftools filter -i 'QUAL>=${params.qual} & FORMAT/DP>=${params.depth} & FORMAT/GQ>=${params.gq}' $bcf  -o $filterbcf
   """
}

process bcf_filter_p {
   input:
      val(flags)
      file(bcf)
   output:
      file(filterbcf)
   script:
      filterbcf="${bcf.simpleName}-filter.bcf"
   """
     hostname
     bcftools filter -i 'QUAL>=${flags['qual']} & FORMAT/DP>=${flags['depth']} & FORMAT/GQ>=${flags['gq']}' $bcf  -o $filterbcf
   """
}

process norm {
  cpus 4
    input:
      tuple path(refdir), path(bcf)
   output:
      file(normbcf)
   script:
      normbcf="${bcf.simpleName}".replace("-filter","-norm")+".bcf"
   """
     bcftools norm --threads 3 --fasta-ref  $refdir/${params.ref}  $bcf  -o $normbcf
   """
}

process index {
   cpus 3
   input:
      file(bcf)
   output:
      file(index)
   publishDir params.output_dir
   script:
      index="${bcf}.csi"
   """
     bcftools index --threads 3  $bcf
   """
}
   