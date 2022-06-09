

process splitVCF{
  time '8d'
  memory '64g'
  input:
    path (all_vcf)
  output:
    path params.split_dir
  script:
     """
     bcftools  +split   -k CHROM,POS,REF,ALT,INFO,FORMAT $all_vcf \
               -Oz -o indivs
     """
}
