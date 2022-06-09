// (c) University of the Witwatersrand, Johannesburg, 2022
// Scott Hazelhurst
// MIT Licence





process compressTuple {
    input:
        path (vcf)
    output:
        tuple val(base), path(vcfgz), path(tbi)
    script:
        base=vcf.simpleName
        vcfgz="${vcf}.gz"
        tbi  ="${vcfgz}.tbi"
        """
        bgzip $vcf
        tabix $vcfgz
        """
}


process tabixTuple {
    input:
        path (vcf)
    output:
        tuple val(base), path(vcf), path(tbi)
    script:
        base=vcf.simpleName
        tbi  ="${vcf}.tbi"
        """
        tabix $vcf
        """
}


process splitVCF{
  time '8d'
  memory '64g'
  input:
    path (all_vcf)
  output:
    path "*indivs/*vcf.gz"
  script:
     """
     bcftools  +split   -k CHROM,POS,REF,ALT,INFO,FORMAT $all_vcf \
               -Oz -o indivs
     """
}
