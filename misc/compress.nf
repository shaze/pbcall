// (c) University of the Witwatersrand, Johannesburg, 2022
// Scott Hazelhurst
// MIT Licence





process compress {
    input:
        file(vcf)
    output:
        tuple val(base), file(vcfgz), file(tbi)
    script:
        base=vcf.simpleName
        vcfgz="${vcf}.gz"
        tbi  ="${vcfgz}.tbi"
        """
        bgzip $vcf
        tabix $vcfgz
        """
}
