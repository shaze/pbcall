// (c) University of the Witwatersrand, Johannesburg, 2022
// Scott Hazelhurst
// MIT Licence


params.dist=1
params.support=1
params.sv_type=1
params.strand=0
params.enabled=0
params.min_size=30



process survivor {
    input:
       path vcfs
    output:
      tuple path("${combined_vcf}.gz"), path(index)
    publishDir  params.out_dir, enabled: params.publish
    script:
       combined_vcf="${params.vcf_name}.vcf"
       index = "${combined_vcf}.tbi"
       """
       ls *vcf > sample_files
       SURVIVOR merge sample_files ${params.dist} ${params.support} ${params.sv_type} \
	         ${params.strand}  ${params.enabled} ${params.min_size} $combined_vcf
       bcftools sort $combined_vcf -Oz -o ${combined_vcf}.gz    
       bcftools index -t ${combined_vcf}.gz
       """
}
    
