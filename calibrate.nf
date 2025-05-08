

process do_filter {
   maxForks 10
   input:
     path(vcf)
     each gq
     each depth
     each qual
   output:
       tuple path(out), path("${out}.csi")
   script:
       base=vcf.simpleName
       out="${base}-${gq}-${depth}-${qual}.bcf"
       """
       bcftools filter -i \
            'QUAL>=${qual} & FORMAT/DP>=${depth} & FORMAT/GQ>=${gq}' $vcf -Ob -o $out
       bcftools index $out
       """
}


process check_mendelian {
   input:
      tuple path(bcf), path(csi)
   output:
      path ("${base}.rpt")
    script:
       base = bcf.simpleName
       """
       bcftools +mendelian2 $bcf -p ${params.trio_ids}  -m c -o ${base}.rpt
       """
}

process extract_stats {
   input:
     path(report)
   output:
     path(stats)
   script:
     stats = report.simpleName+".stats"
     """
       get_calibrate_stats.py $report $stats
     """
}




process report_stats {
    input:
       path(reports)
    output:
       path(out)
    publishDir params.results_dir, mode : 'copy'
    script:
        out = params.result_name+".tsv"
        """
        echo "GQ\tDP\tQUAL\tGOOD\tBAD\tTOTAL\tMERR%\tALL_SITES" > $out
        sort -n -k 2  *stats >> $out
        """
}


process passthrough {
  input:
    tuple path(bcf), path(csi)
  output:
    tuple path(bcf), path(csi), path("${report}*{pdf,txt}")
  publishDir params.results_dir
  script:
     report = params.result_name
     """
       plot_merr.py $bcf 4_000_000 $report 
     """
}
workflow  {
	 gqs   = Channel.of(*params.gqs)
	 dps   = Channel.of(*params.depth_vals)
	 quals = Channel.of(*params.quals)
	 trio  = Channel.fromPath(params.trio)
	 main:
	   do_filter(trio,gqs,dps,quals) \
	   | check_mendelian \
	   | extract_stats \
	   | toList        \
	   | report_stats
	   do_filter.out | filter { fn ->  fn[0].name.contains("-30-10-35") }  | passthrough

}