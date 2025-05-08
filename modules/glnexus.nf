


chrom_size_b38 = ['01':248956422, '02':242193529, '03':198295559, '04':190214555,
	      '05':181538259, '06':170805979, '07':159345973, '08':145138636,
	      '09':138394717, '10':133797422, '11':135086622, '12':133275039,
	      '13':114364328, '14':107043718, '15':101991180, '16':90338345,
	      '17':83257441,  '18':80373285,  '19':50617616,  '20':64444167,
	      '21':47709983,  '22':50818468,  'X' :156040896, 'Y' :57227415,
	      'M':17000]

chrom_size_t2t = ['01': 248387328, '02': 242696752, '03': 201105948, '04': 193574945,
	      '05': 182045439, '06': 172126628, '07': 160567428, '08': 146259331,
	      '09': 150617247, '10': 134758134, '11': 135127769, '12': 133324548,
	      '13': 113566686, '14': 101161492, '15': 99753195, '16': 96330374,
	      '17': 84276897, '18': 80542538, '19': 61707364, '20': 66210255,
	      '21': 45090682, '22': 51324926, 'X': 154259566, 'Y': 62460029, 'M': 16569]


chrom_size_b37= ['01' :  249250621, '02' : 243199373, '03' : 198022430, '04' : 191154276,
		 '05' : 180915260, '06' : 171115067, '07' : 159138663, '08' : 146364022,
		 '09' : 141213431, '10' : 135534747, '11' : 135006516, '12' : 133851895,
		 '13' : 115169878, '14' : 107349540, '15' : 102531392, '16' : 90354753 ,
		 '17' : 81195210 , '18' : 78077248 , '19' : 59128983 , '20' : 63025520 ,
		 '21' : 48129895 , '22' : 51304566 , 'X' : 155270560, 'Y' : 59373566 ]


if (params.build=='t2t') {
   chrom_size = chrom_size_t2t
} else if (params.build='b37') {
  chrom_size = chrom_size_b37 
} else if (params.build='b38') {
  chrom_size = chrom_size_b38
} else {
   println "No such build : ${params.build}"
   System.exit(10)
}

/* for v large scale
def memoryForChrom ( the_chrom ) {
    chrom =  the_chrom.getName()
    req   = Math.round(chrom_size[chrom]*2.96/1000000+246)
    return "${req}.G"
}

*/

if (params.scratch=="false")
    scratch=false;
else if (params.scratch=="true")
    scratch==true;
else
    scratch=params.scratch;


process GLnexus {
    scratch scratch
    memory "64G"
    cpus   24
    input:
      tuple val(chrom), path(vcfs), path(tbis)
    output:
       path(bcf)
    script:
       bcf="${params.joint_name}-${chrom}.bcf"
       res="dv${chrom}.t"
       """
       hostname 
       set limit descriptors 3600
       /usr/bin/time -f "Elapsed time: %e; MAXRSS=%M" -o $res\
            glnexus_cli  --config DeepVariantWGS $vcfs > $bcf
       """
}

