
# Simple pipeline to call PacBio Data

This pipeline uses PacBio's smrtlink _pbmm2_ tool to align reads to a reference genome and then uses Google's _deepvariant_ tool to call, _glnexus_ to do joint calling,  and the PacBio `pbsv` tools to call SV.

## Preparation 

You need to have a FASTA file that to represent the reference genome. This is used by both of the above tools. This needs to be indexed -- and the two tools use different indexes so you have to index twice

* Index this first using _pbmm2_ e.g. `pbmm2 index ref38.fasta`   This produces a file with an _mmi_ suffix. 

* Index using samtools e.g. `samtools faidx ref38.fasta`. This produces a file with an `fai` suffix

The index files should be in the same directory as the fasta file. If necessary using symbolic links.

## Parameters

To be specified on the command line or in the config file

* `--ref_dir`: The full path to the reference genome _directory_.
* `--ref`: the name of the reference genome -- it should be in the directory. The index files should be there too
* `--fqs`: a glob with the bgzipped FastQ files -- one file per sample. You don't need to provide this is your input is a BAM file. (see `no_align` below)
* `--fast_qc`: directory where the QCed FastQ files will be placed. You don't need to provide this is your input is a BAM file. (see `no_align` below).
* `--exclude_regions`: BED files of regions which should be ignored
* `--output_vcf` directory where per sample VCFs should be placed
* `--output_gvcf` directory where gVCF files should go
* `--joint_dir  where the joinly called BCFs go
* `--joint_name` what the name of the file should be
* `--no_align` If this set to true, the BAM files provided are used as input        
* `--tandem_example`  Needed for pbsv
* `--par_regions_bed` PAR Regions
* `--bam`: Where the BAM files should be placed. Usually an output directory but see `no_align`
* `--bamify_cpus`: how many cores does creating the BAM file use (default is 16)
* `--bamify_mem`: how much memory BAM creation needs (default 32GB)
* `--call_cpus`: how many cores calling requires (default is 16)
* `--call_mem`: how much memory calling requires (default 48GB)
* `--output`: The name of the jointly called VCF file output by `pbsv`
* `--chrom_prefix`. The default is `chr`. BAM/VCF files can refer to a chromosome as chr7 or just as 7. The various tools in the pipeline need to know which. For build 38, `chr7` is more common and this is the default but if your data is different you need to set this.
* `--skip-sv`. Default is true (don't call SV)



## Wits users ##



Google's DeepVariant pipeline relies on TensorFlow which in turn relies on computers with AVX instruction support. A few of the older nodes on the 
cluster do not support AVX instructions so you need to make sure that SLURM gives you the nodes you need. Add the following options. (If you see jobs failing with a 252 error you may have overlooked this)

* `--constraint=avx2`

# calibrate.nf

Used to take a trio, run bcftools through a parameter sweep to find Mendelian error rates



