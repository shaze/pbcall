
Simplistic pipeline to call PacBio Data

This pipeline uses PacBio's smrtlink _pbmm2_ tool to align reads to a reference genome and then uses Google's _deepvariant_ tool to call and the PacBio `pbsv` tools to call.

## Preparation 

You need to have a FASTA file that to represent the reference genome. This is used by both of the above tools. This needs to be indexed --
and the two tools use different indexes so you have to index twice

* Index this first using _pbmm2_ e.g. `pbmm2 index ref38.fasta`   This produces a file with an _mmi_ suffix. 

* Index using samtools e.g. `samtools faidx ref38.fasta`. This produces a file with an `fai` suffix

The index files should be in the same directory as the fasta file. If necessary using symbolic links.

## Parameters

To be specified on the command line or in the config file

* `--input`: A comma separated list of input directories. Each directory should contain the FASTQ files used as input. No spaces in list
* `--ref_seq`: The full path to the reference genome. This should be a FASTA file
* `--ref_mmi`: The full path to an index file (_mmi_)  of the reference genome. This is produced by `pbmm2`
* `--ref_fai`: The full path to an index file (_fai_)  of the reference genome. This is produced by `samtools faidx`
* `--bam`: Where the BAM files should be placed. Usually an output directory but see `has_bam`
*  `--has_bam` (default `false`), set to to `true` if the `bam` directory above should be treated as an input
*  `--tandem_example`. A BED file with tandem repeat annotation to help discovery. Download from here https://github.com/PacificBiosciences/pbsv/blob/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed (the build 37 annotation also available)
* `--bamify_cpus`: how many cores does creating the BAM file use (default is 16)
* `--bamify_mem`: how much memory BAM creation needs (default 32GB)
* `--call_cpus`: how many cores calling requires (default is 16)
* `--call_mem`: how much memory calling requires (default 48GB)
* `--output`: The name of the jointly called VCF file output by `pbsv`
* `--chrom_prefix`. The default is `chr`. BAM/VCF files can refer to a chromosome as chr7 or just as 7. The various tools in the pipeline need to know which. For build 38, `chr7` is more common and this is the default but if your data is different you need to set this.



## Wits users ##

Please put `module load smrtlink` into your `.bashrc` file.

Google's DeepVariant pipeline relies on TensorFlow which in turn relies on computers with AVX instruction support. A few of the older nodes on the 
cluster do not support AVX instructions so you need to make sure that SLURM gives you the nodes you need. Add the following options. (If you see jobs failing with a 252 error you may have overlooked this)

* `--constraint=avx2`




