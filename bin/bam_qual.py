#!/usr/bin/env python3

import sys


import pysam

input = sys.argv[1]
min_mapq= int(sys.argv[2])
output  = sys.argv[3]

def mark_low_mapq_unmapped(input_bam, output_bam, min_mapq=30):
    """Mark reads with low MAPQ as unmapped and sort the output"""
    
    # Read input BAM
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        # Get header for output
        header = infile.header
        
        # Temporary unsorted BAM
        temp_bam = output_bam + ".tmp"
        
        with pysam.AlignmentFile(temp_bam, "wb", header=header) as outfile:
            for read in infile:
                if read.mapping_quality < min_mapq:
                    # Mark as unmapped
                    read.is_unmapped = True
                    read.reference_id = -1
                    read.reference_start = -1
                    read.mapping_quality = 0
                    
                    # Optionally preserve original data in tags
                    if not read.is_unmapped:  # If it was mapped before
                        read.set_tag("XM", read.reference_name)  # Original reference
                        read.set_tag("XP", read.reference_start)  # Original position
                        read.set_tag("XQ", read.mapping_quality)  # Original MAPQ
                
                outfile.write(read)
    
    # Sort the output BAM
    pysam.sort("-@","12", "-o", output_bam, temp_bam)
    
    # Index the sorted BAM
    pysam.index(output_bam)
    
    # Clean up temporary file
    import os
    os.remove(temp_bam)

# Usage
mark_low_mapq_unmapped(input, output, min_mapq=min_mapq)
