#!/usr/bin/env python3
"""
BCF MERR Analyzer

This program analyzes a BCF file (binary VCF) to determine the percentage of positions
with MERR>0 among positions where F_MISSING=0. It generates plots by chromosome
using a kernel density estimation approach and a summary report.

Usage:
    python bcf_merr_analyzer.py <bcf_file> <width> <output_base>

Arguments:
    bcf_file: Path to the BCF file
    width: Bandwidth parameter for the kernel density estimation (integer)
    output_base: Base name for output files
"""

import sys
import os
import re
import matplotlib.pyplot as plt
import numpy as np
import pysam
import math
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker

import argparse

centromere_38= {
    "chr1": (122026459,124932724),
    "chr10": (39686682,41593521),
    "chr11": (51078348,54425074),
    "chr12": (34769407,37185252),
    "chr13": (16000000,18051248),
    "chr14": (16000000,18173523),
    "chr15": (17083673,19725254),
    "chr16": (36311158,38265669),
    "chr17": (22813679,26616164),
    "chr18": (15460899,20861206),
    "chr19": (24498980,27190874),
    "chr2": (92188145,94090557),
    "chr20": (26436232,30038348),
    "chr21": (10864560,12915808),
    "chr22": (12954788,15054318),
    "chr3": (90772458,93655574),
    "chr4": (49712061,51743951),
    "chr5": (46485900,50059807),
    "chr6": (58553888,59829934),
    "chr7": (58169653,61528020),
    "chr8": (44033744,45877265),
    "chr9": (43389635,45518558),
    "chrX": (58605579,62412542),
    "chrY": (10316944,10544039)}

centromere_t2t= {
    "chr1": (122026459,124932724),
    "chr10": (39686682,41593521),
    "chr11": (51078348,54425074),
    "chr12": (34769407,37185252),
    "chr13": (16000000,18051248),
    "chr14": (16000000,18173523),
    "chr15": (17083673,19725254),
    "chr16": (36311158,38265669),
    "chr17": (22813679,26616164),
    "chr18": (15460899,20861206),
    "chr19": (24498980,27190874),
    "chr2": (92188145,94090557),
    "chr20": (26436232,30038348),
    "chr21": (10864560,12915808),
    "chr22": (12954788,15054318),
    "chr3": (90772458,93655574),
    "chr4": (49712061,51743951),
    "chr5": (46485900,50059807),
    "chr6": (58553888,59829934),
    "chr7": (58169653,61528020),
    "chr8": (44033744,45877265),
    "chr9": (43389635,45518558),
    "chrX": (58605579,62412542),
    "chrY": (10316944,10544039)}





def parse_arguments():
    # Create the parser
    parser = argparse.ArgumentParser(description='Process input files and parameters.')
    
    # Add the optional argument
    parser.add_argument('--compare', type=str, default="", help='Optional file name for comparison')
    
    # Add the three compulsory arguments
    parser.add_argument('bcf', type=str, help='Input file name')
    parser.add_argument('width', type=int, help='Integer value')
    parser.add_argument('base', type=str, help='String value')
    
    # Parse the arguments
    args = parser.parse_args()
    return args

def process_bcf_file(bcf_file, width, output_base):
    """
    Process the BCF file, analyze MERR values, and generate reports.
    
    Args:
        bcf_file: Path to the BCF file
        width: Interval size for binning
        output_base: Base name for output files
    """
    # Open the BCF file
    try:
        vcf_reader = pysam.VariantFile(bcf_file)
    except Exception as e:
        print(f"Error opening BCF file: {e}")
        sys.exit(1)
    
    # Dictionary to store data by chromosome
    chroms_data = {}

    old_chrom="00"
    # Process each record in the BCF file
    for record in vcf_reader:
        chrom = record.chrom
        if chrom != old_chrom:
            old_chrom=chrom
            #if chrom=="chr2":break
        pos = record.pos
        
        # Extract INFO fields
        info = record.info

        # Skip if F_MISSING is not 0
        if 'F_MISSING' not in info or info['F_MISSING'][0] != 0:
            continue
        
        # Get MERR value
        merr_value = info.get('MERR', 0)

        # Initialize chromosome data if not already present
        if chrom not in chroms_data:
            chroms_data[chrom] = {
                'positions': [],
                'merr_status': []
            }
        
        # Add position and MERR status
        chroms_data[chrom]['positions'].append(pos)
        chroms_data[chrom]['merr_status'].append(1 if merr_value > 0 else 0)
    
    vcf_reader.close()
    
    # Create reports and plots
    return create_reports_and_plots(chroms_data, width, output_base, width)


def create_reports_and_plots(chroms_data, width, output_base, window_width=None):
    """
    Create reports and plots from the processed data.
    
    Args:
        chroms_data: Dictionary containing data by chromosome
        width: Interval size for binning
        output_base: Base name for output files
        window_width: Width of the sliding window used for plots
    """
    # Prepare report data
    report_data = []
    total_positions = 0
    total_merr_positions = 0
    
    # Sort chromosomes naturally
    def natural_key(s):
        return [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', s)]
    
    chromosomes = sorted(chroms_data.keys(), key=natural_key)

    result=[]
    # Process each chromosome
    for chrom in chromosomes:
        positions = chroms_data[chrom]['positions']
        merr_status = chroms_data[chrom]['merr_status']
        
        if not positions:
            continue
        
        # Count statistics
        num_positions = len(positions)
        num_merr_positions = sum(merr_status)
        merr_percentage = (num_merr_positions / num_positions) * 100 if num_positions > 0 else 0
        
        # Add to totals
        total_positions += num_positions
        total_merr_positions += num_merr_positions
        
        # Add to report data
        report_data.append({
            'chrom': chrom,
            'total_positions': num_positions,
            'merr_positions': num_merr_positions,
            'merr_percentage': merr_percentage
        })
        
        # Create plot for this chromosome
        result.append(create_chromosome_plot(chrom, positions, merr_status, width, output_base))
    
    # Calculate overall percentage
    overall_percentage = (total_merr_positions / total_positions) * 100 if total_positions > 0 else 0
    
    # Write summary report
    write_summary_report(report_data, total_positions, total_merr_positions, overall_percentage, output_base, window_width)
    return result


def create_chromosome_plot(chrom, positions, merr_status, width, output_base):
    """
    Create a plot for a chromosome showing the percentage of MERR>0 by position
    using a kernel density estimation approach.
    
    Args:
        chrom: Chromosome name
        positions: List of positions
        merr_status: List of MERR status (1 for MERR>0, 0 otherwise)
        width: Width parameter for the smoothing kernel
        output_base: Base name for output files
    """
    if not positions:
        return
    
    # Convert to numpy arrays for better performance
    positions = np.array(positions)
    merr_status = np.array(merr_status)
    
    # Sort positions and corresponding MERR status


    x_pos=[]
    curr=prev=0
    num_pos=num_errs=0
    pct_err=[]
    step=width//32
    for bp in range(step,positions[-1],step):
        while curr<len(positions) and positions[curr]<bp:
            num_pos=num_pos+1
            if merr_status[curr]: num_errs=num_errs+1
            curr=curr+1
        while bp>width and bp-positions[prev]>width:
            num_pos=num_pos-1
            if merr_status[prev]: num_errs=num_errs-1
            prev=prev+1
        if num_pos>0:
            x_pos.append(bp/1_000_000)
            pct_err.append(num_errs/num_pos*100)
    plt.figure(figsize=(12, 6))
    plt.yscale('log')
    formatter = ScalarFormatter(useOffset=False)
    formatter.set_scientific(False)
    ax = plt.gca()
    ax.axhline(y=1, color='lightgreen', linestyle='-', alpha=0.7)
    yticks=[0.1,0.2,0.4,0.6,1,2,5,10]
    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{y:.1f}" for y in yticks])
    ax.set_xlim(left=0, right=positions[-1]//1_000_000+1)
    ax.set_ylim(bottom=None, top=11)
    ax.yaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    formatter_x = ScalarFormatter(useOffset=False)
    formatter_x.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter_x)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
    plt.plot(x_pos,pct_err, 'o',markersize=1.8)
    median=np.median(pct_err)
    plt.xlabel('Position (Mb)')
    plt.ylabel('Percentage with MERR > 0 (%)')
    plt.title(f'{output_base} Chromosome {chrom}: Percentage of Positions with MERR > 0 : Median MERR is {median:.2f}')
    plt.grid(True, alpha=0.3)
    
    # Format the x-axis with appropriate scale
    #plt.ticklabel_format(axis='x', style='scientific', scilimits=(0, 0))
    
    # Pad the chromosome number with leading zero if needed
    chrom_num = re.sub(r'\D', '', chrom)  # Extract numeric part
    if chrom_num:
        padded_chrom = chrom_num.zfill(2)
    else:
        padded_chrom = chrom  # Keep original if no numeric part
    
    # Save the plot
    plt.tight_layout()
    plt.savefig(f"{output_base}-{padded_chrom}.pdf")
    plt.close()
    return (x_pos,pct_err)

def write_summary_report(report_data, total_positions, total_merr_positions, overall_percentage, output_base, width=None):
    """
    Write a summary report with statistics by chromosome.
    
    Args:
        report_data: List of dictionaries with chromosome statistics
        total_positions: Total number of positions across all chromosomes
        total_merr_positions: Total number of positions with MERR>0
        overall_percentage: Overall percentage of positions with MERR>0
        output_base: Base name for output files
        width: Width of the sliding window used for analysis
    """
    with open(f"{output_base}.report.txt", "w") as f:
        f.write("BCF MERR Analysis Report\n")
        f.write("======================\n\n")
        
        f.write("Summary by Chromosome:\n")
        f.write("-" * 60 + "\n")
        f.write(f"{'Chromosome':<12} {'Total Positions':<18} {'MERR>0 Positions':<18} {'Percentage (%)':<15}\n")
        f.write("-" * 60 + "\n")
        
        for data in report_data:
            f.write(f"{data['chrom']:<12} {data['total_positions']:<18} {data['merr_positions']:<18} {data['merr_percentage']:.2f}\n")
        
        f.write("-" * 60 + "\n")
        f.write(f"{'TOTAL':<12} {total_positions:<18} {total_merr_positions:<18} {overall_percentage:.2f}\n")
        f.write("-" * 60 + "\n\n")
        
        f.write("Analysis Parameters:\n")
        f.write(f"- Only positions with F_MISSING=0 were considered\n")
        f.write(f"- Percentage represents the proportion of positions with MERR>0\n")
        if width:
            f.write(f"- Window of={width} was used for the plots\n")

def show_comparison(first,second):
    chrom=0
    tot_num_a=tot_errs_a=tot_num_b=tot_errs_b=0
    for (a,b) in zip(first,second):
        chrom=chrom+1
        pos_a,a_val=a
        pos_b,b_val=b
        a_only=[]
        b_only=[]
        common=[]
        ratio =[]
        i_a=i_b=0
        while i_a<len(pos_a) or i_b<len(pos_b):
           if i_b>=len(pos_b) or i_a<len(pos_a) and pos_a[i_a]<pos_b[i_b]:
               a_only.append(pos_a[i_a])
               i_a=i_a+1
           elif i_a>=len(pos_a) or i_b<len(pos_b) and  pos_b[i_b]<pos_a[i_a]:
               b_only.append(pos_b[i_b])
               i_b=i_b+1
           else:
               diff=a_val[i_a]-b_val[i_b]       
               ratio.append(diff)
               common.append(pos_a[i_a])
               i_a=i_a+1
               i_b=i_b+1
        plt.figure(figsize=(12, 6))
        plt.plot(common, ratio, 'o',markersize=1)
        plt.plot(a_only, [-8]*len(a_only), 'o',markersize=1,color='r')
        plt.plot(b_only, [8]*len(b_only), 'o',markersize=8,color='g')                
        plt.xlabel('Position')
        plt.ylabel('MERR diff (Negative T2T better, Positive GRCh38 better, Red: T2T only')
        plt.title(f'Chromosome {chrom}: T2T versus GrCh38')
        plt.grid(True, alpha=0.3)

        # Format the x-axis with appropriate scale
        plt.ticklabel_format(axis='x', style='scientific', scilimits=(0, 0))
        # Pad the chromosome number with leading zero if needed
        chrom_num = str(chrom)
        if chrom_num:
            padded_chrom = chrom_num.zfill(2)
        else:
            padded_chrom = chrom  # Keep original if no numeric part
        plt.tight_layout()
        plt.savefig(f"comparison-{padded_chrom}.pdf")
        plt.close()


            
def main():

    args = parse_arguments()
    
    bcf_file = args.bcf
    
    try:
        width = args.width
        if width <= 0:
            raise ValueError("Bandwidth parameter must be a positive integer")
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    
    output_base = args.base
    
    # Process the BCF file
    first=process_bcf_file(bcf_file, width, output_base)

    if args.compare:
       second=process_bcf_file(args.compare, width, "38")
       show_comparison(first,second)
       


if __name__ == "__main__":
    main()
    

    
