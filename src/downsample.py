"""
Helper file with functions to downsample a VCF to a specified number of sites
"""

import subprocess
import random
import os.path

def get_keep_lines(vcf_path, num):
    """
    A helper function that determines how many data lines are in a vcf file (excluding header),
    downsamples to a specified number of lines, and outputs a list of the subsampled
    lines to keep 
    """
    # find out how many lines are in the file
    all_lines = subprocess.run(["wc", "-l", vcf_path], capture_output = True, text = True).stdout
    # also find out how many header rows there are
    header_lines = subprocess.run(["grep", "^#", "-c", vcf_path], capture_output = True, text = True).stdout
    header_lines = int(header_lines.strip().split(" ")[0])
    # calculate how many data lines there are in the file:
    all_lines = int(all_lines.strip().split(" ")[0]) - header_lines
    if all_lines < num:
        raise Exception(f"Error: requested downsampling to {num} lines, but only {all_lines} data lines are available in the file {vcf_path}. Please choose a smaller number of lines to downsample to.")
    print(f"all lines: {all_lines}")
    print(f"downsample to: {num} lines")
    # select num rows from the all_lines available data rows
    to_keep = random.sample(range(header_lines, all_lines), num)
    # keep all of the header lines
    to_keep = to_keep + list(range(header_lines))
    to_keep.sort()
    # to_keep is the sorted list of the lines that we're keeping
    return(to_keep)

def downsample(vcf_path, new_vcf, num):
    """
    A function that only saves a specified number 'num' of positions to the new vcf 
    These positions are randomly distributed across the input file
    """
    to_keep = get_keep_lines(vcf_path, num)
    with open(vcf_path, "r") as file1, open(new_vcf, "w") as outfile:
        line_count = 0
        for line in file1:
            # once we've written out all of the lines we're keeping, terminate
            if not to_keep:
                break
            # the first index in to_keep corresponds to the next line to keep
            elif line_count == to_keep[0]:
                line = line.strip().split("\t")
                # once we've written this line, remove it from to_keep
                to_keep.remove(line_count)
                outfile.write("\t".join(line) + "\n")
            line_count += 1   
