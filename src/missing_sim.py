"""
Helper file with functions to add missing genotypes to a VCF at a specified rate
"""

from helper_functions import parse_header, multiallelic
import random

def missing(gt, rate):
    """
    A helper function that returns a missing genotype at a specified rate
    """
    if random.random() <= rate:
        return("./.")
    else:
        return(gt)

def add_missingness(vcf_path, new_vcf, sample_list, rate):
    """
    A function that converts a specified rate of genotypes to missing in the 
    specified samples
    """
    with open(vcf_path, "r") as file1, open(new_vcf, "w") as outfile:
        for line in file1:
            line = line.strip().split("\t")
            if line[0].startswith("#"): 
                outfile.write("\t".join(line) + "\n")
                if "#CHROM" in line[0]:
                    header_ix, include, name = parse_header(line, sample_list)
                    print(f"adding missingness to: {name}")
            else:
                # convert genotypes in the specified samples to missing at the given rate
                line = [missing(line[i], rate) if i in include else line[i] for i in range(len(line))]
                outfile.write("\t".join(line) + "\n")

