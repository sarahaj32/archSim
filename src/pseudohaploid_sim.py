"""
Helper file with functions to randomly sample an allele
from heterozygous positions, outputting pseudo-haploid genotype calls
"""
from helper_functions import parse_header
import random

def ph_geno_call(base):
    """
    helper function that randomly selects an allele
    from the given genotypes
    """
    if base[0] == base[2]:
        return base
    else:
        # randomly select one of the two alleles, and make homozygous for that allele
        call = random.choice([base[0], base[2]])
        return(str(call) + base[1] + str(call) + base[3:])

def make_pseudohaploid(vcf_path, new_vcf, sample_list):
    """
    Finds the individuals to be pseudohaplotyped,
    at heterozygous sites, randomly selects an allele
    and converts genotype to homozygous for that allele
    """
    with open(vcf_path, "r") as file1, open(new_vcf, "w") as outfile:
        for line in file1:
            line = line.strip().split("\t")
            if line[0].startswith("#"): 
                outfile.write("\t".join(line) + "\n")
                if "#CHROM" in line[0]:
                    # include: the indices of the samples for the simulation to be applied to
                    _, include = parse_header(line, sample_list)
                    print(f"making {sample_list} samples pseudohaploid")
            else:
                # convert to pseudohaploid if the individual is specified to be pseudohaploid
                line = [ph_geno_call(line[i]) if i in include else line[i] for i in range(len(line))]
                outfile.write("\t".join(line) + "\n")
