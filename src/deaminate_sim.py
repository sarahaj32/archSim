"""
Helper file with functions to deaminate genotypes of a 
target individual. 
"""

from helper_functions import parse_header
import random

# function to create the heterozygous deaminated genome call
def deam_geno_call(base, rate):
    """
    Helper function that adds contamination to homozygous
    reference positions at a specified rate
    base is a biallelic genotype 
    rate is the rate of deamination
    simulation retains phasing 
    """
    # only change homozygous reference positions
    if base[0] != '0' or base[2] != '0':
        return base
    # if NOT deaminated (by chance), don't alter
    elif random.random() >= rate:
        return base
    # randomly select which allele is deaminated
    if base[1] == "/":
        return(f"0/1{base[3:]}")
    elif round(random.random()):
        return(f"1{base[1]}0{base[3:]}")
    else:
        return(f"0{base[1]}1{base[3:]}")

def add_deam(vcf_path, new_vcf, sample_list, rate):
    """
    converts homozygous reference transition genotypes 
    to heterozygous at specified deamination rate. 
    homozgyous alternative positions become heterozygous
    """
    with open(vcf_path, "r") as file1, open(new_vcf, "w") as outfile:
        for line in file1:
            line = line.strip().split("\t")
            if line[0].startswith("#"): 
                outfile.write("\t".join(line) + "\n")
                if "#CHROM" in line[0]:
                    header_ix, include = parse_header(line, sample_list)
                    print(f"deaminating {sample_list}")
                    print(f"include:{include}")
            else:
                # find the transition sites and add deamination
                # at all transition sites, convert homozygous reference to heterozygous 
                # what we should actually do otherwise is set a proportion of positions as transitions (ie: 77.6) and ignore the alleles themselves
                if (line[header_ix["ref_ix"]] == "C" and line[header_ix["alt_ix"]] == "T") or (line[header_ix["ref_ix"]] == "G" and line[header_ix["alt_ix"]] == "A"):
                    line = [deam_geno_call(line[i], rate) if i in include else line[i] for i in range(len(line))]
                # could remove this line to more closely mimic what happens biologically (probably not necessary if we're just trying to induce biased error)
                elif (line[header_ix["ref_ix"]] == "T" and line[header_ix["alt_ix"]] == "C") or (line[header_ix["ref_ix"]] == "A" and line[header_ix["alt_ix"]] == "G"):
                    line = [deam_geno_call(line[i], rate) if i in include else line[i] for i in range(len(line))]
                # the line may be updated by the above 2 steps or not. Regardless, write it out
                outfile.write("\t".join(line) + "\n")
