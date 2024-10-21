"""
Helper file with functions that randomly samples read depths from a given distribution 
for specified individuals. Does not alter the genotypes of other individuals.
Also contains an option to remove
"""

from helper_functions import parse_header
import numpy as np
import random

def get_depth(mean, var):
    """
    A helper function that samples depth from a given normal distribution,
    with a given mean and variance.
    This is the total reads for this site, is distributed across alleles in later functions
    I should change so that it's not only normal
    """

    max(round(np.random.normal(mean, var)), 0)
    depth1 = depth2 = depth
    return depth1, depth2

def pos_depth_all(pos, mean, var):
    """
    A function that adds depth AND
    1. outputs missing genotypes if depth is 0
    2. outputs homozygous genotype if depth is 3 or less
        randomly selects allele
    3. randomly disperses 3 or more reads across alelels  
    """

    depth1, depth2 = get_depth(mean, var)

    # add missing genotype if there are no reads
    if depth == 0:
        gt = "./."
    
    # homozygous: all reads support 1 allele 
    elif pos[0] == pos[2]:
        gt = pos[0:3]
        if pos[0] == '0':
            depth2 = 0
        else:
            depth1 = 0

    # heterozygous: induce false homozygous if there are 3 or fewer reads
    elif depth <= 3:
        # randomly select the allele to be homozygous
        allele = random.choice(pos[0], pos[2]) # could adjust this line to induce a reference bias
        gt = f"{allele}{pos[1]}{allele}"
        if allele == 0:
            depth2 = 0
        else:
            depth1 = 0
    # if there are more than 3 reads, randomly distribute the total reads between alleles
    else:
        gt = pos[0:3]
        # ensure that there is at least 1 read supporting each allele - but randomly distribute
        depth1 = random.randint(1, (depth - 1))
        depth2 = depth - depth1 
    # output genotype, with depth information
    return(f"{gt}{pos[3:]}:{str(depth1)},{str(depth2)}:{str(depth)}")

def pos_depth_only(pos, mean, var):
    """
    A function that adds depth ONLY
    *does not alter any genotypes*
    randomly disperses 3 or more reads across alleles  
    """

    depth1, depth2 = get_depth(mean, var)

    # homozygous: all reads support 1 allele     
    if pos[0] == pos[2]:
        if pos[0] == '0':
            depth2 = 0
        else:
            depth1 = 0

    # randomly distribute the total reads between alleles at heterozygous sites
    else:
        # ensure that there is at least 1 read supporting each allele - but randomly distribute
        depth1 = random.randint(1, max(1, (depth - 1)))
        depth2 = depth - depth1 

    return(f"{pos}:{str(depth1)},{str(depth2)}:{str(depth)}")

def pos_depth_homo(pos, mean, var):
    """
    A function that adds depth AND
    1. outputs homozygous genotype if depth is 3 or less
        randomly selects allele
    2. randomly disperses 3 or more reads across alleles
    *does not output missing data, even if depth = 0*
    """

    depth1, depth2 = get_depth(mean, var)

    # homozygous: all reads support 1 allele 
    if pos[0] == pos[2]:
        gt = pos[0:3]
        if pos[0] == '0':
            depth2 = 0
        else:
            depth1 = 0

    # heterozygous: induce false homozygous if there are 3 or fewer reads
    elif depth <= 3:
        # randomly select the allele to be homozygous
        allele = random.choice(pos[0], pos[2]) # could adjust this line to induce a reference bias
        gt = f"{allele}{pos[1]}{allele}"
        if allele == 0:
            depth2 = 0
        else:
            depth1 = 0

    # if there are more than 3 reads, randomly distribute the total reads between alleles
    else:
        gt = pos[0:3]
        # ensure that there is at least 1 read supporting each allele - but randomly distribute
        depth1 = random.randint(1, (depth - 1))
        depth2 = depth - depth1 

    return(f"{gt}{pos[3:]}:{str(depth1)},{str(depth2)}:{str(depth)}")


def add_depth(vcf_path, new_vcf, sample_list, mean, var, distribution, fun_name):
    """
    function that 
    """
    exclude = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

    with open(vcf_path, "r") as file1, open(new_vcf, "w") as outfile:

        # This is the function to apply: either pos_depth_all, pos_depth_only, pos_depth_homo
        dp_fun = globals()[fun_name]
        print(dp_fun)
        for line in file1:
            line = line.strip().split("\t")
            if line[0].startswith("#"):
                if "#CHROM" in line[0]:
                    # add on header with FORMAT information (so that the depth we add can still be parsed)
                    outfile.write('##FORMAT=<ID=AD,Number=1,Type=String,Description="allele_depth">\n')
                    outfile.write('##FORMAT=<ID=DP,Number=1,Type=String,Description="total depth">\n')
                    header_ix, include = parse_header(line, sample_list)
                    exclude_ix = [i for i,name in enumerate(line) if name in exclude]
                    outfile.write("\t".join(line) + "\n")

                else:
                    outfile.write("\t".join(line) + "\n")

            else:
                # this allows us to add AD and DP to all lines:
                metrics = [line[i] for i in range(len(line)) if i in exclude_ix]
                metrics[header_ix["format_ix"]] += ":AD:DP"
                # if line[header_ix["alt_ix"]] != ".": # this would allow us to only focus on variable sites - probably not a good idea

                # add in depth for every line - but only to the specified individuals
                line = [dp_fun(line[i], mean, var) if i in include else (line[i] + "::") for i in range(len(line))]
                # remove "metrix" from the line (those don't need any depth annotations)
                line = [line[i] for i in range(len(line)) if i not in exclude_ix]
                # then recombine
                line = metrics + line
                
                outfile.write("\t".join(line) + "\n")
