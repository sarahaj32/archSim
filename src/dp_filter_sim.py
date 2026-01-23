"""
Helper file with functions that randomly samples read depths for a sample from a negative binomial distribution,
then distributes those reads across alleles following binomial distribtuion, creating missing data for sites
that no reads support, and inducing homozygous genotypes when depth is below a given threshold with the allele
being selected with the given reference bias. 
"""

from helper_functions import parse_header
import numpy as np
import random

# pnina's updated function that allows for negative binomial sampling
def get_sample_depth(mean, var, dist="nbinom"):
    """
    A helper function that samples DP for a sample from 
    a nbinomial distribution with a given mean and variance.
    This is the total reads for this site, and is distributed 
    across alleles in later functions.
    """
    m = max(0.0, float(mean))
    v = float(var)
    if dist == "nbinom":
        # v = m + m^2/size  → size = m^2 / (v - m)
        if v <= m or m <= 0:
            dist = "poisson"
        else:
            size = (m*m) / max(v - m, 1e-9)
            p = size / (size + m)
            return int(np.random.negative_binomial(size, p))
    if dist == "poisson":
        return int(max(0, np.random.poisson(lam=m)))
    if dist == "normal":
        sd = math.sqrt(v if v is not None else m + 1.0)
        return max(0, int(round(np.random.normal(loc=m, scale=sd))))
    raise ValueError("depth_dist must be 'nbinom','poisson','normal'")

def pos_depth(pos, mean, var, ref_bias, dropout_depth, dist):
    """
    A function that adds depth AND
    1. adds missing data if no reads are sampled for a given position and individual
    2. outputs homozygous genotype if depth is below the dropout depth or less
        randomly selects allele with the provided reference bias
    3. at heterozygous sites disperses 3 or more reads across alleles following a binomial distribution.
    If no reads are assigned to an allele, converts to homozygous
    """

    dp = get_sample_depth(mean, var, dist)
    a1=pos[0]
    a2=pos[2]
    sep=pos[1]
    # if no reads, output missing genotype
    if dp <= 0:
        return "./.:.:0"

    # homozygote: all reads support 1 allele 
    if a1 == a2:
        depth1 = depth2 = dp
        gt = pos[0:3]
        if a1 == '0':
            depth2 = 0
        else:
            depth1 = 0
        return(f"{gt}:{str(depth1)},{str(depth2)}:{str(dp)}")
    #heterozygote: low-DP dropout → miscall as homo (REF-tilted)
    if dp <= dropout_depth:
        homo = '0' if random.random() < ref_bias else '1'
        depth1, depth2 = (dp, 0) if homo == '0' else (0, dp)
        return f"{homo}{sep}{homo}:{depth1},{depth2}:{dp}"

    # if there are more than 3 reads, distribute the total reads between alleles
    depth1 = np.random.binomial(dp, ref_bias)
    depth2 = dp - depth1
   
    # if all reads go to one allele, also change GT to homozygous
    if depth1 == 0:
        # all ALT reads → 1/1
        return f"1{sep}1:0,{dp}:{dp}"
    elif depth2 == 0:
        # all REF reads → 0/0
        return f"0{sep}0:{dp},0:{dp}"

    # otherwise it's a true heterozygote
    return(f"{pos[:3]}:{str(depth1)},{str(depth2)}:{str(dp)}")

def add_depth(vcf_path, new_vcf, sample_list, mean, var, ref_bias, dropout, dist):
    """
    Function that simulates depth in the specified individuals, from the specified distribution 
    with specified mean and variance. Then distributes the reads across alleles with the given 
    dropout threshold and reference bias. Outputs :: or depth annotations to all positions, so that 
    downstream analysis with bcftools etc. can be conducted. Removes multiallelic positions.
    """
    exclude = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

    with open(vcf_path, "r") as file1, open(new_vcf, "w") as outfile:
        for line in file1:
            line = line.strip().split("\t")
            if line[0].startswith("#"):
                if "#CHROM" in line[0]:
                    # add on header with FORMAT information (so that the depth we add can still be parsed by bcftools etc.)
                    outfile.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n')
                    outfile.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
                    header_ix, include, name = parse_header(line, sample_list)
                    if include == []:
                        raise Exception("Provided population names not in file. Breaking")
                    else:
                        print(f"Simulation depth in: {name}")

                    exclude_ix = [i for i,name in enumerate(line) if name in exclude]
                    outfile.write("\t".join(line) + "\n")
                # write out all header lines
                else:
                    outfile.write("\t".join(line) + "\n")

            else:
                # this allows us to add AD and DP to all lines:
                metrics = [line[i] for i in range(len(line)) if i in exclude_ix]
                metrics[header_ix["format_ix"]] += ":AD:DP"

                # add in depth for every line - but only to the specified individuals
                # all other individuals get a "::" added to their FORMAT field
                if not multiallelic(line, header_ix):
                    line = [pos_depth(line[i], mean, var, ref_bias, dropout, dist) if i in include else (line[i] + "::") for i in range(len(line))]
                    # remove "metrics" from the line (those don't need any depth annotations)
                    line = [line[i] for i in range(len(line)) if i not in exclude_ix]
                    # then recombine
                    line = metrics + line

                    outfile.write("\t".join(line) + "\n")
