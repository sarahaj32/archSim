from helper_functions import parse_header
import numpy as np
import random

def pos_depth(pos, mean, var, apply_filter):
    depth = max(round(np.random.normal(mean, var)), 0)
    depth1 = depth2 = depth
    # if not apply_filter:
    #     return(f"{pos}:{str(depth1)},{str(depth2)}:{str(depth)}")
    # else:
    # missing data
    if depth == 0:
        gt = "./."
    
    # homozygous
    elif pos[0] == pos[2]:
        gt = pos[0:3]
        if pos[0] == '0':
            depth2 = 0
        else:
            depth1 = 0
    # heterozygous: induce false homozygous
    elif depth <= 3:
        allele = random.randint(0, 1) # could adjust this line to induce a reference bias
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

def add_depth(vcf_path, new_vcf, sample_list, mean, var, distribution, apply_filter):
    exclude = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

    with open(vcf_path, "r") as file1, open(new_vcf, "w") as outfile:
#        counter = 1
        for line in file1:
            # counter += 1
            # if counter > 20:
            #     break
            line = line.strip().split("\t")
            if line[0].startswith("#"):  # looking at if every time
                if "#CHROM" in line[0]:
                    header_ix, include = parse_header(line, sample_list)
                    exclude_ix = [i for i,name in enumerate(line) if name in exclude]
                    outfile.write("\t".join(line) + "\n")
                else:
                    outfile.write("\t".join(line) + "\n")
            else:
                metrics = [line[i] for i in range(len(line)) if i in exclude_ix]
                if line[header_ix["alt_ix"]] != ".":
                    line = [pos_depth(line[i], mean, var, apply_filter) if i in include else (line[i] + "::") for i in range(len(line))]
                line = [line[i] for i in range(len(line)) if i not in exclude_ix]
                metrics[header_ix["format_ix"]] += ":AD:DP"
                line = metrics + line
                outfile.write("\t".join(line) + "\n")
