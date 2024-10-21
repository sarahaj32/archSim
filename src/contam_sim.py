"""
Helper file with functions that adds ancestral contamination 
or modern human contamination from a given set of modern populations
at a given contamination rate
"""

from helper_functions import parse_header
import random

def anc_geno_call(base, rate):
    """
    Helper function that adds ancestral contamination to a position
    base is a biallelic genotype 
    rate is the rate of contamination
    simulation retains phasing if included in input file
    """
    # change any derived position
    if base[0] == '0' and base[2] == '0':
        return base
    elif random.random() >= rate:
        return base
    # convert to the ancestral position at a set rate
    # maintain the phasing of the initial file
    else:
        if (base[0] == '0' and base[2] == '1') or (base[0] == '1' and base[2] == '0'):
            return(f"0{base[1]}0{base[3:]}")
        elif random.random() >= 0.5:
            if base[1] == "/":
                return("0/1")
            elif random.random() >= 0.5: # this should be phased data
                return(f"0{base[1]}1{base[3:]}")
            else:
                return(f"1{base[1]}0{base[3:]}")
        else:
            return(f"0{base[1]}0{base[3:]}")

def add_anc_contamination(vcf_path, new_vcf, sample_list, contamination):
    """
    converts derived positions to the reference at the given contamination rate
    heterozygous positions become homozygous reference
    homozgyous alternative positions become heterozygous
    """
    with open(vcf_path, "r") as file1, open(new_vcf, "w") as outfile:
        for line in file1:
            line = line.strip().split("\t")
            if line[0].startswith("#"): 
                outfile.write("\t".join(line) + "\n")
                if "#CHROM" in line[0]:
                    _, include = parse_header(line, sample_list)
                    if include == []:
                        print("provided population names not in file. Breaking")
                        break
                    else:
                        print(f"making {sample_list} samples contaminated")
                    print(f"include:{include}")
            else:
                # add contamination to the populations of interest 
                line = [anc_geno_call(line[i], contamination) if i in include else line[i] for i in range(len(line))]
                outfile.write("\t".join(line) + "\n")

def add_contam_object(dict, pos, cont_range, receive, contam):
    """
    helper function that adds modern human contaminating chunks 
    to a dictionary to keep track of. start position, lenght of 
    contaminating fragment, individual receiving, individual contaminating
    """
    
    dict[receive] = {"start": pos, "end": pos + cont_range, "contam_ix": contam}

def add_mh_contamination(vcf_path, new_vcf, sample_list, modern_human, rate, cont_range):
    """
    A function that will simulate modern human contamination
    by replacing the "receiving" individuals genotypes with the "contaminating" 
    individuals genotypes. Randomly selects a contaminating individual from the modern human
    populations. Converts genotypes in genomic chunks
    """
    contam_dict={}
    exclude = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    with open(vcf_path, "r") as file1, open(new_vcf, "w") as outfile:
        for line in file1:
            line = line.strip().split("\t")
            if line[0].startswith("#"): 
                outfile.write("\t".join(line) + "\n")
                if "#CHROM" in line[0]:
                    header_ix, include = parse_header(line, sample_list)
                    modern = [i for i,name in enumerate(line) if any([f"{s}_" in name for s in modern_human])]
                    if modern == []:
                        print("list of modern populations is required")
                        break
                    print(f"indices of modern pops to contaminate with: {modern}")
                    if sample_list == []:
                        # contaminate all non-modern samples (so just remove VCF specific rows)
                        to_contam = [i for i,name in enumerate(line) if name not in exclude and name not in modern_human]
                    else:
                        to_contam = [i for i,name in enumerate(line) if name in include]
                        if to_contam == []:
                            print("provided population names not in file. Breaking")
                            print(f"select from: {line}")
                            break
            else:
                pos = int(line[header_ix["pos_ix"]])
                # loop through samples to receive contamination:
                for i in to_contam:
                    # see if the position is in an existing modern human contamination chunk 
                    if i in contam_dict and contam_dict[i]["start"] <= pos <= contam_dict[i]["end"]:
                        # replace receiving genotype with contaminating genotype
                        line[i] = line[contam_dict[i]["contam_ix"]]
                    else:
                        # remove the modern human chunk if we have moved past it
                        if i in contam_dict: 
                            del contam_dict[i]
                        # at a given rate, create another modern human contamination chunk
                        # save the start and end of this chunk so we can contaminate any other positions
                        # within the chunk
                        if random.random() <= rate:
                            add_contam_object(contam_dict, pos, cont_range, i, random.choice(modern))
                            line[i] = line[contam_dict[i]["contam_ix"]]
                outfile.write("\t".join(line) + "\n")