from helper_functions import parse_header
import random

def anc_geno_call(base, rate):
    # change any derived position
    if base[0] == '0' and base[2] == '0':
        return base
    elif random.random() >= rate:
        return base
    # convert to the ancestral position at a set rate
    # maintain the phasing of the initial file
    else:
        if (base[0] == '0' and base[2] == '1') or (base[0] == '1' and base[2] == '0'):
            # print("haplo")
            return(f"0{base[1]}0{base[3:]}")
        elif random.random() >= 0.5:
            # print("homo")
            if base[1] == "/":
                return("0/1")
            elif random.random() >= 0.5: # this should be phased data
                return(f"0{base[1]}1{base[3:]}")
            else:
                return(f"1{base[1]}0{base[3:]}")
        else:
            # print("homo2")
            return(f"0{base[1]}0{base[3:]}")

def add_anc_contamination(vcf_path, new_vcf, sample_list, contamination):
    # seems to be 1 order of magnitude off of what the rate should be?!?!
    with open(vcf_path, "r") as file1, open(new_vcf, "w") as outfile:
        for line in file1:
            line = line.strip().split("\t")
            if line[0].startswith("#"):  # looking at if every time
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
                line = [anc_geno_call(line[i], contamination) if i in include else line[i] for i in range(len(line))]
                outfile.write("\t".join(line) + "\n")

def add_contam_object(dict, pos, cont_range, receive, contam):
    dict[receive] = {"start": pos, "end": pos + cont_range, "contam_ix": contam}

def add_mh_contamination(vcf_path, new_vcf, sample_list, modern_human, rate, cont_range):
    # modern_human and rate are required
    contam_dict={}
    exclude = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    print("HERE")
    print(modern_human)
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
                        to_contam = [i for i,name in enumerate(line) if any([f"{s}_" in name for s in sample_list])]
                        if to_contam == []:
                            print("provided population names not in file. Breaking")
                            print(f"select from: {line}")
                            break
            else:
                pos = int(line[header_ix["pos_ix"]])
                # get index of specified samples
                for i in to_contam:
                    if i in contam_dict and contam_dict[i]["start"] <= pos <= contam_dict[i]["end"]:
                        line[i] = line[contam_dict[i]["contam_ix"]]
                    else:
                        if i in contam_dict: 
                            del contam_dict[i]
                        if random.random() <= rate:
                            # updates the dictionary 
                            add_contam_object(contam_dict, pos, cont_range, i, random.choice(modern))
                            line[i] = line[contam_dict[i]["contam_ix"]]
                outfile.write("\t".join(line) + "\n")