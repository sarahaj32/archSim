from helper_functions import parse_header
import random

def ph_geno_call(base):
    if base[0] == base[2]:
        return base
    else:
        call = round(random.random())
        return(str(call) + base[1] + str(call) + base[3:])

def make_pseudohaploid(vcf_path, new_vcf, sample_list):
    with open(vcf_path, "r") as file1, open(new_vcf, "w") as outfile:
        for line in file1:
            line = line.strip().split("\t")
            if line[0].startswith("#"):  # looking at if every time
                outfile.write("\t".join(line) + "\n")
                if "#CHROM" in line[0]:
                    # include: the indices of the samples for the simulation to be applied to
                    _, include = parse_header(line, sample_list)
                    # if include == []:
                    #     print("provided population names not in file. Breaking")
                    #     break
                    # else:
                        # update this to say all samples if no sample_list is provided
                    print(f"making {sample_list} samples pseudohaploid")
                    #print(f"include:{include}")
            else:
                line = [ph_geno_call(line[i]) if i in include else line[i] for i in range(len(line))]
                outfile.write("\t".join(line) + "\n")
