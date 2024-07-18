import subprocess
import random
from helper_functions import parse_header

def get_keep_lines(vcf_path, num):
    # i should use bcftools - that would allow this to be gzipped
    all_lines = subprocess.run(["wc", "-l", vcf_path], capture_output = True, text = True).stdout
    # also find out how many header rows there are
    header_lines = subprocess.run(["grep", "^#", "-c", vcf_path], capture_output = True, text = True).stdout
    header_lines = int(header_lines.strip().split(" ")[0])

    # this is how many data lines there are in the file:
    all_lines = int(all_lines.strip().split(" ")[0]) - header_lines

    # select num rows from the all_lines available data rows
    to_keep = random.sample(list(range(all_lines)), num)
    to_keep.sort()
    return(to_keep)

def downsample(vcf_path, new_vcf, num):
    line_count = 0
    to_keep = get_keep_lines(vcf_path, num)
    with open(vcf_path, "r") as file1, open(new_vcf, "w") as outfile:
        for line in file1:
            line = line.strip().split("\t")
            if line[0].startswith("#"): 
                outfile.write("\t".join(line) + "\n")
            else:
                if line_count in to_keep:
                    outfile.write("\t".join(line) + "\n")
                line_count += 1

def missing(gt, rate):
    if random.random() <= rate:
        return("./.")
    else:
        return(gt)

def add_missingness(vcf_path, new_vcf, sample_list, rate):
    with open(vcf_path, "r") as file1, open(new_vcf, "w") as outfile:
        for line in file1:
            line = line.strip().split("\t")
            if line[0].startswith("#"): 
                outfile.write("\t".join(line) + "\n")
                if "#CHROM" in line[0]:
                    header_ix, include = parse_header(line, sample_list)
                    print(f"adding missingness to: {sample_list}")
                    print(f"include:{include}")
            else:
                line = [missing(line[i], rate) if i in include else line[i] for i in range(len(line))]
               # the line may be updated by the above step or not. Regardless, write it out
                outfile.write("\t".join(line) + "\n")

