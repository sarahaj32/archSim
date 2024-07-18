def parse_header(line, sample_list):
    ref_ix = line.index("REF")
    alt_ix = line.index("ALT")
    pos_ix = line.index("POS")
    format_ix = line.index("FORMAT")

    # write a dictionary that stores all of the header index information
    header_ix = {"ref_ix":ref_ix, "alt_ix":alt_ix, "pos_ix":pos_ix, "format_ix":format_ix} 
    
    if sample_list == []:
        # use all samples (so just remove VCF specific rows)
        exclude = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        include = [i for i,name in enumerate(line) if name not in exclude]
        sample_list = [name for i,name in enumerate(line) if name not in exclude]
    else:
        # get index of specified samples
        include = [i for i,name in enumerate(line) if any([f"{s}_" in name for s in sample_list])]
    if include == []:
        # print("provided population names not in file. Breaking")
        raise Exception("provided population names not in file. Breaking")
    return header_ix, include

def parse_target_indivs(input):
    sample_list = input.strip()
    if sample_list == "":
        sample_list = []
    else:
        sample_list = [s.strip() for s in sample_list.split(",")]  
    return(sample_list)

