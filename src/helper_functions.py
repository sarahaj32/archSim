import json

def parse_indivs(samp_input, ind_type):
    """
    A helper function that returns a list of individuals parsed from either a json file or a comma-separated string
    """
    if ".json" in samp_input:
        return(parse_indivs_from_json(samp_input, ind_type))
    # if the input is a comma-separated string, parse the individuals
    sample_list = samp_input.strip()
    if sample_list == "":
        sample_list = []
    else:
        sample_list = [s.strip() for s in sample_list.split(",")]  
    return(sample_list)

def parse_indivs_from_json(samp_json, ind_type):
    """
    A helper function that parses a json file, and returns a list of individuals matching the ind_type
    ind_type should be either "target" or "contam" to match the json file
    """
    with open(samp_json, "r") as f:
        data = json.load(f)
    # make all keys lowercase in case they were typed differently in the json file
    data = {key.lower():value for key,value in data.items()}
    # select the target or contam individuals
    samples = data[ind_type]
    sample_list = [s.strip() for s in samples]  
    return(sample_list)

def parse_header(line, sample_list):
    """
    A helper function that parses the header line of a VCF file
    and outputs a dictionary with the index of the:
    -reference allele, -alternative allele, -position, -format
    Additionally creates and outputs a list of indices for the specified samples
    """

    ref_ix = line.index("REF")
    alt_ix = line.index("ALT")
    pos_ix = line.index("POS")
    format_ix = line.index("FORMAT")

    # write a dictionary that stores all of the header index information
    header_ix = {"ref_ix":ref_ix, "alt_ix":alt_ix, "pos_ix":pos_ix, "format_ix":format_ix} 
    
    if sample_list == []:
        # if there are no specified samples, use all samples (so just remove VCF specific rows)
        exclude = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        include_idx = [i for i,name in enumerate(line) if name not in exclude]
    else:
        # get index of specified samples
        # this is matching if the only population name was given
        include_idx = [i for i,name in enumerate(line) if any([f"{s}_" in name for s in sample_list])]
        # if that doesn't give any matches, then look for exact matches
        if include_idx == []:
            include_idx = [i for i,name in enumerate(line) if any([f"{s}" == name for s in sample_list])]

        missing = [s for s in sample_list if not any([s in name for name in line])]
        if missing != []:
            print(f"\nWarning: the following specified samples were not found in the VCF and will be ignored: {missing}\n")
    if include_idx == []:
        raise Exception("None of the provided sample names were found in file. Breaking")
    include_names = [line[i] for i in include_idx]
    return header_ix, include_idx, include_names

def multiallelic(line, header_idx):
    if "," in line[header_idx["alt_ix"]]:
        print(f"Skipping multiallelic site at position {line[header_idx['pos_ix']]}")
        return True
    else:
        return False