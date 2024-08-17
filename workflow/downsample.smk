rule make_downsampled:
    input:
        vcf = "data/{timescale}/{param}/full/vcf/{param}_full.vcf",
        script = "src/main.py"
    output:
        "data/{timescale}/{param}/{d}/vcf/{param}_{d}-{n}.vcf"
    run: 
        # somehow this isn't working - it's not skipping the full
        if wildcards.d != "full":
            shell("python {input.script} downsample -v {input.vcf} -n {wildcards.d} -o {output}")

rule make_eigenstrat:
# we only need the geno and snp files, will use the same .ind file for all files (since the populations are unchanged)
    input:
        "data/{timescale}/{param}/{d}/vcf/{param}_{d}-{n}.vcf"
    output:
        geno = "data/{timescale}/{param}/{d}/eig/{param}_{d}-{n}.geno", 
        snp = "data/{timescale}/{param}/{d}/eig/{param}_{d}-{n}.snp"#, 
        #ind = "data/{timescale}/{param}/{d}/eig/{param}_{d}-{n}.ind"
    shell:
        """ 
        python2 {lab_scripts}/vcf2eigenstrat.py -v {input} -o data/{wildcards.timescale}/{wildcards.param}/{wildcards.d}/eig/{wildcards.param}_{wildcards.d}-{wildcards.n}
        """ 

rule make_selection_ind:
    # create .ind files that in some cases only have a select number (1) of individuals per population
    input:
        ind = expand("data/{timescale}/original/full/eig/original_all.ind", timescale = config["timescale"]), 
        script = "helper_files/make_selection_ind.sh"
    output:
        expand("data/{timescale}/original/full/eig/original_{select}.ind", select = [i for i in config["selection"] if i != "all"], timescale = config["timescale"])
    run:
        for s in config["selection"]: # ["all", "T", "AEO"]
            new_file = f"data/{config['timescale']}/original/full/eig/original_{s}.ind"
            if s != "all":
                # populations to downsize have n individuals
                pops = f"\"{' '.join(config['selection_dict'][s])}\""
                n = "1"
                # script to create new .ind file with fewer samples per selected population
                shell("helper_files/make_selection_ind.sh {input.ind} {new_file} {pops} {n}")
