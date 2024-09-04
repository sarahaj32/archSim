rule make_full_eigenstrat:
    input:
        "data/{timescale}/{param}/full/vcf/{param}_full.vcf"
    output:
        geno = "data/{timescale}/{param}/full/eig/{param}_full.geno", 
        snp = "data/{timescale}/{param}/full/eig/{param}_full.snp", 
        # output from the make eigenstrat script
        ind_full = "data/{timescale}/{param}/full/eig/{param}_full.ind",
        # newly created file with population info
        ind_all = "data/{timescale}/{param}/full/eig/{param}_all.ind"

    params:
        tmp_file = "data/{timescale}/{param}/full/eig/tmp_full"
    shell:
        """ 
        python2 {lab_scripts}/vcf2eigenstrat.py -v {input} -o data/{wildcards.timescale}/{wildcards.param}/full/eig/{wildcards.param}_full
        awk '$3=$1' {output.ind_full} | awk 'BEGIN{{FS=OFS="_"}}NF--' > {params.tmp_file}
        mv {params.tmp_file} {output.ind_all}
        """ 

rule calc_full_fstatistics:
    input:
        geno = "data/{timescale}/{param}/full/eig/{param}_full.geno",
        snp = "data/{timescale}/{param}/full/eig/{param}_full.snp", 
        full_ind = "data/{timescale}/original/full/eig/original_all.ind",
        select_ind=expand("data/{timescale}/original/full/eig/original_{s}.ind", s = config["selection"], timescale = config["timescale"]),
        script = config["calc_script"]       
    output:
        "data_output/{timescale}/{param}/full/{stat}/{param}_full_{stat}_{s}.tsv"
    params:
        ind_basepath = "data/{timescale}/original/full/eig/original_",
        time='03:00:00'
    resources: 
        cpus=12
    shell:
        """
        Rscript {input.script} {input.geno} {input.snp} {params.ind_basepath} {wildcards.s} data_output/{wildcards.timescale}/{wildcards.param}/full/ {wildcards.stat}
        """

rule calc_fst:
# calculate FST and run pca on the full simulation data
    input:
        params = "helper_files/pca_parameters",
        geno = "data/{timescale}/original/full/eig/original_full.geno",
        snp = "data/{timescale}/original/full/eig/original_full.snp",
        ind = "data/{timescale}/original/full/eig/original_all.ind"
    
    output:
        evec = "output/{timescale}/pca/pca.evec",
        eval_file = "output/{timescale}/pca/pca.eval",
        phyl = "output/{timescale}/pca/pca.phyl"
    shell:
        """
        module load gsl
        cut -d' ' -f3 {input.ind} | uniq > helper_files/pca_popfile
        /global/scratch/users/sarahj32/software/miniconda3/envs/r_python/bin/smartpca -p {input.params} > output/{wildcards.timescale}/pca/pca_output
        sed -i -e 's/\-0.00[0-9]/ 0/g' {output.phyl}
        """
