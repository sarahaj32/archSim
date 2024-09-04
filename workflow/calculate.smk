rule calc_fstatistics:
    input:
        geno = "data/{timescale}/{param}/{d}/eig/{param}_{d}-{n}.geno",
        snp = "data/{timescale}/{param}/{d}/eig/{param}_{d}-{n}.snp", 
        # needs all of the ind files
        ind = "data/{timescale}/original/full/eig/original_all.ind",
        select_ind=expand("data/{timescale}/original/full/eig/original_{select}.ind", select = config["selection"], timescale = config["timescale"]),
        script = config["calc_script"]       
    output:
        "data_output/{timescale}/{param}/{d}/{stat}/{param}_{d}-{n}_{stat}_{s}.tsv" # s = [all, T, AEO]
    conda:
        "r_python"
    params:
        ind_basepath = "data/{timescale}/original/full/eig/original_"
    shell:
        """
        Rscript {input.script} {input.geno} {input.snp} {params.ind_basepath} {wildcards.s} data_output/{wildcards.timescale}/{wildcards.param}/{wildcards.d}/ {wildcards.stat}
        """   

rule calc_summary_statistics:
    input:
        "data/{timescale}/{param}/full/vcf/{param}_full.vcf"
    output:
        "results/{timescale}/full_summary_statistics/{timescale}_{param}_stats.tsv"
    shell:
        """
        bcftools stats -s - {input} > {output}
        """

rule summarize_results:
# call an R function that reads in and combines all of the output files with f-statistics 
    input:
        expand("data_output/{timescale}/{param}/full/{stat}/{param}_full_{stat}_{s}.tsv", param = parameters, test = tests, stat = statistics, s = config["selection"], timescale = config["timescale"]),
        expand("data_output/{timescale}/{param}/{d}/{stat}/{param}_{d}-{n}_{stat}_{s}.tsv", param = parameters, d = [i for i in ds if i != "full"], n = reps, stat = statistics, test = tests, s = config["selection"], timescale = config["timescale"])
    output:
        expand("analysis/{timescale}/{timescale}_fstats_results.tsv", timescale = config['timescale'])
    params:
        ts = config["timescale"],
        selection = config["selection"],
    conda:
        "r_python"
    shell:
        """
        echo "Rscript helper_files/fstat_results.R data_output/ analysis/{params.ts} '{parameters}' '{statistics}' '{ds}' '{reps}' '{params.selection}' '{params.ts}'"
        Rscript helper_files/fstat_results.R data_output/ analysis/{params.ts} '{parameters}' '{statistics}' '{ds}' '{reps}' '{params.selection}' '{params.ts}'
        """
# Rscript helper_files/fstat_results.R data_output/ analysis/human_sediment 'admix relna' 'original ph-deam-contam pseudohaploid deaminated coverage contam05T' 'dstat f3' 'full 1000000 100000 50000 30000 10000 5000 2500 1500 500 250' '0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99' 'T all' 'human_sediment'