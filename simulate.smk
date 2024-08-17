rule simulate: # takes ~10 minutes for 100Mb
# simulate specified populations using specified demography
# produce a vcf of the populations of interest, as well as the tree file from the simulation (so we can access and compute on it later)
    input: 
        deme = config["demography"]["demo_yaml"],
        json = config["demography"]["samples_json"],
        script = "../superarchaic_sims/workflow/simulate_from_deme_only.py"
    output: 
        vcf = "dat/raw/super_human_sim_{rep}.vcf",
        #bed = "dat/{target}/frags/introgressed_fragments_{rep}.bed",
        tree = "dat/raw/super_human_sim_{rep}.tree",
        biallelic = "dat/raw/super_human_sim_{rep}_bi.vcf"
        #super_bed = "dat/{target}/frags/introgressed_super_fragments_{rep}.bed"
    conda: 
        "sims"
    resources:
        # increasing the cpus is the way to increase memory available when submitting with savio_htc
        # these need a lot of memory
        cpus=3,
        time='3:00:00'
    #params:
       # source = lambda wildcards: config["targets"][wildcards.target],
       # time = config["demography"]["split_time"]
    shell: 
        """
        python {input.script} -j {input.json} -c {wildcards.rep} -y {input.deme} -o {output.vcf} -g 29 -l 100000000 -m 1.29e-8 -r 1e-8 -t {output.tree}
        bcftools view -m2 -M2 -v snps -o {output.biallelic} {output.vcf}
        """
#        python {input.script} -j {input.json} -c {wildcards.rep} -y {input.deme} -o {output.vcf} -b {output.bed} -t {wildcards.target} -s {params.source} -e {params.time} -g 29 -l 100000000 -m 1.29e-8 -r 1e-8 -f {output.tree}
#        sed -i 's/chr//g' {output.bed}

rule compress:
    input: 
        "dat/raw/super_human_sim_{rep}_bi.vcf"
    output: 
        bcf = "dat/raw/super_human_sim_{rep}.bcf",
        index = "dat/raw/super_human_sim_{rep}.bcf.csi"
    shell:
        """
        bcftools view -Ob {input} -o {output.bcf}
        bcftools index {output.bcf}
        """

rule calculate_window_stats:
    input:
        intro_vcf="dat/{target}/intro/introgressed_fragments_{rep}.vcf",
        not_intro_vcf="dat/{target}/not_intro/not_introgressed_fragments_{rep}.vcf"

        #den_snps = "dat/ArtificialGenome/fragment_bcfs/{name}/{name}_den_SNPs_chr{chr}.vcf"
    output:
        intro = "dat/{target}/results/intro/window/results_intro_den_{rep}.csv",
        not_intro = "dat/{target}/results/not_intro/window/results_not_intro_den_{rep}.csv",
        intro_ancden = "dat/{target}/results/intro_ancden/window/results_intro_ancden_{rep}.csv",
        not_intro_ancden = "dat/{target}/results/not_intro_ancden/window/results_not_intro_ancden_{rep}.csv"
    shell:
        """
        python ../superarchaic_analysis/helper_scripts/calc_statistics_precalc.py -r {wildcards.rep} -c window -w 10000000 -o {output.intro} -n Altai_1 -a AFR -m 50 -f {input.intro_vcf} -b 0 -d {wildcards.target}_1
        python ../superarchaic_analysis/helper_scripts/calc_statistics_precalc.py -r {wildcards.rep} -c window -w 10000000 -o {output.intro_ancden} -n Altai_1 -a AFR -m 50 -f {input.intro_vcf} -b 0 -d anc_den_1

        python ../superarchaic_analysis/helper_scripts/calc_statistics_precalc.py -r {wildcards.rep} -c window -w 10000000 -o {output.not_intro} -n Altai_1 -a AFR -m 50 -f {input.not_intro_vcf} -b 0 -d {wildcards.target}_1
        python ../superarchaic_analysis/helper_scripts/calc_statistics_precalc.py -r {wildcards.rep} -c window -w 10000000 -o {output.not_intro_ancden} -n Altai_1 -a AFR -m 50 -f {input.not_intro_vcf} -b 0 -d anc_den_1
        """

rule calculate_full_window_stats:
    input:
        "dat/{target}/raw/super_human_sim_{rep}_bi.vcf",
    output:
        "dat/{target}/results/full/window/results_full_den_{rep}.csv",
    shell:
        """
        python ../superarchaic_analysis/helper_scripts/calc_statistics_precalc.py -r {wildcards.rep} -c window -w 10000000 -o {output} -n Altai_1 -a AFR -m 50 -f {input} -b 0 -d {wildcards.target}_1
        """

