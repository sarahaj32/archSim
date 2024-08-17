import json

# read in the number of individuals - useful for removing singletons
json_file = config["demography"]["samples_json"]
with open(json_file, "r") as f:
    data = json.load(f)
indivs_list = [data[pop]["num_samples"] for pop in data.keys()]
indivs = sum(indivs_list)
indiv_chroms = indivs * 2
max_haplo = int(indiv_chroms) - 2
# when pseudohaploid - we need to remove 2 more chromosomes (since they are the same within an individual)
max_haplo_ph = max_haplo - 2

if "pops" in config.keys() and "target" in config["pops"].keys():
    target_pops = ",".join(config["pops"]["target"])
    print(target_pops)
else:
    target_pops = ""

if "mh_contam" in config["parameters"]:
    modern_pops = ",".join(config["pops"]["modern"])


rule simulate:
# takes ~ 10 minutes
    input:
        deme = config["demography"]["demo_yaml"],
        json = config["demography"]["samples_json"],
        script = "helper_files/simulate_from_deme.py"
    output:
        "data/{timescale}/original/full/vcf/original_filtered_{chrom}.vcf"
    conda:
        # figure out why i can no longer specify the yml file, that would really be ideal
        # "envs/sims_environment.yml"
        "sims"
    resources: 
        cpus = 2
    shell: 
     # simulate populations using msprime, from yaml file and json file of individuals
        """
        echo {input.script}
        python {input.script} -y {input.deme} -j {input.json} -c {wildcards.chrom} -rl -o {output}
        """ 

rule make_pseudohaploid: 
    # randomly select an allele for all haploid positions in select individuals
    input:
        vcf = "data/{timescale}/original/full/vcf/original_filtered_{chrom}.vcf",
        script = "src/main.py"
    output:
        "data/{timescale}/pseudohaploid/full/vcf/pseudohaploid_filtered_{chrom}.vcf"
    params:
        pops = target_pops
    shell:
        """
        python {input.script} pseudohaploid -v {input.vcf} -t {params.pops} -o {output}
        """

# add deamination at each specified amount
for deam_rate in deam_rates:
    rule:
        name: f"deaminate-{deam_rate}"
        input:
            vcf = "data/{timescale}/original/full/vcf/original_filtered_{chrom}.vcf",
            script = "src/main.py"
        output:
            f"data/{{timescale}}/deaminate-{deam_rate}/full/vcf/deaminate-{{deam_rate}}_filtered_{{chrom}}.vcf"
        params:
            pops = target_pops,
            rate = deam_rate
        shell:
            "python {input.script} deaminate -t {params.pops} -v {input.vcf} -o {output} -r {params.rate}"

# add contamination at each specified amount
for anc_contam_rate in contam_rates:
    rule:
        name: f"ancContam-{anc_contam_rate}"
        input:
            vcf = "data/{timescale}/original/full/vcf/original_filtered_{chrom}.vcf",
            script = "src/main.py"
        output:
            f"data/{{timescale}}/ancContam-{anc_contam_rate}/full/vcf/ancContam-{anc_contam_rate}_filtered_{{chrom}}.vcf"
        params:
            pops = target_pops,
            rate = anc_contam_rate
        shell:
            "python {input.script} contaminate -t {params.pops} -v {input.vcf} -o {output} -anc -r {params.rate}"

for mh_contam_rate in contam_rates:
    rule:
        name: f"mhContam-{mh_contam_rate}"
        input:
            vcf = "data/{timescale}/original/full/vcf/original_filtered_{chrom}.vcf",
            script = "src/main.py"
        output:
            f"data/{{timescale}}/mhContam-{mh_contam_rate}/full/vcf/mhContam-{mh_contam_rate}_filtered_{{chrom}}.vcf"
        params:
            pops = target_pops,
            modern_pops = config["pops"]["modern"],
            rate = mh_contam_rate
        shell:
            "python {input.script} contaminate -t {params.pops} -v {input.vcf} -o {output} -mh -r {params.rate} -modern {params.modern_pops}"

rule add_coverage:
    input:
        vcf = "data/{timescale}/original/full/vcf/original_filtered_{chrom}.vcf",
        script = "src/main.py"
    output:
        "data/{timescale}/dpFilter/full/vcf/dpFilter_filtered_{chrom}.vcf"
    params:
        pops = target_pops
    shell:
        "python {input.script} dpFilter -t {params.pops} -v {input.vcf} -o {output}"

rule add_coverage_and_missing:
    input:
        vcf = "data/{timescale}/original/full/vcf/original_filtered_{chrom}.vcf",
        script = "src/main.py"
    output:
        "data/{timescale}/dpFilterMissing/full/vcf/dpFilterMissing_filtered_{chrom}.vcf"
    params:
        pops = target_pops
    shell:
        "python {input.script} dpFilter -t {params.pops} -v {input.vcf} -o {output} -missing"

#  rule add_missingness:
#     input:
#         vcf = "data/{timescale}/original/full/vcf/original_filtered_{chrom}.vcf",
#         script = "src/main.py"
#     output:
#         "data/{timescale}/dpFilter/full/vcf/dpFilter_filtered_{chrom}.vcf"
#     params:
#         pops = target_pops
#     shell:
#         "python {input.script} dpFilter -t {params.pops} -v {input.vcf} -o {output} -r "
   
for missing_rate in missing_rates:
    rule:
        name: f"missing-{missing_rate}"
        input:
            vcf = "data/{timescale}/original/full/vcf/original_filtered_{chrom}.vcf",
            script = "src/main.py"
        output:
            f"data/{{timescale}}/missing-{missing_rate}/full/vcf/missing-{missing_rate}_filtered_{{chrom}}.vcf"
        params:
            pops = target_pops,
            rate = missing_rate
        shell:
            "python {input.script} missing -t {params.pops} -v {input.vcf} -o {output} -r {params.rate}"
