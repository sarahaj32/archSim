import json
json_file = config["demography"]["samples_json"]
with open(json_file, "r") as f:
    data = json.load(f)
indivs_list = [data[pop]["num_samples"] for pop in data.keys()]
indivs = sum(indivs_list)
indiv_chroms = indivs * 2
max_haplo = int(indiv_chroms) - 2
# when pseudohaploid - we need to 
max_haplo_ph = max_haplo - 2

target_pops = ",".join(config["target_pops"])
print(target_pops)

rule simulate:
# takes ~ 10 minutes
    input:
        deme = config["demography"]["demo_yaml"],
        json = config["demography"]["samples_json"],
        script = "helper_files/simulate_from_deme.py"
    output:
        "data/{timescale}/original/full/vcf/sim1_original_full_{chrom}.vcf"
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

rule filter:
    # remove singletons
    # currently also filters to biallelic (that might not actually be necessary though)
    input:
        "data/{timescale}/original/full/vcf/sim1_original_full_{chrom}.vcf"
    output:
        "data/{timescale}/original/full/vcf/sim1_original_filtered_{chrom}.vcf"
    shell: 
        """
        module load bcftools
        wc -l {input}
        bcftools view -m2 -M2 --min-ac=2 --max-ac={max_haplo} -v snps {input} > {output}
        wc -l {output}
        """ 

# Define a rule for each sample using a for loop
print(config["deamination_rates"])
for deam_rate in config["deamination_rates"]:
    rule:
        name: f"deaminate-{deam_rate}"
        input:
            vcf = "data/{timescale}/original/full/vcf/sim1_original_filtered_{chrom}.vcf",
            script = "src/main.py"
        output:
            f"data/{{timescale}}/deaminated-{deam_rate}/full/vcf/sim1_deaminated-{{deam_rate}}_filtered_{{chrom}}.vcf"
        params:
            pops = target_pops,
            rate = deam_rate
        shell:
            "python {input.script} deaminate -t {params.pops} -v {input.vcf} -o {output} -r {params.rate}"

for anc_contam_rate in config["contamination_rates"]:
    rule:
        name: f"ancContam-{anc_contam_rate}"
        input:
            vcf = "data/{timescale}/original/full/vcf/sim1_original_filtered_{chrom}.vcf",
            script = "src/main.py"
        output:
            f"data/{{timescale}}/ancContam-{anc_contam_rate}/full/vcf/sim1_ancContam-{anc_contam_rate}_filtered_{{chrom}}.vcf"
        params:
            pops = target_pops,
            rate = anc_contam_rate
        shell:
            "python {input.script} contaminate -t {params.pops} -v {input.vcf} -o {output} -anc -r {params.rate}"

for mh_contam_rate in config["contamination_rates"]:
    rule:
        name: f"mhContam-{mh_contam_rate}"
        input:
            vcf = "data/{timescale}/original/full/vcf/sim1_original_filtered_{chrom}.vcf",
            script = "src/main.py"
        output:
            f"data/{{timescale}}/mhContam-{mh_contam_rate}/full/vcf/sim1_mhContam-{mh_contam_rate}_filtered_{{chrom}}.vcf"
        params:
            pops = target_pops,
            modern_pops = config["modern_pops"],
            rate = mh_contam_rate
        shell:
            "python {input.script} contaminate -t {params.pops} -v {input.vcf} -o {output} -mh -r {params.rate} -modern {params.modern_pops}"

rule make_pseudohaploid: 
    # randomly select an allele for all haploid positions in select individuals
    # again filter to remove singletons
    input:
        vcf = "data/{timescale}/original/full/vcf/sim1_original_filtered_{chrom}.vcf",
        script = "src/main.py"
    output:
        "data/{timescale}/pseudohaploid/full/vcf/sim1_pseudohaploid_filtered_{chrom}.vcf"
    params:
        tmp = "data/{timescale}/pseudohaploid/full/vcf/sim1_pseudohaploid_filtered_{chrom}_tmp.vcf",
        pops = target_pops
    shell:
        """
        module load bcftools
        python {input.script} pseudohaploid -v {input.vcf} -t {params.pops} -o {output}
        bcftools view -m2 -M2 --min-ac=4 --max-ac={max_haplo_ph} -v snps {output} > {params.tmp}
        mv {params.tmp}  {output}
        """
        
rule add_coverage:
# lets just do this to the original file intiially
    input:
        "data/{timescale}/original/full/vcf/sim1_original_full_filtered_{chrom}.vcf"
    output:
        "data/{timescale}/coverage/full/vcf/sim1_coverage_full_filtered_{chrom}.vcf"
    params:
        pops = target_pops,
        rate = mh_contam_rate
    shell:
            "python {input.script} contaminate -t {params.pops} -v {input.vcf} -o {output} -mh -r {params.rate} -modern {params.modern_pops}"