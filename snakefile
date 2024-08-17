# create visual of workflow:
# snakemake --forceall --rulegraph | dot -Tpdf > output/workflow_dag.pdf
# snakemake --dag /global/scratch/users/sarahj32/ancient_individuals_sims/output/{A,B}.ps > dag.svg


# set a variable to be used in rule files
lab_scripts = "/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/bin/"

# read in simulation parameters from the config file
reps = [i for i in range(config["reps"])] # 100
ds = config["ds"]
chromosomes = config["chromosomes"]
tests = config["tests"]
statistics = config["statistics"]
chroms = config["chromosomes"]

# set default rates if they're not specified
if "rates" in config.keys():
    if "deamination" in config["rates"]:
        deam_rates = config["rates"]["deamination"]
    else:
        deam_rates = [0.05]
    if "contamination" in config["rates"]:
        contam_rates = config["rates"]["contamination"]
    else:
        contam_rates = [0.05]
    if "missing" in config["rates"]:
        missing_rates = config["rates"]["missing"]
    else:
        missing_rates = [0.1]
else:
    contam_rates = deam_rates = missing_rates = [0.05]

parameters = config["parameters"]
if "deaminate" in parameters:
    parameters = [i for i in parameters if i != "deaminate" ]
    parameters.extend([f"deaminate-{r}" for r in deam_rates])

if "ancContam" in parameters:
    parameters = [i for i in parameters if i != "ancContam"]
    parameters.extend([f"ancContam-{r}" for r in contam_rates])

if "mhContam" in parameters:
    parameters = [i for i in parameters if i != "mhContam"]
    parameters.extend([f"mhContam-{r}" for r in contam_rates])
if "missing" in parameters:
    parameters = [i for i in parameters if i != "missing"]
    parameters.extend([f"missing-{r}" for r in missing_rates])

#configfile: "workflow/configs/archaic.yaml"
print("HERE")
print(parameters)

# if "cumulative" in config.keys() and config["cumulative"] == "yes":
#     include: "workflow/simulate_cumulatively.smk"
# else:
print("NEW")
#"
include: "workflow/combine_vcfs.smk"
include: "workflow/downsample.smk"
include: "workflow/full.smk"
include: "workflow/calculate.smk"
include: "workflow/simulate_individually.smk"

print("NEXT")
rule all:
    input:
        #expand("data/{timescale}/original/full/vcf/original_full_{chrom}.vcf", timescale = config['timescale'], chrom = chroms ),
        expand("results/{timescale}/full_summary_statistics/{timescale}_{param}_stats.tsv", timescale = config["timescale"], param = parameters),
        expand("data/{timescale}/{param}/{d}/eig/{param}_{d}-{n}.snp", timescale = config["timescale"], param = parameters, d = [i for i in ds if i != "full" ], n = reps)
        #expand("analysis/{timescale}/{timescale}_fstats_results.tsv", timescale = config['timescale'])



