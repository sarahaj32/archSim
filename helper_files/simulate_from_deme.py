import msprime
import argparse
import demes 
import json

def get_SampleSet_from_json(json_file, gen_years):
    with open(json_file, "r") as f:
        data = json.load(f)
    samples = []
    for key in data.keys():
        print(key)
        num_samples = data[key]["num_samples"]
        population = data[key]["population"]
        if "years" in data[key].keys():
            gens = data[key]["years"] / gen_years
        elif "generations" in data[key].keys():
            gens = data[key]["generations"]
        else: 
            print("using default sampling time: 0")
            gens = 0
        if "ploidy" in data[key].keys():
            ploidy = data[key]["ploidy"]
        else:
            print("using default ploidy: 2")
            ploidy = 2
        samples.append(msprime.SampleSet(num_samples = num_samples, population = population, time = gens, ploidy = ploidy))
    return(samples)

def get_SampleNames_from_json(json_file):
    with open(json_file, "r") as f:
        data = json.load(f)
    indiv_names = []
    for key in data.keys():
        indiv_names.extend([key + '_' + str(j + 1) for j in range(data[key]["num_samples"])])
    return(indiv_names)
    
def write_vcf(mts, indiv_names, vcf_name, rep):
    # create names and write out vcf:
    with open(vcf_name, "w") as vcf_file:
            mts.write_vcf(vcf_file, individual_names = indiv_names, contig_id = rep)
            print("Done")

def simulate(graph, samples, length, mu, rr, rep):
    print("starting simulation")
    ts = msprime.sim_ancestry(
        samples=samples,
        demography=graph,
        sequence_length=length,
        recombination_rate=rr,
        random_seed=rep+1,
        record_migrations=True
    )
    print("done, adding mutations")
    # Overlay mutations.
    mts = msprime.sim_mutations(
        tree_sequence=ts, rate=mu,
        random_seed=rep+1
        #discrete_genome=False
    )
    print("done, writing out")
    return(mts)

parser = argparse.ArgumentParser("Process input")
parser.add_argument("-c", "--chrom", help = "chromosome or repetition number", type = int, default = 22)
parser.add_argument("-y", "--yaml", help = "path to yaml file with demography", type = str, required = True)
parser.add_argument("-o", "--out_file", help = "path and name output vcf file", type = str, default = "simulation.vcf")
parser.add_argument("-l", "--length", help = "length of chromsome", type = int, default = 10**7)
parser.add_argument("-rl", "--real_lengths", help = "simulate the true lengths of chromosomes", action = 'store_true')
parser.add_argument("-m", "--mu", help = "mutation rate", type = float, default = 1.5e-8)
parser.add_argument("-r", "--rr", help = "recombination rate", type = float, default = 1e-8)
parser.add_argument("-g", "--gens", help = "generation time", type = int, default = 29)
parser.add_argument("-j", "--json", help = "json file with individuals info", type = str, required = True)
parser.add_argument("-t", "--tree_file", help = "path and name of output tree sequence", type = str)

args = parser.parse_args()

out_file = args.out_file
gens = args.gens
json_file = args.json
tree_file = args.tree_file

if args.real_lengths:
    real_lengths_dict = {1:2.49e8, 2:2.42e8, 3:1.98e8, 4:1.90e8, 5:1.82e8, 6:1.71e8, 7:1.59e8, 8:1.45e8, 9:1.38e8, 10:1.34e8, 11:1.35e8, 12:1.33e8, 13:1.14e8,14:1.07e8, 15:1.02e8, 16:9.03e7, 17:8.33e7, 18:8.04e7, 19:5.86e7, 20:6.44e7, 21:4.67e7, 22:5.08e07}
    length = real_lengths_dict[args.chrom]
    print(f"simulating chromosome {args.chrom} of length {length}")
else:
    length = args.length


if __name__ == "__main__":
    print("Beginning")

    print(f"mu= {args.mu}")
    print(f"rr= {args.rr}")
    # set up the demography, currently with defaults 
    demo_deme = demes.load(args.yaml)

    demo_graph = msprime.Demography.from_demes(demo_deme)

    # read in sampling parameters from json file
    samples = get_SampleSet_from_json(json_file, gens)
    indiv_names = get_SampleNames_from_json(json_file)

    # set parameters
    print("simulating")
    # simulate data
    mts = simulate(graph = demo_graph, samples = samples, length = length, mu = args.mu, rr = args.rr, rep = args.chrom)

    # write to vcf
    write_vcf(mts = mts, indiv_names = indiv_names, vcf_name = out_file, rep = args.chrom)
    
    # also save the ts so that we can use it later to look at coalescence times and snp ages, etc.
    if tree_file != None:
        print(f"writing out tree file: {tree_file}")
        mts.dump(tree_file)


