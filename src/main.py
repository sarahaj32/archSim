import argparse
from deaminate_sim import add_deam
from pseudohaploid_sim import make_pseudohaploid
from contam_sim import add_anc_contamination, add_mh_contamination
from downsample import downsample, add_missingness
from helper_functions import parse_target_indivs
from dp_filter_sim import add_depth

def main():
    print("main beginning")
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    subparser = parser.add_subparsers(dest = 'mode')

    # Pseudohaploid
    ph_subparser = subparser.add_parser('pseudohaploid', help='convert diploid data to pseudohaploid')
    ph_subparser.add_argument("-vcf",help="path to vcf to simulate damage in", type=str, required = True)
    ph_subparser.add_argument("-out",help="path to output simulated vcf", type=str, default = "out.vcf")
    ph_subparser.add_argument("-targets", metavar='',help="Target individuals to simulate features on", type=str, default = "")

    # Deamination
    deamination_subparser = subparser.add_parser('deaminate', help='Simulate deamination at given rate')
    deamination_subparser.add_argument("-vcf",help="path to vcf to simulate damage in", type=str, required = True)
    deamination_subparser.add_argument("-out",help="[required] path to output simulated vcf", type=str, default = "out.vcf")
    deamination_subparser.add_argument("-targets", metavar='',help="Target individuals to simulate features on", type=str, default = "")
    deamination_subparser.add_argument("-rate", metavar='',help="rate at which to convert to deminated", type=float, default=0.05)

    # Contamination
    contam_subparser = subparser.add_parser('contaminate', help='Simulate modern human contamination')
    contam_subparser.add_argument("-vcf",help="path to vcf to simulate damage in", type=str, required = True)
    contam_subparser.add_argument("-mh",help="modern human flag", action = 'store_true')
    contam_subparser.add_argument("-anc",help="ancestral flag human flag", action = 'store_true')
    contam_subparser.add_argument("-out",help="[required] path to output simulated vcf", type=str, default = "out.vcf") 
    contam_subparser.add_argument("-targets", metavar='',help="Target individuals to simulate features on", type=str, default = "")
    contam_subparser.add_argument("-rate", metavar='',help="rate at which to convert to deminated", type=float, default=0.05)
    contam_subparser.add_argument("-length", metavar='',help="length of contaminating modern human fragments", type=int, default=1000)
    contam_subparser.add_argument("-modern", metavar='',help="Modern human individuals to simulate features on", type=str, default = "")

    # downsample VCF
    downsample_subparser = subparser.add_parser('downsample', help='Downsample VCF')
    downsample_subparser.add_argument("-vcf",help="path to vcf to simulate damage in", type=str, required = True)
    downsample_subparser.add_argument("-out",help="[required] path to output simulated vcf", type=str, default = "out.vcf")
    downsample_subparser.add_argument("-num", metavar='',help="number of snps to downsample to", type=int, default=30_000)

    # add missingness
    missing_subparser = subparser.add_parser('missing', help='add missingness to VCF')
    missing_subparser.add_argument("-vcf",help="path to vcf to simulate damage in", type=str, required = True)
    missing_subparser.add_argument("-targets", metavar='',help="Target individuals to simulate features on", type=str, default = "")
    missing_subparser.add_argument("-out",help="[required] path to output simulated vcf", type=str, default = "out.vcf")
    missing_subparser.add_argument("-rate", metavar='',help="rate of missingness to induce", type=int, default=.1)

    # add depth and filter
    depth_subparser = subparser.add_parser('dp_filter', help='Downsample VCF')
    depth_subparser.add_argument("-vcf",help="path to vcf to simulate damage in", type=str, required = True)
    depth_subparser.add_argument("-out",help="[required] path to output simulated vcf", type=str, default = "out.vcf")
    depth_subparser.add_argument("-mean", metavar='',help="mean depth to simulate", type=int, default=5)
    depth_subparser.add_argument("-variance", metavar='',help="variance of depth to simulate", type=int, default=2)
    depth_subparser.add_argument("-distribution", metavar='',help="distribution to use to simulate depth from", type=str, default="normal")
    #depth_subparser.add_argument("-filter", help="flag that indiciates if filtering should be done", action = "store_true")
    depth_subparser.add_argument("-targets", metavar='',help="Target individuals to simulate features on", type=str, default = "")


    # could make target and sources files? 
    args = parser.parse_args()

    
    # run the features (not only one can be called at a time)
    if args.mode == "pseudohaploid":
        sample_list = parse_target_indivs(args.targets)
        print("Pseudohaplotyping")
        make_pseudohaploid(args.vcf, args.out, sample_list)
        print("finished pseudohaploid calls")
    
    if args.mode == "contaminate":
        sample_list = parse_target_indivs(args.targets)
        if args.anc:
            add_anc_contamination(args.vcf, args.out, sample_list, args.rate)
        elif args.mh:
            if args.modern == "":
                raise Exception("modern human population is necessary for adding modern human contamination")
            add_mh_contamination(args.vcf, args.out, sample_list, args.modern, args.rate, args.length)
        else:
            raise Exception("must specify either anc for ancestral contamination or mh for modern human contamination")

    if args.mode == "deaminate":
        sample_list = parse_target_indivs(args.targets)
        if not (0 < args.rate < 1):
            raise Exception("deamination rate must be above 0 and less than 1")
        else:
            add_deam(args.vcf, args.out, sample_list, args.rate)

    if args.mode == "dp_filter":
        sample_list = parse_target_indivs(args.targets)
        if args.distribution == "poisson":
            args.variance = args.mean
        else:
            add_depth(args.vcf, args.out, sample_list, args.mean, args.variance, args.distribution, args.filter)

    if args.mode == "downsample":
        if args.num <=0:
            raise Exception("num of positions must be > 0")
        else:
            downsample(args.vcf, args.out, args.num)

    if args.mode == "missing":
        sample_list = parse_target_indivs(args.targets)
        if not (0 < args.rate < 1):
            raise Exception("rate must be between 0 and 1")
        else:
            add_missingness(args.vcf, args.out, sample_list, args.rate)

    # print(args)
    # sample_list = args.target_populations.strip()
    # out_dir = args.out_dir
    # vcf_path = args.vcf_path
    # name = args.name
    # pseudohaploid = args.pseudo_haplotype
    # deaminate = args.deaminate
    # ds = args.downsample

    # modern_contamination = args.mh_contamination
    # modern_pops = args.modern_human.strip()
    # print(modern_pops)
    # mh_contam_range = args.mh_range

    # assert(deaminate >= 0 and deaminate < 1)
    # contamination = args.contamination
    # assert(contamination >= 0 and contamination < 1)
    # depth = args.reduce_depth
    # # ascertain = args.ascertain

    # # no longer using this option:
    # # rename files using the parameters in the output directory name
    # basename = vcf_path.split("/")[-1].replace(".vcf", "")
    # print(f"basename:{basename}")
    # #if name != "":
    # #name = f"_{name}"
    # #else:
    # #    name = "_simulated"

    # new_vcf = f"{out_dir}/{basename}{name}.vcf"
    # if "data" in out_dir:
    #     replace = str.split(out_dir, sep = "/")
    #     replace_name = replace[replace.index("data") + 2]
    #     new_vcf = re.sub("sim1_\w*_full", f"sim1_{replace_name}_full", new_vcf)

    # print(f"writing new file:{new_vcf}")
    # if sample_list == "":
    #     sample_list = []
    # else:
    #     sample_list = [s.strip() for s in sample_list.split(",")]

    # if modern_pops == "":
    #     modern_pops = []
    # else:
    #     modern_pops = [s.strip() for s in modern_pops.split(",")]

    # print(f"sample list: {sample_list}")


# print('reading in file')

if __name__ == "__main__":
    print("starting")
    main()

#     if deaminate > 0:
#         print("Deaminating")
#         add_deam(vcf_path, new_vcf, sample_list, deaminate)
#         print("finished adding deamination")
#     if contamination > 0:
#         print("Adding Ancestral Contamination")
#         add_contamination(vcf_path, new_vcf, sample_list, contamination)
#         print("finished adding ancestral contamination")
#     if depth > 0:
#         print("adding depths")
#         add_depth(vcf_path, new_vcf, sample_list, depth, 2)
#         print("finished simulating depths")
#     if ds > 0:
#         print("downsampling")
#         downsample(vcf_path, new_vcf, ds)
#     if modern_contamination > 0:
#         print('adding modern human contamination')
#         add_modern_human_cont(vcf_path, new_vcf, sample_list, modern_pops, modern_contamination, mh_contam_range)
    # if ascertain != "":
    #     ascertain_positions(vcf_path, new_vcf, ascertain)
    # # if i > 0: # i found a better way to do this - just change the ind files
    #     print("selecting samples")
    #     select_samples(vcf_path, new_vcf, sample_list, select, i)
    #     print("finished selecting samples")


# variables for testing:
# vcf_path = "/global/scratch/users/sarahj32/sim-fstats/data/medium/original/full/vcf/sim1_original_full_filtered_17.vcf"
# sample_list = ["pop_1", "pop_2", "pop_3"]
# new_vcf = "test.vcf"
# sample_list = []

# other example vcfs:
# bcftools view /global/scratch/p2p3/pl1_moorjani/sarahj32/superarchaic_introgression/superarchaic_analysis/data/VCFs/hg38/AltaiDen_merged_hg38_chr21.v 10000 > test_multi_geno.vcf 
# vcf_path = "test_multi_geno.vcf"
# sample_list = ["AltaiNea"]

# python helper_files/simulator_rf.py -v <vcf_path> -s <sample_list>
