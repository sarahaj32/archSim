import argparse

from deaminate_sim import add_deam
from pseudohaploid_sim import make_pseudohaploid
from contam_sim import add_anc_contamination, add_mh_contamination
from downsample import downsample
from missing_sim import add_missingness
from helper_functions import parse_indivs
from dp_filter_sim import add_depth

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    subparser = parser.add_subparsers(dest = 'mode')

    # Pseudohaploid
    ph_subparser = subparser.add_parser('pseudohaploid', help='convert diploid data to pseudohaploid')
    ph_subparser.add_argument("-vcf",help="[Required] path to vcf to simulate damage in", type=str, required = True)
    ph_subparser.add_argument("-out",help="path to output simulated vcf", type=str, default = "pseudohaploid.vcf")
    ph_subparser.add_argument("-targets", metavar='',help="Target individuals to simulate features on, comma separated list or json file", type=str, default = "")

    # Deamination
    deamination_subparser = subparser.add_parser('deaminate', help='Simulate deamination at given rate')
    deamination_subparser.add_argument("-vcf",help="[Required] path to vcf to simulate damage in", type=str, required = True)
    deamination_subparser.add_argument("-out",help="path to output simulated vcf", type=str, default = "deaminate.vcf")
    deamination_subparser.add_argument("-targets", metavar='',help="Target individuals to simulate features on, comma separated list or json file", type=str, default = "")
    deamination_subparser.add_argument("-rate", metavar='',help="rate at which to convert to deminated", type=float, default=0.05)
    deamination_subparser.add_argument("-proportion", metavar='',help="proporion of sites that will be specified as transitions. If not specified, transitions will be determined by VCF alleles", type=float, default=None)

    # Contamination
    contam_subparser = subparser.add_parser('contaminate', help='Simulate contamination (either ancestral or modern human)')
    contam_subparser.add_argument("-vcf",help="[Required] path to vcf to simulate damage in", type=str, required = True)
    contam_subparser.add_argument("-out",help="path to output simulated vcf", type=str, default = "contaminate.vcf") 
    contam_subparser.add_argument("-mh",help="modern human flag", action = 'store_true')
    contam_subparser.add_argument("-ancestral",help="ancestral flag human flag", action = 'store_true')
    contam_subparser.add_argument("-targets", metavar='',help="Target individuals to simulate features on, comma separated list or json file", type=str, default = "")
    contam_subparser.add_argument("-rate", metavar='',help="rate at which to contaminate sites of target individuals", type=float, default=0.05)
    contam_subparser.add_argument("-length", metavar='',help="length of contaminating modern human fragments", type=int, default=1)
    contam_subparser.add_argument("-modern", metavar='',help="comma-separated list of individuals to contaminate with. If not provided, will look for 'contam' individuals in a json file provide with targets", type=str, default = "")

    # downsample VCF
    downsample_subparser = subparser.add_parser('downsample', help='Downsample VCF')
    downsample_subparser.add_argument("-vcf",help="[Required] path to vcf to simulate damage in", type=str, required = True)
    downsample_subparser.add_argument("-out",help="path to output simulated vcf", type=str, default = "downsample.vcf")
    downsample_subparser.add_argument("-num", metavar='',help="number of snps to downsample to", type=int, default=30_000)

    # add missingness
    missing_subparser = subparser.add_parser('missing', help='add missingness to VCF')
    missing_subparser.add_argument("-vcf",help="[Required] path to vcf to simulate damage in", type=str, required = True)
    missing_subparser.add_argument("-targets", metavar='',help="Target individuals to simulate features on", type=str, default = "")
    missing_subparser.add_argument("-out",help="path to output simulated vcf", type=str, default = "missing.vcf")
    missing_subparser.add_argument("-rate", metavar='',help="rate of missingness to simulate", type=float, default=0.1)

    # add depth and filter
    depth_subparser = subparser.add_parser('dpFilter', help='Downsample VCF')
    depth_subparser.add_argument("-vcf",help="path to vcf to simulate damage in", type=str, required = True)
    depth_subparser.add_argument("-out",help="path to output simulated vcf", type=str, default = "out.vcf")
    depth_subparser.add_argument("-mean", metavar='',help="mean depth to simulate", type=int, default=5)
    depth_subparser.add_argument("-variance", metavar='',help="variance of depth to simulate", type=int, default=2)
    # allowing varying disbributions has not been implemented yet
    # depth_subparser.add_argument("-distribution", metavar='',help="distribution to use to simulate depth from", type=str, default="normal")
    depth_subparser.add_argument("-missing", help="flag that indiciates to remove sites with 0 reads", action = "store_true")
    depth_subparser.add_argument("-annotate", help="flag that indiciates to only annotate with depth", action = "store_true")
    depth_subparser.add_argument("-targets", metavar='',help="Target individuals to simulate features on", type=str, default = "")

    args = parser.parse_args()

    # check if a simulation mode was provided
    if not args.mode:
        print("\nA simulation mode is required. Options are:")
        print("pseudohaploid, deaminate, contaminate, downsample, missing, dpFilter\n")

    else:
        # run the features (for now, only one can be called at a time)
        if args.mode == "pseudohaploid":
            sample_list = parse_indivs(args.targets, "target")
            print("Pseudohaplotyping")
            make_pseudohaploid(args.vcf, args.out, sample_list)
            print("finished pseudohaploid calls")
        
        if args.mode == "deaminate":
            sample_list = parse_indivs(args.targets, "target")
            if not (0 < args.rate < 1):
                raise Exception("deamination rate must be above 0 and less than 1")
            if args.proportion:
                if not (0 <= args.proportion <= 1):
                    raise Exception("If specified, the proportion of transitions sites must be above 0 and less than 1")
            print("Adding deamination")
            add_deam(args.vcf, args.out, sample_list, args.rate, args.proportion)
            print("Finished adding deamination")

        if args.mode == "contaminate":
            sample_list = parse_indivs(args.targets, "target")
            if args.ancestral and args.mh:
                raise Exception("Cannot specify both ancestral and modern human contamination flags")
            if args.ancestral:
                print(f"Adding ancestral contamination at {args.rate} rate")
                add_anc_contamination(args.vcf, args.out, sample_list, args.rate)
                print("Finished adding ancestral contamination")
            # check if a contaminating modern human population was provided as an argument or in the 
            elif args.mh:
                if args.modern:
                    modern_list = parse_indivs(args.modern, "contam")
                else:
                    raise Exception("Modern human population is necessary for adding modern human contamination")

                print(f"Adding modern contamination from {modern_list} at {args.rate} rate")
                add_mh_contamination(args.vcf, args.out, sample_list, modern_list, args.rate, args.length)
                print("Finished adding modern contamination")
            else:
                raise Exception("must specify either anc for ancestral contamination or mh for modern human contamination")

        if args.mode == "downsample":
            if args.num <=0:
                raise Exception("num of positions must be > 0")
            else:
                print(f"Beginnning downsampling to {args.num} positions")
                downsample(args.vcf, args.out, args.num)
                print("Finished downsampling")

        if args.mode == "missing":
            sample_list = parse_indivs(args.targets, "target")
            if not (0 < args.rate < 1):
                raise Exception("rate must be between 0 and 1")
            else:
                print(f"Adding missingness at rate: {args.rate}")
                add_missingness(args.vcf, args.out, sample_list, args.rate)
                print("Finished adding missingness")

        if args.mode == "dpFilter":
            sample_list = parse_indivs(args.targets)
            if args.missing:
                print("Adding depth annotions, including missingness, and creating false homozygotes")
                fun = "pos_depth_all"
            elif args.annotate:
                print("Adding depth annotations")
                fun = "pos_depth_only"
            else: # by default - annotate and add false homozygotes but don't induce missingness
                print("Adding depth annotations and creating false homozygotes")
                fun = "pos_depth_homo"
        
            add_depth(args.vcf, args.out, sample_list, args.mean, args.variance, fun) # args.distribution, 
            print("Finished adding depth")


if __name__ == "__main__":
    print("starting")
    main()

# variables for testing:
# vcf_path = "/global/scratch/users/sarahj32/sim-fstats/data/medium/original/full/vcf/sim1_original_full_filtered_17.vcf"
# sample_list = ["pop_1", "pop_2", "pop_3"]
# new_vcf = "test.vcf"
# sample_list = []

# other example vcfs:
# bcftools view /global/scratch/p2p3/pl1_moorjani/sarahj32/superarchaic_introgression/superarchaic_analysis/data/VCFs/hg38/AltaiDen_merged_hg38_chr21.v 10000 > test_multi_geno.vcf 
# vcf_path = "test_multi_geno.vcf"
# sample_list = ["AltaiNea"]
