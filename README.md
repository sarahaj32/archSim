## ArchSim

Dependencies: Python3
Required arguments:
simulation mode
-vcf <vcf_path>

### Optional arguments for all simulation modes:
-out [name and path to output file] default = out.vcf

-targets [comma separated list of target populations] default = all

### Module specific arguments:

#### deaminate:
-rate rate at which deamination will be simulated in the sites of the target population. Default = 0.05
(at this rate, sites in the selected population will convert to ancestral (50% to het, 50% to homozygous ancestral)

#### contamination:
-rate rate at which to contaminate sites in target population. Default = 0.05
-ancestral flag that indicates that contaminated sites will be converted to ancestral alleles
-mh flag that indicates that contaminated sites will be contributed from modern human populations
-modern [required if -mh flag] name of modern human population to contaminate from
-length length of contaminating modern human fragments. Default = 1000

#### downsample
-num number of snps to downsample to. Default = 30,000

#### missing
-rate rate of missingness to simulate. Default = 0.1

#### dpFilter 
-mean mean depth to simulate. Default = 5
-variance variance of depth to simulate. Default = 2
-missing flag that to convert genotypes at sites with 0 reads to missing (./.)
-annotate flag that indiciates to only annotate with depth (do not add missing genotypes, do not change genotypes to simulate homozygous bias)

*** Currently, only one of these can be done at a time (adding multiple flags will run but will overwrite each file with each newly created file)

## Examples:
vcf = data/archaic_admix_MHtarget/original/full/vcf/original_filtered_21.vcf

Example adding deamination to all individuals
python src/main.py deaminate -vcf data/archaic_admix_MHtarget/original/full/vcf/original_filtered_21.vcf -o ./test/deaminated_21.vcf -r 0.5

Example pseudohaplotyping all individuals
python src/main.py pseudohaploid -vcf data/archaic_admix_MHtarget/original/full/vcf/original_filtered_21.vcf -o ./test/pseudohaploid_21.vcf 

Example pseudohaplotyping only populations 0 and 10
python src/main.py pseudohaploid -vcf data/archaic_admix_MHtarget/original/full/vcf/original_filtered_21.vcf -o ./test/pseudohaploid_pop0pop10_21.vcf -target pop_0,pop_10

Example adding 1% modern human contamination in blocks of 500 bp to populations 0 and 10 from populations 2 and 3
python src/main.py contaminate -mh -vcf data/archaic_admix_MHtarget/original/full/vcf/original_filtered_21.vcf -o ./test/mhcontam_pop0pop10 -rate 0.01 -length 500 -target pop_0,pop_10 -modern pop_2,pop_3