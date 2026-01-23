# ArchSim
## A package for simulating archaic sediment and skeletal DNA features in VCF's

---
1. [Installation](#Installation)
2. [Usage](#Usage)
3. [Simulation Details](#Simulation Details)
3. [Examples](#Examples)
---

## Installation

Clone the repository:
```bash
git clone https://github.com/sarahaj32/archSim.git
```

Dependencies: Python3, json
Ensure that python3 is installed

## Usage 
```note

The package is run on command line 
### Required arguments for all simulation modes:
**simulation mode** (options: psuedohaploid, deaminate, contaminate, downsample, missing, dpFilter)

-**vcf** <vcf_path>

### Optional arguments for all simulation modes:
-**out** [name and path to output file] default = out.vcf

-**targets** [comma separated list of target populations] default = all
or a json file with "targets": 

### Module specific arguments:

> psuedohaploid         
    -vcf                [required] path to the vcf to simulate data in
    -out                outputfile (defaults to pseudohaploid.vcf)
    -targets            a comma separated list of targets sample names, or path to a json file containing target sample names. These target samples are the ones that will become psuedohaploid

> deaminate:
    -vcf                [required] path to the vcf to simulate data in
    -out                outputfile (defaults to deaminate.vcf)
    -targets            a comma separated list of targets sample names, or path to a json file containing target sample names. These target samples are the ones that will become psuedohaploid
    -rate               Rate at which deamination will be simulated in the sites of the target population, between 0 and 1. Default = 0.05. (at this rate, homozygous reference gentoypes at transition sites in the selected individuals will convert to heterozygous)
    -proportion         The proportion of sites that will be specified as transitions. If this argument is not included, transition sites will be identified from VCF alleles

> contaminate
    -vcf                [required] path to the vcf to simulate data in
    -out                outputfile (defaults to contaminate.vcf)
    -targets            a comma separated list of targets sample names, or path to a json file 
    -ancestral          flag that indicates that contaminated derived sites will be converted to heterozygous
    -rate               the rate at which contamination is induced, between 0 and 1. Default = 0.05.(at this rate, homozygous alternative gentoypes in the selected individuals will convert to heterozygous)
    -mh                 flag that indicates that contaminated sites will be contributed from a given set of modern human individuals
    -modern             comma separated list of individuals to contaminate with, or path to a json file
    -length             length of contaminating modern human fragments. Default = 1

> downsample
    -vcf                [required] path to the vcf to simulate data in
    -out                outputfile (defaults to downsample.vcf)
    -num                number of snps to downsample to. Default = 30,000

#### missing
-**rate** rate of missingness to simulate. Default = 0.1

#### dpFilter 
-**mean** mean depth to simulate. Default = 5

-**variance** variance of depth to simulate. Default = 2

-**missing** flag that to convert genotypes at sites with 0 reads to missing (./.)

-**annotate** flag that indiciates to only annotate with depth (do not add missing genotypes, do not change genotypes to simulate homozygous bias)

*** Currently, only one mode can be run at a time (adding multiple flags will run but will overwrite each file with each newly created file)
```

## Example Data:
To provide examples of archSim's usage, we simulated human data under a simple demography with an outgroup (representing African modern humans), a deeply divergent lineage (representing archaic hominins), and an admixed population between the outgroup and diverged lineage (representing non-African modern humans). The admixed population was sampled shortly after admixture and at present day, representing an ancient individual and modern human.

vcf=simulated_data_21.vcf

There is also an example json file, which identifies all of the "admixed" individuals as targets, and all of the "modern human" indivduals as contamination sources.


Example adding 1% modern human contamination in blocks of 500 bp to populations Vindija and Denisova from populations AFR and mh_contam:

python src/main.py contaminate -mh -vcf data/archaic_admix_MHtarget/original/full/vcf/original_filtered_21.vcf -out ./test/mhcontam_21.vcf -rate 0.01 -length 500 -target Vindija,Denisova -modern AFR,mh_contam

Example adding deamination to all individuals:

python src/main.py deaminate -vcf $vcf -out ./test/deaminated_21.vcf -r 0.5

Example making all populations pseudohaploid:

python src/main.py pseudohaploid -vcf $vcf -out ./test/pseudohaploid_21.vcf 

Example making only populations 0 and 10 pseudohaploid:

python src/main.py pseudohaploid -vcf $vcf -out ./test/pseudohaploid_pop0pop10_21.vcf -target pop_0,pop_10

## Simulation details
# pseudohaploid:
At every position, an allele will be randomly selected from each target individual, and the genotype will become homozygous for that allele. Homozygous positions will be unaffected, but heterozygous positions will be converted to homozygous reference or homozygous alternative with an even probability. If the input is phased VCF, the output pseudohaploid VCF will still be phased. 

# Example:
Let's simulate the case where only the target individuals are converted to psuedohaploid:
```note
python src/main.py pseudohaploid -vcf $vcf -targets test/individuals_all.json -out ./test/simulated_human_pseudohaploid_21.vcf 
```

We could also simulate a case where all individuals are converted to pseudohaploid:
```note
python src/main.py pseudohaploid -vcf $vcf -out ./test/simulated_human_pseudohaploid_all_21.vcf 
```

If we only have 2 admixed individuals that we want to be pseudohaploid, we can specify those in the command line without needing a json file:
We could also simulate a case where all individuals are converted to pseudohaploid:
```note
python src/main.py pseudohaploid -vcf $vcf -targets admix_1,admix_2 -out ./test/simulated_human_pseudohaploid_2Admix_21.vcf 
```


# deaminate:
We induce genotype errors at transition sites, representative of errors likely to be seen from deamination, by converting homozygous reference calls to heterozygous calls at the specified rate (0.05 by default). By default, transition sites are determined from the VCF alleles as C/T, T/C, A/G, G/A sites. However, if the "-proportion" argument is specified, transition sites will be determined as that proprotion of all sites. If the input genotypes are phased, the output genotypes will still be phased, and any new alternative alleles will be randomly assigned to either chromosome.  

# Example:
Let's simulate deamination-related errors in the target individuals only:
```note
python src/main.py deaminate -vcf $vcf -targets test/individuals_all.json -out ./test/simulated_human_deaminated_21.vcf 
```

# contaminate (ancestral)
We simulate faunal contamination by converting derived alleles to the ancestral at the specified rate (0.05 by default). If a site is contaminated and converted to heterozygous, and the orignal genotype was phased, the ancestral allele will be randomly assinged between chromosomes. Multiallelic positions are skipped in this step and removed from the output file. 

Note - this assumes that the VCF is polarized so that the reference allele is ancestral and the alternative allele is derived. If the VCF is not polarized, ancestral contamination can be simulated using the "modern human" simulation mode (see next section).

# Example:
Let's simulate ancestral contamination in the target individuals only, at the default rate:
```note
python src/main.py contaminate -vcf $vcf -targets test/individuals_all.json -ancestral -out ./test/simulated_human_ancContam_21.vcf 
```

Let's simulate a large amount of ancestral contamination in the target individuals only, at a rate of 20%:
```note
python src/main.py contaminate -vcf $vcf -targets test/individuals_all.json -ancestral -rate 0.2 -out ./test/simulated_human_highAncContam_21.vcf 
```

# contaminate (modern human)
We simulate modern human contamination by replacing the genotypes of the target individuals with the genotypes of a randomly selected contaminating individual, at at specified rate (default = 0.05). At each position in each individual, the contaminating genotype will be randomly selected from genotypes at that position of all the contaminating individuals. 

Contamination from the same individual can be simulated into the same target individual in chunks by specifying the length argument. If the `-length 1000` argument is added, then at the given rate the target individual's next 1000 genotypes will be replaced by the next 1000 genotypes of one of the contaminating individuals. Note that if all VCF snps are at least 1000 bp apart, the `-length 1000` will yield the same results as not specifying the length (or `-length 1`). If VCF snps are closer than 1000 bp, then the resulting contamination amount will be slightly higher than the provided rate, since the rate specifies how often the chunk of 1000 bp are replaced, but each contamination event contributes more than 1 contaminating genotype

While this module was created to simulate modern human contamination, any individual can be used for the contaminating individual. For example, if the VCF is not polarized, an ancestral individual could be added to the VCF and specified as the "contaminating" individual.

# Example:
Let's simulate contamination in the target individuals from the specified contaminating individuals, at the default rate:
```note
python src/main.py contaminate -vcf $vcf -mh -targets test/individuals_all.json -modern test/individuals_all.json -out ./test/simulated_human_mhContam_21.vcf 
```

Let's simulate a contamination from the specified individuals into the target individuals, in chunks of 1000 bp, at a rate of 0.01:
```note
python src/main.py contaminate -vcf $vcf -mh -targets test/individuals_all.json -modern test/individuals_all.json -rate 0.01 -length 1000 -out ./test/simulated_human_chunkMhContam_21.vcf 
```
Again, the actual fraction of contaminating genotypes will be slightly higher than 1% if any of the snps are within 1000 bp of another snp.

Let's simulate a contamination from a single individual (Afr_8) into the target individuals, at a rate of 8%:
```note
python src/main.py contaminate -vcf $vcf -mh -targets test/individuals_all.json -modern Afr_8 -rate 0.08 -out ./test/simulated_human_Afr8MhContam_21.vcf 
```

At every position in all specified individuals, at a specified rate contamination may be added. When contamination is added:
1. a modern human individual is randomly selected
2. the genotype of the contaminatED individual is swapped for the genotype of the contaminatING individual
3. The genotypes are similarly swapped for the next 'length' basepairs (ie: if length = 1000, any SNPs in the next 1000 basepairs will receive the genotype of the contaminating individual). This is tracked separately for each contaminatED individual, so all contaminating "chunks" are independent
4. After this contaminating "chunk" the process repeats, and at the same rate again at the next position, that individual may be contaminated

# Example:
Let's simulate deamination-related errors in the target individuals only:
```note
python src/main.py deaminate -vcf $vcf -targets test/individuals_all.json -out ./test/simulated_human_deaminated_21.vcf 
```

# contaminate
At every position in all specified individuals, at a specified rate contamination may be added. When contamination is added:
1. a modern human individual is randomly selected
2. the genotype of the contaminatED individual is swapped for the genotype of the contaminatING individual
3. The genotypes are similarly swapped for the next 'length' basepairs (ie: if length = 1000, any SNPs in the next 1000 basepairs will receive the genotype of the contaminating individual). This is tracked separately for each contaminatED individual, so all contaminating "chunks" are independent
4. After this contaminating "chunk" the process repeats, and at the same rate again at the next position, that individual may be contaminated

> modern human contamination
At every position in all specified individuals, at a specified rate contamination may be added. When contamination is added:
1. a modern human individual is randomly selected
2. the genotype of the contaminatED individual is swapped for the genotype of the contaminatING individual
3. The genotypes are similarly swapped for the next 'length' basepairs (ie: if length = 1000, any SNPs in the next 1000 basepairs will receive the genotype of the contaminating individual). This is tracked separately for each contaminatED individual, so all contaminating "chunks" are independent
4. After this contaminating "chunk" the process repeats, and at the same rate again at the next position, that individual may be contaminated
```