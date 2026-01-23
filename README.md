# ArchSim
## A package for simulating archaic sediment and skeletal DNA features in VCF's

### Contact:
sarahj32@berkeley.edu

---
1. [Installation](#installation)
2. [Usage](#usage)
3. [Simulation Details and Examples](#simulation-details-and-examples)
      - [Psuedohaploid](#pseudohaploid)
      - [Deaminate](#deaminate)
      - [Contaminate](#contaminate)
      - [Downsample](#downsample)
      - [Missing](#missing)
      - [dpFilter](#dpfilter)
4. [Example Data](#example-data)
5. [Combining Features](#combining-features)
---

## Installation

Clone the repository:
```bash
git clone https://github.com/sarahaj32/archSim.git
```

Dependencies: Python3

## Usage 
archSim is run on command line by calling: `python src/main.py`

Parameter details: (this can also be shown for each method with the help flag: `python src/main.py deaminate -h`)

```note
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

> mising
    -vcf                [required] path to the vcf to simulate data in
    -targets            a comma separated list of targets sample names, or path to a json file 
    -out                outputfile (defaults to missing.vcf)
    -rate               rate at which genotypes will be converted to missing (./.). Default = 0.1

> dpFilter
    -vcf                [required] path to the vcf to simulate data in
    -targets            a comma separated list of targets sample names, or path to a json file 
    -out                outputfile (defaults to dpFilter.vcf)
    -mean               the mean depth to simulate. Default = 5
    -variance           the variance of depth to simulate. Default = 2
    -dist               distribution to sample reads around. Default = nbinom. Other options = poisson, normal
    -bias               Reference bias threshold. Default = 0.55
    -dropout            Minimum number of reads necessary to have a heterozygous call. Default = 3

```

## Example Data:
To provide examples of archSim's usage, we simulated human data under a simple demography with an outgroup (representing African modern humans), a deeply divergent lineage (representing archaic hominins), and an admixed population between the outgroup and diverged lineage (representing non-African modern humans). The admixed population was sampled shortly after admixture and at present day, representing an ancient individual and modern human.

`vcf=test/simulated_data.vcf`
In this dataset, the following names correspond to the following populations:
Afr - African outgroup
admix - modern humans with archaic admixture, sampled shortly after admixture
AMH - anatomically modern humans, sampled prior to archaic admixture
mh_contam - modern humans with archaic admixture, sampled at present time

There is also an example json file `test/individuals_all.json`, which identifies all of the "admixed" individuals as targets, and all of the "modern human" indivduals as contamination sources.

With these files, you can replicate all of the below examples.

## Simulation Details and Examples

## pseudohaploid:
At every position, an allele will be randomly selected from each target individual, and the genotype will become homozygous for that allele. Homozygous positions will be unaffected, but heterozygous positions will be converted to homozygous reference or homozygous alternative with an even probability. If the input is phased VCF, the output pseudohaploid VCF will still be phased. 

### Example:
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

## deaminate:
We induce genotype errors at transition sites, representative of errors likely to be seen from deamination, by converting homozygous reference calls to heterozygous calls at the specified rate (0.05 by default). By default, transition sites are determined from the VCF alleles as C/T, T/C, A/G, G/A sites. However, if the "-proportion" argument is specified, transition sites will be determined as that proprotion of all sites. If the input genotypes are phased, the output genotypes will still be phased, and any new alternative alleles will be randomly assigned to either chromosome.  

### Example:
Let's simulate deamination-related errors in the target individuals only:
```note
python src/main.py deaminate -vcf $vcf -targets test/individuals_all.json -out ./test/simulated_human_deaminated_21.vcf 
```

## contaminate 
## (ancestral):
We simulate faunal contamination by converting derived alleles to the ancestral at the specified rate (0.05 by default). If a site is contaminated and converted to heterozygous, and the orignal genotype was phased, the ancestral allele will be randomly assinged between chromosomes. Multiallelic positions are skipped in this step and removed from the output file. 

Note - this assumes that the VCF is polarized so that the reference allele is ancestral and the alternative allele is derived. If the VCF is not polarized, ancestral contamination can be simulated using the "modern human" simulation mode (see next section).

### Example:
Let's simulate ancestral contamination in the target individuals only, at the default rate:
```note
python src/main.py contaminate -vcf $vcf -targets test/individuals_all.json -ancestral -out ./test/simulated_human_ancContam_21.vcf 
```

Let's simulate a large amount of ancestral contamination in the target individuals only, at a rate of 20%:
```note
python src/main.py contaminate -vcf $vcf -targets test/individuals_all.json -ancestral -rate 0.2 -out ./test/simulated_human_highAncContam_21.vcf 
```

## (modern human):
We simulate modern human contamination by replacing the genotypes of the target individuals with the genotypes of a randomly selected contaminating individual, at at specified rate (default = 0.05). At each position in each individual, the contaminating genotype will be randomly selected from genotypes at that position of all the contaminating individuals. 

Contamination from the same individual can be simulated into the same target individual in chunks by specifying the length argument. If the `-length 1000` argument is added, then at the given rate the target individual's next 1000 genotypes will be replaced by the next 1000 genotypes of one of the contaminating individuals. Note that if all VCF snps are at least 1000 bp apart, the `-length 1000` will yield the same results as not specifying the length (or `-length 1`). If VCF snps are closer than 1000 bp, then the resulting contamination amount will be slightly higher than the provided rate, since the rate specifies how often the chunk of 1000 bp are replaced, but each contamination event contributes more than 1 contaminating genotype

While this module was created to simulate modern human contamination, any individual can be used for the contaminating individual. For example, if the VCF is not polarized, an ancestral individual could be added to the VCF and specified as the "contaminating" individual.

### Example:
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

## downsample:
We downsample the vcf by first identifying the number of data lines in the file (lines that do not start with #), and randomly select the specified "num" number of lines (30,000 by default). The output file will contain these downsampled lines, as well as all header lines.

### Example:
Let's downsample our VCF to 10,000 SNPs:
```note
python src/main.py downsample -vcf $vcf -num -out ./test/simulated_human_downsampled_21.vcf 
```
The resulting VCF has 10,008 rows, the 8 header rows from the original VCF and the 10,000 downsampled data rows. 

## missing
We simulate missing data by converting genotypes of the specified target individuals to missing (./.) at the given rate (default = 0.1)

### Example:
Let's add 10% missingness (the default) to all target individuals:
```note
python src/main.py missing -vcf $vcf -targets test/individuals_all.json -out ./test/simulated_human_missing_21.vcf 
```

Let's add 5% missingness to two target individuals only: admix_1, and admix_3:
```note
python src/main.py missing -vcf $vcf -targets admix_1,admix_3 -rate 0.05 -out ./test/simulated_human_fewMissing_21.vcf 
```

## dpFilter:
We simulate depth for each target individual and position from a distribution with a provided mean and variance. Multiallelic positions are skipped in this step and removed from the output file. By default, depth (DP) is sampled from a negative binomial distribution wtih mean 5 and variance 3, though these parameters can be updated with input arguments. The gentoype is set to missing (./.) if there are no reads. At homozygous sites with reads, all reads are assigned to the homozygous allele. At heterozygous sites, the reads (DP) are distributed across alleles in the following way:
- If the number of reads is less than the dropout threshold (default 3, but can be specified with the `dropout` argument), then the genotype is converted to homozygous. The allele is selected by random, following the reference bias (default 0.55 but cant be specified with the `bias` argument). For example, if bias = 0.55, 55% of the time the position will become homozygous reference, and 45% of the time the position will become homozygous alternative
- If the number of reads are gereater than the dropout threshold, they are distributed following a binomial distribution. if either allele has no reads assigned then the genotype is homozygous for the other allele. Otherwise the reads supporting each allele are output

### Example:
Let's simulate depth in all target individuals, with the default mean, variance, distribution, dropuout, and bias
```note
python src/main.py dpFilter -vcf $vcf -targets test/individuals_all.json  
```

Let's simulate depth in two target individuals only (admix_2, admix_4), with mean depth of 10 and variance of 2, a dropout depth of 4 and reference bias of 0.7
```note
python src/main.py dpFilter -vcf $vcf -targets admix_2,admix_4 -mean 10 -variance 2 -dropout 4 -bias 0.7  
```

## Combining Simulated Features:
Each run of archSim only will simulate one feature. However, features can be combined by applying each method to the output of the previous method. Be sure to specify the output file names.

Let's simulate a dataset where the target individuals have 5% deamination and 10% contamination from two modern human individuals: Afr_1 and Afr_9. This dataset will also be downsampled to 10,000 positions, and all individuals will be pseudohaploid.

```note
python src/main.py deaminate -vcf $vcf -targets test/individuals_all.json -out ./test/simulated_human_deaminated_21.vcf

python src/main.py contaminate -vcf ./test/simulated_human_deaminated_21.vcf -targets test/individuals_all.json -mh -modern Afr_1,Afr_9 -r 0.1 -out ./test/simulated_human_deamContam_21.vcf

python src/main.py downsample -vcf ./test/simulated_human_deamContam_21.vcf -num 10000 -out ./test/simulated_human_deamContamDs_21.vcf

python src/main.py pseudohaploid -vcf ./test/simulated_human_deamContam_21.vcf -out ./test/simulated_human_deamContamDsPsuedohap_21.vcf
```

The output dataset should be a valid VCF that is compatible with other downstream tools, like bcftools


