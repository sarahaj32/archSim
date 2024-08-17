# rule make_selection_files:
#     input:
#     output:
#         expand("data/{timescale}/selection_files/{select}.tsv", select = selection)
#     run:
#         for s in selection:
#             file = f"data/{timescale}/selection_files/{s}.tsv"
#             samples = selection_dict[s]
#             samples = [f"pop_{x}_{i}" for x in samples for i in range(10)]
#             samples_out = "\n".join(samples)
#             with open(file, "w") as f:
#                 f.write(samples_out)

# rule make_select_vcfs:
#     input:
#         vcf = "data/{timescale}/contaminated/full/vcf/contaminated_filtered_{chrom}.vcf",
#         sample_list = "data/{timescale}/selection_files/{select}.tsv"
#     output:
#         "data/{timescale}/contaminated/full/vcf/contaminated_{select}_filtered_{chrom}.vcf"
#     shell:
#         """
#         bgzip {input.vcf} > {input.vcf}.gz
#         tabix -p vcf {input.vcf}.gz
#         bcftools view {input.vcf}.gz -S {select}.tsv | head -n 1000 > filtered_file
#         """

rule combine_chromosomes:
    # combine all non-downsampled chromosome vcfs
    # create a file with the SNP positions so we can quickly downsample
    input:
        expand("data/{timescale}/{{param}}/full/vcf/{{param}}_filtered_{chrom}.vcf", chrom = chromosomes, timescale = config["timescale"])
    output:
        data = "data/{timescale}/{param}/full/vcf/{param}_full.vcf",
        indices = "data/{timescale}/{param}/full/vcf/chr_pos_{param}.vcf"
    run: 
        input_list =  [f"data/{config['timescale']}/{wildcards.param}/full/vcf/{wildcards.param}_filtered_{chrom}.vcf" for chrom in chromosomes]
        out_file = f"data/{config['timescale']}/{wildcards.param}/full/vcf/{wildcards.param}_full.vcf"
        out_indices = f"data/{config['timescale']}/{wildcards.param}/full/vcf/chr_pos_{wildcards.param}.vcf"
        shell("module load bcftools; bcftools concat -o {out_file} {input_list}")
        shell("cut -f1-2 {out_file} > {out_indices}")



