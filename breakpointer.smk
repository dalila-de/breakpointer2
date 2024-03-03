"""
Breakpointer2
"""
import os
configfile= "config.yaml"
config["tool"] = "Breakpointer"

filepath = os.path.dirname(os.path.realpath(workflow.snakefile))

rule all:
    input:
        #

# Define a function to generate input paths dynamically
def get_bed_file(config, sp1, sp2):
    return f"/input/bed_files/{sp1}_against_{sp2}.bed"

"""
Only chromosome-scale genomes .fasta files prepared like on the link https://www.metagenomics.wiki/tools/fastq/subset
"""

rule chr_id:
    input:
        genome = config["species"][wildcards.sample]["genome"]
        num_chr = config["species"][wildcards.sample]["num_chr"]
    output:
        chr_id = "/input/chr_id_genomes/{wildcards.sample}.txt"
    shell:
        "zcat {input.genome} | grep '>'  | cut -f 1 -d ' ' |  sed 's/>//g' | head -n {input.num_chr} > {output.chr_id}"

rule chr_scale_only:
    input:
        genome = config["species"][wildcards.sample]["genome"],
        chr_id = "/input/chr_id_genomes/{wildcards.sample}.txt"
    output:
        chr_scale = "/input/chr_scale_genomes/{wildcards.sample}.fasta"
    shell:
        "seqtk subseq {input.genome} {input.chr_id} > {output.chr_scale}"

rule minimap2:
    input:
        chr_scale_1 = "/input/chr_scale_genomes/{sample1}.fasta",
        chr_scale_2 = "/input/chr_scale_genomes/{sample2}.fasta"
    output:
        paf = "/input/paf_files/{sample1}_vs_{sample2}.paf"
    params:
        seq_arg = config["seq_arg"]
    shell:
        "minimap2 {params.seq_arg} {input.chr_scale_1} {input.chr_scale_2} > {output.paf}"
        
rule breakpointer:
    input:
        input_paf = config["input_paf"],
        sp1 = config["sp1"],
        sp2 = config["sp2"],
        len_co = config["len_co"],
        len_bp = config["len_bp"],
        q = config["q"],
        bed_file_1 = lambda wildcards: get_bed_file(config, wildcards.sp1, wildcards.sp2),
        bed_file_2 = lambda wildcards: get_bed_file(config, wildcards.sp2, wildcards.sp1)
    output:
        dynamic("{plot}.png", plot=get_plot_names()),
        overall_table = "/output/breakpointer2/{sp1}+'_vs_'+{sp2}+'_breakpoints.tsv'"
    shell:
        "python {filepath}/breakpointer2.py {input.input_paf} {input.sp1} {input.sp2} {input.len_co} {input.len_bp} {input.q}"

rule clean_bed:
    input:
        bed_file_1 = "/input/bed_files/{sp1}_against_{sp2}.bed",
        bed_file_2 = "/input/bed_files/{sp2}_against_{sp1}.bed"
    output:
        cleaned_bed_1 = "/output/breakpointer2/{sp1}_against_{sp2}.bed",
        cleaned_bed_2 = "/output/breakpointer2/{sp2}_against_{sp1}.bed"
    shell:
        """
        cut -f2- <{input.bed_file_1} | sort -k1,1 -k2,2n >{output.cleaned_bed_1}
        cut -f2- <{input.bed_file_2} | sort -k1,1 -k2,2n >{output.cleaned_bed_2}
        """

rule bedtools_intersect_pairs:
    input:
        cleaned_bed_1 = "/output/breakpointer2/{sp1}_against_{sp2}.bed",
        cleaned_bed_2 = "/output/breakpointer2/{sp2}_against_{sp1}.bed"
    output:
        intersected_bed = "/output/breakpointer2/intersected_{sp1}_{sp2}.bed"
    shell:
        "bedtools intersect -a {input.bed_file_1} -b {input.bed_file_2} -wa -wb >{output.intersected_bed}"

        