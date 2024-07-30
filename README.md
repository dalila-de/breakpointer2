# Breakpointer2

Breakpointer2 is a program that parses insulation score data relative to chromosomal breakpoints between two closely-related species. This  allows us to explore how fixed rearrangements throughout evolution relate to topological units of the genome and their boundaries (topologically associated domains, TADs) in non-model organisms.

Breakpointer2 is a Python3-based program that uses:
1. genome-genome alignments of chromosome-scale genomes in .paf format
2. information on insulation scores (topological information- as described here https://vaquerizaslab.github.io/fanc/fanc-executable/fanc-analyse-hic/domains.html)

# Summary
The script *Breakpointer2.py* was used to generate data for a Master Thesis available here https://utheses.univie.ac.at/detail/67483. The program is explained in great detail in Chapter 4 *Scale-free investigation of rearrangements in the context of genome topology*.

If you want to generate the files that the script can execute, please run *minimap2* under default parameters (https://github.com/lh3/minimap2), and for the generation of insulation scores refer to the bash script available in the repo (*Generation_of_fanc.sh*).

To run the script *Breakpointer2.py*, please see the following parameters that should be optimized based on the genomes of the lineage you want to test:\
`--input_paf` is the pairwise alignment file between two chromosome-scale genomes you want to test\
`--input_ins_score_sp1` insulation score file in *.bed* format for species 1 (the first species in the *.paf* file)\
`--input_ins_score_sp2` insulation score file in *.bed* format for species 2 (the second species in the *.paf* file)\
`--sp1` is species 1 in *.paf* file\
`--sp2` is species 2 in *.paf* file\
`--len_co` is an argument specifying the alignment length cutoff in the *.paf* file\
`--len_bp` is a parameter that specifies how distant two points in the alignment have to be to call it a breakpoint\
`--q` is the quality cutoff of the alignments\
`--num_rounds` is the number of rounds for the one-tailed-permutation test\
`--output_fin_tab` output file

The snakemake script (*breakpointer2_v5.smk*) is the latest workflow incorporating this program.
The goal is to generalize the code and have a workflow generating data for multiple pairwise comparisons simultaneously. The current workflow outputs only the number of inversion events in the pairwise comparison. More than two input species can be used. Future efforts will incorporate insulation score data to infer topological information more broadly across the tree of life.

# System Requirements

The Python script *Breakpointer2.py* was tested on MacOS and the Life Science Compute Cluster (LiSC) of the University of Vienna. The script itself can be run on a personal computer, but to generate the necessary files, one needs a computing cluster.

Python packages that are required are Pandas v1.5.2. (Team, 2023), Matplotlib v3.6.3.(Hunter, 2007), numpy v1.24.1. (Harris et al., 2020), os, and statistics.

# Installation Guide

`git clone https://github.com/dalila-de/breakpointer2.git`

The script in the repository is ready to use as is after cloning.

# Detailed Instructions

The program was tested and performed analyses for the following genomes: GCF_006345805.1, GCF_001194135.2, GCA_951406725.2, and https://figshare.com/s/fa09f5dadcd966f020f3.

## Generation of genome-genome alignments:

For that analysis, we picked a reference genome- *Octopus vulgaris* GCA_951406725.2, and aligned it against the other three genomes using Minimap2 v2.24 (H. Li, 2018, 2021).

`minimap2 species1_genome.fasta species2_genome.fasta > species1_species2.paf`
`minimap2 species1_genome.fasta species3_genome.fasta > species1_species3.paf`

## Generation of insulation score files:

Necessary data:
- Genomes used in the alignments (*.fasta*)
- Hi-C reads (*.fastq*)

Each step should be done for each analyzed species.

Steps:

1. Use Chromap v0.2.3 (H. Zhang et al., 2021) for genome indexing and mapping of Hi-C reads.

`chromap -i -r species1_genome.fasta -o species1_genome.index`
`chromap -t 20 --preset hic -x species1_genome.index -r species1_genome.fasta -1 forward_hic_reads.fastq.gz -2 reverse_hic_reads.fastq.gz -q 0 -o species1_0.pairs`

2. Compress the *.pairs* file using Pairix v0.3.7 (Lee et al., 2022).

`pairix/bin/bgzip -f species1_0.pairs`

3. Generate a *.txt* file containing scaffold names and lengths using bioawk v1.0 (https://github.com/lh3/bioawk).

`bioawk -cfastx '{printf("%s\t%d\n", $name, length($seq))}' species1_genome.fasta > species1_chrom_size.txt`

4. Generate a *.cool* matrix using Cooler v0.9.1 (Abdennur & Mirny, 2020). A bin of 5000 worked well for our analyses.

`cooler cload pairix --nproc 4 species1_chrom_size.txt:5000 species1_0.pairs.gz species1.cool`

5. This is an optional step. If you wish to visualize the .cool file, you can use HiCExplorer toolkit v3.7.2 (Ramírez et al., 2018; Wolff et al., 2018, 2020), and also generate a normalized and balanced .cool file.

`hicCorrectMatrix diagnostic_plot -m species1.cool -o species1.png`
`hicNormalize -m species1.cool --normalize norm_range -o species_normalized.cool`
`#Here, you use Cooler again to create a balanced matrix`
`cp species1_normalized.cool species1_normalized_balanced.cool`
`cooler balance --force species1_normalized_balanced.cool`

6. Use FAN-C Toolkit (Kruse et al., 2020) to convert the *.cool* files to the compatible format for the downstream analyses.

`fanc from-cooler species1.cool species1.hic`

7. Generate insulation score files using FAN-C Toolkit (Kruse et al., 2020). Multiple windows are used, and optimization might be necessary depending on the genome size.

`fanc insulation species1.hic species1.interactions -g -w 100000 250000 500000 750000 1000000 -o bed`

## Running Breakpointer2

Now that the alignment (*.paf*) and insulation score files (*.bed*) have been generated, Breakpointer2 can be run on the pairwise data.

Example command:

`python /path/breakpointer2/Breakpointer2.py --input-paf species1_species2.paf --input_ins_score_sp1 species1.bed --input_ins_score_sp2 species2.bed --sp1 Ovu --sp2 Obi --len-co 100000 --len-bp 5000000 --q 30 --num-rounds 100000 --output_fin_tab species1_species2.tsv`

There are additional output files, including a statistical summary for the colocalization of breakpoints and borders of topologically associated domains.
A *.bed* file is also generated for the identified breakpoints in each species. This means that if you have a reference genome, you can view breaks in different species against this reference in IGV (https://igv.org/doc/desktop/#), or if you want to view the breakpoints relative to the Hi-C heatmap, one can visualize the *.bed* files on HiGlass (Kerpedjiev et al., 2018).

Future efforts are towards the automatization of quantification of shared and unique breakpoints in the lineages, as viewing the *.bed* files against a reference genome requires manual counting of shared or unique breakpoints.

