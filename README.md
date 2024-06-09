# Breakpointer2

Breakpointer2 is a program that parses insulation score data relative to chromosomal breakpoints between two closely-related species. This  allows us to explore how fixed rearrangements throughout evolution relate to topological units of the genome and their boundaries (topologically associated domains, TADs) in non-model organisms.

Breakpointer2 is a Python3-based program that uses:
1. genome-genome alignments of chromosome-scale genomes in .paf format
2. information on insulation scores (topological information- as described here https://vaquerizaslab.github.io/fanc/fanc-executable/fanc-analyse-hic/domains.html)

# More details on the program
The script *Breakpointer2.py* was used to generate data for a Master Thesis available here https://utheses.univie.ac.at/detail/67483.

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
