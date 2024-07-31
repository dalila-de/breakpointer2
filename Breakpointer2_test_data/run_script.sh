#!/bin/bash

python Breakpointer2.py species1_to_species2.paf species1_500kb_window_size.bed species2_500kb_window_size.bed sp1 sp2 100000 5000000 30 100000 sp1_sp2

#sp1 and sp2 should be used to differengtiate species in the final results
#100000 is the minimum length of the alignment to be processed downstream
#5000000 is the distance required between two genomic regions within a chromosome needed to call a breakpoint in the alignment against the chromosome of the other species
#30 is the minimum quality score of the alignment to be kept in downstream analyses
#100000 is the number of trials in one-tailed permutation test
#sp1_sp2 is the name of the file generated containing all the detected breakpoints
