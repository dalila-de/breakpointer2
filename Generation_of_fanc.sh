#!/bin/bash
#SBATCH --job-name=the_pairs
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL

#First the variables
PREFIX=
THREAD=
ASSEM=
R1=
R2=
MYBIN=

#Then load Chromap to index the genome and map the Hi-C reads
conda activate chromap

chromap -i -r ${ASSEM} -o ${PREFIX}.index
for qual in {0,2,10}
do
chromap -t ${THREAD} --preset hic -x ${PREFIX}.index -r ${ASSEM} -1 ${R1} -2 ${R2} -q ${qual} -o ${PREFIX}_${qual}.pairs
done

#Use bioawk to get a .txt file with all the chromosome sizes
conda deactivate
conda activate base

bioawk -cfastx '{printf("%s\t%d\n", $name, length($seq))}' ${ASSEM} > ${PREFIX}_chrom_size.txt

#Use pairix to compress the .pairs file and create a .gz file
${MYBIN}pairix/bin/bgzip -f ${PREFIX}_0.pairs

#Generate the .cool file
conda activate cooler
cooler cload pairix --nproc 4 ${PREFIX}_chrom_size.txt:5000 ${PAIRS}.gz ${PREFIX}.cool

#If you wish to visualize the .cool file, you can use hicexplorer, and also generate a normalized and balanced .cool file, but it is not necessary
conda activate hicexplorer
hicCorrectMatrix diagnostic_plot -m ${PREFIX}.cool -o ${PREFIX}.png
hicNormalize -m ${PREFIX}.cool --normalize norm_range -o ${PREFIX}_normalized.cool

conda activate cooler
cp ${PREFIX}_normalized.cool ${PREFIX}_normalized_balanced.cool
cooler balance --force ${PREFIX}_normalized_balanced.cool

#Now generate a .fanc file from the .cool file (fanc file has the suffix .hic, but is not the same as Juicer .hic files?)
module load fanc
fanc from-cooler ${PREFIX}.cool ${PREFIX}.hic
fanc from-cooler ${PREFIX}_normalized_balanced.cool ${PREFIX}_normalized_balanced.hic

#Now, generate insulation score files from the .fanc files
fanc insulation ${PREFIX}.hic ${PREFIX}.interactions -g -w 100000 250000 500000 750000 1000000 -o bed
fanc insulation ${PREFIX}.hic ${PREFIX}_normalized_balanced.interactions -g -w 100000 250000 500000 750000 1000000 -o bed