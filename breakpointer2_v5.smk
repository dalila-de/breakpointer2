import os
from Bio import SeqIO
import gzip
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
configfile: "config.yaml"
config["tool"] = "Breakpointer2"
config["analyses"] = {}

#to disable the warning when adding the breakpoints column to the original dataframe
pd.options.mode.chained_assignment = None

"""
Setting up the snakemake workflow: analyses part to create pairs in the workflow (similar to what is used in odp https://github.com/conchoecia/odp
"""

all_species = list(sorted(list(config["species"].keys())))
# could use permutations
for i in range(len(all_species)):
    for j in range(i+1, len(all_species)):
        species1 = all_species[i]
        species2 = all_species[j]
        analysis_name =  f"{species1}_{species2}"
        config["analyses"][analysis_name] = {"species1": species1,
                                             "species2": species2}

file_targets = []
for this_analysis in config["analyses"]:
    sp1 = config["analyses"][this_analysis]["species1"]
    sp2 = config["analyses"][this_analysis]["species2"]
    file_targets.append(config["tool"] + f"/input/paf_frag/{sp1}_to_{sp2}.raw.paf")
    file_targets.append(config["tool"] + f"/input/paf_frag/{sp2}_to_{sp1}.raw.paf")
    file_targets.append(config["tool"] + f"/input/paf_compl/{sp1}_to_{sp2}.paf")
    file_targets.append(config["tool"] + f"/input/paf_compl/{sp2}_to_{sp1}.paf")
    file_targets.append(config["tool"] + f"/input/processed/{sp1}_to_{sp2}/{sp1}_to_{sp2}_filtered.paf")
    file_targets.append(config["tool"] + f"/input/processed/{sp2}_to_{sp1}/{sp2}_to_{sp1}_filtered.paf")
    file_targets.append(config["tool"] + f"/input/processed/{sp1}_to_{sp2}/{sp1}_to_{sp2}_FET_test.tsv")
    file_targets.append(config["tool"] + f"/input/processed/{sp2}_to_{sp1}/{sp2}_to_{sp1}_FET_test.tsv")
    file_targets.append(config["tool"] + f"/output/{sp1}_to_{sp2}_run/{sp1}_vs_{sp2}_breakpoints.tsv")
    file_targets.append(config["tool"] + f"/output/{sp2}_to_{sp1}_run/{sp2}_vs_{sp1}_breakpoints.tsv")
    file_targets.append(config["tool"] + f"/input/raw_bedlike_files/{sp1}_to_{sp2}_run/{sp1}_{sp2}.bed")
    file_targets.append(config["tool"] + f"/input/raw_bedlike_files/{sp1}_to_{sp2}_run/{sp2}_{sp1}.bed")
    file_targets.append(config["tool"] + f"/input/raw_bedlike_files/{sp2}_to_{sp1}_run/{sp1}_{sp2}.bed")
    file_targets.append(config["tool"] + f"/input/raw_bedlike_files/{sp2}_to_{sp1}_run/{sp2}_{sp1}.bed")

wildcard_constraints:
    sp="[A-Za-z0-9-.]+",
    sp1="[A-Za-z0-9-.]+",
    sp2="[A-Za-z0-9-.]+",

rule all:
    input:
        expand(config["tool"] + "/input/chr_scale_genomes/{sp}_chr_scale.fasta", sp=config["species"]),
        expand(config["tool"] + "/input/fragmented_genomes/{sp}_fragmented.fasta", sp=config["species"]),
        # expand(config["tool"] + "/output/{sp1}_to_{sp2}_to_{sp}_run/{sp1}_vs_{sp2}_vs_{sp}_intersect.bed", sp=config["species"], 
        # sp1=config["analyses"][this_analysis]["species1"], sp2=config["analyses"][this_analysis]["species2"]),
        file_targets

""""
This will be a part with a bunch of functions used in Breakpointer2
"""
#This function finds the homologous chromosome pairs in two species
def fisher_exact_func(summed_df_sp1, sp1, sp2, table_of_homologous_pairs, n):
    dfs = []
    sp1_chr = summed_df_sp1[sp1+'_chr'].unique()
    sp2_chr = summed_df_sp1[sp2+'_chr'].unique()
    
    # Iterate over each unique value for sp1+'_chr' and sp2+'_chr'
    for chr1 in sp1_chr:
        for chr2 in sp2_chr:
            in1_in2 = summed_df_sp1.loc[(summed_df_sp1[sp1+'_chr'] == chr1) & (summed_df_sp1[sp2+'_chr'] == chr2), 'alen'].sum()
            in1_notin2 = summed_df_sp1.loc[(summed_df_sp1[sp1+'_chr'] == chr1) & (summed_df_sp1[sp2+'_chr'] != chr2), 'alen'].sum()
            notin1_in2 = summed_df_sp1.loc[(summed_df_sp1[sp1+'_chr'] != chr1) & (summed_df_sp1[sp2+'_chr'] == chr2), 'alen'].sum()
            notin1_notin2 = summed_df_sp1.loc[(summed_df_sp1[sp1+'_chr'] != chr1) & (summed_df_sp1[sp2+'_chr'] != chr2), 'alen'].sum()
        
            contingency_matrix = [[int(in1_in2/n), int(in1_notin2/n)], [int(notin1_in2/n), int(notin1_notin2/n)]]
            odds_ratio, p_value = fisher_exact(contingency_matrix, alternative='greater')
            dict1 = {'Chr1': chr1, 'Chr2': chr2, 'Odds ratio': odds_ratio, 'P-value': p_value}
            dfs.append(dict1)

    # Concatenate dictionaries into a DataFrame
    table_of_homologous_pairs = pd.DataFrame(dfs)

    # Calculate Bonferroni corrected p-value
    p_adj = table_of_homologous_pairs['P-value'] * (len(sp1_chr) * len(sp2_chr))
    table_of_homologous_pairs['P-adj'] = p_adj
    #table_of_homologous_pairs.to_csv('FET_test.tsv', index=False, sep='\t')
    return table_of_homologous_pairs

#This function finds the rows that have sudden "jumps" in the alignment and calls breakpoints
# def process_species(pairs, species_list, len_bp, sp1, sp2):
#     the_list = []
    
#     for z in species_list:
#         pairs = pairs.sort_values([z+'_chr', z+'_start', z+'_stop']).reset_index()
#         for j in range(len(pairs) - 1):
#             row = pairs.iloc[j]
#             next_row = pairs.iloc[j+1]
#             first_row_index = pairs.index.tolist()[0]

#             if next_row[sp2+'_start'] < (row[sp2+'_stop'] - len_bp) or \
#                next_row[sp2+'_start'] > (row[sp2+'_stop'] + len_bp) or \
#                next_row[sp1+'_start'] < (row[sp1+'_stop'] - len_bp) or \
#                next_row[sp1+'_start'] > (row[sp1+'_stop'] + len_bp):
                
#                 pairs.at[first_row_index+j, 'breakpoint_in_species'] = row['breakpoint_in_species'] = z
#                 pairs.at[first_row_index+j+1, 'breakpoint_in_species'] = next_row['breakpoint_in_species'] = z
#                 the_list.append(row)
#                 the_list.append(next_row)

#         pairs = pairs.set_index('index')
    
#     return pairs, the_list

#This function subsets the dataframe to include only rows that have breakpoints
# def breakp_range(table_i,sp):
#     sublist_fr_ins_sc = []
#     table_i = table_i.sort_values([sp+'_chr', sp+'_start', sp+'_stop']).reset_index()
#     table_i = table_i.drop(['index'], axis=1)
#     for n in range(len(table_i)-1):
#         if n!=0:
#             if table_i.iloc[n-1]['breakpoints']=='True':
#                 continue
#         row = table_i.iloc[n]
#         next_row = table_i.iloc[n+1]
#         if row['breakpoints'] == 'True':
#             sublist_fr_ins_sc.append(row)
#             sublist_fr_ins_sc.append(next_row)
#     table_bp = pd.concat(sublist_fr_ins_sc, axis=1).transpose().drop_duplicates()
#     table_bp = table_bp.set_index('index').reset_index()
#     one_row_data = pd.DataFrame()
#     table_x = pd.DataFrame()
#     for n in range(0,(len(table_bp)-1),2):
#         row = table_bp.iloc[n]
#         next_row = table_bp.iloc[n+1]
#         one_row_data = pd.concat([row[[sp+'_chr', sp+'_stop']], next_row[[sp+'_start']]], axis=1).transpose()
#         one_row_data.at[n, sp + '_start'] = next_row[sp+'_start']
#         table_x = pd.concat([table_x, one_row_data])
#     return table_x

#This function plots the breakpoints on the chromosomes
def plot_chromosome_comparison(pairs, sp1, sp2, chr1, chr2):
    list_pos_rel_ori = []
    list_neg_rel_ori = []
    fig, ax = plt.subplots(figsize=(10, 10))
    
    for l in range(len(pairs)):
        row = pairs.iloc[l]
        if row['relative_orientation'] == '+':
            list_pos_rel_ori.append(row)
        else:
            list_neg_rel_ori.append(row)
    if len(list_pos_rel_ori)>0:
        positive_table = pd.concat(list_pos_rel_ori, axis=1).transpose()
        x1_values = [positive_table[sp1+'_start'], positive_table[sp1+'_stop']]
        y1_values = [positive_table[sp2+'_start'], positive_table[sp2+'_stop']]
        plt.plot(x1_values, y1_values, 'b', linestyle="solid")
    if len(list_neg_rel_ori)>0:
        negative_table = pd.concat(list_neg_rel_ori, axis=1).transpose()
        x2_values = [negative_table[sp1+'_stop'], negative_table[sp1+'_start']]
        y2_values = [negative_table[sp2+'_start'], negative_table[sp2+'_stop']]
        plt.plot(x2_values, y2_values, 'b', linestyle="solid")

    for m in range(len(pairs)-1):
        row = pairs.iloc[m]
        if row['breakpoints'] == 'True':
            ax.scatter(row[sp1+'_start'], row[sp2+'_stop'], color='r')
        
    ax.ticklabel_format(style='plain', scilimits=(0, 0))
    ax.set_xlabel(sp1 + '_' + chr1)
    ax.set_ylabel(sp2 + '_' + chr2)
    ax.set_title('Chromosome ' + chr1 + ' vs chromosome ' + chr2)
    plt.savefig('comparison_' + chr1 + '_vs_' + chr2 + '.png')

"""
The rules start here
"""
#TODO
#len of the genome
#1%+of the genome to extract (probably a better way)
#or assume that the user knows what they're doing and are giving the chr scale scaffs
rule extract_first_n_sequences:
    input:
        genome = lambda wildcards: config["species"][wildcards.sp]["genome"]
    output:
        output_genome = config["tool"] + "/input/chr_scale_genomes/{sp}_chr_scale.fasta"
    params:
        num_sequences = lambda wildcards: config["species"][wildcards.sp]["num_chr"]
    threads: 1
    run:
        # Check if the input genome file is gzipped
        if input.genome.endswith(".gz"):
            with gzip.open(input.genome, "rt") as handle:  # Open the gzipped file for reading in text mode ("rt")
                records = list(SeqIO.parse(handle, "fasta"))
        else:
            with open(input.genome, "r") as handle:  # Open the file normally
                records = list(SeqIO.parse(handle, "fasta"))

        with open(output.output_genome, "w") as output_handle:
            SeqIO.write(records[:params.num_sequences], output_handle, "fasta")

rule split_query_fasta:
    """
    rule from https://github.com/conchoecia/genome_assembly_pipelines/blob/master/snakefiles/GAP_dgenies_prep
    """
    input:
        assem = config["tool"] + "/input/chr_scale_genomes/{sp}_chr_scale.fasta"
    output:
        assem = config["tool"] + "/input/fragmented_genomes/{sp}_fragmented.fasta"
    threads: 1
    run:
        outhandle = open(output.assem, "w")
        stepsize = 1000000
        for record in SeqIO.parse(input.assem, "fasta"):
            i = 0
            done = False
            record.id = "{}_.-._chunk{}".format(record.id, i)
            while not done:
                SeqIO.write(record[i*stepsize:(i*stepsize)+stepsize], outhandle, 'fasta')
                i += 1
                record.id = record.id.split("_.-._")[0] + "_.-._chunk" + str(i)
                if (i*stepsize)-1 >= len(record.seq):
                    done = True
        outhandle.close()

rule align_genomes:
    input:
        genome1 = config["tool"] + "/input/chr_scale_genomes/{sp2}_chr_scale.fasta",
        genome2 = config["tool"] + "/input/fragmented_genomes/{sp1}_fragmented.fasta"
    output:
        paf = config["tool"] + "/input/paf_frag/{sp1}_to_{sp2}.raw.paf"
    threads: 10
  #  params:
        #seq_arg = config["seq_arg"]
    shell:
        "minimap2 {input.genome1} {input.genome2} > {output.paf}"

rule cleanup_paf:
    """
    rule from https://github.com/conchoecia/genome_assembly_pipelines/blob/master/snakefiles/GAP_dgenies_prep
    """
    input:
        assem  = config["tool"] + "/input/chr_scale_genomes/{sp1}_chr_scale.fasta",
        paf    = config["tool"] + "/input/paf_frag/{sp1}_to_{sp2}.raw.paf"
    output:
        paf    = config["tool"] + "/input/paf_compl/{sp1}_to_{sp2}.paf"
    run:
        # first we get the scaffold size
        scaf_to_size = {}
        for record in SeqIO.parse(input.assem, "fasta"):
            scaf_to_size[record.id] = len(record.seq)

        # now we fix the coordinates of the query sequences
        outhandle = open(output.paf, "w")
        #fields for paf are
        #Col 	Type 	Description
        #1 	string 	Query sequence name
        #2 	int 	Query sequence length
        #3 	int 	Query start (0-based; BED-like; closed)
        #4 	int 	Query end (0-based; BED-like; open)
        #5 	char 	Relative strand: "+" or "-"
        #6 	string 	Target sequence name
        #7 	int 	Target sequence length
        #8 	int 	Target start on original strand (0-based)
        #9 	int 	Target end on original strand (0-based)
        #10 	int 	Number of residue matches
        #11 	int 	Alignment block length
        #12 	int 	Mapping quality (0-255; 255 for missing)
        with open(input.paf, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    fields = line.split("\t")
                    multiplier = int(fields[0].split("_.-._chunk")[1])
                    scaf       = fields[0].split("_.-._chunk")[0]
                    fields[0] = scaf
                    fields[1] = str(scaf_to_size[scaf])
                    fields[2] = str(int(fields[2]) + (1000000 * multiplier))
                    fields[3] = str(int(fields[3]) + (1000000 * multiplier))
                    print("\t".join(fields), file = outhandle)
        outhandle.close()

rule clean_paf:
    input:
        paf = config["tool"] + "/input/paf_compl/{sp1}_to_{sp2}.paf"
    output:
        filtered_tab = config["tool"] + "/input/processed/{sp1}_to_{sp2}/{sp1}_to_{sp2}_filtered.paf",
        FET_test = config["tool"] + "/input/processed/{sp1}_to_{sp2}/{sp1}_to_{sp2}_FET_test.tsv"
    params:
        sp1 = lambda wildcards: wildcards.sp1,
        sp2 = lambda wildcards: wildcards.sp2,
        len_co = config["len_co"],
        q = config["q"],
        FET= config["FET"]
    run: 
        #Input .paf file
        input_paf = input.paf
        #write the name of the query species in the alignment
        sp1 = params.sp1
        #write the name of the reference species in the alignment
        sp2 = params.sp2
        #Write the alignment length cutoff for the table filtering- this step might need some optimization depending on the data
        len_co= int(params.len_co)
        #Write the quality cutoff for the alignments
        q = int(params.q)
        #normalize
        normalization =int(params.FET)
        #write the output file
        filtered_tab = output.filtered_tab
        #write out the FET test results
        FET_test = output.FET_test
        #read in only the necessary columns of a .paf file and give them names
        columns_to_read = list(range(12))  # First 12 columns (0 to 11)
        column_names = [f"{sp1}_chr", f"{sp1}_chr_size", f"{sp1}_start", f"{sp1}_stop","relative_orientation",f"{sp2}_chr", f"{sp2}_chr_size", f"{sp2}_start", f"{sp2}_stop", "nmatch", "alen", "qual_score"]
        sp1_to_sp2 = pd.read_csv(input_paf, delimiter='\t', header=None, usecols=columns_to_read, names=column_names)
        #This one sorts the by the values of 6 columns
        sp1_to_sp2= sp1_to_sp2.sort_values([sp1+'_chr',sp1+'_start',sp1+'_stop'])
        #This filters the dataframe to include only relevant alignments- the length of the alignment must be higher than the values given in len_co and wuality of the alignment must be higher than value q
        sp1_to_sp2 = sp1_to_sp2[(sp1_to_sp2.alen > len_co) & (sp1_to_sp2.qual_score > q)]
        sp1_to_sp2 = sp1_to_sp2.reset_index()
        sp1_to_sp2 = sp1_to_sp2.drop(['index'], axis=1)
        sp1_to_sp2.to_csv(filtered_tab, sep='\t', index=False)
        sp1_to_sp2_hm = sp1_to_sp2[[sp1+'_chr', sp2+'_chr', 'alen']]
        # Group by pairs in the alignment and sum up the alignment lengths
        summed_df_sp1 = sp1_to_sp2_hm.groupby([sp1+'_chr', sp2+'_chr'])['alen'].sum().reset_index()
        #generate a new dataframe to be able to save the Fisher's exact test results
        result_df1 = pd.DataFrame()
        result_df1 = fisher_exact_func(summed_df_sp1, sp1, sp2, result_df1, normalization)
        result_df1 = result_df1[result_df1['P-adj'] < 0.01]
        result_df1.to_csv(FET_test, index=False, sep='\t')

rule breakpointer2:
    input:
        paf = config["tool"] + "/input/processed/{sp1}_to_{sp2}/{sp1}_to_{sp2}_filtered.paf",
        FET_test = config["tool"] + "/input/processed/{sp1}_to_{sp2}/{sp1}_to_{sp2}_FET_test.tsv"
    output:
        breakpoints = config["tool"] + "/output/{sp1}_to_{sp2}_run/{sp1}_vs_{sp2}_breakpoints.tsv",
        bed_file_filt_sp1 = config["tool"] + "/input/raw_bedlike_files/{sp1}_to_{sp2}_run/{sp1}_{sp2}.bed",
        bed_file_filt_sp2 = config["tool"] + "/input/raw_bedlike_files/{sp1}_to_{sp2}_run/{sp2}_{sp1}.bed",
        bed_file_sp1 = config["tool"] + "/input/raw_bedlike_files/{sp1}_to_{sp2}_run/{sp1}_{sp2}_all.bed",
        bed_file_sp2 = config["tool"] + "/input/raw_bedlike_files/{sp1}_to_{sp2}_run/{sp2}_{sp1}_all.bed"
    params:
        sp1 = lambda wildcards: wildcards.sp1,
        sp2 = lambda wildcards: wildcards.sp2,
        len_bp = config["len_bp"],
        bp_co = config["bp_co"]
    run:
        #Write what's the distance between two ends of the aligned fragments to call it a breakpoint
        len_bp = int(params.len_bp)
        bp_co = int(params.bp_co)
        #write the name of the query species in the alignment
        sp1 = params.sp1
        #write the name of the reference species in the alignment
        sp2 = params.sp2
        #read in the FET test results and take the homologous pairs
        sp1_to_sp2 = pd.read_csv(input.paf, delimiter='\t')    
        result_df1 = pd.read_csv(input.FET_test, delimiter='\t')
        list_sp1 = result_df1['Chr1'].tolist()
        list_sp2 = result_df1['Chr2'].tolist()

        chromosome_pairs = dict(zip(list_sp1, list_sp2))

        pairs_dataframe_fin_sp1 = pd.DataFrame()
        pairs_dataframe_fin_sp2 = pd.DataFrame()
        breakpoints_fin_filt_sp1 = pd.DataFrame()
        breakpoints_fin_filt_sp2 = pd.DataFrame()
        breakpoints_fin_sp1 = pd.DataFrame()
        breakpoints_fin_sp2 = pd.DataFrame()
        #The for loop
        species_list = [sp1,sp2]

        for k in chromosome_pairs:
            pairs= format(f"table_{k}")
            globals()[pairs]=sp1_to_sp2.loc[(sp1_to_sp2[sp1+'_chr']==k) & (sp1_to_sp2[sp2+'_chr']==chromosome_pairs[k]),]
            the_list = []
            globals()[pairs]['breakpoint_in_species'] = np.nan
            globals()[pairs]['breakpoint_in_species']=globals()[pairs]['breakpoint_in_species'].astype(str)
            #globals()[pairs], the_list = process_species(globals()[pairs], species_list, len_bp, sp1, sp2)
            for z in species_list:
                globals()[pairs] = globals()[pairs].sort_values([z+'_chr',z+'_start',z+'_stop'])
                globals()[pairs] = globals()[pairs].reset_index()
                for j in range(len(globals()[pairs])-1):
                    row = globals()[pairs].iloc[j]
                    next_row = globals()[pairs].iloc[j+1]
                    first_row_index = globals()[pairs].index.tolist()[0]
                    if next_row[sp2+'_start'] < (row[sp2+'_stop'] - len_bp):
                        #if not any(x.isin(row).all() for x in the_list):
                        globals()[pairs].at[first_row_index+j, 'breakpoint_in_species']=row['breakpoint_in_species']=z
                        globals()[pairs].at[first_row_index+j+1, 'breakpoint_in_species']=next_row['breakpoint_in_species']=z
                        the_list.append(row)
                        the_list.append(next_row)
                    
                    if next_row[sp2+'_start'] > (row[sp2+'_stop'] + len_bp):
                        #if not any(x.isin(row).all() for x in the_list):
                        globals()[pairs].at[first_row_index+j, 'breakpoint_in_species']=row['breakpoint_in_species']=z
                        globals()[pairs].at[first_row_index+j+1, 'breakpoint_in_species']=next_row['breakpoint_in_species']=z
                        the_list.append(row)
                        the_list.append(next_row)
                        
                    if next_row[sp1+'_start'] < (row[sp1+'_stop'] - len_bp):
                        #if not any(x.isin(row).all() for x in the_list):
                        globals()[pairs].at[first_row_index+j, 'breakpoint_in_species']=row['breakpoint_in_species']=z
                        globals()[pairs].at[first_row_index+j+1, 'breakpoint_in_species']=next_row['breakpoint_in_species']=z
                        the_list.append(row)
                        the_list.append(next_row)
                                    
                    if next_row[sp1+'_start'] > (row[sp1+'_stop'] + len_bp):
                        #if not any(x.isin(row).all() for x in the_list):
                        globals()[pairs].at[first_row_index+j, 'breakpoint_in_species']=row['breakpoint_in_species']=z
                        globals()[pairs].at[first_row_index+j+1, 'breakpoint_in_species']=next_row['breakpoint_in_species']=z
                        the_list.append(row)
                        the_list.append(next_row)
                globals()[pairs] = globals()[pairs].set_index('index')
            if len(the_list)>0:
                globals()[pairs] = globals()[pairs].sort_values([sp1+'_chr',sp1+'_start',sp1+'_stop'])
                breakp_df = format(f"concat_list_{k}")
                globals()[breakp_df]= pd.concat(the_list, axis=1)
                globals()[breakp_df] = globals()[breakp_df].transpose()#.drop_duplicates()
                globals()[breakp_df] = globals()[breakp_df].set_index('index')
                globals()[breakp_df] = globals()[breakp_df][~globals()[breakp_df].index.duplicated(keep='first')].sort_values([sp1+'_chr',sp1+'_start',sp1+'_stop'])
                globals()[pairs] = globals()[pairs].drop(['breakpoint_in_species'], axis=1)
                globals()[breakp_df]['breakpoints']= 'True'
                subset = globals()[breakp_df][['breakpoints','breakpoint_in_species']]
                globals()[pairs] = globals()[pairs].join(subset, how='left')
                #now plot all the datapoints in the dataframe to visualized called breakpoints
                #plot_chromosome_comparison(globals()[pairs], sp1, sp2, k, chromosome_pairs[k])
                bed_file = format(f"bed_file_{k}")
                bed_file_all = format(f"bed_file_all_{k}")
                #globals()[bed_file] = breakp_range(globals()[pairs],i)
                table_i = pd.DataFrame()
                table_all = pd.DataFrame()
                for i in species_list:
                    table_i = globals()[pairs].sort_values([i+'_chr', i+'_start', i+'_stop']).reset_index()
                    table_i = table_i.drop(['index'], axis=1)
                    sublist_fr_ins_sc = []
                    for n in range(len(table_i)-1):
                        if n != 0 and table_i.iloc[n-1]['breakpoints'] == 'True':
                            continue
                        row = table_i.iloc[n]
                        next_row = table_i.iloc[n+1]
                        if row['breakpoints'] == 'True':
                            sublist_fr_ins_sc.append(row)
                            sublist_fr_ins_sc.append(next_row)
                    table_bp = pd.DataFrame(sublist_fr_ins_sc)
                    table_bp = table_bp.sort_values([i+'_chr', i+'_start', i+'_stop']).reset_index()
                    table_bp = table_bp.drop(['index'], axis=1)
                    one_row_data = pd.DataFrame()   
                    table_x = pd.DataFrame()
                    table_all = pd.DataFrame()
                    for n in range(0,(len(table_bp)-1),2):
                        row = table_bp.iloc[n]
                        next_row = table_bp.iloc[n+1]
                        one_row_data = pd.concat([row[[i+'_chr', i+'_stop']], next_row[[i+'_start']]], axis=1).transpose()
                        one_row_data.at[n, i + '_start'] = next_row[i+'_start']
                        table_all = pd.concat([table_all, one_row_data])
                        if (one_row_data.at[n,i+'_start']-one_row_data.at[n,i+'_stop'])<bp_co:
                            table_x = pd.concat([table_x, one_row_data])
                        else:
                            continue
                    table_x = table_x.dropna()
                    table_all = table_all.dropna()
                    globals()[bed_file] = table_x
                    globals()[bed_file_all] = table_all
                    if i==sp1:
                        pairs_dataframe_fin_sp1 = pd.concat([pairs_dataframe_fin_sp1, globals()[pairs]])
                        breakpoints_fin_filt_sp1 = pd.concat([breakpoints_fin_filt_sp1, globals()[bed_file]]).dropna()
                        breakpoints_fin_sp1 = pd.concat([breakpoints_fin_sp1, globals()[bed_file_all]]).dropna()
                    if i==sp2:
                        pairs_dataframe_fin_sp2 = pd.concat([pairs_dataframe_fin_sp2, globals()[pairs]])
                        breakpoints_fin_filt_sp2 = pd.concat([breakpoints_fin_filt_sp2, globals()[bed_file]]).dropna()
                        breakpoints_fin_sp2 = pd.concat([breakpoints_fin_sp2, globals()[bed_file_all]]).dropna()
                    else:
                        continue
        pairs_dataframe_fin_sp1.to_csv(output.breakpoints, sep='\t', index=False)
        breakpoints_fin_sp1.to_csv(output.bed_file_sp1, sep='\t', header=False, index=False)
        breakpoints_fin_sp2.to_csv(output.bed_file_sp2, sep='\t', header=False, index=False)
        breakpoints_fin_filt_sp1.to_csv(output.bed_file_filt_sp1, sep='\t', header=False, index=False)
        breakpoints_fin_filt_sp2.to_csv(output.bed_file_filt_sp2, sep='\t', header=False, index=False)