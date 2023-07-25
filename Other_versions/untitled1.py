#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 11:11:58 2023

@author: dalilad
"""

# -*- coding: utf-8 -*-
"""
To be written
use this code to get the scaffold names for the assemblies: grep '^>' fin_oct_vulg_30_scaff.fa | sed -E 's/^[[:space:]]*[^[:space:]]+/&,/' | sed 's/>//g' | sed 's/[^[:space:],]\+/"&"/g' | tr '\n' ' '
"""
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import statistics
#from mlxtend.evaluate import permutation_test
#from scipy.stats import ttest_ind

#write the name of the query species in the alignment
sp1 = 'Obi'
#write the name of the reference species in the alignment
sp2 = 'Ovu'
#Write the alignment length cutoff for the table filtering
len_co= 100000
#Write the quality cutoff for the alignments
q = 30
#Write what's the distance between two ends of the aligned fragments to call it a breakpoint
len_bp = 5000000
#set the number of permutations
num_rounds=1000
# define a custom function to perform the permutation test


#EBR=evolutionary breakpoints
#make sure to know where you're working, can be modified later as an input
os.chdir(r'/Users/dalilad/Desktop/The_octopus_code/breakpointer2-main')
#os.chdir(r'/Users/dalilad/Desktop/The_octopus_code')
#print (os.getcwd())

#to disable the warning when adding the breakpoints column to the original dataframe
pd.options.mode.chained_assignment = None

#The .paf file part, aka the analysis of the breakpoints in genome-genome alignments

#This part reads in the .paf file as Pandas DataFrame
sp1_to_sp2= pd.read_csv('bimaculoides_to_vulgaris.paf', delimiter='\t', header=None, names= [sp1+'_chr',sp1+'_chr_size',sp1+'_start',sp1+'_stop','relative_orientation',sp2+'_chr',sp2+'_chr_size',sp2+'_start',sp2+'_stop','nmatch','alen','qual_score','col12','col13','col14','col15','col16','col17','col18'])

#The .paf file has some extra columns which are unnecessary, so this line removes them
sp1_to_sp2= sp1_to_sp2.drop(['col12','col13','col14','col15','col16','col17','col18'], axis =1)

#This one sorts the by the values of 6 columns
sp1_to_sp2= sp1_to_sp2.sort_values([sp1+'_chr',sp1+'_start',sp1+'_stop'])

#This filters the dataframe to include only relevant alignments- the length of the alignment must be higher than the values given in len_co and wuality of the alignment must be higher than value q
sp1_to_sp2 = sp1_to_sp2[(sp1_to_sp2.alen > len_co) & (sp1_to_sp2.qual_score > q)]

sp1_to_sp2 = sp1_to_sp2.reset_index()

sp1_to_sp2 = sp1_to_sp2.drop(['index'], axis=1)
#Now I want to make a dictionary of homologous chromosomes
#since the vulg has been named according to it's homology with bimac, i can just use that, however, it would be better to have it automated and usable for wider type of datasets
#got an idea from here: https://www.geeksforgeeks.org/python-convert-two-lists-into-a-dictionary/
dict_key=["CM046601.1", "CM046602.1", "CM046603.1", "CM046604.1", "CM046605.1", "CM046606.1", "CM046607.1", "CM046608.1", "CM046609.1", "CM046610.1", "CM046611.1", "CM046612.1", "CM046613.1", "CM046614.1", "CM046615.1", "CM046616.1", "CM046617.1", "CM046618.1", "CM046619.1", "CM046620.1"]#, "CM046621.1", "CM046622.1", "CM046623.1", "CM046624.1", "CM046625.1", "CM046626.1", "CM046627.1", "CM046628.1", "CM046629.1", "CM046630.1"]
dict_value=["Ovu01", "Ovu02", "Ovu03", "Ovu04", "Ovu05", "Ovu06", "Ovu07", "Ovu08", "Ovu09", "Ovu10", "Ovu11", "Ovu12", "Ovu13", "Ovu14", "Ovu15", "Ovu16", "Ovu17", "Ovu18", "Ovu19", "Ovu20"]#, "Ovu21", "Ovu22", "Ovu23", "Ovu24", "Ovu25", "Ovu26", "Ovu27", "Ovu28", "Ovu29", "Ovu30"]
chromosome_pairs = dict(zip(dict_key, dict_value))

#.gff or .bed file part
# 100kb windows for .bed
#try to plot 5k thing on the HiGlass for example
bimac_insul_sc_500kb = pd.read_csv('oct_bim_0_normalized_balanced_diff_wind_500kb.bed', delimiter='\t', header= None, names=[sp1+'_chr', sp1+'_start',sp1+'_stop', 'irrel1', 'insul_score','irrel2'])
bimac_insul_sc_500kb = bimac_insul_sc_500kb.dropna().drop(['irrel1','irrel2'], axis=1)

vulg_insul_sc_500kb = pd.read_csv('fin_oct_vulg_500kb.bed', delimiter='\t', header= None, names=[sp2+'_chr', sp2+'_start',sp2+'_stop', 'irrel1', 'insul_score','irrel2'])
vulg_insul_sc_500kb = vulg_insul_sc_500kb.dropna().drop(['irrel1','irrel2'], axis=1)

#This function takes the insulation score table and species according to which the table was sorted and adds information on insulation scores of the breakpoints
def insul_scores_sp(breakp_df, sp_insul_k,i):
    breakp_df = breakp_df.reset_index()
    for n in range(0,(len(breakp_df)-1),2):
        row = breakp_df.iloc[n]
        next_row = breakp_df.iloc[n+1]
        midpoint_diff_i = abs(row[i+'_stop']-next_row[i+'_start'])/2
        midpoint_value = row[i+'_stop']+midpoint_diff_i
        breakp_df.at[n, 'midpoint_'+i] = row['midpoint_'+i] = midpoint_value
        for o in range(len(sp_insul_k)):
            row_ins_i = sp_insul_k.iloc[o]
            if row_ins_i[i+'_start'] <= midpoint_value <= row_ins_i[i+'_stop']:
                breakp_df.at[n, 'insul_score_'+i] = row['insul_score_'+i] = row_ins_i['insul_score']
    if len(breakp_df) % 2 == 1:
        last_row = breakp_df.iloc[-1]
        midpoint_diff_i = abs(next_row[i+'_stop']-last_row[i+'_start'])
        midpoint_value = next_row[i+'_stop'] + midpoint_diff_i
        breakp_df.at[len(breakp_df)-1, 'midpoint_'+i] = last_row['midpoint_'+i] = midpoint_value
        for o in range(len(sp_insul_k)):
            row_ins_i = sp_insul_k.iloc[o]
            if row_ins_i[i+'_start'] <= midpoint_value <= row_ins_i[i+'_stop']:
                breakp_df.at[len(breakp_df)-1, 'insul_score_'+i] = row['insul_score_'+i] = row_ins_i['insul_score']
    breakp_df = breakp_df.set_index('index')
    return breakp_df
"""
because of the gaps in the breakp dataframe, it is not possible to calculate potential midpoints for the opposite species, therefore, a different 
approach is needed
table_i=the merged table for a chromosome, including both breakpoints and non-breakpoints
i_insul_sc= insulation score table for the opposite species
sp = the opposite species
"""
def insul_score_oppo_sp(table_i, i_insul_k,sp):
    sublist_fr_ins_sc = []
    table_i = table_i.sort_values([sp+'_chr', sp+'_start', sp+'_stop']).reset_index()
    for n in range(len(table_i)-1):
        if n!=0:
            if table_i.iloc[n-1]['breakpoints']=='True':
                continue
        row = table_i.iloc[n]
        next_row = table_i.iloc[n+1]
        if row['breakpoints'] == 'True':
            sublist_fr_ins_sc.append(row)
            sublist_fr_ins_sc.append(next_row)
    table_opp_sp = pd.concat(sublist_fr_ins_sc, axis=1).transpose().drop_duplicates()
    table_opp_sp = table_opp_sp.set_index('index').reset_index()
    for n in range(0,(len(table_opp_sp)-1),2):
        row = table_opp_sp.iloc[n]
        next_row = table_opp_sp.iloc[n+1]
        midpoint_diff_i = abs(row[sp+'_stop']-next_row[sp+'_start'])/2
        midpoint_value = row[sp+'_stop']+midpoint_diff_i
        table_opp_sp.at[n, 'midpoint_'+sp] = row['midpoint_'+sp] = midpoint_value
        for o in range(len(i_insul_k)):
            row_ins_i = i_insul_k.iloc[o]
            if row_ins_i[sp+'_start'] <= midpoint_value <= row_ins_i[sp+'_stop']:
                table_opp_sp.at[n, 'insul_score_'+sp] = row['insul_score_'+sp] = row_ins_i['insul_score']
    if len(table_opp_sp) % 2 == 1:
        last_row = table_opp_sp.iloc[-1]
        midpoint_diff_i = abs(next_row[sp+'_stop']-last_row[sp+'_start'])
        midpoint_value = next_row[sp+'_stop'] + midpoint_diff_i
        table_opp_sp.at[len(table_opp_sp)-1, 'midpoint_'+sp] = last_row['midpoint_'+sp] = midpoint_value
        for o in range(len(i_insul_k)):
            row_ins_i = i_insul_k.iloc[o]
            if row_ins_i[sp+'_start'] <= midpoint_value <= row_ins_i[sp+'_stop']:
                table_opp_sp.at[len(table_opp_sp)-1, 'insul_score_'+sp] = row['insul_score_'+sp] = row_ins_i['insul_score']
    table_opp_sp = table_opp_sp.set_index('index')
    table_i = table_i.drop(['midpoint_'+sp, 'insul_score_'+sp], axis=1).set_index('index')
    table_opp_sp = table_opp_sp[['midpoint_'+sp, 'insul_score_'+sp]]
    table_i = table_i.join(table_opp_sp, how='left')
    return table_i


def flatten(l):
    return [item for sublist in l for item in sublist]

def permut_test_histogram (insul_sc_sp_bp, insul_sc_sp, sp):
    median_sp_bp = statistics.median(insul_sc_sp_bp)
    median_list = []
    permut_samples = []
    for p in range(num_rounds):
        permut_set = np.random.choice(insul_sc_sp, size= len(insul_sc_sp_bp),replace=False)
        permut_median = statistics.median(permut_set)
        median_list.append(permut_median)
        permut_samples.append(permut_set)
    #write a line here to save the permut samples    
    permut_samples = flatten(permut_samples)
    fig, ax = plt.subplots(figsize=(10, 10))
    # plot the data points as dots
    ax.plot(insul_sc_sp_bp, np.zeros_like(insul_sc_sp_bp), 'ro', alpha=0.5, label='Data Points')
    # plot the median of the data as a vertical line
    ax.axvline(median_sp_bp, color='r', label='Data Median')
    # plot the medians of the permutations as a histogram
    ax.hist(median_list, bins=10, alpha=0.5, label='Permutation Medians')
    if sp == sp1 and len(insul_sc_sp_bp)<80:
        plt.gca().set(title='Histogram for chr ' + k, ylabel='Frequency')
    if sp == sp2 and len(insul_sc_sp_bp)<80:
        plt.gca().set(title='Histogram for chr ' + chromosome_pairs[k], ylabel='Frequency')
    if len(insul_sc_sp_bp)>=80:
        plt.gca().set(title='Histogram for overall chromosomes of ' + sp, ylabel='Frequency')
    plt.legend()
    bigger_than= np.sum(median_list>median_sp_bp)
    smaller_than = np.sum(median_list<median_sp_bp)
    fdr = smaller_than/num_rounds
    if fdr == 0:
        corr= 1/num_rounds
        fdr= format(f"<{corr}")
    print(fdr)
    index = list(range(5))
    if sp == sp1:
        stats_df = {
            'Dataset' : df_name,
            'Chromosome' : k,
            'Insulation_breakpoint_median' : median_sp_bp,
            'Flase_detection_rate' : fdr,
            'Number_of_smaller_permutation_medians' : smaller_than,
            'Number_of_bigger_permutation_medians' : bigger_than
            }
    else:
        stats_df = {
            'Dataset' : df_name,
            'Chromosome' : chromosome_pairs[k],
            'Insulation_breakpoint_median' : median_sp_bp,
            'Flase_detection_rate' : fdr,
            'Number_of_smaller_permutation_medians' : smaller_than,
            'Number_of_bigger_permutation_medians' : bigger_than
            }
    stats_df = pd.DataFrame(data = stats_df, index = index)
    stats_df = stats_df.drop_duplicates()
    return fig, ax, stats_df

pairs_dataframe_fin_sp1 = pd.DataFrame()
pairs_dataframe_fin_sp2 = pd.DataFrame()
breakpoint_df_sp1 = pd.DataFrame()
breakpoint_df_sp2= pd.DataFrame()
stats_df_sp2 = pd.DataFrame()
stats_df_sp1 = pd.DataFrame()
overall_res = pd.DataFrame()
#The for loop
species_list = [sp1,sp2]
#This one sorts the by the values of 6 columns
for i in species_list:
    df_name = format(f"sorted_by_{i}")
    globals()[df_name] = sp1_to_sp2.sort_values([i+'_chr',i+'_start',i+'_stop'])
    globals()[df_name] = globals()[df_name].reset_index().drop(['index'], axis=1)
    final_insul_list_sp1 = []
    final_insul_list_sp2 = []
    insul_table_sp2 = []
    insul_table_sp1 = []
    for k in chromosome_pairs:
        pairs= format(f"table_{i}")
        globals()[pairs]=globals()[df_name].loc[(globals()[df_name][sp1+'_chr']==k) & (globals()[df_name][sp2+'_chr']==chromosome_pairs[k]),]
        sp1_insul_k = bimac_insul_sc_500kb.loc[(bimac_insul_sc_500kb[sp1+'_chr'] == k),]
        sp2_insul_k = vulg_insul_sc_500kb.loc[(vulg_insul_sc_500kb[sp2+'_chr']== chromosome_pairs[k]),]
        #globals()[pairs] = globals()[pairs].reset_index()
        #I open this list as it has been shown as the easiest-least-crashiest way to do what I want to do
        the_list = []
    
        #bed file starts and stops: the start is the least num and the stop is higher (regardless of the strand)
        #write an if statement and plot them vice versa for - direction
        globals()[pairs]['breakpoint_in_species'] = np.nan
        globals()[pairs]['breakpoint_in_species']=globals()[pairs]['breakpoint_in_species'].astype(str)
        #bed file starts and stops: the start is the least num and the stop is higher (regardless of the strand)
        #write an if statement and plot them vice versa for - direction

        for j in range(len(globals()[pairs])-1):
            row = globals()[pairs].iloc[j]
            next_row = globals()[pairs].iloc[j+1]
            first_row_index = globals()[pairs].index.tolist()[0]
            if i==sp1:    
                if next_row[sp2+'_start'] < (row[sp2+'_stop'] - len_bp):
                    #if not any(x.isin(row).all() for x in the_list):
                    globals()[pairs].at[first_row_index+j, 'breakpoint_in_species']=row['breakpoint_in_species']=i
                    globals()[pairs].at[first_row_index+j+1, 'breakpoint_in_species']=next_row['breakpoint_in_species']=i
                    the_list.append(row)
                    the_list.append(next_row)
                
                if next_row[sp2+'_start'] > (row[sp2+'_stop'] + len_bp):
                    #if not any(x.isin(row).all() for x in the_list):
                    globals()[pairs].at[first_row_index+j, 'breakpoint_in_species']=row['breakpoint_in_species']=i
                    globals()[pairs].at[first_row_index+j+1, 'breakpoint_in_species']=next_row['breakpoint_in_species']=i
                    the_list.append(row)
                    the_list.append(next_row)
            else:    
                if next_row[sp1+'_start'] < (row[sp1+'_stop'] - len_bp):
                    #if not any(x.isin(row).all() for x in the_list):
                    globals()[pairs].at[first_row_index+j, 'breakpoint_in_species']=row['breakpoint_in_species']=i
                    globals()[pairs].at[first_row_index+j+1, 'breakpoint_in_species']=next_row['breakpoint_in_species']=i
                    the_list.append(row)
                    the_list.append(next_row)
                
                if next_row[sp1+'_start'] > (row[sp1+'_stop'] + len_bp):
                    #if not any(x.isin(row).all() for x in the_list):
                    globals()[pairs].at[first_row_index+j, 'breakpoint_in_species']=row['breakpoint_in_species']=i
                    globals()[pairs].at[first_row_index+j+1, 'breakpoint_in_species']=next_row['breakpoint_in_species']=i
                    the_list.append(row)
                    the_list.append(next_row)
                    
        #Seems like breakpoint can be detected in both, should be aware of that
        if len(the_list)>0:
            breakp_df = format(f"concat_list_{i}")
            #the_list = pd.Series(the_list).drop_duplicates(keep='last').tolist()
            globals()[breakp_df]= pd.concat(the_list, axis=1)
            
            globals()[breakp_df] = globals()[breakp_df].transpose().drop_duplicates(keep='last')
            globals()[pairs] = globals()[pairs].drop(['breakpoint_in_species'], axis=1)
            globals()[breakp_df]['breakpoints']= 'True'
            
            #globals()[breakp_df] = globals()[breakp_df].reset_index()
            globals()[breakp_df][['midpoint_'+sp1,'midpoint_'+sp2, 'insul_score_'+sp1, 'insul_score_'+sp2]] = np.nan
            
            #result_table = format(f"checkpoint_insul_sc_{i}_{k}")
            result_table = format(f'{i}_midpoint')
            
            if i == sp1:
                globals()[result_table] = insul_scores_sp(globals()[breakp_df], sp1_insul_k,i)
            if i == sp2:
                globals()[result_table] = insul_scores_sp(globals()[breakp_df], sp2_insul_k,i)
            globals()[breakp_df]=globals()[result_table]
            subset = globals()[breakp_df][['breakpoints','breakpoint_in_species', 'midpoint_'+sp1,'midpoint_'+sp2, 'insul_score_'+sp1, 'insul_score_'+sp2]]
            globals()[pairs] = globals()[pairs].join(subset, how='left')
           
            if i ==sp1:
               globals()[pairs] = insul_score_oppo_sp(globals()[pairs], sp2_insul_k,sp2)
            else:
               globals()[pairs] = insul_score_oppo_sp(globals()[pairs], sp1_insul_k,sp1)
                
            globals()[pairs] = globals()[pairs].sort_values([i+'_chr',i+'_start',i+'_stop'])
            table = globals()[pairs]
            pairs= format(f"table_{i}_{k}")
            globals()[pairs] = table
            insul_sc_sp1_bp = globals()[pairs]['insul_score_'+sp1].to_numpy()
            insul_sc_sp1_bp = [x for x in insul_sc_sp1_bp if not np.isnan(x)]
            insul_sc_sp2_bp = globals()[pairs]['insul_score_'+sp2].to_numpy()
            insul_sc_sp2_bp = [x for x in insul_sc_sp2_bp if not np.isnan(x)]
            sp1_insul_k = sp1_insul_k['insul_score'].tolist()
            sp2_insul_k = sp2_insul_k['insul_score'].tolist()
            
            fig, ax, stats_df = permut_test_histogram(insul_sc_sp2_bp, sp2_insul_k, sp2)
            stats_df_sp2 = pd.concat([stats_df,stats_df_sp2]).sort_values('Chromosome')
            final_insul_list_sp2.append(insul_sc_sp2_bp)
            insul_table_sp2.append(sp2_insul_k)
            
            fig, ax, stats_df = permut_test_histogram(insul_sc_sp1_bp, sp1_insul_k, sp1)
            stats_df_sp1 = pd.concat([stats_df,stats_df_sp1]).sort_values('Chromosome')
            final_insul_list_sp1.append(insul_sc_sp1_bp)
            insul_table_sp1.append(sp1_insul_k)
                        
            list_pos_rel_ori=[]
            list_neg_rel_ori=[]
            #if i == sp1:
            fig, ax = plt.subplots(figsize=(10, 10))
            for l in range(len(globals()[pairs])):
                row = globals()[pairs].iloc[l]
                if row['relative_orientation'] == '+':
                    list_pos_rel_ori.append(row)
                else:
                    list_neg_rel_ori.append(row)
            positive_table=pd.concat(list_pos_rel_ori, axis=1)
            positive_table= positive_table.transpose()
            x1_values=[positive_table[sp1+'_start'], positive_table[sp1+'_stop']]
            y1_values=[positive_table[sp2+'_start'], positive_table[sp2+'_stop']]
            negative_table=pd.concat(list_neg_rel_ori, axis=1)
            negative_table= negative_table.transpose()
            x2_values=[negative_table[sp1+'_stop'], negative_table[sp1+'_start']]
            y2_values=[negative_table[sp2+'_start'], negative_table[sp2+'_stop']]
            plt.plot(x1_values, y1_values, 'b', linestyle="solid")
            plt.plot(x2_values, y2_values, 'b', linestyle="solid")
            # # connect adjacent points with lines
            # for i in range(len(globals()[pairs])-1):
            #     row=globals()[pairs].iloc[i]
            #     ax.plot(row[sp1+'_start'], row[sp2+'_stop'], 'k-', alpha=0.5)
            breakp=[]

            for m in range(len(globals()[pairs])-1):
                row = globals()[pairs].iloc[m] 
                #z=(row[sp1+'_start'])/(row[sp1+'_chr_size'])
                if row['breakpoints']=='True':
                    ax.scatter(row[sp1+'_start'], row[sp2+'_stop'], color='r')
                    breakp.append('yes')
            
            ax.ticklabel_format(style='plain', scilimits=(0,0))
            ax.set_xlabel('Octopus bimaculoides chromosome '+k)
            ax.set_ylabel('Octopus vulgaris chromosome '+chromosome_pairs[k])
            ax.set_title('Chromosome '+k+' vs chromosome '+chromosome_pairs[k])
            
            
        if i==sp1:
            pairs_dataframe_fin_sp1 = pd.concat([pairs_dataframe_fin_sp1, globals()[pairs]])
            breakpoint_df_sp1 = pd.concat([breakpoint_df_sp1, globals()[breakp_df]])
        if i==sp2:
            pairs_dataframe_fin_sp2 = pd.concat([pairs_dataframe_fin_sp2, globals()[pairs]])
            breakpoint_df_sp2 = pd.concat([breakpoint_df_sp2, globals()[breakp_df]])
        else:
         continue
    k = i
    chromosome_pairs[k] = i
    final_insul_list_sp1= flatten(final_insul_list_sp1)
    insul_table_sp1 = flatten(insul_table_sp1)
    fig, ax, stats_df = permut_test_histogram(final_insul_list_sp1, insul_table_sp1, sp1)
    overall_res = pd.concat([stats_df, overall_res])
    final_insul_list_sp2= flatten(final_insul_list_sp2)
    insul_table_sp2 = flatten(insul_table_sp2)
    fig, ax, stats_df = permut_test_histogram(final_insul_list_sp2, insul_table_sp2, sp2)
    overall_res = pd.concat([stats_df, overall_res])




