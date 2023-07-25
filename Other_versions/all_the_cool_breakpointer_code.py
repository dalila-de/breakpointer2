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
from scipy.stats import ttest_ind

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
dict_key=["CM046601.1", "CM046602.1", "CM046603.1", "CM046604.1", "CM046605.1", "CM046606.1", "CM046607.1", "CM046608.1", "CM046609.1", "CM046610.1", "CM046611.1", "CM046612.1", "CM046613.1", "CM046614.1", "CM046615.1", "CM046616.1", "CM046617.1", "CM046618.1", "CM046619.1", "CM046620.1", "CM046621.1", "CM046622.1", "CM046623.1", "CM046624.1", "CM046625.1", "CM046626.1", "CM046627.1", "CM046628.1", "CM046629.1", "CM046630.1"]
dict_value=["Ovu01", "Ovu02", "Ovu03", "Ovu04", "Ovu05", "Ovu06", "Ovu07", "Ovu08", "Ovu09", "Ovu10", "Ovu11", "Ovu12", "Ovu13", "Ovu14", "Ovu15", "Ovu16", "Ovu17", "Ovu18", "Ovu19", "Ovu20", "Ovu21", "Ovu22", "Ovu23", "Ovu24", "Ovu25", "Ovu26", "Ovu27", "Ovu28", "Ovu29", "Ovu30"]
chromosome_pairs = dict(zip(dict_key, dict_value))

#.gff or .bed file part
# 100kb windows for .bed
#try to plot 5k thing on the HiGlass for example
bimac_insul_sc_500kb = pd.read_csv('oct_bim_0_normalized_balanced_diff_wind_500kb.bed', delimiter='\t', header= None, names=[sp1+'_chr', sp1+'_start',sp1+'_stop', 'irrel1', 'insul_score','irrel2'])
bimac_insul_sc_500kb = bimac_insul_sc_500kb.dropna().drop(['irrel1','irrel2'], axis=1)

vulg_insul_sc_500kb = pd.read_csv('fin_oct_vulg_500kb.bed', delimiter='\t', header= None, names=[sp2+'_chr', sp2+'_start',sp2+'_stop', 'irrel1', 'insul_score','irrel2'])
vulg_insul_sc_500kb = vulg_insul_sc_500kb.dropna().drop(['irrel1','irrel2'], axis=1)

pairs_dataframe_fin_sp1 = pd.DataFrame()
pairs_dataframe_fin_sp2 = pd.DataFrame()
breakpoint_df_sp1 = pd.DataFrame()
breakpoint_df_sp2= pd.DataFrame()
stat_test_df_sp1 = pd.DataFrame()
stat_test_df_sp2 = pd.DataFrame()
final_insul_list_sp1 = []
final_insul_list_sp2 = []
insul_table = []
#The for loop
species_list = [sp1,sp2]
#This one sorts the by the values of 6 columns
for i in species_list:
    df_name = format(f"sorted_by_{i}")
    globals()[df_name] = sp1_to_sp2.sort_values([i+'_chr',i+'_start',i+'_stop'])
    globals()[df_name] = globals()[df_name].reset_index().drop(['index'], axis=1)
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
            
            globals()[breakp_df] = globals()[breakp_df].reset_index()
            globals()[breakp_df][['midpoint_'+sp1,'midpoint_'+sp2, 'insul_score_'+sp1, 'insul_score_'+sp2]] = np.nan
            
            if i ==sp1:
                for n in range(0,(len(globals()[breakp_df])-1),2):
                    row= globals()[breakp_df].iloc[n]
                    next_row = globals()[breakp_df].iloc[n+1]
                    midpoint_diff_sp1=abs(row[sp1+'_stop']-next_row[sp1+'_start'])/2
                    midpoint_value = row[sp1+'_stop']+midpoint_diff_sp1
                    globals()[breakp_df].at[n, 'midpoint_'+sp1]=row['midpoint_'+sp1]=midpoint_value
                    for o in range(len(sp1_insul_k)):
                        row_ins_sp1= sp1_insul_k.iloc[o]
                        if row_ins_sp1[sp1+'_start'] <= midpoint_value <= row_ins_sp1[sp1+'_stop']:
                            globals()[breakp_df].at[n,'insul_score_'+sp1]=row['insul_score_'+sp1]=row_ins_sp1['insul_score']
                if len(globals()[breakp_df])%2==1:
                    last_row = globals()[breakp_df].iloc[-1]
                    midpoint_diff_sp1= abs(next_row[sp1+'_stop']-last_row[sp1+'_start'])
                    midpoint_value = next_row[sp1+'_stop']+midpoint_diff_sp1
                    globals()[breakp_df].at[len(globals()[breakp_df]), 'midpoint_'+sp1]=last_row['midpoint_'+sp1]=midpoint_value
                    for o in range(len(sp1_insul_k)):
                        row_ins_sp1= sp1_insul_k.iloc[o]
                        if row_ins_sp1[sp1+'_start'] <= midpoint_value <= row_ins_sp1[sp1+'_stop']:
                            globals()[breakp_df].at[n,'insul_score_'+sp1]=row['insul_score_'+sp1]=row_ins_sp1['insul_score']
                globals()[breakp_df] = globals()[breakp_df].set_index('index')
                globals()[breakp_df] = globals()[breakp_df].sort_values([sp2+'_chr',sp2+'_start',sp2+'_stop']).reset_index()
                for m in range(0,(len(globals()[breakp_df])-1),2):
                    row= globals()[breakp_df].iloc[m]
                    next_row = globals()[breakp_df].iloc[m+1]
                    #row_index = globals()[breakp_df].index.tolist()[m]
                    midpoint_diff_sp2=abs(row[sp2+'_stop']-next_row[sp2+'_start'])/2
                    midpoint_value = row[sp2+'_stop']+midpoint_diff_sp2
                    globals()[breakp_df].at[m, 'midpoint_'+sp2]=row['midpoint_'+sp2]=midpoint_value
                    for o in range(len(sp2_insul_k)):
                        row_ins_sp2= sp2_insul_k.iloc[o]
                        if row_ins_sp2[sp2+'_start'] <= midpoint_value <= row_ins_sp2[sp2+'_stop']:
                            globals()[breakp_df].at[m,'insul_score_'+sp2]=row['insul_score_'+sp2]=row_ins_sp2['insul_score']
                if len(globals()[breakp_df])%2==1:
                    last_row = globals()[breakp_df].iloc[-1]
                    midpoint_diff_sp2= abs(next_row[sp2+'_stop']-last_row[sp2+'_start'])
                    midpoint_value = next_row[sp2+'_stop']+midpoint_diff_sp2
                    globals()[breakp_df].at[len(globals()[breakp_df]), 'midpoint_'+sp2]=last_row['midpoint_'+sp2]=midpoint_value
                    for o in range(len(sp2_insul_k)):
                        row_ins_sp2= sp2_insul_k.iloc[o]
                        if row_ins_sp2[sp2+'_start'] <= midpoint_value <= row_ins_sp2[sp2+'_stop']:
                            globals()[breakp_df].at[len(globals()[breakp_df]),'insul_score_'+sp2]=row['insul_score_'+sp2]=row_ins_sp2['insul_score']
                globals()[breakp_df] = globals()[breakp_df].sort_values([sp1+'_chr',sp1+'_start',sp1+'_stop']).set_index('index')
            else:
                for n in range(0,(len(globals()[breakp_df])-1),2):
                    row= globals()[breakp_df].iloc[n]
                    next_row = globals()[breakp_df].iloc[n+1]
                    midpoint_diff_sp2=abs(row[sp2+'_stop']-next_row[sp2+'_start'])/2
                    midpoint_value = row[sp2+'_stop']+midpoint_diff_sp2
                    globals()[breakp_df].at[n, 'midpoint_'+sp2]=row['midpoint_'+sp2]=midpoint_value
                    for o in range(len(sp2_insul_k)):
                        row_ins_sp2= sp2_insul_k.iloc[o]
                        if row_ins_sp2[sp2+'_start'] <= midpoint_value <= row_ins_sp2[sp2+'_stop']:
                            globals()[breakp_df].at[n,'insul_score_'+sp2]=row['insul_score_'+sp2]=row_ins_sp2['insul_score']
                if len(globals()[breakp_df])%2==1:
                    last_row = globals()[breakp_df].iloc[-1]
                    midpoint_diff_sp2= abs(next_row[sp2+'_stop']-last_row[sp2+'_start'])
                    midpoint_value = next_row[sp2+'_stop']+midpoint_diff_sp2
                    globals()[breakp_df].at[len(globals()[breakp_df]), 'midpoint_'+sp2]=row['midpoint_'+sp2]=midpoint_value
                    for o in range(len(sp2_insul_k)):
                        row_ins_sp2= sp2_insul_k.iloc[o]
                        if row_ins_sp2[sp2+'_start'] <= midpoint_value <= row_ins_sp2[sp2+'_stop']:
                            globals()[breakp_df].at[len(globals()[breakp_df]),'insul_score_'+sp2]=row['insul_score_'+sp2]=row_ins_sp2['insul_score']

                globals()[breakp_df] = globals()[breakp_df].set_index('index')
                globals()[breakp_df] = globals()[breakp_df].sort_values([sp1+'_chr',sp1+'_start',sp1+'_stop']).reset_index()
                for m in range(0,(len(globals()[breakp_df])-1),2):
                    row= globals()[breakp_df].iloc[m]
                    next_row = globals()[breakp_df].iloc[m+1]
                    #row_index = globals()[breakp_df].index.tolist()[m]
                    midpoint_diff_sp1=abs(row[sp1+'_stop']-next_row[sp1+'_start'])/2
                    midpoint_value = row[sp1+'_stop']+midpoint_diff_sp1
                    globals()[breakp_df].at[m, 'midpoint_'+sp1]=row['midpoint_'+sp1]=midpoint_value
                    for o in range(len(sp1_insul_k)):
                        row_ins_sp1= sp1_insul_k.iloc[o]
                        if row_ins_sp1[sp1+'_start'] <= midpoint_value <= row_ins_sp1[sp1+'_stop']:
                            globals()[breakp_df].at[n,'insul_score_'+sp1]=row['insul_score_'+sp1]=row_ins_sp1['insul_score']
                if len(globals()[breakp_df])%2==1:
                    last_row = globals()[breakp_df].iloc[-1]
                    midpoint_diff_sp1= abs(next_row[sp1+'_stop']-last_row[sp1+'_start'])
                    midpoint_value = next_row[sp1+'_stop']+midpoint_diff_sp1
                    globals()[breakp_df].at[len(globals()[breakp_df]), 'midpoint_'+sp1]=row['midpoint_'+sp1]=midpoint_value
                    for o in range(len(sp1_insul_k)):
                        row_ins_sp1= sp1_insul_k.iloc[o]
                        if row_ins_sp1[sp1+'_start'] <= midpoint_value <= row_ins_sp1[sp1+'_stop']:
                            globals()[breakp_df].at[len(globals()[breakp_df]),'insul_score_'+sp1]=row['insul_score_'+sp1]=row_ins_sp1['insul_score']
                    
                globals()[breakp_df] = globals()[breakp_df].sort_values([sp2+'_chr',sp2+'_start',sp2+'_stop']).set_index('index')
                

    
            #globals()[breakp_df] = globals()[breakp_df].set_index('index')          
            subset = globals()[breakp_df][['breakpoints','breakpoint_in_species', 'midpoint_'+sp1,'midpoint_'+sp2, 'insul_score_'+sp1, 'insul_score_'+sp2]]
            globals()[pairs] = globals()[pairs].join(subset, how='left')
            
            # insul_sc_sp1_bp = globals()[breakp_df]['insul_score_'+sp1].to_numpy()
            # insul_sc_sp1_bp = [x for x in insul_sc_sp1_bp if not np.isnan(x)]
            # insul_sc_sp2_bp = globals()[breakp_df]['insul_score_'+sp2].to_numpy()
            # insul_sc_sp2_bp = [x for x in insul_sc_sp2_bp if not np.isnan(x)]
            # median_sp1_bp = statistics.median(insul_sc_sp1_bp)
            # sp1_insul_k = sp1_insul_k['insul_score'].to_numpy()
            # sp2_insul_k = sp2_insul_k['insul_score'].to_numpy()
            
            # median_list = []
            # for p in range(num_rounds):
            #     permut_set = np.random.choice(sp1_insul_k, size= len(insul_sc_sp1_bp),replace=False)
            #     permut_median = statistics.median(permut_set)
            #     median_list.append(permut_median)
            # fig, ax = plt.subplots(figsize=(10, 10))
            # plt.hist([median_list, median_sp1_bp], bins=10, density= False)
            # plt.gca().set(title='Histogram for chr '+k, ylabel='Frequency')
            # bigger_than= np.sum(median_list<median_sp1_bp)
            # smaller_than = np.sum(median_list>median_sp1_bp)
            # fdr = smaller_than/num_rounds
            # print(fdr)
            # final_insul_list_sp1.append(insul_sc_sp1_bp)
            # insul_table.append(sp1_insul_k)
            
            insul_sc_sp2_bp = globals()[breakp_df]['insul_score_'+sp2].to_numpy()
            insul_sc_sp2_bp = [x for x in insul_sc_sp2_bp if not np.isnan(x)]
            insul_sc_sp2_bp = globals()[breakp_df]['insul_score_'+sp2].to_numpy()
            insul_sc_sp2_bp = [x for x in insul_sc_sp2_bp if not np.isnan(x)]
            median_sp2_bp = statistics.median(insul_sc_sp2_bp)
            sp2_insul_k = sp2_insul_k['insul_score'].to_numpy()
            sp2_insul_k = sp2_insul_k['insul_score'].to_numpy()
            
            median_list = []
            for p in range(num_rounds):
                permut_set = np.random.choice(sp2_insul_k, size= len(insul_sc_sp2_bp),replace=False)
                permut_median = statistics.median(permut_set)
                median_list.append(permut_median)
            fig, ax = plt.subplots(figsize=(10, 10))
            plt.hist([median_list, median_sp2_bp], bins=10, density= False)
            plt.gca().set(title='Histogram for chr '+k, ylabel='Frequency')
            bigger_than= np.sum(median_list<median_sp2_bp)
            smaller_than = np.sum(median_list>median_sp2_bp)
            fdr = smaller_than/num_rounds
            print(fdr)
            final_insul_list_sp2.append(insul_sc_sp2_bp)
            insul_table.append(sp2_insul_k)
            # for n in range(0,(len(globals()[breakp_df])-1),2):
            #     print(n)
            #     row= globals()[breakp_df].iloc[n]
            #     next_row = globals()[breakp_df].iloc[n+1]
            #     first_row_index = globals()[breakp_df].index.tolist()[0]
            #     for o in range(len(sp1_insul_k)):
            #         row_ins_sp1= sp1_insul_k.iloc[o]
            #         if row_ins_sp1[sp1+'_start'] <= row[sp1+'_stop'] and row[sp1+'_stop'] <= row_ins_sp1[sp1+'_stop']:
            #             globals()[breakp_df].at[first_row_index+n, 'insul_stop_'+sp1]=row['insul_stop_'+sp1]=row_ins_sp1['insul_score']
            #         if row_ins_sp1[sp1+'_start'] <= next_row[sp1+'_start'] and next_row[sp1+'_start'] <= row_ins_sp1[sp1+'_stop']:
            #             globals()[breakp_df].at[first_row_index+n+1, 'insul_start_'+sp1]=next_row['insul_start_'+sp1]=row_ins_sp1['insul_score']      
            #     for p in range(len(sp2_insul_k)):
            #         row_ins_sp2= sp2_insul_k.iloc[p]
            #         if row_ins_sp2[sp2+'_start'] <= row[sp2+'_stop'] and row[sp2+'_stop'] <= row_ins_sp2[sp2+'_stop']:
            #             globals()[breakp_df].at[first_row_index+n, 'insul_stop_'+sp2]=row['insul_stop_'+sp2]=row_ins_sp2['insul_score']
            #         if row_ins_sp2[sp2+'_start'] <= next_row[sp2+'_start'] and next_row[sp2+'_start'] <= row_ins_sp2[sp2+'_stop']:
            #             globals()[breakp_df].at[first_row_index+n+1, 'insul_start_'+sp2]=next_row['insul_start_'+sp2]=row_ins_sp2['insul_score']
                        
            # globals()[breakp_df] = globals()[breakp_df].set_index('index')          
            # subset = globals()[breakp_df][['breakpoints','breakpoint_in_species']]
            #globals()[pairs] = globals()[pairs].join(subset, how='left')
            #globals()[pairs] = globals()[pairs].reset_index().round()
            #stat_test_subset_1 = globals()[pairs].query("breakpoints != 'True'")
            
            # for t in range(0,(len(globals()[pairs])-1),2):
            #     row= globals()[pairs].iloc[t]
            #     next_row = globals()[pairs].iloc[t+1]
            #     first_row_index = globals()[pairs].index.tolist()[0]
            #     for o in range(len(sp1_insul_k)):
            #         row_ins_sp1= sp1_insul_k.iloc[o]
            #         if row_ins_sp1[sp1+'_start'] <= row[sp1+'_stop'] and row[sp1+'_stop'] <= row_ins_sp1[sp1+'_stop']:
            #             globals()[pairs].at[first_row_index+t, 'insul_stop_'+sp1]=row['insul_stop_'+sp1]=row_ins_sp1['insul_score']
            #         if row_ins_sp1[sp1+'_start'] <= row[sp1+'_start'] and row[sp1+'_start'] <= row_ins_sp1[sp1+'_stop']:
            #             globals()[pairs].at[first_row_index+t, 'insul_start_'+sp1]=row['insul_start_'+sp1]=row_ins_sp1['insul_score']
            #         if row_ins_sp1[sp1+'_start'] <= next_row[sp1+'_start'] and next_row[sp1+'_start'] <= row_ins_sp1[sp1+'_stop']:
            #             globals()[pairs].at[first_row_index+t+1, 'insul_start_'+sp1]=next_row['insul_start_'+sp1]=row_ins_sp1['insul_score']    
            #         if row_ins_sp1[sp1+'_start'] <= next_row[sp1+'_stop'] and next_row[sp1+'_stop'] <= row_ins_sp1[sp1+'_stop']:
            #             globals()[pairs].at[first_row_index+t+1, 'insul_stop_'+sp1]=next_row['insul_stop_'+sp1]=row_ins_sp1['insul_score']
            #     for p in range(len(sp2_insul_k)):
            #         row_ins_sp2= sp2_insul_k.iloc[p]
            #         if row_ins_sp2[sp2+'_start'] <= row[sp2+'_stop'] and row[sp2+'_stop'] <= row_ins_sp2[sp2+'_stop']:
            #             globals()[pairs].at[first_row_index+t, 'insul_stop_'+sp2]=row['insul_stop_'+sp2]=row_ins_sp2['insul_score']
            #         if row_ins_sp2[sp2+'_start'] <= row[sp2+'_start'] and row[sp2+'_start'] <= row_ins_sp2[sp2+'_stop']:
            #             globals()[pairs].at[first_row_index+t, 'insul_start_'+sp2]=row['insul_start_'+sp2]=row_ins_sp2['insul_score']
            #         if row_ins_sp2[sp2+'_start'] <= next_row[sp2+'_start'] and next_row[sp2+'_start'] <= row_ins_sp2[sp2+'_stop']:
            #             globals()[pairs].at[first_row_index+t+1, 'insul_start_'+sp2]=next_row['insul_start_'+sp2]=row_ins_sp2['insul_score']    
            #         if row_ins_sp2[sp2+'_start'] <= next_row[sp2+'_stop'] and next_row[sp2+'_stop'] <= row_ins_sp2[sp2+'_stop']:
            #             globals()[pairs].at[first_row_index+t+1, 'insul_stop_'+sp2]=next_row['insul_stop_'+sp2]=row_ins_sp2['insul_score']
            #     if t==(len(globals()[pairs])-1) and (len(globals()[pairs])-1)%2==1:
            #         last_row = globals()[pairs].iloc[-1]
            #         for o in range(len(sp1_insul_k)):
            #             row_ins_sp1= sp1_insul_k.iloc[o]
            #             if row_ins_sp1[sp1+'_start'] <= last_row[sp1+'_stop'] and last_row[sp1+'_stop'] <= row_ins_sp1[sp1+'_stop']:
            #                 globals()[pairs].at[first_row_index+t, 'insul_stop_'+sp1]=last_row['insul_stop_'+sp1]=row_ins_sp1['insul_score']
            #             if row_ins_sp1[sp1+'_start'] <= last_row[sp1+'_start'] and last_row[sp1+'_start'] <= row_ins_sp1[sp1+'_stop']:
            #                 globals()[pairs].at[first_row_index+t, 'insul_start_'+sp1]=last_row['insul_start_'+sp1]=row_ins_sp1['insul_score']
            #         for p in range(len(sp2_insul_k)):
            #             row_ins_sp2= sp2_insul_k.iloc[p]
            #             if row_ins_sp2[sp2+'_start'] <= last_row[sp2+'_stop'] and last_row[sp2+'_stop'] <= row_ins_sp2[sp2+'_stop']:
            #                 globals()[pairs].at[first_row_index+t, 'insul_stop_'+sp2]=last_row['insul_stop_'+sp2]=row_ins_sp2['insul_score']
            #             if row_ins_sp2[sp2+'_start'] <= last_row[sp2+'_start'] and last_row[sp2+'_start'] <= row_ins_sp2[sp2+'_stop']:
            #                 globals()[pairs].at[first_row_index+t, 'insul_start_'+sp2]=last_row['insul_start_'+sp2]=row_ins_sp2['insul_score']
            #     else:
            #         continue
            # globals()[pairs] = globals()[pairs].set_index('index')
            # select random two subsequent rows from the subset
            # if (len(globals()[breakp_df]))%2==0:
            #     num_samples = int(len(globals()[breakp_df])/2)  # set the desired number of random row pairs
            #     sample_indices = np.random.choice(range(len(stat_test_subset_1) - 1), size=num_samples, replace=False)
            #     stat_test_subset = [(stat_test_subset_1.iloc[t], stat_test_subset_1.iloc[t+1]) for t in sample_indices]
            # else:
            #     num_samples = int(len(globals()[breakp_df])//2)
            #     sample_indices = np.random.choice(range(len(stat_test_subset_1) - 1), size=num_samples, replace=False)
            #     if i==(len(globals()[breakp_df])-1):
            #         stat_test_subset = [(stat_test_subset_1.iloc[t], stat_test_subset_1.iloc[t+1], stat_test_subset_1.sample(n=1, random_state=42)) for t in sample_indices]
            #     else:
            #         stat_test_subset = [(stat_test_subset_1.iloc[t], stat_test_subset_1.iloc[t+1]) for t in sample_indices]
            # list_of_tuples = []
            # for s in stat_test_subset:
            #     list_of_tuples.append(pd.DataFrame(s))
            # stat_test_subset = pd.concat(list_of_tuples, axis=0)
            # list_of_test_rows = []
            # #stat_test_subset = stat_test_subset.sample(len(globals()[breakp_df]), random_state=42)
            # #stat_test_subset = pd.DataFrame(stat_test_subset, columns=["Series1", "Series2"])
            # for r in range(0, (len(stat_test_subset))-1, 2):
            #     row= stat_test_subset.iloc[r]
            #     next_row = stat_test_subset.iloc[r+1]
            #     row_index = stat_test_subset.index.tolist()[r]
            #     print(row_index)
            #     for o in range(len(sp1_insul_k)):
            #         row_ins_sp1= sp1_insul_k.iloc[o]
            #         if row_ins_sp1[sp1+'_start'] <= row[sp1+'_stop'] and row[sp1+'_stop'] <= row_ins_sp1[sp1+'_stop']:
            #             globals()[pairs].at[row_index, 'insul_stop_'+sp1]=row['insul_stop_'+sp1]=row_ins_sp1['insul_score']
            #             list_of_test_rows.append(row)
            #         if row_ins_sp1[sp1+'_start'] <= next_row[sp1+'_start'] and next_row[sp1+'_start'] <= row_ins_sp1[sp1+'_stop']:
            #             globals()[pairs].at[row_index+1, 'insul_start_'+sp1]=next_row['insul_start_'+sp1]=row_ins_sp1['insul_score']
            #             list_of_test_rows.append(next_row)
            #     for p in range(len(sp2_insul_k)):
            #         row_ins_sp2= sp2_insul_k.iloc[p]
            #         if row_ins_sp2[sp2+'_start'] <= row[sp2+'_stop'] and row[sp2+'_stop'] <= row_ins_sp2[sp2+'_stop']:
            #             globals()[pairs].at[row_index, 'insul_stop_'+sp2]=row['insul_stop_'+sp2]=row_ins_sp2['insul_score']
            #             list_of_test_rows.append(row)
            #         if row_ins_sp2[sp2+'_start'] <= next_row[sp2+'_start'] and next_row[sp2+'_start'] <= row_ins_sp2[sp2+'_stop']:
            #             globals()[pairs].at[row_index+1, 'insul_start_'+sp2]=next_row['insul_start_'+sp2]=row_ins_sp2['insul_score']
            #             list_of_test_rows.append(next_row)
            #     subset_test = pd.concat(list_of_test_rows, axis=1).transpose()
                        
        #     list_pos_rel_ori=[]
        #     list_neg_rel_ori=[]
        #     #if i == sp1:
        #     fig, ax = plt.subplots(figsize=(10, 10))
        #     for l in range(len(globals()[pairs])):
        #         row = globals()[pairs].iloc[l]
        #         if row['relative_orientation'] == '+':
        #             list_pos_rel_ori.append(row)
        #         else:
        #             list_neg_rel_ori.append(row)
        #     positive_table=pd.concat(list_pos_rel_ori, axis=1)
        #     positive_table= positive_table.transpose()
        #     x1_values=[positive_table[sp1+'_start'], positive_table[sp1+'_stop']]
        #     y1_values=[positive_table[sp2+'_start'], positive_table[sp2+'_stop']]
        #     negative_table=pd.concat(list_neg_rel_ori, axis=1)
        #     negative_table= negative_table.transpose()
        #     x2_values=[negative_table[sp1+'_stop'], negative_table[sp1+'_start']]
        #     y2_values=[negative_table[sp2+'_start'], negative_table[sp2+'_stop']]
        #     plt.plot(x1_values, y1_values, 'b', linestyle="solid")
        #     plt.plot(x2_values, y2_values, 'b', linestyle="solid")
        #     # # connect adjacent points with lines
        #     # for i in range(len(globals()[pairs])-1):
        #     #     row=globals()[pairs].iloc[i]
        #     #     ax.plot(row[sp1+'_start'], row[sp2+'_stop'], 'k-', alpha=0.5)
        #     breakp=[]           
        #     for m in range(len(globals()[pairs])):
        #         row = globals()[pairs].iloc[m]
        #         #z=(row[sp1+'_start'])/(row[sp1+'_chr_size'])
        #         if row['breakpoints']=='True':
        #             ax.scatter(row[sp1+'_start'], row[sp2+'_stop'], color='r')
        #             breakp.append('yes')
        #     #print(k+' vs '+chromosome_pairs[k])
        #     #print(len(breakp)/2+1)
            
        #     ax.ticklabel_format(style='plain', scilimits=(0,0))
        #     # formatter = ticker.ScalarFormatter(useMathText=True)
        #     # formatter.set_powerlimits((-3, 4))  # Set the exponent range to show
        #     # ax.xaxis.set_major_formatter(formatter)
        #     ax.set_xlabel('Octopus bimaculoides chromosome '+k)
        #     ax.set_ylabel('Octopus vulgaris chromosome '+chromosome_pairs[k])
        #     ax.set_title('Chromosome '+k+' vs chromosome '+chromosome_pairs[k])
        # if i==sp1:
        #     pairs_dataframe_fin_sp1 = pd.concat([pairs_dataframe_fin_sp1, globals()[pairs]])
        #     breakpoint_df_sp1 = pd.concat([breakpoint_df_sp1, globals()[breakp_df]])
        #     # stat_test_df_sp1 = pd.concat([stat_test_df_sp1, subset_test]).drop_duplicates()
        # if i==sp2:
        #     pairs_dataframe_fin_sp2 = pd.concat([pairs_dataframe_fin_sp2, globals()[pairs]])
        #     breakpoint_df_sp2 = pd.concat([breakpoint_df_sp2, globals()[breakp_df]])
        #     # stat_test_df_sp2 = pd.concat([stat_test_df_sp2, subset_test]).drop_duplicates()
        #else:
         #continue
    else:
        continue

def flatten(l):
    return [item for sublist in l for item in sublist]
final_insul_list_sp2= flatten(final_insul_list_sp2)
final_insul_list_sp2_median = statistics.median(final_insul_list_sp2)
median_list = []
insul_table = flatten(insul_table)
for p in range(num_rounds):
    permut_set = np.random.choice(insul_table, size= len(final_insul_list_sp2),replace=False)
    permut_median = statistics.median(permut_set)
    median_list.append(permut_median)
fig, ax = plt.subplots(figsize=(10, 10))
plt.hist([median_list, median_sp2_bp], bins=10, density= False)
plt.gca().set(title='Histogram for all', ylabel='Frequency')
bigger_than= np.sum(median_list<median_sp2_bp)
smaller_than = np.sum(median_list>median_sp2_bp)
fdr = smaller_than/num_rounds
print(fdr)
# median_value_1_stop = stat_test_df_sp1['insul_stop_'+sp1].median()
# median_value_1_stop_2 = stat_test_df_sp1['insul_stop_'+sp2].median()
# median_value_bp1 = breakpoint_df_sp1['insul_stop_'+sp1].median()
# median_value_bp1_stop_2 = breakpoint_df_sp1['insul_stop_'+sp2].median()
# # median_value_2 = stat_test_df_sp2['insul_stop_'+sp2].median()
# median_values_bp2 = breakpoint_df_sp2['insul_stop_'+sp1].median()
# median_sp1 = pairs_dataframe_fin_sp1[['insul_stop_'+sp1, 'insul_stop_'+sp2, 'insul_start_'+sp1, 'insul_start_'+sp2]].median()
# median_sp2 = pairs_dataframe_fin_sp2[['insul_stop_'+sp1, 'insul_stop_'+sp2, 'insul_start_'+sp1, 'insul_start_'+sp2]].median()
# pairs_dataframe_fin_sp1.to_csv('table_with_insul_sc')
# plt.hist(pairs_dataframe_fin_sp1[['insul_stop_'+sp1, 'insul_stop_'+sp2, 'insul_start_'+sp1, 'insul_start_'+sp2]], density=False, bins=10)
#plt.hist(pairs_dataframe_fin_sp2[['insul_stop_'+sp1, 'insul_stop_'+sp2, 'insul_start_'+sp1, 'insul_start_'+sp2]], density=False, bins=10)
# plt.hist(stat_test_df_sp1[['insul_start_'+sp2,'insul_stop_'+sp2]], bins=10, density=True, alpha=0.5)
# plt.hist(stat_test_df_sp1[['insul_start_'+sp2,'insul_stop_'+sp2]], bins=10, density=True, alpha=0.5)

# # Create a figure with two subplots
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

# # Plot the first histogram in the first subplot
# ax1.hist(stat_test_df_sp1[['insul_start_'+sp2,'insul_stop_'+sp2]], bins=10, density=False, alpha=0.5)
# ax1.set_xlabel('X Label for subplot 1')
# ax1.set_ylabel('Y Label for subplot 1')
# ax1.set_title('Title for subplot 1')

# # Plot the second histogram in the second subplot
# ax2.hist(breakpoint_df_sp1[['insul_start_'+sp2,'insul_stop_'+sp2]], bins=10, density=False, alpha=0.5)
# ax2.set_xlabel('X Label for subplot 2')
# ax2.set_ylabel('Y Label for subplot 2')
# ax2.set_title('Title for subplot 2')

# # Display the figure
# plt.show()

# # Create a figure with two subplots
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

# # Plot the first histogram in the first subplot
# ax1.hist(stat_test_df_sp2[['insul_start_'+sp1,'insul_stop_'+sp1]], bins=10, density=True, alpha=0.5)
# ax1.set_xlabel('X Label for subplot 1')
# ax1.set_ylabel('Y Label for subplot 1')
# ax1.set_title('Title for subplot 1')

# # Plot the second histogram in the second subplot
# ax2.hist(breakpoint_df_sp2[['insul_start_'+sp1,'insul_stop_'+sp1]], bins=10, density=True, alpha=0.5)
# ax2.set_xlabel('X Label for subplot 2')
# ax2.set_ylabel('Y Label for subplot 2')
# ax2.set_title('Title for subplot 2')

# # Display the figure
# plt.show()

# # Create a figure with two subplots
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

# # Plot the first histogram in the first subplot
# ax1.hist(stat_test_df_sp2[['insul_start_'+sp2,'insul_stop_'+sp2]], bins=10, density=True, alpha=0.5)
# ax1.set_xlabel('X Label for subplot 1')
# ax1.set_ylabel('Y Label for subplot 1')
# ax1.set_title('Title for subplot 1')

# # Plot the second histogram in the second subplot
# ax2.hist(breakpoint_df_sp2[['insul_start_'+sp2,'insul_stop_'+sp2]], bins=10, density=True, alpha=0.5)
# ax2.set_xlabel('X Label for subplot 2')
# ax2.set_ylabel('Y Label for subplot 2')
# ax2.set_title('Title for subplot 2')

# # Display the figure
# plt.show()

# # Create a figure with two subplots
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

# # Plot the first histogram in the first subplot
# ax1.hist(stat_test_df_sp1[['insul_start_'+sp1,'insul_stop_'+sp1]], bins=10, density=True, alpha=0.5)
# ax1.set_xlabel('X Label for subplot 1')
# ax1.set_ylabel('Y Label for subplot 1')
# ax1.set_title('Title for subplot 1')

# # Plot the second histogram in the second subplot
# ax2.hist(breakpoint_df_sp1[['insul_start_'+sp1,'insul_stop_'+sp1]], bins=10, density=True, alpha=0.5)
# ax2.set_xlabel('X Label for subplot 2')
# ax2.set_ylabel('Y Label for subplot 2')
# ax2.set_title('Title for subplot 2')

# Display the figure
plt.show()