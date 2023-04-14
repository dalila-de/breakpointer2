# -*- coding: utf-8 -*-
"""
To be written
use this code to get the scaffold names for the assemblies: grep '^>' fin_oct_vulg_30_scaff.fa | sed -E 's/^[[:space:]]*[^[:space:]]+/&,/' | sed 's/>//g' | sed 's/[^[:space:],]\+/"&"/g' | tr '\n' ' '
"""
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

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
sp1_to_sp2= sp1_to_sp2.sort_values([sp1+'_chr',sp1+'_start',sp1+'_stop',sp2+'_chr',sp2+'_start',sp2+'_stop'])

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
bimac_insul_sc_500kb = pd.read_csv('oct_bim_0_diff_wind_500kb.bed', delimiter='\t', header= None, names=[sp1+'_chr', sp1+'_start',sp1+'_stop', 'irrel1', 'insul_score','irrel2'])
bimac_insul_sc_500kb = bimac_insul_sc_500kb.dropna().drop(['irrel1','irrel2'], axis=1)

vulg_insul_sc_500kb = pd.read_csv('fin_oct_vulg.bed_500kb.bed', delimiter='\t', header= None, names=[sp2+'_chr', sp2+'_start',sp2+'_stop', 'irrel1', 'insul_score','irrel2'])
vulg_insul_sc_500kb = vulg_insul_sc_500kb.dropna().drop(['irrel1','irrel2'], axis=1)

pair_dataframe_fin = pd.DataFrame()

#The for loop
for k in chromosome_pairs:
    pair_dataframe=sp1_to_sp2.loc[(sp1_to_sp2[sp1+'_chr']==k) & (sp1_to_sp2[sp2+'_chr']==chromosome_pairs[k]),]

    #I open this list as it has been shown as the easiest-least-crashiest way to do what I want to do
    the_list = []
    #bed file starts and stops: the start is the least num and the stop is higher (regardless of the strand)
    #write an if statement and plot them vice versa for - direction
    pair_dataframe['breakpoint_in_species'] = np.nan
    pair_dataframe['breakpoint_in_species']=pair_dataframe['breakpoint_in_species'].astype(str)
    #bed file starts and stops: the start is the least num and the stop is higher (regardless of the strand)
    #write an if statement and plot them vice versa for - direction
    for i in range(len(pair_dataframe)-1):
        row = pair_dataframe.iloc[i]
        next_row=pair_dataframe.iloc[i+1]
        if next_row[sp2+'_start'] > (row[sp2+'_stop'] + len_bp):
            pair_dataframe.at[i,'breakpoint_in_species']==sp2
            pair_dataframe.at[i+1,'breakpoint_in_species']==sp2
            the_list.append(row)
            the_list.append(next_row)
        if next_row[sp2+'_start'] > (row[sp2+'_start'] + len_bp):
            pair_dataframe.at[i,'breakpoint_in_species']==sp2
            pair_dataframe.at[i+1,'breakpoint_in_species']==sp2
            the_list.append(row)
            the_list.append(next_row)
        if next_row[sp2+'_start'] < (row[sp2+'_stop'] - len_bp):
            pair_dataframe.at[i,'breakpoint_in_species']==sp2
            pair_dataframe.at[i+1,'breakpoint_in_species']==sp2
            the_list.append(row)
            the_list.append(next_row)
        if next_row[sp2+'_stop'] > (row[sp2+'_start'] + len_bp):
            pair_dataframe.at[i,'breakpoint_in_species']==sp2
            pair_dataframe.at[i+1,'breakpoint_in_species']==sp2
            the_list.append(row)
            the_list.append(next_row)
        if next_row[sp1+'_start'] < (row[sp1+'_stop'] - len_bp):
            pair_dataframe.at[i,'breakpoint_in_species']==sp1
            pair_dataframe.at[i+1,'breakpoint_in_species']==sp1
            the_list.append(row)
            the_list.append(next_row)
        if next_row[sp1+'_start'] > (row[sp1+'_start'] + len_bp):
            pair_dataframe.at[i,'breakpoint_in_species']==sp1
            pair_dataframe.at[i+1,'breakpoint_in_species']==sp1
            the_list.append(row)
            the_list.append(next_row)
        if next_row[sp1+'_stop'] > (row[sp1+'_start'] + len_bp):
            pair_dataframe.at[i,'breakpoint_in_species']==sp1
            pair_dataframe.at[i+1,'breakpoint_in_species']==sp1
            the_list.append(row)
            the_list.append(next_row)
        if next_row[sp1+'_start'] > (row[sp1+'_stop'] + len_bp):
            pair_dataframe.at[i,'breakpoint_in_species']==sp1
            pair_dataframe.at[i+1,'breakpoint_in_species']==sp1
            the_list.append(row)
            the_list.append(next_row)
    
    #Seems like breakpoint can be detected in both, should be aware of that            
    sp1_to_sp2_final_table= pd.concat(the_list, axis=1)
    sp1_to_sp2_final_table = sp1_to_sp2_final_table.transpose().drop_duplicates()
    
    sp1_to_sp2_breakpoints=sp1_to_sp2.isin(sp1_to_sp2_final_table).all(axis=1)
    pair_dataframe['breakpoints'] = sp1_to_sp2_breakpoints
    #for every row in test_final get the ymin and ymax (red bars) and then plot using those ymin and ymax in plot function.
    
    #make it with lines
    #draw vertical/horizontal axis with breakpoints
    #generate a figure and axis object
    #I wrote this incredibly complicated part just to make sure that
    #it's plotted properly according to the relative orientation- without it looks strange
    fig, ax = plt.subplots(figsize=(10, 10))
    list_pos_rel_ori=[]
    list_neg_rel_ori=[]
    
    for i in range(len(pair_dataframe)):
        row = pair_dataframe.iloc[i]
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
    # for i in range(len(pair_dataframe)-1):
    #     row=pair_dataframe.iloc[i]
    #     ax.plot(row[sp1+'_start'], row[sp2+'_stop'], 'k-', alpha=0.5)
    breakp=[]           
    for i in range(len(pair_dataframe)):
        row = pair_dataframe.iloc[i]
        #z=(row[sp1+'_start'])/(row[sp1+'_chr_size'])
        if row['breakpoints']:
            ax.scatter(row[sp1+'_start'], row[sp2+'_stop'], color='r')
            breakp.append('yes')
    print(k+' vs '+chromosome_pairs[k])
    print(len(breakp)/2+1)
    
    ax.ticklabel_format(style='plain', scilimits=(0,0))
    # formatter = ticker.ScalarFormatter(useMathText=True)
    # formatter.set_powerlimits((-3, 4))  # Set the exponent range to show
    # ax.xaxis.set_major_formatter(formatter)
    ax.set_xlabel('Octopus bimaculoides chromosome '+k)
    ax.set_ylabel('Octopus vulgaris chromosome '+chromosome_pairs[k])
    ax.set_title('Chromosome '+k+' vs chromosome '+chromosome_pairs[k])
    
    plt.show()
    
    #the .bed file part
    #pandas.query function if this value is this or that----filters for you, doesnt iterate, dont need a for loop
    for i in range(len(pair_dataframe)):
        row = pair_dataframe.iloc[i]
        if row['breakpoints']:
            for j in range(len(bimac_insul_sc_500kb)):
                insul_row=bimac_insul_sc_500kb.iloc[j]
                if (insul_row[sp1+'_start'] + 5000) < row[sp1+'_stop'] and (insul_row[sp1+'_start'] - 5000) > row[sp1+'_stop']:
                    print(insul_row['insult_score'])
    
    pair_dataframe_fin = pd.concat([pair_dataframe_fin, pair_dataframe])
    






#list of ranges-----each breakpoint left side and right side--- both species
#plot horizontal or vertical line between the breakpoints(midpoint)




#plt.save can pick pdf
#goal to have a D-Genies plot I can zoom in on


