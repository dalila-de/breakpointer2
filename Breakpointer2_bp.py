import os
import argparse
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt

parser =argparse.ArgumentParser(description='This script shall process the .paf file and find breakpoints)')
parser.add_argument('input_paf', help='This is the alignment file')
parser.add_argument('sp1', help = 'species 1 in .paf file')
parser.add_argument('sp2', help = 'species 2 in .paf file')
parser.add_argument('len_co', help='minimum alignment length to be kept', default=100000)
parser.add_argument('len_bp', help = 'distance needed to call a break', default=5000000)
parser.add_argument('q', help = 'quality cutoff value', default=30)

args = parser.parse_args()

#set values like this because why not

#Input .paf file
input_paf = args.input_paf
#write the name of the query species in the alignment
sp1 = args.sp1
#write the name of the reference species in the alignment
sp2 = args.sp2
#Write the alignment length cutoff for the table filtering- this step might need some optimization depending on the data
len_co= int(args.len_co)
#Write the quality cutoff for the alignments
q = int(args.q)
#Write what's the distance between two ends of the aligned fragments to call it a breakpoint
len_bp = int(args.len_bp)
#write the output file
output_fin_tab = sp1+'_vs_'+sp2+'_breakpoints.tsv'

os.chdir(os.path.dirname(os.path.realpath('Breakpointer2.py')))

#to disable the warning when adding the breakpoints column to the original dataframe
pd.options.mode.chained_assignment = None

#The .paf file part, aka the analysis of the breakpoints in genome-genome alignments

#This part reads in the .paf file as Pandas DataFrame
sp1_to_sp2= pd.read_csv(input_paf, delimiter='\t', header=None, names= [sp1+'_chr',sp1+'_chr_size',sp1+'_start',sp1+'_stop','relative_orientation',sp2+'_chr',sp2+'_chr_size',sp2+'_start',sp2+'_stop','nmatch','alen','qual_score','col12','col13','col14','col15','col16','col17','col18'])

#The .paf file has some extra columns which are unnecessary, so this line removes them
sp1_to_sp2= sp1_to_sp2.drop(['col12','col13','col14','col15','col16','col17','col18'], axis =1)

#This one sorts the by the values of 6 columns
sp1_to_sp2= sp1_to_sp2.sort_values([sp1+'_chr',sp1+'_start',sp1+'_stop'])

#This filters the dataframe to include only relevant alignments- the length of the alignment must be higher than the values given in len_co and wuality of the alignment must be higher than value q
sp1_to_sp2 = sp1_to_sp2[(sp1_to_sp2.alen > len_co) & (sp1_to_sp2.qual_score > q)]

sp1_to_sp2 = sp1_to_sp2.reset_index()

sp1_to_sp2 = sp1_to_sp2.drop(['index'], axis=1)

sp1_to_sp2_hm = sp1_to_sp2[[sp1+'_chr', sp2+'_chr', 'alen']]
# Group by pairs in the alignment and sum up the alignment lengths
summed_df_sp1 = sp1_to_sp2_hm.groupby([sp1+'_chr', sp2+'_chr'])['alen'].sum().reset_index()

#This part finds the homologous chromosome pairs in two species
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
    table_of_homologous_pairs.to_csv('FET_test.tsv', index=False, sep='\t')
    return table_of_homologous_pairs

def breakp_range(table_i,sp):
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
    table_bp = pd.concat(sublist_fr_ins_sc, axis=1).transpose().drop_duplicates()
    table_bp = table_bp.set_index('index').reset_index()
    one_row_data = pd.DataFrame()
    table_x = pd.DataFrame()
    for n in range(0,(len(table_bp)-1),2):
        row = table_bp.iloc[n]
        next_row = table_bp.iloc[n+1]
        one_row_data = pd.concat([row[[sp+'_chr', sp+'_stop']], next_row[[sp+'_start']]], axis=1).transpose()
        one_row_data.at[n, sp + '_start'] = next_row[sp+'_start']
        table_x = pd.concat([table_x, one_row_data])
    return table_x

# USe the function to get the homologous pairs

result_df1 = pd.DataFrame()
result_df1 = fisher_exact_func(summed_df_sp1, sp1, sp2, result_df1, 1000000)
result_df1 = result_df1[result_df1['P-adj'] < 0.01]
list_sp1 = result_df1['Chr1'].tolist()
list_sp2 = result_df1['Chr2'].tolist()

chromosome_pairs = dict(zip(list_sp1, list_sp2))

def process_species(pairs, species_list, len_bp, sp1, sp2):
    the_list = []
    
    for z in species_list:
        pairs = pairs.sort_values([z+'_chr', z+'_start', z+'_stop']).reset_index()
        for j in range(len(pairs) - 1):
            row = pairs.iloc[j]
            next_row = pairs.iloc[j+1]
            first_row_index = pairs.index.tolist()[0]

            if next_row[sp2+'_start'] < (row[sp2+'_stop'] - len_bp) or \
               next_row[sp2+'_start'] > (row[sp2+'_stop'] + len_bp) or \
               next_row[sp1+'_start'] < (row[sp1+'_stop'] - len_bp) or \
               next_row[sp1+'_start'] > (row[sp1+'_stop'] + len_bp):
                
                pairs.at[first_row_index+j, 'breakpoint_in_species'] = row['breakpoint_in_species'] = z
                pairs.at[first_row_index+j+1, 'breakpoint_in_species'] = next_row['breakpoint_in_species'] = z
                the_list.append(row)
                the_list.append(next_row)

        pairs = pairs.set_index('index')
    
    return pairs, the_list

pairs_dataframe_fin_sp1 = pd.DataFrame()
pairs_dataframe_fin_sp2 = pd.DataFrame()
stats_df_sp2 = pd.DataFrame()
stats_df_sp1 = pd.DataFrame()
overall_res = pd.DataFrame()
breakpoints_fin_sp1 = pd.DataFrame()
breakpoints_fin_sp2 = pd.DataFrame()
#The for loop
species_list = [sp1,sp2]
#This one sorts the by the values of 6 columns
for i in species_list:
    df_name = format(f"sorted_by_{i}")
    globals()[df_name] = sp1_to_sp2.sort_values([i+'_chr',i+'_start',i+'_stop'])
    final_insul_list_sp1 = []
    final_insul_list_sp2 = []
    for k in chromosome_pairs:
        pairs= format(f"table_{i}")
        globals()[pairs]=globals()[df_name].loc[(globals()[df_name][sp1+'_chr']==k) & (globals()[df_name][sp2+'_chr']==chromosome_pairs[k]),]
        #globals()[pairs] = globals()[pairs].reset_index()
        #I open this list as it has been shown as the easiest-least-crashiest way to do what I want to do
        the_list = []
        #bed file starts and stops: the start is the least num and the stop is higher (regardless of the strand)
        #write an if statement and plot them vice versa for - direction
        globals()[pairs]['breakpoint_in_species'] = np.nan
        globals()[pairs]['breakpoint_in_species']=globals()[pairs]['breakpoint_in_species'].astype(str)
        #bed file starts and stops: the start is the least num and the stop is higher (regardless of the strand)
        #write an if statement and plot them vice versa for - direction
        #k = "CM046602.1"
        globals()[pairs], the_list = process_species(globals()[pairs], species_list, len_bp, sp1, sp2)       
        #Seems like breakpoint can be detected in both, should be aware of that
        if len(the_list)>0:
            globals()[pairs] = globals()[pairs].sort_values([i+'_chr',i+'_start',i+'_stop'])
            breakp_df = format(f"concat_list_{i}")
            #the_list = pd.Series(the_list).drop_duplicates(keep='last').tolist()
            globals()[breakp_df]= pd.concat(the_list, axis=1)
            globals()[breakp_df] = globals()[breakp_df].transpose().drop_duplicates().set_index('index')
            globals()[breakp_df] = globals()[breakp_df][~globals()[breakp_df].index.duplicated(keep='first')].sort_values([i+'_chr',i+'_start',i+'_stop'])
            globals()[pairs] = globals()[pairs].drop(['breakpoint_in_species'], axis=1)
            globals()[breakp_df]['breakpoints']= 'True'
        
            #globals()[breakp_df] = globals()[breakp_df].reset_index()
            subset = globals()[breakp_df][['breakpoints','breakpoint_in_species']]
            globals()[pairs] = globals()[pairs].join(subset, how='left')
            bed_file = format(f"bed_file_{i}")
            globals()[bed_file] = breakp_range(globals()[pairs],i)
        
            list_pos_rel_ori=[]
            list_neg_rel_ori=[]
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

            for m in range(len(globals()[pairs])-1):
                row = globals()[pairs].iloc[m] 
                #z=(row[sp1+'_start'])/(row[sp1+'_chr_size'])
                if row['breakpoints']=='True':
                    ax.scatter(row[sp1+'_start'], row[sp2+'_stop'], color='r')
        
            ax.ticklabel_format(style='plain', scilimits=(0,0))
            ax.set_xlabel(sp1+'_'+k)
            ax.set_ylabel(sp2+'_'+chromosome_pairs[k])
            ax.set_title('Chromosome '+k+' vs chromosome '+chromosome_pairs[k])
            plt.savefig('comparison_'+k+'.png')
        if i==sp1:
            pairs_dataframe_fin_sp1 = pd.concat([pairs_dataframe_fin_sp1, globals()[pairs]])
            breakpoints_fin_sp1 = pd.concat([breakpoints_fin_sp1, globals()[bed_file]]).dropna()
        if i==sp2:
            pairs_dataframe_fin_sp2 = pd.concat([pairs_dataframe_fin_sp2, globals()[pairs]])
            breakpoints_fin_sp2 = pd.concat([breakpoints_fin_sp2, globals()[bed_file]]).dropna()
        else:
            continue
pairs_dataframe_fin_sp1.to_csv(output_fin_tab, sep='\t', index=False)
breakpoints_fin_sp1.reset_index(drop=True, inplace=True)
breakpoints_fin_sp2.reset_index(drop=True, inplace=True)
breakpoints_fin_sp1.to_csv(sp1+'_breakp.bed', sep = '\t', header = False)
breakpoints_fin_sp2.to_csv(sp2+'_breakp.bed', sep = '\t', header = False)