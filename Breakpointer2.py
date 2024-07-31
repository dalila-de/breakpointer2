"""
Breakpointer2.py is a script that uses the output of the .paf alignments and the FAN-C insulation score files
to detect breakpoints in the genome-genome alignments (species-species comparison).
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import statistics
import argparse
import matplotlib
from scipy.stats import fisher_exact

parser =argparse.ArgumentParser(description='This script uses the output of the .paf alignments and the FAN-C insulation score files to detect breakpoints in the genome-genome alignments (species-species comparison')
parser.add_argument('input_paf', help='This is the alignment file')
parser.add_argument('input_ins_score_sp1', help= 'FAN-C insulation score for sp1')
parser.add_argument('input_ins_score_sp2', help= 'FAN-C insulation score for sp2')
parser.add_argument('sp1', help = 'species 1 in .paf file')
parser.add_argument('sp2', help = 'species 2 in .paf file')
parser.add_argument('len_co', help='minimum alignment length to be kept', default=100000)
parser.add_argument('len_bp', help = 'distance needed to call a break', default=5000000)
parser.add_argument('q', help = 'quality cutoff value', default=30)
parser.add_argument('num_rounds', help='number of permutations for one-tailed permutation test', default=100000)
parser.add_argument('output_fin_tab', help='Overall table')
args = parser.parse_args()
matplotlib.use('Agg')
#write the name of the query species in the alignment
sp1 = args.sp1
#write the name of the reference species in the alignment
sp2 = args.sp2
#Write the alignment length cutoff for the table filtering
len_co= int(args.len_co)
#Write the quality cutoff for the alignments
q = int(args.q)
#Write what's the distance between two ends of the aligned fragments to call it a breakpoint
len_bp = int(args.len_bp)
#set the number of permutations
num_rounds=int(args.num_rounds)

os.chdir(os.path.dirname(os.path.realpath('Breakpointer2.py')))

#to disable the warning when adding the breakpoints column to the original dataframe
pd.options.mode.chained_assignment = None

#the analysis of the breakpoints in genome-genome alignments

#This part reads in the .paf file as Pandas DataFrame
sp1_to_sp2= pd.read_csv(args.input_paf, delimiter='\t', header=None, names= [sp1+'_chr',sp1+'_chr_size',sp1+'_start',sp1+'_stop','relative_orientation',sp2+'_chr',sp2+'_chr_size',sp2+'_start',sp2+'_stop','nmatch','alen','qual_score','col12','col13','col14','col15','col16','col17','col18'])

#The .paf file has some extra columns that are unnecessary, so this line removes them
sp1_to_sp2= sp1_to_sp2.drop(['col12','col13','col14','col15','col16','col17','col18'], axis =1)

#This one sorts them by the values of 6 columns
sp1_to_sp2= sp1_to_sp2.sort_values([sp1+'_chr',sp1+'_start',sp1+'_stop'])

#This filters the dataframe to include only relevant alignments- the length of the alignment must be higher than the values given in len_co and quality of the alignment must be higher than value q
sp1_to_sp2 = sp1_to_sp2[(sp1_to_sp2.alen > len_co) & (sp1_to_sp2.qual_score > q)]

sp1_to_sp2 = sp1_to_sp2.reset_index()

sp1_to_sp2 = sp1_to_sp2.drop(['index'], axis=1)

sp1_to_sp2_hm = sp1_to_sp2[[sp1+'_chr', sp2+'_chr', 'alen']]

#Group by pairs in the alignment and sum up the alignment lengths
summed_df_sp1 = sp1_to_sp2_hm.groupby([sp1+'_chr', sp2+'_chr'])['alen'].sum().reset_index()
summed_df_sp2 = sp1_to_sp2_hm.groupby([sp2+'_chr', sp1+'_chr'])['alen'].sum().reset_index()


#This part finds the homologous chromosome pairs in two species
def homologous_pairs(summed_df_sp1, sp1, sp2, table_of_homologous_pairs):
    
    dfs = []

    # Iterate over each unique value for sp1+'_chr'
    for chr1_value in summed_df_sp1[sp1+'_chr'].unique():

        # Filter for the current sp1+chr1 value
        chr1_df = summed_df_sp1[summed_df_sp1[sp1+'_chr'] == chr1_value]

        # Find the row with the maximum alignment length
        max_pair = chr1_df.loc[chr1_df['alen'].idxmax()]

        # Sum alignment lengths for all other pairs
        other_entries_sum = chr1_df[chr1_df[sp2+'_chr'] != max_pair[sp2+'_chr']]['alen'].sum()

        # Append a new DataFrame to the list
        dfs.append(pd.DataFrame({
            'Chr1': [max_pair[sp1+'_chr']],
            'Chr2': [max_pair[sp2+'_chr']],
            'Max alignment length': [max_pair['alen']],
            'Sum of other entries': [other_entries_sum]
        }))

    # Concatenate DataFrames in the list
    table_of_homologous_pairs = pd.concat(dfs, ignore_index=True)

    return table_of_homologous_pairs


# Use the function
result_df1 = pd.DataFrame()
result_df2 = pd.DataFrame()
result_df1 = homologous_pairs(summed_df_sp1, sp1, sp2, result_df1)
result_df2 = homologous_pairs(summed_df_sp2, sp2, sp1, result_df2)

#Now get statistical support for the pairs found above
#Fisher's exact test    
#We have to normalize the counts by dividing the number of bases by 10^x (huge numbers)
count_normalization = [10000,100000,1000000]
    
for x in count_normalization:
    FET_x = format("FET_test_normalized_%s" % x)
    df = []
    globals()[FET_x] = pd.DataFrame()
    for index,row in result_df1.iterrows():
        chr1 = row['Chr1']
        chr2 = row['Chr2']
        max_alignment_df1 = row['Max alignment length']
        sum_other_df1 = row['Sum of other entries']
        max_alignment_df2 = result_df2.loc[(result_df2['Chr1'] == chr2) & (result_df2['Chr2'] == chr1), 'Max alignment length'].values[0]
        sum_other_df2 = result_df2.loc[(result_df2['Chr1'] == chr2) & (result_df2['Chr2'] == chr1), 'Sum of other entries'].values[0]
        # Calculate the sum of all other entries in 'Max alignment length' column
        sum_max_alignment_other = result_df1[result_df1['Chr1'] != chr1]['Max alignment length'].sum()
        # Calculate the sum of all other entries in 'Sum of other entries' column
        sum_sum_other_other = result_df1[result_df1['Chr1'] != chr1]['Sum of other entries'].sum()
        #As the number of total bases becomes insane, maybe divide it by the number of chromosomes?
        other_to_other= (sum_max_alignment_other + sum_sum_other_other)/len(result_df1)

        # Construct the contingency matrix
        contingency_matrix = [[int(max_alignment_df1/x+x/x), int(sum_other_df1/x +x/x)], [int(sum_other_df2/x+x/x), int(other_to_other/x+x/x)]]
        df_contingency = pd.DataFrame(contingency_matrix, columns = [chr2, 'Other chromosomes'], index = [chr1, 'Other chromosomes'])
        # Perform Fisher's exact test
        odds_ratio, p_value = fisher_exact(contingency_matrix)
        #Bonferonni correction for FET, 4 tests are performed, because 2x2 contingency matrix is used
        p_adj = p_value*4
        df.append(pd.DataFrame({'Chr1': [chr1], 'Chr2': [chr2], 'Odds ratio': [odds_ratio], 'P-value': [p_value], 'P-value adj': [p_adj], 'Normalization factor': [x]}))
        globals()[FET_x] = pd.concat(df, ignore_index=True)
    globals()[FET_x].to_csv('FET_normalized_%s.tsv' % x, index=False, sep='\t')

list_sp1 = result_df1['Chr1'].tolist()
list_sp2 = result_df1['Chr2'].tolist()

chromosome_pairs = dict(zip(list_sp1, list_sp2))

#.bed file part
bimac_insul_sc_500kb = pd.read_csv(args.input_ins_score_sp1, delimiter='\t', header= None, names=[sp1+'_chr', sp1+'_start',sp1+'_stop', 'irrel1', 'insul_score','irrel2'])
bimac_insul_sc_500kb = bimac_insul_sc_500kb.dropna().drop(['irrel1','irrel2'], axis=1)

vulg_insul_sc_500kb = pd.read_csv(args.input_ins_score_sp2, delimiter='\t', header= None, names=[sp2+'_chr', sp2+'_start',sp2+'_stop', 'irrel1', 'insul_score','irrel2'])
vulg_insul_sc_500kb = vulg_insul_sc_500kb.dropna().drop(['irrel1','irrel2'], axis=1)

"""
because of the gaps in the breakp dataframe, it is not possible to calculate potential midpoints for the opposite species, therefore, a different 
approach is needed
table_i=the merged table for a chromosome, including both breakpoints and non-breakpoints
i_insul_sc= insulation score table for the opposite species
sp = the opposite species
"""
def insul_score_sp(table_i, i_insul_k,sp):
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
    for n in range(0,(len(table_bp)-1),2):
        row = table_bp.iloc[n]
        next_row = table_bp.iloc[n+1]
        midpoint_diff_i = abs(row[sp+'_stop']-next_row[sp+'_start'])/2
        midpoint_value = row[sp+'_stop']+midpoint_diff_i
        table_bp.at[n, 'midpoint_'+sp] = row['midpoint_'+sp] = midpoint_value
        for o in range(len(i_insul_k)):
            row_ins_i = i_insul_k.iloc[o]
            if row_ins_i[sp+'_start'] <= midpoint_value <= row_ins_i[sp+'_stop']:
                table_bp.at[n, 'insul_score_'+sp] = row['insul_score_'+sp] = row_ins_i['insul_score']
    if len(table_bp) % 2 == 1:
        last_row = table_bp.iloc[-1]
        midpoint_diff_i = abs(next_row[sp+'_stop']-last_row[sp+'_start'])
        midpoint_value = next_row[sp+'_stop'] + midpoint_diff_i
        table_bp.at[len(table_bp)-1, 'midpoint_'+sp] = last_row['midpoint_'+sp] = midpoint_value
        for o in range(len(i_insul_k)):
            row_ins_i = i_insul_k.iloc[o]
            if row_ins_i[sp+'_start'] <= midpoint_value <= row_ins_i[sp+'_stop']:
                table_bp.at[len(table_bp)-1, 'insul_score_'+sp] = row['insul_score_'+sp] = row_ins_i['insul_score']
    table_bp = table_bp.set_index('index')
    table_i = table_i.drop(['midpoint_'+sp, 'insul_score_'+sp], axis=1).set_index('index')
    table_bp = table_bp[['midpoint_'+sp, 'insul_score_'+sp]]
    table_i = table_i.join(table_bp, how='left')
    return table_x, table_i

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
    ax.hist(median_list, bins=20, alpha=0.7, label='Permutation Medians')
    if sp == sp1 and len(insul_sc_sp_bp)<80:
        plt.gca().set(title='Histogram for chr ' + k, ylabel='Frequency')
        plt.legend()
        plt.savefig('Histogram_'+k+'.png')
    if sp == sp2 and len(insul_sc_sp_bp)<80:
        plt.gca().set(title='Histogram for chr ' + chromosome_pairs[k], ylabel='Frequency')
        plt.legend()
        plt.savefig('Histogram_'+k+'.png')
    if len(insul_sc_sp_bp)>=80:
        plt.gca().set(title='Histogram for overall chromosomes of ' + sp, ylabel='Frequency')
        plt.legend()
        plt.savefig('Histogram_'+sp+'.png')
    bigger_than= np.sum(median_list>median_sp_bp)
    smaller_than = np.sum(median_list<median_sp_bp)
    fdr = smaller_than/num_rounds
    if fdr == 0:
        corr= 1/num_rounds
        fdr= format(f"<{corr}")
    index = list(range(5))
    if sp == sp1 and len(insul_sc_sp_bp)<=80:
        stats_df = {
            'Dataset' : df_name,
            'Chromosome' : k,
            'Insulation_breakpoint_median' : median_sp_bp,
            'Num_of_Breakpoint_insul_scores' : len(insul_sc_sp_bp),
            'False_detection_rate' : fdr,
            'Number_of_smaller_permutation_medians' : smaller_than,
            'Number_of_bigger_permutation_medians' : bigger_than,
            'Estimated_number of breakpoints': breakp
            }
    if sp == sp2 and len(insul_sc_sp_bp)<=80:
        stats_df = {
            'Dataset' : df_name,
            'Chromosome' : chromosome_pairs[k],
            'Insulation_breakpoint_median' : median_sp_bp,
            'Num_of_Breakpoint_insul_scores' : len(insul_sc_sp_bp),
            'False_detection_rate' : fdr,
            'Number_of_smaller_permutation_medians' : smaller_than,
            'Number_of_bigger_permutation_medians' : bigger_than,
            'Estimated_number of breakpoints': breakp
            }
    if len(insul_sc_sp_bp)>80:
        stats_df = {
            'Species' : sp,
            'Insulation_breakpoint_median' : median_sp_bp,
            'Num_of_Breakpoint_insul_scores' : len(insul_sc_sp_bp),
            'False_detection_rate' : fdr,
            'Number_of_smaller_permutation_medians' : smaller_than,
            'Number_of_bigger_permutation_medians' : bigger_than
            }
    stats_df = pd.DataFrame(data = stats_df, index = index)
    stats_df = stats_df.drop_duplicates()
    return fig, ax, stats_df

pairs_dataframe_fin_sp1 = pd.DataFrame()
pairs_dataframe_fin_sp2 = pd.DataFrame()
stats_df_sp2 = pd.DataFrame()
stats_df_sp1 = pd.DataFrame()
overall_res = pd.DataFrame()
breakpoints_fin_sp1 = pd.DataFrame()
breakpoints_fin_sp2 = pd.DataFrame()
#The for loop
species_list = [sp1,sp2]
#This one sorts them by the values of 6 columns
for i in species_list:
    df_name = format(f"sorted_by_{i}")
    globals()[df_name] = sp1_to_sp2.sort_values([i+'_chr',i+'_start',i+'_stop'])
    final_insul_list_sp1 = []
    final_insul_list_sp2 = []
    insul_table_sp2 = []
    insul_table_sp1 = []
    for k in chromosome_pairs:
        pairs= format(f"table_{i}")
        globals()[pairs]=globals()[df_name].loc[(globals()[df_name][sp1+'_chr']==k) & (globals()[df_name][sp2+'_chr']==chromosome_pairs[k]),]
        sp1_insul_k = bimac_insul_sc_500kb.loc[(bimac_insul_sc_500kb[sp1+'_chr'] == k),]
        sp2_insul_k = vulg_insul_sc_500kb.loc[(vulg_insul_sc_500kb[sp2+'_chr']== chromosome_pairs[k]),]
        #I open this list as it has been shown as the easiest-least-crashiest way to do what I want to do
        the_list = []
        #write an if statement and plot them vice versa for - direction
        globals()[pairs]['breakpoint_in_species'] = np.nan
        globals()[pairs]['breakpoint_in_species']=globals()[pairs]['breakpoint_in_species'].astype(str)
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
            globals()[pairs] = globals()[pairs].sort_values([i+'_chr',i+'_start',i+'_stop'])
            breakp_df = format(f"concat_list_{i}")
            #the_list = pd.Series(the_list).drop_duplicates(keep='last').tolist()
            globals()[breakp_df]= pd.concat(the_list, axis=1)
            
            globals()[breakp_df] = globals()[breakp_df].transpose().drop_duplicates().set_index('index')
            globals()[breakp_df] = globals()[breakp_df][~globals()[breakp_df].index.duplicated(keep='first')].sort_values([i+'_chr',i+'_start',i+'_stop'])
            globals()[pairs] = globals()[pairs].drop(['breakpoint_in_species'], axis=1)
            globals()[breakp_df]['breakpoints']= 'True'
            globals()[breakp_df][['midpoint_'+sp1,'midpoint_'+sp2, 'insul_score_'+sp1, 'insul_score_'+sp2]] = np.nan
            breakp = len(globals()[breakp_df])
            result_table = format(f'{i}_midpoint')
            subset = globals()[breakp_df][['breakpoints','breakpoint_in_species', 'midpoint_'+sp1,'midpoint_'+sp2, 'insul_score_'+sp1, 'insul_score_'+sp2]]
            globals()[pairs] = globals()[pairs].join(subset, how='left')
            
            if i == sp1:
                breakpoints_sp1, globals()[result_table] = insul_score_sp(globals()[pairs], sp1_insul_k,i)
            if i == sp2:
                breakpoints_sp2, globals()[result_table] = insul_score_sp(globals()[pairs], sp2_insul_k,i)
            globals()[pairs]=globals()[result_table]
            
            globals()[pairs] = globals()[pairs].sort_values([i+'_chr',i+'_start',i+'_stop'])
            if i ==sp1:
               breakpoints_sp1, globals()[pairs] = insul_score_sp(globals()[pairs], sp2_insul_k,sp2)
            else:
               breakpoints_sp2, globals()[pairs] = insul_score_sp(globals()[pairs], sp1_insul_k,sp1)
               
            globals()[pairs] = globals()[pairs].sort_values([i+'_chr',i+'_start',i+'_stop'])
            
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
        if i==sp1:
            pairs_dataframe_fin_sp1 = pd.concat([pairs_dataframe_fin_sp1, globals()[pairs]])
            breakpoints_fin_sp1 = pd.concat ([breakpoints_fin_sp1, breakpoints_sp1]).dropna()
            
            #breakpoint_df_sp1 = pd.concat([breakpoint_df_sp1, globals()[breakp_df]])
        if i==sp2:
            pairs_dataframe_fin_sp2 = pd.concat([pairs_dataframe_fin_sp2, globals()[pairs]])
            breakpoints_fin_sp2 = pd.concat ([breakpoints_fin_sp2, breakpoints_sp2]).dropna()
           #breakpoint_df_sp2 = pd.concat([breakpoint_df_sp2, globals()[breakp_df]])
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
pairs_dataframe_fin_sp1.to_csv(args.output_fin_tab)
stats_df_sp1.drop_duplicates(subset=['Chromosome'], keep='first').to_csv('sp1_stats')
stats_df_sp2.drop_duplicates(subset=['Chromosome'], keep='first').to_csv('sp2_stats')
overall_res.to_csv('overall_res')
breakpoints_fin_sp1.to_csv('sp1_breakp.bed', sep = '\t')
breakpoints_fin_sp2.to_csv('sp2_breakp.bed', sep = '\t')
