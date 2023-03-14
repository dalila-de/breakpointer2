# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import os
import matplotlib.pyplot as plt

os.chdir(r'/Users/dalilad/Desktop/The_octopus_code')
#os.chdir(r'/Users/dalilad/Desktop/The_octopus_code')
#print (os.getcwd())

#to disable the warning when adding the breakpoints column to the original dataframe
pd.options.mode.chained_assignment = None

#The .paf file part
#rename the columns- find and replace- sp.name_something
bimac_to_vulg= pd.read_csv('bimaculoides_to_vulgaris.paf', delimiter='\t', header=None, names= ['col1','col2','col3','col4','col5','col6','col7','col8','col9','col10','col11','col12','col13','col14','col15','col16','col17','col18'])

#don't rename the dfs all the time

bimac_to_vulg_sorted=bimac_to_vulg.sort_values(['col1','col3','col4','col6','col8','col9'])

#write magic number at the beginning and name them properly, also qscore or length

test1 = bimac_to_vulg_sorted[(bimac_to_vulg_sorted.col11 > 100000) & (bimac_to_vulg_sorted.col12 > 30)]
#test2 = bimac_to_vulg_sorted[(bimac_to_vulg_sorted.col11 > 500000) & (bimac_to_vulg_sorted.col12 > 30)]

#try to plot as lines (plot each row as a line)
the_list = []
#max_col11 = test1['col11'].max()
#make a column for both axes breakpoints, x start x stop, y start and y stop
n = 1000000
for i in range(len(test1)-1):
    row=test1.iloc[i]
    next_row=test1.iloc[i+1]
    if row['col5'] == '+':
        if next_row['col8'] > row['col9'] + n:
            the_list.append(row)
            the_list.append(next_row)
        #elif next_row['col8'] < row['col9'] + 100000:
         #   the_list.append(row)
          #  the_list.append(next_row)
    if row['col5']== '-':
        if next_row['col9'] > row['col8'] + n:
            the_list.append(row)
            the_list.append(next_row)
        #elif next_row['col8'] < row['col9'] + 100000:
         #   the_list.append(row)
          #  the_list.append(next_row)    
    if row['col5'] != next_row['col5']:
        the_list.append(row)
        the_list.append(next_row)
            
test1_final_table= pd.concat(the_list, axis=1)
test1_final_table = test1_final_table.transpose().drop_duplicates()

test1_breakpoints=test1.isin(test1_final_table).all(axis=1)
test1['breakpoints'] = test1_breakpoints
#for every row in test_final get the ymin and ymax (red bars) and then plot using those ymin and ymax in plot function.

#make it with lines
#draw vertical/horizontal axis with breakpoints
#generate a figure and axis object
fig, ax = plt.subplots(figsize=(25, 25))
test_grouped=test1.groupby(['col1', 'col6'])
#filters row based on the two col names
#a=unique(test1['col1'])
#b=unique(test1['col6'])
plotdf=test1.loc[(test1['col1']=='CM046601.1') & (test1['col6']=='Ovu01'),]
#for i, group in test_grouped:
ax.scatter(plotdf['col3'], plotdf['col9'], alpha=0.5)


#plt.axvline(x=39786979, ymin = 0.20, ymax = 0.30, color = "red", linewidth=0.5)
# z= (39756707/225692979)
# plt.axvline(x=39756707, ymin = (z-0.08), ymax = (z+0.04), color = "red", linewidth=0.5)
# for l in range(len(test1_final_table)):
#     row = test1_final_table.iloc[l]
#     z=(row['col4'])/(row['col2'])
#     plt.axvline(x=row['col4'], ymin = (z-0.04), ymax = (z+0.03), color = "red", linewidth=0.5)
    
for i in range(len(test1_final_table)):
    row=test1_final_table.iloc[i]
        #ax.stem(row['col3'], row['col9'], markerfmt ='')
    #normalized_data_row=(row['col4']-test1_final_table['col4'].min())/(plotdf['col4'].max()-plotdf['col4'].min())
    #normalized_data_row_bim=(row['col9']-plotdf['col9'].min())/(plotdf['col9'].max()-plotdf['col9'].min())
    #plt.axvline(x=row['col3'] & row['col4'], ymin = 0, ymax = row['col4'], color = "red", linewidth=0.5)
    plt.axvline(x=row['col8'], ymin = 0.25, ymax = 0.75, color = "red", linewidth=0.5)
        #plt.vlines(row['col3'] & row['col4'], ymin=0, ymax= row['col4'].max(), colors=None, linestyles='solid', label='')
# ax.scatter(test1['col4'], test1['col9'], alpha=0.5)
# ax.set_xlim(0, plotdf['col8'].max())
# ax.set_ylim(0, plotdf['col9'].max())
# for i in range(len(plotdf)-1):
#     row=plotdf.iloc[i]
#     if plotdf[row,'breakpoints']==True:
#       plt.axvline(x=plotdf['col3'] & plotdf['col9'], ymin = 0, ymax = plotdf['col8'].max(), color = "red")
#set axis labels and title
ax.set_xlabel('Bimac')
ax.set_ylabel('Vulgaris')
ax.set_title('Homology')

# add a legend
ax.legend()

# show the plot
plt.show()
#plt.save can pick pdf
#goal to have a D-Genies plot I can zoom in on

#.gff or .bed file part

# bimac_insul_sc_500gb= pd.read_csv('oct_bim_0_diff_wind_500kb.bed', delimiter='\t', header= None, names= ['col1','col2','col3','col4','col5','col6'])
# test3= bimac_insul_sc_500gb.dropna()

# vulg_insul_sc_500gb= pd.read_csv('fin_oct_vulg.bed_500kb.bed', delimiter='\t', header= None, names= ['col1','col2','col3','col4','col5','col6'])
# test4= vulg_insul_sc_500gb.dropna()

# # for i in range(len(test1)-1):
# #     row=test1.iloc[i] 
# #     if test1.loc[i,'breakpoints'] == True:
# #         test1['place']=test1['col3'].between(test3['col2'], test3['col3'])
# #         if test1['place'] == True:
# #             test1['insul_score']=test3['col5']

# for index, row in test1.iterrows():
#     if row['breakpoints'] == True:
#         insul_score_query = test3.query(f"{row['col3']} >= col2 and {row['col3']} <= col3")['col5']
#         if not insul_score_query.empty:
#             insul_score = insul_score_query.iloc[0]
#             test1.loc[index, 'insul_score_b'] = insul_score
#         else:
#             test1.loc[index, 'insul_score_b'] = 'NaN'
            
# for index, row in test1.iterrows():
#     if row['breakpoints'] == True:
#         insul_score_query = test4.query(f"{row['col8']} >= col2 and {row['col3']} <= col3")['col5']
#         if not insul_score_query.empty:
#             insul_score = insul_score_query.iloc[0]
#             test1.loc[index, 'insul_score_v'] = insul_score
#         else:
#             test1.loc[index, 'insul_score_v'] = 'NaN'


            