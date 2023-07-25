#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 13:01:26 2023

@author: dalilad
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
#from scipy.stats import ttest_ind
#from statsmodels.stats.multitest import fdrcorrection

sp1='Obi'
sp2='Ovu'
#make sure to know where you're working, can be modified later as an input
os.chdir(r'/Users/dalilad/Desktop/The_octopus_code/breakpointer2-main')
dataframe = pd.read_csv('table_with_insul_sc')

dataframe_1 = dataframe.query('breakpoints==True')
dataframe_2 = dataframe.query('breakpoints!=True')
stat_test_dataframe = pd.DataFrame()
the_list=[]
dataframe_1 = dataframe_1.reset_index()
dataframe_2 = dataframe_2.reset_index()
for i in range(0,(len(dataframe_1))-1,2):
    row = dataframe_1.iloc[i]
    next_row = dataframe_1.iloc[i+1]
    the_list.append(row[['insul_stop_'+sp1, 'insul_stop_'+sp2]].tolist() + next_row[['insul_start_'+sp1, 'insul_start_'+sp2]].tolist())
the_other_list = []
for j in range(0, len(dataframe_2)-1):
    row = dataframe_2.iloc[j]
    next_row = dataframe_2.iloc[j+1]
    the_other_list.append(row[['insul_stop_'+sp1, 'insul_stop_'+sp2]].tolist() + next_row[['insul_start_'+sp1, 'insul_start_'+sp2]].tolist())

plt.hist(dataframe[['insul_stop_'+sp1, 'insul_stop_'+sp2,'insul_start_'+sp1, 'insul_start_'+sp2 ]], bins=10, density=True, alpha=0.5)
plt.hist(dataframe_1[['insul_stop_'+sp1, 'insul_stop_'+sp2,'insul_start_'+sp1, 'insul_start_'+sp2 ]], bins=10, density=True, alpha=0.5)
plt.hist(dataframe_2[['insul_stop_'+sp1, 'insul_stop_'+sp2,'insul_start_'+sp1, 'insul_start_'+sp2 ]], bins=10, density=True, alpha=0.5)
plt.hist(the_other_list,bins=10, density=True, alpha=0.5)

alpha = 0.05

flat_list=[]

for i in the_list:
    for item in i:
        flat_list.append(item)
plt.hist(flat_list, bins=10, alpha=0.5)

the_other_flat_list = []

for k in the_other_list:
    for item in k:
        the_other_flat_list.append(item)
plt.hist(the_other_flat_list, bins=10, alpha=0.5)

 
#list1_subset = random.sample(the_other_list, len(the_list))

# perform a two-sample t-test on each pair of sublists
# t_values = []
# for i in range(len():
#     t, p = ttest_ind(list1[i], list2[i])
#     t_values.append(t)

# # perform permutation test with 1000 repeats
# t_permuted = []
# for i in range(1000):
#     permuted_list1 = np.random.permutation(list1)[:len(list2)]
#     permuted_t = []
#     for j in range(len(permuted_list1)):
#         t, p = ttest_ind(permuted_list1[j], list2[j])
#         permuted_t.append(t)
#     t_permuted.append(permuted_t)

# # calculate p-values from the permuted t-values
# p_values_permuted = []
# for i in range(len(list1)):
#     permuted_t_values = [t[i] for t in t_permuted]
#     p = (sum(permuted_t_values >= t_values[i]) + 1) / (len(permuted_t_values) + 1)
#     p_values_permuted.append(p)

# # apply FDR correction to the p-values
# reject, p_values_fdr = fdrcorrection(p_values_permuted, alpha=alpha)

# print("t-values:", t_values)
# print("p-values (permuted):", p_values_permuted)
# print("p-values (FDR corrected):", p_values_fdr)
# print("significant sublists:", [i for i in range(len(list1)) if reject[i]])