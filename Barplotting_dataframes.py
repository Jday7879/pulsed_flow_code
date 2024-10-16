# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 14:37:35 2022

@author: Jordan
"""

import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd
import seaborn as sns
sns.set_style("whitegrid")
sns.set_context("talk")
df = pd.read_csv('dataframe_pulse_cont_summary.csv')
df['Removal Rate'] = (df['Influent Substrate'] - df['Effluent Substrate'])/df['Influent Substrate'] *100
case = []
for i in np.arange(df.shape[0]):
    case.append('S Influent: ' + str(df['Influent Substrate'][i]) + '\n Baffles: ' + str(df['Number of Baffles'][i]))
    
df['Case'] = case
df['Removal per Anode'] = df['Removal Rate']/df['Number of Baffles']

plt.close('all')
plt.figure()
sns.barplot(data = df, x = 'Flow', y = 'Current Density', hue = 'Number of Baffles')
plt.figure()
sns.barplot(data = df, x = 'Flow', y = 'Effluent Substrate', hue = 'Number of Baffles')

plt.figure()
sns.barplot(data = df, x = 'Flow', y = 'Current Density', hue = 'Influent Substrate')
plt.figure()
sns.barplot(data = df, x = 'Flow', y = 'Effluent Substrate', hue = 'Influent Substrate')

plt.figure(figsize = (15,14))
sns.set(style="darkgrid")
sns.barplot(data = df, x = 'Case', y = 'Current Density', hue = 'Flow')
plt.ylabel('Current Density (A/m$^2$)',fontsize = 36)
plt.xlabel('Case',fontsize = 36)
plt.tick_params(axis='both', labelsize=28)
plt.legend(fontsize = 34)
plt.tight_layout()
    


plt.figure(figsize = (15,14))
sns.set(style="darkgrid")
sns.barplot(data = df, x = 'Case', y = 'Removal Rate', hue = 'Flow')
plt.ylabel('Substrate Removal (%)',fontsize = 36)
plt.xlabel('Case',fontsize = 36)
plt.tick_params(axis='both', labelsize=28)
plt.tight_layout()
plt.legend(fontsize = 34)
    