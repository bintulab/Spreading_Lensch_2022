import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('pdf') # do this because environment does not have GUI backend
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['pdf.fonttype'] = 42

#read in csv with signal 
df = pd.read_csv("bedgraphAUC.csv")

#split dataframe by dox/no dox
df_nodox = df[df['Condition']=="nodox"]
print(df_nodox)
df_dox = df[df['Condition']=="dox"]
print(df_dox)


#calculate signal difference 
sig_diff = np.array(df_dox['Signal'].values.tolist())-np.array(df_nodox['Signal'].values.tolist())
print(sig_diff)
df_dox["Signal_Difference"]=sig_diff
#normalize signal difference by length of element 
df_dox['Length'] = df_dox['End'].astype(int) - df_dox['Start'].astype(int)
df_dox['Signal Difference per kb'] = 1000*df_dox['Signal_Difference']/df_dox['Length']

#subset dataframe to antibodies to be plotted
df_subset = df_dox[df_dox.Antibody != 'H3K4me3']
df_subset = df_subset[df_subset.Antibody != 'IgG']

#configure figure size 
fig,  ax = plt.subplots(figsize = (2.5,2.5))

#rename elements 
# df_dox = df_dox.replace(to_replace = '5kb-spacer', value ='5kb')

#plot signal difference for different antibodies on the same plot 
f = sns.barplot(data = df_dox[df_dox.Antibody == 'H3Kac'], ax=ax, x='Element', y='Signal Difference per kb', color='#F8981C')
g = sns.barplot(data = df_dox[df_dox.Antibody == 'H3K9me3'], ax=ax, x='Element', y='Signal Difference per kb', color = '#325DDB')

#configure axes 
#f.set_ylim(-28,20)
#g.set_ylim(-28,20)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

#configure legend
legend_elements = [Line2D([0],[0], marker='$\u25AC$', color = 'white', label='H3Kac', markerfacecolor='#F8981C', markersize =13),
					Line2D([0],[0], marker='$\u25AC$', color = 'white', label='H3K9me3', markerfacecolor='#325DDB', markersize =13),]
ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(0.3, 1), frameon=False) 
plt.axhline(y=0, linewidth=1, c='black')
plt.xticks(rotation=45)

#save plot 
plt.tight_layout()
plt.savefig('barplot_K562_5kb_sig_diff_norm_Ac_K9.png')
plt.savefig('barplot_K562_5kb_sig_diff_norm_Ac_K9.pdf')
plt.close()