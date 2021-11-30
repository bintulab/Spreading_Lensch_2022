import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('pdf') # do this because environment does not have GUI backend
import matplotlib.pyplot as plt
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

#plot signal difference for all antibodies 
pal = ['#EC2890', '#325DDB','#F8981C', '#777777']
g = sns.catplot(data=df_dox, x='Element', y='Signal Difference per kb', hue='Antibody', kind='bar',  height=2.8, aspect=1, palette=pal, legend=False)
plt.tight_layout()
plt.xticks(rotation=45)

#save plot
plt.savefig('catplot_K562_5kb_sig_diff_norm_Ac_K4_k9.png')
plt.savefig('catplot_K562_5kb_sig_diff_norm_Ac_K4_k9.pdf')
plt.close()