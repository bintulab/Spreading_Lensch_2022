import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('pdf') # do this because environment does not have GUI backend
import matplotlib.pyplot as plt
import seaborn as sns

#coordinates of elements to be quanitified 
elemDict = {'PuroR':[119, 718], 'TRE':[1035, 1358], 'mCit':[3085, 3801], '5kb' :[4290,9305],
			 'pRSV':[9324, 9586], 'mCh':[9996, 10706]}

#save input directory and initialize dictionary of elements
inputDir = sys.argv[1]
initDict = {'PuroR':0, 'TRE':0, 'mCit':0, '5kb':0,
			'pRSV':0, 'mCh':0}

#initialize dictionary of sample conditions (antibody, dox)
sigDict = {'H3K9me3_5kb_dox':initDict.copy(),
		   'H3K9me3_5kb_nodox':initDict.copy(),
		   'H3K4me3_5kb_dox':initDict.copy(),
		   'H3K4me3_5kb_nodox':initDict.copy(),
		   'IgG_5kb_dox':initDict.copy(),
		   'IgG_5kb_nodox':initDict.copy(),
		   'H3Kac_5kb_dox':initDict.copy(),
		   'H3Kac_5kb_nodox':initDict.copy()}

#read in bedgraph files and sum signal at different elements in various sample conditions 
for file in os.listdir(inputDir):
	if file[-8:] != 'bedgraph':
		continue
	#if ('IgG' in file) | ('H3K4me3' in file):
	#	continue
	filepath = os.path.join(inputDir, file)
	with open(filepath, 'r') as bedgraph:
		print(file)
		file_contents = file.split('_')
		spacer = file_contents[3]
		condition = file_contents[5]
		antibody = file_contents[6]
		sigDictkey = '_'.join([antibody, spacer, condition])
		for line in bedgraph:
			start = int(line.rstrip('\n').split('\t')[1])
			end = int(line.rstrip('\n').split('\t')[2])
			signal = float(line.rstrip('\n').split('\t')[3])
			for key in elemDict:
				if (start >= elemDict[key][0]) & ((end - 1) <= (elemDict[key][1] + 1)):
					sigDict[sigDictkey][key] += signal
	bedgraph.close()

#format signal into lists to make dataframe
AbList = []
spacerList = []
conditionList= []
elemList = []
sigList = []
for sample in sigDict:
	Abs, spacer, condition = sample.split('_')
	for element in sigDict[sample]:
		AbList.append(Abs)
		spacerList.append(spacer)
		conditionList.append(condition)
		elemList.append(element)
		sigList.append(sigDict[sample][element])
#construct dataframe
df = pd.DataFrame({'Antibody':AbList, 'Spacer':spacerList, 'Condition':conditionList,
				   'Element':elemList, 'Signal':sigList})
#define order of elements and conditions
elemOrder = ['PuroR', 'TRE', 'mCit', '5kb', 'pRSV', 'mCh']
condOrder = ['nodox', 'dox']
#reorder and sort dataframe
df.Element = pd.Categorical(df.Element, categories=elemOrder, ordered=True)
df.Condition = pd.Categorical(df.Condition, categories=condOrder, ordered=True)
df = df.sort_values(by=['Condition', 'Spacer', 'Antibody', 'Element'])

#calculate length of element and signal/kb
elemStart = {}
elemEnd = {}
for key in elemDict:
	elemStart[key] = elemDict[key][0]
	elemEnd[key] = elemDict[key][1]
df['Start'] = df['Element'].map(elemStart)
df['End'] = df['Element'].map(elemEnd)
df['Length'] = df['End'].astype(int) - df['Start'].astype(int)
df['Signal per kb'] = 1000*df['Signal']/df['Length']

#save dataframe for plotting signal 
df.to_csv('bedgraphAUC.csv', index=False)
print(df)

#plot signal for all conditions 
pal = ['#999999', '#333333']
# g = sns.FacetGrid(data=df, col='Antibody', hue='Condition', aspect=3, height=5)
g = sns.catplot(data=df, x='Element', y='Signal per kb', hue='Condition', col='Antibody', kind='bar', height=4, aspect=1, palette=pal)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.tight_layout()
plt.savefig('bedgraphAUC_plot_col-antibody_sigperkb_5kb.png')
plt.close()