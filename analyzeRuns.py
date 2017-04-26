# import necessary python packages
import random as rand
import math as math
import numpy as numpy
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats, integrate
import matplotlib.pyplot as plt
import itertools as it
import pickle as pickle

# make a figure comparing distance of closes ligand at protein binding over Lfactor values at a particular ligand number
def rbyL(dataDict,Lfactors,ligandNum):
	results=[]
	# loop over all runs with a single ligand
	for Lfactor in Lfactors: 
		# get data
		tempResult=[]
		proteinStates=dataDict[ligandNum+Lfactor][4]
		proteinStates=filter(lambda a: len(a)>0, proteinStates)
		# iterate over trials, find min ligand distance when protein was activated
		for run in proteinStates:
			runTemp=list(run)
			runTemp.pop(0) # throw out protein state
			Rs=[]
			for ligand in runTemp:
				Rs.append(numpy.linalg.norm(ligand))
			tempResult.append(float(min(Rs)))
		results.append(tempResult)
	sns.set_style("ticks")

	sns.set(font_scale=2) 
	ax=plt.boxplot(results)
	sns.set_palette(sns.light_palette("purple/blue", input="xkcd", reverse=True))
	sns.set_style("white")
	#sns.set_style("whitegrid")
	#plt.tight_layout()
	plt.xticks(range(len(Lfactors)),Lfactors)
	plt.ylabel('Closest ligand distance')
	plt.xlabel('Lfactor')
	plt.title('R at protein activation- '+str(ligandNum)+' ligands')
	plt.tight_layout()
	plt.savefig('RvsLfactor_ligand_'+str(ligandNum)+'.png', bbox_inches='tight')
	#plt.show()


def IFfracPlot(dataDict,Lfactors,ligands):
	results=[]
	# loop over all ligand numbers
	for ligand in ligands:
		tempResult=[]
		for Lfactor in Lfactors: # loop over all Lfactors 
			# get data
			IF=sum(dataDict[ligand+Lfactor][1]) # how many IF
			CS=sum(dataDict[ligand+Lfactor][2]) # how many CS
			# find IF fraction
			if CS>0:
				ratio=1.*IF/(IF+CS)
			else:
				ratio=1
			tempResult.append(ratio)
		results.append(tempResult)
	#plot IF fractions
	sns.set_style("ticks")
	#df = pd.DataFrame.from_items(results, ligands)
	sns.set_style("whitegrid")
	sns.set(font_scale=2) 
	sns.set_palette(sns.light_palette("purple/blue", input="xkcd", reverse=True))

	y_pos = 3*np.arange(len(Lfactors))
	width = 0.45       # the width of the bars
	fig, ax = plt.subplots()
	colors=['b','c','g','y','r']
	# make a bar plot with legend
	recter=[]
	for i in range(len(results)):
		recter.append(ax.bar(y_pos+i*width, results[i], width, color=colors[i]))
	rects=[rect[0] for rect in recter]
	ax.legend(tuple(rects), tuple(ligands),bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	#plt.bar(y_pos, results[0], align='center', alpha=0.5, capsize=20)
	plt.xticks(y_pos+2*width, Lfactors)
	plt.ylabel('Proportion IF')
	plt.xlabel('L_factor')
	plt.title('Fraction IF by L_factor and Ligand Number')
	plt.tight_layout()
	sns.despine()
	plt.tight_layout()
	plt.savefig('IFfracPlot.png', bbox_inches='tight')
	#plt.show()


if __name__ == '__main__':

	Lfactors=[.001, .0025, .005, .01, .025 ] # Lfactors that were tested
	ligands=[ 1, 2, 5, 10, 25] # ligands that were tested

	dataDict={}
	# read in data from runs into a dictionary that specifies starting state
	for ligand in ligands:
		for Lfactor in Lfactors:
			data=pickle.Unpickler(open( "Lfactor_"+str(Lfactor)+"_ligand_"+str(ligand)+"_output", "rb" )).load()
			dataDict[ligand+Lfactor]=data
	# the dict is saving a list with the following variables: [reacts- did the reaction proceed?, IFs- via IF?, CSs- Visa CS?, proteinActivateStates- where were ligands when receptor activated, ligandReactStates- where were ligands when the reaction happened?]
	IFfracPlot(dataDict,Lfactors,ligands)
	# generage r vs Lfactor plots
	for ligand in ligands:
		rbyL(dataDict,Lfactors,ligand)