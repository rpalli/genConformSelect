
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

# Energy function extracted from Grieves and Zhou
def U(r, R1, L, U0):
	return -U0*(math.tanh((r-R1)/L)/2)

class ligPotential:
	# class for calculation of the Boltzmann factors for each individual ligand... initialize with U0, L, and R... 
	def __init__(self, R, U0, Lfactor):
		self.U0=U0
		self.L=Lfactor*R  # set sharpness of switch
		self.R=R
		self.R1=1.1*R # set the outer radius of the protein shell as in Grieves and Zhou
	# Boltzmann factor calculation
	def BFcalc(self, r):
		return U(r, self.R1, self.L, self.U0)
	def BF(self, r,active):
		if active:
			return self.BFcalc(r)
		else:
			return 0

class recPotential:
	# class for calculation of the Boltzmann factors for the receptor... initialize with U0, L, and R... 
	def __init__(self, R, U0, Lfactor):
		self.U0=U0
		self.L=Lfactor*R  # set sharpness of switch
		self.R=R
		self.R1=1.1*R # set the outer radius of the protein shell as in Grieves and Zhou
	# Boltzmann factor calculation
	def BF(self,states):
		return numpy.sum([U(numpy.linalg.norm(states[i]), self.R1, self.L, self.U0) for i in range(1,len(states)) ])

def getInitialStates(numLigands, dims, R2, R):
	states=2.*R2*numpy.random.random_sample((numLigands, dims))-1.*R2
	for i in range(len(states)):
		while not R<numpy.linalg.norm(states[i])<R2:
			states[i]=2.*R2*numpy.random.random_sample((dims))-1.*R2
		# print(states[i])
		# print(numpy.linalg.norm(states[i]))
	return [False]+list(states)

def updateLigState(state, potential, R2, delta):
	x=1+math.floor(1.*(len(state)-1)*rand.random())
	possible=False
	while(not possible):
		deltas=delta*(.5-numpy.random.random_sample(tuple([len(state[x])])))
		newState=numpy.add(state[x],deltas)
		# print(newState)
		r1=numpy.linalg.norm(state[x])
		r=numpy.linalg.norm(newState)
		if potential.R<r<R2:
			possible=True
	energy=potential.BF(float(r),bool(state[0]))-potential.BF(float(r1),bool(state[0]))
	# print(energy)
	# print(math.exp(energy))
	if rand.random()< math.exp(energy):
		state[x]=newState
	return state

# function to run MC sim and return trajectory
# our formalism lets R=1
def MCSim(steps, numLigands, ligPotential, pActivateStep, R2, dims, delta, omegaNeg, gamma):
	rec=recPotential(ligPotential.R, ligPotential.U0, ligPotential.L)
	react=False 
	trajectory=[]
	truthList=[]
	state=getInitialStates(numLigands, dims, R2, ligPotential.R) # set up initial ligand positions
	# trajectory.append(list(state)) # add initial state to the trajectory
	IF=False
	CS=False
	# print(state)
	proteinActivateState=[]
	ligandReactState=[]
	for step in range(steps): # iterate over MC steps
		x=rand.random() #generate random number to pick what type of step to do
		if x<pActivateStep: # check if step to be performed is protein activation, if so, see if receptor should be set to active state
			y=rand.random()
			if state[0]:
				if y<omegaNeg:
					state[0]=False
					proteinActivateState=[]
					CS=False
					IF=False
			else:
				if y< omegaNeg*math.exp(rec.BF(state)):
					state[0]=True
					proteinActivateState=list(state)
					# print('protein active')
					for k in range(1,len(state)):
						if numpy.linalg.norm(state[k])< ligPotential.R1:
							IF=True
					if not IF:
						CS=True
		else: # if not a receptor step, move a ligand
			state=updateLigState(state, ligPotential, R2, delta)
		truth=[]
		if state[0]==True:
			truth.append(True)
		else:
			truth.append(False)
		bound=False
		for ligand in state[1:]:
			if numpy.linalg.norm(ligand)<ligPotential.R1:
				bound=True
		truth.append(bound)
		if(truth[0] and bound):
			truth.append(True)
		else:
			truth.append(False)
		# trajectory.append(list(state))
		checkReact=False
		if state[0]==True:
			for k in range(1,len(state)):
				if numpy.linalg.norm(state[k])< ligPotential.R1:
					checkReact=True
		if checkReact:
			x=rand.random()
			if x < (1-math.exp(-1.*gamma)): # check to see if reaction happens based on initial and final rates of reaction completion
				react=True
				ligandReactState=list(state)
				# print('reacted')
				break
	return [react,IF, CS, proteinActivateState, ligandReactState]
	#return trajectory

def runSim(numLigands, Lfactor,dims, omegaNeg):
	steps=1000000
	trials=1000
	R=1.
	U0=math.log(100.)
	stepsBeforeCheck=100
	R2=11.
	delta=.25
	pActivateStep=1/(stepsBeforeCheck*numLigands)
	gamma=10
	numLigands=2
	omegaNeg=.5
	Lfactor=.005
	dims=2
	
	pActivateStep=1/(stepsBeforeCheck*numLigands)
	ligPot=ligPotential( R, U0, Lfactor)

	ligPot=ligPotential( R, U0, Lfactor)
	outStates=[]
	IFs=[]
	CSs=[]
	reacts=[]
	proteinActivateStates=[]
	ligandReactStates=[]
	#MCSim(steps, numLigands, ligPot, pActivateStep, R2, dims, delta, omegaNeg)
	for k in range(trials):
		[react,IF, CS, proteinActivateState, ligandReactState]=MCSim(steps, numLigands, ligPot, pActivateStep, R2, dims, delta, omegaNeg, gamma)
		print([react,IF, CS, proteinActivateState, ligandReactState])
		reacts.append(react)
		IFs.append(IF)
		CSs.append(CS)
		proteinActivateStates.append(proteinActivateState)
		ligandReactStates.append(ligandReactState)
	print(sum(reacts))
	print(sum(CSs))
	print(sum(IFs))
	return [reacts, IFs, CSs, proteinActivateStates, ligandReactStates]

def testLigMoves():
	steps=100000
	trials=10
	numLigands=2
	R=1.
	U0=100.
	Lfactor=.005
	stepsBeforeCheck=100
	R2=11.
	dims=2
	delta=.25
	pActivateStep=1/(stepsBeforeCheck*numLigands)
	ligPot=ligPotential( R, U0, Lfactor)
	omegaNeg=.5
	
	outStates=[]
	for k in range(trials):
		outStates.extend(MCSim(steps, numLigands, ligPot, pActivateStep, R2, dims, delta, omegaNeg))
	#rawX=[x[1] if x[0]==False else 100 for x in outStates ]
	rawX=filter(lambda a: a[0]==False, outStates)
	rawX=[x[1] for x in rawX]
	filteredX=[rawX[i] for i in range(0,len(rawX),10)]
	labels = ['x', 'y']
	df = pd.DataFrame.from_records(filteredX, columns=labels)
	sns.jointplot( x='x',y='y',data=df, kind="kde")
	plt.savefig('off'+str(delta*100)+'.png', bbox_inches='tight')
	rawX=filter(lambda a: a[0]==True, outStates)
	rawX=[x[1] for x in rawX]
	filteredX=[rawX[i] for i in range(0,len(rawX),10)]
	labels = ['x', 'y']
	df = pd.DataFrame.from_records(filteredX, columns=labels)
	sns.jointplot( x='x',y='y',data=df, kind="kde")
	plt.savefig('on'+str(delta*100)+'.png', bbox_inches='tight')

#function to estimate avgs and errors by repeat calls to MCSim
def estimateMCerror(steps, potentialName, temp, trials):
	PA=[]# storage var for estimates
	for i in range(trials): #iterate over num trials
		PA.append(MCSim(steps,potentialName,temp)) # run a trial
	mean=1.*sum(PA)/trials # find mean
	variance=sum([(mean-sample)**2 for sample in PA])/len(PA) #find variance
	return mean, variance
# testLigMoves()
def LfactorTest():
	for ligands in [2, 1, 5, 10]:
		print("ligand")
		print(ligands)
		for Lfactor in [.005, .0025, .01, .001, .025]:
			print("Lfactor")
			print(Lfactor)
			output=runSim(ligands, Lfactor,2, .5)
			pickle.dump(output , open( "Lfactor_"+str(Lfactor)+"_ligand_"+str(ligands)+"_output", "wb" ) )
LfactorTest()