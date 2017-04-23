
# import necessary python packages
import random as rand
import math as math
import numpy as numpy
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats, integrate
import matplotlib.pyplot as plt



# Energy function extracted from Grieves and Zhou
def U(r, R1, L, U0):
	return -U0*(tanh((r-R1)/L)/2)

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
	def BF(states):
		return numpy.sum([U(numpy.linalg.norm(states[i]), self.R1, self.L, self.U0) for i in range(1,len(states)) ])

def rateCalc(state):
	return 0
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
		r=numpy.linalg.norm(newState)
		if potential.R<r<R2:
			possible=True
	energy=potential.BF(float(r),bool(state[0]))
	# print(energy)
	# print(math.exp(energy))
	if rand.random()< math.exp(energy):
		state[x]=newState
	return state

# function to run MC sim and return trajectory
# our formalism lets R=1
def MCSim(steps, numLigands, ligPotential, pActivateStep, R2, dims, delta):
	rec=recPotential(ligPotential.R, ligPotential.U0, ligPotential.L)
	react=False 
	trajectory=[]
	state=getInitialStates(numLigands, dims, R2, ligPotential.R) # set up initial ligand positions
	# initRate=rateCalc(state) # set up initial rate of reaction completion given starting states
	trajectory.append(list(state)) # add initial state to the trajectory
	print(state)
	for step in range(steps): # iterate over MC steps
		x=rand.random() #generate random number to pick what type of step to do
		if x<pActivateStep: # check if step to be performed is protein activation, if so, see if receptor should be set to active state
			if state[0]:
				state[0]=False
			else:
				y=rand.random()
				if y<math.exp(recPotential.BF(state)):
					state[0]=True
		else: # if not a receptor step, move a ligand
			state=updateLigState(state, ligPotential, R2, delta)
		# finalRate=rateCalc(state) # calculate rate of reaction completion at the current state
		trajectory.append(list(state))
		x=rand.random()
		# if x < (1-math.exp(-.5*(initRate+finalRate))): # check to see if reaction happens based on initial and final rates of reaction completion
		# 	react=True
		# 	break
	return trajectory

def testLigMoves():
	steps=1000000
	numLigands=1
	R=1.
	U0=100.
	Lfactor=.005
	stepsBeforeCheck=100
	R2=11.
	dims=2
	delta=.25
	pActivateStep=1/(stepsBeforeCheck*numLigands)
	ligPot=ligPotential( R, U0, Lfactor)

	for delta in [.25,.5,1]:
		outStates=MCSim(steps, numLigands, ligPot, pActivateStep, R2, dims, delta)
		rawX=[x[1] for x in outStates if x[0]==False]
		filteredX=[rawX[i] for i in range(0,len(rawX),10)]
		labels = ['x', 'y']
		df = pd.DataFrame.from_records(filteredX, columns=labels)

		sns.jointplot( x='x',y='y',data=df, kind="kde")
		plt.savefig('off'+str(delta*100)+'.png', bbox_inches='tight')

		rawX=[x[1] if x[0]==True for x in outStates]
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

testLigMoves()