
# import necessary python packages
import random as rand
import math
import numpy as numpy

# Energy function extracted from Grieves and Zhou
def U(r, R1, L, U0):
	return -U0*(tanh((r-R1)/L)/2)

class ligPotential:
	# class for calculation of the Boltzmann factors for each individual ligand... initialize with U0, L, and R... 
	def __init__(self, R, U0, Lfactor):
		self.U0=U0
		self.L=L_factor*R  # set sharpness of switch
		self.R=R
		self.R1=1.1*R # set the outer radius of the protein shell as in Grieves and Zhou
	# Boltzmann factor calculation
	def BF(r):
		return U(r, self.R1, self.L, self.U0)

class recPotential:
	# class for calculation of the Boltzmann factors for the receptor... initialize with U0, L, and R... 
	def __init__(self, R, U0, Lfactor):
		self.U0=U0
		self.L=L_factor*R  # set sharpness of switch
		self.R=R
		self.R1=1.1*R # set the outer radius of the protein shell as in Grieves and Zhou
	# Boltzmann factor calculation
	def BF(states):
		return numpy.sum([U(numpy.linalg.norm(states[i]), self.R1, self.L, self.U0) for i in range(1,len(states)) ])

def rateCalc(state):
	return 0

def getInitialStates(numLigands, ligPotential):


def updateLigState(state, potential):



# function to run MC sim and return trajectory
# our formalism lets R=1
def MCSim(steps, numLigands, ligPotential, pActivateStep):
	rec=recPotential(ligPotential.R, ligPotential.U0, ligPotential.L)
	react=False 
	trajectory=[]
	state=getInitialStates(numLigands, ligPotential) # set up initial ligand positions
	# initRate=rateCalc(state) # set up initial rate of reaction completion given starting states
	trajectory.append(state) # add initial state to the trajectory

	for step in range(steps): # iterate over MC steps
		x=rand.random() #generate random number to pick what type of step to do
		if x<pActivateStep: # check if step to be performed is protein activation, if so, see if receptor should be set to active state
			if state[0]:
				state[0]=False
			else:
				y=rand.random()
				if y<exp(recPotential.BF(state)):
					state[0]=True
		else: # if not a receptor step, move a ligand
			state=updateLigState(state, potential)
		# finalRate=rateCalc(state) # calculate rate of reaction completion at the current state
		trajectory.append(state)
		x=rand.random()
		# if x < (1-exp(-.5*(initRate+finalRate))): # check to see if reaction happens based on initial and final rates of reaction completion
		# 	react=True
		# 	break
	return trajectory

#function to estimate avgs and errors by repeat calls to MCSim
def estimateMCerror(steps, potentialName, temp, trials):
	PA=[]# storage var for estimates
	for i in range(trials): #iterate over num trials
		PA.append(MCSim(steps,potentialName,temp)) # run a trial
	mean=1.*sum(PA)/trials # find mean
	variance=sum([(mean-sample)**2 for sample in PA])/len(PA) #find variance
	return mean, variance