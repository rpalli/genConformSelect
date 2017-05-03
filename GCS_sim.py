# import necessary python packages
import random as rand # random number generator
import math as math # math utility
import numpy as numpy # fast numeric calculations
import pandas as pd # data management wrapper
import seaborn as sns # plotting software that uses above pandas data wrappers
from scipy import stats # statistical utility
import matplotlib.pyplot as plt # finer tuning of plots
import itertools as it # utility for simpler loops
import pickle as pickle # package to save data in easy to re-read python format

# Energy function extracted from Grieves and Zhou for when protein is active
def U(r, R1, L, U0):
	return -U0*(math.tanh((r-R1)/L)/2)

# class for calculation of the Boltzmann factors for each individual ligand... initialize with U0, L, and R... 
class ligPotential:
	def __init__(self, R, U0, Lfactor):
		self.U0=U0
		self.L=Lfactor*R  # set sharpness of switch
		self.R=R
		self.R1=1.1*R # set the outer radius of the protein shell as in Grieves and Zhou
	# Boltzmann factor calculation if protein is active
	def BFcalc(self, r):
		return U(r, self.R1, self.L, self.U0)
	# return Bolzmann factor that is correct
	def BF(self, r,active):
		if active:
			return self.BFcalc(r)
		else:
			return 0

# class for calculation of the Boltzmann factors for the receptor... initialize with U0, L, and R... 
class recPotential:
	def __init__(self, R, U0, Lfactor):
		self.U0=U0
		self.L=Lfactor*R  # set sharpness of switch
		self.R=R
		self.R1=1.1*R # set the outer radius of the protein shell as in Grieves and Zhou
	# Boltzmann factor calculation for the receptor at a particular state
	def BF(self,states):
		return numpy.sum([U(numpy.linalg.norm(states[i]), self.R1, self.L, self.U0) for i in range(1,len(states)) ])

# generate a set of randomly initialized initial states for numLigands ligands in dims dimensions with distace between R and R2. 
def getInitialStates(numLigands, dims, R2, R):
	# we do this by generating states between -R2 and R2 in each dimension then throwing them out if they are too far away. 
	states=2.*R2*numpy.random.random_sample((numLigands, dims))-1.*R2
	for i in range(len(states)):
		while not R<numpy.linalg.norm(states[i])<R2:
			states[i]=2.*R2*numpy.random.random_sample((dims))-1.*R2
	return [False]+list(states) # return a "state" which includes locations of all the ligands and the state of the protein. We always start the protein off

# this is function to perform a ligand move
def updateLigState(state, potential, R2, delta):
	# we pick a random ligand, generate a move, calculate the energy, and see if the move should be accepted
	# pick a ligand, x, to attempt to move
	x=1+math.floor(1.*(len(state)-1)*rand.random())
	possible=False # keep generating moves until one that is possible is generated
	while(not possible):
		deltas=delta*(.5-numpy.random.random_sample(tuple([len(state[x])])))
		newState=numpy.add(state[x],deltas)
		# print(newState)
		r1=numpy.linalg.norm(state[x])
		r=numpy.linalg.norm(newState)
		if potential.R<r<R2:
			possible=True
	# find energy of generated possible move
	energy=potential.BF(float(r),bool(state[0]))-potential.BF(float(r1),bool(state[0]))
	# decide whether to accept move
	if rand.random()< math.exp(energy):
		state[x]=newState
	return state

# function to perform MC and return trajectory or important characteristics of the simulation together
# our formalism lets R=1
def MCSim(steps, numLigands, ligPotential, pActivateStep, R2, dims, delta, omegaNeg, gamma):
	
	# steps gives the number of MC steps
	# numLigands gives the number of ligands
	# ligPotential gives the ligand potential
	# pActivateStep gives the probability of a move to attempt to activate or inactivate the ligand
	# R2 gives the radius that ligands are allowed to diffuse in
	# dims gives the number of dimensions
	# delta gives the diffision constant i.e. how much a ligand is allowed to move on each step
	# omegaNeg is the rate of reversion of receptor to the inactive state at each MC step
	# gamma is the rate constant for the reaction after the ligand is in an active receptor

	# set up receptor potential, initial states of outputlists, system
	rec=recPotential(ligPotential.R, ligPotential.U0, ligPotential.L)
	react=False 
	trajectory=[] # stores successive list of states
	truthList=[]
	state=getInitialStates(numLigands, dims, R2, ligPotential.R) # set up initial ligand positions
	IF=False
	CS=False
	proteinActivateState=[] # stores states at which the protein was activated before final reaction
	ligandReactState=[] # stores locations of ligands when the protein became activated before final reaction
	
	trajectory.append(list(state)) # add initial state to the trajectory

	# iterate over n MC steps
	for step in range(steps): # iterate over MC steps
		x=rand.random() #generate random number to pick what type of step to do
		if x<pActivateStep: # check if step to be performed is protein activation, if so, see if receptor should be set to active state
			y=rand.random()
			# if the protein is active, then try to inactivate
			if state[0]:
				if y<omegaNeg:
					state[0]=False
					proteinActivateState=[]
					CS=False
					IF=False
			else:
				# if protein is inactive, try to activate. 
				if y< omegaNeg*math.exp(rec.BF(state)):
					state[0]=True
					# update states of various flags to indicate whether IF or CS mechanism, ligand locations as of this activation
					proteinActivateState=list(state)
					for k in range(1,len(state)):
						if numpy.linalg.norm(state[k])< ligPotential.R1:
							IF=True
					if not IF:
						CS=True
		else: # if not a receptor step, move a ligand
			state=updateLigState(state, ligPotential, R2, delta)

		trajectory.append(list(state)) # add step to the trajectory
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
	return [react,IF, CS, proteinActivateState, ligandReactState] # return important var: did we react? was it IF? was it CS? where were ligands when the protein activated? where were ligands when the protein reacted? 
	#return trajectory # switch out return statements if we want to see the trajectory

# function to run an experiment by running many MC sims and collecting results
def runSim(numLigands, Lfactor,dims, omegaNeg):
	# set up default parameters
	steps=1000000
	trials=1000
	R=1.
	U0=math.log(100.)
	stepsBeforeCheck=100
	R2=11.
	delta=.25
	pActivateStep=1/(stepsBeforeCheck*numLigands)
	gamma=10

	# set up ligand potential, output arrays
	ligPot=ligPotential( R, U0, Lfactor)
	outStates=[]
	IFs=[]
	CSs=[]
	reacts=[]
	proteinActivateStates=[]
	ligandReactStates=[]
	# run MC sim over number of trials, print results, compile them, then print summary of all results and return result arrays
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

# make a density plot of the trajectory
# requires changing the return statement to trajectory in MCSim before running
def testLigMoves():
	# set up parameters for ligand move test
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
	
	# do 10 trials of 100000 steps
	outStates=[]
	for k in range(trials):
		outStates.extend(MCSim(steps, numLigands, ligPot, pActivateStep, R2, dims, delta, omegaNeg))
	# filter out just the protein off states
	rawX=filter(lambda a: a[0]==False, outStates)
	rawX=[x[1] for x in rawX]
	filteredX=[rawX[i] for i in range(0,len(rawX),10)]
	# plot and save
	labels = ['x', 'y']
	df = pd.DataFrame.from_records(filteredX, columns=labels)
	sns.jointplot( x='x',y='y',data=df, kind="kde")
	plt.savefig('off'+str(delta*100)+'.png', bbox_inches='tight')
	# filter out just the protein on states
	rawX=filter(lambda a: a[0]==True, outStates)
	rawX=[x[1] for x in rawX]
	filteredX=[rawX[i] for i in range(0,len(rawX),10)]
	# plot and save
	labels = ['x', 'y']
	df = pd.DataFrame.from_records(filteredX, columns=labels)
	sns.jointplot( x='x',y='y',data=df, kind="kde")
	plt.savefig('on'+str(delta*100)+'.png', bbox_inches='tight')

# function to run runSim with different numbers of ligands and lfactors for D=.5 and dims=2...
# it will save each parameter set seperately in a picklefile with four lists. 
def LfactorTest():
	for ligands in [2, 5]:
		print("ligand")
		print(ligands)
		for Lfactor in [1, .5, .1 ]:
			print("Lfactor")
			print(Lfactor)
			output=runSim(ligands, Lfactor,2, .5)
			pickle.dump(output , open( "Lfactor_"+str(Lfactor)+"_ligand_"+str(ligands)+"_output", "wb" ) )
LfactorTest()