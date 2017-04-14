
# import stuff I need
import random as rand
import math as math
import numpy as np
kb=.0019872041 # define kb so it can be used easily later

#define v1 and v2 as functions
def V1(x):
	return .5*(x**4-3*x**2+.5*x)
def V2(x):
	return V1(x)+.2*math.cos(20*x)
class potential:
	# class for calculation of the Boltzmann factors... using a given potential. 
	# if statements below map functions to the potential so potential.potentialName gives the name
	# and potential.energy(x) returns the actual value of the potential at x
	# and potential.BF(x, T) returns the boltzmann factor for x at T given the potential initialized
	def __init__(self,potentialName):
		self.potentialName=potentialName
		# identify the energy function with the correct v1 or v2 function from above
		if potentialName=='v1':
			self.energy=V1
		if potentialName=='v2':
			self.energy=V2
	# function of initialized potential that returns Boltzman factor
	def BF(self,x,T):
		return self.energy(x)/(kb*T)

# function to run MC sim and return P(A) 
def MCSim(steps, potentialName, temp):
	x=rand.gauss(0,1) #set initial value
	utility=potential(potentialName) # initialize potential class with correct potential
	currentBF=utility.BF(x, temp)
	states=[] #list of states taken on by the simulation (i.e. the trajectory)
	for i in range(steps): # iterate over number of steps
		y=x+rand.gauss(0,1) #generate random gaussian step from current x value
		newBF= utility.BF(y,temp) #calc BF for step
		if newBF-currentBF<0:
			x=y
			currentBF=newBF
		elif rand.random() < math.exp(-newBF+currentBF): # accept step with probability prop to BF difference
			x=y
			currentBF=newBF
		states.append(x)
	return sum([1 for state in states if state<0])/(1.*len(states))

#function to estimate avgs and errors by repeat calls to MCSim
def estimateMCerror(steps, potentialName, temp, trials):
	PA=[]# storage var for estimates
	for i in range(trials): #iterate over num trials
		PA.append(MCSim(steps,potentialName,temp)) # run a trial
	mean=1.*sum(PA)/trials # find mean
	variance=sum([(mean-sample)**2 for sample in PA])/len(PA) #find variance
	return mean, variance


# numerically solve p
def numericalIntegrator(potentialName, temp):
	positiveSide=0
	negativeSide=0
	utility=potential(potentialName)
	for i in np.arange(.0001,100,.00001):
		positiveSide+= math.exp(-utility.BF(i,temp)) #add positive numbers by BF
		negativeSide+= math.exp(-utility.BF(-i,temp)) # add negative numbers by BF
	return negativeSide/(negativeSide+positiveSide)

print(estimateMCerror(10**5,'v1',300,10))
print(numericalIntegrator('v1',300))
print(estimateMCerror(10**6,'v1',300,10))
print(estimateMCerror(10**5,'v1',400,10))
print(estimateMCerror(10**6,'v1',400,10))
print(estimateMCerror(10**6,'v2',300,10))

print('extra credit section')



#extra credit section. 
def extraCreditFunction(steps, potentialName, temp):
	
	#run one trajectory

	x=rand.gauss(0,1) #set initial value
	utility=potential(potentialName) # initialize potential class with correct potential
	currentBF=utility.BF(x, temp)
	states=[] #list of states taken on by the simulation (i.e. the trajectory)
	for i in range(steps): # iterate over number of steps
		y=x+rand.gauss(0,1) #generate random gaussian step from current x value
		newBF= utility.BF(y,temp) #calc BF for step
		if newBF-currentBF<0:
			x=y
			currentBF=newBF
		elif rand.random() < math.exp(-newBF+currentBF): # accept step with probability prop to BF difference
			x=y
			currentBF=newBF
		states.append(x)
	cs=[]

	states= [1 if state<0 else 0 for state in states] 
	# time correlation approach

	#find average and variance
	average=sum(states)/(1.*len(states)) 
	variance=sum([(state-average)**2 for state in states])/(1.*len(states))
	# go through each step length to find correlation coefficient at each
	for j in range(steps):
		c=0
		# start to calculate c by stepping over trajectory
		for i in range(0,len(states)-j):
			c+=(states[i]-average)*(states[i+j]-average)
		cs.append(c/((steps-j)*variance)) #  normalize c so it can be simply added to get tau
	# find first negative number, where the long unstable tail of c starts
	firstneg=[index for index in range(len(cs)) if cs[index]<0].pop(0)
	#print(firstneg) #print the first negative
	tau=sum([cs[num] for num in range(0,firstneg)])# integrate cs to get tau by going up to when the long unstable tail starts
	print(tau) 
	print(variance*tau/steps)


	# block averaging approach

	variances=[] # store stdevs calculated as we vary the size of block
	for n in range(1,(steps-1)/10): # loop over size of block
		newblocks=[]
		#divide into blocks
		for numblock in range(int(math.floor(steps/n))):
			newblocks.append(states[(numblock*n):(numblock*(n+1))])
		# calculate stdev by blocksize
		blockavgs=[]
		# calculate average over each block 
		for block in newblocks:
			if sum(block)>0:
				blockavgs.append(sum(block)/(1.*len(block)))
			else:
				blockavgs.append(0)
		# calculate average of blocks of size n
		nAvgs= sum(blockavgs)/len(blockavgs)
		# append stdev of blocks of size n
		variances.append(sum([(state-nAvgs)**2 for state in blockavgs])/(1.*len(blockavgs)**2))
	print(variances)
print('10,000 steps')
extraCreditFunction(10000, 'v1', 300)
print('1,000 steps')
extraCreditFunction(1000, 'v1', 300)