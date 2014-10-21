import random
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#attemtps to implement a gillespie direct method on an SEIR model
#with metapopulations
#
#TODO 
#Also: make sure the population selection + "popprob" probabilities are doing
#what they're supposed to be doing
#
#Add import event

#first define all parameters
mu = 1.00/70                             #death rate
gamma = 365.00/10                        #recovery rate
R_0 = 15.00                              #basic reproductive ratio
t = 0.0                                  #initial time (years)
tmax = 1                                 # end time
t_append = 0.0                           #when to save data
nu = mu                                  #birth rate
alpha = 365.00/10                        #latent period
beta = R_0*(gamma + mu)*(alpha+mu)/alpha #transmission rate
metapopsize = 5                          #number of populations in metapopulation

class Population(object):
    """A population has all the variables N, S, E, I, R"""
    def __init__(self, index, N=10000, S=1500, E=0, I=10, R=10000-1500-10):
        self.index = index
        self.N = N
        self.S = S
        self.E = E
        self.I = I
        self.R = R
    
    def birth(self):
        self.N += 1
	self.S += 1

    def sdeath(self):
        self.N -= 1
	self.S -= 1

    def expose(self):
        self.S -= 1
	self.E += 1

    def edeath(self):
        self.N -= 1
	self.E -= 1

    def symptom(self):
	self.E -= 1
	self.I += 1

    def ideath(self):
	self.N -= 1
        self.I -= 1

    def recover(self):
	self.I -= 1
	self.R += 1

    def rdeath(self):
	self.N -= 1
	self.R -= 1

def pbirth(n):
    p_birth = nu*n
    return p_birth

def psdeath(s):
    p_sdeath = mu*s
    return p_sdeath

def pexpose(s, i, n):
    p_infect = beta*s*i/(n)    #warning, adding 1 to n to avoid /0 error (never mind)
    return p_infect

def pedeath(e):
    p_edeath = mu*e
    return p_edeath

def psymptom(e):
    p_symptom = alpha*e
    return p_symptom

def pideath(i):
    p_ideath = mu*i
    return p_ideath

def precover(i):
    p_recover = gamma*i
    return p_recover

def prdeath(r):
    p_rdeath = mu*r
    return p_rdeath

pop = []
psum = []
# probably make this into a function or something so I can clean up the big
# loop
for x in range(0, metapopsize):
    pop.append(Population(x, 10000, 1500, 0, 10, 8490))
    psum.append((pbirth(pop[x].N) + psdeath(pop[x].S) +
              pexpose(pop[x].S, pop[x].I, pop[x].N) +
              pedeath(pop[x].E) + psymptom(pop[x].E) + pideath(pop[x].I) + 
	      precover(pop[x].I) + prdeath(pop[x].R)))
	
pmetasum = sum(psum)
popprob = []
for x in range(0, metapopsize): #this needs to go in the loop :(
    popprob.append(psum[x]/pmetasum)

#print precover(pop[2].I)
#print pexpose(pop[3].S, pop[3].I, pop[3].N)

dataArray = np.array([0, 10000, 1500, 0, 10, 8490]) #default values for now :(

while t < tmax:
    randomtime, randomevent, randompopulation = random.random(), random.random(), random.random()
    dt = -math.log(randomtime)/pmetasum
    popchoice = 0
    popselect = 0
    #make sure (prove?) the following selection of which population to perform
    #the next event is actually choosing properly (it's not)
    while (popchoice < randompopulation) and (popselect<(metapopsize)):
        popchoice = popchoice + popprob[popselect]
	popselect += 1
    popselect -=1
    #print popselect
    if randomevent < pbirth(pop[popselect].N)/psum[popselect]:
        pop[popselect].birth()
    elif randomevent < (pbirth(pop[popselect].N)+psdeath(pop[popselect].S))/psum[popselect]:
	pop[popselect].sdeath()
    elif randomevent < ((pbirth(pop[popselect].N)+psdeath(pop[popselect].S)
	      + pexpose(pop[popselect].S, pop[popselect].I, pop[popselect].N))
              /psum[popselect]):
	pop[popselect].expose()
    elif randomevent < ((pbirth(pop[popselect].N)+psdeath(pop[popselect].S)
	      + pexpose(pop[popselect].S, pop[popselect].I, pop[popselect].N)+
	      pedeath(pop[popselect].E))/psum[popselect]):
	pop[popselect].edeath()
    elif randomevent < ((pbirth(pop[popselect].N)+psdeath(pop[popselect].S)
	      + pexpose(pop[popselect].S, pop[popselect].I, pop[popselect].N)+
	      pedeath(pop[popselect].E)+psymptom(pop[popselect].E))/psum[popselect]):
	pop[popselect].symptom()
    elif randomevent < ((pbirth(pop[popselect].N)+psdeath(pop[popselect].S)
	      + pexpose(pop[popselect].S, pop[popselect].I, pop[popselect].N)+
	      pedeath(pop[popselect].E)+psymptom(pop[popselect].E)+
	      pideath(pop[popselect].I))/psum[popselect]):
	pop[popselect].ideath()
    elif randomevent < ((pbirth(pop[popselect].N)+psdeath(pop[popselect].S)
	      + pexpose(pop[popselect].S, pop[popselect].I, pop[popselect].N)+
	      pedeath(pop[popselect].E)+psymptom(pop[popselect].E)+
	      pideath(pop[popselect].I)+precover(pop[popselect].I))/psum[popselect]):
	pop[popselect].recover()
    else:
	pop[popselect].rdeath()
	#now update all the probabilities qq
    psum[popselect] = ((pbirth(pop[popselect].N) + psdeath(pop[popselect].S) +
              pexpose(pop[popselect].S, pop[popselect].I, pop[popselect].N) +
              pedeath(pop[popselect].E) + psymptom(pop[popselect].E) +
	      pideath(pop[popselect].I) + precover(pop[popselect].I) +
	      prdeath(pop[popselect].R)))
    pmetasum = sum(psum)
    popprob[popselect]=(psum[popselect]/pmetasum)
    if t > t_append:
	t_append = t_append + 1.0/365.0
	t += dt
	#dataArray stuff save data for whatever
	dataArray = np.vstack((dataArray,[t, pop[1].N, pop[1].S, pop[1].E, pop[1].I,
	                     pop[1].R]))
	#print popselect
    else:
        t += dt
	
plt.figure()
plt.ioff()
plt.plot(dataArray[:,0],dataArray[:,2])
plt.hold(True)
plt.plot(dataArray[:,0],dataArray[:,3])
plt.plot(dataArray[:,0],dataArray[:,4])
#plt.plot(dataArray[:,0],dataArray[:,1])
plt.show()
#print psum
#print pmetasum

#print popprob
#print sum(popprob)
