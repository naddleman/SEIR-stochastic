import random
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#implements the Gillespie Direct Method algorithm for a basic SIR model
#
# Now uses pandas to store the data as a Data Frame
# Can output the data from one simulation to a .csv (comma separated values) file
# 2014-10-12 rev 4
# implementing E class and metapopulations DONE!


def runSim():
    """simulation function"""
    #print "Test2"
    tmax = 1                            #simulation time (years)
    N_init = 10000                      #initial population size
    S_init = 1500                       #initial susceptible population
    E_init = 0                          #initial exposed population
    I_init = 10                         #initial infected population
    mu = 1.00/70                        #death rate
    gamma = 365.00/10                   #recovery rate
    R_0 = 15.00                         #R_0, basic reproductive ratio
    nsims = 10                          #number of simulations to run
    R_init = N_init - S_init - I_init - E_init   # initial recovered populataion
    nu = mu                             # birth rate
    alpha = 365.00/10                   #Latent period
    beta = R_0*(gamma + mu)             #transmission rate
    #n                                  #which simulation, needed when looping
    time_list = []                      #list of times for each event
    S_list = []
    E_list = []
    I_list = []
    R_list = []
    N_list = []
    t = 0                               #initial time
    j = 1
    S, E, I, R = S_init, E_init, I_init, R_init
    N = N_init
    p_birth = nu*N                          #probability parameters (initial values)
    p_sdeath = mu*S
    p_infection = beta*S*I/N
    p_edeath = mu*E
    p_symptom = alpha*E
    p_ideath = mu*I
    p_recovery = gamma*I
    p_rdeath = mu*R
    while t < tmax:
        A = p_birth+p_sdeath+p_infection+p_edeath+p_symptom+p_ideath+p_recovery+p_rdeath
	#sum of the transmission rates
        y_1, y_2 = random.random(), random.random() # generate 2 random numbers in [0,1)
	dt = -math.log(y_1)/A           #time step: time until next event
	#simulate each evenet individually, update SEIRN and parameters
	if y_2 < p_birth/A:                                         #Birth
            S = S+1                    
            N = N+1
            p_birth = mu*N
            p_sdeath = mu*S
        elif y_2 < (p_birth+p_sdeath)/A:                            #S death
            S = S-1
            N = N-1
            p_birth = mu*N
            p_sdeath = mu*S
        elif y_2 < (p_birth+p_sdeath+p_infection)/A:                #infection
            S = S-1
            E = E+1
            p_sdeath = mu*S
            p_edeath = mu*E
	    p_symptom = alpha*E
        elif y_2 < (p_birth+p_sdeath+p_infection+p_edeath)/A:       #E death
            N = N-1
            E = E-1
	    p_birth = mu*N
            p_edeath = mu*E
            p_symptom = alpha*E
        elif y_2 < (p_birth+p_sdeath+p_infection+p_edeath+p_symptom)/A: #onset
            E = E-1
	    I = I+1
	    p_edeath = mu*E
	    p_symptom = alpha*E
	    p_ideath = mu*I
        elif y_2 < (p_birth+p_sdeath+p_infection+p_edeath+p_symptom+p_ideath)/A: #ideath  
            I = I-1
	    N = N-1
	    p_birth = mu*N
	    p_ideath = mu*I
	    p_recovery = gamma*I
        elif y_2 < (p_birth+p_sdeath+p_infection+p_edeath+p_symptom+p_ideath+p_recovery)/A:   #recovery
            R = R+1
            I = I-1
            p_ideath = mu*I
            p_recovery = gamma*I
            p_rdeath = mu*R
        else:                                                           #R death 
	    R = R-1
	    N = N-1
	    p_rdeath = mu*R
        p_infection=beta*S*I/N
        t = t+dt
        j =  j+1
        time_list.append(t)
        S_list.append(S)
	E_list.append(E)
        I_list.append(I)
        R_list.append(R)
        N_list.append(N)
    #create a DataFrame from S,I,R,N lists which can be graphed, or saved as .csv
    #TODO: only save values for variables once per day instead of each dt
    dfaa = {'Time': time_list,'Susceptibles': S_list,'Exposed': E_list, 'Infected': I_list, 'Recovered': R_list, 'Population': N_list} 
    dfa = pd.DataFrame(data=dfaa)
    dfa.replace(',','.')
    dfa = dfa.replace('-','NaN')
    #print(dfa)
    #dfa.to_csv('test.csv')
    plt.figure() 
    plt.ioff()
    dfa.plot(x='Time',y='Susceptibles')
    plt.hold(True)
    dfa.plot(x='Time',y='Exposed')
    dfa.plot(x='Time',y='Infected')
    plt.show()

runSim()
#print(dfa) 
#plt.ioff()  
#plt.plot(time_list,S_list)
#plt.hold(True)
#plt.plot(time_list,I_list,'r')
#plt.plot(time_list,R_list,'g')
#plt.show()
