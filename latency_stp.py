from numpy import random
from brian2 import *
import resource
from stp_syn import *

###to supress logger warnings
get_logger()
logger.suppress_hierarchy('brian2')

start_scope()

prefs.codegen.target = 'numpy'


seed=742360
np.random.RandomState(seed)
np.random.seed(seed)

loc = ['p']

#define plasticity locus-------------
if loc == ['p']:
    eqwss = eqws
    eqpotss = eqpotP#for changing p
    eqdepss = eqdepP#for changing p
    valmin = 0.01
    valmax = 1.
elif loc == ['q']: 
    eqwss = eqwsq
    eqpotss = eqpotq #for changing q
    eqdepss = eqdepq #for changing q
    valmin = 0.01
    valmax = 1.
elif loc == ['pq']:
    eqwss = eqws
    eqpotss = eqpotB #for changing both
    eqdepss = eqdepB #for changing both
    valmin = np.sqrt(0.15)
    valmax = np.sqrt(0.5)

#simulation parameters------------------
nruns = 51

init=100
plasticity_on = 0
stim=25 
period=350
freq=100*Hz

timestep = 0.1*ms

#phase reference
T=np.zeros(nruns)
for k in range(nruns):
    T[k]=k*(period+stim)
cicle = np.zeros(nruns)

#output neuron parameters
tauv = 20*ms
taug = 5*ms
E_v = -74*mV

#synaptic parameters
tau_stdp = 20*ms
initialvalue = 0.5
qmax = 0.02
Ap = 0.005
Ad = -Ap*1.05

#stochastic STP
nv = 2 #min
tau_r=200*ms
tau_p=50*ms

#output-----------------------------
eqs = '''
dv/dt = (-(ge+1)*v+E_v)/tauv : volt  
dge/dt = -ge/taug : 1
'''

P = NeuronGroup(1, eqs, threshold='v>-54*mV', reset='v=-74*mV', refractory=1*ms, order=2, method='euler', dt=timestep)

P.v = -74*mV
P.ge = 0
#-----------------------------------


#Input------------------------------
N=1000
I = NeuronGroup(N,'rates : Hz', threshold='rand()<rates*dt', order=0, method='euler', dt=timestep)


#Synapses------------------

S = Synapses(I, P, eqwss, on_pre={'pre_t': eqpotss}, on_post=eqdepss, order=1, method = 'euler', dt=timestep)
#more contacts:
#S1 = Synapses(I, P, eqwss, on_pre={'pre_t': eqpotss}, on_post=eqdepss, order=1, method = 'euler', dt=timestep)
#S2 = Synapses(I, P, eqwss, on_pre={'pre_t': eqpotss}, on_post=eqdepss, order=1, method = 'euler', dt=timestep)
#S3 = Synapses(I, P, eqwss, on_pre={'pre_t': eqpotss}, on_post=eqdepss, order=1, method = 'euler', dt=timestep)
#S4 = Synapses(I, P, eqwss, on_pre={'pre_t': eqpotss}, on_post=eqdepss, order=1, method = 'euler', dt=timestep)

S.connect()

S.U=initialvalue
S.p = S.U
S.R = 1

#define transmission delays, centered around 100 ms to avoid negative delays
npos = 0
npre = 0
Gw = np.zeros(N, dtype = int)
for k in range(len(S.i)):
    a = np.random.normal(100,15)
    S.pre_t.delay[k]=a*ms
    if(a>100):
        Gw[k] = 1
        npos = npos +1
    else:
        Gw[k] = 0
        npre = npre +1


for nrun in range(nruns):
    P.v = E_v
    P.ge = 0

    if nrun > plasticity_on: #to look at interval before plasticity
        Ap = 0.05
        Ad = -Ap*1.05
    else:
        Ad=0.
        Ap=0.

    Mf = SpikeMonitor(P,record=True)
    I.rates=freq
    run(stim*ms)

    I.rates=0*Hz
    run(period*ms)
    print(Mf.num_spikes)#S.Apre,S.Apost)
    if Mf.num_spikes>0:
        cicle[nrun]=Mf.t[0]/ms
        if Mf.num_spikes>1:
        #    for k in range(1,Ms.num_spikes):
             with open('output_lat_stp_duration','a') as nfile:
                nfile.write("%f %d \n" % ((Mf.t[Mf.num_spikes-1]-Mf.t[0])/ms,nrun))
             with open('output_lat_stp_frequency','a') as nfile:
                nfile.write("%f %d \n" % (Mf.num_spikes*ms/(Mf.t[Mf.num_spikes-1]-Mf.t[0])*1000,nrun))
        else:
             with open('output_lat_stp_duration','a') as nfile:
                 nfile.write("0 %d \n" % (nrun))
             with open('output_lat_stp_frequency','a') as nfile:
                 nfile.write("%f %d \n" % (1./36,nrun))
    else:
        #if there is no spike the latency is registered as 0
        cicle[nrun]=T[nrun]
        with open('output_lat_stp_duration','a') as nfile:
            nfile.write("0 %d \n" % (nrun))
        with open('output_lat_stp_frequency','a') as nfile:
            nfile.write("0 %d \n" % (nrun))
    del Mf

    ###to print average weights
    acum = np.zeros(2)
    for k in range(N):
        acum[Gw[k]] =acum[Gw[k]]+S.U[k]
    print(acum[0]/npre,acum[1]/npos,nrun)
    #with open('output_lat_sto_w','a') as nfile:
    #    nfile.write("%f %f %d\n" % (acum[0]/npre,acum[1]/npos,nrun))

###to print latency
with open('output_lat_stp','a') as nfile:
    for k in range(nruns):
        nfile.write("%f %d\n" % (cicle[k]-T[k],k))





























