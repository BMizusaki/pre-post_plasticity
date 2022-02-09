import numpy as np
from brian2 import *
import resource

###to supress logger warnings
#get_logger()
#logger.suppress_hierarchy('brian2')

start_scope()

prefs.codegen.target = 'numpy'

U_on = 1.
A_on = 1.
###to block presynaptic plasticity
#A_on = 0
###to block postsynaptic plasticity 
#U_on = 0

seed=742360
np.random.RandomState(seed)
np.random.seed(seed)

#simulation parameters------------------
nruns = 51

init=100
plasticity_on = 0
spread=7
stim=35 
period=250
freq=100*Hz

timestep = 0.1*ms

#phase reference
T=np.zeros(nruns)
cicle = np.zeros(nruns)

#output neuron parameters
tauv = 20*ms
taug = 5*ms
E_v = -74*mV
C = 281 * pF
gL = 30 * nS
taum = C / gL
EL = -70.6 * mV
DeltaT = 2 * mV
vt = -50.4 * mV
tauvt = 50 * ms

#synaptic parameters
tau_stdp = 20*ms
initialvalue = 0.4
s = 0.25*gL#100. * nS
Amax = 0.5
Amin = 0
Ainit = 0.1
Umax = 1.
Umin = 0
Uinit = 0.1

tau_r=200*ms
tau_u=50*ms

#output-----------------------------
eqs = '''
        dvm/dt=(gL*(EL-vm)+gL*DeltaT*exp((vm-vt)/DeltaT)-ge*vm-x)/C : volt
        dx/dt=(c*(vm-EL)-x)/tauw : amp  
        dge/dt = -ge/taug : siemens   
        NA : 1
'''
tauw, c, b, Vr = 144 * ms, 4 * nS, 0.0805 * nA, -70.6 * mV
P = NeuronGroup(1, eqs, threshold='vm>vt', reset="vm=Vr;x+=b", refractory=1*ms, order=2, method='euler', dt=timestep)

P.vm = EL
P.x = 0
P.ge = 0
P.NA = 0
#-----------------------------------


#Input------------------------------
Nsyn=100
I = NeuronGroup(Nsyn,'rates : Hz', threshold='rand()<rates*dt', order=0, method='euler', refractory=5*ms, dt=timestep)

params = [0.1771,    0.0327,    0.1548,    0.2302,    0.0618,    0.0666];
AFBn = params[0]
tau_FBn = params[1]*1e3 * ms
AFBp = params[2]
tau_FBp = params[3]*1e3 * ms
AFFp = params[4]
tau_FFp = params[5]*1e3 * ms
dFBn = 0
dFBp = 0
dFFp = 0

#Synapses------------------

eqw='''
w : 1
U : 1
A : 1         
NA_post=AFFp*FFp*FBn : 1 (summed)
dFFp/dt=-FFp/tau_FFp : 1 (event-driven)
dFBp/dt=-FBp/tau_FBp : 1 (event-driven)
dFBn/dt=-FBn/tau_FBn : 1 (event-driven)
dR/dt=(1-R)/tau_r : 1 (clock-driven)
du/dt=(U-u)/tau_u : 1 (clock-driven) 
'''
eqdep='''
A=A+etaA*(AFFp*FFp*FBn)
A=A-etaA*0.5*NA_post/Nsyn
A=clip(A,Amin,Amax)
w=U*A
FBp+=1.
FBn+=1.
'''
eqpot='''
ge_post +=  s*A*R*u 
U=clip(U+etaU*(-AFBn*FBn*FBp + AFBp*FBp*FFp),Umin,Umax)
w=U*A
FFp+=1
R-=R*u
u+=U*(1-u)
'''

S = Synapses(I, P, eqw, on_pre={'pre_t':eqpot}, on_post=eqdep, order=1, method = 'euler', dt=timestep)
#more contacts:
#S1 = Synapses(I, P, eqwss, on_pre={'pre_t': eqpotss}, on_post=eqdepss, order=1, method = 'euler', dt=timestep)
#S2 = Synapses(I, P, eqwss, on_pre={'pre_t': eqpotss}, on_post=eqdepss, order=1, method = 'euler', dt=timestep)
#S3 = Synapses(I, P, eqwss, on_pre={'pre_t': eqpotss}, on_post=eqdepss, order=1, method = 'euler', dt=timestep)
#S4 = Synapses(I, P, eqwss, on_pre={'pre_t': eqpotss}, on_post=eqdepss, order=1, method = 'euler', dt=timestep)

S.connect()

S.U=initialvalue
S.FBp=0
S.FBn=0
S.FFp=0
S.R=1
S.A = initialvalue*Amax
S.u = S.U

#define transmission delays, centered around 100 ms to avoid negative delays
npos = 0
npre = 0
Gw = np.zeros(Nsyn, dtype = int)
for k in range(len(S.i)):
    a = np.random.normal(init,spread)
    S.pre_t.delay[k]=a*ms
    if(a>init):
        Gw[k] = 1
        npos = npos +1
    else:
        Gw[k] = 0
        npre = npre +1

for nrun in range(nruns):
    P.ge = 0
    P.NA = 0
    P.vm = EL
    P.x = 0
    S.FBp=0
    S.FBn=0
    S.FFp=0
    S.R = 1
    S.u = S.U
    if nrun > plasticity_on: #to look at interval before plasticity
        etaA = A_on*0.15
        etaU = U_on*0.15
    else:
        etaA = 0.
        etaU = 0.

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
    acum = np.zeros(4)
    for k in range(Nsyn):
        acum[Gw[k]] =acum[Gw[k]]+S.U[k]/Umax
        acum[Gw[k]+2] =acum[Gw[k]+2]+S.A[k]/Amax
    print(acum[0]/npre,acum[1]/npos,acum[2]/npre,acum[3]/npos,nrun)
    #with open('output_lat_sto_w','a') as nfile:
    #    nfile.write("%f %f %d\n" % (acum[0]/npre,acum[1]/npos,nrun))

###to print latency
with open('output_lat_stp','a') as nfile:
    for k in range(nruns):
        nfile.write("%f %d\n" % (cicle[k]-T[k],k))





























