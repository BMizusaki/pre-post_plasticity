from numpy import *
from brian2 import *
import resource

start_scope()


U_on = 1.
A_on = 1.
###to block presynaptic plasticity
#A_on = 0
###to block postsynaptic plasticity 
#U_on = 0



prefs.codegen.target = 'numpy'
seed=860872
np.random.seed(seed)

#simulation parameters
nruns = 1500
corr_time = 20*ms
timestep = 1*ms

#output neuron parameters

taug = 5*ms
C = 281 * pF
gL = 30 * nS
taum = C / gL
EL = -70.6 * mV
DeltaT = 2 * mV
vt = -50.4 * mV
tauvt = 50 * ms


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


#Input parameters
initialvalue = 0.5
plasticity_on = 0
Nsyn = 200
I = NeuronGroup(Nsyn,'rates : Hz', threshold='rand()<rates*dt', order=0, method='euler', dt=timestep)

freq = 15

N1 = int(Nsyn/2)
I1 = I[:N1]
I2 = I[N1:]

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

tau_u = 50 * ms
tau_r = 200 * ms
s = 0.25*gL
Amax = 1.2
Amin = 0.01
Ainit = 0.1
Umax = 1.
Umin = 0.01
Uinit = 0.1
s2 = sqrt(2.)

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
R-=R*u
u+=U*(1-u)
FFp+=1
'''

S = Synapses(I, P, eqw, on_pre={'pre_t':eqpot}, on_post=eqdep, order=1, method = 'euler', dt=timestep)
S.connect()

S.A = initialvalue*Amax
S.U = initialvalue
S.u = S.U
S.FBp=0
S.FBn=0
S.FFp=0
S.R=1


for nrun in range(nruns):
    if nrun > plasticity_on: #to look at interval before plasticity
        etaA = A_on*0.1
        etaU = U_on*0.1
    else:
        etaA = 0.
        etaU = 0.

    Mf = SpikeMonitor(P,record=True)

    y=np.random.normal(0,1)
    xa=np.random.normal(0,1,Nsyn)
    for i in range(N1):
        I1.rates[i]=(max(0,(freq*(1+0.3*xa[i]+2*y))))*Hz
        I2.rates[i]=(max(0,(freq*(1+0.3*s2*xa[i+N1]))))*Hz
    duration=max(1*ms,corr_time+np.random.normal(0,1)*ms)
    run(duration)

    del Mf

    ###to print average weights
    acum = np.zeros(4)
    for k in range(N1):
        acum[0] =acum[0]+S.U[k]/Umax
        acum[1] =acum[1]+S.U[k+N1]/Umax
        acum[2] =acum[2]+S.U[k]/Umax
        acum[3] =acum[3]+S.A[k+N1]/Amax
    print(acum[0]/N1,acum[1]/N1,acum[2]/N1,acum[3]/N1,nrun)










