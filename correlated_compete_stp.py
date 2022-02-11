from numpy import *
from brian2 import *
import resource
from scipy.special import erf

start_scope()
prefs.codegen.target = 'numpy'


#to suppress log messages
get_logger()
logger.suppress_hierarchy('brian2')

#Function to rectify Gaussian (to generate corr inputs)
def gaus(rater):
#    np.seterr(invalid='warn')
    a = 1.+np.real(erf(rater/(np.sqrt(2))+0j))
    mu2=np.array((1./np.sqrt(2.*np.pi))*np.exp(-0.5*rater*rater)+.5*rater*a,dtype=np.float64)
    ns = np.zeros(2)
    ns[0] = (rater-mu2)*mu2
    ns[1] = 0.5*a
    sigma2=np.array(np.sqrt(ns[0]+ns[1]),dtype=np.float64)
    return mu2,sigma2 


seed=860872
np.random.seed(seed)


#simulation parameters
nruns = 200
step_time = 50*ms
defaultclock.dt = 0.1*ms

#Correlation time scale
rate = 25*Hz #average input spiking frequency
cr = 0.9 #correlation coefficient
tauc = 20*ms
delta = defaultclock.dt/tauc
olambda=np.exp(-delta)

#output neuron parameters
tauv = 20*ms
taug = 5*ms
E_v = -74*mV


#synaptic parameters
tau_stdp = 20*ms
initialvalue = 0.5
qmax = 0.25
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
NA1 : 1
NA2 : 1
NA = NA1 + NA2 : 1
'''

P = NeuronGroup(1, eqs, threshold='v>-54*mV', reset='v=-74*mV', refractory=1*ms, order=2)
P.v = -74*mV
P.ge = 0

#Input------------------------------
Nsyn = 200 
N1 = int(Nsyn/2)
eqwss1 = '''
U : 1
vs : 1
initi : 1
NA1_post = U : 1 (summed)
dR/dt=(1-R)/tau_r : 1 (clock-driven)
dp/dt=(U-p)/tau_p : 1 (clock-driven)
dApre/dt = -Apre/tau_stdp : 1 (event-driven) 
dApost/dt = -Apost/tau_stdp : 1 (event-driven)
'''
eqwss2 = '''
U : 1
vs : 1
initi : 1
NA2_post = U : 1 (summed)
dR/dt=(1-R)/tau_r : 1 (clock-driven)
dp/dt=(0.5-p)/tau_p : 1 (clock-driven)
dApre/dt = -Apre/tau_stdp : 1 (event-driven) 
dApost/dt = -Apost/tau_stdp : 1 (event-driven)
'''

#changing P
eqpotss1 = '''
vs = clip(int(R*nv),0,1)
w = int(rand()+p)*vs
ge_post += 0.5 * qmax * w
R-=R*w/nv
p+=U*(1-p)
dw = 0.5*Apost*w
Apre += Ap*w
U = clip(U+dw, valmin, valmax)
'''
eqdepss1 = '''
Apost += Ad
dw = 0.5*Apre
U = clip(U+dw, valmin, valmax)-clip(0.5,NA_post/Nsyn,10.) +0.5
'''

#changing q
eqpotss2 = '''
vs = clip(int(R*nv),0,1)
w = int(rand()+p)*vs
ge_post += U * qmax * w
R-=R*w/nv
p+=0.5*(1-p)
dw = 0.5*Apost*w
Apre += Ap*w
U = clip(U+dw, valmin, valmax)
'''
eqdepss2 = '''
Apost += Ad
dw = 0.5*Apre
U = clip(U+dw, valmin, valmax)-clip(0.5,NA_post/Nsyn,10.) +0.5
'''

valmin = 0.01
valmax = 0.99

I1 = NeuronGroup(N1,'v : 1 (shared)',threshold = 'rand()<(v+mu)')
I1.run_regularly('v = v*olambda+randn()*sg ')
I1.v = 0

I2 = NeuronGroup(N1,'v : 1 (shared)',threshold = 'rand()<(v+mu)')
I2.run_regularly('v = v*olambda+randn()*sg ')
I2.v = 0

S1 = Synapses(I1, P, eqwss1, on_pre={'pre_t': eqpotss1}, on_post=eqdepss1, order=1)
S2 = Synapses(I2, P, eqwss2, on_pre={'pre_t': eqpotss2}, on_post=eqdepss2, order=1)

S1.connect()
S2.connect()

S1.p = initialvalue
S2.p = initialvalue
S1.U=initialvalue
S1.R = 1
S2.U=initialvalue
S2.R = 1

sigmar=np.sqrt(cr*rate/(2.*tauc))
xa = np.array(rate/sigmar,dtype=np.float64)
mur,sig = gaus(xa)

#Transform to inverse Gaussian
fx = mur/sig-xa
x1 = np.array(1.1*rate/sigmar,dtype=np.float64)
step = 0
#Find solution
while True: 
    mur,sig = gaus(x1)
    fx0 = fx
    fx = np.array(mur/sig-rate/sigmar,dtype=np.float64)
    xb = x1
    x1 = x1 - ((x1-xa)/(fx-fx0))*fx
    xa = xb
    step+=1
    if(np.fabs(fx)<1e-4):
        break

sigmav = rate/(np.exp(-0.5*x1*x1)/(np.sqrt(2.*np.pi))+.5*x1*(1.+erf(x1/np.sqrt(2.))))
ratev = x1*sigmav

#Dispersion from average rate
sg = sigmav*np.sqrt(1.-np.exp(-2*delta))*defaultclock.dt
#Shared average rate 
mu = ratev*defaultclock.dt



#run simulation
t0 = 0
for nrun in range(nruns):            
    #Ms = StateMonitor(P,'v',record=True,dt=1*ms)

    duration=max(1*ms,step_time+np.random.normal(0,1)*ms)
    run(duration)
    t0 = t0 + duration/ms


    #del Ms
    with open('output_w','wb') as nfile:
        Nw1=sum(S1.U)/(N1)#average weight changing P
        Nw2=sum(S2.U)/(N1)#average weight changing q
        print(nrun,Nw1,Nw2)



