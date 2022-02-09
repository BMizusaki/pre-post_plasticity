eqws = '''
U : 1
vs : 1
dR/dt=(1-R)/tau_r : 1 (clock-driven)
dp/dt=(U-p)/tau_p : 1 (clock-driven)
dApre/dt = -Apre/tau_stdp : 1 (event-driven) 
dApost/dt = -Apost/tau_stdp : 1 (event-driven)
'''
eqwsq = '''
U : 1
vs : 1
dR/dt=(1-R)/tau_r : 1 (clock-driven)
dp/dt=(initialvalue-p)/tau_p : 1 (clock-driven)
dApre/dt = -Apre/tau_stdp : 1 (event-driven) 
dApost/dt = -Apost/tau_stdp : 1 (event-driven)
'''
#equations for stochastic synapse--- change p
eqpotP = '''
vs = clip(int(R*nv),0,1)
w = int(rand()+p)*vs
ge_post += initialvalue * qmax * w
R-=R*w/nv
p+=U*(1-p)
dw = initialvalue*Apost*w
Apre += Ap*w
U = clip(U+dw, valmin, valmax)
'''
eqdepP = '''
Apost += Ad
dw = initialvalue*Apre
U = clip(U+dw, valmin, valmax)
'''
#-----------------------------------

#equations for stochastic synapse--- change q
eqpotq = '''
vs = clip(int(R*nv),0,1)
w = int(rand()+p)*vs
ge_post += U * qmax * w
R-=R*w/nv
p+=initialvalue*(1-p)
dw = initialvalue*Apost*w
Apre += Ap*w
U = clip(U+dw, valmin, valmax)
'''
eqdepq = '''
Apost += Ad
dw = initialvalue*Apre
U = clip(U+dw, valmin, valmax)
'''
#-----------------------------------


#equations for stochastic synapse--- change both
eqpotB = '''
vs = clip(int(R*nv),0,1)
w = int(rand()+p)*vs
ge_post += U * qmax * w
R-=R*w/nv
p+=U*(1-p)
Apre += Ap*w
U = clip(U-(U-sqrt(U*U+initialvalue*Apost))*w, valmin, valmax)
'''
eqdepB = '''
Apost += Ad
U = clip(sqrt(U*U+initialvalue*Apre)), valmin, valmax)
'''
#-----------------------------------
