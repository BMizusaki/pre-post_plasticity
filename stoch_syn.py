
#equations for stochastic synapse--- change p
eqpotP = '''
w = int(rand()+U)
ge_post += 0.5 * qmax * w
dw = 0.5*Apost*w
Apre += Ap*w
U = clip(U+dw, valmin, valmax)
'''
eqdepP = '''
Apost += Ad
dw = 0.5*Apre
U = clip(U+dw, valmin, valmax) 
'''
#-----------------------------------

#equations for stochastic synapse--- change q
eqpotq = '''
w = int(rand()+0.5)
ge_post += U * qmax * w
dw = 0.5*Apost*w
Apre += Ap*w
U = clip(U+dw, valmin, valmax)
'''
eqdepq = '''
Apost += Ad
dw = 0.5*Apre
U = clip(U+dw, valmin, valmax) 
'''
#-----------------------------------

#equations for stochastic synapse--- change both
eqpotB = '''
w = int(rand()+U) 
ge_post += U * qmax * w 
Apre += Ap * w
U = clip(U-(U-sqrt(U*U+0.5*Apost))*w, valmin, valmax) 
'''
#p = clip(p+Apost*w, 0, 1)
eqdepB = '''
Apost += Ad
U = clip(sqrt(U*U+0.5*Apre), valmin, valmax)
'''
#-----------------------------------
