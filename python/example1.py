#!/usr/bin/env python2

import numpy as np
from FBT import FBT
import pylab as py
from  matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)

test = lambda b: b*np.exp(-b)
exact = lambda qT: (1+qT**2)**(-1.5)

test1 = lambda b:   b*np.exp(-b**2)
exact1 = lambda qT: np.exp(-qT**2/4.)/2.

N=10
Q=10.0 # inverse of where b*test(b) peaks in b space
q=1.
nu=0.0


fbt = FBT(nu)

wexact = exact(q)
wfbt  = fbt.fbt(test,q,N,Q,0)
wfbtu  = fbt.fbt(test,q,N,Q,1)

print( "Relative error estimate of doubling N = ", N, "at q = ", q," is ",fbt.fbterror(test1,q,N,Q))

print (" b*np.exp(-b) :")
print ("Exact=", wexact)
print ("Fbt=", wfbt)
print ("Fbtu=", wfbtu)

wexact1 = exact1(q)
wfbt1  = fbt.fbt(test1,q,N,Q,0)
wfbtu1  = fbt.fbt(test1,q,N,Q,1)

print (" b*np.exp(-b**2) :")
print ("Exact=", wexact1)
print ("Fbt=", wfbt1)
print ("Fbtu=", wfbtu1)
