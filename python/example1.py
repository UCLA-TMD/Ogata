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
Q=10.0 # inverse of where test(b) peaks in bt space
q=1.
nu=0.0


fbt = FBT(nu)

wexact = exact(q)
wfbt  = fbt.fbt(test,q,N,Q)
wfbtu  = fbt.fbtu(test,q,N,Q)

print " b*np.exp(-b) :"
print "Exact=", wexact
print "Fbt=", wfbt
print "Fbtu=", wfbtu

wexact1 = exact1(q)
wfbt1  = fbt.fbt(test1,q,N,Q)
wfbtu1  = fbt.fbtu(test1,q,N,Q)

print " b*np.exp(-b**2) :"
print "Exact=", wexact1
print "Fbt=", wfbt1
print "Fbtu=", wfbtu1
