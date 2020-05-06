#!/usr/bin/env python2

import numpy as np
from FBT import FBT
import pylab as py
from  matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)

test = lambda b: b*np.exp(-b)

N=10
Q=1.0 # inverse of where test(b) peaks in bt space
q=0.1
nu=0.0


fbt = FBT(nu)

exact = lambda qT: (1+qT**2)**(-1.5)
wexact = exact(q)
wfbt  = fbt.fbt(test,q,N,Q,nu)
ratios = wfbt/wexact

print "Exact=", wexact
print "Fbt=", wfbt
