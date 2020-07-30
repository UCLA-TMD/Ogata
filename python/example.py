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
q=np.linspace(0.01,1,101)
nu=0

fbt = FBT(nu)

exact = lambda qT: (1+qT**2)**(-1.5)/2./np.pi
wexact = [exact(_q) for _q in q]
wfbt  = [fbt.fbt(test,_q,N,Q) for _q in q]
ratios = [wfbt[i]/wexact[i] for i in range(len(q))]

ax=py.subplot(121)
ax.plot(q,wexact,'-',label='exact')
ax.plot(q,wfbt,'-.',label='numerical')
ax.set_xlabel(r'$q_{\perp}\; \rm (GeV)$',fontsize=20)
ax.set_ylabel(r'$W(q_{\perp})$',fontsize=20)
ax.legend(fontsize=12)
ax=py.subplot(122)
ax.plot(q,ratios)
ax.set_xlabel(r'$q_{\perp}\; \rm (GeV)$',fontsize=20)
ax.set_ylabel(r'\rm FBT/Exact',fontsize=20)
py.tight_layout()
py.show()
