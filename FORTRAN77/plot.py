#!/usr/bin/env python2

import numpy as np
import pylab as py
from  matplotlib import rc
import pandas as pd
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)

df = pd.read_csv("output.dat", delimiter = ',',delim_whitespace=True)

ax=py.subplot(121)
ax.plot(df.qT,df.exact,'-',label='exact')
ax.plot(df.qT,df.fbt, '-.',label='numerical')
ax.set_xlabel(r'$q_{\perp}\; \rm (GeV)$',fontsize=20)
ax.set_ylabel(r'$W(q_{\perp})$',fontsize=20)
ax.legend(fontsize=12)
ax=py.subplot(122)
ax.plot(df.qT,df.fbt/df.exact)
ax.set_xlabel(r'$q_{\perp}\; \rm (GeV)$',fontsize=20)
ax.set_ylabel(r'\rm fbt/Exact',fontsize=20)
py.tight_layout()
py.show()
