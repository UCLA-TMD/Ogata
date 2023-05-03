#!/usr/bin/env python

###############################################################################
#                                                                             #
#                Fast Bessel Transform (FBT) for TMDs                         #
#     Zhongbo Kang, Alexei Prokudin, Nobuo Sato, John Terry                   #
#                   Please cite ArXiv:1906.05949                              #
#                      h is spacing parameter                                 #
#                  N is number of function calls                              #
#                  nu is Bessel function order                                #
#                                                                             #
###############################################################################

import numpy as np
from scipy.special import jv, jn_zeros, yv
from scipy.optimize import fsolve,minimize_scalar

class FBT:

    def __init__(self,nu=0.0):
        """ Constructor, sets nu of Jnu()                          """
        self.nu = nu
        self.setup()
        
    def setup(self,N = 100):
        self.jn_zeros0 = jn_zeros(self.nu,N)
        self.xi = self.jn_zeros0/np.pi
        Jp1=jv(self.nu+1,np.pi*self.xi)
        self.w=yv(self.nu, np.pi * self.xi) / Jp1
        
    def get_psi(self, t):
        return t*np.tanh(np.pi/2*np.sinh(t))
     
    def get_psip(self, t):
        return np.pi*t*(-np.tanh(np.pi*np.sinh(t)/2)**2 + 1)*np.cosh(t)/2 + np.tanh(np.pi*np.sinh(t)/2)

    def _ogatat(self,f,h,N,nu):
        """Transformed Ogata quadrature sum. Equation 8 in the reference."""
        N = int(N)
        if N > 100:
           self.setup(N)
        knots=np.pi/h*self.get_psi(h*self.xi[:N])
        Jnu=jv(nu,knots)
        psip=self.get_psip(h*self.xi[:N])
        # AP 5/3/23 This line does not when the functions does not accept a list as the argument
        #F=f(knots)
        # the solution:
        F=[f(knot) for knot in knots]
        psip[np.isnan(psip)]=1.0
        val=0.5*np.sum(self.w[:N]*F*Jnu*psip)
        return val

    def _ogatau(self,f,h,N,nu):
        """Untransformed Ogata quadrature sum. Equation 7 in the reference."""
        if N > 100:
           self.setup(N)
        knots = self.xi[:N]*h
        g=lambda x: f(x)*jv(nu,x)
        # AP 5/3/23 This line does not when the functions does not accept a list as the argument
        #F=g(knots)
        # the solution:
        F=[g(knot) for knot in knots]
        val=h*np.sum(self.w[:N]*F)/2./np.pi
        return val

    def _get_hu(self,f,q,Q):
        """Determines the untransformed hu by maximizing contribution to first node. Equation 11 in ref."""
        zero1 = self.jn_zeros0[0]
        h = lambda x: -abs(x*f(x/q))
        """Use brent method to maximize."""
        hu = minimize_scalar(h, bracket=None, bounds=(Q/10,10*Q), args=(), method='brent', tol=0.01, options=None).x/zero1*np.pi
        if hu>2.:
            hu = 2.
        return hu

    def _get_ht(self,hu,N):
        """Determine transformed ht from untransformed hu. Equation 13 in ref."""
        zeroN = self.jn_zeros0[int(N-1)]
        ht = np.pi/zeroN*np.arcsinh(2/np.pi*np.arctanh(hu/np.pi))
        return ht


    def fbt(self,g,q,N=10,Q=10.,option=0):
        """ Transformed optimized Ogata of a function f.                              """
        """            /Infty                                                         """
        """  result = |  d x g(x) J_nu(q*x)                                           """
        """           /0                                                              """
        """ Parameters                                                                """
        """ g: function of a single argument that has a single maximum in [0,Infinity)"""
        """ q: double precision:                                                      """
        """    the position where Fourier transform to be avaluated                   """
        """ N: int optional:                                                          """
        """    number of nodes to use                                                 """
        """ Q: double optional:                                                       """
        """    user's guess of the position of the maximum of function f              """
        """ option: int optional:                                                     """
        """    0 transformed Ogata  optimized h                                       """
        """    1 untransformed Ogata  optimized h                                     """
        """    2 transformed Ogata   with default h = 0.05                          """
        nu = self.nu
        f = lambda x: g(x/q)/q
        if option == 0:   # Transformed Ogata optimized h
            hu = self._get_hu(g,q,Q)
            ht = self._get_ht(hu,N)
            result = self._ogatat(f,ht,N,nu)
        elif option == 1: # Untransformed Ogata optimized h
            hu = self._get_hu(g,q,Q)
            result=self._ogatau(f,hu,N,nu)
        elif option == 2: # Transformed Ogata with h = 0.05
            hu = 0.05
            result=self._ogatat(f,hu,N,nu)
        return result

    def fbterror(self,g,q,N=10,Q=10.,option=0):
        """Transformed optimized Ogata error estimate."""
        N1 = N
        resultN =  self.fbt(g,q,N1,Q,option);
        N2 = 2*N
        result2N =  self.fbt(g,q,N2,Q,option);
        return np.abs((result2N-resultN)/result2N)

    def fbt_findN(self,g,q,err,Q=10,option=0):
        _N = 2
        _Q= Q
        _option = option
        W1 = self.fbt(g,q,N=  _N,Q=_Q,option=_option)
        W2 = self.fbt(g,q,N=2*_N,Q=_Q,option=_option)
        while abs(W1-W2)/W2>err:
            _N = 2*_N
            W1 = self.fbt(g,q,N=  _N,Q=_Q,option=_option)
            W2 = self.fbt(g,q,N=2*_N,Q=_Q,option=_option)
        return _N
