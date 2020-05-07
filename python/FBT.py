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

    def __init__(self,nu):
        self.nu = nu
        """Sets maximum number of nodes to about 2^15."""
        self.maxN = 32769
        """Imports zeros of the Bessel function. Initializing this way speeds up calls"""
        self.jn_zeros0 = jn_zeros(nu,self.maxN)

    """Transformed Ogata quadrature sum. Equation 8 in the reference."""
    def ogatat(self,f,h,N,nu):
        N = int(N)
        zeros =  self.jn_zeros0[:N]
        xi=zeros/np.pi
        Jp1=jv(nu+1,np.pi*xi)
        w=yv(nu, np.pi * xi) / Jp1
        get_psi=lambda t: t*np.tanh(np.pi/2*np.sinh(t))
        get_psip=lambda t:np.pi*t*(-np.tanh(np.pi*np.sinh(t)/2)**2 + 1)*np.cosh(t)/2 + np.tanh(np.pi*np.sinh(t)/2)
        knots=np.pi/h*get_psi(h*xi)
        Jnu=jv(nu,knots)
        psip=get_psip(h*xi)
        F=f(knots)
        psip[np.isnan(psip)]=1.0
        val=np.pi*np.sum(w*F*Jnu*psip)
        return val

    """Untransformed Ogata quadrature sum. Equation 7 in the reference."""
    def ogatau(self,f,h,N,nu):
        zeros=self.jn_zeros0[:N]
        xi=zeros/np.pi
        Jp1=jv(nu+1,np.pi*xi)
        w=yv(nu, np.pi * xi) / Jp1
        knots = xi*h
        g=lambda x: f(x)*jv(nu,x)
        F=g(knots)
        val=h*np.sum(w*F)
        return val#,h*w*F

    """Determines the untransformed hu by maximizing contribution to first node. Equation 11 in ref."""
    def get_hu(self,f,nu,q,Q):
        zero1 = self.jn_zeros0[0]
        h = lambda x: -abs(x*f(x/q))
        """Use brent method to maximize."""
        hu = minimize_scalar(h, bracket=None, bounds=(Q/10,10*Q), args=(), method='brent', tol=0.01, options=None).x/zero1
        if hu>3.:
            hu = 3.
            print 'Warning: Number of nodes is too small.'
        return hu

    """Determine transformed ht from untransformed hu. Equation 13 in ref."""
    def get_ht(self,hu,nu,N):
        zeroN = self.jn_zeros0[int(N-1)]
        #ht = fsolve(lambda h: hu-np.pi*np.tanh(np.pi/2*np.sinh(h*zeroN/np.pi)),2*hu/np.pi/zeroN)[0]
        ht = np.pi/zeroN*np.arcsinh(2/np.pi*np.arctanh(hu/np.pi))
        return ht

    """Untransformed optimized Ogata."""
    def fbtu(self,g,q,N,Q=10.,nu=None):
        if nu is None:
          nu = self.nu
        hu = self.get_hu(g,nu,q,Q)
        f = lambda x: g(x/q)/q
        return self.ogatau(f,hu,N,nu)

    """Transformed optimized Ogata."""
    def fbt(self,g,q,N,Q=10.,nu=None):
        if nu is None:
          nu = self.nu
        hu = self.get_hu(g,nu,q,Q)
        ht = self.get_ht(hu,nu,N)
        f = lambda x: g(x/q)/q
        return self.ogatat(f,ht,N,nu)

if __name__ == "__main__":
    f = lambda x: x*np.exp(-x)
    fbt = FBT(0)
    print fbt.fbt(f,1.,20,1.)
