#!/usr/bin/env python

###############################################################################
#                                                                             #
#                Adaptive Ogata (AdOg) for TMDs                               #
#     Zhongbo Kang, Alexei Prokudin, Nobuo Sato, John Terry                   #
#                   Please cite ArXiv.........                                #
#   Untransformed Ogata (adogu) and transformed (Ogata adogt)                 #
#                      h is spacing parameter                                 #
#                  N is number of function calls                              #
#                  nu is Bessel function order                                #
#                                                                             #
###############################################################################

import numpy as np
from scipy.special import jv, jn_zeros, yv
from scipy.optimize import fsolve,minimize_scalar

class AdOg:

    def __init__(self,nu):
        """Sets maximum number of nodes to about 2^15."""
        self.maxN = 32769
        """Imports zeros of the Bessel function."""
        self.jn_zeros0 = jn_zeros(nu,self.maxN)

    """Transformed Ogata quadrature sum."""
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

    """Untransformed Ogata quadrature sum."""
    def ogatau(self,f,h,N,nu):
        zeros=self.jn_zeros0[:N]
        xi=zeros/np.pi
        Jp1=jv(nu+1,np.pi*xi)
        w=yv(nu, np.pi * xi) / Jp1
        knots = xi*h
        g=lambda x: f(x)*jv(0,x)
        F=g(knots)
        val=h*np.sum(w*F)
        return val#,h*w*F

    """Determines the untransformed hu by maximizing contribution to first node."""
    def get_hu(self,f,nu,q,Q):
        zero1 = self.jn_zeros0[0]
        h = lambda x: -abs(x*f(x/q))
        """Use brent method to maximize."""
        hu = minimize_scalar(h, bracket=None, bounds=(Q,10*Q), args=(), method='brent', tol=0.01, options=None).x/zero1
        return hu

    """Determine transformed ht from untransformed hu."""
    def get_ht(self,hu,nu,N):
        zeroN = self.jn_zeros0[int(N-1)]
        ht = fsolve(lambda h: hu-np.pi*np.tanh(np.pi/2*np.sinh(h*zeroN/np.pi)),2*hu/np.pi/zeroN)[0]
        return ht

    """Adaptive untransformed Ogata."""
    def adogu(self,g,q,N,Q=1.,nu=None):
        if nu is None:
          nu = self.nu
        hu = self.get_hu(g,nu,q,Q)
        f = lambda x: g(x/q)/q
        return self.ogatau(f,hu,N,nu)

    """Adaptive transformed Ogata."""
    def adogt(self,g,q,N,Q=1.,nu=None):
        if nu is None:
          nu = self.nu
        hu = self.get_hu(g,nu,q,Q)
        ht = self.get_ht(hu,nu,N)
        f = lambda x: g(x/q)/q
        return self.ogatat(f,ht,N,nu)


