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
        """Sets maximum number of nodes to about 2^15."""
        self.maxN = 32769
        """Imports zeros of the Bessel function. Initializing this way speeds up calls"""
        self.jn_zeros0 = jn_zeros(nu,self.maxN)

    def _ogatat(self,f,h,N,nu):
        """Transformed Ogata quadrature sum. Equation 8 in the reference."""
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

    def _ogatau(self,f,h,N,nu):
        """Untransformed Ogata quadrature sum. Equation 7 in the reference."""
        zeros=self.jn_zeros0[:N]
        xi=zeros/np.pi
        Jp1=jv(nu+1,np.pi*xi)
        w=yv(nu, np.pi * xi) / Jp1
        knots = xi*h
        g=lambda x: f(x)*jv(nu,x)
        F=g(knots)
        val=h*np.sum(w*F)
        return val#,h*w*F

    def _get_hu(self,f,q,Q):
        """Determines the untransformed hu by maximizing contribution to first node. Equation 11 in ref."""
        zero1 = self.jn_zeros0[0]
        h = lambda x: -abs(x*f(x/q))
        """Use brent method to maximize."""
        hu = minimize_scalar(h, bracket=None, bounds=(Q/10,10*Q), args=(), method='brent', tol=0.01, options=None).x/zero1
        if hu>3.:
            hu = 3.
            print ('Warning: Number of nodes is too small.')
        return hu

    def _get_ht(self,hu,N):
        """Determine transformed ht from untransformed hu. Equation 13 in ref."""
        zeroN = self.jn_zeros0[int(N-1)]
        #ht = fsolve(lambda h: hu-np.pi*np.tanh(np.pi/2*np.sinh(h*zeroN/np.pi)),2*hu/np.pi/zeroN)[0]
        ht = np.pi/zeroN*np.arcsinh(2/np.pi*np.arctanh(hu/np.pi))
        return ht


    def fbt(self,g,q,N=10,Q=10.,option=0):
        """ Transformed optimized Ogata of a function f.                              """
        """            /Infty                                                         """
        """  result = |  d x g(x) J_nu(q*x)                                           """
        """           /0                                                              """
        """ Parameters                                                                """
        """ g: function of single argument that has a single maximum in [0,Infinity)  """
        """ q: double precision:                                                      """
        """    the position where Fourier transform to be avaluated                   """
        """ N: int optional:                                                          """
        """    number of nodes to use                                                 """
        """ Q: double optional:                                                       """
        """    user's guess of the position of the maximum of function f              """
        """ option: int optional:                                                     """
        """    0 transformed Ogata  optimized h                                       """
        """    1 untransformed Ogata  optimized h                                     """
        """    2 untransformed Ogata   with deafult h = 0.05                          """
        nu = self.nu
        f = lambda x: g(x/q)/q
        if option == 0: # optimized Ogata
            hu = self._get_hu(g,q,Q)
            ht = self._get_ht(hu,N)
            result = self._ogatat(f,ht,N,nu)
        elif option == 1: # untransformed Ogata  optimized h
            hu = self._get_hu(g,q,Q)
            result=self._ogatau(f,hu,N,nu)
        elif option == 2: # untransformed Ogata   with deafult h = 0.05
            hu = 0.05
            result=self._ogatau(f,hu,N,nu)
        return result

    def fbterror(self,g,q,N=10,Q=10.,option=0):
        """Transformed optimized Ogata error estimate."""
        N1 = N
        resultN =  self.fbt(g,q,N1,Q,option);
        N2 = 2*N
        result2N =  self.fbt(g,q,N2,Q,option);
        return np.abs((result2N-resultN)/result2N)



if __name__ == "__main__":
    f = lambda x: x*np.exp(-x)
    fbt = FBT(0)
    print (fbt.fbt(f,1.,20,1.))
