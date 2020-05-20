#!/usr/bin/env python2

import numpy as np
from FBT import FBT
import pylab as py
from  matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)

testexp  = lambda b: b*np.exp(-b)          #Input exponential test function
exactexp = lambda qT: (1+qT**2)**(-1.5)    #Analytical result for FT of test function

testgau  = lambda b: b*np.exp(-b**2)       #Input Gaussian test function
exactgau = lambda qT: np.exp(-qT**2/4.)/2. #Analytical result for FT of test function

N=10
Qexp=10.0 # The initial guess for the inverse of where b*exp(-b) peaks in b space 
Qgau=10.0 # The initial guess for the inverse of where b*gau(-b**2) peaks in b space 
#Q=10.0 # inverse of where b*test(b) peaks in b space
q=1.
nu=0.0


fbt = FBT(nu)

print (" b*np.exp(-b) :")
wexpexact  = exactexp(q)
wexpt = fbt.fbt(testexp,q,N,Qexp,0) #FBT using transformed method
wexpu = fbt.fbt(testexp,q,N,Qexp,1) #FBT using untransformed method
wexpo = fbt.fbt(testexp,q,N,Qexp,2) #FBT using traditional Ogata method

erexpt = fbt.fbterror(testexp,q,N,Qexp,0) #FBT error estimate using transformed method
erexpu = fbt.fbterror(testexp,q,N,Qexp,1) #FBT error estimate using untransformed method
erexpo = fbt.fbterror(testexp,q,N,Qexp,2) #FBT error estimate using traditional Ogata method

print ("Exact        =", wexpexact)
print ("Transformed  =", wexpt, ", Relative error estimate = ", erexpt)
print ("Untransformed=", wexpu, ", Relative error estimate = ", erexpu)
print ("Traditional  =", wexpo, ", Relative error estimate = ", erexpo)

print (" b*np.exp(-b**2) :")
wgauexact  = exactgau(q)
wgaut = fbt.fbt(testgau,q,N,Qgau,0) 
wgauu = fbt.fbt(testgau,q,N,Qgau,1) 
wgauo = fbt.fbt(testgau,q,N,Qgau,2) 

ergaut = fbt.fbterror(testgau,q,N,Qgau,0) 
ergauu = fbt.fbterror(testgau,q,N,Qgau,1) 
ergauo = fbt.fbterror(testgau,q,N,Qgau,2) 

print ("Exact        =", wgauexact)
print ("Transformed  =", wgaut, ", Relative error estimate = ", ergaut)
print ("Untransformed=", wgauu, " , Relative error estimate = ", ergauu)
print ("Traditional  =", wgauo, ", Relative error estimate = ", ergauo)
