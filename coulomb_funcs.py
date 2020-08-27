# A compilation of Coulomb wave functions for F, G, and Whittaker function, and their derivatives against k*r
# Note the mpmath's fp coulomb wave functions are not as accurate as believed. Here we use mp. coulomb functions and 
# use fp.mpf to change the results to float numbers. 
# To speed up the calculation, the Coulomb functions from a Fortran code embedded in the r-matrix code is also 
# included here, imported as mycoulfg . 
# So I have code to return mycoulfy calculations in the region the Fortran code is reliable and mpmmath calculations 
# outside that region. These functions are labeled by _mix 
# Also the rescaled Coulomb functions are also provided here, to mitigate the problem of F and G becoming too small 
# or too big in the near zero energy region.
  

import numpy as np
from mpmath import fp,mp                                           
from scipy import special as ss 
from coulomb_Barnett   import mycoulfg
mp.dps=20

coulombf_ufunc=np.frompyfunc(lambda angL, eta, rho: fp.mpf(mp.coulombf(angL, eta, rho ) ), 3, 1) # argment: angL, eta, rho=k*r
coulombg_ufunc=np.frompyfunc(lambda angL, eta, rho: fp.mpf(mp.coulombg(angL, eta, rho ) ), 3, 1)   # argment: angL, eta, rho=k*r
coulombw_ufunc=np.frompyfunc(lambda angL, eta, rho: fp.mpf(mp.whitw(eta,angL+0.5, 2*rho ) ), 3, 1) # argment: angL, eta, rho=k*r 
# argment: angL, eta, rho=k*r, note this order is different from whitw's argument, i
#and there is a factor of 2 for the 3rd argument
drho_coulombf_ufunc=np.frompyfunc(lambda  angL, eta, rho, hvalue, dirvalue: 
        fp.mpf(mp.diff(lambda x: mp.coulombf(angL, eta, x), rho, 1, h=hvalue, direction=dirvalue)) , 5,1)
drho_coulombg_ufunc=np.frompyfunc(lambda  angL, eta, rho, hvalue, dirvalue: 
        fp.mpf(mp.diff(lambda x: mp.coulombg(angL, eta, x), rho, 1, h=hvalue, direction=dirvalue)) , 5,1)  
drho_coulombw_ufunc=np.frompyfunc(lambda  angL, eta, rho, hvalue, dirvalue: 
        fp.mpf(mp.diff(lambda x: mp.whitw(eta,angL+0.5, 2*x), rho, 1, h=hvalue, direction=dirvalue)) , 5,1)

coulombc_ufunc=np.frompyfunc(lambda angL, eta : fp.mpf(mp.coulombc(angL, eta) ), 2 , 1)

coulomb_heta_ufunc = lambda eta : ( ss.psi(1j*eta)  + 0.5/(1j*eta) - np.log(1j*eta) )

def coulombf_rescaled(angL, eta, rho):
    return fp.mpf(mp.coulombf(angL, eta, rho))/fp.mpf(mp.coulombc(angL, eta))*eta**(angL+1)
def coulombg_rescaled(angL, eta, rho):                                                                                                                                                                      return fp.mpf(mp.coulombg(angL, eta, rho))*fp.mpf(mp.coulombc(angL, eta))/eta**angL

def drho_coulombf_rescaled(angL, eta, rho, hvalue, dirvalue): 
    return fp.mpf(mp.diff(lambda x: mp.coulombf(angL, eta, x), rho, 1, h=hvalue, direction=dirvalue))/fp.mpf(mp.coulombc(angL, eta))*eta**(angL+1)                                                                                                                             
def drho_coulombg_rescaled(angL, eta, rho, hvalue, dirvalue):
    return fp.mpf(mp.diff(lambda x: mp.coulombg(angL, eta, x), rho, 1, h=hvalue, direction=dirvalue))*fp.mpf(mp.coulombc(angL, eta))/eta**angL 

coulombf_rescaled_ufunc=np.frompyfunc(coulombf_rescaled, 3, 1) # argment: angL, eta, rho=k*r
coulombg_rescaled_ufunc=np.frompyfunc(coulombg_rescaled, 3, 1) # argment: angL, eta, rho=k*r

drho_coulombf_rescaled_ufunc=np.frompyfunc(drho_coulombf_rescaled, 5, 1)  # argment: angL, eta, rho=k*r, diff step size h, diff direction (0:central, -1 left, 1 right)
drho_coulombg_rescaled_ufunc=np.frompyfunc(drho_coulombg_rescaled, 5, 1)

coulombfdfgdg_ufunc=np.frompyfunc(mycoulfg, 3, 5)

def mycoulfg_rescaled(angL,eta,rho):
    f, df, g, dg, IFAIL = mycoulfg(angL,eta,rho)
    f=f/fp.mpf(mp.coulombc(angL, eta))*eta**(angL+1)
    df=df/fp.mpf(mp.coulombc(angL, eta))*eta**(angL+1)
    g=g*fp.mpf(mp.coulombc(angL, eta))/eta**angL
    dg=dg*fp.mpf(mp.coulombc(angL, eta))/eta**angL
    return f, df, g, dg, IFAIL 

coulombfdfgdg_rescaled_ufunc=np.frompyfunc(mycoulfg_rescaled, 3, 5)

def mycoulfg_mix(angL,eta,rho, hvalue, dirvalue):   
    if rho> (eta + np.sqrt(eta**2+angL*(angL+1))) :
        f,df,g,dg,IFAIL=mycoulfg(angL,eta,rho)
        if(IFAIL!=0):
            print(IFAIL)
    else:
        f=fp.mpf(mp.coulombf(angL, eta, rho ) )
        g=fp.mpf(mp.coulombg(angL, eta, rho ) )
        df=fp.mpf(mp.diff(lambda x: mp.coulombf(angL, eta, x), rho, 1, h=hvalue, direction=dirvalue))
        dg=fp.mpf(mp.diff(lambda x: mp.coulombg(angL, eta, x), rho, 1, h=hvalue, direction=dirvalue))
#    print(f,df,g,dg)    
    return f, df, g, dg 

def mycoulfg_mix_rescaled(angL,eta,rho, hvalue, dirvalue):                                                    
    f,df,g,dg=mycoulfg_mix(angL,eta,rho, hvalue, dirvalue)
    f=f/fp.mpf(mp.coulombc(angL, eta))*eta**(angL+1)
    df=df/fp.mpf(mp.coulombc(angL, eta))*eta**(angL+1)  
    g=g*fp.mpf(mp.coulombc(angL, eta))/eta**angL
    dg=dg*fp.mpf(mp.coulombc(angL, eta))/eta**angL
    return f, df, g, dg 

mycoulfg_mix_ufunc=np.frompyfunc(lambda  angL, eta, rho, hvalue, dirvalue: mycoulfg_mix(angL,eta,rho, hvalue, dirvalue), 5,4)
mycoulfg_mix_rescaled_ufunc=np.frompyfunc(lambda  angL, eta, rho, hvalue, dirvalue: mycoulfg_mix_rescaled(angL,eta,rho, hvalue, dirvalue), 5,4)















