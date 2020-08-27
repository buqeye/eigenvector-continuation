import numpy as np
from scipy.integrate import quadrature 
from scipy import special as ss
from Constants import * 
from  coulomb_funcs import  mycoulfg_mix_rescaled_ufunc, coulombw_ufunc, drho_coulombw_ufunc, \
coulombc_ufunc, coulomb_heta_ufunc, coulombf_ufunc, mycoulfg_mix_ufunc 
from  rmatrix_f2py import rmat_ini,pot_output_rmatrix, gauleg 

#np.seterr('raise')

class two_body_pot:
    """
    Define a two-body potenial object with tuneable parameters

    Radial wave functions u(r) = r*R(r) are used. 

    Issues and to-do list:

    Parameters
    ----------
    pot_func : non-local potential function that takes a vector np.array([r, rprime])
               default 0 . In the unit of MeV fm^-1
    local_pot_func : local potential function takes argument of r 
                     In the unit of MeV
    mu : float
        Reduced mass.  mu = 1 by default but explicit in class.
    hbar : float
        Planck's constant divided by 2\pi. hbar = 1 (and not explicit).
    E_mesh : mesh points for energy where we will compute scattering
            a vector
    r_c : the channel radius for r-matrix analysis, beyond which the 
          strong interaction becomes negligable .  In the unit of fm
    r_max : the largest radius that we want to compute wave function. 
           in the unit of fm
    r_in_mesh_pts : int 
        Number of r mesh points up to r_c 
    r_out_mesh_pts: int
        Number of r mesh between r_c and r_max 
    angL : partial wave 

    z_t : charge of target 

    z_p : charge of projectile 

    compute_ps: True to compute phase shift, False for not computing it.     
    
    kc: the Coulomb momentum  

    test_not_coulomb: logical, True for non-Coulomb, False for Coulomb case

    r_in_mesh: the mesh for r < r_c 
    r_in_weight: the weight for r_in_mesh 

    r_out_mesh: the mesh for r > r_c
    r_out_weight: the weight for r_out_mesh 

    pot_mesh: the non-local potential computed at r and rp in the 
              r_in_mesh  \otimes r_in_mesh 
 
    local_pot_mesh : the local potential computed at r on the r_in_mesh 
                     Note in the r_matrix calculation, a point Coulomb potential 
                     is already included. If the Coulomb potetnail in the 
                     problem is not point, the difference should be included 
                     either here if the difference is local potential, or 
                     included in the pot_mesh, if the difference is non-local. 
                     Of course the potetnials should be computed from the pot_func 
                     and or local_pot_func.    
 
    E_mesh_pts: Number of meshes in the E_mesh 

    k_mesh : the k mesh corresponding to the E_mesh 
 
    delta_mesh : phase shift on the E_mesh

    scatt_wf_in_re_mesh: real part of the scattering wave function at r_in_mesh 
    scatt_wf_in_im_mesh: imaginary part of the scattering wave function at r_in_mesh 
    
    tau_mesh: the tau on the E_mesh 

    Methods
    asym_WF_sc(r, E, rescale=False): computes asymptotic scattering wave functions
              and dervatives wrt to k*r, rescale=True would rescale the 
              wave function such that the Coulomb case can be better handled. 
    asym_WF_b(self, r, E): compute the asymptotic wave function for bound state 
    eta_func(k) : compute the Sommerfeld parameter eta 
    get_ere_mesh() : compute ere function at the energy mesh points. 
                   for the Coulomb case, the first return is the full ERE 
                   and the 2nd is the good ERE, while for the non-Coulomb case, 
                   both are the same. 
                   It returns ere_mesh  

    get_scatt_wf_out_mesh() : compute scattering wave function at r_out_mesh
                              it returns scatt_wf_out_mesh  
    -------
    
    """
    def __init__(self, pot_func = lambda rvec: 0,  
                 local_pot_func= lambda r : 0,  
                 mu=M_N, hbar=hbarc, E_mesh=np.array([0.1,1]), r_c=10., r_max=20., r_in_mesh_pts=48, 
                 r_out_mesh_pts=48, angL=0, z_t=0, z_p = 0, compute_ps=True, compute_wf_out=False ):
        self.pot_func = pot_func 
        self.local_pot_func = local_pot_func
        self.mu = mu
        self.hbar = hbar
        self.E_mesh=E_mesh
        self.r_c=r_c
        self.r_max = r_max 
        self.r_in_mesh_pts=int(r_in_mesh_pts)
        self.r_out_mesh_pts=int(r_out_mesh_pts)
        self.angL=int(angL) # partial wave
        self.z_t=z_t  # target's EM charge
        self.z_p=z_p  # projectile's charge
        self.kc= z_t*z_p*alpha_EM*mu/hbar # kc momentum in the unit of wave number
        if int(z_t*z_p)==0:
            self.test_not_coulomb=True
        else :
            self.test_not_coulomb=False
        zrma,r_c0,wle,xle, tc, blo0, blo1, blo2, q1, q2 = rmat_ini(self.r_in_mesh_pts,self.r_c)
        self.r_in_mesh=zrma 
        self.r_in_weight=wle * self.r_c  
        self.r_out_mesh, self.rout_weight = gauleg(self.r_c,self.r_max,self.r_out_mesh_pts)   
        r_in_mesh_2= np.array([ [ [r,rp]  for rp in self.r_in_mesh]  for r in self.r_in_mesh] )
        self.pot_mesh= np.apply_along_axis(self.pot_func, 2,  r_in_mesh_2)
        self.local_pot_mesh = np.vectorize(self.local_pot_func)(self.r_in_mesh)
        E_mesh_pts=np.shape(self.E_mesh)[0]
        self.k_mesh=np.sqrt(2*self.mu*self.E_mesh)/self.hbar
        if compute_ps : 
            self.delta_mesh, self.scatt_wf_in_re_mesh, self.scatt_wf_in_im_mesh  =  \
                  pot_output_rmatrix(self.hbar,alpha_EM,self.mu,self.z_t,self.z_p,
                                E_mesh_pts,self.E_mesh,self.r_in_mesh_pts,self.r_c, 
                               self.r_in_mesh,self.pot_mesh,self.local_pot_mesh, 
                               r_c0,wle,xle,tc,blo0,blo1,blo2,q1,q2,self.angL)
            self.tau_mesh= np.tan(self.delta_mesh)/self.k_mesh 
            self.get_ere_mesh() 
        if compute_wf_out:
            self.get_scatt_wf_out_mesh() 

 
    def asym_WF_sc(self, r, E, rescale=False):                                                                         
        """
        The asymptotic wave functions and their derivatives (d/drho!) 
        for both neutral and charged particle scattering cases.
        Note the g and dg function have a minus sign in front of it. 
        """
        k=np.sqrt(2*self.mu*E)/self.hbar
        rho=k *r 
        eta=self.eta_func(k)
#        print('in asym_WF_sc', "E,  k, eta, r,  rho" )
#        print( E,  k, eta, r, rho )
        if np.sum(np.shape(rho))==0 : 
            if np.absolute(rho)<1.e-16: rho=1.e-16
        else:
            rho[np.isclose(rho,0, atol=1.e-16)]=1.e-16
# in order to deal with the case rho=0
        if self.test_not_coulomb:
            tp=ss.spherical_jn(self.angL, rho)
            f = rho* tp  
            df = tp + rho*ss.spherical_jn(self.angL, rho, derivative=True)
            tp=ss.spherical_yn(self.angL, rho)
            g = -rho* tp
            dg = -tp - rho*ss.spherical_yn(self.angL, rho, derivative=True) 
        else:
            hvalue=None
            dirvalue=None
            if rescale:
                f,df,g,dg = mycoulfg_mix_rescaled_ufunc(self.angL,eta,rho, hvalue, dirvalue) 
                f=np.array(f).astype(float)
                g=np.array(g).astype(float)
                df=np.array(df).astype(float)
                dg=np.array(dg).astype(float)
            else :    
                f,df,g,dg = mycoulfg_mix_ufunc(self.angL,eta,rho, hvalue, dirvalue)
                f=np.array(f).astype(float) 
                g=np.array(g).astype(float) 
                df=np.array(df).astype(float) 
                dg=np.array(dg).astype(float) 
#        print("f,df,g,dg")
#        print(f,df,g,dg)
        return (f, df), (g, dg)

    def asym_WF_b(self, r, E):                                                                         
        """
        The asymptotic wave function and its derivative (d/drho!) 
        for both neutral and charged particle bound states. 
        """
        gammaB=np.sqrt(-2*self.mu*E)/self.hbar  # binding momentum
        rho=gammaB *r 
        eta=self.eta_func(gammaB) # corresponding to binding momentum gammaB
        if np.sum(np.shape(rho))==0 :
            if np.absolute(rho)<1.e-16: rho=1.e-16
        else:
            rho[np.isclose(rho,0, atol=1.e-16)]=1.e-16

        if self.test_not_coulomb:
            f=ss.spherical_kn(self.angL, rho)*2./np.pi*rho # note 2/pi factor and the rho factor
            df=f/rho+ss.spherical_kn(self.angL, rho, derivative=True)*2./np.pi*rho
        else:
            hvalue=None
            dirvalue=None
            f=np.array(coulombw_ufunc(self.angL, -eta,  rho)).astype(float)
            df=np.array(drho_coulombw_ufunc(self.angL, -eta,  rho, hvalue, dirvalue)).astype(float)
        return f, df

    def eta_func(self,k):
        """
        compute eta variable
        """
        if self.test_not_coulomb:
#            print(self.kc,k)
            return 0
        else:
            return  self.kc/k
   
    def get_ere_mesh(self):
        """                                                                                                            
        return w(angL)*(c_{eta,angL})^2 k^{2angL+1} (\cot\delta(k) + [0,  2 eta*k^(2 angL +1) *v*Re[Heta] ] ), 
        and phase-shift, with w(angL)=(Gamma[2 l + 2]/Gamma[l + 1]/2^l)^2; 
        v=\prod_{j=1}^angL(1+eta^2/j^2) (note v=1 for angL=0);  
        c_{eta,l} defined in coulombc_ufunc in Coulomb_funcs.py.
        This ERE definition smoothly transitions to the ere functin of the eta=0 case for neutral particle scattering
        In order to save computation time, the delta phase shift is put as the 3rd entry. 
        """
        L=self.angL
        if self.test_not_coulomb:
            goodere=self.k_mesh**(2*L+1)/np.tan(self.delta_mesh) 
            self.ere_mesh = goodere, goodere
        else : 
            eta=self.kc/self.k_mesh
            cetaLsq=(np.array(coulombc_ufunc(L,eta)).astype(float))**2
            w=(ss.factorial(2*L+1)/ss.factorial(L)/2**L)**2
            if L==0:
                v=1.
            else:
                v=np.prod(np.array([1+eta**2/j**2 for j in range(1,L+1) ]), axis=0)
            heta=np.array(coulomb_heta_ufunc(eta)).astype(complex)
            fullere= self.k_mesh**(2*L+1) * w * cetaLsq / np.tan(self.delta_mesh)
            goodere=fullere+self.k_mesh**(2*L+1)*2*eta* v * np.real(heta)
            self.ere_mesh = fullere, goodere
#        return self.ere_mesh 
        
    def get_scatt_wf_out_mesh(self):
        """
        Exterior scattering wave function with energy E at radius r, normalized
         with 1/k factor with k = np.sqrt(2*mu*E).
        """       
        (f,df), (g,dg)=  self.asym_WF_sc(self.r_out_mesh, np.array([self.E_mesh]).T  )
        self.scatt_wf_out_mesh = np.array([ f[i]/k + np.tan(self.delta_mesh[i])*g[i]/k   
                                              for i, k in enumerate(self.k_mesh ) ] ) 
#        return self.scatt_wf_out_mesh 

#################################################################################

