import numpy as np
import scipy.linalg as la
from scipy.linalg import eig, eigh
from scipy.integrate import simps, quadrature
from scipy import special as ss
from  coulomb_funcs import  mycoulfg_mix_rescaled_ufunc, coulombw_ufunc, drho_coulombw_ufunc, \
coulombc_ufunc, coulomb_heta_ufunc, coulombf_ufunc, mycoulfg_mix_ufunc
from two_body_comp_pot  import two_body_pot



def nonlocal_matrix_element_gauss(wf1, wf2, r_nodes, r_weights, op=None):

    """
    Calculate the matrix element of NONLOCAL operator op,  using quadratures
    from the two_body_pot objects.
    op when None, op = delta(r-rp) on mesh, which is delta_{r,rp}/weight. 
    """
    if op is None:
        op = np.diag(1/r_weights)

    return  (wf1 * r_weights)  @ op  @ (wf2 * r_weights)   

def local_matrix_element_gauss(wf1, wf2, r_nodes, r_weights, op=None):

    """
    Calculate the matrix element of op (which is taken to be
     the identity if omitted) using Gaussian quadratures from the two_body_pot
    objects. 
    """
    if op is None:
        op = np.ones(len(wf1)) 

    return (wf1 * op * wf2) @ r_weights


def ere_from_tau(angL=0, kc=0, k=10, tau=0.1, test_not_coulomb=True):
        """
        based on tau values, compute good ERE:
        w(angL)*(c_{eta,angL})^2 k^{2angL+1} (\cot\delta(k) +  2 eta*k^(2 angL +1) *v*Re[Heta] ] ,
        with w(angL)=(Gamma[2 l + 2]/Gamma[l + 1]/2^l)^2;
        v=\prod_{j=1}^angL(1+eta^2/j^2) (note v=1 for angL=0);
        c_{eta,l} defined in coulombc_ufunc in Coulomb_funcs.py.
        This ERE definition smoothly transitions to the ere functin of the eta=0 case for neutral particle scattering
        """
        L=angL
        tandelta=tau*k
        if test_not_coulomb:
            goodere=k**(2*L+1)/tandelta
        else :
            eta=kc/k
            cetaLsq=(np.array(coulombc_ufunc(L,eta)).astype(float))**2
            w=(ss.factorial(2*L+1)/ss.factorial(L)/2**L)**2
            if L==0:
                v=1.
            else:
                v=np.prod(np.array([1+eta**2/j**2 for j in range(1,L+1) ]), axis=0)
            heta=np.array(coulomb_heta_ufunc(eta)).astype(complex)
            fullere= k**(2*L+1) * w * cetaLsq / tandelta
            goodere=fullere+k**(2*L+1)*2*eta* v * np.real(heta)
        return goodere


class EigenvectorContinuationScattering:
    """
    Set up and carry out eigenvectorcontinuation for the two_body_comp_pot objects
      but for scattering instead of bound states. 
      NOTE THIS IS FOR COMPLEX POTENTIAL.
      We'll assume that hbar and mu and r_c and r_max and the mesh are the 
      same for all potentials.
     
    Should calculate the matrix element of the potential by using matrix
      elements pre-calculated rather than doing it again.
    
      
    Parameters
    ----------
    pc_array : array of potential objects for the appropriate Hamiltonian
    
    Methods
    ---------
    find_EVC_scattering: takes in nonlocal potential function, new_pot_func, and
                         local potential function, new_local_pot_func, and use
                         EVC to compute the full potential's scatterings at the
                         E_mesh from the pc_array objects.
                         Three ways to handle ill-conditioned matrices:
                         direct inverse using la.inv,
                         or the pseduo inverse np.linalg.pinv
                         or adding nuggest.
                         the cond would be input for the latter two
                         It returns EVC-estimated tau, tav_var_mesh, on the E_mesh,
                         and the corresponding ere_var_mesh, c_vec_mesh, lag_mesh,
                         and delta_tilde_U_condition_mesh.
    """
    
    def __init__(self, pc_array):
        self.pc_array = pc_array    
         
    
    def find_EVC_scattering(self, new_pot_func = lambda rvec: 0 ,
                 new_local_pot_func= lambda r : 0 , pinv=True, nugget=False, cond=1.e-11):
        """Try out EVC for scattering
        """
        mu = self.pc_array[0].mu
        hbar = self.pc_array[0].hbar
        E_mesh=self.pc_array[0].E_mesh
        r_c= self.pc_array[0].r_c 
        r_max= self.pc_array[0].r_max
        r_in_mesh_pts=self.pc_array[0].r_in_mesh_pts
        r_out_mesh_pts=self.pc_array[0].r_out_mesh_pts
        angL=self.pc_array[0].angL
        z_t=self.pc_array[0].z_t
        z_p=self.pc_array[0].z_p
         
        new_pc = two_body_pot(pot_func=new_pot_func, local_pot_func=new_local_pot_func, 
                              mu=mu, hbar=hbar, E_mesh=E_mesh, r_c=r_c, r_max=r_max,
                              r_in_mesh_pts=r_in_mesh_pts, r_out_mesh_pts= r_out_mesh_pts, 
                              angL=angL,z_t=z_t,z_p=z_p, compute_ps=False)
## In the two_body_pot, we need to set it NOT to compute phase-shift using r-matrix, 
# otherwise, EVC looses the point. 
        num_terms= np.shape(self.pc_array)[0]
        k_mesh=self.pc_array[0].k_mesh
        r_in_mesh = self.pc_array[0].r_in_mesh
        r_in_weight = self.pc_array[0].r_in_weight  
        tau_mesh_vec =np.array( [pc.tau_mesh for pc in self.pc_array] )
        d_pot_mesh_vec = [-pc.pot_mesh + new_pc.pot_mesh for pc in self.pc_array]   
        d_local_pot_mesh_vec= [- pc.local_pot_mesh + new_pc.local_pot_mesh  for pc in self.pc_array]     
        d_full_pot_mesh_vec=[ d_pot_mesh_vec[i] + np.diag(d_local_pot_mesh_vec[i]/r_in_weight)  
                                         for i, pc in enumerate(self.pc_array)]  

        delta_U = np.zeros( (num_terms, num_terms) )*1j
        delta_tilde_U_condition_mesh =np.array([])
        tau_var_mesh=np.array([])
        ere_var_mesh=np.array([]) 
        lag_mesh=np.array([]) 
        for iE, E in enumerate(E_mesh): 
            for i, pc_i in enumerate(self.pc_array):
                for j, pc_j in enumerate(self.pc_array):
                    delta_U[i,j] =  (2*mu/hbar**2) \
                        * nonlocal_matrix_element_gauss( pc_i.scatt_wf_in_mesh[iE], 
                                                  pc_j.scatt_wf_in_mesh[iE], 
                                                  r_in_mesh, r_in_weight,
                                                  d_full_pot_mesh_vec[j] )
            
            delta_tilde_U = delta_U + delta_U.T
            if pinv:
                delta_tilde_U_inv = np.linalg.pinv(delta_tilde_U, rcond=cond)
            elif nugget:
                delta_tilde_U = delta_tilde_U + np.eye(len(delta_tilde_U)) * cond
                delta_tilde_U_inv = la.inv(delta_tilde_U)
            else :
                delta_tilde_U_inv = la.inv(delta_tilde_U)
            delta_tilde_U_condition_mesh =np.append(delta_tilde_U_condition_mesh, np.linalg.cond(delta_tilde_U))
            tau_vec=tau_mesh_vec[:,iE]
        
            lagrange = ( np.sum(delta_tilde_U_inv @ tau_vec) - 1) \
                      / np.sum(delta_tilde_U_inv)
            c_vec = delta_tilde_U_inv @ (tau_vec - lagrange)
#            c_vec_solve = la.solve(delta_tilde_U, tau_vec - lagrange)
        #print(delta_tilde_U.shape, delta_tilde_U_inv.shape, lagrange.shape, c_vec.shape)
       
            tau_var = c_vec @ tau_vec - c_vec.T @ delta_U @ c_vec
            ere_var = ere_from_tau(angL=angL, kc=self.pc_array[0].kc, k=k_mesh[iE], 
                                  tau=tau_var, test_not_coulomb=self.pc_array[0].test_not_coulomb)
#            tau_var_solve = c_vec_solve @ tau_vec - c_vec_solve.T @ delta_U @ c_vec_solve
#            tau_var_t = c_vec @ tau_vec
#            tau_var_alt = tau_var_t - c_vec.T @ delta_tilde_U @ c_vec / 2.
            tau_var_mesh=np.append(tau_var_mesh, tau_var)
            ere_var_mesh=np.append(ere_var_mesh, ere_var)
            lag_mesh=np.append(lag_mesh, lagrange)
            if iE ==0 :
                c_vec_mesh=np.array([c_vec])
            else:
                c_vec_mesh=np.append(c_vec_mesh, [c_vec], axis=0)

#            print(delta_U)
        return tau_var_mesh, ere_var_mesh, c_vec_mesh, lag_mesh, delta_tilde_U_condition_mesh
        # , lagrange, delta_tilde_U, delta_tilde_U_inv
      
