Part of the calculations, i.e., solving the scattering Schrodinger equation directly, relie on the R-matrix package. The package is documented at Computer Physics Communications 200 (2016) 199–219, by P. Descouvemont.  

In order to import this package in F90 langauge into Python code, we apply f2py to wrapp the F90 code. For example, to wrap the rmatrix_f2py.f code, run the following command in Linux:
python3.8 -m numpy.f2py -c  -m rmatrix_f2py   rmatrix_f2py.f  -lliblapack 

Note the lapack library is linked in the compilation. The same way can be used to wrap other F90 code in this directory.   

There are three jupyter notebooks for NN, p-alpha, and alpha-Pb scatterings, with their names self-explained. These scatterings have been discussed in the details in the paper. 

Other python files:

Constants.py include the constants used in the calculations. 

coulomb_funcs.py compute coulomb functions needed in the calculations. Two types of calculations are hybrided together, including a F77 code from the R-matrix package, and mpmath python package. The latter offers arbitrary precision calculation, which is needed for computing low-energy processes with very large Sommerfeld parameter. 

two_body_pot.py computes various results for a given strong interaction potential. The Coulomb potential in the point charge form is wired in the rmatrix_f2py.f and rmatrix_f2py_complex_potential.f. Therefore, if the Coulomb potential in a particular calculation is NOT point charge, i.e., the potential at short distance is different from the point charge form, then the difference between the two needs to be included in the strong interaction potential. 

two_body_comp_pot.py works for the complex potential. 

evc_two_body.py provdies EC calculation. Its main input is the array of two_body_pot class instances constructed based on a selected potential parameter sets. 

evc_two_body_comp_pot.py provides EC calculation for complex potential. 

   