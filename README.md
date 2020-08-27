# Notebooks for Eigenvalue Continuation (EC) for scattering 

The Jupyter notebooks and associated codes can be used to reproduce and extend the figures in [arXiv:2007.03635](https://arxiv.org/abs/2007.03635), "Efficient emulators for scattering using eigenvector continuation" by R.J. Furnstahl, A.J. Garcia, P.J. Millican, and Xilin Zhang.

Parts of the calculations, i.e., solving the scattering Schrodinger equation directly, rely on the R-matrix package documented at [Computer Physics Communications 200 (2016) 199–219](https://www.sciencedirect.com/science/article/pii/S0010465515003951?via%3Dihub), by P. Descouvemont. The needed Fortran codes are already included in the current directory. They have been modfified here so that they can be further processed as discussed below. The original package itself can be downloaded from [http://cpc.cs.qub.ac.uk/summaries/AEYP_v1_0.html](http://cpc.cs.qub.ac.uk/summaries/AEYP_v1_0.html).  

In order to import this package, written in F90, into the Python code, we apply `f2py` to wrap the F90 code. 

* To wrap the `rmatrix_f2py.f` code, run the following command in a terminal window: <br>
`python -m numpy.f2py -c  -m rmatrix_f2py   rmatrix_f2py.f  -lliblapack` <br>

* To wrap the `rmatrix_f2py_complex_potential.f` code, run the following command in a terminal window: <br>
`python -m numpy.f2py -c  -m rmatrix_f2py_complex_potential    rmatrix_f2py_complex_potential.f  -lliblapack` <br>

* To wrap the `coulomb_Barnett.f` code, run the following command in a terminal window: <br>
`python -m numpy.f2py -c  -m coulomb_Barnett    coulomb_Barnett.f ` <br>

Note that the lapack library is linked in the compilations of `rmatrix_f2py.f` and `rmatrix_f2py_complex_potential.f`, which is required for the current set up of `cminv_nsym()` and `cminv_sym()` subroutines in the `rmatrix_f2py.f` and `rmatrix_f2py_complex_potential.f`. See the comments in those subroutines in the two Fortran codes; also see "section 3. Description of the package" in the aforementioned paper by P. Descouvemont.   

There are three Jupyter notebooks: for NN (nucleon-nucleon), p-alpha (proton-alpha particle), and alpha-Pb (alpha particle-lead) scatterings. Details are discussed in the paper. 

Other Python files: 

* `Constants.py` includes the constants used in the calculations. 

* `coulomb_funcs.py` computes the Coulomb functions needed in the calculations. Two types of calculations are hybrided together, including a F77 code from the R-matrix package and the mpmath python package. The latter enables arbitrary precision calculations, which are needed for computing low-energy processes with a very large Sommerfeld parameter. 

* `two_body_pot.py` computes various results for a given strong-interaction potential. The Coulomb potential in the point-charge form is wired into the `rmatrix_f2py.f` and `rmatrix_f2py_complex_potential.f`. Therefore, if the Coulomb potential in a particular calculation is NOT point charge, i.e., the potential at short distance is different from the point-charge form, then the difference between the two needs to be included in the strong interaction potential. 

* `two_body_comp_pot.py` works for the complex potential. 

* `evc_two_body.py` provides an EC calculation. Its main input is the array of two_body_pot class instances constructed based on a selected potential parameter sets. 

* `evc_two_body_comp_pot.py` provides an EC calculation for complex potentials. 

   
