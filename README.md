CorrectBoxDriftMain.f90 and CorrectBoxDriftmodule.f90 are the main program and module program, respectively, which could be used to eliminate position deviation between the reference lattice system and the simulation systems with periodic boundary conditions for the purpose of correct identification of vacancies, self-interstitial atoms, and their clusters when using the defect location methods, such as Wigner-Seitz cell and Lindemann spheres methods. They can be added to Fortran editors (such as Intel Visual Fortran Composer, Compaq Visual Fortran 6.5/6.6, and Fortran PoweStation4.0) for compilation, linking, and running. 

For the BCC crystal structure, we provide two examples of input files, namely Config_W_LA20LB20LC20_19SIA_1500K and ReferLatt_W_LA20LB20LC20. The former contains the positions of all atoms in the W box that undergoes significant drift, while the latter contains the positions of all sites in the reference lattice system. 

Similarly, for the FCC structure, two input files are provided, namely Config_Cu_LA10LB10LC10_3SIA_700K and ReferLatt_Cu_LA10LB10LC10.

The format of input files, the settings of input parameters, and the functions of each subroutine can be found in the comments section of these programs.
