program CorrctBoxDrift
!***This program could be used to eliminate position deviation between the reference lattice system and 
!   the simulation systems with periodic boundary conditions for the purpose of correct identification 
!   of vacancies, self-interstitial atoms, and their clusters when using the defect location methods, 
!   such as Wigner-Seitz cell and Lindemann spheres methods.
!*** All length units are in lattice length units

use CorrctBoxDriftModule
implicit none
 
!***Input_AtomsPosi_name specifies the name and path of the input file. This input file contains information about the types and
!   postions of all atoms in the system.In this input file, the first column represents the type of each atom, and the second,
!   third, and fourth columns represent the position coordinate components of each atom, respectively. 
!   The default value for the type of substrate atom is 1 .    
 
!***Input_ReferLatt_name specifies the name and path of the input file.This input file contaning information about the postions 
!   of all sites in a perfect lattice system.In this input file, the first,second,and third columns represent the position 
!   coordinate components of each lattice site, respectively. 
 
!***Output_AtomsPosi_name specifies the name and path of the output file.The output file containing information about the types 
!   and the postions of all the atoms in the system after correction.

!***The input file "Config_W-LA20LB20LC20_19SIA_1500K" here is an example file.It contains the types and positions of all atoms in 
!   the W crystal that undergoes significant drift. The example file "ReferLatt_W_LA20LB20LC20" contains the positions of all sites 
!   in the reference lattice system.The following parameters are all set based on this example.

character*256::Input_AtomsPosi_name='C:\test\Config_W_LA20LB20LC20_19SIA_1500K' 
character*256::Input_ReferLatt_name='C:\test\ReferLatt_W_LA20LB20LC20' 
character*256::Output_AtomsPosi_name='C:\test\correct_Config_W_LA20LB20LC20_19SIA_1500K'

!*** Box sizes with periodic boundary conditions in x, y, z direction
real*8::LA     =20                     
real*8::LB     =20        
real*8::LC     =20  

!***Crystal structural type. 1: BCC structure;   2: FCC structure   3: HCP structure 
integer::CrystalType=1

!*** The type of substrate atom crosspondng to the input file. The default value is set to 1. 
integer::subatomtype=1

!*** An adjustable parameter representing the range of atomic vibrations
real*8::Radius=0.06D0

!***The type of boundary.  The default value is 1.
!*** BoundType=1: The boundary positions in the X direction are set at -LA/2 and LA/2, respectively
!***              The boundary positions in the Y direction are set at -LB/2 and LB/2, respectively
!***              The boundary positions in the Z direction are set at -LC/2 and LC/2, respectively
!*** BoundType=2: The boundary positions in the X direction are set at 0 and LA, respectively
!***              The boundary positions in the Y direction are set at 0 and LB, respectively
!***              The boundary positions in the Z direction are set at 0 and LC, respectively
integer::BoundType=1  

!***If the crystal structure is HCP (CrystalType=3), please provide the axial ratio(Hcpca=c/a) corresponding to the material.
!   For example, Hcpca=1.586 for Ti, and for Zn, Hcpca=1.861. When CrystalType is not equal to 3, this parameter is not working.
Hcpca=1.587D0

!***To read in the data from the input files
call Initialize(Input_AtomsPosi_name,Input_ReferLatt_name) 

!***To Search for atoms in the system, which are able to maintain the crystal structure(BCC,FCC or HCP).
call Search_For_Atoms_CrystalType(CrystalType,Radius,subatomtype) 

!*** To determine the position deviation between the reference lattice system and the simulation system.
call Determine_The_Drift_of_Box()

!***To correct the postions of all the atoms in box according the deviation.
call Correct_ALL_Atoms_Position(LA,LB,LC,BoundType)

!***To write the data into the output file.
call Write_Into_Output_File(Output_AtomsPosi_name)

end program CorrctBoxDrift
