# AMBER2LAMMPS
python utility code to create LAMMPS data file



 simply run convert.py. Its that simple 

You need to have antechamber, mdanalysis and RdKit (if you want use smiles) installed. 

For installing antechmaber check 

(https://ambermd.org/AmberTools.php)

and for MDanalysis check 

https://www.mdanalysis.org

and for RdKit check 

https://www.rdkit.org/docs/Install.html



Usage:

lmp_png < in.lammps

Code check with InterMOL:


For installing Intermol check 

https://github.com/shirtsgroup/InterMol


InterMOL:


Per MPI rank memory allocation (min/avg/max) = 132.1 | 132.1 | 132.1 Mbytes

   E_bond        E_angle        E_dihed        E_impro         E_pair         E_vdwl         E_coul         E_long         E_tail         PotEng   

 2.3161274      6.0940384      12.475809      0             -8.8739004      10.824738      97.869973     -117.56861     -0.0044166818   12.012074     

Loop time of 7.5e-07 on 1 procs for 0 steps with 53 atoms

My code:

  E_bond        E_angle        E_dihed        E_impro         E_pair         E_vdwl         E_coul         E_long         E_tail         PotEng    

 2.3161274      6.0940126      12.475827      0             -9.7069268      9.9917113      97.869973     -117.56861     -0.0043506337   11.17904   


Before check make sure the box size is the same in the data file 

