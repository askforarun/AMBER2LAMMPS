# AMBER2LAMMPS
python utility code to create LAMMPS data file

 pmd2lmp is the function that generates data file and input file 

 simply run main.py. Its that simple 

You need to have antechamber, mdanalysis and RdKit (if you want use smiles) installed. 

Example is for short PVA (polymer) chain.

Usage:

pmd2lmp(name)

lmp_png < in.name

#input: name.mol2, name.frcmod and name.top.

#output: data.name, parm.name and in.name  

#name is the name of mol2, frcmod and top file and these three files should be in the same folder

#main.py generates the three files 
