# AMBER2LAMMPS
python utility code to create LAMMPS data file

chek main.py which calls pm2lmp. 

You need to have antechamber, mdanalysis and RdKit (if you want use smiles) installed. 

Example is for short PVA (polymer) chain.

Usage:

pmd2lmp(name)

lmp_png < in.name

#input: name.mol2, name.frcmod and name.top. These should be in the same folder or current directory
#output: data.name, parm.name and in.name  
#name is the name of mol2, frcmod and top file
