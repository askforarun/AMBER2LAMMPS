
import MDAnalysis as mda
import subprocess 
import os
from amber2lammps import amber2lammps


name='PVA' # name of the compound 

# if you know the smiles string if not comment the following three lines and make sure name.pdb is in the current folder 
smiles='CCC(O)CC(O)CC(C)O'
u = mda.Universe.from_smiles('{}'.format(smiles))
u.atoms.write('{}.pdb'.format(name))



ff="gaff2"
file_path = './leap.log'
# Check if the file exists before attempting to delete it
if os.path.exists(file_path):
    os.remove(file_path)


#This checks the pdb file for any errors 
cmd1="antechamber -i {}.pdb -o {}.pdb -fi pdb -fo pdb -dr y -rn LIG".format(name,name)
subprocess.run(cmd1, shell=True)
cmd1="antechamber -j 4 -at {} -dr no -fi pdb -fo mol2 -i {}.pdb -o {}.mol2 -c bcc".format(ff,name,name)
subprocess.run(cmd1, shell=True)
cmd1="parmchk2 -i {}.mol2 -o {}.frcmod -f mol2 -a Y".format(name,name)
subprocess.run(cmd1, shell=True)

with open("tleap.in","w") as f:
    f.write("source leaprc.{}\n".format(ff))
    f.write("SUS = loadmol2 {}.mol2\n".format(name)) 
    f.write("check SUS\n")
    f.write("loadamberparams {}.frcmod\n".format(name))
    f.write("saveamberparm SUS {}.top {}.crd\n".format(name,name))
    f.write("quit")

cmd1="tleap -f tleap.in"
subprocess.run(cmd1, shell=True)

#Calling the function 
amber2lammps(name)















































































































































