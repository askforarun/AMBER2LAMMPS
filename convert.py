
import MDAnalysis as mda
import subprocess 
import os
from amber_to_lammps import amber2lammps

#This file generates mol2, frcmod and top file

name='epon' # name of the pdb file

# if you know the smiles strings 
# smiles='CCC(O)CC(O)CC(C)O'
# u = mda.Universe.from_smiles('{}'.format(smiles))
# u.atoms.write('{}.pdb'.format(name))



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
    f.write("saveamberparm SUS {}.prmtop {}.crd\n".format(name,name))
    f.write("quit")

cmd1="tleap -f tleap.in"
subprocess.run(cmd1, shell=True)

#Calling the function 

amber2lammps("epon.prmtop","epon.mol2","epon.frcmod")


#Write lammps input file
with open("in.lammps", "w") as f:
    f.write("units real\n")
    f.write("dimension 3\n")
    f.write("boundary p p p\n")
    f.write("atom_style full\n")
    f.write("read_data data.lammps\n")
    f.write("pair_style      lj/cut/coul/long 9 9\n")
    f.write("kspace_style    pppm 1.0e-8\n")
    f.write("pair_modify     tail yes\n")
    f.write("bond_style      harmonic\n")
    f.write("angle_style      harmonic\n")
    f.write("dihedral_style    fourier\n")
    f.write("special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 0.83333333\n")
    f.write("include parm.lammps\n")
    f.write("thermo_style custom ebond eangle edihed eimp epair evdwl ecoul elong etail pe\n")
    f.write("run 0")












































































































































