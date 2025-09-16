
import parmed as pmd
import numpy as np

#The main code that writes LAMMPS input file and data file for AMBER forcefield
def amber2lammps(ars):
    AmberParm=pmd.amber.AmberParm
    printBonds=pmd.tools.actions.printBonds
    printAngles=pmd.tools.actions.printAngles
    printDihedrals=pmd.tools.actions.printDihedrals

    #input: name.mol2, name.frcmod and name.top. These should be in the same folder or current directory
    #output: data.name, parm.name and in.name  
    #name is the name of mol2, frcmod and top file

    with open("data.{}".format(name), "w") as f:
        f.write("LAMMS data file from parmed\n\n")

    with open("parm.{}".format(name), "w") as f:
        f.write("# force field parameters for {}\n\n".format(name))

    
    parm = AmberParm('{}.top'.format(name))
    natoms=len(parm.atoms)

    Bonds=printBonds(parm,"@1-{}".format(natoms))
    with open("bonds.txt", "w") as f:
        print(Bonds,file=f)

    Angles=printAngles(parm,"@1-{}".format(natoms))
    with open("angles.txt", "w") as f:
        print(Angles,file=f)

    Dihedrals=printDihedrals(parm,"@1-{}".format(natoms))
    with open("dihedrals.txt", "w") as f:
        print(Dihedrals,file=f)


    atype={}
    mass=[]
    with open("{}.frcmod".format(name)) as f:
        for line in f:
            if "MASS" in line:
                for i,line_1 in enumerate(f):
                    if "BOND" in line_1:
                        break
                    if line_1.split():
                        #atype.append({line_1.split()[0]:i+1})
                        atype[line_1.split()[0]]=i+1
                        mass.append(line_1.split()[1])
                break

   
    buffer=3.8
    x=[]
    y=[]
    z=[]

    charges=[]
    with open("{}.mol2".format(name)) as f: 
        for line in f:
            if "ATOM" in line:
                for i,line_1 in enumerate(f):
                    if "BOND" in line_1:
                        break
                    if line_1.split():
                        charges.append(float(line_1.split()[8]))
                        x.append(float(line_1.split()[2]))
                        y.append(float(line_1.split()[3]))
                        z.append(float(line_1.split()[4]))
                        pass
                break

    xlo=np.min(x)-buffer
    xhi=np.max(x)+buffer
    ylo=np.min(y)-buffer
    yhi=np.max(y)+buffer
    zlo=np.min(z)-buffer
    zhi=np.max(z)+buffer

    with open("data.{}".format(name), "a") as f:
        f.write("{} atoms \n".format(natoms))
        f.write("{} atom types \n".format(len(atype)))
        f.write("{} bonds \n".format(len(parm.bonds)))
        f.write("{} bond types \n".format(len(parm.bonds)))
        f.write("{} angles \n".format(len(parm.angles)))
        f.write("{} angle types \n".format(len(parm.angles)))
        f.write("{} dihedrals \n".format(len(parm.dihedrals)))
        f.write("{} dihedral types \n\n".format(len(parm.dihedrals)))
        f.write("{} {} xlo xhi \n".format(xlo,xhi))
        f.write("{} {} ylo yhi \n".format(ylo,yhi))
        f.write("{} {} zlo zhi \n\n".format(zlo,zhi))
        f.write("Masses \n\n")
        for i in range(len(atype)):
            f.write("{} {} \n".format(i+1,mass[i]))

        f.write("\nAtoms\n\n")
                
    if np.sum(charges) >= 0:
        charges=charges - (abs(np.sum(charges))/len(charges))
    elif np.sum(charges) < 0:
        charges=charges + (abs(np.sum(charges))/len(charges))

    with open("{}.frcmod".format(name)) as f:
        for line in f:
            if "NONBON" in line:
                for i,line_1 in enumerate(f):
                    if line_1.split():
                        with open("parm.{}".format(name), "a") as flammps:
                            flammps.write("pair_coeff {} {} {} {} # {}\n".format(i+1,i+1,line_1.split()[2],0.890898718140339*2*float(line_1.split()[1]),line_1.split()[0])) 


    with open("{}.mol2".format(name)) as f:  
        for line in f:
            if "ATOM" in line:
                with open("data.{}".format(name), "a") as flammps:
                    for i,line_1 in enumerate(f):
                        if "BOND" in line_1:
                            break
                        if line_1.split():
                            flammps.write("{} {} {} {} {} {} {} {} {} {}\n".format(i+1,1,atype.get(str(line_1.split()[5])),charges[i],line_1.split()[2],line_1.split()[3],line_1.split()[4],0,0,0))
                break

    i=0
    with open("data.{}".format(name), "a") as flammps:
        flammps.write("\nBonds \n\n")
        with open("bonds.txt") as f:
            for line in f:
                if "Atom" not in line and line.split():
                    flammps.write("{} {} {} {}\n".format(i+1,i+1,int(line.split()[0]),int(line.split()[4])))
                    with open("parm.{}".format(name), "a") as flammpsparm:
                        flammpsparm.write("bond_coeff {} {} {}\n".format(i+1,line.split()[9],line.split()[8]))
                    i=i+1


    i=0
    with open("data.{}".format(name), "a") as flammps:
        flammps.write("\nAngles \n\n")
        with open("angles.txt") as f:
            for line in f:
                if "Atom" not in line and line.split():
                    flammps.write("{} {} {} {} {}\n".format(i+1,i+1,int(line.split()[0]),int(line.split()[4]),int(line.split()[8])))
                    with open("parm.{}".format(name), "a") as flammpsparm:
                        flammpsparm.write("angle_coeff {} {} {}\n".format(i+1,line.split()[12],line.split()[13]))
                    i=i+1

    i=0
    with open("data.{}".format(name), "a") as flammps:
        flammps.write("\nDihedrals \n\n")
        with open("dihedrals.txt") as f:
            for line in f:
                if "Atom" not in line and line.split():
                    if str(line.split()[0])=='M' or str(line.split()[0])=='I':
                        flammps.write("{} {} {} {} {} {}\n".format(i+1,i+1,int(line.split()[1]),int(line.split()[5]),int(line.split()[9]),int(line.split()[13])))
                        with open("parm.{}".format(name), "a") as flammpsparm:
                            flammpsparm.write("dihedral_coeff {} 1 {} {} {}\n".format(i+1,line.split()[17],int(float(line.split()[18])),line.split()[19]))
                        i=i+1
                    else:
                        flammps.write("{} {} {} {} {} {}\n".format(i+1,i+1,int(line.split()[0]),int(line.split()[4]),int(line.split()[8]),int(line.split()[12])))
                        with open("parm.{}".format(name), "a") as flammpsparm:
                            flammpsparm.write("dihedral_coeff {} 1 {} {} {}\n".format(i+1,line.split()[16],int(float(line.split()[17])),line.split()[18]))
                        i=i+1

    with open("in.{}".format(name), "w") as f:
        f.write("units real\n")
        f.write("dimension 3\n")
        f.write("boundary p p p\n")
        f.write("atom_style full\n")
        f.write("read_data data.{}\n".format(name))
        f.write("pair_style      lj/cut/coul/long 10\n")
        f.write("kspace_style    pppm 1.0e-4\n")
        f.write("pair_modify     mix arithmetic tail yes\n")
        f.write("bond_style      harmonic\n")
        f.write("angle_style      harmonic\n")
        f.write("dihedral_style    fourier\n")
        f.write("special_bonds   amber\n")
        f.write("include parm.{}\n".format(name))
    









