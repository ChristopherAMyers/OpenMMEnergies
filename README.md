# OpenMMEnergies
Nifty tool for analyzing Force field energies using OpenMM. The main file, energy.py, will loop over any frames in a PDB file and print out its energy components. An .xyz file can be used instead, but the atoms in both the .pdb and .xyz files must align with eachother. If no topology file is provided, the program will use the default AMBER99 force field. 

usage: energy.py [-h] -pdb PDB [-psf PSF] [-xyz XYZ] [-qcin QCIN] [-chg CHG]
                 [-top TOP] [-ipt IPT]

optional arguments:\
  -h, --help  show this help message and exit\
  -pdb PDB    PDB file to base the calcualtions off of\
  -psf PSF    CHARMM PSF topology file\
  -xyz XYZ    XYZ file of coordinate frames to loop over\
  -qcin QCIN  Q-Chem input file with coordantes to use\
  -chg CHG    Supplimental column of charges to use\
  -top TOP    Gromacs topology file with force field info\
  -ipt IPT    Input file with options to controll program behavior\
