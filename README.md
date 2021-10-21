# OpenMMEnergies
Nifty tool for analyzing Force field energies using OpenMM. The main file, energy.py, will loop over any frames in a PDB file and print out its energy components. An .xyz file can be used instead, but the atoms in both the .pdb and .xyz files must align with eachother. If no topology file is provided, the program will use the default AMBER99 force field. 

## Usage

```
python3 energy.py [-h] -pdb PDB [-psf PSF] [-xyz XYZ] [-qcin QCIN] [-chg CHG]
                 [-top TOP] [-ipt IPT]

optional arguments:
  -h, --help  show this help message and exit
  -pdb PDB    PDB file to base the calcualtions off of
  -psf PSF    CHARMM PSF topology file
  -xyz XYZ    XYZ file of coordinate frames to loop over
  -qcin QCIN  Q-Chem input file with coordantes to use
  -chg CHG    Supplimental column of charges to use
  -top TOP    Gromacs topology file with force field info
  -xml XML    OpenMM Force Field XML file to use
  -ipt IPT    Input file with options to controll program behavior
  ```

## Input File
The optional input file uses a Q-Chem format where each input section starts with `$section` and ends with `$end`. For example, the main options section , `$rem`, with their default values applied would look something like this:

```
$rem
  optimize            False     ! optimize structure first; uses only the first frame of PDB coords
  print_eda           True      ! print energy decomposition
  print_forces        False     ! print forces too (not implimented yet)
  opt_mode            openmm    ! optimization mode: openmm or build_in BFGS solver
  opt_freeze_main     False     ! freeze main atoms (not drude particles)
  opt_freeze_drude    False     ! freeze drude particles
  density_chg         True      ! replace charges with point nuclei and Slater electron densities
$end
```

Comment lines begin with a `!` character.
