from openmm.app import *
import parmed

# pylint: disable=no-member
import openmm.unit
picoseconds = openmm.unit.picoseconds
picosecond = picoseconds
nanometer = openmm.unit.nanometer
femtoseconds = openmm.unit.femtoseconds
# pylint: enable=no-member

def create_system(prmtop_file):
    #top = AmberPrmtopFile(prmtop_file)
    #system = top.createSystem(nonbondedCutoff=NoCutoff)

    parm = parmed.amber.AmberParm(prmtop_file)
    system = parm.createSystem()
    
    return system