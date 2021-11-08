from re import split
from scipy import optimize
from openmm.app import *
from openmm.unit import *
import numpy as np
import os


# pylint: disable=no-member
import openmm
picoseconds = openmm.unit.picoseconds
picosecond = picoseconds
nanometer = openmm.unit.nanometer
femtoseconds = openmm.unit.femtoseconds
# pylint: enable=no-member

def create_system(psf_file):

    psf = CharmmPsfFile(psf_file)
    psf.setBox(6.24737*nanometer,6.24737*nanometer,6.24737*nanometer)
    ff_dir = os.path.dirname(__file__)
    str_files = ()
    for file in os.listdir(ff_dir):
        if os.path.splitext(file)[-1] == '.str':
            str_files += (os.path.abspath(os.path.join(ff_dir, file)),)
    params = CharmmParameterSet(*str_files)
    system = psf.createSystem(params,  nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=1*nanometer, verbose=True)

    return system