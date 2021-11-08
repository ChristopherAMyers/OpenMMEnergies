#!/usr/bin/python 
 
import openmm as mm 
from openmm.openmm import app
from openmm.app import forcefield as ff
from openmm.unit import * 
from sys import stdout, exit, stderr, argv
from os import path
#from drudeutility import * 

# pylint: disable=no-member
import openmm
picoseconds = openmm.unit.picoseconds
picosecond = picoseconds
nanometer = openmm.unit.nanometer
femtoseconds = openmm.unit.femtoseconds
# pylint: enable=no-member
 
jobname=argv[1]
 # read PSF 
#psf = app.CharmmPsfFile('ecor1_drude.xplor.psf')
psf = app.CharmmPsfFile(jobname + '.psf')

# load equilibrated coordinates 
#crd = app.PDBFile('ecor1.drude.0.pdb')
crd = app.PDBFile(jobname + '.pdb') 
 
# set initial box vectors 
psf.setBox(6.24737*nanometer,6.24737*nanometer,6.24737*nanometer)  
# read Drude force field 
ff_dir = path.dirname(__file__)
master_str = path.join(ff_dir, 'toppar_drude_master_protein_2013e.str')
drude_str = path.join(ff_dir, 'toppar_drude_nucleic_acid_2017b.str')
rna_str = path.join(ff_dir, 'toppar_all36_na_rna_modified.str')
params = app.CharmmParameterSet(master_str, drude_str, rna_str)
 
# create the system and set attributes 
system = psf.createSystem(params,  nonbondedMethod=ff.PME, nonbondedCutoff=1.2*nanometer, switchDistance=1.0*nanometer, ewaldErrorTolerance = 0.0001, constraints=ff.HBonds) 
 
nbforce = [system.getForce(i) for i in range(system.getNumForces()) if isinstance(system.getForce(i), mm.NonbondedForce)][0] 

nbfix = [system.getForce(i) for i in range(system.getNumForces()) if isinstance(system.getForce(i), mm.CustomNonbondedForce)][0] 
 
nbforce.setNonbondedMethod(mm.NonbondedForce.PME) 
nbforce.setEwaldErrorTolerance(0.0001) 
nbforce.setCutoffDistance(1.2*nanometer) 
nbforce.setUseSwitchingFunction(True) 
nbforce.setSwitchingDistance(1.0*nanometer) 
 
nbfix.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic) 
nbfix.setCutoffDistance(1.2*nanometer) 
nbfix.setUseSwitchingFunction(True) 
nbfix.setSwitchingDistance(1.0*nanometer) 
 
# set up barostat 
system.addForce(mm.MonteCarloBarostat(1*bar, 298*kelvin)) 
 
# set up integrator and dual Langevin thermostat 
integrator = mm.DrudeLangevinIntegrator(298*kelvin, 5/picosecond, 
1*kelvin, 20/picosecond, 0.001*picoseconds) 
# Drude hard wall 
integrator.setMaxDrudeDistance(0.02) 
simulation = app.Simulation(psf.topology, system, integrator) 
simulation.context.setPositions(crd.positions) 
simulation.context.computeVirtualSites() 
 
# print out energy of initial configuration 
state=simulation.context.getState(getPositions=True, getForces=True, getEnergy=True) 
print(state.getPotentialEnergy().value_in_unit(kilocalories_per_mole) )
exit()
 
# start from the equilibrated positions and velocities 
with open(jobname+'.0.rst', 'r') as f: 	 
	simulation.context.setState(mm.XmlSerializer.deserialize(f.read())) 
 
nsavcrd=10000   # save coordinates every 10 ps 
nstep=10000000  # simulate for 10 ns 
nprint=10000    # report every 10 ps 
 
# run the simulation 
dcd=app.DCDReporter(jobname+'.dcd', nsavcrd)
firstdcdstep = nsavcrd 
dcd._dcd = app.DCDFile(dcd._out, simulation.topology, simulation.integrator.getStepSize(), firstdcdstep, nsavcrd)
simulation.reporters.append(dcd) 
simulation.reporters.append(app.StateDataReporter(jobname+'.out', nprint, step=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True))  
 
simulation.step(nstep)
simulation.reporters.pop()
simulation.reporters.pop() 
 
# write restart file at the end of the run 
state = simulation.context.getState( getPositions=True, getVelocities=True ) 
with open(jobname+'.rst', 'w') as f: 
 	f.write(mm.XmlSerializer.serialize(state)) 
