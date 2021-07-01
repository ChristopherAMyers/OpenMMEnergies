import re
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.openmm.openmm import CustomNonbondedForce, Integrator, VerletIntegrator, NonbondedForce
#from simtk.openmm.openmm import *
from simtk.unit import *
import argparse
import numpy as np
from copy import copy
import altered_forces
from os import environ, path

from EnergyReporter import EnergyReporter
from molFileReader import molFileReader
from InputFileParser import InputFile

# pylint: disable=no-member
import simtk
picoseconds = simtk.unit.picoseconds
picosecond = picoseconds
nanometer = simtk.unit.nanometer
femtoseconds = simtk.unit.femtoseconds
# pylint: enable=no-member

def assign_charges(chg_file, system, topology, simulation):
    """ Change the charges in a forcefield

        Parameters
        ----------
        system: openmm.System
            System that contains either NonbondedForce or CustomNonbondedForce force
        topology: openmm.Topology
            Topology object (typically from PDB) to construct system from
        simulation: openmm.Simulation
            The Simulation object to update parameters in context
            
        Returns
        -------
        System object to use with a Simulation object
    """
    charges = np.loadtxt(chg_file)
    nb_force = None
    for n, force in enumerate(system.getForces()):
        if isinstance(force, CustomNonbondedForce):
            nb_force = force
            break
        elif isinstance(force, NonbondedForce):
            nb_force = force
            break

    for atom in topology.atoms():
        params = list(nb_force.getParticleParameters(atom.index))
        params[0] = charges[atom.index]*elementary_charge*1.0
        if isinstance(nb_force, CustomNonbondedForce):
            nb_force.setParticleParameters(atom.index, params)
        else:
            nb_force.setParticleParameters(atom.index, *params)

    nb_force.updateParametersInContext(simulation.context)

def create_system(args, topol):
    """ Create the system of forces

        Parameters
        ----------
        args: args
            Script arguments parsed from ArgumentParser
        topol: Topology
            Topology object (typically from PDB) to construct system from
            
        Returns
        -------
        System object to use with a Simulation object
    """
    system = None
    if args.top:
        home = environ['HOME']
        inc_dir = path.join(home, 'ChenRNALab/bin/gromacs-2019.4/bin/share/gromacs/top')
        top = GromacsTopFile(args.top, includeDir=inc_dir)
        system = top.createSystem(nonbondedCutoff=2*nanometer)
    else:
        print("Using Amber99 as default forcefield")
        amber = ForceField('amber99sb.xml')
        [templates, residues] = amber.generateTemplatesForUnmatchedResidues(topol)
        for n, template in enumerate(templates):
            amber_template = amber._templates[template.name]
            for atom in template.atoms:
                for amber_atom in amber_template.atoms:
                    if amber_atom.name == atom.name:
                        atom.type = amber_atom.type
                        break
                    if atom.type == None:
                        atom.type = '1581' # RA-H8
            template.name = str(n) + template.name   
            amber.registerResidueTemplate(template)

        #   create system
        system = amber.createSystem(topol, nonbondedCutoff=2*nanometer)

    return system
    
def get_options(input_file):
    parser = InputFile()
    parser.add_argument('dens', default=False, help='Replace charges with Gaussian electron densities')
    parser.add_argument('eda', default=False, help='Print the various energy terms for each frame')
    parser.add_argument('')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-pdb', help='PDB file to base the calcualtions off of', required=True)
    parser.add_argument('-xyz', help='XYZ file of coordinate frames to loop over')
    parser.add_argument('-chg', help='Supplimental column of charges to use')
    parser.add_argument('-top', help='Gromacs topology file with force field info')
    parser.add_argument('-dens', help='Replace charges with Gaussian electron densities', action='store_true')
    args = parser.parse_args()

    #   load in pdb and assign atom types
    pdb = PDBFile(args.pdb)
    topol = pdb.getTopology()
    for residue in topol.residues():
        if residue.name == 'UNK':
            print('WARNING: Residue name "UNK" might cause bonding issues')

    #   create system object. This holds the forces used
    system = create_system(args, topol)

    #   replace charges with point nuclei and gaussian electron densities
    if args.dens:
        altered_forces.gaussian_density(system, topol)

    #   for DEBUG only
    while system.getNumForces() > 1 and False:
        for i in range(system.getNumForces()):
            force = system.getForce(i)
            if not isinstance(force, PeriodicTorsionForce):
                system.removeForce(i)
                break

    #   set up custom energy reporter
    eng_report = EnergyReporter(1, system)

    #   set up integrator and simulation objects
    #   integrator is not actually used
    integrator = VerletIntegrator(2*femtoseconds)
    simulation = Simulation(topol, system, integrator)
    simulation.context.setPositions(pdb.getPositions())

    #   replace charges in force field with provided charge list
    if args.chg:
        print("Replacing force field charges with provided charge file")
        assign_charges(args.chg, system, topol, simulation)

    #   extract coordinates to loop energies over
    coords_to_use = []
    if args.xyz:
        print(" Program will use XYZ frames.")
        mol = molFileReader.XYZ()
        mol.import_xyz(args.xyz)
        for frame in mol.frames:
            coords_to_use.append(np.copy(frame.coords)*angstroms)
    else:
        print(" Program will use PDB frames.")
        for n in range(pdb.getNumFrames()):
            coords_to_use.append(pdb.getPositions(asNumpy=True, frame=n))
    print(" There are {:d} frames to loop over".format(len(coords_to_use)))


    #   get self_energy
    if True:
        self_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()*0.0
        if topol.getNumResidues() == 2:
            pos = copy(pdb.getPositions())
            for atom in list(topol.residues())[1].atoms():
                pos[atom.index] += Vec3(10, 0, 0)*nanometer
            simulation.context.setPositions(pos)
            self_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    
    if False:
        for x in np.arange(-0.5, 2.5, 0.1):
            pos = np.copy(pdb.getPositions(True)/angstrom)*angstrom
            pos[15:] += np.array([x, 0, 0])*angstroms
            simulation.context.setPositions(pos)
            energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
            print("{:5.2f}  {:12.3f} {:s}".format(x, (energy - self_energy)/energy.unit, str(energy.unit)))
        exit()

    
    for n, coords in enumerate(coords_to_use):
        print(" \n Frame {:d}: ".format(n))
        simulation.context.setPositions(coords)
        state = simulation.context.getState(getEnergy=True)
        eng_report.report(simulation, state, total_only=False)

