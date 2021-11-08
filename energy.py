#!/usr/bin/env python3
from openmm.app import *
from openmm import *
from openmm.openmm import CustomNonbondedForce, Integrator, RBTorsionForce, VerletIntegrator, NonbondedForce, VirtualSite, HarmonicBondForce
from openmm.unit import *
from openmm.app import element
import argparse
import numpy as np
from copy import copy, deepcopy
import altered_forces
from os import environ, path

from EnergyReporter import EnergyReporter
from molFileReader import molFileReader
from minimize import BFGS
from Drude import drude
from InputOptions import InputOptions
from Nonbonded import *


# pylint: disable=no-member
import openmm
picoseconds = openmm.unit.picoseconds
picosecond = picoseconds
nanometer = openmm.unit.nanometer
femtoseconds = openmm.unit.femtoseconds
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
        if args.xml is not None:
            xml_file = args.xml
            print(" Using included .xml file for the forcefield")
        else:
            xml_file = 'amber99sb.xml'
            print(" Using Amber99 as default forcefield")
        amber = ForceField(xml_file)
        # [templates, residues] = amber.generateTemplatesForUnmatchedResidues(topol)
        # for n, template in enumerate(templates):
        #     amber_template = amber._templates[template.name]
        #     for atom in template.atoms:
        #         print(atom.name)
        #     for atom in template.atoms:
        #         for amber_atom in amber_template.atoms:
        #             if amber_atom.name == atom.name:
        #                 atom.type = amber_atom.type
        #                 break
        #             if atom.type == None:
        #                 atom.type = '1581' # RA-H8

        #     template.name = str(n) + template.name   
        #     amber.registerResidueTemplate(template)

        
        system = amber.createSystem(topol, nonbondedCutoff=4*nanometer, constraints=None)
        
        #   make sure that all h-bonds have a bond force. This is needed when using bond constraints
        #force = system.getForce(0)
        force = [f for f in system.getForces() if isinstance(f, HarmonicBondForce)][0]
        force_bonds = []
        for n in range(force.getNumBonds()):
            force_bonds.append(tuple(sorted(force.getBondParameters(n)[0:2])))
        for bond in topol.bonds():
            a1 = bond.atom1
            a2 = bond.atom2
            top_bond = tuple(sorted([a1.index, a2.index]))
            if top_bond not in force_bonds and \
            (a1.element is element.hydrogen or a2.element is element.hydrogen):
                print("     * Adding H-bond: ", bond)
                force.addBond(a1.index, a2.index, 1*angstrom, 400000)
                force_bonds.append(top_bond)

    return system

def set_self_energies(eng_report, mol):
    if not isinstance(mol, molFileReader.QC):
        raise NotImplementedError("Self energies are only supported for -qcin molecules")
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-pdb', help='PDB file to base the calcualtions off of', required=True)
    parser.add_argument('-psf', help='CHARMM PSF topology file')
    parser.add_argument('-xyz', help='XYZ file of coordinate frames to loop over')
    parser.add_argument('-qcin', help='Q-Chem input file with coordantes to use')
    parser.add_argument('-chg', help='Supplimental column of charges to use')
    parser.add_argument('-top', help='Gromacs topology file with force field info')
    parser.add_argument('-xml', help='OpenMM Force Field XML file to use')
    parser.add_argument('-ipt', help='Input file with options to controll program behavior')
    args = parser.parse_args()
    
    #   get program options
    opts = InputOptions(args.ipt)

    #   load in pdb and assign atom types
    print(" Loading PDB File")
    pdb = PDBFile(args.pdb)
    print(" Done Loading PDB")
    topol = pdb.getTopology()

    for residue in topol.residues():
        if residue.name == 'UNK':
            print('WARNING: Residue name "UNK" might cause bonding issues')

    #   create system object. This holds the forces used
    if args.psf:
        print(" Creating Drude System")
        system = drude.create_system(args.psf)
    else:
        print(" Creating System")
        system = create_system(args, topol)

    #   replace charges with point nuclei and gaussian electron densities
    if opts.density_chg:
        print(" Using density charges")
        altered_forces.gaussian_density(system, topol)

    #   for DEBUG only
    while system.getNumForces() > 1 and False:
        for i in range(system.getNumForces()):
            force = system.getForce(i)
            if not isinstance(force, PeriodicTorsionForce):
                system.removeForce(i)
                break
            

    #   set up custom energy reporter
    if opts.nonbonded_eda:
        create_decomposed_forces(system, topol=topol, exclude_self=opts.nonbonded_res_only)
    eng_report = EnergyReporter(1, system, nonbonded_eda=opts.nonbonded_eda)

    #   set up integrator and simulation objects
    #   integrator is not actually used
    integrator = VerletIntegrator(2*femtoseconds)
    simulation = Simulation(topol, system, integrator)
    simulation.context.setPositions(pdb.getPositions())
    #exit()

    #   replace charges in force field with provided charge list
    if args.chg:
        print(" Replacing force field charges with provided charge file")
        assign_charges(args.chg, system, topol, simulation)

    #   extract coordinates to loop energies over
    coords_to_use = []
    if args.xyz:
        print(" Program will use XYZ frames.")
        mol = molFileReader.XYZ()
        mol.import_xyz(args.xyz)
        for frame in mol.frames:
            coords_to_use.append(np.copy(frame.coords)*angstroms)
    elif args.qcin:
        print(" Program will use Q-Chem input file for coordinates.")
        mol = molFileReader.QC()
        mol.import_qc(args.qcin)
        all_coords = []
        for frag in mol.fragments:
            for coord in frag.coords:
                all_coords.append(list(coord))
        coords_to_use.append(np.array(all_coords)*angstroms)
    else:
        print(" Program will use PDB frames.")
        for n in range(pdb.getNumFrames()):
            coords_to_use.append(pdb.getPositions(asNumpy=True, frame=n))

    if False:
        try:
            import debug      
        except:
            print("FAIL")
            pass
        else:
            psf = CharmmPsfFile(args.psf)
            debug.write_psf_xml(system, topol, psf)
            exit()
            #debug.drude(system, simulation, topol)
            #debug.pull(topol, simulation, pdb)
    
    if opts.optimize:
        print(" Minimizing structure")
        if opts.opt_freeze_main:
            print(" Freezing main particles")
            for n, atom in enumerate(topol.atoms()):
                #if 'POT' in atom.name: continue
                if atom.name[0] != 'D':
                    system.setParticleMass(atom.index, 0*dalton)
        if opts.opt_freeze_drude:
            print(" Freezing Drude particles")
            for n, atom in enumerate(topol.atoms()):
                #if 'POT' in atom.name: continue
                if atom.name[0] == 'D':
                    print(atom)
                    system.setParticleMass(atom.index, 0*dalton)


        if opts.opt_mode == 'bfgs':
            bfgs = BFGS(simulation.context, out_pdb='opt_bfgs.pdb', topology=topol)
            bfgs.minimize()
        elif opts.opt_mode == 'openmm':
            with open('opt_openmm.pdb', 'w') as opt_file:
                for n in range(100):
                    PDBFile.writeModel(topol, simulation.context.getState(getPositions=True).getPositions(), file=opt_file, modelIndex=n)
                    simulation.minimizeEnergy(maxIterations=5)
        opt_state = simulation.context.getState(getPositions=True)
        print(" Minimized Energy:")
        eng_report.report(simulation, opt_state)
    else:
        print(" There are {:d} frames to loop over".format(len(coords_to_use)))
        for n, coords in enumerate(coords_to_use):
            print(" \n Frame {:d}: ".format(n))
            simulation.context.setPositions(coords)
            state = simulation.context.getState(getEnergy=True)
            eng_report.report(simulation, state, total_only=False)


