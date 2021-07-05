import sys
from simtk.openmm.app.forcefield import CMAPTorsion, CustomTorsion, PME
from simtk.openmm import openmm as mm
import simtk.unit as unit

class EnergyReporter(object):
    PME
    """Custom state reporter to print out energy components of a force field"""
    def __init__(self, interval, system, file=sys.stdout):
        """
            Create the EnergyReporter

            Parameters
            ----------
            interval: int
                The interval (in time steps) at which to write energies
            system: openmm.system
                Simulation system to grab forces from
            file: string or file
                The file to write to, specified as a file name or file object
        """
        if isinstance(file, str):
            self._out = open(file, 'w')
        else:
            self._out = file
        self._interval = int(interval)
        self._force_gorups = self._get_force_groups(system)
        self._type2name = { 
            mm.HarmonicBondForce: "Bond",
            mm.HarmonicAngleForce: "Angle",
            mm.PeriodicTorsionForce: "Torsion",
            mm.NonbondedForce: "Nonbonded",
            mm.CustomNonbondedForce: "Custom Nonbonded",
            mm.CustomBondForce: "Custom Bond",
            mm.CMMotionRemover: "C.O.M. Remover",
            mm.RBTorsionForce: "R.B. Torsion",
            mm.CMAPTorsionForce: "CMAP Torsion",
            mm.DrudeForce: "Drude",
            mm.CustomTorsionForce: "Custom Torsion"
        }
        self._energy_terms = {}

    def __del__(self):
        self._out.close()

    def _get_force_groups(self, system):
        """ Assign group numbers to all forces in the system

            Parameters
            ----------
            system: openmm.system
                Simulation system to grab forces from

            Returns
            -------
            A dictionary with keys as the force pointers and values of the 
            force group numbers
        """
        groups = {}
        for i in range(system.getNumForces()):
            force = system.getForce(i)
            force.setForceGroup(i)
            groups[force] = i
        return groups

    def get_energy_terms(self):
        return self._energy_terms

    def describeNextReport(self, simulation):
        """ Get information about the next report

            Parameters
            ----------
            simulation: openmm.Simulation
                The Simulation object to generate the report for
        """
        steps = self._interval - simulation.currentStep % self._interval
        return (steps, False, False, True, True, None)


    def report(self, simulation, state=None, total_only=False):
        """ Get information about the next report

            Parameters
            ----------
            simulation: openmm.Simulation
                The Simulation object to generate the report for
            state: State
                The current state of the system generated from Context
            total_only: bool
                Whether or not to print only the total energy and skip decomposition
        """
        if not total_only:
            print("")
            print(" --------------------------------------- ")
            print("          Energy Decomposition           ")
            print(" --------------------------------------- ")
            total = 0 * unit.kilojoule_per_mole
            for f, i in self._force_gorups.items():
                energy = simulation.context.getState(getEnergy=True, groups=2**i).getPotentialEnergy()
                energy_type = self._type2name.get(type(f), str(f))
                self._energy_terms[energy_type] = energy
                print(" {:16} {:15.3f} kJ/mol".format(energy_type, energy/unit.kilojoule_per_mole))
                total += energy


                # forces = simulation.context.getState(getForces=True, groups=2**i).getForces()
                # for n, force in enumerate(forces):
                #     print(n, force)

            print(" --------------------------------------- ")
        try:
            total = state.getPotentialEnergy()
        except:
            total = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        
        print(" {:16} {:15.3f} kJ/mol".format("Total energy", total/unit.kilojoule_per_mole))


