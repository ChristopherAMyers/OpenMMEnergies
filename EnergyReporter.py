import sys
from simtk.openmm.openmm import HarmonicAngleForce, \
        HarmonicBondForce, NonbondedForce, PeriodicTorsionForce, \
        CustomBondForce, CustomNonbondedForce, CMMotionRemover
import simtk.unit as unit

class EnergyReporter(object):
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
            HarmonicBondForce: "Bond",
            HarmonicAngleForce: "Angle",
            PeriodicTorsionForce: "Torsion",
            NonbondedForce: "Nonbonded",
            CustomNonbondedForce: "Custom Nonbonded",
            CustomBondForce: "Custom Bond",
            CMMotionRemover: "C.O.M. Remover"
        }

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
                print(" {:16} {:15.3f} kJ/mol".format(self._type2name.get(type(f), str(f)), energy/unit.kilojoule_per_mole))
                total += energy
            print(" --------------------------------------- ")
        if state:
            total = state.getPotentialEnergy()
        else:
            total = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        
        print(" {:16} {:15.3f} kJ/mol".format("Total energy", total/unit.kilojoule_per_mole))


