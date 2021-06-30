import sys
from simtk.openmm.openmm import HarmonicAngleForce, \
        HarmonicBondForce, NonbondedForce, PeriodicTorsionForce, \
        CustomBondForce, CustomNonbondedForce, CMMotionRemover
import simtk.unit as unit

class EnergyReporter(object):
    def __init__(self, interval, system, file=sys.stdout):
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

    def describeNextReport(self, simulation):
            steps = self._interval - simulation.currentStep % self._interval
            return (steps, False, False, True, True, None)

    def report(self, simulation, state):
        print("")
        print("          Energy Decomposition           ")
        print(" --------------------------------------- ")
        total = 0 * unit.kilojoule_per_mole
        for f, i in self._force_gorups.items():
            energy = simulation.context.getState(getEnergy=True, groups=2**i).getPotentialEnergy()
            print(" {:16} {:15.3f} kJ/mol".format(self._type2name[type(f)], energy/unit.kilojoule_per_mole))
            total += energy 
        print(" {:16} {:15.3f} kJ/mol".format("Total energy", total/unit.kilojoule_per_mole))
        print(" --------------------------------------- ")


    def _get_force_groups(self, system):
        groups = {}
        for i in range(system.getNumForces()):
            force = system.getForce(i)
            force.setForceGroup(i)
            groups[force] = i
        return groups
