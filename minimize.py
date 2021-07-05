from scipy import optimize
from simtk.openmm.openmm import Context_swigregister
from simtk.unit import *
import numpy as np

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

constr_2_idx = {
    'X': [0], 'Y': [1], 'Z': [2],
    'XY': [0, 1], 'XZ': [1,2], 'YZ': [1,2],
    'YX': [0, 1], 'ZX': [1,2], 'ZY': [1,2],
    'XYZ': [0, 1, 2]
}

def BFGS(context, constraints=None):

    #constraints = dict(zip(np.arange(64), ['Z']*64))

    init_state = context.getState(getForces=True, getEnergy=True, getPositions=True)
    init_pos = init_state.getPositions(True).value_in_unit(nanometer)
    init_energy, init_forces = target_func(init_pos, context, constraints)
    force_norms = [np.linalg.norm(f) for f in init_forces]
    print(" Initial max. force: {:15.3f} kJ/mol".format(np.max(force_norms)))
    print(" Initial energy:     {:15.3f} kJ/mol/nm".format(init_energy))

    args = (context, constraints)
    res = optimize.minimize(target_func, init_pos, args=args, method='L-BFGS-B', jac=True, 
    options=dict(maxiter=500, disp=False, gtol=5))
    final_pos = res.x.reshape(-1,3)

    final_energy, final_forces = target_func(final_pos, context, constraints)
    force_norms = [np.linalg.norm(f) for f in final_forces]
    print(" Final max. force:   {:15.3f} kJ/mol".format(np.max(force_norms)))
    print(" Final energy:       {:15.3f} kJ/mol/nm".format(final_energy))

def target_func(pos, context, constraints=None):
    context.setPositions(pos.reshape(-1,3))
    state = context.getState(getEnergy=True, getForces=True)
    forces = state.getForces(asNumpy=True)
    energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    forces = forces.value_in_unit(kilojoule_per_mole/nanometer)

    if constraints is not None:
        for n, constr in constraints.items():
            for idx in constr_2_idx[constr.upper()]:
                forces[n][idx] *= 0

    return energy, -forces.flatten()

