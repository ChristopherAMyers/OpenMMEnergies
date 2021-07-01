import numpy as np
from simtk.openmm.app import *
from simtk.unit import *
from simtk.openmm.openmm import CustomNonbondedForce, NonbondedForce, CustomBondForce

# pylint: disable=no-member
import simtk
picoseconds = simtk.unit.picoseconds
picosecond = picoseconds
nanometer = simtk.unit.nanometer
femtoseconds = simtk.unit.femtoseconds
# pylint: enable=no-member

Z2alpha = {
    1: 1.2050,
    6: 0.7146,
    7: 0.7606,
    8: 0.9593,
    15: 3.0127
}

def gaussian_density(system, topology, remove_nbf=True):
        #   find the nonbonded force
    nb_force = None
    nb_force_idx = -1
    for n, force in enumerate(system.getForces()):
        if isinstance(force, NonbondedForce):
            nb_force = force
            nb_force_idx = n
            break

    forceString = "4*epsilon*((sigma/r)^12 - (sigma/r)^6) "
    forceString += " + 138.935458*q1*q2/r; "
    #forceString += " + 138.935458*q1*q2*erf(alpha*r)/r; "
    #forceString += " + 138.935458*( (Z1-q1)*(Z2-q2)*erf(alpha*r)/r - (Z1-q1)*Z2*erf(sqrt(a1)*r)/r - (Z2-q2)*Z1*erf(sqrt(a2)*r)/r + Z1*Z2/r ); "
    forceString += "sigma=0.5*(sigma1+sigma2); "
    forceString += "epsilon=sqrt(epsilon1*epsilon2); "
    forceString += "alpha = sqrt(a1*a2/(a1+a2)); "
    custom_force = CustomNonbondedForce(forceString)
    custom_force.addPerParticleParameter("q")
    custom_force.addPerParticleParameter("sigma")
    custom_force.addPerParticleParameter("epsilon")
    custom_force.addPerParticleParameter("Z")
    custom_force.addPerParticleParameter("a")

    #   extract and scale parameters
    #   also separate out atoms from each fragment
    n_part = nb_force.getNumParticles()
    sigma     = []
    epsilon   = []
    ff_charge = []
    total_charge = 0 * elementary_charge
    for n, atom in enumerate(topology.atoms()):
        params = nb_force.getParticleParameters(n)
        if atom.element.atomic_number > 1:
            params[1] = 1.00 * params[1]
            params[2] = 1.00 * params[2]
        ff_charge.append(params[0])
        sigma.append(params[1])
        epsilon.append(params[2])
        z = int(atom.element.atomic_number)
        alpha = Z2alpha.get(z, 1.00)*(18.8973**2)*200 #   convert to 1/nm^2
        total_charge += ff_charge[n]
        custom_force.addParticle([ff_charge[n], sigma[n], epsilon[n], z, alpha])

    bonds = []
    for bond in topology.bonds():
        bond_idx = (bond.atom1.index, bond.atom2.index)
        if -1 in bond_idx or 15 in bond_idx:
            print(bond_idx, bond) 
        bonds.append(bond_idx)
    custom_force.createExclusionsFromBonds(bonds, 3)
    custom_force.setCutoffDistance(2*nanometer)

    print('adding exceptions:')
    print('Number of exceptions before:', custom_force.getNumExclusions())
    print('Number of nb exceptions to add:', nb_force.getNumExceptions())

    
    bond_expression  = "4*epsilon*((sigma/r)^12 - (sigma/r)^6) + 138.935458*chargeprod/r;"
    custom_bond_force = CustomBondForce(bond_expression)
    custom_bond_force.addPerBondParameter('chargeprod')
    custom_bond_force.addPerBondParameter('sigma')
    custom_bond_force.addPerBondParameter('epsilon')
    for idx in range(nb_force.getNumExceptions()):
        idx, jdx, chg, sigma, epsilon = nb_force.getExceptionParameters(idx)
        if chg/chg.unit != 0 or epsilon/epsilon.unit != 0:
            custom_bond_force.addBond(idx, jdx, [chg, sigma, epsilon])

    if remove_nbf:
        system.removeForce(nb_force_idx)

    system.addForce(custom_force)
    system.addForce(custom_bond_force)
    return custom_force