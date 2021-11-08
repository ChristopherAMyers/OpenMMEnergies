import openmm.unit as unit
from openmm.openmm import CustomBondForce, CustomNonbondedForce, NonbondedForce
from copy import deepcopy
from openmm.unit.quantity import Quantity
from math import sqrt

def create_decomposed_forces(system, topol=None, exclude_self=False):
    ''' Replace Nonbonded Force with separate forces for electrostatics
        and Lennard-Jones. This will add global parameter derivatives for 
        each interaction type. Each global parameter starts with the "EDA_"
        and multiplies each energy term.

        Parameters
        ----------
        system: OpenMM System object
            the system being calculated
    '''
    nbforce = None
    nbd_force_idx = None
    for n, force in enumerate(system.getForces()):
        if isinstance(force, NonbondedForce):
            nbforce = force
            nbd_force_idx = n

    #   excract original parameters
    orig_params = [nbforce.getParticleParameters(n) for n in range(nbforce.getNumParticles())]
    exceptions = [nbforce.getExceptionParameters(n) for n in range(nbforce.getNumExceptions())]
    exception_pairs = [(min(x[0:2]), max(x[0:2])) for x in exceptions]



    #   add new force using same interaction terms
    chg_energy_str = "EDA_chg*138.935458*chg1*chg2/r;"
    chg_force = CustomNonbondedForce(chg_energy_str)
    chg_force.addPerParticleParameter('chg')
    chg_force.addGlobalParameter('EDA_chg', 1.0)
    chg_force.addEnergyParameterDerivative('EDA_chg')
    disp_energy_str = "4*epsilon*(EDA_pauli*(sigma/r)^12-EDA_disp*(sigma/r)^6); \
                                    sigma=0.5*(sigma1+sigma2); \
                                    epsilon=sqrt(epsilon1*epsilon2)"
    disp_force = CustomNonbondedForce(disp_energy_str)
    disp_force.addPerParticleParameter('sigma')
    disp_force.addPerParticleParameter('epsilon')
    disp_force.addGlobalParameter('EDA_pauli', 1.0)
    disp_force.addGlobalParameter('EDA_disp', 1.0)
    disp_force.addEnergyParameterDerivative('EDA_pauli')
    disp_force.addEnergyParameterDerivative('EDA_disp')

    #   used for particle pair exceptions
    chg_except_force = CustomBondForce("EDA_chg_14*138.935458*q/r")
    chg_except_force.addPerBondParameter('q')
    chg_except_force.addGlobalParameter('EDA_chg_14', 1.0)
    chg_except_force.addEnergyParameterDerivative('EDA_chg_14')
    disp_except_force = CustomBondForce("4*epsilon*(EDA_pauli_14*(sigma/r)^12-EDA_disp_14*(sigma/r)^6)")
    disp_except_force.addPerBondParameter('sigma')
    disp_except_force.addPerBondParameter('epsilon')
    disp_except_force.addGlobalParameter('EDA_disp_14', 1.0)
    disp_except_force.addGlobalParameter('EDA_pauli_14', 1.0)
    disp_except_force.addEnergyParameterDerivative('EDA_disp_14')
    disp_except_force.addEnergyParameterDerivative('EDA_pauli_14')
    
    #   add all particles to the nonbonded forces
    total_charge = 0
    for n in range(nbforce.getNumParticles()):
        chg, sigma, eps = orig_params[n]
        disp_force.addParticle([sigma, eps])
        chg_force.addParticle([chg])
        total_charge += chg/unit.elementary_charge

    if topol is not None and exclude_self:
        print(" Excluding inter-residue nonbonded energies")
        for res in topol.residues():
            for atom1 in res.atoms():
                for atom2 in res.atoms():
                    if atom2.index <= atom1.index: continue
                    p1, p2, chg_prod, sigma, eps = atom1.index, atom2.index, 0.0, 0.0, 0.0
                    exception = (p1, p2, chg_prod, sigma, eps)
                    pair = (min(p1, p2), max(p1, p2))
                    if pair not in exception_pairs:
                        exceptions.append(exception)
                        exception_pairs.append(pair)
                    else:
                        except_idx = exception_pairs.index(pair)
                        exceptions[except_idx] = exception
                        
                    #chg_force.addExclusion(atom1.index, atom2.index)
                    #disp_force.addExclusion(atom1.index, atom2.index)

    #   add the exceptions
    for exception in exceptions:
        p1, p2, chg_prod, sigma, eps = exception
        chg_except_force.addBond(p1, p2, [chg_prod])
        disp_except_force.addBond(p1, p2, [sigma, eps])

        chg_force.addExclusion(p1, p2)
        disp_force.addExclusion(p1, p2)

    system.removeForce(nbd_force_idx)
    system.addForce(chg_force)
    system.addForce(disp_force)
    system.addForce(chg_except_force)
    system.addForce(disp_except_force)
    print(" Total charge: ", total_charge)