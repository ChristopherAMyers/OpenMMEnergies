import openmm.unit as unit
from openmm.openmm import CustomBondForce, CustomNonbondedForce, NonbondedForce
from itertools import combinations

def determine_fragment_index(topol, res_ids):
    fragment_index_list = []

    #   first make a dictionary to map residue id to it's Residue object
    residues = {}
    for res in topol.residues():
        residues[int(res.id)] = res

    #   now grab the atom index for all atoms in each fragment
    for name, id_list in res_ids.items():
        index_list = []
        for res_id in id_list:
            for atom in residues[res_id].atoms():
                index_list.append(atom.index)
        fragment_index_list.append(tuple(index_list))
    
    #   print summary
    print("\n A fragment interaction calculation has been requested for {:d} fragments".format(len(res_ids)))
    print(" ---------------------------------------------")
    for i, (name, res_list) in enumerate(res_ids.items()):
        print(" Fragment {:d}".format(i + 1))
        print("     Name:          {:s}".format(name))
        print("     Num. Residues: {:d}".format(len(res_list)))
        print("     Num. Atoms:    {:d}".format(len(fragment_index_list[i])))
    print(" ---------------------------------------------\n")

    return fragment_index_list

def create_decomposed_forces(system, topol=None, exclude_self=False, fragments=None):
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

    #   determine if fragments should be used
    fragment_index_list = {}
    if topol is not None and fragments is not None:
        fragment_index_list = determine_fragment_index(topol, fragments)

    print("\n Interaction energy decomposition is requested")

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

    #   add the exceptions
    for exception in exceptions:
        p1, p2, chg_prod, sigma, eps = exception

        #   only update forces if there are no fragments being used
        #   or if each exception bemong to different fragments
        add_to_force = True
        if len(fragment_index_list) != 0:
            p1_frag = -1
            p2_frag = -1
            for i, index_list in enumerate(fragment_index_list):
                if p1 in index_list:
                    p1_frag = i
                if p2 in index_list:
                    p2_frag = i
            if p1_frag == p2_frag:
                add_to_force = False

        if add_to_force:
            chg_except_force.addBond(p1, p2, [chg_prod])
            disp_except_force.addBond(p1, p2, [sigma, eps])
            chg_force.addExclusion(p1, p2)
            disp_force.addExclusion(p1, p2)

    if len(fragment_index_list) != 0:
        for p1, p2 in combinations(range(len(fragment_index_list)), 2):
            chg_force.addInteractionGroup(fragment_index_list[p1], fragment_index_list[p2])
            disp_force.addInteractionGroup(fragment_index_list[p1], fragment_index_list[p2])

    system.removeForce(nbd_force_idx)
    system.addForce(chg_force)
    system.addForce(disp_force)
    system.addForce(chg_except_force)
    system.addForce(disp_except_force)
    print("     Total charge of system: {:.5f}".format(total_charge))