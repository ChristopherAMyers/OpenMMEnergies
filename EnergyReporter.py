import sys
from openmm.app.forcefield import CMAPTorsion, CustomTorsion, PME
from openmm import openmm as mm
import openmm.unit as unit
import json
import os
import numpy as np


class EnergyReporter(object):
    
    """Custom state reporter to print out energy components of a force field"""
    def __init__(self, interval, system, file=sys.stdout, print_eda=True, nonbonded_eda=False, compute_forces=0, keep_results=False):
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
        self._force_groups = self._get_force_groups(system)
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
        self._group2name = {}
        self._energy_terms = {}
        self._print_eda = print_eda
        self._nonbonded_eda = nonbonded_eda
        self._compute_forces = int(compute_forces)

        self._glob2name = {
            'EDA_chg': 'Electrostatic',
            'EDA_pauli': 'Pauli 12 power',
            'EDA_disp': 'Dispersion',
            'EDA_chg_14': 'Electrostatic 1-4',
            'EDA_pauli_14': 'Pauli 12 power 1-4',
            'EDA_disp_14': 'Dispersion 1-4',
        }

        #   saves all data to file
        self._data_files = {'energies': {}, 'forces': {}}
        self._opened_files = False
        self._all_saved_data = []
        self._keep_results = keep_results
        if self._keep_results:
            os.makedirs('energy_report', exist_ok=True)
            
    def __del__(self):
        self._out.close()
        if self._data_files['energies']:
            self._data_files['energies'].close()
        for file in self._data_files['forces'].values():
            file.close()

    def _create_data_files(self, results):
        #   create files if they haven't been created already
        if self._keep_results and not self._opened_files:
            os.makedirs('energy_report', exist_ok=True)
            file = open('energy_report/{:s}.txt'.format('energies'), 'w')
            self._data_files['energies'] = file
            file.write("#")
            for name in results['energies'].keys():
                file.write('{:21s}  '.format(name.replace(' ', '_')))
            file.write('\n')

            self._data_files['forces'] = {}
            for force_name in results['forces']:
                file = open('energy_report/force_{:s}.txt'.format(force_name), 'w')
                self._data_files['forces'][force_name] = file

            self._opened_files = True

        if self._keep_results:
            for energy_name in results['energies']:
                energy = results['energies'][energy_name]
                self._data_files['energies'].write('{:21.14e}  '.format(energy))
            self._data_files['energies'].write('\n')
            for force_name, file in self._data_files['forces'].items():
                force = results['forces'][force_name]
                np.savetxt(file, force)
                file.write('\n')
                
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

    def report(self, simulation, state_in=None, total_only=False):
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

        #   try to get information from provided state
        saved_data = {'energies': {}, 'forces': {}}
        state_args = {}
        if state_in:
            try:
                total = state_in.getPotentialEnergy()
            except:
                state_args['getEnergy'] = True

            if self._compute_forces >= 1:
                try:
                    total_forces = state_in.getForces(asNumpy=True)
                except:
                    state_args['getForces'] = True
            if len(state_args) > 0:
                state = simulation.context.getState(**state_args)
                if 'getEnergy' in state_args:
                    total = state.getPotentialEnergy()
                if 'getForces' in state_args:
                    total_forces = state.getForces(asNumpy=True)
                
        #   energy decomposition
        if self._print_eda:
            print("")
            print(" --------------------------------------------- ")
            print("             Energy Decomposition              ")
            print(" --------------------------------------------- ")
            total = 0 * unit.kilojoule_per_mole
            for f, i in self._force_groups.items():
                
                state = simulation.context.getState(getEnergy=True, getForces=self._compute_forces, groups={i})
                #state = simulation.context.getState(getEnergy=True, getForces=self._compute_forces, groups=2**i)
                energy = state.getPotentialEnergy()
                energy_name = self._type2name.get(type(f), str(f))

                if energy_name in saved_data['energies']:
                    counter = 2
                    while True:
                        test_name = energy_name + " " + str(counter)
                        if test_name in saved_data:
                            counter += 1
                        else:
                            energy_name = test_name
                            break
                
                #energy_type = self._force_labels[f]
                self._energy_terms[i] = {}
                self._energy_terms[i]['energy'] = energy
                saved_data['energies'][energy_name] = energy/unit.kilojoule_per_mole

                print(" {:22} {:15.3f} kJ/mol".format(energy_name, energy/unit.kilojoule_per_mole))
                total += energy

                if self._compute_forces in [2]:
                    saved_data['forces'][energy_name] = state.getForces(asNumpy=True)


        print(" --------------------------------------------- ")
        print(" {:22} {:15.3f} kJ/mol".format("Total energy", total/unit.kilojoule_per_mole))
        print(" --------------------------------------------- ")

        if self._nonbonded_eda:
            state = simulation.context.getState(getParameterDerivatives=True)
            param_derivs = state.getEnergyParameterDerivatives().asdict()
            total = 0
            print("")
            print(" --------------------------------------------- ")
            print("        Nonbonded Energy Decomposition         ")
            print(" --------------------------------------------- ")
            for key, value in param_derivs.items():
                energy_type = self._glob2name.get(key, key)
                self._energy_terms[energy_type] = {}
                self._energy_terms[energy_type]['energy'] = value
                saved_data['energies'][energy_type] = value
                print(" {:22} {:15.3f} kJ/mol".format(energy_type, value))
                total += value

            print(" --------------------------------------------- ")
            print(" {:22} {:15.3f} kJ/mol".format("Total Nonbonded energy", total))
            print(" --------------------------------------------- ")

            if self._compute_forces == 3:

                #   determine which groups to compute
                groups = set()
                for force, i in self._force_groups.items():
                    if isinstance(force, mm.CustomNonbondedForce):
                        groups.add(i)
                
                #   not compute with all other parameters set to zero
                for main_param in param_derivs:

                    #   skip 1-4 interactions
                    energy_type = self._glob2name.get(main_param, main_param)
                    if '1-4' in energy_type:
                        continue

                    #   set all other parameters to zero except main one
                    simulation.context.setParameter(main_param, 1.0)
                    for other_param in param_derivs:
                        if main_param != main_param:
                            simulation.context.setParameter(other_param, 0.0)

                    #   now compute the forces, which should only with proportional
                    #   to the main_param
                    state = simulation.context.getState(getForces=True, groups=groups)
                    saved_data['forces'][energy_type] = state.getForces(asNumpy=True)
                    print("FORCES: ", energy_type)

        if self._keep_results:
            saved_data['energies']['total'] = total
            if self._compute_forces:
                saved_data['forces']['total'] = total_forces

            self._create_data_files(saved_data)


