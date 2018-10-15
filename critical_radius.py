import numpy as np
from string import Template
from scipy.optimize import minimize_scalar
from subprocess import call, DEVNULL
import math
import os
from mcnp_inputs import HomogeneousInput
import argparse

fp = open('base_input.txt', 'r')
lines = fp.read()
fp.close()

temp = Template(lines)
target_keff = 1.01

global core_mass, refl_mass

save_keff = open('keff_res.txt', 'w')
g_to_kg = 0.001

def calc_keff_error(radius, config):
    """Calculate keff deviation from target.
    """
    frac = config['fuel_frac']
    basename = "{0}_{1}.i".format(round(frac,5), round(radius,5))
    write_inp(radius, basename, config)
    call(["mcnp6", "n= {0} tasks 8".format(basename)], stdout=DEVNULL)
    keff = parse_output(basename)
    os.remove('{0}r'.format(basename))
    os.remove('{0}s'.format(basename))
    os.remove('{0}o'.format(basename))
    save_keff.write('radius: {0} frac: {1} keff: {2}\n'.format(radius, frac, keff))

    return abs(keff - target_keff)

def parse_output(basename):
    """Parse the output for keff value.
    """

    fp = open(basename + 'o', 'r')
    lines = fp.readlines()
    fp.close()

    res_idx = []
    for idx, line in enumerate(lines):
        if 'final result' in line:
            res_idx.append(idx)
    keff = float(lines[res_idx[-1]].split()[2])

    return keff

def write_inp(core_r, basename, configuration):
    """Write mcnp input for keff.
    """
    global core_mass, refl_mass

    configuration['core_r'] = core_r
    input = HomogeneousInput(config=configuration)
    homog_comp = input.homog_core()
    input.write_input(basename)
    
    core_mass = input.core_mass
    refl_mass = input.refl_mass

def find_radius(config):
    """Optimize keff error to determine critical radius.
    """
    res = minimize_scalar(calc_keff_error, method='bounded', bounds=(5, 55),
            args=(config), options={'xatol':1e-4})

    return res

def frac_iterate(target, config, range):
    """Get critical radius = f(fuel_frac)
    """

    global core_mass, refl_mass
    coolant = config['cool']
    fuel = config['fuel']

    resfile = open('{0}_{1}_results.txt'.format(coolant, fuel) , '+a')
    
    resfile.write('{0},r_crit,core_mass,refl_mass,total_mass\n'.format(target))
    resfile.close()
    
    for sample in np.arange(range[0], range[1], range[2]):
        config[target] = sample
        res = find_radius(config)
        resfile = open('{0}_{1}_results.txt'.format(coolant, fuel) , '+a')
        resfile.write("{0:.3f},{1:.3f},{2:.3f},{3:.3f},{4:.3f}\n".format(sample, 
                                                                     res.x, 
                                                                     core_mass*g_to_kg,
                                                                     refl_mass*g_to_kg, 
                                                                     (core_mass+refl_mass)*g_to_kg))
        resfile.close()

def optimize_target(coolant, fuel, clad, matr):
    """Determine the optimal reflector thickness for a given reactor
    configuration.
    """
    rhos = {'CO2' : 233.89e-3, 'H2O' : 123.48e-3}
    fracs = [0.3, 0.6, 0.9]
    config = {'fuel' : fuel,
              'matr' : matr,
              'cool' : coolant,
              'clad' : clad,
              'rho_cool' : rhos[coolant],
              'ref_mult' : 0.165
#              'fuel_frac' : 0.6
             }

    frac_iterate('fuel_frac', config, [0.1, 1, 0.1])

if __name__ == '__main__':
    optimize_target('CO2', 'UO2', 'Inconel-718', None)
    optimize_target('H2O', 'UO2', 'Inconel-718', None)
    optimize_target('CO2', 'UN',  'Inconel-718', 'W')
    optimize_target('H2O', 'UN',  'Inconel-718', 'W')
    save_keff.close()
