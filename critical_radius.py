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

def calc_keff_error(radius, config):
    """Calculate keff deviation from target.
    """
    frac = config['fuel_frac']
    basename = "{0}_{1}.i".format(round(frac,5), round(radius,5))
    mass = write_inp(radius, basename, config)
    call(["mcnp6", "n= {0} tasks 8".format(basename)], stdout=DEVNULL)
    keff = parse_output(basename)
    os.remove('{0}r'.format(basename))
    os.remove('{0}s'.format(basename))
    os.remove('{0}o'.format(basename))
#    os.remove(basename)
    
    print("radius: {0:.4f} keff: {1:.3f} mass: {2:.3f}".format(radius, keff, mass/1000))

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
    AR = 1
    L = core_r * AR
    input = HomogeneousInput(core_r, L, config=configuration)
    homog_comp = input.homog_core()
    input.write_input(basename)
    
    return input.core_mass

def find_radius(config):
    """Optimize keff error to determine critical radius.
    """
    res = minimize_scalar(calc_keff_error, method='bounded', bounds=(5, 150),
            args=(config), options={'xatol':1e-4})

    return res

def frac_iterate(coolant, fuel, clad, matr, refl_mult):
    """Get critical radius = f(fuel_frac)
    """
    rhos = {'CO2' : 233.89e-3, 'H2O' : 123.48e-3}
    
    resfile = open('{0}_{1}_results.txt'.format(coolant, fuel) , '+a')
    
    config = {'fuel' : fuel,
              'matr' : matr, 
              'clad' : clad,
              'rho_cool' : rhos[coolant],
              'ref_mult' : refl_mult,
              'fuel_frac' : 0
             }

    for frac in np.arange(0.1, 1.05, 0.05):
        config['fuel_frac'] = frac
        res = find_radius(config)
        resfile.write("{0} {1}\n".format(frac, round(res.x, 5)))
    
    resfile.close()

if __name__ == '__main__':
    frac_iterate('CO2', 'UO2', 'Inconel-718', None, 1.05)
    frac_iterate('H2O', 'UO2', 'Inconel-718', None, 1.05)
    frac_iterate('CO2', 'UN',  'Inconel-718', 'W',  1.05)
    frac_iterate('H2O', 'UN',  'Inconel-718', 'W',  1.05)
