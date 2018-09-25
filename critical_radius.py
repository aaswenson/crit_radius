import numpy as np
from string import Template
from scipy.optimize import minimize, minimize_scalar
from subprocess import call, DEVNULL
import math
import os
from mcnp_inputs import HomogeneousInput
import argparse

fp = open('base_input.txt', 'r')
lines = fp.read()
fp.close()

temp = Template(lines)


def calc_keff_error(radius, frac, coolant, fuel, clad, matr, cool_rho):
    basename = "{0}_{1}.i".format(round(frac,5), round(radius,5))
    mass = write_inp(frac, radius, basename, coolant, fuel, matr, cool_rho)
    call(["mcnp6", "n= {0} tasks 8".format(basename)], stdout=DEVNULL)
    keff = parse_output(basename)
    os.remove('{0}r'.format(basename))
    os.remove('{0}s'.format(basename))
    os.remove('{0}o'.format(basename))
    os.remove(basename)
    
    print("radius: {0:.4f} keff: {1:.3f} mass: {2:.3f}".format(radius, keff, mass/1000))

    return abs(keff - 1.01)

def parse_output(basename):

    fp = open(basename + 'o', 'r')
    lines = fp.readlines()
    fp.close()

    res_idx = []
    for idx, line in enumerate(lines):
        if 'final result' in line:
            res_idx.append(idx)
    keff = float(lines[res_idx[-1]].split()[2])

    return keff

def write_inp(frac, core_r, basename, coolant, fuel, clad, matr, cool_rho):
    """
    """
    AR = 1
    L = core_r * AR
    input = HomogeneousInput(core_r, L, frac)
    homog_comp = input.homog_core(fuel, coolant, clad, matr, 0.93, 0.5, cool_rho)
    input.write_mat_string(homog_comp)
    input.write_input(basename)
    
    return input.core_mass

def find_radius(frac, coolant, fuel, clad, matr, cool_rho):
    
    res = minimize_scalar(calc_keff_error, method='bounded', bounds=(5, 150),
            args=(frac, coolant, fuel, clad, matr, cool_rho), options={'xatol':1e-4})

    return res

def frac_iterate(coolant, fuel, clad, matr):
    rhos = {'CO2' : 233.89e-3, 'H2O' : 123.48e-3}

    resfile = open('{0}_{1}_results.txt'.format(coolant, fuel) , '+a')
    
    for frac in np.arange(0.1, 1.05, 0.05):
        res = find_radius(frac, coolant, fuel, matr, rhos[coolant])
        resfile.write("{0} {1}\n".format(frac, round(res.x, 5)))
    
    resfile.close()

if __name__ == '__main__':
    frac_iterate('CO2', 'UO2', 'In', None)
    frac_iterate('H2O', 'UO2', 'In', None)
    frac_iterate('CO2', 'UN',  'In', 'W')
    frac_iterate('H2O', 'UN',  'In', 'W')
