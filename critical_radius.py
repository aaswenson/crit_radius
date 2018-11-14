import numpy as np
from string import Template
from scipy import interpolate
from scipy.optimize import minimize_scalar, curve_fit
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

opt_refl_mult = {'UN'  : {'CO2' : 0.344, 'H2O' : 0.296},
                 'UO2' : {'CO2' : 0.165, 'H2O' : 0.158}
                }

def calc_keff(config):
    """Calculate keff deviation from target.
    """
    frac = config['fuel_frac']
    radius = config['core_r']

    basename = "{0}_{1}.i".format(round(frac,5), round(radius,5))
    write_inp(radius, basename, config)
    call(["mcnp6", "n= {0} tasks 8".format(basename)], stdout=DEVNULL)
    keff = parse_output(basename)
    os.remove('{0}r'.format(basename))
    os.remove('{0}s'.format(basename))
    os.remove('{0}o'.format(basename))
    os.remove(basename)

    return keff

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
    stdv = float(lines[res_idx[-1]].split()[3])

    return keff, stdv

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

def cubic(x, a, b, c, d):
    """
    """
    A = a*np.power(x,3)
    B = b*np.power(x,2)
    C = c*x

    return np.add(A,np.add(B,np.add(C,d)))

def load_data_from_file(lines):
    """
    """
    results = {}
    for line in lines:
        if 'fuel fraction' in line:
            frac = float(line.split(':')[1])
            results[frac] = {}
        else:
            data = [float(x) for x in line.split(',')]
            results[frac][data[0]] = data[1]
    
    return results

def crit_radius(config, range, load_data=None):
    """Get critical radius = f(fuel_frac)
    """

    coolant = config['cool']
    fuel = config['fuel']
    name = '{0}_{1}'.format(coolant, fuel)
    
    x = []
    keff = []
    
    if load_data:
        for radius in load_data:
            x.append(radius)
            keff.append(load_data[radius])
    else:
        resfile = open('{0}_fits.csv'.format(name), '+a')
        resfile.write('fuel fraction: {0:.2f}\n'.format(config['fuel_frac']))
        for radius in np.linspace(range[0], range[1], range[2]):
            config['core_r'] = radius
            res, stdv = calc_keff(config)
            x.append(radius)
            keff.append(res)
            resfile.write('{0},{1},{2}\n'.format(radius, res, stdv))
        
    keff_reduced = np.array(keff) - target_keff
    coeffs, copt = curve_fit(cubic, x, keff_reduced)
    
    roots = np.roots(coeffs)
    real_roots = roots[np.isreal(roots)].real
    # return only the real root
    print(coolant, fuel, config['fuel_frac'])
    print(real_roots)
    real_roots = [x for x in real_roots if x >= 0]

    return min(real_roots)
    
def fuel_frac(coolant, fuel, clad, matr, load_data=True):
    """Determine the optimal reflector thickness for a given reactor
    configuration.
    """
    rhos = {'CO2' : 252.638e-3, 'H2O' : 141.236e-3}
    config = {'fuel' : fuel,
              'matr' : matr,
              'cool' : coolant,
              'clad' : clad,
              'rho_cool' : rhos[coolant],
              'ref_mult' : opt_refl_mult[fuel][coolant],
             }
    
    results = open('{0}_{1}_results.txt'.format(coolant, fuel) , '+w')
    results.write('fuel_frac,crit_radius\n') 
    
    data = None

    for frac in np.arange(0.1, 1, 0.1):
        config['fuel_frac'] = frac
        # load saved results
        if load_data:
            filename = '{0}_{1}_fits.csv'.format(coolant, fuel)
            lines = open(filename, 'r').readlines()
            keff_data = load_data_from_file(lines)
            data = keff_data[round(frac, 1)]
        # get critical radius
        r = crit_radius(config, [5, 70, 5], data)
        results.write('{0:.2f},{1}\n'.format(frac, r))

    results.close()

def refl_mult(coolant, fuel, clad, matr, load_data=True):
    """Determine the optimal reflector thickness for a given reactor
    configuration.
    """
    rhos = {'CO2' : 252.638e-3, 'H2O' : 141.236e-3}
    config = {'fuel' : fuel,
              'matr' : matr,
              'cool' : coolant,
              'clad' : clad,
              'rho_cool' : rhos[coolant],
              'fuel_frac' : 0.6
             }
    
    results = open('{0}_{1}_results.txt'.format(coolant, fuel) , '+w')
    results.write('fuel_frac,crit_radius\n') 
    
    data = None

    for mult in np.arange(0.1, 5, 0.1):
        config['ref_mult'] = mult
        # load saved results
        if load_data:
            filename = '{0}_{1}_fits.csv'.format(coolant, fuel)
            lines = open(filename, 'r').readlines()
            keff_data = load_data_from_file(lines)
            data = keff_data[round(mult, 1)]
        # get critical radius
        r = crit_radius(config, [5, 70, 5], data)
        results.write('{0:.2f},{1}\n'.format(mult, r))

    results.close()

if __name__ == '__main__':
#   load = True
#   fuel_frac('CO2', 'UO2', 'Inconel-718', None, load)
#   fuel_frac('H2O', 'UO2', 'Inconel-718', None, load)
#   fuel_frac('CO2', 'UN',  'Inconel-718', 'W', load)
#   fuel_frac('H2O', 'UN',  'Inconel-718', 'W', load)
#   save_keff.close()
    # reflector thickness
    load = True
    refl_mult('CO2', 'UO2', 'Inconel-718', None, load)
    refl_mult('H2O', 'UO2', 'Inconel-718', None, load)
    refl_mult('CO2', 'UN',  'Inconel-718', 'W', load)
    refl_mult('H2O', 'UN',  'Inconel-718', 'W', load)
    save_keff.close()
