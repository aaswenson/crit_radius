import sys
import os
import numpy as np
from scipy.optimize import minimize_scalar, curve_fit
from subprocess import call, DEVNULL
from mcnp_inputs import HomogeneousInput
import fit_data as fd

target_keff = 1.01
# critical radius range to sweep
domain = (5, 28.75)

g_to_kg = 0.001

def calc_keff(config):
    """Calculate keff deviation from target.
    """
    frac = config['fuel_frac']
    radius = config['core_r']

    basename = "{0}_{1}.i".format(round(frac,5), round(radius,5))
    write_inp(basename, config)
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

def write_inp(basename, configuration):
    """Write mcnp input for keff.
    """
    input = HomogeneousInput(config=configuration)
    homog_comp = input.homog_core()
    input.write_input(basename)
    
def cubic(x, a, b, c, d):
    """
    """
    A = a*np.power(x,3)
    B = b*np.power(x,2)
    C = c*x

    return np.add(A,np.add(B,np.add(C,d)))

def quad(x, a, b, c):
    """
    """
    A = a*np.power(x,2)
    B = b*x

    return np.add(A,np.add(B,c))

def crit_radius(config, steps):
    """Get critical radius = f(fuel_frac)
    """

    coolant = config['cool']
    fuel = config['fuel']
    name = '{0}_{1}'.format(coolant, fuel)
    
    x = []
    keff = []
    
    resfile = open('{0}_fits.csv'.format(name), '+a')
    resfile.write('fuel fraction: {0:.2f}\n'.format(config['fuel_frac']))
    
    mid_r = (domain[1] + domain[0]) / 2
    config['core_r'] = mid_r
    kmid, stdvmid = calc_keff(config)
    
    if kmid > target_keff:
        lower_r = domain[0]
        upper_r = mid_r
        radii = np.linspace(lower_r, upper_r, steps)[:-1]
    else:
        lower_r = mid_r
        upper_r = domain[1]
        radii = np.linspace(lower_r, upper_r, steps)[1:]
        x.append(lower_r)
        keff.append(kmid)
        resfile.write('{0},{1},{2}\n'.format(lower_r, kmid, stdvmid))

    for radius in radii:
        config['core_r'] = radius
        res, stdv = calc_keff(config)
        x.append(radius)
        keff.append(res)
        resfile.write('{0},{1},{2}\n'.format(radius, res, stdv))

    if kmid > target_keff:
        x.append(upper_r)
        keff.append(kmid)
        resfile.write('{0},{1},{2}\n'.format(upper_r, kmid, stdvmid))
    
    keff_reduced = np.array(keff) - target_keff
    coeffs, copt = curve_fit(quad, x, keff_reduced)
    
    roots = np.roots(coeffs)
    # return only the real root
    real_roots = roots[np.isreal(roots)].real
    real_roots = [x for x in real_roots if x >= 0]
    

    return min(real_roots)
    
def fuel_frac(coolant, fuel, clad, matr, frac=None):
    """Determine the optimal reflector thickness for a given reactor
    configuration.
    """
    rhos = {'CO2' : 252.638e-3, 'H2O' : 141.236e-3}
    config = {'fuel' : fuel,
              'matr' : matr,
              'cool' : coolant,
              'clad' : clad,
              'rho_cool' : rhos[coolant],
             }
    resname = '{0}_{1}_results.txt'.format(coolant, fuel)
    results = open(resname, '+w')
    results.write('fuel_frac,crit_radius\n') 
    results.close()

    steps = 5
    start = 0.6
    stop = 1
    # allow for 1 frac to speed up calc
    if frac:
        steps = 1
        start = frac
        stop = frac
    for frac in np.linspace(start, stop, steps):
        results = open(resname, 'a')
        config['fuel_frac'] = frac
        config['ref_mult'] = refl_mult(config)
        # get critical radius
        r = crit_radius(config, 5)
        results.write('{0:.2f},{1:.5f},{2:.5f}\n'.format(frac, r,
                                                         config['ref_mult']))
        results.close()

def refl_mult(config):
    """Determine the optimal reflector thickness for a given reactor
    configuration.
    """

    mults = np.linspace(0.001, 0.4, 5)
    data = {'mass' : [], 'mult' : mults}
    refl_res = open('refl_results.txt', 'a')
    for mult in mults:
        config['ref_mult'] = mult
        # get critical radius
        r = crit_radius(config, 5)
        # get the masses
        config['core_r'] = r
        input = HomogeneousInput(config=config)
        homog_comp = input.homog_core()
        tot_mass = input.core_mass + input.refl_mass + input.PV_mass
        data['mass'].append(tot_mass)
        # save results
        refl_res.write('{0},{1},{2},{3},{4},{5}\n'.format(config['fuel_frac'], 
                                                          config['ref_mult'], 
                                                          input.core_mass,
                                                          input.refl_mass,
                                                          input.PV_mass,
                                                          tot_mass))

    popt, cov = fd.fit_data(data, fd.poly, 'mult', 'mass')
    opt_mult = fd.min_mult(popt, fd.poly)
    refl_res.close()
    
    return opt_mult

if __name__ == '__main__':
    cool, fuel, clad, matr = sys.argv[1:5]
    if matr == 'None':
        matr = None
    if len(sys.argv) > 5:
        frac = sys.argv[5]
        fuel_frac(cool, fuel, clad, matr,float(frac))
    else:
        fuel_frac(cool, fuel, clad, matr,frac)

