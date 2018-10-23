import glob
import pandas
import numpy as np
import argparse
from scipy.optimize import curve_fit, minimize_scalar
import matplotlib.pyplot as plt

cm_to_m = 0.01

def min_mult(params):
    res = minimize_scalar(poly, args=((a, b, c)), method='bounded', bounds=(0.1,2))
    
    return res.x

def load_data(file):
    """Load the data from text file
    """

    data = pandas.read_csv(file)

    return data

def power(x, a, b):
    """Exponential fitting function
    """

    return a*np.power(x, b)

def poly(x, a, b, c):

    return np.add(np.add(np.power(x,a), np.multiply(x, b)), c)

def fit_data(data, func, x, y):
    """Fit the data to function
    """    
    popt, pcov = curve_fit(func, data[x], data[y])

    return popt, pcov

def plot_results(data):
    """
    """
    fig, ax = plt.subplots()


    line_formats = {'CO2 UO2 ' : 'r--',
                    'H2O UO2 ' : 'r-',
                    'CO2 UN '  : 'b--',
                    'H2O UN '  : 'b-'}

    for rxtr in data:
        res = data[rxtr]
        ax.plot(res[0], res[1], line_formats[rxtr], label=rxtr)

    plt.legend()

    plt.title('Critical Radius: Buried Reactor on Mars')
    plt.xlabel('Fuel Fraction [-]')
    plt.ylabel('Critical Core Radius [cm]')
    plt.savefig('crit_radius.png', dpi=700)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", '--reflector', action='store_true', default=False, help="reflector mode")
    parser.add_argument("-c", '--crit_radius', action='store_true', default=False,
            help="radius mode")
    args = parser.parse_args()

    
    res_files = glob.glob('./*results.txt')
    plot_data = {}
    for file in res_files:
        data = load_data(file)
        name = " ".join(file.split('results.txt')[0].split('_')).strip('./')
        if args.reflector:
            popt, cov = fit_data(data, poly, 'ref_mult', 'total_mass')
            x = min_mult(popt)
            print(name + ' {0:.3f}'.format(x))
            print(name + ' ' + str(popt[0]) + ' ' + str(popt[1]))

        elif args.crit_radius:
            plot_data[name] = (data['fuel_frac'], data['crit_radius'])
            popt, cov = fit_data(data, power, 'fuel_frac', 'crit_radius')
            print(name + ' ' + str(popt[0]) + ' ' + str(popt[1]))

    if args.crit_radius:
        plot_results(plot_data)
