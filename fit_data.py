import numpy as np
import argparse
from scipy.optimize import curve_fit

cm_to_m = 0.01

def load_data(file):
    """Load the data from text file
    """

    lines = open(file, 'r').readlines()

    x = []
    y = []

    for line in lines:
        x.append(float(line.split()[0]))
        y.append(float(line.split()[1])*cm_to_m)

    return x, y

def power(x, a, b):
    """Exponential fitting function
    """

    return a*np.power(x, b)

def fit_data(x, y, func):
    """Fit the data to function
    """

    popt, pcov = curve_fit(func, x, y)

    return popt, pcov

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", type=str, help="datafile")
    args = parser.parse_args()

    x, y = load_data(args.d)

    popt, cov = fit_data(x, y, power)

    print("Fitted Function: {0:.5f}*fuel_frac^({1:.5f})".format(popt[0],
                                                                popt[1]))
