"""Mass fraction compositions for reactor materials.
"""

# fuel fraction in cermet 

rhos = {'UO2': 10.97, 'W' : 19.3, 'UN' : 11.3, 'In' : 8.2}

mats = {'UO2' : {8016  : 1.5204e-01,
                 92000 : 8.4796e-01
                },

        'UN'  : {92000 : 0.94441, 
                 7015 : 0.05559
                },

        'CO2' : {6000 : 2.7291e-01,
                 8016 : 7.2709e-01
                },

        'W'   : {74180 : 1.1746e-03,
                 74182 : 2.6227e-01,
                 74183 : 1.4241e-01,
                 74184 : 3.0658e-01,
                 74186 : 2.8757e-01
                },
        'In'  : {5010 : 9.2155e-06,
                 5011 : 4.0785e-05,
                 6012 : 7.2154e-04,
                 6013 : 8.4565e-06,
                 13027 : 5.0000e-03,
                 14028 : 2.9214e-03,
                 14029 : 1.5371e-04,
                 14030 : 1.0494e-04,
                 15031 : 1.4000e-04,
                 16032 : 1.3260e-04,
                 16033 : 1.0797e-06,
                 16034 : 6.3031e-06,
                 16036 : 1.5704e-08,
                 22046 : 7.1281e-04,
                 22047 : 6.5680e-04,
                 22048 : 6.6461e-03,
                 22049 : 4.9790e-04,
                 22050 : 4.8644e-04,
                 24050 : 7.9300e-03,
                 24052 : 1.5903e-01,
                 24053 : 1.8380e-02,
                 24054 : 4.6614e-03,
                 25055 : 3.1800e-03,
                 26054 : 9.5974e-03,
                 26056 : 1.5623e-01,
                 26057 : 3.6726e-03,
                 26058 : 4.9733e-04,
                 27059 : 9.1000e-03,
                 28058 : 3.5279e-01,
                 28060 : 1.4057e-01,
                 28061 : 6.2126e-03,
                 28062 : 2.0133e-02,
                 28064 : 5.2922e-03,
                 29063 : 1.8695e-03,
                 29065 : 8.6052e-04,
                 41093 : 5.1250e-02,
                 42092 : 4.2445e-03,
                 42094 : 2.7310e-03,
                 42095 : 4.7781e-03,
                 42096 : 5.0814e-03,
                 42097 : 2.9569e-03,
                 42098 : 7.5898e-03,
                 42100 : 3.1183e-03
                },
        'H2O' : {1001 : 1.1187e-01,
                 1002 : 2.5713e-05,
                 8016 : 8.8811e-01,
                }
       }


def enrich_fuel(enrich, fuel):
    """Enrich Uranium fuel and mix with specified compound.

    Arguments:
    ----------
        enrich (float): U-235 mass conc.
        fuel (dict) (optional): mass composition of fuel
    Returns:
    --------
        fuel_comp (dict): isotopic mass vector of Uranium fuel compound.
    """
    fuel_comp = {}
    # build fuel comp. starting with bounded compound
    for isotope in mats[fuel]:
        if isotope == 92000:
            # add enriched Uranium to fuel composition
            fuel_comp.update({92235 : mats[fuel][92000]*enrich,
                              92238 : mats[fuel][92000]*(1-enrich) })
        else:
            fuel_comp.update({isotope : mats[fuel][isotope]})
    
    mmO = 32
    mmN = 14.0067
    mmU25 = 235.04
    mmU28 = 238.05

    if fuel == 'UN':
        mmfuel = mmN + 1 / ((enrich / mmU25) + ((1-enrich) / mmU28))
    if fuel == 'UO2':
        mmfuel = mmO + 1 / ((enrich / mmU25) + ((1-enrich) / mmU28))
    
    return fuel_comp, mmfuel
