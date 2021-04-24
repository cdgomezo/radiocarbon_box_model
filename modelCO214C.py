#Author: Carlos Gómez-Ortiz
#Affiliation: Department of Physical Geography and Ecosystem Science, Lund University
#2021

from numpy import *


def modelCO214C(initCond, fluxes):

    # Constants
    M = 5.137e21/1000     # mass of the atmosphere [g]
    Ma = 28.964      # molar mass of air [g/mol]
    MC = 12          # molar mass of carbon [g/mol]
    r = MC*(M/2/Ma)  # conversion from molar fraction to kgC
    dt = 1.
    tau = 0.5
    
    D_ff = -1.

    # Create a list containing just c0
    conc = zeros((len(fluxes['time'])+1, 6))
    conc[0,:] = array((initCond['co2nat']['nh'], initCond['co2nat']['sh'], initCond['co2ff']['nh'], initCond['co2ff']['sh'], 
                       initCond['14c']['nh'], initCond['14c']['sh']))

    for i in range(len(fluxes['time'])):
        
        conc[i+1, 0] = conc[i, 0] + (tau * (-conc[i, 0] + conc[i, 1]) + 
                                     (fluxes['co2']['nh']['bio'][i] + fluxes['co2']['nh']['oce'][i]) / r) * dt
        conc[i+1, 1] = conc[i, 1] + (tau * (conc[i, 0] - conc[i, 1]) + 
                                     (fluxes['co2']['sh']['bio'][i] + fluxes['co2']['sh']['oce'][i]) / r) * dt
        conc[i+1, 2] = conc[i, 2] + (tau * (-conc[i, 2] + conc[i, 3]) + 
                                     fluxes['co2']['nh']['ff'][i] / r) * dt
        conc[i+1, 3] = conc[i, 3] + (tau * (conc[i, 2] - conc[i, 3]) + 
                                     fluxes['co2']['sh']['ff'][i] / r) * dt
        conc[i+1, 4] = conc[i, 4] + (tau * (-conc[i, 4] + conc[i, 5]) + 
                                     (fluxes['14c']['nh']['nuc'][i] + fluxes['14c']['nh']['cosm'][i] + 
                                      D_ff * fluxes['co2']['nh']['ff'][i] + 
                                      fluxes['14c']['nh']['D_atm'][i] * (fluxes['co2']['nh']['bio'][i] + fluxes['co2']['nh']['oce'][i]) +
                                      fluxes['14c']['nh']['biodis'][i] + fluxes['14c']['nh']['ocedis'][i]) / r) * dt                       
        conc[i+1, 5] = conc[i, 5] + (tau * (conc[i, 4] - conc[i, 5]) + 
                                     (fluxes['14c']['sh']['nuc'][i] + fluxes['14c']['sh']['cosm'][i] + 
                                      D_ff * fluxes['co2']['sh']['ff'][i] + 
                                      fluxes['14c']['sh']['D_atm'][i] * (fluxes['co2']['sh']['bio'][i] + fluxes['co2']['sh']['oce'][i]) +
                                      fluxes['14c']['sh']['biodis'][i] + fluxes['14c']['sh']['ocedis'][i]) / r) * dt
    
    obs = {
        'co2nat': {
            'nh': {
                'time': fluxes['time'],
                'value': conc[1:,0].reshape(-1,1)
            },
            'sh': {
                'time': fluxes['time'],
                'value': conc[1:,1].reshape(-1,1)
            }
        },
        'co2ff': {
            'nh': {
                'time': fluxes['time'],
                'value': conc[1:,2].reshape(-1,1)
            },
            'sh': {
                'time': fluxes['time'],
                'value': conc[1:,3].reshape(-1,1)
            }
        },
        '14c': {
            'nh': {
                'time': fluxes['time'],
                'value': conc[1:,4].reshape(-1,1)
            },
            'sh': {
                'time': fluxes['time'],
                'value': conc[1:,5].reshape(-1,1)
            }
        }
    }

    return obs
