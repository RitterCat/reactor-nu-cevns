import numpy as np
from scipy.interpolate import CubicSpline
from lib.lib_constants import *

# the energy per fission for the four main fissile isotopes: u235, u238, pu239, pu241
ENERGY_PER_FISSION_I = np.array([201.92, 205.52, 209.99, 213.60])*MeV

def get_spectrum(filename):
    energies = []
    fluxes = []

    for line in open(filename):
        if not line.startswith("#"):
            line=line.strip("\n")
            lineParts=line.split(",")
            energies.append(float(lineParts[0]))
            fluxes.append(float(lineParts[1]))
    
    return (np.array(energies)*MeV, np.array(fluxes)/MeV)

FLUX_ENU_MIN = 12.5*keV
FLUX_ENU_MAX = 12.5*MeV - 12.5*keV

FISSILE_ISOTOPES_TXTFILE_FORMAT = ['u235', 'u238', 'pu239', 'pu241']
spectrum_source = 'bestiole'

fission_spectra_data = []
fission_spectra = []

for fi in FISSILE_ISOTOPES_TXTFILE_FORMAT:
    spec = get_spectrum(f"fluxData/{spectrum_source}_{fi}.txt")
    fission_spectra_data.append(spec)

    fission_spectra.append(CubicSpline(*spec, extrapolate=False))

# def reactor_flux_old(Enu, fuel_fractions, thermal_power):
#     norm_fuel_fractions = np.array(fuel_fractions)/sum(fuel_fractions)

#     mean_energy_per_fission = sum(norm_fuel_fractions*ENERGY_PER_FISSION_I)

#     return (thermal_power/mean_energy_per_fission)*sum([fi*spec(Enu) for fi, spec in zip (norm_fuel_fractions, fission_spectra)])

# def reactor_flux(Enu, fission_rate_per_isotope):
#     return sum([fission_rate*spec(Enu) for fission_rate, spec in zip (fission_rate_per_isotope, fission_spectra)])

def reactor_flux(fission_rate_per_isotope):
    return lambda Enu: sum([fission_rate*spec(Enu) for fission_rate, spec in zip (fission_rate_per_isotope, fission_spectra)])

def get_reactor_flux_Pth_fuel_fractions(thermal_power, fuel_fractions):
    fuel_fractions = np.array(fuel_fractions) # make sure the fuel fractions are an array

    # then calculate true fission rate per isotope, which is what we input into the rate calculation
    mean_energy_per_fission = sum(fuel_fractions*ENERGY_PER_FISSION_I)
    total_fission_rate = thermal_power/mean_energy_per_fission
    fission_rate_per_isotope = total_fission_rate*fuel_fractions

    return reactor_flux(fission_rate_per_isotope)