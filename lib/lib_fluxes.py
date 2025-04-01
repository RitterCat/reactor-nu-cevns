import numpy as np
from scipy.interpolate import CubicSpline
from lib.lib_constants import *

# the energy per fission for the four main fissile isotopes: u235, u238, pu239, pu241
ENERGY_PER_FISSION_I = np.array([201.92, 205.52, 209.99, 213.60])*MeV

# This function reads the spectrum from a data file
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

# HERE I GET INTERPOLATING FUNCTIONS FOR THE SPECTRA FOR THE FOUR FISSILE ISOTOPES

FISSILE_ISOTOPES_TXTFILE_FORMAT = ['u235', 'u238', 'pu239', 'pu241']
spectrum_source = 'bestiole'

fission_spectra_data = []
fission_spectra = []

for fi in FISSILE_ISOTOPES_TXTFILE_FORMAT:
    spec = get_spectrum(f"./fluxData/{spectrum_source}_{fi}.txt") # Locally, use ./ and on spartan, use ../
    fission_spectra_data.append(spec)

    fission_spectra.append(CubicSpline(*spec, extrapolate=False))

# HERE I GET INTERPOLATING FUNCTIONS FOR THE SPECTRA FOR THE FOUR BREEDING ISOTOPES

BREEDING_ISOTOPES_TXTFILE_FORMAT = ['u239', 'np239']

breeding_spectra_data = []
breeding_spectra = []

for bi in BREEDING_ISOTOPES_TXTFILE_FORMAT:
    spec = get_spectrum(f"./fluxData/{spectrum_source}_{bi}.txt") # Locally, use ./ and on spartan, use ../
    breeding_spectra_data.append(spec)

    breeding_spectra.append(CubicSpline(*spec, extrapolate=False))

# REACTOR FLUX CALCULATIONS

# This functions gets the neutrino flux (/s/cm^2) coming from the reactor, given the fission rate per isotope
def reactor_flux(Enu, fission_rate_per_isotope, breeding_rate=0):
    # This gives the neutrino flux (#nu/s) only arising from the four fissile isotopes
    nu_rate = sum([fission_rate*spec(Enu) for fission_rate, spec in zip (fission_rate_per_isotope, fission_spectra)])

    # If there is a non-zero breeding rate, then we add on the neutrino flux arising from the two breeding beta decays
    if breeding_rate != 0:
        breeding_nu_rate = sum([breeding_rate*spec(Enu) if not np.isnan(spec(Enu)) else 0 for spec in breeding_spectra])
        nu_rate += breeding_nu_rate

    return nu_rate

# This functions gets the neutrino flux (/s/cm^2) coming from the reactor, 
# given the total thermal power and the fraction of fissions due to each fissile isotope
def get_reactor_flux_Pth_fuel_fractions(thermal_power, fuel_fractions):
    fuel_fractions = np.array(fuel_fractions) # make sure the fuel fractions are an array

    # then calculate true fission rate per isotope, which is what we input into the rate calculation
    mean_energy_per_fission = sum(fuel_fractions*ENERGY_PER_FISSION_I)
    total_fission_rate = thermal_power/mean_energy_per_fission
    fission_rate_per_isotope = total_fission_rate*fuel_fractions

    return reactor_flux(fission_rate_per_isotope)

##### OLD FUNCTIONS #####

# def reactor_flux_old(Enu, fuel_fractions, thermal_power):
#     norm_fuel_fractions = np.array(fuel_fractions)/sum(fuel_fractions)

#     mean_energy_per_fission = sum(norm_fuel_fractions*ENERGY_PER_FISSION_I)

#     return (thermal_power/mean_energy_per_fission)*sum([fi*spec(Enu) for fi, spec in zip (norm_fuel_fractions, fission_spectra)])

# def reactor_flux(Enu, fission_rate_per_isotope):
#     return sum([fission_rate*spec(Enu) for fission_rate, spec in zip (fission_rate_per_isotope, fission_spectra)])
