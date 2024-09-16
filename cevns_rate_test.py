###########
# IMPORTS #
###########

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
import matplotlib.pyplot as plt
import pymultinest

#############
# CONSTANTS #
#############

# fundamental units and physical constants
GeV = 1
MeV = 1e-3*GeV
keV = 1e-6*GeV
eV = 1e-9*GeV
HBARC = 0.19732698 # GeV fm
SPEED_OF_LIGHT = 2.997925e8 # m/s
G_FERMI = 1.166378e-5/GeV**2
ELECTRON_CHARGE = 1.60217663e-19 # coulombs
AMU = 0.931494

# derived units, all in terms of GeV
METER = 1e15/HBARC
CENTIMETER = 1e-2*METER
SECOND = SPEED_OF_LIGHT*METER
YEAR = 60*60*24*365.25*SECOND
JOULE = eV/ELECTRON_CHARGE
GIGAWATT = 1e9*JOULE/SECOND
KILOGRAM = SPEED_OF_LIGHT**2*JOULE
KILOGRAMDAY = KILOGRAM*(60*60*24*SECOND)
CM2S = CENTIMETER**2*SECOND

# isotope properties for different detector materials
# isotopes of xenon
xe124 = {"Z": 54, "A": 124, "abundance": 0.00095, "mass": 123.906*AMU}
xe126 = {"Z": 54, "A": 126, "abundance": 0.00089, "mass": 125.904*AMU}
xe128 = {"Z": 54, "A": 128, "abundance": 0.0191, "mass": 127.904*AMU}
xe129 = {"Z": 54, "A": 129, "abundance": 0.26401, "mass": 128.905*AMU}
xe130 = {"Z": 54, "A": 130, "abundance": 0.04071, "mass": 129.904*AMU}
xe131 = {"Z": 54, "A": 131, "abundance": 0.21232, "mass": 130.905*AMU}
xe132 = {"Z": 54, "A": 132, "abundance": 0.26909, "mass": 131.904*AMU}
xe134 = {"Z": 54, "A": 134, "abundance": 0.10436, "mass": 133.905*AMU}
xe136 = {"Z": 54, "A": 136, "abundance": 0.08857, "mass": 135.907*AMU}
isotopesXe = [xe124, xe126, xe128, xe129, xe130, xe131, xe132, xe134, xe136]
mXe = sum([iso["abundance"]*iso["mass"] for iso in isotopesXe])

# isotopes of silicon
si28 = {"Z": 14, "A": 28, "abundance": 0.92223, "mass": 27.977*AMU}
si29 = {"Z": 14, "A": 29, "abundance": 0.04685, "mass": 28.977*AMU}
si30 = {"Z": 14, "A": 30, "abundance": 0.03092, "mass": 29.974*AMU}
isotopesSi = [si28, si29, si30]
mSi = sum([iso["abundance"]*iso["mass"] for iso in isotopesSi])

#isotopes of germanium
ge70 = {"Z": 32, "A": 70, "abundance": 0.2052, "mass": 69.924*AMU}
ge72 = {"Z": 32, "A": 72, "abundance": 0.2745, "mass": 71.922*AMU}
ge73 = {"Z": 32, "A": 73, "abundance": 0.0776, "mass": 72.923*AMU}
ge74 = {"Z": 32, "A": 74, "abundance": 0.3652, "mass": 73.921*AMU}
ge76 = {"Z": 32, "A": 76, "abundance": 0.0775, "mass": 75.921*AMU}
isotopesGe = [ge70, ge72, ge73, ge74, ge76]
mGe = sum([iso["abundance"]*iso["mass"] for iso in isotopesGe])

ISOTOPES = {"Xe": isotopesXe, "Si": isotopesSi, "Ge": isotopesGe}

# function to get the target mass (i.e. average mass of isotopes) for a given element
def mTarget(element):
    return sum([iso["abundance"]*iso["mass"] for iso in ISOTOPES[element]])

# weak nuclear charges of protons and neutrons
Qp = 0.0747
Qn = -1.0235

# the energy per fission for the four main fissile isotopes: u235, u238, pu239, pu241
ENERGY_PER_FISSION_I = np.array([201.92, 205.52, 209.99, 213.60])*MeV

##############
# KINEMATICS #
##############

def get_Enu_min(ER, mT):
    return (ER*mT/2)**0.5

def get_ER_max(Enu, mT):
    return 2*Enu**2/(mT + 2*Enu)

#######################
# CEvNS CROSS SECTION #
#######################

def Fhelm(q, A):
    s = 0.9
    r = ( (1.23*A**(1/3) - 0.6)**2 + (7/3)*np.pi**2*0.52**2 - 5*s**2 )**0.5
    return 3 * ( np.sin(q*r/HBARC) - (q*r/HBARC)*np.cos(q*r/HBARC) ) / ( q*r/HBARC )**3 * np.exp( -1*(q*s/HBARC)**2/2 )

def dsigma_dER(ER, Enu, isotope):
    mT = isotope["mass"]
    Z = isotope["Z"]
    A = isotope["A"]

    prefactor = (G_FERMI**2 * isotope["mass"]) / (4*np.pi)
    kinematic_factor = 1 - (ER/Enu) - (ER*mT)/(2*Enu**2)
    weak_charge = Z*Qp + (A-Z)*Qn
    q = (2*mT*ER)**0.5

    return prefactor * kinematic_factor * weak_charge**2 * Fhelm(q, A)

###################
# NEUTRINO FLUXES #
###################

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

fissile_isotopes_txtfile_format = ['u235', 'u238', 'pu239', 'pu241']
spectrum_source = 'bestiole'

fission_spectra_data = []
fission_spectra = []

for fi in fissile_isotopes_txtfile_format:
    spec = get_spectrum(f"fluxData/{spectrum_source}_{fi}.txt")
    fission_spectra_data.append(spec)

    fission_spectra.append(CubicSpline(*spec, extrapolate=False))

def reactor_flux(Enu, fuel_fractions, thermal_power):
    norm_fuel_fractions = np.array(fuel_fractions)/sum(fuel_fractions)

    mean_energy_per_fission = sum(norm_fuel_fractions*ENERGY_PER_FISSION_I)

    return (thermal_power/mean_energy_per_fission)*sum([fi*spec(Enu) for fi, spec in zip (norm_fuel_fractions, fission_spectra)])

##############
# CEvNS RATE #
##############

def dR_dER(ER, detector_material, fuel_fractions, thermal_power, L):
    
    flux = lambda Enu: reactor_flux(Enu, fuel_fractions, thermal_power)

    flux_Enu_min, flux_Enu_max = FLUX_ENU_MIN, FLUX_ENU_MAX

    mT = mTarget(detector_material)

    flux_norm = 1/(4*np.pi*L**2)

    def integrand(Enu, isotope):
        return flux(Enu) * isotope["abundance"] * dsigma_dER(ER, Enu, isotope)

    def unnormalised_rate_per_isotope(isotope):

        min_endpoint = min( max( get_Enu_min(ER, isotope["mass"]), flux_Enu_min ), flux_Enu_max )
        max_endpoint = flux_Enu_max
        
        integral = quad(integrand, min_endpoint, max_endpoint, args=(isotope))
        return integral[0]

    return (flux_norm/mT)*sum([unnormalised_rate_per_isotope(isotope) for isotope in ISOTOPES[detector_material]])

def total_CEvNS_rate(threshold, detector_material, fuel_fractions, thermal_power, L):
    
    flux_Enu_max = FLUX_ENU_MAX

    mT = mTarget(detector_material)

    return quad(dR_dER, threshold, get_ER_max(flux_Enu_max, mT), args=(detector_material, fuel_fractions, thermal_power, L))[0]

##########
# SCRIPT #
##########

if __name__ =='__main__':

    # specify detector material
    detector_material = 'Xe'
    mT = mTarget(detector_material)
    offset = 10*METER

    # specify reactor propeties (including true fuel fractions)
    thermal_power = GIGAWATT
    true_fuel_fractions = [1,0,0,0]

    false_fuel_fractions = [0,0,0,1]

    # specify bin properties (nbins, ERmin, ERmax)
    nbins = 5
    ER_min = 0.1*keV # this is the threshold value
    ER_max = get_ER_max(FLUX_ENU_MAX, mT)

    bin_width = ( ER_max - ER_min ) / nbins

    true_bin_counts = []
    false_bin_counts = []

    for i in range(nbins):
        true_bin_counts.append(quad(dR_dER, ER_min + i*bin_width, ER_min + (i+1)*bin_width, args=(detector_material, true_fuel_fractions, thermal_power, offset))[0]*KILOGRAM*YEAR)
        false_bin_counts.append(quad(dR_dER, ER_min + i*bin_width, ER_min + (i+1)*bin_width, args=(detector_material, false_fuel_fractions, thermal_power, offset))[0]*KILOGRAM*YEAR)

    print(true_bin_counts)
    print(false_bin_counts)