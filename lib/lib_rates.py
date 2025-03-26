from lib.lib_constants import *
from lib.lib_fluxes import *
from lib.lib_detectors import *
from lib.lib_kinematics import *
import numpy as np
from scipy.integrate import quad 

# weak nuclear charges of protons and neutrons
Qp = 0.0747
Qn = -1.0235

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

def dR_dER(ER, detector_material, L, fission_rate_per_isotope, breeding_rate=0): # fuel_fractions, thermal_power
    
    # flux = lambda Enu: reactor_flux(Enu, fission_rate_per_isotope) # fuel_fractions, thermal_power
    
    flux = lambda Enu: reactor_flux(Enu, fission_rate_per_isotope, breeding_rate) # fuel_fractions, thermal_power

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

def total_CEvNS_rate(threshold, detector_material, L, fission_rate_per_isotope, breeding_rate=0): # fuel_fractions, thermal_power
    
    flux_Enu_max = FLUX_ENU_MAX

    mT = mTarget(detector_material)

    return quad(dR_dER, threshold, get_ER_max(flux_Enu_max, mT), args=(detector_material, L, fission_rate_per_isotope, breeding_rate))[0]