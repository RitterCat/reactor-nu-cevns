from dis import dis
from xml.etree.ElementTree import ElementTree
from lib.lib_rates import *

def table_II():
    # replicate Table II

    thresholds = np.array([0, 10, 100, 1000])*eV

    elements = ["Si", "Ge", "Xe"]

    fuel_fractions = [1, 0, 0, 0]
    power = 0.1*GIGAWATT
    distance = 10*METER


    print('Threshold [eV] | 0  10  100  1000')
    for element in elements:
        CEvNS_rates = [total_CEvNS_rate(threshold, element, fuel_fractions, power, distance)*KILOGRAM*YEAR for threshold in thresholds]
        print(f'{element} | ' + ' '.join([f'{rate:.0f}' for rate in CEvNS_rates]))

# CALCULATING THE IBD RATE

# the IBD threshold energy

ENU_IBD_THRESHOLD = 1.806*MeV
IBD_PHASE_SPACE_FRACTION = 1.7152
ELECTRON_MASS = 0.510999*MeV
NEUTRON_LIFETIME = 878.4*SECOND
mCH2 = 14.02658*AMU

NEUTRON_MASS = 1.008665*AMU
PROTON_MASS = 1.00727647*AMU

def total_IBD_rate(fuel_fractions, thermal_power, L):
    # This will be the IBD rate for a CH2 detector (i.e. using the target mass as the mass of a CH2 molecule)
    
    flux = lambda Enu: reactor_flux(Enu, fuel_fractions, thermal_power)

    mT = mCH2

    flux_norm = 1/(4*np.pi*L**2)

    def integrand(Enu):
        Ee = Enu - NEUTRON_MASS + PROTON_MASS

        pe = (Ee**2 - ELECTRON_MASS**2)**0.5

        sigma_IBD = (2*np.pi/(ELECTRON_MASS**5 * IBD_PHASE_SPACE_FRACTION * NEUTRON_LIFETIME)) * Ee * pe

        return 2*(flux_norm/mT)*flux(Enu)*sigma_IBD # the factor of 2 is for the two protons in CH2

    return quad(integrand, ENU_IBD_THRESHOLD, FLUX_ENU_MAX)[0]

if __name__ == "__main__":
    BH_rates = [418, 636, 288, 398]

    for i in range(4):
        fuel_fractions = np.zeros(4)
        fuel_fractions[i] = 1
        IBD_rate = total_IBD_rate(fuel_fractions, 0.1*GIGAWATT, 10*METER)*KILOGRAM*YEAR
        print(IBD_rate, IBD_rate/BH_rates[i])