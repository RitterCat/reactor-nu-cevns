###########
# IMPORTS #
###########

from lib.lib_rates import *
from scipy.stats import poisson, norm
import pymultinest
import configargparse

def get_arguments():
    parser = configargparse.ArgumentParser()

    # Define arguments that can come from command line, config file, or environment variable
    parser.add('-d', '--detector-config', is_config_file=True, 
               help="Path to the configuration file for the properties of the neutrino detector")
    parser.add('--material', type=str, help="Target material of the neutrino detector")
    parser.add('--offset', type=float, help="Offset of the neutrino detector (in meters)")
    parser.add('--threshold', type=float, help="Threshold of the neutrino detector (in keV)")
    parser.add('--exposure', type=float, help="Exposure of the neutrino detector (in kg year)")
    parser.add('--er_max', type=float, default=0, help="The upper recoil energy region of interest (in keV)")
    parser.add('--nbins', type=int, 
               help="""Number of bins between the threshold and the 
               max recoil energy in the neutrino detector""")
    
    parser.add('-r', '--reactor-config', is_config_file=True,
               help="Path to the configuration file for the properties of the nuclear reactor")
    parser.add('--thermal_power', type=float, help="Reactor power (in gigawatt)")
    parser.add('--true_fuel_fractions', type=str, help="A comma-delimited string of the true fuel fractions")

    # Parse the command line arguments
    args = parser.parse_args()

    return args

def main():

    args = get_arguments()

    # detector properties
    detector_material = args.material
    mT = mTarget(detector_material)
    offset = args.offset*METER

    threshold = args.threshold*keV
    exposure = args.exposure*KILOGRAM*YEAR

    # reactor propeties
    thermal_power = args.thermal_power*GIGAWATT
    true_fuel_fractions = np.array([float(f) for f in args.true_fuel_fractions.split(',')])

    # calculate true fission rate per isotope, which is what we input into the rate calculation
    mean_energy_per_fission = sum(true_fuel_fractions*ENERGY_PER_FISSION_I)
    true_total_fission_rate = thermal_power/mean_energy_per_fission
    true_fission_rate_per_isotope = true_total_fission_rate*true_fuel_fractions

    max_fission_rate = 2*true_total_fission_rate # This is going to be used as the upper limit on my fission rate priors

    # specify bin properties (nbins, ERmin, ERmax)
    nbins = args.nbins # this is a stand-in for the energy resolution of the reactor
    ER_min = threshold
    if args.er_max==0:
        ER_max = get_ER_max(FLUX_ENU_MAX, mT)
    else:
        ER_max = args.er_max*keV

    bin_edges = np.logspace(np.log10(ER_min), np.log10(ER_max), nbins+1)

    ### BIN COUNTS ###

    def get_bin_counts(fission_rate_per_isotope, detector_material=detector_material, bin_edges=bin_edges, offset=offset):
        bin_counts = []
        for i in range(nbins):
            bin_counts.append(np.floor(quad(
                dR_dER, bin_edges[i], bin_edges[i+1], args=(detector_material, fission_rate_per_isotope, offset)
                )[0]*exposure))
        return bin_counts
    
    true_bin_counts = get_bin_counts(true_fission_rate_per_isotope)

    ### PRINT TRUE BIN COUNTS ###
    # This lets me test how long it takes to calculate each set of bin_counts,
    # which indicates how long the total job will take to run

    print(true_bin_counts)

    return

if __name__ =='__main__':
    main()
