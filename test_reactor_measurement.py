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
    parser.add('--background_rate', type=float, help="The rate of background events (in events/day/kg/keV)")
    parser.add('--nbins', type=int, 
               help="""Number of bins between the threshold and the 
               max recoil energy in the neutrino detector""")
    
    parser.add('-r', '--reactor-config', is_config_file=True,
               help="Path to the configuration file for the properties of the nuclear reactor")
    parser.add('--thermal_power', type=float, help="Reactor power (in gigawatt)")
    parser.add('--true_fuel_fractions', type=str, help="A comma-delimited string of the true fuel fractions")
    parser.add('--plutonium_rate', type=float, help="The rate of production of plutonium in the breeding blanket (in kg/year)")

    parser.add('-p', '--pymultinest-config', is_config_file=True,
               help="Path to the configuration file that contains options for the pymultinest run")
    parser.add('--sampling_efficiency', type=float, default=0.8)
    parser.add('--evidence_tolerance', type=float, default=0.5)
    parser.add('--n_live_points', type=int, default=250)

    parser.add('-o', '--outputfiles_basename', type=str)

    parser.add('--resume', action='store_true', default=False, help="Make the pymultinest run resume from a previous state")

    parser.add('--power_prior', action='store_true', default=False, help="Use the knowledge of the reactor power as a prior")

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

    background_rate = args.background_rate/KILOGRAMDAY/keV

    # reactor propeties
    thermal_power = args.thermal_power*GIGAWATT
    true_fuel_fractions = np.array([float(f) for f in args.true_fuel_fractions.split(',')])

    plutonium_rate = args.plutonium_rate*KILOGRAM/YEAR
    true_breeding_rate = plutonium_rate/M_PU239

    max_breeding_rate = 2*true_breeding_rate

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

    def get_bin_counts(fission_rate_per_isotope, breeding_rate=0, background_rate=background_rate, detector_material=detector_material, bin_edges=bin_edges, offset=offset):
        bin_counts = []
        for i in range(nbins):
            cevns_counts = quad(
                dR_dER, bin_edges[i], bin_edges[i+1], args=(detector_material, offset, fission_rate_per_isotope, breeding_rate)
                )[0]*exposure
            
            bin_width = bin_edges[i+1] - bin_edges[i]            
            background_counts = background_rate*exposure*bin_width

            print(cevns_counts, background_counts)

            bin_counts.append(np.floor(cevns_counts + background_counts))
        return bin_counts
    
    true_bin_counts = get_bin_counts(true_fission_rate_per_isotope, true_breeding_rate)

    print(true_bin_counts)

    ### PRIOR ###

    def prior(cube, ndim, nparams):
        # The first four cube entries here are just the fission rates for each isotope
        f1 = cube[0]*max_fission_rate
        f2 = cube[1]*max_fission_rate
        f3 = cube[2]*max_fission_rate
        f4 = cube[3]*max_fission_rate
        fission_rates = np.array([f1, f2, f3, f4])
        cube[0] = f1
        cube[1] = f2
        cube[2] = f3
        cube[3] = f4

        # The fifth cube entry is the breeding rate
        breeding_rate = cube[4]*max_breeding_rate
        cube[4] = breeding_rate

        # Now I am storing other parameters in the following cube entries
        cube[5] = sum(fission_rates*ENERGY_PER_FISSION_I) # total power

        cube[6] = f3 + f4 # total fission rate of plutonium

        # the following four parameters are the fuel fractions of each isotope
        total_fission_rate = sum(fission_rates)
        fuel_fractions = fission_rates/total_fission_rate
        cube[7] = fuel_fractions[0]
        cube[8] = fuel_fractions[1]
        cube[9] = fuel_fractions[2]
        cube[10] = fuel_fractions[3]

    ### LOGLIKE ###

    def loglike(cube, ndim, nparams):
        fission_rates = cube[:4]

        breeding_rate = cube[4]

        bin_counts = get_bin_counts(fission_rates, breeding_rate)

        loglikelihood = 0
        for bc, tbc in zip(bin_counts, true_bin_counts):
            logPoisson = poisson.logpmf(bc, tbc)
            loglikelihood += logPoisson

        if args.power_prior:
            logNorm = norm.logpdf(cube[4], loc=thermal_power, scale=0.05*thermal_power)
            loglikelihood += logNorm
        
        return loglikelihood

    ndims = 5
    nparams = 11
    
    outputfiles_basename = args.outputfiles_basename

    sampling_efficiency = args.sampling_efficiency
    evidence_tolerance = args.evidence_tolerance
    n_live_points = args.n_live_points
            

    # pymultinest.run( loglike, prior, n_dims=ndims, n_params = nparams,
    #         outputfiles_basename=f"/data/gpfs/projects/punim0011/ritter/reactor_nu/out/{outputfiles_basename}", 
    #             verbose=True, importance_nested_sampling = False, resume = args.resume, 
    #             n_live_points=n_live_points, sampling_efficiency=sampling_efficiency, evidence_tolerance=evidence_tolerance)

if __name__ =='__main__':
    main()
