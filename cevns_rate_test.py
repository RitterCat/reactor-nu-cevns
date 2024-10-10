###########
# IMPORTS #
###########

from lib.lib_rates import *
from scipy.stats import poisson, norm
import pymultinest

##########
# SCRIPT #
##########

if __name__ =='__main__':

    ### FUNCTION TO GET BINNED EVENTS ###

    # specify detector material
    detector_material = 'Xe'
    mT = mTarget(detector_material)
    offset = 10*METER

    # specify reactor propeties (including true fuel fractions)
    mean_thermal_power = GIGAWATT
    thermal_power_var = 0.1*mean_thermal_power
    true_fuel_fractions = [0.2, 0.3, 0.2, 0.3]

    # specify bin properties (nbins, ERmin, ERmax)
    nbins = 5
    ER_min = 0.1*keV # this is the threshold value
    ER_max = get_ER_max(FLUX_ENU_MAX, mT)

    bin_edges = np.logspace(np.log10(ER_min), np.log10(ER_max), nbins+1)

    def get_bin_counts(fuel_fractions, thermal_power):
        bin_counts = []
        for i in range(nbins):
            bin_counts.append(np.floor(quad(dR_dER, bin_edges[i], bin_edges[i+1], args=(detector_material, fuel_fractions, thermal_power, offset))[0]*KILOGRAM*YEAR))
        return bin_counts

    ### CUBE ###

    def prior(cube, ndim, nparams):
        # f_U235: can be anything between zero and 1
        cube[0] = cube[0]
        # f_U238: cannot be greater than 1 - f_U235
        cube[1] = cube[1] * (1 - cube[0])
        # f_Pu239: cannot be greater than 1 - f_U235 - f_U238
        cube[2] = cube[2] * (1 - cube[0] - cube[1])
        # f_Pu241: is now specified by the choices of the other three fractions
        cube[3] = 1 - cube[0] - cube[1] - cube[2]
        # P is selected from a normal distribution with mean thermal_power
        cube[4] = norm.ppf(cube[4], loc=mean_thermal_power, scale = thermal_power_var)

        return

    ### LOGLIKE ###

    def loglike(cube, ndim, nparams):

        fuel_fractions = cube[:4]
        thermal_power = cube[4]


        true_bin_counts = get_bin_counts(true_fuel_fractions, thermal_power)

        bin_counts = get_bin_counts(fuel_fractions, thermal_power)
        
        loglikelihood = 0
        for bc, tbc in zip(bin_counts, true_bin_counts):
            logPoisson = poisson.logpmf(bc, tbc)
            loglikelihood += logPoisson
        
        return loglikelihood

    ndims = 3 
    nparams = 4

    # pymultinest.run( loglike, prior, n_dims=ndims, n_params = nparams,
    #             outputfiles_basename="out/ff2323/2323", verbose=False,
    #             importance_nested_sampling = False, resume = False, n_live_points = 100,
    #             sampling_efficiency=0.8, evidence_tolerance=0.5)