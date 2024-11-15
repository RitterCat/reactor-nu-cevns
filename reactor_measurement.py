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

    def get_bin_counts(fuel_fractions, thermal_power, detector_material=detector_material, bin_edges=bin_edges, offset=offset):
        bin_counts = []
        for i in range(nbins):
            bin_counts.append(np.floor(quad(dR_dER, bin_edges[i], bin_edges[i+1], args=(detector_material, fuel_fractions, thermal_power, offset))[0]*KILOGRAM*YEAR))
        return bin_counts

    ### CUBE ###

    def prior(cube, ndim, nparams):
        # The first dimension is the reactor power, P
        # P is selected from a normal distribution with mean thermal_power
        cube[0] = norm.ppf(cube[0], loc=mean_thermal_power, scale = thermal_power_var)
        
        # I then select four fractions that add up to 1 by sorting the next three cube entries and using them to 'partition' the interval [0,1]
        partitions = sorted(cube[1:4])
        p1, p2, p3 = partitions
        f1 = p1
        f2 = p2 - p1
        f3 = p3 - p2
        f4 = 1 - p3
        # The four fractions are then the first four cube entries: f_U235, f_U238, f_Pu239, f_Pu241
        cube[1] = f1
        cube[2] = f2
        cube[3] = f3
        cube[4] = f4
        
        return

    ### LOGLIKE ###

    def loglike(cube, ndim, nparams):

        thermal_power = cube[0]
        fuel_fractions = cube[1:5]

        true_bin_counts = get_bin_counts(true_fuel_fractions, thermal_power)

        bin_counts = get_bin_counts(fuel_fractions, thermal_power)
        
        loglikelihood = 0
        for bc, tbc in zip(bin_counts, true_bin_counts):
            logPoisson = poisson.logpmf(bc, tbc)
            loglikelihood += logPoisson
        
        return loglikelihood

    ndims = 4
    nparams = 5

    pymultinest.run( loglike, prior, n_dims=ndims, n_params = nparams,
                outputfiles_basename="out/ff2323/partitioned2323", verbose=True,
                importance_nested_sampling = False, resume = False, n_live_points = 100,
                sampling_efficiency=0.8, evidence_tolerance=0.5)