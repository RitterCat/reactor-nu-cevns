###########
# IMPORTS #
###########

from lib.lib_rates import *
from scipy.stats import poisson, norm
import pymultinest
import matplotlib.pyplot as plt
import time

##########
# SCRIPT #
##########

if __name__ =='__main__':

    ### SET-UP OF TEST CONDITIONS ###

    # specify detector material and properties
    detector_material = 'Xe'
    mT = mTarget(detector_material)
    offset = 25*METER

    threshold = 0.3*keV
    exposure=10*KILOGRAM*YEAR

    # specify reactor propeties (including true fuel fractions)
    mean_thermal_power = 3*GIGAWATT
    thermal_power_var = 0.1*mean_thermal_power # only need this if thermal power is a prior
    true_fuel_fractions = np.array([0.561, 0.076, 0.307, 0.056])

    # then calculate true fission rate per isotope, which is what we input into the rate calculation
    mean_energy_per_fission = sum(true_fuel_fractions*ENERGY_PER_FISSION_I)
    true_total_fission_rate = mean_thermal_power/mean_energy_per_fission
    true_fission_rate_per_isotope = true_total_fission_rate*true_fuel_fractions

    max_fission_rate = 2*true_total_fission_rate # This is going to be used as the upper limit on my fission rate priors

    # specify bin properties (nbins, ERmin, ERmax)
    nbins = 5 # this is a stand-in for the energy resolution of the reactor
    ER_min = threshold
    ER_max = get_ER_max(FLUX_ENU_MAX, mT)

    bin_edges = np.logspace(np.log10(ER_min), np.log10(ER_max), nbins+1)

    ### BIN COUNTS ###

    def get_bin_counts(fission_rate_per_isotope, offset=offset, detector_material=detector_material, bin_edges=bin_edges):
        bin_counts = []
        for i in range(nbins):
            bin_counts.append(np.floor(quad(
                dR_dER, bin_edges[i], bin_edges[i+1], args=(detector_material, fission_rate_per_isotope, offset)
                )[0]*exposure))
        return bin_counts
    
    true_bin_counts = get_bin_counts(true_fission_rate_per_isotope)

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

        # Now I am storing other parameters in the following cube entries
        cube[4] = sum(fission_rates*ENERGY_PER_FISSION_I) # total power

        cube[5] = f3 + f4 # total fission rate of plutonium

        # the following four parameters are the fuel fractions of each isotope
        total_fission_rate = sum(fission_rates)
        fuel_fractions = fission_rates/total_fission_rate
        cube[6] = fuel_fractions[0]
        cube[7] = fuel_fractions[1]
        cube[8] = fuel_fractions[2]
        cube[9] = fuel_fractions[3]

    def prior_power_four_fractions(cube, ndim, nparams):
        # This prior is for the previous implementation, where we provided the total reactor power and the fuel fractions

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
        # The four fractions are then the next four cube entries: f_U235, f_U238, f_Pu239, f_Pu241
        cube[1] = f1
        cube[2] = f2
        cube[3] = f3
        cube[4] = f4
        
        return

    ### LOGLIKE ###

    def loglike(cube, ndim, nparams):
        fission_rates = cube[:4]

        bin_counts = get_bin_counts(fission_rates)

        loglikelihood = 0
        for bc, tbc in zip(bin_counts, true_bin_counts):
            logPoisson = poisson.logpmf(bc, tbc)
            loglikelihood += logPoisson
        
        return loglikelihood

    def loglike_power_four_fractions(cube, ndim, nparams):

        thermal_power = cube[0]
        fuel_fractions = cube[1:5]

        true_bin_counts = get_bin_counts(true_fuel_fractions, thermal_power)

        bin_counts = get_bin_counts(fuel_fractions, thermal_power)
        
        loglikelihood = 0
        for bc, tbc in zip(bin_counts, true_bin_counts):
            logPoisson = poisson.logpmf(bc, tbc)
            loglikelihood += logPoisson
        
        return loglikelihood
    
    # the point of this loglikelihood is to return nothing, and so test if the priors are being treated properly, or if there is any weird correlation introduced by our sampling method
    def zero_loglike(cube, ndim, nparams):
        return 0

    ndims = 4
    nparams = 10

    pymultinest.run( loglike, prior, n_dims=ndims, n_params = nparams,
                outputfiles_basename="out/RELICS/run1", verbose=True,
                importance_nested_sampling = False, resume = False, n_live_points = 250,
                sampling_efficiency=0.8, evidence_tolerance=0.5)

    #################
    ### TIME TEST ###
    #################

    # The point of this function is to test how long each sampling step will take,
    # which will give an indication of how long the overall run will take 
    # (at least ~5000 times as long)

    def time_test(nbins, threshold, exposure):
        ER_min = threshold
        ER_max = get_ER_max(FLUX_ENU_MAX, mT)

        bin_edges = np.logspace(np.log10(ER_min), np.log10(ER_max), nbins+1)
        true_bin_counts = []
        test_bin_counts = []

        test_fission_rate_per_isotope = true_total_fission_rate*np.array([0,0,0.5,0.5])

        start = time.time()
        for i in range(nbins):
            true_bin_counts.append(np.floor(quad(
                dR_dER, bin_edges[i], bin_edges[i+1], args=(detector_material, true_fission_rate_per_isotope, offset)
                )[0]*exposure))
            
            test_bin_counts.append(np.floor(quad(
                dR_dER, bin_edges[i], bin_edges[i+1], args=(detector_material, test_fission_rate_per_isotope, offset)
                )[0]*exposure))
            
        end = time.time()

        print(true_bin_counts)
        print(test_bin_counts)

        fig, ax = plt.subplots(figsize=(5,7))

        ax.step(bin_edges[:-1], true_bin_counts)
        ax.step(bin_edges[:-1], test_bin_counts)

        ax.set_yscale('log')
        ax.set_xscale('log')

        plt.show()

        return(end - start)

    # print(time_test(20, 0.1*eV, 1000*KILOGRAM*YEAR))