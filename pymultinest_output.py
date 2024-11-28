import pymultinest
import matplotlib.pyplot as plt

def corner_plot(n_pars, parameters, outputfiles_basename, par_range=None):

    if not par_range:
        par_min = 0
        par_max = n_pars - 1
        n_plots = n_pars
    else:
        par_min, par_max = par_range
        n_plots = par_max - par_min + 1

    a = pymultinest.Analyzer(n_pars, outputfiles_basename=outputfiles_basename)

    p = pymultinest.PlotMarginalModes(a)
    plt.figure(figsize=(5*n_plots, 5*n_plots))

    for i in range(n_plots):
        plt.subplot(n_plots, n_plots, n_plots * i + i + 1)
        p.plot_marginal(i + par_min, with_ellipses = True, with_points = False, grid_points=50)
        plt.ylabel("Probability")
        plt.xlabel(parameters[i])

        for j in range(i):
            plt.subplot(n_plots, n_plots, n_plots * i + j + 1)
            p.plot_marginal(j + par_min, i + par_min, with_ellipses = False, with_points = False, grid_points=30)
            plt.xlabel(parameters[j])
            plt.ylabel(parameters[i])

    plt.subplots_adjust(wspace=.5)

    plt.show()

# HERE IS WHERE I STORE VARIOUS RUNS THAT I HAVE DONE
# The variables specified are n_pars, a list of parameter names, and the basename of the output files

THREE_FRACTIONS = ['f_U235', 'f_U238', 'f_Pu239']
POWER_AND_FOUR_FRACTIONS = ['P_th', 'f_U235', 'f_U238', 'f_Pu239', 'f_Pu241']

ff2323 = (3, THREE_FRACTIONS, "out/ff2323/2323")

partitioned2323 = (5, POWER_AND_FOUR_FRACTIONS, "out/ff2323/partitioned2323")

partitioned2323_4bins = (5, POWER_AND_FOUR_FRACTIONS, "out/ff2323/partitioned2323_4bins", [1, 4])

ff2323_run1 = (5, POWER_AND_FOUR_FRACTIONS, "out/ff2323/run1", [1, 4])

prior_test = (5, POWER_AND_FOUR_FRACTIONS, "out/prior_test/test")

if __name__ == "__main__":
    corner_plot(*ff2323_run1)