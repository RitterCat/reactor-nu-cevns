import pymultinest
import matplotlib.pyplot as plt

def plot2233():

    n_pars=3
    parameters= ['f_U235', 'f_U238', 'f_Pu239']

    a = pymultinest.Analyzer(n_pars, outputfiles_basename="out/ff2323/2323")

    p = pymultinest.PlotMarginalModes(a)
    plt.figure(figsize=(5*n_pars, 5*n_pars))

    for i in range(n_pars):
        plt.subplot(n_pars, n_pars, n_pars * i + i + 1)
        p.plot_marginal(i, with_ellipses = True, with_points = False, grid_points=50)
        plt.ylabel("Probability")
        plt.xlabel(parameters[i])

        for j in range(i):
            plt.subplot(n_pars, n_pars, n_pars * i + j + 1)
            p.plot_marginal(j, i, with_ellipses = False, with_points = False, grid_points=30)
            plt.xlabel(parameters[j])
            plt.ylabel(parameters[i])

    plt.subplots_adjust(wspace=.3)

    plt.show()

    return

def plot_partitioned2233():

    n_pars = 5
    parameters = ['P_th', 'f_U235', 'f_U238', 'f_Pu239', 'f_Pu241']

    a = pymultinest.Analyzer(n_pars, outputfiles_basename="out/ff2323/partitioned2323")

    p = pymultinest.PlotMarginalModes(a)
    plt.figure(figsize=(5*n_pars, 5*n_pars))

    for i in range(n_pars):
        plt.subplot(n_pars, n_pars, n_pars * i + i + 1)
        p.plot_marginal(i, with_ellipses = True, with_points = False, grid_points=50)
        plt.ylabel("Probability")
        plt.xlabel(parameters[i])

        for j in range(i):
            plt.subplot(n_pars, n_pars, n_pars * i + j + 1)
            p.plot_marginal(j, i, with_ellipses = False, with_points = False, grid_points=30)
            plt.xlabel(parameters[j])
            plt.ylabel(parameters[i])

    plt.subplots_adjust(wspace=.3)

    plt.show()

    return

if __name__ == "__main__":
    plot_partitioned2233()