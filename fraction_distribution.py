import fractions
import numpy as np
import matplotlib.pyplot as plt

# This code is to test what my prior distributions are for the fuel fractions in the reactor
# I obtain these priors from transforming the cube, which is an n-dimensional hypercube, i.e. each dimension is an independent uniform random variable between 0 and 1

# Let us call the cube dimension d as Cd

### FIRST ITERATION ###
# Initially I just took the first dimension C1 to be the first fraction (f1 = U235), then took f2 = (1 - f1)*C2, so that f2 + f1 < 1 (i.e. that the total sum of fractions is not greater than 1)
# However, it is likely that this means that the distribution of f2 is warped. This code will check that

def plot_first_iteration():
    cube = np.random.rand(10000, 3)

    fractions = []
    for c1, c2, c3 in cube:
        f1 = c1
        f2 = c2*(1-f1)
        f3 = c3*(1- f1 - f2)
        f4 = 1 - f1 - f2 - f3
        fractions.append([f1, f2, f3, f4])

    fraction_array = np.array(fractions).T

    fig, axes = plt.subplots(2,4, figsize=(10,5))

    for i, (ax, cube_entry) in enumerate(zip(axes[0], cube.T)):
        ax.hist(cube_entry)
        ax.set_xlabel(f'cube[{i}]')

    for i, (ax, fraction_list) in enumerate(zip(axes[1], fraction_array)):
        ax.hist(fraction_list)
        ax.set_xlabel(f'fraction {i}')

    z = np.linspace(0.05, 1, 50)
    minus_logz = -1000*np.log(z) # 1000 = npoints/nbins (10000/10)
    logz_sq = 1000/2*np.log(z)**2
    axes[1][1].plot(z, minus_logz, 'k--')
    axes[1][2].plot(z, logz_sq, 'k--')

    axes[1][1].text(0.4, 1500, r'$-\log(f_2)$', fontsize=14)
    axes[1][2].text(0.3, 2500, r'$\log^2(f_3)/2$', fontsize=14)

    plt.tight_layout()
    plt.show()

    return

### SECOND ITERATION ###

def plot_second_iteration():
    cube = np.random.rand(10000, 3)

    fractions = []
    for c_list in cube:
        c1, c2, c3 = sorted(c_list)
        f1 = c1
        f2 = c2 - c1
        f3 = c3 - c2
        f4 = 1 - c3
        fractions.append([f1, f2, f3, f4])

    fraction_array = np.array(fractions).T

    fig, axes = plt.subplots(2,4, figsize=(10,5))

    for i, (ax, cube_entry) in enumerate(zip(axes[0], cube.T)):
        ax.hist(cube_entry)
        ax.set_xlabel(f'cube[{i}]')

    for i, (ax, fraction_list) in enumerate(zip(axes[1], fraction_array)):
        ax.hist(fraction_list)
        ax.set_xlabel(f'fraction {i}')

    plt.tight_layout()
    plt.show()


    return


### CODE

if __name__ == "__main__":
    plot_second_iteration()