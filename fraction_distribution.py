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
# In this iteration, I generate three random numbers and use them to partition the unit interval [0, 1] into four sections, which I then take as my four fractions

def plot_second_iteration():
    cube = np.random.rand(10000, 3)

    fraction_list = []
    for c_list in cube:
        c1, c2, c3 = sorted(c_list)
        f1 = c1
        f2 = c2 - c1
        f3 = c3 - c2
        f4 = 1 - c3
        fraction_list.append([f1, f2, f3, f4])

    fraction_array = np.array(fraction_list).T

    fig, axes = plt.subplots(2,4, figsize=(10,5))

    for i, (ax, cube_entry) in enumerate(zip(axes[0], cube.T)):
        ax.hist(cube_entry)
        ax.set_xlabel(f'cube[{i}]')

    for i, (ax, fraction_list) in enumerate(zip(axes[1], fraction_array)):
        ax.hist(fraction_list)
        ax.set_xlabel(f'fraction {i}')

    z = np.linspace(0.05, 1, 50)
    minus_logz = -1000*np.log(z) # 1000 = npoints/nbins (10000/10)
    for ax in axes[1]:
        ax.plot(z, minus_logz, 'k--')

    plt.tight_layout()
    plt.show()

    return

### THIRD ITERATION ###
# Here, I am testing the simple idea that Jay suggested: randomly generating four uniform random values and discarding those where the sum is greater than one. 
# ACTUALLY, picking four random doesn't work, as they will probably sum to less than one.
# I'm going to try picking three random, and if their sum is less than one, calculate the fourth

def plot_third_iteration():

    cube = np.random.rand(100000, 3)

    fraction_list = []
    for c_list in cube:
        if sum(c_list) < 1:
            f4 = 1 - sum(c_list)
            fractions = np.append(c_list, f4)
            fraction_list.append(fractions)

    fraction_array = np.array(fraction_list).T

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

### COMPARING TWO AND THREE

def get_partitioned_fractions():
    cube = np.random.rand(50000, 3)

    fraction_list = []
    for c_list in cube:
        c1, c2, c3 = sorted(c_list)
        f1 = c1
        f2 = c2 - c1
        f3 = c3 - c2
        f4 = 1 - c3
        fraction_list.append([f1, f2, f3, f4])

    fraction_array = np.array(fraction_list).T

    return fraction_array

def get_random_fractions():

    cube = np.random.rand(100000, 3)

    fraction_list = []
    for c_list in cube:
        if sum(c_list) < 1:
            f4 = 1 - sum(c_list)
            fractions = np.append(c_list, f4)
            fraction_list.append(fractions)

    fraction_array = np.array(fraction_list).T

    return fraction_array

def compare_two_and_three():
    part_frac = get_partitioned_fractions()
    rand_frac = get_random_fractions()

    fig, axes = plt.subplots(3,4, figsize=(10,8))

    part_counts = []
    for i, (ax, fraction_list) in enumerate(zip(axes[0], part_frac)):
        counts, bins, _bars = ax.hist(fraction_list)
        part_counts.append(counts)

    rand_counts = []
    for i, (ax, fraction_list) in enumerate(zip(axes[1], rand_frac)):
        counts, _bins, _bars = ax.hist(fraction_list)
        rand_counts.append(counts)
        ax.set_xlabel(f'fraction {i}')

    for i, (ax, p_counts, r_counts) in enumerate(zip(axes[2], part_counts, rand_counts)):
        bin_ratios = np.array(p_counts)/np.array(r_counts)
        ax.bar(bins[:-1], bin_ratios)

    # z = np.linspace(0.05, 1, 50)
    # minus_logz = -1*np.log(z)

    # #Getting n_points for each sample allows us to rescale the trial pdf functions correctly
    # nbins = 10
    # scale_pf = len(part_frac[0])/nbins
    # scale_rf = len(rand_frac[0])/nbins
    # scaling_factors = [scale_pf, scale_rf]
    # for ax_row, scaling_factor in zip(axes, scaling_factors):
    #     for ax in ax_row:
    #         ax.plot(z, scaling_factor*minus_logz, 'k--')

    axes[0][0].set_ylabel('Partitioning [0, 1]')
    axes[1][0].set_ylabel('Discarding points \n with sum(f) > 1')

    # axes[0][0].text(0.5, 1500, r'$-\log(f)$', fontsize=14)
    # axes[1][0].text(0.5, 3000, r'$-\log(f)$', fontsize=14)

    plt.tight_layout()
    plt.show()

    return

# They produce the same output! Partitioning [0, 1] is the right idea.


### CODE

if __name__ == "__main__":
    compare_two_and_three()