import pymultinest

import numpy as np
import pymultinest
from scipy.stats import norm

# Generate synthetic data (Gaussian with known parameters)
np.random.seed(42)
true_mean = 5
true_std = 2
n_data = 1000
data = np.random.normal(loc=true_mean, scale=true_std, size=n_data)

# Log-likelihood function
def log_likelihood(cube, ndim, nparams):
    mean, std = cube
    # Ensure std is positive (because standard deviation cannot be negative)
    if std <= 0:
        return -np.inf
    # Compute log-likelihood (sum of log of Gaussian PDF values for all data points)
    log_like = np.sum(np.log(norm.pdf(data, loc=mean, scale=std)))
    return log_like

# Prior function (uniform priors for mean and standard deviation)
def prior(cube, ndim, nparams):
    # Uniform prior between -10 and 10 for mean
    cube[0] = cube[0] * 20 - 10  # mean between -10 and 10
    # Uniform prior between 0 and 10 for std (standard deviation must be positive)
    cube[1] = cube[1] * 10  # std between 0 and 10

# Set up the MultiNest sampler
n_dims = 2  # mean and std are the two parameters we're estimating
n_params = 2
output_dir = 'test_multi_output/test'

# Run the sampler
pymultinest.run(log_likelihood, prior, n_dims, outputfiles_basename=output_dir, verbose=True,
                importance_nested_sampling = False, resume = False, n_live_points=100)