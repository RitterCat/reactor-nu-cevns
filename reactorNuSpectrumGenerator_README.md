# reactorNuSpectrumGenerator
Written originally by Samuel Hedges. Edited by Alex Ritter.
Allows user to specify different sources of reactor neutrino isotope fluxes calculations, weights, and energy range, and generates an output neutrino spectrum. Still a work in progress. Things to change are:
* Improve or give the user extrapolation options
* Include errors authors use in calculations. More difficult when data extracted from curves.

# Usage:
python reactorNuSpectrumGeneratory.py < config file >

# Explanation of config file
## spectrum_settings
* nbins: number of bins in output spectrum
* normalized: either "GW" for neutrinos/sec using supplied power, "pdf" for normalization to integral of 1, or anything else for neutrinos/fission
* emin: minimum energy of output spectrum
* emax: maximum energy of output spectrum

## data_sources
* For each isotope, filename to use from fluxData folder. Format of data is lines of "energy [in MeV], antineutrinos/fission". Comments start with a "#". If data does not cover requested (emin,emax) range, exponential extrapolation is done using five closest data points.
  
## reactor_data
* type: "manual"
* power: power in GW
* fraction_u235: fraction of u235 in core, does not have to be normalized
* fraction_u238: fraction of u238 in core, does not have to be normalized
* fraction_pu239: fraction of pu239 in core, does not have to be normalized
* fraction_pu241: fraction of pu241 in core, does not have to be normalized

## output_settings:
* format: either "text" or "root"
* output_name: name of the output file, with appropriate extension

# Data references
* Hayes calculations: https://arxiv.org/abs/1605.02047
* Huber calculations: https://arxiv.org/pdf/1106.0687.pdf
* Klapdor Pu239 calculations: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.48.127
* Klapdor U235 calculations: https://www.sciencedirect.com/science/article/pii/0370269382908978
* Mueller calculations: https://arxiv.org/pdf/1101.2663.pdf
* Avignone calculation: https://journals.aps.org/prc/pdf/10.1103/PhysRevC.22.594