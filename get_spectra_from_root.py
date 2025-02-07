# This script opens the ROOT file provided by the supplementary material of arXiv:2034.14992,
# which contains the neutrino and beta fission spectrum data as calcuated by the software
# BESTIOLE (Beta Energy Spectrum Tool for an Improved Optimal List of Elements)

# This script then reads out the spectra for the specified fissile isotopes and saves them to text files

import ROOT
from scipy.interpolate import CubicSpline
import numpy as np

fissileIsotopeSpectraFile = ROOT.TFile('./fluxData/BESTIOLE_2023.root')

isotopes_ROOT_format = ["239U", "239Np"] # These are the two relevant isotopes when breeding 239Pu

Enu_min, Enu_max, nbins = 0, 12.5, 500 # spectra go from 0 MeV to 12.5 MeV
bin_width = (Enu_max - Enu_min)/nbins
nu_energies_BESTIOLE = np.arange(bin_width/2, 12.5, bin_width) # Using the centre of each bin to locate the data point (so have data points between 0.0125 MeV and 12.4875 MeV)

# Get the fluxes from the ROOT file
fissionSpectra = fissileIsotopeSpectraFile.Get("activation")
fission_fluxes_BESTIOLE = []
for isotope in isotopes_ROOT_format:
    fission_fluxes_BESTIOLE.append([flux for flux in fissionSpectra.Get(f"{isotope}_nspec")][1:-1])

# Write the spectra out to text files
isotopes_txtfile_format = ['u239', 'np239']
for isotope, fission_flux in zip(isotopes_txtfile_format, fission_fluxes_BESTIOLE):
    with open(f'./fluxData/bestiole_{isotope}.txt', 'w') as outFile:
        outFile.write('#Data from https://arxiv.org/abs/2304.14992 ancillary materials\n')
        for i, (energy, flux) in enumerate(zip(nu_energies_BESTIOLE, fission_flux)):
            line = f'{energy}, {flux}'
            if i != nbins - 1:
                line += '\n'
            outFile.write(line)