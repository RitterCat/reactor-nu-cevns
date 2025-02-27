{
    "spectrum_settings": {
      "nbins": 100,
      "normalized": "GW",
      "emin": 1.0,
      "emax": 9.0
    },
    
    "data_sources": {
       "u235": "hayes_u235.txt",
       "u238": "hayes_u238.txt",
       "pu239": "hayes_pu239.txt",
       "pu241": "hayes_pu241.txt"
     },
     
     "reactor_data": {
       "type": "manual",
       "power": 0.001,
       "fraction_u235":  0.967,
       "fraction_u238":  0.013,
       "fraction_pu239": 0.02,
       "fraction_pu241": 0
     },
     
     "output_settings": {
       "format": "text",
       "output_name": "spectra/TRIGA_flux.text"
     }
  }
  