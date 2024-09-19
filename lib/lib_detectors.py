from lib.lib_constants import *

# isotope properties for different detector materials
# isotopes of xenon
xe124 = {"Z": 54, "A": 124, "abundance": 0.00095, "mass": 123.906*AMU}
xe126 = {"Z": 54, "A": 126, "abundance": 0.00089, "mass": 125.904*AMU}
xe128 = {"Z": 54, "A": 128, "abundance": 0.0191, "mass": 127.904*AMU}
xe129 = {"Z": 54, "A": 129, "abundance": 0.26401, "mass": 128.905*AMU}
xe130 = {"Z": 54, "A": 130, "abundance": 0.04071, "mass": 129.904*AMU}
xe131 = {"Z": 54, "A": 131, "abundance": 0.21232, "mass": 130.905*AMU}
xe132 = {"Z": 54, "A": 132, "abundance": 0.26909, "mass": 131.904*AMU}
xe134 = {"Z": 54, "A": 134, "abundance": 0.10436, "mass": 133.905*AMU}
xe136 = {"Z": 54, "A": 136, "abundance": 0.08857, "mass": 135.907*AMU}
isotopesXe = [xe124, xe126, xe128, xe129, xe130, xe131, xe132, xe134, xe136]
mXe = sum([iso["abundance"]*iso["mass"] for iso in isotopesXe])

# isotopes of silicon
si28 = {"Z": 14, "A": 28, "abundance": 0.92223, "mass": 27.977*AMU}
si29 = {"Z": 14, "A": 29, "abundance": 0.04685, "mass": 28.977*AMU}
si30 = {"Z": 14, "A": 30, "abundance": 0.03092, "mass": 29.974*AMU}
isotopesSi = [si28, si29, si30]
mSi = sum([iso["abundance"]*iso["mass"] for iso in isotopesSi])

#isotopes of germanium
ge70 = {"Z": 32, "A": 70, "abundance": 0.2052, "mass": 69.924*AMU}
ge72 = {"Z": 32, "A": 72, "abundance": 0.2745, "mass": 71.922*AMU}
ge73 = {"Z": 32, "A": 73, "abundance": 0.0776, "mass": 72.923*AMU}
ge74 = {"Z": 32, "A": 74, "abundance": 0.3652, "mass": 73.921*AMU}
ge76 = {"Z": 32, "A": 76, "abundance": 0.0775, "mass": 75.921*AMU}
isotopesGe = [ge70, ge72, ge73, ge74, ge76]
mGe = sum([iso["abundance"]*iso["mass"] for iso in isotopesGe])

ISOTOPES = {"Xe": isotopesXe, "Si": isotopesSi, "Ge": isotopesGe}

# function to get the target mass (i.e. average mass of isotopes) for a given element
def mTarget(element):
    return sum([iso["abundance"]*iso["mass"] for iso in ISOTOPES[element]])