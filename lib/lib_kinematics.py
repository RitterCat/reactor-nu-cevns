def get_Enu_min(ER, mT):
    return (ER*mT/2)**0.5

def get_ER_max(Enu, mT):
    return 2*Enu**2/(mT + 2*Enu)