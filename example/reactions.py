import kinetics
import numpy as np


esterase_params = {"afest_km_ester": (1500, 200),
                   "afest_kcat": (6, 1)}

def esterase_r1(y, s_names, params):
    # Substrate values
    ester = y[s_names.index("Ester")]
    esterase = y[s_names.index("afEst2")]

    # Parameters
    km_ester = params["afest_km_ester"]
    kcat_esterase = params["afest_kcat"]

    rate = kinetics.one_substrate_mm(kcat_esterase, km_ester, esterase, ester)

    substrates = ["Ester"]
    products = ["Acid", "Methanol"]
    y_prime = kinetics.calculate_yprime(y, rate, substrates, products, s_names)

    return y_prime
