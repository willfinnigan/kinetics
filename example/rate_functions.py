import Kinetics
import numpy as np


def esterase_r1(y, s_names, params):
    # Substrate values
    ester = y[s_names.index("Ester")]
    esterase = y[s_names.index("afEst2")]

    # Parameters
    km_ester = params["afest_km_ester"]
    kcat_esterase = params["afest_kcat"]

    # Calculate the rate
    rate = Kinetics.one_substrate_mm(kcat_esterase, esterase, ester, km_ester)

    # Return the change in substrate concentration
    # This section += or -= yprime at the correct index by the rate for the substrates and products listed
    substrates = ["Ester"]
    products = ["Acid", "Methanol"]
    yprime = np.zeros(len(y))
    yprime = Kinetics.yprime_minus(yprime, rate, substrates, s_names)
    yprime = Kinetics.yprime_plus(yprime, rate, products, s_names)

    return yprime
