import kinetics
import numpy as np


class Esterase_r1():
    parameters = {"afest_km_ester": (1500, 200),
                  "afest_kcat": (6, 1)}

    parameter_defaults = kinetics.set_parameter_defaults(parameters)
    parameter_bounds = kinetics.get_bounds_from_std_error(parameters)

    def reaction(self, y, substrate_names, parameters):
        # Substrates
        ester = y[substrate_names.index("Ester")]
        esterase = y[substrate_names.index("afEst2")]

        # Rate equation
        rate = kinetics.one_substrate_mm(kcat=parameters["afest_kcat"],
                                         km_a=parameters["afest_km_ester"],
                                         enz=esterase,
                                         a=ester)

        # Calculate change in substrate concentrations
        substrates = ["Ester"]
        products = ["Acid", "Methanol"]
        y_prime = kinetics.calculate_yprime(y, rate, substrates, products, substrate_names)

        return y_prime