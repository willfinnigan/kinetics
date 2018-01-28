from Kinetics import *

""" ------------------- SUBSTRATE AND ENZYME STARTING CONDITIONS  ------------------- """

species_with_pc_error = {"Ester": (4000, 0.05),
                         "Acid": (0, 0),
                         "Aldehyde": (0, 0),
                         "Alcohol": (0, 0),
                         "Methanol" : (0, 0),

                         'afEst2': (10, 0.05),

                         }

# This makes species_default a dictionary without the error.  eg {"Ester" : 1000, "Acid" : 0}
species_defaults = set_species_defaults(species_with_pc_error)

""" ------------------- PARAMETERS  ------------------- """

parameters_with_error = {"afest_km_ester": (1500, 200),
                         "afest_kcat": (6, 1),

                         }

# This makes parameters_default a dictionary without the error.  eg {"kcat_est" : 1, "km_est" : 0.5}
parameter_defaults = set_parameter_defaults(parameters_with_error)

""" ------------------- REACTIONS  ------------------- """

def esterase_r1(y, s_names, params):
    # Substrate values
    ester = y[s_names.index("Ester")]
    esterase = y[s_names.index("afEst2")]

    # Parameters
    km_ester = params["afest_km_ester"]
    kcat_esterase = params["afest_kcat"]

    # Calculate the rate
    rate = one_substrate_mm(kcat_esterase, esterase, ester, km_ester)

    # Return the change in substrate concentration
    # This section += or -= yprime at the correct index by the rate for the substrates and products listed
    substrates = ["Ester"]
    products = ["Acid", "Methanol"]
    yprime = np.zeros(len(y))
    yprime = yprime_minus(yprime, rate, substrates, s_names)
    yprime = yprime_plus(yprime, rate, products, s_names)

    return yprime



""" ------------------- MODEL TIME SETTINGS  ------------------- """
start = 0
end = 240
number_steps = 241
mx_step = 10000

""" ------------------- MODEL  ------------------- """
# Make an instance of the Model class containing the ordered species names and the parameters
model = Model()
model.set_parameters(parameter_defaults)
model.set_species_names_and_starting_values(species_with_pc_error)
model.set_time(start, end, number_steps, mxsteps=mx_step)

# Add the reaction functions
model.append(esterase_r1)


