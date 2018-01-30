import Kinetics
from rate_functions import *

""" ------------------- SUBSTRATE AND ENZYME STARTING CONDITIONS  ------------------- """

species_with_pc_error = {"Ester": (2800, 0.05),
                         "Acid": (0, 0),
                         "Aldehyde": (0, 0),
                         "Alcohol": (0, 0),
                         "Methanol" : (0, 0),

                         'afEst2': (10, 0.05),

                         }

# This makes species_default a dictionary without the error.  eg {"Ester" : 1000, "Acid" : 0}
species_defaults = Kinetics.set_species_defaults(species_with_pc_error)
species_bounds = Kinetics.get_bounds_from_pc_error(species_with_pc_error)

""" ------------------- PARAMETERS  ------------------- """

parameters_with_error = {"afest_km_ester": (1500, 200),
                         "afest_kcat": (6, 1),

                         }

# This makes parameters_default a dictionary without the error.  eg {"kcat_est" : 1, "km_est" : 0.5}
parameter_defaults = Kinetics.set_parameter_defaults(parameters_with_error)
parameter_bounds = Kinetics.get_bounds_from_std_error(parameters_with_error)

""" ------------------- MODEL TIME SETTINGS  ------------------- """
start = 0
end = 240
number_steps = 241
max_steps = 10000

""" ------------------- MODEL  ------------------- """
# Make an instance of the Model class containing the ordered species names and the parameters
model = Kinetics.Model()
model.set_parameters(parameter_defaults)
model.set_species_names_and_starting_values(species_with_pc_error)
model.set_time(start, end, number_steps, mxsteps=max_steps)

# Add the reaction functions
model.append(esterase_r1)



""" --- Run the model as an uncertainty analysis and sensitivity analysis --- """
substrates_to_show_ua = ["Ester", "Acid"]
substrate_for_sa_analysis = 'Acid'
timepoint_for_sa_analysis = model.time[-1]

print("---Uncertainty Analysis---")
ua = Kinetics.UA(parameter_bounds,
                 species_bounds,
                 model,
                 num_samples=1000,
                 quartile_range=95)
ua.make_lhc_samples()
ua.run_models()
ua.calculate_quartiles()
ua.print_ua_quartile_output(substrates_to_show_ua)
ua.save_ua_quartile_output(substrates_to_show_ua, filename='uncertainty_analysis.txt')


print("---Sensitivity Analysis---")
sa = Kinetics.SA(parameter_bounds,
        species_bounds,
        model,
        timepoint_for_sa_analysis,
        substrate_for_sa_analysis,
        number_samples=500,
        second_order=False,
        conf_level=0.95,
        num_resample=100)

sa.make_saltelli_samples()
sa.run_models()
sa.analyse_sobal_sensitivity()
sa.save_sa_quartile_output(filename='sensitivity_analysis.txt')




