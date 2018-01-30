import kinetics

""" ------------------- SUBSTRATE AND ENZYME STARTING CONDITIONS  ------------------- """

species_with_pc_error = {"Ester": (2800, 0.05),
                         "Acid": (0, 0),
                         "Aldehyde": (0, 0),
                         "Alcohol": (0, 0),
                         "Methanol": (0, 0),
                         'afEst2': (10, 0.05),
                         }

# From species_with_pc_error make two other dictionaries
# species_defaults to set the default values for the model (ie species that wont change)
# species_bounds to make a dict which contains the upper and lower bounds (calculated from the percentage error)
species_defaults = kinetics.set_species_defaults(species_with_pc_error)
species_bounds = kinetics.get_bounds_from_pc_error(species_with_pc_error)

""" ------------------- PARAMETERS  ------------------- """

parameters_with_std_error = {"afest_km_ester": (1500, 200),
                             "afest_kcat": (6, 1),
                             }

# From parameters_with_std_error make two other dictionaries
# parameter_defaults to set the default values for the model (ie species that wont change)
# parameter_bounds to make a dict which contains the upper and lower bounds (calculated from the std error)
parameter_defaults = kinetics.set_parameter_defaults(parameters_with_std_error)
parameter_bounds = kinetics.get_bounds_from_std_error(parameters_with_std_error)

""" ------------------- MODEL TIME SETTINGS  ------------------- """
start = 0
end = 240
number_steps = 241
max_steps = 10000

""" ------------------- MODEL  ------------------- """
# Make an instance of the Model class.  Set the default parameters, species and time.
model = kinetics.Model()
model.set_parameters(parameter_defaults)
model.set_species(species_defaults)
model.set_time(start, end, number_steps, mxsteps=max_steps)

# Add the reaction functions
from example.reaction_functions import *
model.append(esterase_r1)

""" --- Run the model as an uncertainty analysis and sensitivity analysis --- """
substrates_to_show_ua = ["Ester", "Acid"]
substrate_for_sa_analysis = 'Acid'
timepoint_for_sa_analysis = model.time[-1]

print("---Uncertainty Analysis---")
ua = kinetics.UA(parameter_bounds,
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
sa = kinetics.SA(parameter_bounds,
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
