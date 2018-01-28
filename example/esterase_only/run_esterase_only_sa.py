import esterase_only.model_script
import Kinetics

model = esterase_only.model_script.model
species = esterase_only.model_script.species_with_pc_error
parameters = esterase_only.model_script.parameters_with_error

substrate_for_analysis = 'Acid'
timepoint_for_analysis = model.time[-1]






""" Sensitivity analysis"""
print("---Sensitivity Analysis---")
species_bounds = Kinetics.get_bounds_from_pc_error(species)
parameter_bounds = Kinetics.get_bounds_from_std_error(parameters)

sa = Kinetics.SA(parameter_bounds,
        species_bounds,
        model,
        timepoint_for_analysis,
        substrate_for_analysis,
        number_samples=500,
        second_order=False,
        conf_level=0.95,
        num_resample=100)

sa.make_saltelli_samples()
sa.run_models()
sa.analyse_sobal_sensitivity()

sa.save_sa_quartile_output(filename='')



