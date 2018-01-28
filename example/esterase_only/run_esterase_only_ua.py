import esterase_only.model_script
import Kinetics

model = esterase_only.model_script.model
species = esterase_only.model_script.species_with_pc_error
parameters = esterase_only.model_script.parameters_with_error

substrates_to_show = ["Ester", "Acid"]





""" Uncertainty analysis"""
print("---Uncertainty Analysis---")

species_bounds = Kinetics.get_bounds_from_pc_error(species)
parameter_bounds = Kinetics.get_bounds_from_std_error(parameters)

ua = Kinetics.UA(parameter_bounds,
        species_bounds,
        model,
        latin_hypercube_samples=100,
        quartile_range=95)

ua.make_lhc_samples()
ua.run_models()
ua.calculate_quartiles()

ua.print_ua_quartile_output(substrates_to_show)
ua.save_ua_quartile_output(substrates_to_show, filename='')
