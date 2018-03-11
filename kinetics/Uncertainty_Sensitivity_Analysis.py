from SALib.sample import latin, saltelli
from SALib.analyze import sobol
from kinetics.Model import *
import numpy as np
import pandas as pd
from tqdm import tqdm
import time



""" Setup the bounds for uncertainty or sensitivity analsyis"""
def setup_bounds_lists(dict_with_bounds):
    names_list = []
    bounds_list = []

    for name in dict_with_bounds:
        if dict_with_bounds[name] != 0:
            names_list.append(name)
            bounds_list.append(dict_with_bounds[name])

    return names_list, bounds_list


def get_bounds_from_std_error(dict_with_std_error):
    dict_with_bounds = {}

    for name in dict_with_std_error:
        if (dict_with_std_error[name][0] != 0) and (dict_with_std_error[name][1] != 0):
            upper_bound = dict_with_std_error[name][0] + dict_with_std_error[name][1]
            lower_bound = dict_with_std_error[name][0] - dict_with_std_error[name][1]

            if lower_bound < 0:
                lower_bound = 0

            dict_with_bounds[name] = (lower_bound, upper_bound)

        else:
            dict_with_bounds[name] = 0

    return dict_with_bounds


def get_bounds_from_pc_error(dict_with_pc_error):
    dict_with_bounds = {}

    for name in dict_with_pc_error:
        if (dict_with_pc_error[name][0] != 0) and (dict_with_pc_error[name][1] != 0):
            pc_error = dict_with_pc_error[name][1] * dict_with_pc_error[name][0]
            upper_bound = dict_with_pc_error[name][0] + pc_error
            lower_bound = dict_with_pc_error[name][0] - pc_error

            if lower_bound < 0:
                lower_bound = 0

            dict_with_bounds[name] = (lower_bound, upper_bound)

        else:
            dict_with_bounds[name] = 0

    return dict_with_bounds


""" Set up the problem dict which will be used for sampling"""
def setup_problem(parameter_names, parameter_bounds,
                  species_names, species_bounds):
    names_list = []
    names_list.extend(parameter_names)
    names_list.extend(species_names)

    bounds = []
    bounds.extend(parameter_bounds)
    bounds.extend(species_bounds)

    problem = {'num_vars': len(names_list),
               'names': names_list,
               'bounds': bounds}

    return problem


""" Run the models for the uncertainty analysis or sensitivity analysis"""
def parse_samples_to_run(samples, parameter_names, species_names):
    parsed_samples = []

    def unpack_sampled_params_and_species(sample_list, parameter_names_list, species_names_list):
        parameters = {}
        count = 0
        for i in range(len(parameter_names_list)):
            name = parameter_names_list[i]
            parameters[name] = sample_list[i]
            count += 1

        species = {}
        for i in range(count, len(species_names_list) + count):
            name = species_names_list[i - count]
            species[name] = sample_list[i]

        return parameters, species

    for sampled in samples:
        parameters, species = unpack_sampled_params_and_species(sampled, parameter_names, species_names)
        parsed_samples.append([parameters, species])

    # returns a list of tuples containing [(parameter_dict, species_dict), ] for each sample
    return parsed_samples

def run_all_models(parsed_samples, model):
    ua_sa_output = []

    for ua_parameters, ua_species in tqdm(parsed_samples):
        model.update_species(ua_species)
        model.update_parameters(ua_parameters)

        y = model.run_model()

        ua_sa_output.append(y)


    model.reset_model()

    # ua_output will be a list like [y1, y2, y3, ect...]
    return ua_sa_output



""" Process multiple model outputs for uncertainty analysis"""
def return_ys_for_a_single_substrate(time, list_y, species_names, name):
    """
    Collect all the runs for a single substrate (name)

    collected_output = [ [t0, r1, r2, r3],
                         [t1, r1, r2, r3]....
    """

    collected_output = []

    for i in range(len(time)):
        timepoint = []
        timepoint.append(time[i])

        for y in list_y:
            timepoint.append(y[i][species_names.index(name)])

        collected_output.append(timepoint)

    return collected_output

def return_quartiles(ys_for_single_substrate, name, quartile=95):
    quartiles = {"Time": [], "High": [], "Low": [], "Mean": []}

    for i in range(len(ys_for_single_substrate)):
        # output_at_t will be a array.  i=0 is time, after than the substrate values at that time.
        output_at_t = ys_for_single_substrate[i]

        quartiles["Time"].append(output_at_t[0])
        quartiles["High"].append(np.percentile(output_at_t[1:], quartile))
        quartiles["Low"].append(np.percentile(output_at_t[1:], 100 - quartile))
        quartiles["Mean"].append(np.mean(output_at_t[1:]))

    return quartiles




"""Classes to do UA or SA"""
class UA():
    def __init__(self,
                 model,
                 num_samples=1000,
                 quartile_range=95):

        self.parameters_with_bounds = model.parameter_bounds
        self.species_with_bounds = model.species_bounds

        self.parameter_names, self.parameter_bounds = setup_bounds_lists(self.parameters_with_bounds)
        self.species_names, self.species_bounds = setup_bounds_lists(self.species_with_bounds)

        self.model = model
        self.num_samples = num_samples
        self.quartile_range = quartile_range

        self.samples = []
        self.parsed_samples = []
        self.problem = 0

        self.ua_output = []
        self.quartile_output = {}

        self.substrate_dataframes = {}

    def analyse_number_of_samples_vs_parameters(self):
        # Code for text output
        num_params = len(self.parameter_names)
        num_specs = len(self.species_names)
        total = num_params + num_specs
        total_pw = total * total
        total_pw3 = total * total * total
        print("")
        print(num_params, 'parameters and', num_specs, "species in uncertainty analysis")
        print(total, "variables in total")
        print(total, "^2 =", total_pw)
        print(total, "^3 =", total_pw3)
        print(str(self.num_samples), "samples made by lhc")
        print("")

    def make_lhc_samples(self):
        self.problem = setup_problem(self.parameter_names, self.parameter_bounds,
                                     self.species_names, self.species_bounds)

        self.samples = latin.sample(self.problem, self.num_samples)

        self.parsed_samples = parse_samples_to_run(self.samples, self.parameter_names, self.species_names)

        return self.parsed_samples

    def run_models(self):
        print("running all models")
        time.sleep(0.5)
        self.ua_output = run_all_models(self.parsed_samples, self.model)
        print("samples run")
        self.model.reset_model()

        return self.ua_output

    def quartiles_to_dataframe(self):
        self.substrate_dataframes = {}
        for substrate in self.model.species_names:
            df = pd.DataFrame(self.quartile_output[substrate], columns=["Time", "High", "Low", "Mean"])
            self.substrate_dataframes[substrate] = df

        return self.substrate_dataframes

    def calculate_quartiles(self, quartile_range=95):

        for name in self.model.species_names:
            collected_runs = return_ys_for_a_single_substrate(self.model.time, self.ua_output, self.model.species_names, name)
            quartiles = return_quartiles(collected_runs, name, quartile=quartile_range)

            self.quartile_output[name] = quartiles

        print("quartiles calculated")

        self.quartiles_to_dataframe()

        print('saved as dataframes in a dict self.substrate_dataframes:  {"substrate" : dateframe}')

        return self.substrate_dataframes

    def return_ua_info(self):
        info = "--- Species Defaults in Reaction --- \n"
        for name in self.model.species_names:
            if self.model.species_defaults[name] > 0:
                info += (str(name) + " : " + str(self.model.species_defaults[name]) + "\n")

        info += "\n ---Species Concentration bounds--- \n"
        for i in range(len(self.species_names)):
            info += (str(self.species_names[i]) + " : " + str(self.species_bounds[i]) + "\n")

        info += ("\n")

        info += ("--- Parameters Defaults--- \n")
        for name in self.model.parameters:
            info += (str(name) + " : " + str(self.model.parameters[name]) + "\n")

        info += ("\n")

        info += ("--- Parameters Bounds--- \n")
        for i in range(len(self.parameter_bounds)):
            info += (str(self.parameter_names[i]) + " : " + str(self.parameter_bounds[i]) + "\n")

        info += ("\n")

        info += ("Quartile range = " + str(self.quartile_range) + "%" + "\n")
        info += ("Number of samples = " + str(self.num_samples) + "\n")
        info += ("\n")
        num_params = len(self.parameter_names)
        num_specs = len(self.species_names)
        total = num_params + num_specs
        info += (str(num_params) + ' parameters and ' + str(num_specs) + " species in uncertainty analysis" + "\n")
        info += (str(total) + " variables in total" + "\n")
        info += (str(total) + "^2 = " + str(total * total) + "\n")
        info += (str(total) + "^3 = " + str(total * total * total) + "\n")
        info += (str(self.num_samples) + " samples made by lhc" + "\n")

        return info

class SA():
    def __init__(self,
                 model,
                 number_samples=500,
                 second_order=False,
                 conf_level=0.95,
                 num_resample=100):



        self.number_samples = number_samples
        self.second_order = second_order
        self.conf_level = conf_level
        self.num_resample = num_resample

        self.parameters_with_bounds = model.parameter_bounds
        self.species_with_bounds = model.species_bounds

        self.parameter_names, self.parameter_bounds = setup_bounds_lists(self.parameters_with_bounds)
        self.species_names, self.species_bounds = setup_bounds_lists(self.species_with_bounds)

        self.model = model

        self.samples = []
        self.parsed_samples = []
        self.problem = {}

        self.output = []
        self.output_for_analysis = []
        self.analysis = []

        self.dataframe_output = 0

        self.analysis_info = ''

    def make_saltelli_samples(self):
        self.problem = setup_problem(self.parameter_names, self.parameter_bounds,
                                     self.species_names, self.species_bounds)

        self.samples = saltelli.sample(self.problem, self.number_samples, calc_second_order=self.second_order)

        self.parsed_samples = parse_samples_to_run(self.samples,
                                                   self.parameter_names,
                                                   self.species_names)

        return self.parsed_samples

    def run_models(self):
        self.output = run_all_models(self.parsed_samples,
                                     self.model)

        self.model.reset_model()

        return self.output

    def get_outputs_at_timepoint(self, t, substrate):

        closest_timepoint = min(self.model.time, key=lambda x: abs(x - t))
        index = list(self.model.time).index(closest_timepoint)

        outputs = []
        for y in self.output:
            output = y[index][self.model.species_names.index(substrate)]
            outputs.append(output)

        self.output_for_analysis = np.array(outputs)

        return self.output_for_analysis

    def analyse_sobal_sensitivity_substrate_concentration_at_t(self, t, substrate):

        self.output_for_analysis = self.get_outputs_at_timepoint(t, substrate)

        self.analysis = sobol.analyze(self.problem,
                                      self.output_for_analysis,
                                      calc_second_order=self.second_order,
                                      num_resamples=self.num_resample,
                                      conf_level=self.conf_level,
                                      print_to_console=False,
                                      parallel=False,
                                      n_processors=None)

        rows = self.problem['names']
        self.dataframe_output = pd.DataFrame(self.analysis, index=rows)



        print("Saved as dataframe in self.dataframe_output")

        self.analysis_info = " --- Analysis mode = Uncertainty in substrate at timepoint --- \n"
        self.analysis_info += "Timepoint = " + str(t) + "\n"
        self.analysis_info += "Substrate = " + str(substrate) + "\n"

        lower, upper, max_out, min_out = self.get_interquartile_of_output_at_timepoint()
        self.analysis_info += ("\n 95 % Quartiles at t = " + str(lower) + " : " + str(upper) + "\n")
        self.analysis_info += ("Min at t = " + str(min_out) + ",  Max at t = " + str(max_out) + "\n")

        return self.dataframe_output

    def get_interquartile_of_output_at_timepoint(self, quartile=95):

        upper = np.percentile(self.output_for_analysis[1:], quartile)
        lower = np.percentile(self.output_for_analysis[1:], 100 - quartile)
        max_out = np.max(self.output_for_analysis[1:])
        min_out = np.min(self.output_for_analysis[1:])

        return lower, upper, max_out, min_out

    def return_sa_info(self):
        info = ""

        info += ("---Species concentration bounds--- \n")
        for i in range(len(self.species_names)):
            info += (str(self.species_names[i]) + " : " + str(self.species_bounds[i]) + "\n")

        info += ("\n")

        info += "--- Parameters --- \n"
        for i in range(len(self.parameter_bounds)):
            info += (str(self.parameter_names[i]) + " : " + str(self.parameter_bounds[i]) + "\n")

        info += "\n"

        info += self.analysis_info

        return info
