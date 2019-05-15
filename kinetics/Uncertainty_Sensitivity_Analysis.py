from SALib.sample import latin, saltelli
from SALib.analyze import sobol
from kinetics.Model import *
import numpy as np
import pandas as pd
from tqdm import tqdm
import time
import random

""" Setup the bounds for uncertainty or sensitivity analsyis"""
def setup_bounds_lists(dict_with_bounds):
    """
    Converts a dict with bounds to two ordered lists

    Takes a dictinary of format {"name" : (lower_bound, upper_bound), ..}
    Returns two ordered lists of the same order:
      A list of the names - ['name1', 'name2'..]
      A list of the bounds - [ (lower_bound1, upper_bound1), (lower_bound2, upper_bound2), ..]

    :param dict_with_bounds:  {"name" : (lower_bound, upper_bound), ..}
    :return: names_list, bounds_list
    """

    names_list = []
    bounds_list = []

    for name in dict_with_bounds:
        if dict_with_bounds[name] != 0:
            names_list.append(name)
            bounds_list.append(dict_with_bounds[name])

    return names_list, bounds_list

def get_bounds_from_std_error(dict_with_std_error):
    """
    Converts dict with std to dict with lower and upper bounds.

    Takes a dictionary with the format {"name" : (value, std_error_of_value), ..}

    Returns a dictinary with the format {"name" : (lower_bound, upper_bound), ..}

    :param dict_with_std_error: {"name" : (value, std_error_of_value), ..}
    :return: dict_with_bounds: {"name" : (lower_bound, upper_bound), ..}
    """

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
    """
    Converts dict with percentage error to dict with lower and upper bounds.

    Takes a dictionary with the format {"name" : (value, percentage_error_of_value), ..}

    Returns a dictinary with the format {"name" : (lower_bound, upper_bound), ..}

    :param dict_with_pc_error: {"name" : (value, percentage_error_of_value), ..}
    :return: dict_with_bounds: {"name" : (lower_bound, upper_bound), ..}
    """

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
    """
    This function returns a 'problem' dictionary which SALib uses.

    Takes ordered lists for parameters and species, for the names and the bounds.
    Combines these to make a single pair of ordered lists
    Uses these to create the problem dictionary.


    :param parameter_names:  An ordered list of parameter names. ['name1', name2', ..]
                             Matches order of parameter_bounds

    :param parameter_bounds: An ordered list of parameter_bounds. [(lower1, upper1), (lower2, upper2), ..]
                             Matches order of parameter_names

    :param species_names:  An ordered list of species names. ['name1', name2', ..]
                           Matches order of species_bounds

    :param species_bounds: n ordered list of species_bounds. [(lower1, upper1), (lower2, upper2), ..]
                           Matches order of species_names

    :return: a dictinary with the format:
                                        {'num_vars': len(names_list),
                                         'names': names_list,
                                         'bounds': bounds}

    """

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
    """
    Takes a list of samples and converts this to a list of tuples of dictionaries
     which can update the Model class.

    Samples are created using the relevent function (eg LHC or Saltelli)
    Samples will be a list.
      Each entry in the list, is a list of values for the species and parameters,
      in the order they were defined in the 'problem' dictionary.
      This will be the parameters first followed by the species.

    For example, Samples = [ [value1, value2, value3..],  [value1, value2, value3], ...]

    The ordered lists of parameter and species names are used to unpack these lists back into dictionaries.
    These dictionaries can be used to update the model.
    The dictionaries are saved in a new list called parsed samples

    Parsed samples = [ (parameter_dict1, species_dict1), (parameter_dict2, species_dict2) ..]


    :param samples:  a list of lists containing the samples values, in the order they were defined in problem dict.
    :param parameter_names: ordered list of parameter names
    :param species_names: ordered list of species names
    :return: Parsed samples - [ (parameter_dict1, species_dict1), (parameter_dict2, species_dict2) ..]
    """

    parsed_samples = []

    def unpack_sampled_params_and_species(sample_list, parameter_names_list, species_names_list):
        """
        Returns a dictionaries for parameters and species

        :param sample_list: a list of sample values,
                            parameters first then species, of the same order as the names lists
        :param parameter_names_list: an ordered list of parameter_names which matches the first part of sample_list
        :param species_names_list:  an ordered list of species_names which matches the first part of sample_list
        :return: (parameters_dictionary, species_dictionary)   - dictionary format is {'name' : value, ..}
        """

        parameters_dict = {}
        count = 0
        for i in range(len(parameter_names_list)):
            name = parameter_names_list[i]
            parameters_dict[name] = sample_list[i]
            count += 1

        species_dict = {}
        for i in range(count, len(species_names_list) + count):
            name = species_names_list[i - count]
            species_dict[name] = sample_list[i]

        return parameters_dict, species_dict

    for sampled in samples:
        parameters, species = unpack_sampled_params_and_species(sampled, parameter_names, species_names)
        parsed_samples.append([parameters, species])

    # returns a list of tuples containing [(parameter_dict, species_dict), ] for each sample
    return parsed_samples


def run_all_models(parsed_samples, model, logging=True):
    """
    Run all the models for a set of parsed samples. Return the outputs

    Takes the parsed samples, runs each model in turn.  Saves each model output to a list of outputs.
    Returns the list of outputs  - [ys1, ys2, ys3] - where ys is the numpy array of y values from a model run)


    :param parsed_samples: a list of [ (parameter_dict1, species_dict1), (parameter_dict2, species_dict2), ..]
    :param model:  an instance of the Model class
    :return: list of outputs  - [ys1, ys2, ys3] - where ys is a numpy array of y values from a model run)
    """

    output = []

    if logging==True:
        for parameters, species in tqdm(parsed_samples):
            model.update_species(species)
            model.update_parameters(parameters)
            y = model.run_model()
            output.append(y)

    elif logging==False:
        for parameters, species in parsed_samples:
            model.update_species(species)
            model.update_parameters(parameters)

            y = model.run_model()
            output.append(y)

    # Reset the model back to the default values
    model.reset_model()

    # ua_output will be a list like [y1, y2, y3, ect...]
    return output


""" Process multiple model outputs for uncertainty analysis"""
def return_ys_for_a_single_substrate(model, run_all_models_ouput, substrate_name):
    """
    From run_all_models_ouput, return only the outputs for substrate_name

    collected_output = [ [t0, r1, r2, r3],
                         [t1, r1, r2, r3].... ]


    :param model: an instance of the Model class
    :param run_all_models_ouput: the output from the run_all_models function - [ys1, ys2, ys3 .. ]
    :param substrate_name: he name of the substrate to return the ys for
    :return: collected_output = [ [t0, r1, r2, r3],   [t1, r1, r2, r3].... ]
    """

    collected_output = []

    for i in range(len(model.time)):
        timepoint = [model.time[i]]

        for y in run_all_models_ouput:
            timepoint.append(y[i][model.species_names.index(substrate_name)])

        collected_output.append(timepoint)

    return collected_output

def return_quartiles(ys_for_single_substrate, quartile=95):
    """
    Takes ys_for_single_substrate and returns only the quartiles

    ys_for_single_substrate from the function return_ys_for_a_single_substrate

    Runs through each timepoint in turn, and takes the high, low and mean.
    High and Low determined by quartile=95

    Returns quartiles, which is a dictinary containing the time, high, low and mean outputs as lists.


    :param ys_for_single_substrate:  from the function return_ys_for_a_single_substrate
                                     has the format: [ [t0, r1, r2, r3], [t1, r1, r2, r3].. ]

    :param quartile: default to 95.  Is the quartile which is taken at each timepoint (ie 5th and 95th)
    :return: quartiles dictionary - {"Time": [], "High": [], "Low": [], "Mean": []}  Outputs saved as lists.
    """

    quartiles = {"Time": [], "High": [], "Low": [], "Mean": []}

    for i in range(len(ys_for_single_substrate)):
        # output_at_t will be a array.  i=0 is time, after than the substrate values at that time.
        output_at_t = ys_for_single_substrate[i]

        quartiles["Time"].append(output_at_t[0])
        quartiles["High"].append(np.percentile(output_at_t[1:], quartile))
        quartiles["Low"].append(np.percentile(output_at_t[1:], 100 - quartile))
        quartiles["Mean"].append(np.mean(output_at_t[1:]))

    return quartiles

def random_sample(problem, number_of_samples):
    """
    Generate model inputs using random sampling

    Returns a NumPy matrix containing the model inputs generated by random sampling.
    The resulting matrix contains N rows and D columns, where D is the number of parameters.

    :param problem:  {'num_vars': len(names_list),
                      'names': names_list,
                      'bounds': bounds}

    :param number_of_samples: a integer, the number of samples to make

    :return: a list of lists.  Each list contains sampled values with the same order as the problem.
    """

    samples = []
    for i in range(number_of_samples):
        new_sample = []
        for pair_of_bounds in problem['bounds']:
            new_sample.append(random.uniform(pair_of_bounds[0], pair_of_bounds[1]))

        samples.append(new_sample)

    samples = np.array(samples)

    return samples

def check_bounds(bounds_names, bounds_tuples):
    for name, tuple in zip(bounds_names, bounds_tuples):
        try:
            assert tuple[0] < tuple[1]
        except AssertionError as error:
            print(str(name) + ' bounds are incorrect, check smaller bound is first')






"""Classes to do UA or SA"""


class UA(object):

    def __init__(self, model, num_samples=1000, quartile_range=95, logging=True):

        self.parameters_with_bounds = model.parameter_bounds
        self.species_with_bounds = model.species_bounds

        self.parameter_names, self.parameter_bounds = setup_bounds_lists(self.parameters_with_bounds)
        self.species_names, self.species_bounds = setup_bounds_lists(self.species_with_bounds)

        self.model = model
        self.num_samples = num_samples
        self.quartile_range = quartile_range

        self.samples = []
        self.parsed_samples = []
        self.problem = {}

        self.output = []
        self.quartile_output = {}

        self.substrate_dataframes = {}

        self.all_runs_substrate_dataframes = {}

        self.logging = logging

    def load_species_and_parameters_from_model(self):
        self.parameters_with_bounds = self.model.parameter_bounds
        self.species_with_bounds = self.model.species_bounds

        self.parameter_names, self.parameter_bounds = setup_bounds_lists(self.parameters_with_bounds)
        self.species_names, self.species_bounds = setup_bounds_lists(self.species_with_bounds)


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
        check_bounds(self.parameter_names, self.parameter_bounds)
        check_bounds(self.species_names, self.species_bounds)

        self.problem = setup_problem(self.parameter_names, self.parameter_bounds,
                                     self.species_names, self.species_bounds)

        self.samples = latin.sample(self.problem, self.num_samples)

        self.parsed_samples = parse_samples_to_run(self.samples, self.parameter_names, self.species_names)

        if self.logging==True:
            print("self.problem, self.samples and self.parsed_samples set by lhc")

        return self.parsed_samples

    def make_random_samples(self):
        self.problem = setup_problem(self.parameter_names, self.parameter_bounds,
                                     self.species_names, self.species_bounds)

        self.samples = random_sample(self.problem, self.num_samples)

        self.parsed_samples = parse_samples_to_run(self.samples, self.parameter_names, self.species_names)

    def run_models(self):
        if self.logging == True:
            print("running all models")

        time.sleep(0.5)
        self.output = run_all_models(self.parsed_samples, self.model, logging=self.logging)

        if self.logging == True:
            print("samples run, model outputs saved in self.output")

        self.model.reset_model()

        return self.output

    def quartiles_to_dataframe(self):
        self.substrate_dataframes = {}
        for substrate in self.model.species_names:
            df = pd.DataFrame(self.quartile_output[substrate], columns=["Time", "High", "Low", "Mean"])
            self.substrate_dataframes[substrate] = df

        if self.logging == True:
            print("Quartiles for each substrate saved to self.substrate_dataframes")
        return self.substrate_dataframes

    def calculate_quartiles(self):

        for name in self.model.species_names:
            collected_runs = return_ys_for_a_single_substrate(self.model, self.output, name)
            quartiles = return_quartiles(collected_runs, quartile=self.quartile_range)

            self.quartile_output[name] = quartiles

        if self.logging == True:
            print("quartiles calculated, saved in self.quartile_output")

        self.quartiles_to_dataframe()

        return self.substrate_dataframes

    def calculate_df_quartiles_single_substrate(self, substrate):
        collected_runs = return_ys_for_a_single_substrate(self.model, self.output, substrate)
        quartiles = return_quartiles(collected_runs, quartile=self.quartile_range)

        df = pd.DataFrame(quartiles, columns=["Time", "High", "Low", "Mean"])

        return df



    def calculate_all_runs_substrate_dataframes(self):

        for name in self.model.species_names:
            # format: [[t0, r1, r2, r3], [t1, r1, r2, r3]..]
            collected_runs = return_ys_for_a_single_substrate(self.model, self.output, name)
            all_runs = {"Time": []}
            column_titles = ['Time']

            for i in range(1,len(collected_runs[0])):
                all_runs[str(i)] = []
                column_titles.append(str(i))

            for timepoint in collected_runs:
                all_runs['Time'].append(timepoint[0])

                for i in range(1,len(timepoint)):
                    all_runs[str(i)].append(timepoint[i])

            df = pd.DataFrame(all_runs, columns=column_titles)

            self.all_runs_substrate_dataframes[name] = df

        return self.all_runs_substrate_dataframes

    def run_standard_ua(self, output=None):   #output = 'all' or 'quartiles' or 'None'
        self.make_lhc_samples()
        self.run_models()
        self.calculate_all_runs_substrate_dataframes()
        self.calculate_quartiles()

        if output == 'all':
            return self.all_runs_substrate_dataframes
        elif output == 'quartiles':
            return self.substrate_dataframes
        else:
            return None

    def return_ua_info(self):
        info = "--- Species Defaults in Reaction --- \n"
        for name in self.model.species_names:
            if self.model.species_defaults[name] > 0:
                info += (str(name) + " : " + str(self.model.species_defaults[name]) + "\n")

        info += "\n ---Species Concentration bounds--- \n"
        for i in range(len(self.species_names)):
            info += (str(self.species_names[i]) + " : " + str(self.species_bounds[i]) + "\n \n")

        info += "--- Parameters Bounds--- \n"
        for i in range(len(self.parameter_bounds)):
            info += (str(self.parameter_names[i]) + " : " + str(self.parameter_bounds[i]) + "\n \n")

        info += ("Quartile range = " + str(self.quartile_range) + "%" + "\n")
        info += ("Number of samples = " + str(self.num_samples) + "\n \n")
        num_params = len(self.parameter_names)
        num_specs = len(self.species_names)
        total = num_params + num_specs
        info += (str(num_params) + ' parameters and ' + str(num_specs) + " species in uncertainty analysis" + "\n")
        info += (str(total) + " variables in total" + "\n")
        info += (str(total) + "^2 = " + str(total * total) + "\n")
        info += (str(total) + "^3 = " + str(total * total * total) + "\n")
        info += (str(self.num_samples) + " samples made by lhc" + "\n")

        return info

class SA(object):
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

        print('self.problem, self.samples and self.parsed_samples set by make_saltelli_samples')

        return self.parsed_samples

    def run_models(self):
        print("running all models")
        time.sleep(0.5)
        self.output = run_all_models(self.parsed_samples,
                                     self.model)
        print("samples run, model outputs saved in self.output")

        self.model.reset_model()

        return self.output

    def get_outputs_at_timepoint(self, timepoint, substrate):

        closest_timepoint = min(self.model.time, key=lambda x: abs(x - timepoint))
        index = list(self.model.time).index(closest_timepoint)

        outputs = []
        for y in self.output:
            output = y[index][self.model.species_names.index(substrate)]
            outputs.append(output)

        self.output_for_analysis = np.array(outputs)
        print("self.output_for_analysis updated with outputs for substrate",
              str(substrate), "at timepoint", str(timepoint))

        return self.output_for_analysis

    def get_times_to_reach_substrate_concentration(self, substrate_name, substrate_concentration, mode='<='):



        list_of_times = []

        for y in self.output:
            y = np.transpose(y)

            substrate_index = self.model.species_names.index(substrate_name)
            y_substrate = y[substrate_index]

            if mode == '<=':
                index_for_conc = np.where(y_substrate <= substrate_concentration)
            elif mode == '>=':
                index_for_conc = np.where(y_substrate >= substrate_concentration)

            if len(index_for_conc[0]) == 0:
                time = self.model.time[-1]
            else:
                index_for_conc = index_for_conc[0][0]
                time = self.model.time[index_for_conc]

            list_of_times.append(time)


        list_of_times = np.array(list_of_times)
        return list_of_times

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

        print("self.analysis updated with sobal sensitivity analysis output")

        self.dataframe_output = pd.DataFrame(self.analysis, index=rows)
        print("Sobal sensitivity analysis saved as dataframe in self.dataframe_output")

        self.analysis_info = " --- Analysis mode = Uncertainty in substrate at timepoint --- \n"
        self.analysis_info += "Timepoint = " + str(t) + "\n"
        self.analysis_info += "Substrate = " + str(substrate) + "\n"

        lower, upper, max_out, min_out = self.get_interquartile_of_output_at_timepoint()
        self.analysis_info += ("\n 95 % Quartiles at t = " + str(lower) + " : " + str(upper) + "\n")
        self.analysis_info += ("Min at t = " + str(min_out) + ",  Max at t = " + str(max_out) + "\n")

        return self.dataframe_output

    def analysis_sobal_sensitivity_time_to_concentration(self, substrate_name, substrate_concentration, mode='<='):

        self.output_for_analysis = self.get_times_to_reach_substrate_concentration(substrate_name,
                                                                                   substrate_concentration,
                                                                                   mode=mode)

        self.analysis = sobol.analyze(self.problem,
                                      self.output_for_analysis,
                                      calc_second_order=self.second_order,
                                      num_resamples=self.num_resample,
                                      conf_level=self.conf_level,
                                      print_to_console=False,
                                      parallel=False,
                                      n_processors=None)

        rows = self.problem['names']

        print("self.analysis updated with sobal sensitivity analysis output")

        self.dataframe_output = pd.DataFrame(self.analysis, index=rows)
        print("Sobal sensitivity analysis saved as dataframe in self.dataframe_output")

        self.analysis_info = " --- Analysis mode = Uncertainty in time to target substrate conc --- \n"

        return self.dataframe_output

    def get_interquartile_of_output_at_timepoint(self, quartile=95):

        upper = np.percentile(self.output_for_analysis[1:], quartile)
        lower = np.percentile(self.output_for_analysis[1:], 100 - quartile)
        max_out = np.max(self.output_for_analysis[1:])
        min_out = np.min(self.output_for_analysis[1:])

        return lower, upper, max_out, min_out

    def return_sa_info(self):
        info = ""

        info += "---Species concentration bounds--- \n"
        for i in range(len(self.species_names)):
            info += (str(self.species_names[i]) + " : " + str(self.species_bounds[i]) + "\n \n")

        info += "--- Parameters --- \n"
        for i in range(len(self.parameter_bounds)):
            info += (str(self.parameter_names[i]) + " : " + str(self.parameter_bounds[i]) + "\n")

        info += "\n"

        info += self.analysis_info

        return info
