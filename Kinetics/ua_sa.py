from SALib.sample import latin, saltelli
from SALib.analyze import sobol
from Kinetics.kinetics import *
import numpy as np

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


""" Sampling - returns a list of ordered lists containing sampled parameters and sampled species concentrations"""


def compile_params_species(parameter_names_list, parameter_bounds_list, species_names_list, species_bounds_list):
    names_list = []
    names_list.extend(parameter_names_list)
    names_list.extend(species_names_list)

    bounds = []
    bounds.extend(parameter_bounds_list)
    bounds.extend(species_bounds_list)

    return names_list, bounds

def make_samples_latin_hypercube(parameter_names_list, parameter_bounds_list, species_names_list, species_bounds_list,
                                 number):
    names_list, bounds = compile_params_species(parameter_names_list, parameter_bounds_list,
                                                species_names_list, species_bounds_list)

    problem = {'num_vars': len(names_list),
               'names': names_list,
               'bounds': bounds}

    sampled_inputs = latin.sample(problem, number)

    return sampled_inputs, problem

def make_samples_for_sobal(parameter_names_list, parameter_bounds_list, species_names_list, species_bounds_list,
                           number=500, second_order=False):
    names_list, bounds = compile_params_species(parameter_names_list, parameter_bounds_list,
                                                species_names_list, species_bounds_list)

    problem = {'num_vars': len(names_list),
               'names': names_list,
               'bounds': bounds}

    sampled_inputs = saltelli.sample(problem, number, calc_second_order=second_order)

    return sampled_inputs, problem


""" Run uncertainty analysis or sensitivity analysis"""


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

def parse_samples_to_run(samples, parameter_names, species_names):
    param_species_samples = []
    for sampled in samples:
        parameters, species = unpack_sampled_params_and_species(sampled, parameter_names, species_names)
        param_species_samples.append([parameters, species])

    # returns a list of tuples containing [(parameter_dict, species_dict), ] for each sample
    return param_species_samples

def run_all_models(parsed_samples, model):
    ua_sa_output = []
    count = 0
    for ua_parameters, ua_species in parsed_samples:
        model.update_species_names_and_starting_values(ua_species)
        model.update_parameters(ua_parameters)
        count += 1

        if count % 50 == 0:
            print(count, "out of", len(parsed_samples), "complete")

        y = model.run_model()

        ua_sa_output.append(y)

    # ua_output will be a list like [y1, y2, y3, ect...]

    return ua_sa_output


""" Process uncertainty analysis"""


def find_quartiles(ys_for_single_substrate, name, quartile=95):
    quartiles = [["Time"], [str(name) + " High"], [str(name) + " Low"], [str(name) + " Mean"]]

    for i in range(len(ys_for_single_substrate)):
        # output_at_t will be a array.  i=0 is time, after than the substrate values at that time.
        output_at_t = np.array(ys_for_single_substrate[i])

        quartiles[0].append(output_at_t[0])
        quartiles[1].append(np.percentile(output_at_t[1:], quartile))
        quartiles[2].append(np.percentile(output_at_t[1:], 100 - quartile))
        quartiles[3].append(np.mean(output_at_t[1:]))

    final_output = quartiles

    return final_output

def collect_runs_for_substrate(time, list_y, species_names, name):
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

def process_ua_output(ua_output, time, species_names, quartile_range=95):
    output = {}

    for name in species_names:
        collected_runs = collect_runs_for_substrate(time, ua_output, species_names, name)
        quartiles = find_quartiles(collected_runs, name, quartile=quartile_range)

        output[name] = quartiles

    return output


""" Process sensitivity analysis"""


def get_outputs_for_sobal(ua_sa_output, species_names, name, time_point, time):
    closest_timepoint = min(time, key=lambda x: abs(x - time_point))
    index = list(time).index(closest_timepoint)

    outputs = []
    for y in ua_sa_output:
        output = y[index][species_names.index(name)]
        outputs.append(output)

    return np.array(outputs)

def analyse_satelli_sobol(problem, outputs,
                          second_order=False, conf_level=0.95, num_resample=100):
    return sobol.analyze(problem,
                         outputs,
                         calc_second_order=second_order,
                         num_resamples=num_resample,
                         conf_level=conf_level,
                         print_to_console=True,
                         parallel=False,
                         n_processors=None)








"""Classes to do UA or SA"""

class UA():
    def __init__(self, parameters_with_bounds, species_with_bounds, model,
                 latin_hypercube_samples=1000, quartile_range=95):

        self.parameters_with_bounds = parameters_with_bounds
        self.species_with_bounds = species_with_bounds

        self.parameter_names, self.parameter_bounds = setup_bounds_lists(parameters_with_bounds)
        self.species_names, self.species_bounds = setup_bounds_lists(species_with_bounds)

        self.model = model
        self.latin_hypercube_samples = latin_hypercube_samples
        self.quartile_range = quartile_range

        self.samples = []
        self.parsed_samples = []
        self.problem = 0

        self.ua_output = []

    def make_lhc_samples(self):

        # Code for text output
        num_params = len(self.parameter_names)
        num_specs = len(self.species_names)
        total = num_params + num_specs
        total_pw = total * total
        print("")
        print(num_params, 'parameters and', num_specs, "species in uncertainty analysis")
        print(total, "variables in total")
        print(total, "^2 =", total_pw)
        print(str(self.latin_hypercube_samples), "samples made by lhc")
        print("")

        # Code for making samples
        self.samples, self.problem = make_samples_latin_hypercube(self.parameter_names, self.parameter_bounds,
                                                                  self.species_names, self.species_bounds,
                                                                  self.latin_hypercube_samples)
        self.parsed_samples = parse_samples_to_run(self.samples, self.parameter_names, self.species_names)

        return self.parsed_samples

    def run_models(self):
        print("running all models")
        self.ua_output = run_all_models(self.parsed_samples, self.model)
        print("samples run")

        return self.ua_output

    def calculate_quartiles(self, quartile_range=0):
        if quartile_range != 0:
            self.quartile_range=quartile_range

        self.quartile_output = process_ua_output(self.ua_output, self.model.time, self.model.species_names,
                                                 quartile_range=self.quartile_range)
        print("quartiles calculated")

        return self.quartile_output

    def print_ua_quartile_output(self, names_to_show):

        for name in names_to_show:
            print()
            print("--- " + str(name) + " Uncertainty Analysis Quartiles ---")
            for i in range(len(self.quartile_output[name][0])):
                for j in range(len(output[name])):
                    print(str(output[name][j][i]) + ", ", end="")
                print()
            print()


class SA():
    def __init__(self,
                 parameters_with_bounds,
                 species_with_bounds,
                 model,
                 timepoint_for_analysis,
                 substrate_for_analysis,
                 number_samples=500,
                 second_order=False,
                 conf_level=0.95,
                 num_resample=100):

        self.timepoint_for_analysis = timepoint_for_analysis
        self.substrate_for_analysis = substrate_for_analysis

        self.number_samples = number_samples
        self.second_order = second_order
        self.conf_level = conf_level
        self.num_resample = num_resample

        self.parameters_with_bounds = parameters_with_bounds
        self.species_with_bounds = species_with_bounds

        self.parameter_names, self.parameter_bounds = setup_bounds_lists(parameters_with_bounds)
        self.species_names, self.species_bounds = setup_bounds_lists(species_with_bounds)

        self.model = model

        self.samples = []
        self.parsed_samples = []
        self.problem = 0

        self.output = []
        self.analysis = []

    def make_saltelli_samples(self):
        self.samples, self.problem = make_samples_for_sobal(self.parameter_names,
                                                            self.parameter_bounds,
                                                            self.species_names,
                                                            self.species_bounds,
                                                            number=self.number_samples,
                                                            second_order=self.second_order)

        self.parsed_samples = parse_samples_to_run(self.samples,
                                                   self.parameter_names,
                                                   self.species_names)

        return self.parsed_samples

    def run_models(self):
        self.output = run_all_models(self.parsed_samples,
                                     self.model)
        return self.output

    def analyse_sobal_sensitivity(self):
        self.analysis = analyse_satelli_sobol(self.problem,
                                              self.output,
                                              second_order=self.second_order,
                                              conf_level=self.conf_level,
                                              num_resample=self.num_resample)

        return self.analysis

