from SALib.sample import latin, saltelli
from SALib.analyze import sobol
from Kinetics.Model import *
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



""" Process multiple model outputs for uncertainty analysis"""
def return_quartiles(ys_for_single_substrate, name, quartile=95):
    quartiles = [["Time"], [str(name) + " High"], [str(name) + " Low"], [str(name) + " Mean"]]

    for i in range(len(ys_for_single_substrate)):
        # output_at_t will be a array.  i=0 is time, after than the substrate values at that time.
        output_at_t = np.array(ys_for_single_substrate[i])

        quartiles[0].append(output_at_t[0])
        quartiles[1].append(np.percentile(output_at_t[1:], quartile))
        quartiles[2].append(np.percentile(output_at_t[1:], 100 - quartile))
        quartiles[3].append(np.mean(output_at_t[1:]))

    return quartiles

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





"""Classes to do UA or SA"""

class UA():
    def __init__(self, parameters_with_bounds, species_with_bounds, model,
                 num_samples=1000, quartile_range=95):

        self.parameters_with_bounds = parameters_with_bounds
        self.species_with_bounds = species_with_bounds

        self.parameter_names, self.parameter_bounds = setup_bounds_lists(parameters_with_bounds)
        self.species_names, self.species_bounds = setup_bounds_lists(species_with_bounds)

        self.model = model
        self.num_samples = num_samples
        self.quartile_range = quartile_range

        self.samples = []
        self.parsed_samples = []
        self.problem = 0

        self.ua_output = []
        self.quartile_output = {}


    def analyse_number_of_samples_vs_parameters(self):
        # Code for text output
        num_params = len(self.parameter_names)
        num_specs = len(self.species_names)
        total = num_params + num_specs
        total_pw = total * total
        print("")
        print(num_params, 'parameters and', num_specs, "species in uncertainty analysis")
        print(total, "variables in total")
        print(total, "^2 =", total_pw)
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
        self.ua_output = run_all_models(self.parsed_samples, self.model)
        print("samples run")

        return self.ua_output

    def calculate_quartiles(self, quartile_range=95):

        for name in self.model.species_names:
            collected_runs = collect_runs_for_substrate(self.model.time, self.ua_output, self.model.species_names, name)
            quartiles = return_quartiles(collected_runs, name, quartile=quartile_range)

            self.quartile_output[name] = quartiles

        print("quartiles calculated")

        return self.quartile_output

    def print_ua_quartile_output(self, names_to_show):

        for name in names_to_show:
            print()
            print("--- " + str(name) + " Uncertainty Analysis Quartiles ---")
            for i in range(len(self.quartile_output[name][0])):
                for j in range(len(self.quartile_output[name])):
                    print(str(self.quartile_output[name][j][i]) + ", ", end="")
                print()
            print()

    def save_ua_quartile_output(self, names_to_show, filename=''):

        import os
        directory = os.path.dirname('output/')

        try:
            os.stat(directory)
        except:
            os.mkdir(directory)

        if filename == '':
            filename = os.path.basename(__file__)
            filename = filename[0:-3] + '.txt'

        filename = directory + "/" + filename

        file = open(filename, "w")

        file.write("Date: " + str(datetime.datetime.now()) + str("\n"))
        file.write("\n")

        file.write("---Species concentration bounds--- \n")
        for i in range(len(self.species_names)):
            file.write(str(self.species_names[i]) + " : " + str(self.species_bounds[i]) + "\n")

        file.write("\n")

        file.write("--- Parameters --- \n")
        for i in range(len(self.parameter_bounds)):
            file.write(str(self.parameter_names[i]) + " : " + str(self.parameter_bounds[i]) + "\n")

        file.write("\n")

        file.write("Quartile range = " + str(self.quartile_range) + "\n")
        file.write("Number of samples = " + str(self.num_samples) + "\n")
        file.write("\n")
        num_params = len(self.parameter_names)
        num_specs = len(self.species_names)
        total = num_params + num_specs
        file.write(str(num_params) + ' parameters and ' + str(num_specs) + " species in uncertainty analysis" + "\n")
        file.write(str(total) + " variables in total" + "\n")
        file.write(str(total) + "^2 = " + str(total * total) + "\n")
        file.write(str(total) + "^3 = " + str(total * total * total) + "\n")
        file.write(str(self.num_samples) + " samples made by lhc" + "\n")


        file.write("\n")

        for name in names_to_show:
            file.write("\n")
            file.write("--- " + str(name) + " Uncertainty Analysis Quartiles ---")
            file.write("\n")
            for i in range(len(self.quartile_output[name][0])):
                for j in range(len(self.quartile_output[name])):
                    file.write(str(self.quartile_output[name][j][i]) + ", ")
                file.write("\n")
            file.write("\n")

        file.close()

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
        self.problem = {}

        self.output = []
        self.output_for_analysis = []
        self.analysis = []

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
        return self.output

    def get_outputs_at_timepoint(self):

        closest_timepoint = min(self.model.time, key=lambda x: abs(x - self.timepoint_for_analysis))
        index = list(self.model.time).index(closest_timepoint)

        outputs = []
        for y in self.output:
            output = y[index][self.model.species_names.index(self.substrate_for_analysis)]
            outputs.append(output)

        self.output_for_analysis = np.array(outputs)

        return self.output_for_analysis

    def analyse_sobal_sensitivity(self):

        self.output_for_analysis = self.get_outputs_at_timepoint()

        self.analysis = sobol.analyze(self.problem,
                         self.output_for_analysis,
                         calc_second_order=self.second_order,
                         num_resamples=self.num_resample,
                         conf_level=self.conf_level,
                         print_to_console=True,
                         parallel=False,
                         n_processors=None)

        return self.analysis

    def get_interquartile_of_output_at_timepoint(self, quartile=95):

        upper = np.percentile(self.output_for_analysis[1:], quartile)
        lower = np.percentile(self.output_for_analysis[1:], 100 - quartile)

        return lower, upper

    def save_sa_quartile_output(self, filename=''):

        import os
        directory = os.path.dirname('output/')

        try:
            os.stat(directory)
        except:
            os.mkdir(directory)

        if filename == '':
            filename = os.path.basename(__file__)
            filename = filename[0:-3] + '.txt'

        filename = directory + "/" + filename

        file = open(filename, "w")

        file.write("Date: " + str(datetime.datetime.now()) + str("\n"))
        file.write("\n")

        file.write("---Species concentration bounds--- \n")
        for i in range(len(self.species_names)):
            file.write(str(self.species_names[i]) + " : " + str(self.species_bounds[i]) + "\n")

        file.write("\n")

        file.write("--- Parameters --- \n")
        for i in range(len(self.parameter_bounds)):
            file.write(str(self.parameter_names[i]) + " : " + str(self.parameter_bounds[i]) + "\n")

        file.write("\n")

        file.write("--- Sensitivity Analysis Parameters --- \n")
        file.write("Timepoint for analysis: " + str(self.timepoint_for_analysis) + "\n")
        file.write("Substrate for analysis: " + str(self.substrate_for_analysis) + "\n")

        file.write("\n")

        lower, upper = self.get_interquartile_of_output_at_timepoint()
        file.write("5th and 95th quartiles = " + str(lower) + " : " + str(upper))

        file.write("\n")
        file.write("\n")

        file.write("Parameter, S1, S1_conf, ST, ST conf")
        file.write("\n")
        for i in range(len(self.problem['names'])):
            file.write(self.problem['names'][i] + ", ")
            file.write(str(self.analysis['S1'][i]) + ", ")
            file.write(str(self.analysis['S1_conf'][i]) + ", ")
            file.write(str(self.analysis['ST'][i]) + ", ")
            file.write(str(self.analysis['ST_conf'][i]))
            file.write("\n")

        file.close()


