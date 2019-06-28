from kinetics.Uncertainty import check_not_neg, run_all_models
from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_samples(samples, parameter_names, species_names):
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

def make_salib_problem_with_bounds(model, negative_allowed=[], ppf=(0,1)):
    names = list(model.parameter_distributions.keys()) + list(model.species_distributions.keys())

    bounds = []

    for name, distribution in model.parameter_distributions.items():
        upper = distribution.ppf(ppf[1])
        lower = None
        lower_ppf = ppf[0]
        while not check_not_neg(lower, name, negative_allowed):
            lower = distribution.ppf(lower_ppf)
            lower_ppf += 0.01
        bounds.append([lower, upper])

    for name, distribution in model.species_distributions.items():
        upper = distribution.ppf(ppf[1])
        lower = None
        lower_ppf = ppf[0]
        while not check_not_neg(lower, name, negative_allowed):
            lower = distribution.ppf(lower_ppf)
            lower_ppf += 0.01
        bounds.append([lower, upper])

    problem = {'num_vars': len(names),
               'names': names,
               'bounds': bounds}

    return problem

def make_saltelli_samples(model, salib_problem, num_samples, second_order=False):

    saltelli_samples = saltelli.sample(salib_problem, num_samples, calc_second_order=second_order)

    samples = parse_samples(saltelli_samples, list(model.parameter_distributions.keys()), list(model.species_distributions.keys()))

    return samples

def get_concentrations_at_timepoint(model, output, timepoint, substrate):
    closest_timepoint = min(model.time, key=lambda x: abs(x - timepoint))
    index = list(model.time).index(closest_timepoint)

    outputs_for_analysis = []
    for y in output:
        output_at_t = y[index][model.run_model_species_names.index(substrate)]
        outputs_for_analysis.append(output_at_t)

    outputs_for_analysis = np.array(outputs_for_analysis)

    return outputs_for_analysis

def get_time_to_concentration(model, output, concentration, substrate, mode='>='):
    list_of_times = []

    for y in output:
        y = np.transpose(y)

        substrate_index = model.run_model_species_names.index(substrate)
        y_substrate = y[substrate_index]

        if mode == '<=':
            index_for_conc = np.where(y_substrate <= concentration)
        elif mode == '>=':
            index_for_conc = np.where(y_substrate >= concentration)

        if len(index_for_conc[0]) == 0:
            time = model.time[-1]
        else:
            index_for_conc = index_for_conc[0][0]
            time = model.time[index_for_conc]

        list_of_times.append(time)

    list_of_times = np.array(list_of_times)
    return list_of_times

def analyse_sobal_sensitivity(salib_problem, output_to_analyse,
                              second_order=False, num_resample=100, conf_level=0.95):

    analysis = sobol.analyze(salib_problem,
                             output_to_analyse,
                             calc_second_order=second_order,
                             num_resamples=num_resample,
                             conf_level=conf_level,
                             print_to_console=False,
                             parallel=False,
                             n_processors=None)

    rows = salib_problem['names']
    dataframe_output = pd.DataFrame(analysis, index=rows)

    return dataframe_output

def remove_st_less_than(dataframe, column='ST', less_than=0.001):
    new_df = dataframe[dataframe[column] > less_than]

    return new_df

def plot_sa_total_sensitivity(df):
    df.sort_values("ST", inplace=True, ascending=False)

    x_names = df.index.values
    x = np.arange(len(x_names))
    st = df['ST']
    st_err = df['ST_conf']

    plt.bar(x, st, align='center', yerr=st_err, edgecolor='black', color='#000090')
    plt.xticks(x, x_names, rotation=90)
    plt.ylabel("ST")



