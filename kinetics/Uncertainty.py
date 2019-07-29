from tqdm import tqdm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from SALib.sample import latin, saltelli
from SALib.analyze import sobol

""" -- Helper Functions, not for use by user -- """
def check_not_neg(sample, name, negative_allowed):
    if sample == None:
        return False

    if (sample < 0) and (name not in negative_allowed):
        return False

    return True

def return_ys_for_a_single_substrate(model, output, substrate_name):

    collected_output = []
    species_names = list(model.species.keys())

    for i in range(len(model.time)):
        timepoint = [model.time[i]]

        for y in output:
            timepoint.append(y[i][species_names.index(substrate_name)])

        collected_output.append(timepoint)

    return collected_output

def parse_samples(samples, parameter_names, species_names):
    """
    Takes a list of samples and converts this to a list of tuples of dictionaries which can update the Model class.

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


""" -- Generate Samples --"""
def sample_distributions(model, num_samples=1000, negative_allowed=[]):
    """
    Makes a set of samples from the species and parameter distributions in the model.

    Args:
        model (kinetics.Model): A model object
        num_samples (int): Number of samples to make (default 1000)
        negative_allowed (list): A list of any distributions that can be negative.

    Returns:
        A list of samples.  Each entry in the list is a tuple containing (parameter_dict, species_dict) for the samples.
    """

    samples = []
    for i in range(num_samples):
        parameter_dict = {}
        species_dict = {}

        for name, distribution in model.parameter_distributions.items():
            sample = None
            while not check_not_neg(sample, name, negative_allowed):
                sample = distribution.rvs()
            parameter_dict[name] = sample

        for name, distribution in model.species_distributions.items():
            sample=None
            while not check_not_neg(sample, name, negative_allowed):
                sample = distribution.rvs()
            species_dict[name] = sample

        samples.append([parameter_dict, species_dict])

    return samples # samples will = [ (parameter_dict1, species_dict1), (parameter_dict2, species_dict2) ..]

def sample_uniforms(model, num_samples=1000):

    bounds = []
    names = []
    for name, tuple in model.parameter_distributions.items():
        bounds.append(tuple)
        names.append(names)

    for name, tuple in  model.species_distributions.items():
        bounds.append(tuple)
        names.append(names)

    problem = {'num_vars': len(names),
               'names': names,
               'bounds': bounds}

    lhc_samples = latin.sample(problem, num_samples)

    samples = parse_samples(lhc_samples, list(model.parameter_distributions.keys()), list(model.species_distributions.keys()))

    return samples

def salib_problem_with_bounds(model, negative_allowed=[], ppf=(0,1)):
    """
    Make a salib problem by specifying bounds using ppf of scipy distributions

    Args:
        model (Model): The model object
        negative_allowed (list): Any distributions which are allowed negative samples
        ppf (Tuple): Percent point functions to take.  Default is 0 and 1 which will give the absolute upper and lower bounds of a distribution.

    Returns:
        An SALib problem
    """

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
    """
    Use SALib to make saltelli samples

    Args:
        model (Model): The model object
        salib_problem (dict): An SALib Problem
        num_samples (int): number of samples to take
        second_order (bool): look at second order interactions

    Returns:
        Samples from salib which have been parsed into a set of samples which run_all_models can take.

    """

    saltelli_samples = saltelli.sample(salib_problem, num_samples, calc_second_order=second_order)

    samples = parse_samples(saltelli_samples, list(model.parameter_distributions.keys()), list(model.species_distributions.keys()))

    return samples


""" -- Run models"""
def run_all_models(model, samples, logging=True):
    """
    Run all the models for a set of samples.

    Args:
        model (kinetics.Model): A model object
        samples (list): A list of samples in the form [(param_dict1, species_dict1), (param_dict2.... ect}
        logging (bool): Show logging and progress bar.  Default = True

    Returns (list): [y1, y2, y3, y4, ect..]

    """
    output = []

    if logging==True:
        for parameters, species in tqdm(samples):
            model.update_species(species)
            model.run_model_parameters.update(parameters)
            y = model.run_model()
            output.append(y)

    elif logging==False:
        for parameters, species in samples:
            model.update_species(species)
            model.run_model_parameters.update(parameters)

            y = model.run_model()
            output.append(y)

    # Reset the model back to the default values
    model.reset_model_to_defaults()

    # ua_output will be a list like [y1, y2, y3, ect...]
    return output

def dataframes_all_runs(model, output, substrates=[]):
    """
    Gives a dictionary of dataframes - {'Substrate' : dataframe'}
    Each dataframe has time in the first column, followed by every model run in the subsequent columns

    [[t0, r1, r2, r3], [t1, r1, r2, r3]]

    Args:
        model (Model): Model object
        output (list): The output from run_all_models. [y1, y2, y3 ect]
        substrates (list): Substrate names to include. If empty returns all (default).

    Returns:
        Dictionary of dataframes containing all model runs - {'Substrate' : dataframe'}
    """

    all_runs_substrate_dataframes = {}

    if substrates == []:
        substrates = list(model.species.keys())

    for name in substrates:
        # format: [[t0, r1, r2, r3], [t1, r1, r2, r3]..]
        collected_runs = return_ys_for_a_single_substrate(model, output, name)
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

        all_runs_substrate_dataframes[name] = df

    return all_runs_substrate_dataframes

def dataframes_quartiles(model, output, substrates=[], quartile=95):
    """
    Gives a dictionary of dataframes - {'Substrate' : dataframe'}
    Each dataframe has columns ['Time', 'High', 'Low', 'Mean']

    Args:
        model (Model): Model object
        output (list): The output from run_all_models. [y1, y2, y3 ect]
        substrates (list): Substrate names to include. If empty returns all (default).
        quartile (int): The percentile to take.  Default is 95 which gives with 95% and 5% quartiles.

    Returns:
        Dictionary of dataframes containing confidence intervals from the uncertainty analysis.
    """

    dataframes = {}

    if substrates == []:
        substrates = list(model.species.keys())

    for name in substrates:

        quartiles = {"Time": [], "High": [], "Low": [], "Mean": []}
        ys_for_single_substrate = return_ys_for_a_single_substrate(model, output, name)

        for i in range(len(ys_for_single_substrate)):
            # output_at_t will be a array.  i=0 is time, after than the substrate values at that time.
            output_at_t = ys_for_single_substrate[i]

            quartiles["Time"].append(output_at_t[0])
            quartiles["High"].append(np.percentile(output_at_t[1:], quartile))
            quartiles["Low"].append(np.percentile(output_at_t[1:], 100 - quartile))
            quartiles["Mean"].append(np.mean(output_at_t[1:]))

        dataframes[name] = quartiles

    return dataframes


""" -- Plotting uncertainty analysis -- """
def plot_substrate(substrate, dataframes,
                   colour='blue', xlabel="Time (mins)", ylabel="μM",
                   alpha=0.1, linewidth=0.1, y_min=True, plot=False):
    """
    Plot every model run for a single substrate.

    Args:
        substrate (str): Substrate name
        dataframes (dict): A dictionary of dataframes made using dataframes_all_runs
        colour: Colour argument for matplotlib, default = 'blue'
        xlabel (str): Label for x axis, default = 'Time (mins)'
        ylabel (str): Label for y axis, default = 'μM'
        alpha: Alpha argument for matplotlib, default = 0.1
        linewidth: Linewidth argument for matplotlib, defualt = 0.1
        y_min (int): If a number sets the bottom of the axis to this. Default is True
        plot (bool):  If true plots the graph using plt.plot()

    """

    df = dataframes[substrate]
    for i in range(1, len(df.columns)):
        plt.plot(df['Time'], df[str(i)],
                 color=colour, alpha=alpha, linewidth=linewidth)
    print(str(substrate) + ' - ' + str(colour))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if y_min != True:
        plt.ylim(bottom=y_min)

    if plot == True:
        plt.show()

def plot_ci_intervals(substrates_to_add, dataframes, plot=False,
                      colours=['darkred', 'darkblue', 'darkgreen', 'darkorange'],
                      alpha=0.1, units=['$\mu M$', 'Time (mins)']):
    """
    Plot every model run for a single substrate.

    Args:
        substrates_to_add (list): List of substrate names
        dataframes (dict): A dictionary of dataframes made using dataframes_quartiles
        colours (list): Colour arguments for matplotlib, each substrate will cycle through this list.
        alpha (int): Alpha argument for matplotlib, default = 0.1
        units (list): Units for the axis [yaxis_lable, xaxis_lable]
        plot (bool):  If true plots the graph using plt.plot()
    """

    for i, substrate in enumerate(substrates_to_add):
        color = colours.pop(0)
        colours.append(color)

        df = dataframes[substrate]

        time = df['Time']
        high = df['High']
        low = df['Low']
        mean = df['Mean']

        # high_line = plt.plot(time, high, color=color, linestyle="--", linewidth = 0.5)
        # low_line = plt.plot(time, low,  color=color, linestyle="--", linewidth = 0.5)
        plt.plot(time, mean, color=color, linewidth=0.5, label=substrate)
        plt.fill_between(time, high, y2=low, color=color, alpha=alpha, linewidth=0)

        plt.ylabel(units[0])
        plt.xlabel(units[1])

    if plot == True:
        plt.show()

def plot_data(substrates, data_df, plot=False,
              alpha=0.5, size=35, colours=['black'], symbols=["o", "s", '^', 'v']):
    """
    Add experimental data to a plot.

    Args:
        substrates (list): A list of substrate names
        data_df (Dataframe): A pandas dataframe containing experimental data
        plot (bool): If True calls plt.plot()
        alpha (int): Alpha argument for matplotlib, default = 0.1
        size (int): Size argument for matplotlib, default = 35
        colours (list): List of colours for matplotlib. Default=['black']
        symbols (list): Symbols for matpltlib. Defaut=["o", "s", '^', 'v']
    """

    time_data = data_df["Time"]

    for substrate in substrates:
        color = colours.pop(0)
        colours.append(color)

        symbol = symbols.pop(0)
        symbols.append(symbol)

        for column in data_df:
            if substrate in column:
                data_to_plot = data_df[column]
                plt.scatter(time_data, data_to_plot,
                            c=color, alpha=alpha, s=size,
                            marker=symbol)

    if plot==True:
        plt.show()


""" --- Sensitivity analysis --- """
def get_concentrations_at_timepoint(model, output, timepoint, substrate):
    """
    Return a np.array of concentrations at the specified timepoint (or closest timepoint)

    Args:
        model (Model): The model object
        output (list): Output from run_all_models
        timepoint (int): Timepoint of interest
        substrate (str): Substrate name of interest

    Returns:
        A np.array containing the concentrations from run_all_models at the specified timepoint
        [c1, c2, c3...]
    """
    closest_timepoint = min(model.time, key=lambda x: abs(x - timepoint))
    index = list(model.time).index(closest_timepoint)

    outputs_for_analysis = []
    for y in output:
        output_at_t = y[index][model.run_model_species_names.index(substrate)]
        outputs_for_analysis.append(output_at_t)

    outputs_for_analysis = np.array(outputs_for_analysis)

    return outputs_for_analysis

def get_time_to_concentration(model, output, concentration, substrate, mode='>='):
    """
    Return a np.array containing the time it takes to reach a certain concentration for all the models run.

    Args:
        model (Model): A model object
        output (list): Output from run_all_models
        concentration (int): The concentration of interest
        substrate (str): The substrate of interest
        mode (str): Either '>=' or '<=' which looks for more_or_equal or less_or_equal respectively.

    Returns:
        A np.array containing the times taken to reach concentration for all models from run_all_models
    """

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
    """
    Run the SALib sobal sensitivity analysis

    Args:
        salib_problem (dict): The salib problem used to make the samples
        output_to_analyse (np.array): A np.array containing the output of interest.
        second_order (bool): Look at second order interactions. Default=Fa;se
        num_resample(int): salib, number of resamples.  Default=100
        conf_level (float): salib confidence level, default = 0.95

    Returns:
        A dataframe containing the output from the sobal sensitivity analysis
    """

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
    """
    Remove any entry with an ST less than specified

    Args:
        dataframe (pandas.Dataframe): dataframe containing sensitivity analysis output
        column (str): Column name, default is 'ST'
        less_than (float): Remove anything less than this

    Returns:
        New dataframe.
    """

    new_df = dataframe[dataframe[column] > less_than]

    return new_df

def plot_sa_total_sensitivity(df):
    """
    Plot the sensitivity analysis

    Args:
        df: Dataframe containing output of sensitivity analysis.
    """
    df.sort_values("ST", inplace=True, ascending=False)

    x_names = df.index.values
    x = np.arange(len(x_names))
    st = df['ST']
    st_err = df['ST_conf']

    plt.bar(x, st, align='center', yerr=st_err, edgecolor='black', color='#000090')
    plt.xticks(x, x_names, rotation=90)
    plt.ylabel("ST")

