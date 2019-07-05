"""
The Uncertainty Module is for running uncertainty analysis.
There are 2 mains steps.
1.  Make samples
2.  Run models for all the samples.
"""

from tqdm import tqdm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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

def make_samples_from_distributions(model, num_samples=1000, negative_allowed=[]):
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

def dict_of_samples(samples):
    """
    Gives a dictionary containing {'Sample_name' : [all samples]}
    Args:
        samples (list): The output from make_samples. Each entry in the list is a tuple containing (parameter_dict, species_dict)

    Returns:
        A dictionary containing {'Sample_name' : [all samples]}

    """

    samples_dict = {}

    for s in samples:
        # s = (param_dict, species_dict)
        for name in {**s[0], **s[1]}:

            try:
                values = samples_dict[name]
            except KeyError:
                values = samples_dict[name] = []

            sample_to_add = {**s[0], **s[1]}[name]
            values.append(sample_to_add)

    return samples_dict


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
