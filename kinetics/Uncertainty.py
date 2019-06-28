from tqdm import tqdm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns


def check_not_neg(sample, name, negative_allowed):
    if sample == None:
        return False

    if (sample < 0) and (name not in negative_allowed):
        return False

    return True

def return_ys_for_a_single_substrate(model, output, substrate_name):
    """
    From run_all_models_ouput, return only the outputs for substrate_name

    collected_output = [ [t0, r1, r2, r3],
                         [t1, r1, r2, r3].... ]


    :param model: an instance of the Model class
    :param ouput: the output from the run_all_models function - [ys1, ys2, ys3 .. ]
    :param substrate_name: he name of the substrate to return the ys for
    :return: collected_output = [ [t0, r1, r2, r3],   [t1, r1, r2, r3].... ]
    """

    collected_output = []
    species_names = list(model.species.keys())

    for i in range(len(model.time)):
        timepoint = [model.time[i]]

        for y in output:
            timepoint.append(y[i][species_names.index(substrate_name)])

        collected_output.append(timepoint)

    return collected_output


def make_samples_from_distributions(model, num_samples=1000, negative_allowed=[]):

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


def plot_parameters(samples_dict, parameter_names, units={}, plot=True,
                    colour='red', alpha=0.01, size=10, params_to_log=[], num_graphs=6, figsize=(15,5), wspace=3,
                    xaxis_rotation=0):

    fig = plt.figure(figsize=figsize)
    grid = plt.GridSpec(1, num_graphs, wspace=wspace, figure=fig)

    for i, name in enumerate(parameter_names):
        data = samples_dict[name]
        plt.subplot(grid[0, i])

        if name in params_to_log:
            plt.yscale('log')

        sns.stripplot(x=data,orient='v', color=colour, size=size, jitter=False, alpha=alpha)
        plt.xlim(-0.5, 0.5)
        plt.xticks([0], [name], rotation=xaxis_rotation)

        if name in units:
            unit = units[name]
            plt.ylabel(unit)

    if plot==True:
        plt.show()



def plot_distribution(samples_dict, parameter_distributions, param_name, units='',
                      num_points=100, figsize=[8, 5], colour='darkblue', alpha=0.5,
                      log=False, plot=True):

    max_sample = max(samples_dict[param_name])
    max_sample += max_sample*0.05
    min_sample = min(samples_dict[param_name])
    min_sample -= max_sample*0.05

    x = np.linspace(min_sample, max_sample, num_points)

    distribution = parameter_distributions[param_name]
    y = distribution.pdf(x)

    fig, ax = plt.subplots(figsize=figsize)
    plot1 = ax.plot(x, y, color='black')
    plot2 = plt.fill_between(x, y, color=colour, alpha=alpha, linewidth=0)

    # adds a title and axes labels
    ax.set_title(param_name)
    ax.set_xlabel(units)

    # removing top and right borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.yticks([])

    if log == True:
        plt.xscale('log')

    plt.show()

    if plot == True:
        plt.show()




def plot_single_parameter_distribution_and_sample(samples_dict, parameter_distributions, parameter_name,
                                                  units={}, num_points=1000, figsize=[10, 5],
                                                  colour='darkblue', alpha_dist=0.5, alpha_sample=0.01,
                                                  size=10, log=False, plot=False,
                                                  xaxis_rotation=0, width_ratio=[4,1]):

    fig = plt.figure(1)
    grid = gridspec.GridSpec(1, 2, width_ratios=width_ratio)

    max_sample = max(samples_dict[parameter_name])
    max_sample += max_sample*0.05
    min_sample = min(samples_dict[parameter_name])
    min_sample -= max_sample*0.05
    x = np.linspace(min_sample, max_sample, num_points)
    distribution = parameter_distributions[parameter_name]
    y = distribution.pdf(x)

    ax0 = plt.subplot(grid[0])
    ax0.plot(x, y, color='black')
    ax0.fill_between(x, y, color=colour, alpha=alpha_dist, linewidth=0)

    ax0.set_title(parameter_name)

    if parameter_name in units:
        unit = units[parameter_name]
        ax0.set_xlabel(unit)

    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.spines['left'].set_visible(False)
    plt.yticks([])

    ax1 = plt.subplot(grid[1])
    data = samples_dict[parameter_name]
    y = np.random.normal(1, 0.04, size=len(data))
    sns.stripplot(x=data, orient='v', color=colour, size=size, jitter=False, alpha=alpha_sample)
    #ax1.scatter(data, y, color=colour, size=size, alpha=alpha)
    plt.xlim(-0.5, 0.5)
    plt.xticks([0], [parameter_name], rotation=xaxis_rotation)

    if parameter_name in units:
        unit = units[parameter_name]
        plt.ylabel(unit)

    if log == True:
        plt.yscale('log')

    if plot == True:
        plt.show()

def set_units(parameters, kcat='$min^{-1}$', km_ki='$\mu M$', concs=False):

    units = {}
    for name in parameters:
        if 'kcat' in name:
            units[name] = kcat
        elif ('km' in name) or ('ki' in name):
            units[name] = km_ki
        elif concs==True:
            units[name] = km_ki

    return units

def plot_sampling(samples, model, param_logs, folder_name='', units={}, plot=True, save=False, dpi=100):
    samples_dict = dict_of_samples(samples)
    units_dict = set_units(model.parameter_distributions)
    units_dict.update(units)

    for name in model.parameter_distributions:

        log = False
        if name in param_logs:
            log = True

        plot_single_parameter_distribution_and_sample(samples_dict, model.parameter_distributions, name, units=units_dict, log=log)

        if save == True:
            plt.savefig(str(folder_name) + '/' + str(name) + '.png', dpi=dpi)

        if plot == True:
            plt.show()