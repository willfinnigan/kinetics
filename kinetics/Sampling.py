import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns


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

    distribution = parameter_distributions[parameter_name]

    max_sample = distribution.ppf(0.95)
    max_sample += max_sample*0.05
    min_sample = distribution.ppf(0.05)
    min_sample -= max_sample*0.05
    x = np.linspace(min_sample, max_sample, num_points)

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