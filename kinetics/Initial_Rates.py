import pandas as pd
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib as mpl


def calc_initial_rates_ua(model, substrate_name, enzyme_name, substrate_concs,
                       time=1, ua_samples=500, ua_quartile_range=95,
                       logging=False):

    model.set_time(0, time, 2)
    ua = UA(model, num_samples=ua_samples, quartile_range=ua_quartile_range, logging=logging)

    rate_quartiles = pd.DataFrame(columns=["Substrate", "High", "Low", "Mean"])
    rate_all = pd.DataFrame()

    print("Running Initial Rates Uncertainty Analysis")
    for conc in substrate_concs:
        print(str(conc) + ', ', end ="")

        ua.model.reaction_species.update({substrate_name: (conc, 0)})
        ua.model.setup_model()
        ua.run_standard_ua()
        timecourse = ua.substrate_dataframes[substrate_name]

        start = timecourse.iloc[0]
        end = timecourse.iloc[-1]
        difference = start - end
        difference = difference / time  # This will give rates in uM / min
        difference = difference / model.species[enzyme_name]  # This will give rates in uM / min / uM_enz

        rate_quartiles = rate_quartiles.append({'Substrate': conc,
                                                 'Low': difference['High'],
                                                 'High': difference['Low'],
                                                 'Mean': difference['Mean']},
                                                 ignore_index=True)

        timecourse_all = ua.calculate_all_runs_substrate_dataframes()[substrate_name]
        start = timecourse_all.iloc[0]
        end = timecourse_all.iloc[-1]
        difference = start - end
        difference = difference / time  # This will give rates in uM / min
        difference = difference / model.species[enzyme_name]  # This will give rates in uM / min / uM_enz

        rate_all = rate_all.append(difference, ignore_index=True)

    rate_all.drop(columns='Time', inplace=True)
    rate_all.insert(0, 'Substrate', substrate_concs)

    print()

    return rate_quartiles, rate_all


def calc_initial_rates_single(model, substrate_name, enzyme_name, substrate_concs,
                              time=1, verbose=False):

    model.set_time(0, time, 2)
    rates = []

    if verbose==True:
        print('Modelling intial rate experiments with substrate ' + str(substrate_name) + ' at concentrations:')

    for conc in substrate_concs:

        if verbose == True:
            print(str(conc) + ', ', end ="")

        model.reaction_species.update({substrate_name: (conc, 0)})
        model.setup_model()
        model.run_model()
        timecourse = model.results_dataframe()

        start = timecourse.iloc[0][substrate_name]
        end = timecourse.iloc[-1][substrate_name]
        difference = start - end
        difference = difference / time  # This will give rates in uM / min
        difference = difference / model.species[enzyme_name]  # This will give rates in uM / min / uM_enz
        rates.append(difference)

    if verbose == True:
        print()

    return rates

def kcat_to_umolminmg(uM_min_uM_enz, mw_enzyme, volume_ml):

    umol_min_uM_enz = uM_min_uM_enz * (volume_ml/1000)
    umol_enz = 1 * (volume_ml/1000)
    mg_enz = (umol_enz * mw_enzyme)/1000
    umol_min_mg = umol_min_uM_enz / mg_enz

    return umol_min_mg

def concentrations_around_km(reaction, km_param_name,
                             include_zero=True,
                             datapoints=(0, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)):

    km = reaction.parameter_defaults[km_param_name]
    substrate_concs = []

    for point in datapoints:
        substrate_concs.append(point*km)

    if include_zero==False:
        substrate_concs=substrate_concs[1:-1]

    return substrate_concs

def standard_mm_equation(x, km, vmax):
    return vmax * (x / (km + x))

def fit_mm(x_data, y_data, func=standard_mm_equation, param_names=['Km', 'Kcat'], verbose=False):
    popt, pcov = scipy.optimize.curve_fit(func, x_data, y_data)
    perr = np.sqrt(np.diag(pcov))

    parameters = []
    for i in range(len(popt)):
        parameters.append([popt[i], perr[i]])

    x_fit = np.linspace(0,x_data[-1],200)
    y_fit = func(x_fit, *popt)

    to_return = {'x_fit' : x_fit,
                 'y_fit' : y_fit}

    for i in range(len(param_names)):
        name = param_names[i]
        param = round(parameters[i][0], 2)
        error = round(parameters[i][1], 2)
        to_return[name] = (param, error)

        if verbose == True:
            print(str(name) + ' = ' + str(param) + ' ± ' + str(error))

    return to_return


def plot_scatter_all_runs(rates, colour='black', alpha=0.5, size=4):

    concs = rates[1]['Substrate']
    for column in rates[1].columns[1:]:
        y_data = []
        for i in range(len(rates[1][column])):
            v = rates[1][column][i]
            y_data.append(v)
        plt.scatter(concs, y_data, s=size, c=colour, alpha=alpha)

def plot_fit_ua(rates, concs, colour='blue', linewidth=1.5, outer_linewidth=0.5, outer_line_style='--'):
    fit_mean = fit_mm(concs, rates[0]['Mean'])
    fit_high = fit_mm(concs, rates[0]['High'])
    fit_low = fit_mm(concs, rates[0]['Low'])

    plt.plot(fit_high['x_fit'], fit_high['y_fit'], color=colour, linestyle=outer_line_style, linewidth=outer_linewidth)
    plt.plot(fit_low['x_fit'], fit_low['y_fit'], color=colour, linestyle=outer_line_style, linewidth=outer_linewidth)
    plt.plot(fit_mean['x_fit'], fit_mean['y_fit'], color=colour, linewidth=linewidth)