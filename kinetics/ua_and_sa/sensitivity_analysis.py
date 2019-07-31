import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from SALib.analyze import sobol



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
