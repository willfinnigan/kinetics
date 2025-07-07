import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from SALib.analyze import sobol

from kinetics.models.result_classes.multi_result import MultiModelResult


# =============================================================================
# PUBLIC API
# =============================================================================

class SensitivityResult:
    """Simple class to hold sensitivity analysis results and provide plotting"""
    
    def __init__(self, data: pd.DataFrame):
        self.data = data
    
    def plot(self):
        """Plot the sensitivity analysis results"""
        plot_sa_total_sensitivity(self.data)


def analyze_sensitivity_at_timepoint(result: MultiModelResult,
                                   sampler,
                                   species: str,
                                   timepoint: float,
                                   threshold: float = 0.01,
                                   second_order: bool = False,
                                   num_resample: int = 100,
                                   conf_level: float = 0.95) -> SensitivityResult:
    """
    Analyze sensitivity for concentration at a specific timepoint.
    
    Args:
        result: MultiModelResult from model run
        sampler: The sampler used (must have .problem attribute)
        species: Species name to analyze
        timepoint: Time point of interest
        threshold: Filter out sensitivity indices below this value
        second_order: Calculate second order interactions
        num_resample: Number of resamples for confidence intervals
        conf_level: Confidence level for intervals
    
    Returns:
        SensitivityResult object with data and plot method
    """
    # Extract concentrations at timepoint
    concentrations = _get_concentrations_at_timepoint(result, timepoint, species)
    
    # Run Sobol analysis
    df = _run_sobol_analysis(sampler.problem, concentrations, second_order, num_resample, conf_level)
    
    # Filter out low sensitivity parameters
    df_filtered = df[df['ST'] > threshold]
    
    return SensitivityResult(df_filtered)


def analyze_sensitivity_time_to_concentration(result: MultiModelResult,
                                            sampler,
                                            species: str,
                                            concentration: float,
                                            mode: str = '>=',
                                            threshold: float = 0.01,
                                            second_order: bool = False,
                                            num_resample: int = 100,
                                            conf_level: float = 0.95) -> SensitivityResult:
    """
    Analyze sensitivity for time to reach a concentration.
    
    Args:
        result: MultiModelResult from model run
        sampler: The sampler used (must have .problem attribute)
        species: Species name to analyze
        concentration: Concentration threshold
        mode: Either '>=' or '<=' for reaching concentration
        threshold: Filter out sensitivity indices below this value
        second_order: Calculate second order interactions
        num_resample: Number of resamples for confidence intervals
        conf_level: Confidence level for intervals
    
    Returns:
        SensitivityResult object with data and plot method
    """
    # Extract time to concentration
    times = _get_time_to_concentration(result, concentration, species, mode)
    
    # Run Sobol analysis
    df = _run_sobol_analysis(sampler.problem, times, second_order, num_resample, conf_level)
    
    # Filter out low sensitivity parameters
    df_filtered = df[df['ST'] > threshold]
    
    return SensitivityResult(df_filtered)


# =============================================================================
# PRIVATE HELPER FUNCTIONS
# =============================================================================

def _get_concentrations_at_timepoint(result: MultiModelResult,
                                    timepoint: float,
                                    substrate: str):
    """
    Return a np.array of concentrations at the specified timepoint (or closest timepoint)

    Args:
        result (MultiModelResult): The result of a run
        timepoint (float): Timepoint of interest
        substrate (str): Substrate name of interest

    Returns:
        A np.array containing the concentrations from run_all_models at the specified timepoint
        [c1, c2, c3...]
    """
    closest_timepoint = min(result.ts, key=lambda x: abs(x - timepoint))
    index = list(result.ts).index(closest_timepoint)

    outputs_for_analysis = []
    for y in result.multi_ys:
        output_at_t = y[index][result.species_names.index(substrate)]
        outputs_for_analysis.append(output_at_t)

    outputs_for_analysis = np.array(outputs_for_analysis)

    return outputs_for_analysis


def _get_time_to_concentration(result: MultiModelResult,
                              concentration: float,
                              substrate: str,
                              mode='>='):
    """
    Return a np.array containing the time it takes to reach a certain concentration for all the models run.

    Args:
        result (MultiModelResult): The result of a run
        concentration (int): The concentration of interest
        substrate (str): The substrate of interest
        mode (str): Either '>=' or '<=' which looks for more_or_equal or less_or_equal respectively.

    Returns:
        A np.array containing the times taken to reach concentration for all models from run_all_models
    """

    list_of_times = []

    for y in result.multi_ys:
        y = np.transpose(y)

        substrate_index = result.species_names.index(substrate)
        y_substrate = y[substrate_index]

        if mode == '<=':
            index_for_conc = np.where(y_substrate <= concentration)
        elif mode == '>=':
            index_for_conc = np.where(y_substrate >= concentration)

        if len(index_for_conc[0]) == 0:
            time = result.ts[-1]
        else:
            index_for_conc = index_for_conc[0][0]
            time = result.ts[index_for_conc]

        list_of_times.append(time)

    list_of_times = np.array(list_of_times)
    return list_of_times


def _run_sobol_analysis(salib_problem,
                       output_to_analyse,
                       second_order=False, 
                       num_resample=100,
                       conf_level=0.95):
    """
    Run the SALib sobol sensitivity analysis

    Args:
        salib_problem (dict): The salib problem used to make the samples
        output_to_analyse (np.array): A np.array containing the output of interest.
        second_order (bool): Look at second order interactions. Default=False
        num_resample(int): salib, number of resamples.  Default=100
        conf_level (float): salib confidence level, default = 0.95

    Returns:
        A dataframe containing the output from the sobol sensitivity analysis
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


# =============================================================================
# PLOTTING UTILITIES
# =============================================================================

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