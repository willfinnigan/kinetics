from tqdm import tqdm
import pandas as pd
import numpy as np

def run_all_models(model, samples, logging=True):
    """
    Run all the models for a set of samples.

    Args:
        model (kinetics.model_module): A model object
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

def return_ys_for_a_single_substrate(model, output, substrate_name):

    collected_output = []
    species_names = list(model.species.keys())

    for i in range(len(model.time)):
        timepoint = [model.time[i]]

        for y in output:
            timepoint.append(y[i][species_names.index(substrate_name)])

        collected_output.append(timepoint)

    return collected_output

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

def dataframes_quartiles(model, output, substrates=[], quartile=95, logging=False):
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

        if logging == True:
            print(name)

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