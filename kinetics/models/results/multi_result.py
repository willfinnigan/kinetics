from __future__ import annotations

from functools import cache
from typing import List, Dict
import numpy as np
import pandas as pd

from typing import TYPE_CHECKING

from matplotlib import pyplot as plt

if TYPE_CHECKING:
    from kinetics.models.model_class import Model


class MultiModelResult(object):

    def __init__(self, model: Model, multi_ys: List[np.ndarray], species_names):
        """
        This class is used to store the results of multiple model runs where samples have been used.

        Args:
            model (Model): the model that was run
            ys (List[np.ndarray]): the outputs of the model runs
            species_names (List[str]): names of the species
        """
        self.model = model
        self.multi_ys = multi_ys
        self.ts = model.ts
        self.species_names = species_names

    def _return_ys_for_a_single_substrate(self, substrate_name):

        collected_output = []
        for i in range(len(self.ts)):
            timepoint = [self.ts[i]]

            for y in self.multi_ys:
                timepoint.append(y[i][self.species_names.index(substrate_name)])

            collected_output.append(timepoint)

        return collected_output

    @cache
    def dataframe(self) -> Dict[str, pd.DataFrame]:
        """Returns a dictionary of {'species_name': dataframe} for each species in the model."""

        all_runs_substrate_dataframes = {}

        for name in self.species_names:
            # format: [[t0, r1, r2, r3], [t1, r1, r2, r3]..]
            collected_runs = self._return_ys_for_a_single_substrate(name)
            all_runs = {"Time": []}
            column_titles = ['Time']

            for i in range(1, len(collected_runs[0])):
                all_runs[str(i)] = []
                column_titles.append(str(i))

            for timepoint in collected_runs:
                all_runs['Time'].append(timepoint[0])

                for i in range(1, len(timepoint)):
                    all_runs[str(i)].append(timepoint[i])

            df = pd.DataFrame(all_runs, columns=column_titles)
            all_runs_substrate_dataframes[name] = df

        return all_runs_substrate_dataframes

    @cache
    def dataframe_quartiles(self, quartile=95):
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

        for name in self.species_names:

            quartiles = {"Time": [], "High": [], "Low": [], "Mean": []}
            ys_for_single_substrate = self._return_ys_for_a_single_substrate(name)

            for i in range(len(ys_for_single_substrate)):
                # output_at_t will be a array.  i=0 is time, after than the substrate values at that time.
                output_at_t = ys_for_single_substrate[i]

                quartiles["Time"].append(output_at_t[0])
                quartiles["High"].append(np.percentile(output_at_t[1:], quartile))
                quartiles["Low"].append(np.percentile(output_at_t[1:], 100 - quartile))
                quartiles["Mean"].append(np.mean(output_at_t[1:]))

            dataframes[name] = quartiles

        return dataframes

    def plot(self, substrate, units=['', ''],
             colour='blue', alpha=0.1, linewidth=0.1, y_min=True):
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

        """

        dataframes = self.dataframe()

        ylabel = units[0]
        xlabel = units[1]

        df = dataframes[substrate]
        for i in range(1, len(df.columns)):
            plt.plot(df['Time'], df[str(i)],
                     color=colour, alpha=alpha, linewidth=linewidth)
        print(str(substrate) + ' - ' + str(colour))
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        if y_min != True:
            plt.ylim(bottom=y_min)

    def plot_ci(self,
                substrate, ci=95,
                colour='darkblue', alpha=0.1, units=['', '']):
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

        dataframes = self.dataframe_quartiles(ci=ci)

        df = dataframes[substrate]

        time = df['Time']
        high = df['High']
        low = df['Low']
        mean = df['Mean']

        # high_line = plt.plot(time, high, color=color, linestyle="--", linewidth = 0.5)
        # low_line = plt.plot(time, low,  color=color, linestyle="--", linewidth = 0.5)
        plt.plot(time, mean, color=colour, linewidth=0.5, label=substrate)
        plt.fill_between(time, high, y2=low, color=colour, alpha=alpha, linewidth=0)

        plt.ylabel(units[0])
        plt.xlabel(units[1])
