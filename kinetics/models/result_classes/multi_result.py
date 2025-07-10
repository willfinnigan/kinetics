from __future__ import annotations

from functools import cache
from typing import List, Dict
import numpy as np
import pandas as pd

from typing import TYPE_CHECKING

from matplotlib import pyplot as plt

from kinetics.models.result_classes.single_result import plot_data

if TYPE_CHECKING:
    from kinetics.models.model_class import Model


class MultiModelResult:
    """Results container for multiple model runs with parameter uncertainty.
    
    Stores ensemble simulation results and provides methods for statistical
    analysis and uncertainty visualization.
    
    Attributes:
        model: The model that was run
        multi_ys: List of solution arrays from each run
        ts: Time points array
        species_names: Ordered list of species names
        color_palette: Default colors for plotting
        substrate_colors: Color assignments for each species
    """

    def __init__(self, model: Model, multi_ys: List[np.ndarray], species_names: list[str]):
        """Initialize multi-run result container.
        
        Args:
            model: The model that was run
            multi_ys: List of solution arrays from each run
            species_names: Ordered list of species names
        """
        self.model = model
        self.multi_ys = multi_ys
        self.ts = model.ts
        self.species_names = species_names
        self.color_palette = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
        self.substrate_colors = {}
        
    def _get_substrate_color(self, substrate_name: str) -> str:
        """Get or assign a color for a substrate.
        
        Args:
            substrate_name: Name of the substrate
            
        Returns:
            Color string for the substrate
        """
        if substrate_name not in self.substrate_colors:
            color_index = len(self.substrate_colors) % len(self.color_palette)
            self.substrate_colors[substrate_name] = self.color_palette[color_index]
        return self.substrate_colors[substrate_name]

    def _return_ys_for_a_single_substrate(self, substrate_name: str) -> list:
        """Extract time series data for a single substrate across all runs.
        
        Args:
            substrate_name: Name of the substrate
            
        Returns:
            List of [time, run1_value, run2_value, ...] for each timepoint
        """

        collected_output = []
        for i in range(len(self.ts)):
            timepoint = [self.ts[i]]

            for y in self.multi_ys:
                timepoint.append(y[i][self.species_names.index(substrate_name)])

            collected_output.append(timepoint)

        return collected_output

    @cache
    def dataframe(self) -> Dict[str, pd.DataFrame]:
        """Export all run results as DataFrames.
        
        Returns:
            Dictionary mapping species names to DataFrames containing
            time column and one column per simulation run
        """

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
    def dataframe_quartiles(self, quartile: int = 95) -> Dict[str, dict]:
        """Calculate confidence intervals for all species.
        
        Args:
            quartile: Percentile for confidence intervals (default 95%)
            
        Returns:
            Dictionary mapping species names to dictionaries with
            'Time', 'High', 'Low', 'Mean' keys containing arrays
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

    def plot_all(self, substrate: str, units: list[str] = ['', ''],
                 colour: str = None, alpha: float = 0.1, 
                 linewidth: float = 0.1, y_min: bool = True) -> None:
        """Plot every model run for a single substrate.
        
        Args:
            substrate: Substrate name to plot
            units: List of [y_label, x_label] for axis labels
            colour: Color for lines (auto-assigned if None)
            alpha: Line transparency
            linewidth: Width of lines
            y_min: Whether to auto-set y-axis minimum
        """

        dataframes = self.dataframe()

        ylabel = units[0]
        xlabel = units[1]
        
        if colour is None:
            colour = self._get_substrate_color(substrate)

        df = dataframes[substrate]
        for i in range(1, len(df.columns)):
            plt.plot(df['Time'], df[str(i)],
                     color=colour, alpha=alpha, linewidth=linewidth)
        print(str(substrate) + ' - ' + str(colour))
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        if y_min != True:
            plt.ylim(bottom=y_min)

    def plot(self, substrate: str, quartile: int = 95,
             colour: str = None, alpha: float = 0.1, 
             units: list[str] = ['', '']) -> None:
        """Plot confidence intervals for a substrate.
        
        Args:
            substrate: Substrate name to plot
            quartile: Percentile for confidence intervals
            colour: Color for plot (auto-assigned if None)
            alpha: Fill transparency
            units: List of [y_label, x_label] for axis labels
        """

        dataframes = self.dataframe_quartiles(quartile=quartile)

        df = dataframes[substrate]
        
        if colour is None:
            colour = self._get_substrate_color(substrate)

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

    def plot_data(self, substrates: list[str], data_df: pd.DataFrame,
                  alpha: float = 0.5, size: int = 35, 
                  colours: list[str] = ['black'], 
                  symbols: list[str] = ["o", "s", '^', 'v']) -> None:
        """Add experimental data points to current plot.
        
        Args:
            substrates: List of substrate names to plot
            data_df: DataFrame containing experimental data
            alpha: Transparency for data points
            size: Size of data points
            colours: List of colors for different substrates
            symbols: List of symbols for different substrates
        """
        return plot_data(substrates, data_df,
                        alpha=alpha, size=size, colours=colours, symbols=symbols)