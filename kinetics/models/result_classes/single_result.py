from __future__ import annotations
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from kinetics.models.model_class import Model

def plot_data(substrates: list[str], data_df: pd.DataFrame,
              alpha: float = 0.5, size: int = 35, 
              colours: list[str] = ['black'], 
              symbols: list[str] = ["o", "s", '^', 'v']) -> None:
    """Add experimental data points to current matplotlib plot.
    
    Args:
        substrates: List of substrate names to plot
        data_df: DataFrame containing experimental data with 'Time' column
        alpha: Transparency for data points (0-1)
        size: Size of data points
        colours: List of colors for different substrates
        symbols: List of symbols for different substrates
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

class SingleModelResult:
    """Results container for single model run with fixed parameters.
    
    Stores simulation results and provides methods for data export and visualization.
    
    Attributes:
        model: The model that was run
        y: Solution array with shape (n_timepoints, n_species)
        ts: Time points array
        species_names: Ordered list of species names
    """
    
    def __init__(self, model: Model, y: np.ndarray, species_names: list[str]):
        """Initialize result container.
        
        Args:
            model: The model that was run
            y: Solution array from ODE solver
            species_names: Ordered list of species names
        """
        self.model = model
        self.y = y
        self.ts = model.ts
        self.species_names = species_names

    def dataframe(self) -> pd.DataFrame:
        """Export results as a pandas DataFrame.
        
        Returns:
            DataFrame with time column and one column per species
        """
        ys_at_t = {'Time': self.ts}

        for i in range(len(self.species_names)):
            name = self.species_names[i]
            ys_at_t[name] = []

            for t in range(len(self.ts)):
                ys_at_t[name].append(self.y[t][i])

        df = pd.DataFrame(ys_at_t)

        return df

    def plot(self, substrate: str, units: list[str] = ['', '']) -> None:
        """Plot substrate concentration vs time.
        
        Args:
            substrate: Name of substrate to plot
            units: List of [y_label, x_label] for axis labels
        """

        ys_at_t = []
        i = self.species_names.index(substrate)
        for t in range(len(self.ts)):
            ys_at_t.append(self.y[t][i])

        plt.plot(self.ts, ys_at_t, label=substrate)
        plt.ylabel(units[0])
        plt.xlabel(units[1])
        plt.legend()

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


