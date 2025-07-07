from __future__ import annotations
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from kinetics.models.model_class import Model


class SingleModelResult(object):
    def __init__(self, model: Model, y: np.ndarray, species_names):
        """
        This class is used to store the results of a model run.

        Args:
            model (Model): the model that was run
            y (np.ndarray): the output of the model run
        """
        self.model = model
        self.y = y
        self.ts = model.ts
        self.species_names = species_names

    # Export results as dataframe and plot
    def dataframe(self):
        """
        Gives the results of a model run as a dataframe

        Returns:
            Pandas dataframe of results
        """
        ys_at_t = {'Time': self.ts}

        for i in range(len(self.species_names)):
            name = self.species_names[i]
            ys_at_t[name] = []

            for t in range(len(self.ts)):
                ys_at_t[name].append(self.y[t][i])

        df = pd.DataFrame(ys_at_t)

        return df

    def plot(self, substrate, units=['', '']):
        """
        Plot a graph of substrate concentration vs time.

        Need to call plt.show() afterwards

        Args:
            substrate (str): Name of substrate to plot
            plot (bool): Default False.  If True calls plt.show()
        """

        ys_at_t = []
        i = self.species_names.index(substrate)
        for t in range(len(self.ts)):
            ys_at_t.append(self.y[t][i])

        plt.plot(self.ts, ys_at_t, label=substrate)
        plt.ylabel(units[0])
        plt.xlabel(units[1])
        plt.legend()