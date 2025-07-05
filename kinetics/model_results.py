from __future__ import annotations
import numpy as np
import pandas as pd

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from kinetics import Model


class ModelResult(object):

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
    def results_dataframe(self):
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