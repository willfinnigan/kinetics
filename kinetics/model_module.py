from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import integrate
import matplotlib.pyplot as plt

from kinetics.solvers.jax_solver import JaxSolver
from kinetics.solvers.scipy_solver import SciPy_Solver




class TimeSeries(object):

    def __init__(self, start: int, end: int, steps: int, mxsteps=10000, mode='linear'):
        """
        This function sets the time parameters for the model.  This is how long the model will simulate

        Args:
            start (int): the start time - usually 0
            end (int): the end time (default is 100)
            steps (int): the number of timepoints for the output
        """

        self.start = start
        self.end = end
        self.steps = steps
        if mode == 'linear':
            self.t = np.linspace(self.start, self.end, self.steps)
        elif mode == 'log':
            self.t = np.logspace(np.log10(self.start), np.log10(self.end), self.steps)
        self.mxsteps = mxsteps


class ModelResult(object):

    def __init__(self, model: Model, y: np.ndarray):
        """
        This class is used to store the results of a model run.

        Args:
            model (Model): the model that was run
            y (np.ndarray): the output of the model run
        """
        self.model = model
        self.y = y
        self.t = model.timeseries.t
        self.species_names = model.species_names()

    # Export results as dataframe and plot
    def results_dataframe(self):
        """
        Gives the results of a model run as a dataframe

        Returns:
            Pandas dataframe of results
        """
        ys_at_t = {'Time': self.t}

        for i in range(len(self.species_names)):
            name = self.species_names[i]
            ys_at_t[name] = []

            for t in range(len(self.t)):
                ys_at_t[name].append(self.y[t][i])

        df = pd.DataFrame(ys_at_t)

        return df

class Model(object):

    def __init__(self):

        """ Reactions - a list of reaction classes """
        self._reactions = []
        self._species = {}
        self._parameters = {}
        self.timeseries: TimeSeries = TimeSeries(0, 100, 100)
        self.y = []

    def set_time(self, start: int, end: int, steps: int, mxsteps=5000, mode='linear'):
        self.timeseries = TimeSeries(start, end, steps, mxsteps=mxsteps, mode=mode)

    def set_species(self, species):
        self._species.update(species)

    def species_names(self):
        return list(self._species.keys())

    def add_reaction(self, reaction):
        self._reactions.append(reaction)
        reaction.set_parameter_defaults_to_mean()

        # Load the model's parameters from the reaction's parameters
        for name in reaction.parameters:
            if name not in self._parameters:
                self._parameters[name] = reaction.parameters[name]
            else:
                print('Warning - parameter ' + name + ' already set in model, not overwriting with reaction value')

        # Load the model's species from the reaction's species
        for substrate in reaction.substrates + reaction.products + reaction.reaction_substrate_names:
            if substrate not in self._species:
                self._species[substrate] = 0

    def setup_model(self):
        species_names = self.species_names()
        for reaction in self._reactions:
            reaction.get_indexes(species_names)
            reaction.get_parameters(self._parameters)
            for modifier in reaction.modifiers:
                modifier.get_substrate_indexes(reaction.reaction_substrate_names)
                modifier.get_parameter_indexes(reaction.parameter_names)

    def run_model(self, mode='scipy') -> ModelResult:
        self.setup_model()

        if mode == 'scipy':
            solver = SciPy_Solver(mxsteps=self.timeseries.mxsteps)
        elif mode == 'jax':
            solver = JaxSolver(use_gpu=False)
        elif mode == 'jax_gpu':
            solver = JaxSolver(use_gpu=True)
        else:
            raise ValueError(f"Unknown solver mode: {mode}. Use 'scipy', 'jax', or 'jax_gpu'.")

        y = solver.run(self._reactions, self._species, self._parameters, self.timeseries.t)
        result = ModelResult(self, y)
        return result






