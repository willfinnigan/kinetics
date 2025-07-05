from __future__ import annotations
import numpy as np

from kinetics.model_results import ModelResult
from kinetics.solvers.jax_solver import JaxSolver
from kinetics.solvers.scipy_solver import SciPy_Solver

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from kinetics.model_results import ModelResult

class Model(object):

    def __init__(self):

        """ Reactions - a list of reaction classes """
        self._reactions = []
        self.timeseries = None
        self.y = []
        self.ts = None

    def set_time(self, start: int, end: int, steps: int, mode='linear'):
        if mode == 'linear':
            self.ts = np.linspace(start, end, steps)
        elif mode == 'log':
            self.ts = np.logspace(np.log10(start), np.log10(end), steps)
        else:
            raise ValueError(f"Unknown time series mode: {mode}. Use 'linear' or 'log'.")

    def add_reaction(self, reaction):
        self._reactions.append(reaction)

    def _parameters_and_species_from_reactions(self):
        """This loads the default parameters from the model, and sets all required species to 0"""

        species, parameters = {}, {}

        for reaction in self._reactions:
            reaction.set_parameter_defaults_to_mean()

            # Load the model's parameters from the reaction's parameters
            for name in reaction.parameters:
                if name not in parameters:
                    parameters[name] = reaction.parameters[name]
                else:
                    print('Warning - parameter ' + name + ' already set in model, not overwriting with reaction value')

            # Load the model's species from the reaction's species
            for substrate in reaction.substrates + reaction.products + reaction.reaction_substrate_names:
                if substrate not in species:
                    species[substrate] = 0

        return species, parameters

    def _setup_model(self, species_names, parameter_names):
        """ This loads the indexes of species and parameters into each reaction
        so that these can be accessed when the model is run.
        """
        for reaction in self._reactions:
            reaction.setup_reaction(species_names, parameter_names)

    def run_model(self, starting_concentrations, mode='scipy') -> ModelResult:

        # Get the default values for all parameters and species (which is 0)
        species, parameters = self._parameters_and_species_from_reactions()

        # Set the starting concentrations
        species.update(starting_concentrations)

        # Extract species and parameters into ordered lists
        species_names, species_values = zip(*species.items())
        parameter_names, parameter_values = zip(*parameters.items())

        # Set the indexes in the reactions
        self._setup_model(species_names, parameter_names)

        if mode == 'scipy':
            solver = SciPy_Solver(mxsteps=5000)
        elif mode == 'jax':
            solver = JaxSolver()
        else:
            raise ValueError(f"Unknown solver mode: {mode}. Use 'scipy', 'jax', or 'jax_gpu'.")

        y = solver.run(self._reactions, species_names, species_values, parameter_values, self.ts)
        result = ModelResult(self, y, species_names)
        return result






