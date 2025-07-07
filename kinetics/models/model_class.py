from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from tqdm import tqdm

from kinetics.sampling.scipy_sampling import ScipyDist_Sampler
from kinetics.solvers.scipy_solver import SciPySolver

from kinetics.models.result_classes.single_result import SingleModelResult
from kinetics.models.result_classes.multi_result import MultiModelResult

if TYPE_CHECKING:
    from kinetics.sampling.sampling_interface import Sampler
    from kinetics.solvers.solver_interface import ODESolver

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
                    raise Exception(f"Parameter {name} already set in model, can't overwrite")

            # Load the model's species from the reaction's species
            for substrate in reaction.substrates + reaction.products + reaction.reaction_substrate_names:
                if substrate not in species:
                    species[substrate] = 0

        return species, parameters

    def get_parameter_distributions(self):
        # Get parameter distributions from the reactions - update parameters with these
        parameter_distributions = {}
        for reaction in self._reactions:
            # check we're not overwriting existing parameters
            for name in reaction.parameter_distributions:
                if name in parameter_distributions:
                    raise Exception(f"Parameter {name} already set in model, can't overwrite")

            # update the model's parameter distributions with the reaction's
            parameter_distributions.update(reaction.parameter_distributions)

        return parameter_distributions

    def _setup_model(self, species_names, parameter_names):
        """This loads the indexes of species and parameters into each reaction
        so that these can be accessed when the model is run.
        """
        for reaction in self._reactions:
            reaction.setup_reaction(species_names, parameter_names)

    def _set_default_species(self, species):
        """Any species which are distributions, convert to the mean"""
        for name, value in species.items():
            if hasattr(value, 'rvs'):
                species[name] = value.mean()
            if isinstance(value, tuple) and len(value) == 2:
                # If its a tuple, take the value directly in the middle of the two values
                species[name] = (value[0] + value[1]) / 2
        return species

    def run_single(self,
                   starting_concentrations: dict,
                   solver: ODESolver = SciPySolver()) -> SingleModelResult:

        # Get the default values for all parameters and species (which is 0)
        species, parameters = self._parameters_and_species_from_reactions()


        # Set the starting concentrations
        starting_concentrations = self._set_default_species(starting_concentrations)
        species.update(starting_concentrations)

        # Extract species and parameters into ordered lists
        species_names, species_values = zip(*species.items())
        parameter_names, parameter_values = zip(*parameters.items())

        # Set the indexes in the reactions
        self._setup_model(species_names, parameter_names)

        y = solver.run(self._reactions, species_names, species_values, parameter_values, self.ts)
        result = SingleModelResult(self, y, species_names)
        return result

    def run_multi(self,
                  starting_concentrations: dict,
                  sampler: Sampler = ScipyDist_Sampler(num_samples=1000),
                  solver: ODESolver = SciPySolver()) -> MultiModelResult:
        """
        Run the model sampling parameter distributions
        """

        # Get the default values for all parameters and species (which is 0)
        species, parameters = self._parameters_and_species_from_reactions()

        # Set the indexes in the reactions
        self._setup_model(list(species.keys()), list(parameters.keys()))

        # Set the starting concentrations - can contain distributions or fixed values
        species.update(starting_concentrations)

        # Get the parameter distributions from the reactions
        parameter_distributions = self.get_parameter_distributions()

        # Now we sample
        samples = sampler.sample(parameter_distributions, species)

        # Ok here is where I think we can vectorize
        # But for now we're just run one at a time (which is the existing implementation)
        results = []
        for parameter_dict, species_dict in tqdm(samples):
            # Update the species and parameters with the sampled values
            # Note - this is important to maintain the correct ordering
            species.update(species_dict)
            parameters.update(parameter_dict)

            # Extract species and parameters into ordered lists
            species_names, species_values = zip(*species.items())
            parameter_names, parameter_values = zip(*parameters.items())

            y = solver.run(self._reactions, species_names, species_values, parameter_values, self.ts)
            results.append(y)

        return MultiModelResult(self, results, list(species.keys()))


    # TODO
    ## Multi output
    ## Look into vectorisation and gpu acceleration
















