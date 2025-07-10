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
    """Main class for modeling enzyme reactions using ordinary differential equations.
    
    The Model class orchestrates reactions, manages time series, and coordinates
    the execution of kinetic simulations. It supports both single parameter values
    and probability distributions for uncertainty analysis.
    """

    def __init__(self):
        """Initialize empty model with no reactions."""
        self._reactions = []
        self.timeseries = None
        self.y = []
        self.ts = None

    def set_time(self, start: int, end: int, steps: int, mode: str = 'linear') -> None:
        """Set the time points for integration.
        
        Args:
            start: Start time
            end: End time
            steps: Number of time steps
            mode: Time series mode ('linear' or 'log')
            
        Raises:
            ValueError: If mode is not 'linear' or 'log'
        """
        if mode == 'linear':
            self.ts = np.linspace(start, end, steps)
        elif mode == 'log':
            self.ts = np.logspace(np.log10(start), np.log10(end), steps)
        else:
            raise ValueError(f"Unknown time series mode: {mode}. Use 'linear' or 'log'.")

    def add_reaction(self, reaction: 'Reaction') -> None:
        """Add a reaction to the model.
        
        Args:
            reaction: Reaction object to add
        """
        self._reactions.append(reaction)

    def _parameters_and_species_from_reactions(self) -> tuple[dict, dict]:
        """Extract parameters and species from all reactions.
        
        Loads default parameters from reactions and initializes all species to 0.
        
        Returns:
            Tuple of (species_dict, parameters_dict)
        """

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

    def get_parameter_distributions(self) -> dict:
        """Get parameter distributions from all reactions.
        
        Returns:
            Dictionary mapping parameter names to distribution objects
        """
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

    def _setup_model(self, species_names: list[str], parameter_names: list[str]) -> None:
        """Set up reaction indexes for efficient access during simulation.
        
        Loads species and parameter indexes into each reaction for fast lookup
        during ODE solving.
        
        Args:
            species_names: Ordered list of species names
            parameter_names: Ordered list of parameter names
        """
        for reaction in self._reactions:
            reaction.setup_reaction(species_names, parameter_names)

    def _set_default_species(self, species: dict) -> dict:
        """Convert species distributions to mean values.
        
        Args:
            species: Dictionary of species with values or distributions
            
        Returns:
            Dictionary with distributions converted to mean values
        """
        for name, value in species.items():
            if hasattr(value, 'rvs'):
                species[name] = value.mean()
            if isinstance(value, tuple) and len(value) == 2:
                # If its a tuple, take the value directly in the middle of the two values
                species[name] = (value[0] + value[1]) / 2
        return species

    def run_single(self,
                   starting_concentrations: dict,
                   solver: 'ODESolver' = None) -> SingleModelResult:
        """Run model with single parameter values.
        
        Args:
            starting_concentrations: Initial species concentrations
            solver: ODE solver to use (default: SciPySolver)
            
        Returns:
            SingleModelResult containing simulation results
        """
        if solver is None:
            solver = SciPySolver()

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
                  sampler: 'Sampler' = None,
                  solver: 'ODESolver' = None) -> MultiModelResult:
        """Run model with parameter distributions for uncertainty analysis.
        
        Args:
            starting_concentrations: Initial species concentrations (can include distributions)
            sampler: Sampler for parameter distributions (default: ScipyDist_Sampler)
            solver: ODE solver to use (default: SciPySolver)
            
        Returns:
            MultiModelResult containing ensemble simulation results
        """
        if sampler is None:
            sampler = ScipyDist_Sampler(num_samples=1000)
        if solver is None:
            solver = SciPySolver()

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

    def run_multi_vectorized(self,
                           starting_concentrations: dict,
                           sampler: 'Sampler' = None,
                           solver: 'ODESolver' = None) -> MultiModelResult:
        """Run model with parameter distributions using vectorized solver.
        
        This method uses vectorized solvers to process multiple parameter sets
        simultaneously, which can be significantly faster than run_multi.
        
        Args:
            starting_concentrations: Initial species concentrations (can include distributions)
            sampler: Sampler for parameter distributions (default: ScipyDist_Sampler)
            solver: Vectorized ODE solver to use (default: SciPyVectorizedSolver)
            
        Returns:
            MultiModelResult containing ensemble simulation results
        """
        from kinetics.solvers.scipy_vectorized_solver import SciPyVectorizedSolver
        
        if sampler is None:
            sampler = ScipyDist_Sampler(num_samples=1000)
        if solver is None:
            solver = SciPyVectorizedSolver()

        # Get the default values for all parameters and species
        species, parameters = self._parameters_and_species_from_reactions()

        # Set the indexes in the reactions
        self._setup_model(list(species.keys()), list(parameters.keys()))

        # Set the starting concentrations - can contain distributions or fixed values
        species.update(starting_concentrations)

        # Get the parameter distributions from the reactions
        parameter_distributions = self.get_parameter_distributions()

        # Sample parameters and species
        samples = sampler.sample(parameter_distributions, species)

        # Convert samples to batch format
        parameter_samples = []
        species_samples = []
        for parameter_dict, species_dict in samples:
            parameters.update(parameter_dict)
            species.update(species_dict)
            
            # Extract in the same order as the model
            parameter_names, parameter_values = zip(*parameters.items())
            species_names, species_values = zip(*species.items())
            
            parameter_samples.append(list(parameter_values))
            species_samples.append(list(species_values))

        parameter_samples = np.array(parameter_samples)
        species_samples = np.array(species_samples)

        # Use the first sample to get the species names and values for the solver
        species_names, _ = zip(*species.items())
        species_values = species_samples[0]  # Initial conditions from first sample

        # Run vectorized solver
        y_batch = solver.run(self._reactions, species_names, species_values, 
                           parameter_samples, self.ts)

        # Convert batch results to individual results for MultiModelResult
        results = []
        for i in range(parameter_samples.shape[0]):
            results.append(y_batch[:, i, :])

        return MultiModelResult(self, results, list(species.keys()))

    def run_multi_multiprocessing(self,
                                 starting_concentrations: dict,
                                 sampler: 'Sampler' = None,
                                 solver: 'ODESolver' = None) -> MultiModelResult:
        """Run model with parameter distributions using multiprocessing solver.
        
        This method uses multiprocessing to solve multiple ODE systems in parallel
        across CPU cores for better performance on multi-core systems.
        
        Args:
            starting_concentrations: Initial species concentrations (can include distributions)
            sampler: Sampler for parameter distributions (default: ScipyDist_Sampler)
            solver: Multiprocessing ODE solver to use (default: SciPyMultiprocessingSolver)
            
        Returns:
            MultiModelResult containing ensemble simulation results
        """
        from kinetics.solvers.scipy_multiprocessing_solver import SciPyMultiprocessingSolver
        
        if sampler is None:
            sampler = ScipyDist_Sampler(num_samples=1000)
        if solver is None:
            solver = SciPyMultiprocessingSolver()

        # Get the default values for all parameters and species
        species, parameters = self._parameters_and_species_from_reactions()

        # Set the indexes in the reactions
        self._setup_model(list(species.keys()), list(parameters.keys()))

        # Set the starting concentrations - can contain distributions or fixed values
        species.update(starting_concentrations)

        # Get the parameter distributions from the reactions
        parameter_distributions = self.get_parameter_distributions()

        # Sample parameters and species
        samples = sampler.sample(parameter_distributions, species)

        # Convert samples to batch format
        parameter_samples = []
        species_samples = []
        for parameter_dict, species_dict in samples:
            parameters.update(parameter_dict)
            species.update(species_dict)
            
            # Extract in the same order as the model
            parameter_names, parameter_values = zip(*parameters.items())
            species_names, species_values = zip(*species.items())
            
            parameter_samples.append(list(parameter_values))
            species_samples.append(list(species_values))

        parameter_samples = np.array(parameter_samples)
        species_samples = np.array(species_samples)

        # Use the first sample to get the species names and values for the solver
        species_names, _ = zip(*species.items())
        species_values = species_samples[0]  # Initial conditions from first sample

        # Run multiprocessing solver
        y_batch = solver.run(self._reactions, species_names, species_values, 
                           parameter_samples, self.ts)

        # Convert batch results to individual results for MultiModelResult
        results = []
        for i in range(parameter_samples.shape[0]):
            results.append(y_batch[:, i, :])

        return MultiModelResult(self, results, list(species.keys()))

    # TODO
    ## Multi output
    ## Look into vectorisation and gpu acceleration
















