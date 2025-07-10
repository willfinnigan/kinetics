from __future__ import annotations

import numpy as np
import copy
from typing import Any

try:
    import jax.numpy as jnp
    HAS_JAX = True
except ImportError:
    jnp = None
    HAS_JAX = False

def calculate_yprime_old(y: np.ndarray, rate: float, substrates: list[str], products: list[str], substrate_names: list[str]) -> np.ndarray:
    """Calculate derivative array for ODE integration.

    Args:
        y: Current species concentrations array
        rate: Reaction rate value
        substrates: Substrate names (rate subtracted)
        products: Product names (rate added)
        substrate_names: Ordered species names for indexing

    Returns:
        Derivative array with rate contributions
    """

    # Create zeros array compatible with input array type (JAX or NumPy)
    if HAS_JAX and hasattr(y, 'shape') and str(type(y)).startswith('<class \'jax'):
        y_prime = jnp.zeros_like(y)
    else:
        y_prime = np.zeros(len(y))

    for name in substrates:
        idx = substrate_names.index(name)
        if HAS_JAX and hasattr(y, 'shape') and str(type(y)).startswith('<class \'jax'):
            y_prime = y_prime.at[idx].add(-rate)
        else:
            y_prime[idx] -= rate

    for name in products:
        idx = substrate_names.index(name)
        if HAS_JAX and hasattr(y, 'shape') and str(type(y)).startswith('<class \'jax'):
            y_prime = y_prime.at[idx].add(rate)
        else:
            y_prime[idx] += rate

    return y_prime


def calculate_yprime(y: np.ndarray,
                     rate: float,
                     substrate_indices: list[int],
                     product_indices: list[int]) -> np.ndarray:
    """Calculate derivative array using pre-computed indices."""

    # Create zeros array compatible with input array type (JAX or NumPy)
    if HAS_JAX and hasattr(y, 'shape') and str(type(y)).startswith('<class \'jax'):
        y_prime = jnp.zeros_like(y)
        y_prime = y_prime.at[substrate_indices].add(-rate)
        y_prime = y_prime.at[product_indices].add(rate)
    else:
        y_prime = np.zeros(len(y))
        y_prime[substrate_indices] -= rate
        y_prime[product_indices] += rate

    return y_prime

def check_positive(y_prime: np.ndarray) -> np.ndarray:
    """Ensure species concentrations remain non-negative.
    
    Args:
        y_prime: Derivative array
        
    Returns:
        Modified derivative array with negative values set to zero
    """

    for i in range(len(y_prime)):
        if y_prime[i] < 0:
            y_prime[i] = 0

    return y_prime

class Reaction:
    """Base class for all reaction types.
    
    Provides common functionality for parameter management, species handling,
    and rate calculation. Subclasses implement specific kinetic equations.
    
    Attributes:
        parameters (dict): Fixed parameter values
        parameter_distributions (dict): Parameter probability distributions
        substrates (list): Substrate species names
        products (list): Product species names
        reaction_substrate_names (list): All species involved in reaction
        parameter_names (list): All parameter names
        modifiers (list): Reaction modifiers (e.g., inhibitors)
        check_positive (bool): Whether to enforce non-negative concentrations
    """

    def __init__(self):
        """Initialize reaction with empty parameter and species lists."""

        # These are set by the user
        self.parameters = {}
        self.parameter_distributions = {}

        # indexes used to access values during model run
        self.reaction_substrate_indexes = []
        self.parameter_indexes = []


        # These are set when the reaction is set up
        self.reaction_substrate_names = []
        self.parameter_names = []
        self.substrates = []
        self.products = []
        self.substrate_indexes = []
        self.product_indexes = []

        # These are added as needed
        self.modifiers = []
        self.check_positive = False
        self.check_limits_functions = []

    def set_parameter_defaults_to_mean(self) -> None:
        """Set default parameter values from distribution means."""
        for name in self.parameter_distributions:
            if name not in self.parameters:
                if type(self.parameter_distributions[name]) == list or type(self.parameter_distributions[name]) == tuple:
                    self.parameters[name] = (self.parameter_distributions[name][0] + self.parameter_distributions[name][1]) / 2
                else:
                    self.parameters[name] = self.parameter_distributions[name].mean()

    def setup_reaction(self, species_names: list[str], parameter_names: list[str]) -> None:
        """Set up reaction indexes for efficient parameter/species access.
        
        Args:
            species_names: Ordered list of all species in model
            parameter_names: Ordered list of all parameters in model
        """

        # get indexes
        self.reaction_substrate_indexes = [species_names.index(name) for name in self.reaction_substrate_names]
        self.substrate_indexes = [species_names.index(name) for name in self.reaction_substrate_names]
        self.product_indexes = [species_names.index(name) for name in self.products]
        self.parameter_indexes = [parameter_names.index(name) for name in self.parameter_names]

        # set up modifiers
        for modifier in self.modifiers:
            modifier.get_substrate_indexes(self.reaction_substrate_names)
            modifier.get_parameter_indexes(self.parameter_names)

    def add_modifier(self, modifier: Any) -> None:
        """Add a modifier (e.g., inhibitor) to the reaction.
        
        Args:
            modifier: Modifier object with parameter_names and substrate_names
        """
        for name in modifier.parameter_names:
            if name not in self.parameter_names:
                self.parameter_names.append(name)

        for name in modifier.substrate_names:
            if name not in self.reaction_substrate_names:
                self.reaction_substrate_names.append(name)

        self.modifiers.append(modifier)

    def calculate_rate(self, substrates: list[float], parameters: list[float]) -> float:
        """Calculate reaction rate
        
        Args:
            substrates: Current substrate concentrations
            parameters: Parameter values
            
        Returns:
            Reaction rate
        """
        return 0

    def calculate_rate_batch(self, substrates_batch: np.ndarray, parameters_batch: np.ndarray) -> np.ndarray:
        """Calculate reaction rates for multiple parameter sets.
        
        This method provides vectorized rate calculation for batch processing.
        Default implementation processes each parameter set individually.
        Subclasses should override this for true vectorization.
        
        Args:
            substrates_batch: Substrate concentrations (n_samples, n_substrates)
            parameters_batch: Parameter values (n_samples, n_parameters)
            
        Returns:
            Rate array (n_samples,)
        """
        # n_samples = substrates_batch.shape[0]
        # rates = np.zeros(n_samples)
        
        # for i in range(n_samples):
        #     substrates_list = substrates_batch[i].tolist()
        #     parameters_list = parameters_batch[i].tolist()
        #     rates[i] = self.calculate_rate(substrates_list, parameters_list)
        
        # return rates
        raise NotImplementedError("Subclasses should implement calculate_rate_batch for vectorized rate calculation.")

    def reaction(self, y: np.ndarray, substrate_names: list[str], parameter_values: list[float]) -> np.ndarray:
        """Calculate rate and return species concentration derivatives.
        
        Args:
            y: Current species concentrations
            substrate_names: Ordered species names
            parameter_values: Parameter values
            
        Returns:
            Derivative array (dy/dt)
        """

        # Get the substrates from y using the substrate indexes
        substrates = []
        for index in self.reaction_substrate_indexes:
            substrates.append(y[index])

        # Get the parameters using the parameter indexes
        parameters = []
        for index in self.parameter_indexes:
            parameters.append(parameter_values[index])

        # calculate the effects of any modifiers
        for modifier in self.modifiers:
            substrates, parameters = modifier.calc_modifier(substrates, parameters)

        # calculate the rate (this function is modified by the user)
        rate = self.calculate_rate(substrates, parameters)

        # calculate the change in substrate concentrations (y_prime)
        y_prime = calculate_yprime(y, rate, self.substrate_indexes, self.product_indexes)
        y_prime = self.modify_product(y_prime, substrate_names)

        if self.check_positive == True:
            y_prime = check_positive(y_prime)

        return y_prime

    def reaction_batch(self, y_batch: np.ndarray, substrate_names: list[str], 
                      parameter_values_batch: np.ndarray) -> np.ndarray:
        """Calculate rates for multiple parameter sets simultaneously.
        
        This method provides batch processing for vectorized solvers.
        Uses vectorized rate calculation when available.
        
        Args:
            y_batch: Species concentrations (n_samples, n_species)
            substrate_names: Ordered species names
            parameter_values_batch: Parameter values (n_samples, n_parameters)
            
        Returns:
            Derivative array (n_samples, n_species)
        """
        n_samples = y_batch.shape[0]
        n_species = y_batch.shape[1]
        
        # Get substrates for all samples
        substrates_batch = np.zeros((n_samples, len(self.substrate_indexes)))
        for i, index in enumerate(self.substrate_indexes):
            substrates_batch[:, i] = y_batch[:, index]
        
        # Get parameters for all samples
        parameters_batch = np.zeros((n_samples, len(self.parameter_indexes)))
        for i, index in enumerate(self.parameter_indexes):
            parameters_batch[:, i] = parameter_values_batch[:, index]
        
        # Calculate rates for all samples using vectorized method
        rates_batch = self.calculate_rate_batch(substrates_batch, parameters_batch)
        
        # Calculate derivatives for all samples
        y_prime_batch = np.zeros((n_samples, n_species))

        # Compute indices once during initialization or setup
        substrate_indices = [substrate_names.index(name) for name in self.substrates]
        product_indices = [substrate_names.index(name) for name in self.products]

        # Then use the pre-computed indices
        y_prime_batch[:, substrate_indices] -= rates_batch
        y_prime_batch[:, product_indices] += rates_batch
        
        return y_prime_batch

    def modify_product(self, y_prime: np.ndarray, substrate_names: list[str]) -> np.ndarray:
        """Modify product formation (can be overridden by subclasses).
        
        Args:
            y_prime: Current derivative array
            substrate_names: Species names
            
        Returns:
            Modified derivative array
        """
        return y_prime

    def sampling_limits(self, parameter_dict: dict) -> bool:
        """Check if parameter values are within acceptable limits.
        
        Args:
            parameter_dict: Parameter values to check
            
        Returns:
            True if parameters are within limits, False otherwise
        """
        for func in self.check_limits_functions:
            if func(parameter_dict) == False:
                return False

        return True
