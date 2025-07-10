from __future__ import annotations

import numpy as np
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
        # Convert indices to tuple for JAX compatibility
        y_prime = y_prime.at[tuple(substrate_indices)].add(-rate)
        y_prime = y_prime.at[tuple(product_indices)].add(rate)
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
        self.substrate_indexes = []
        self.parameter_indexes = []


        # These are set when the reaction is set up
        self.reaction_substrate_names = []
        self.parameter_names = []
        self.substrates = []
        self.products = []

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
        self.substrate_indexes = [species_names.index(name) for name in self.reaction_substrate_names]
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
        for index in self.substrate_indexes:
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
        substrate_indices = [substrate_names.index(name) for name in self.substrates]
        product_indices = [substrate_names.index(name) for name in self.products]
        y_prime = calculate_yprime(y, rate, substrate_indices, product_indices)
        y_prime = self.modify_product(y_prime, substrate_names)

        if self.check_positive == True:
            y_prime = check_positive(y_prime)

        return y_prime


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
