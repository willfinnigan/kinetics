from __future__ import annotations

import numpy as np
from scipy import integrate
from typing import TYPE_CHECKING

from kinetics.solvers.solver_interface import ODESolver

if TYPE_CHECKING:
    from kinetics.models.reaction_class import Reaction


class SciPyVectorizedSolver(ODESolver):
    """Vectorized SciPy-based ODE solver for multiple parameter sets.
    
    Solves ODEs simultaneously for multiple parameter sets using scipy.integrate.solve_ivp
    with vectorized right-hand side function. This is more efficient than running
    separate integrations for each parameter set.
    
    Args:
        method: Integration method (default: 'LSODA')
        rtol: Relative tolerance
        atol: Absolute tolerance
        max_step: Maximum step size
    """

    def __init__(self, method: str = 'LSODA', rtol: float = 1e-6, 
                 atol: float = 1e-9, max_step: float = np.inf):
        """Initialize vectorized solver.
        
        Args:
            method: Integration method for solve_ivp
            rtol: Relative tolerance
            atol: Absolute tolerance  
            max_step: Maximum step size
        """
        self.method = method
        self.rtol = rtol
        self.atol = atol
        self.max_step = max_step

    def vectorized_deriv(self, t: float, y_flat: np.ndarray, 
                        reactions: list['Reaction'], species_names: list[str], 
                        parameter_values_batch: np.ndarray, 
                        n_species: int, n_samples: int) -> np.ndarray:
        """Calculate derivatives for all parameter sets simultaneously.
        
        Args:
            t: Current time
            y_flat: Flattened state vector (n_species * n_samples,)
            reactions: List of reaction objects
            species_names: Ordered species names
            parameter_values_batch: Parameter values (n_samples, n_parameters)
            n_species: Number of species
            n_samples: Number of parameter sets
            
        Returns:
            Flattened derivative array (n_species * n_samples,)
        """
        # Reshape flat state to (n_samples, n_species)
        y_batch = y_flat.reshape(n_samples, n_species)
        
        # Initialize derivatives
        yprime_batch = np.zeros_like(y_batch)
        
        # Use vectorized reaction calculations
        for reaction in reactions:
            yprime_batch += reaction.reaction_batch(y_batch, species_names, parameter_values_batch)
        
        return yprime_batch.flatten()

    def run(self, reactions: list['Reaction'], species_names: list[str], 
            species_values: list[float], parameter_values_batch: np.ndarray, 
            time: np.ndarray) -> np.ndarray:
        """Solve ODE system for multiple parameter sets.
        
        Args:
            reactions: List of reaction objects
            species_names: Ordered species names
            species_values: Initial species concentrations
            parameter_values_batch: Parameter values (n_samples, n_parameters)
            time: Time points for integration
            
        Returns:
            Solution array with shape (n_timepoints, n_samples, n_species)
        """
        n_samples = parameter_values_batch.shape[0]
        n_species = len(species_values)
        
        # Create initial conditions for all parameter sets
        y0 = np.array(species_values, dtype=float)
        y0_batch = np.tile(y0, n_samples)  # Flatten: [y0_sample0, y0_sample1, ...]
        
        # Solve the vectorized system
        sol = integrate.solve_ivp(
            fun=self.vectorized_deriv,
            t_span=(time[0], time[-1]),
            y0=y0_batch,
            t_eval=time,
            args=(reactions, species_names, parameter_values_batch, n_species, n_samples),
            method=self.method,
            rtol=self.rtol,
            atol=self.atol,
            max_step=self.max_step
        )
        
        # Reshape solution from (n_species*n_samples, n_timepoints) to (n_timepoints, n_samples, n_species)
        solution = sol.y.T.reshape(len(time), n_samples, n_species)
        
        return solution


class SciPyVectorizedSolverOptimized(ODESolver):
    """Optimized vectorized SciPy solver with true vectorization.
    
    This version attempts to vectorize the reaction calculations themselves,
    which requires reactions to support batch processing.
    """

    def __init__(self, method: str = 'LSODA', rtol: float = 1e-6, 
                 atol: float = 1e-9, max_step: float = np.inf):
        self.method = method
        self.rtol = rtol
        self.atol = atol
        self.max_step = max_step

    def vectorized_deriv_optimized(self, t: float, y_flat: np.ndarray, 
                                 reactions: list['Reaction'], species_names: list[str], 
                                 parameter_values_batch: np.ndarray, 
                                 n_species: int, n_samples: int) -> np.ndarray:
        """Optimized derivative calculation with vectorized reactions.
        
        This version assumes reactions can handle batch inputs.
        """
        # Reshape to (n_samples, n_species)
        y_batch = y_flat.reshape(n_samples, n_species)
        
        # Initialize derivatives
        yprime_batch = np.zeros_like(y_batch)
        
        # Try to vectorize reaction calculations
        for reaction in reactions:
            # Check if reaction supports batch processing
            if hasattr(reaction, 'reaction_batch'):
                # Use vectorized reaction method
                yprime_batch += reaction.reaction_batch(y_batch, species_names, parameter_values_batch)
            else:
                # Fall back to individual processing
                for i in range(n_samples):
                    yprime_batch[i] += reaction.reaction(y_batch[i], species_names, parameter_values_batch[i])
        
        return yprime_batch.flatten()

    def run(self, reactions: list['Reaction'], species_names: list[str], 
            species_values: list[float], parameter_values_batch: np.ndarray, 
            time: np.ndarray) -> np.ndarray:
        """Solve using optimized vectorized approach."""
        n_samples = parameter_values_batch.shape[0]
        n_species = len(species_values)
        
        y0 = np.array(species_values, dtype=float)
        y0_batch = np.tile(y0, n_samples)
        
        sol = integrate.solve_ivp(
            fun=self.vectorized_deriv_optimized,
            t_span=(time[0], time[-1]),
            y0=y0_batch,
            t_eval=time,
            args=(reactions, species_names, parameter_values_batch, n_species, n_samples),
            method=self.method,
            rtol=self.rtol,
            atol=self.atol,
            max_step=self.max_step
        )
        
        solution = sol.y.T.reshape(len(time), n_samples, n_species)
        return solution