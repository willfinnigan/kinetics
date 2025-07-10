from __future__ import annotations

import numpy as np
from scipy import integrate
from typing import TYPE_CHECKING
from multiprocessing import Pool, cpu_count
from functools import partial
import os

from kinetics.solvers.solver_interface import ODESolver

if TYPE_CHECKING:
    from kinetics.models.reaction_class import Reaction


def _solve_single_system(args):
    """Solve a single ODE system - used for multiprocessing."""
    (reactions, species_names, species_values, parameter_values, time, 
     method, rtol, atol, max_step, mxsteps) = args
    
    def deriv(y, t, reactions, species_names, parameter_values):
        """Calculate derivatives for single system."""
        yprime = np.zeros(len(y))
        for reaction in reactions:
            yprime += reaction.reaction(y, species_names, parameter_values)
        return yprime
    
    y0 = np.array(species_values, dtype=float)
    
    if method == 'odeint':
        # Use odeint for compatibility with original solver
        y = integrate.odeint(
            deriv, y0, time, 
            args=(reactions, species_names, parameter_values),
            mxstep=mxsteps
        )
    else:
        # Use solve_ivp for more modern approach
        sol = integrate.solve_ivp(
            fun=lambda t, y: deriv(y, t, reactions, species_names, parameter_values),
            t_span=(time[0], time[-1]),
            y0=y0,
            t_eval=time,
            method=method,
            rtol=rtol,
            atol=atol,
            max_step=max_step
        )
        y = sol.y.T
    
    return y


class SciPyMultiprocessingSolver(ODESolver):
    """Multiprocessing-based SciPy solver for parallel ODE solving.
    
    This solver uses the original approach of solving individual ODE systems
    but parallelizes across CPU cores using multiprocessing for better performance.
    """

    def __init__(self, method: str = 'odeint', rtol: float = 1e-6, 
                 atol: float = 1e-9, max_step: float = np.inf, 
                 mxsteps: int = 5000, n_processes: int = None):
        """Initialize multiprocessing solver.
        
        Args:
            method: Integration method ('odeint' or solve_ivp methods)
            rtol: Relative tolerance
            atol: Absolute tolerance
            max_step: Maximum step size
            mxsteps: Maximum number of steps
            n_processes: Number of processes to use (default: CPU count)
        """
        self.method = method
        self.rtol = rtol
        self.atol = atol
        self.max_step = max_step
        self.mxsteps = mxsteps
        self.n_processes = n_processes or cpu_count()

    def run(self, reactions: list['Reaction'], species_names: list[str], 
            species_values: list[float], parameter_values_batch: np.ndarray, 
            time: np.ndarray) -> np.ndarray:
        """Solve ODE systems in parallel using multiprocessing.
        
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
        
        # Prepare arguments for each system
        args_list = []
        for i in range(n_samples):
            args = (
                reactions, species_names, species_values, parameter_values_batch[i], time,
                self.method, self.rtol, self.atol, self.max_step, self.mxsteps
            )
            args_list.append(args)
        
        # Solve systems in parallel
        with Pool(processes=self.n_processes) as pool:
            results = pool.map(_solve_single_system, args_list)
        
        # Convert results to the expected format
        solution = np.zeros((len(time), n_samples, n_species))
        for i, result in enumerate(results):
            solution[:, i, :] = result
        
        return solution


class SciPyMultiprocessingSolverChunked(ODESolver):
    """Chunked multiprocessing solver for better memory efficiency.
    
    This version processes systems in chunks to avoid memory issues with
    very large parameter sets while still using multiprocessing.
    """

    def __init__(self, method: str = 'odeint', rtol: float = 1e-6, 
                 atol: float = 1e-9, max_step: float = np.inf, 
                 mxsteps: int = 5000, n_processes: int = None,
                 chunk_size: int = 1000):
        """Initialize chunked multiprocessing solver.
        
        Args:
            method: Integration method
            rtol: Relative tolerance
            atol: Absolute tolerance
            max_step: Maximum step size
            mxsteps: Maximum number of steps
            n_processes: Number of processes to use (default: CPU count)
            chunk_size: Number of systems to process in each chunk
        """
        self.method = method
        self.rtol = rtol
        self.atol = atol
        self.max_step = max_step
        self.mxsteps = mxsteps
        self.n_processes = n_processes or cpu_count()
        self.chunk_size = chunk_size

    def run(self, reactions: list['Reaction'], species_names: list[str], 
            species_values: list[float], parameter_values_batch: np.ndarray, 
            time: np.ndarray) -> np.ndarray:
        """Solve ODE systems in parallel using chunked processing.
        
        This approach processes systems in chunks to manage memory usage
        while still benefiting from multiprocessing.
        """
        n_samples = parameter_values_batch.shape[0]
        n_species = len(species_values)
        solution = np.zeros((len(time), n_samples, n_species))
        
        # Process in chunks
        for start_idx in range(0, n_samples, self.chunk_size):
            end_idx = min(start_idx + self.chunk_size, n_samples)
            chunk_params = parameter_values_batch[start_idx:end_idx]
            chunk_size = end_idx - start_idx
            
            # Prepare arguments for this chunk
            args_list = []
            for i in range(chunk_size):
                args = (
                    reactions, species_names, species_values, chunk_params[i], time,
                    self.method, self.rtol, self.atol, self.max_step, self.mxsteps
                )
                args_list.append(args)
            
            # Solve chunk in parallel
            with Pool(processes=self.n_processes) as pool:
                chunk_results = pool.map(_solve_single_system, args_list)
            
            # Store results
            for i, result in enumerate(chunk_results):
                solution[:, start_idx + i, :] = result
        
        return solution