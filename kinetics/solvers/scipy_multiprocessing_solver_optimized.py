from __future__ import annotations

import numpy as np
from scipy import integrate
from typing import TYPE_CHECKING
from multiprocessing import Pool, cpu_count
from concurrent.futures import ProcessPoolExecutor
import os

from kinetics.solvers.solver_interface import ODESolver

if TYPE_CHECKING:
    from kinetics.models.reaction_class import Reaction


def _solve_single_system_optimized(args):
    """Optimized single system solver for multiprocessing."""
    (reactions, species_names, species_values, parameter_values, time, 
     method, rtol, atol, max_step, mxsteps) = args
    
    def deriv(y, t, reactions, species_names, parameter_values):
        """Calculate derivatives for single system."""
        yprime = np.zeros(len(y))
        for reaction in reactions:
            yprime += reaction.reaction(y, species_names, parameter_values)
        return yprime
    
    y0 = np.array(species_values, dtype=float)
    
    # Use odeint as it's typically faster for simple systems
    y = integrate.odeint(
        deriv, y0, time, 
        args=(reactions, species_names, parameter_values),
        mxstep=mxsteps,
        rtol=rtol,
        atol=atol
    )
    
    return y


class SciPyMultiprocessingSolverOptimized(ODESolver):
    """Optimized multiprocessing solver with multiple approaches."""

    def __init__(self, method: str = 'odeint', rtol: float = 1e-6, 
                 atol: float = 1e-9, max_step: float = np.inf, 
                 mxsteps: int = 5000, n_processes: int = None, 
                 mp_method: str = 'imap_unordered'):
        """Initialize optimized multiprocessing solver.
        
        Args:
            method: Integration method
            rtol: Relative tolerance
            atol: Absolute tolerance
            max_step: Maximum step size
            mxsteps: Maximum number of steps
            n_processes: Number of processes (default: all available cores)
            mp_method: Multiprocessing method ('map', 'imap', 'imap_unordered', 'futures')
        """
        self.method = method
        self.rtol = rtol
        self.atol = atol
        self.max_step = max_step
        self.mxsteps = mxsteps
        self.mp_method = mp_method
        
        # Use all available processors, including logical cores
        if n_processes is None:
            if hasattr(os, 'sched_getaffinity'):
                self.n_processes = len(os.sched_getaffinity(0))
            else:
                self.n_processes = cpu_count()
        else:
            self.n_processes = n_processes

    def run(self, reactions: list['Reaction'], species_names: list[str], 
            species_values: list[float], parameter_values_batch: np.ndarray, 
            time: np.ndarray) -> np.ndarray:
        """Solve ODE systems in parallel using optimized multiprocessing."""
        n_samples = parameter_values_batch.shape[0]
        n_species = len(species_values)
        
        # Skip multiprocessing for very small batches
        if n_samples <= self.n_processes:
            return self._solve_sequential(reactions, species_names, species_values, 
                                        parameter_values_batch, time)
        
        # Prepare arguments for each system
        args_list = []
        for i in range(n_samples):
            args = (
                reactions, species_names, species_values, parameter_values_batch[i], time,
                self.method, self.rtol, self.atol, self.max_step, self.mxsteps
            )
            args_list.append(args)
        
        # Choose multiprocessing method
        if self.mp_method == 'futures':
            results = self._solve_with_futures(args_list)
        elif self.mp_method == 'imap':
            results = self._solve_with_imap(args_list)
        elif self.mp_method == 'imap_unordered':
            results = self._solve_with_imap_unordered(args_list, n_samples)
        else:  # default to map
            results = self._solve_with_map(args_list)
        
        # Convert results to the expected format
        solution = np.zeros((len(time), n_samples, n_species))
        for i, result in enumerate(results):
            solution[:, i, :] = result
        
        return solution
    
    def _solve_sequential(self, reactions, species_names, species_values, 
                         parameter_values_batch, time):
        """Solve sequentially for small batches."""
        n_samples = parameter_values_batch.shape[0]
        n_species = len(species_values)
        solution = np.zeros((len(time), n_samples, n_species))
        
        for i in range(n_samples):
            args = (
                reactions, species_names, species_values, parameter_values_batch[i], time,
                self.method, self.rtol, self.atol, self.max_step, self.mxsteps
            )
            result = _solve_single_system_optimized(args)
            solution[:, i, :] = result
        
        return solution
    
    def _solve_with_map(self, args_list):
        """Solve using Pool.map()."""
        with Pool(processes=self.n_processes) as pool:
            results = pool.map(_solve_single_system_optimized, args_list)
        return results
    
    def _solve_with_imap(self, args_list):
        """Solve using Pool.imap()."""
        with Pool(processes=self.n_processes) as pool:
            results = list(pool.imap(_solve_single_system_optimized, args_list))
        return results
    
    def _solve_with_imap_unordered(self, args_list, n_samples):
        """Solve using Pool.imap_unordered() and preserve order."""
        with Pool(processes=self.n_processes) as pool:
            # Add index to preserve order
            indexed_args = [(i, args) for i, args in enumerate(args_list)]
            
            def solve_with_index(indexed_args):
                i, args = indexed_args
                result = _solve_single_system_optimized(args)
                return i, result
            
            unordered_results = pool.imap_unordered(solve_with_index, indexed_args)
            
            # Restore order
            results = [None] * n_samples
            for i, result in unordered_results:
                results[i] = result
        
        return results
    
    def _solve_with_futures(self, args_list):
        """Solve using concurrent.futures.ProcessPoolExecutor."""
        with ProcessPoolExecutor(max_workers=self.n_processes) as executor:
            futures = [executor.submit(_solve_single_system_optimized, args) 
                      for args in args_list]
            results = [future.result() for future in futures]
        return results