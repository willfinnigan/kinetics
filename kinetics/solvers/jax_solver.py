from __future__ import annotations

from kinetics.solvers.solver_interface import ODESolver
import numpy as np
import pandas as pd
import jax.numpy as jnp
from jax import jit, devices, device_put
from diffrax import diffeqsolve, ODETerm, Dopri5, SaveAt, PIDController
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kinetics.models.reaction_class import Reaction


class JaxSolver(ODESolver):
    """JAX-based ODE solver using Diffrax for high-performance computation.
    
    Uses JAX for automatic differentiation and GPU acceleration (not yet supported),
    with the Diffrax library for numerical integration.
    
    This solver is optimized for performance-critical applications
    and supports vectorized operations.
    """
    def deriv(self, t: float, y: jnp.ndarray, args: tuple) -> jnp.ndarray:
        """Calculate derivatives for JAX-based ODE integration.
        
        Called by diffrax.diffeqsolve at each integration step.
        Computes the rate of change for each species using JAX arrays.
        
        Args:
            t: Current time
            y: Current species concentrations (JAX array)
            args: Tuple containing (reactions, species_names, parameter_values)
            
        Returns:
            JAX array of derivatives (dy/dt) for each species
        """
        (reactions, species_names, parameter_values) = args

        # y is now a JAX array, and reaction_class.reaction is expected to handle it.
        # Initialize yprime as a JAX array.
        yprime = jnp.zeros_like(y)

        for reaction_class in reactions:
            yprime += reaction_class.reaction(y, species_names, parameter_values)

        return yprime

    def run(self, reactions: list['Reaction'], species_names: list[str], 
            species_values: list[float], parameter_values: list[float], 
            time: np.ndarray) -> np.ndarray:
        """Solve the ODE system using JAX and Diffrax.
        
        Args:
            reactions: List of reaction objects
            species_names: Ordered species names
            species_values: Initial species concentrations
            parameter_values: Parameter values
            time: Time points for integration
            
        Returns:
            Solution array with shape (n_timepoints, n_species)
        """
        
        y0 = jnp.array(species_values)
        time = jnp.asarray(time)

        # Package arguments for the derivative function
        deriv_args = (reactions, species_names, parameter_values)

        term = ODETerm(self.deriv)
        solver = Dopri5()

        saveat = SaveAt(ts=time)

        # Using default rtol/atol for PIDController, can be tuned if needed.
        stepsize_controller = PIDController(rtol=1e-7, atol=1e-7)

        # Determine initial dt0; a small fraction of the first time step or a fixed small value.
        # Ensure dt0 is a concrete Python float.
        if len(time) > 1:
            # Convert JAX array elements to NumPy/Python float for this calculation
            time_np = np.asarray(time)
            dt_val = (time_np[1] - time_np[0]) / 10.0
            dt0 = float(dt_val)
        else:
            # Default dt0 if only one time point or for safety.
            dt0 = 0.1 # This is already a Python float.

        solution = diffeqsolve(term,
                               solver,
                               t0=time[0],
                               t1=time[-1],
                               dt0=dt0,
                               y0=y0,
                               args=deriv_args,
                               saveat=saveat,
                               stepsize_controller=stepsize_controller,
                               max_steps=100_000)  # Increased max_steps for stiff/long problems

        # Convert JAX array solution to NumPy array for compatibility
        y = np.asarray(solution.ys)
        return y
