from __future__ import annotations

import numpy as np
from scipy import integrate
from typing import TYPE_CHECKING

from kinetics.solvers.solver_interface import ODESolver

if TYPE_CHECKING:
    from kinetics.models.reaction_class import Reaction


class SciPySolver(ODESolver):
    """SciPy-based ODE solver using scipy.integrate.odeint.
    
    Uses the LSODA algorithm for efficient integration of stiff and non-stiff
    ordinary differential equations.
    
    Args:
        mxsteps: Maximum number of steps allowed during integration
    """

    def __init__(self, mxsteps: int = 5000):
        """Initialize solver with maximum step limit.
        
        Args:
            mxsteps: Maximum number of integration steps
        """
        self.mxsteps = mxsteps

    def deriv(self, y: np.ndarray, t: float, reactions: list['Reaction'], 
              species_names: list[str], parameter_values: list[float]) -> np.ndarray:
        """Calculate derivatives for ODE integration.
        
        Called by scipy.integrate.odeint at each integration step.
        Computes the rate of change for each species based on all reactions.
        
        Args:
            y: Current species concentrations
            t: Current time (required by odeint)
            reactions: List of reaction objects
            species_names: Ordered species names
            parameter_values: Parameter values
            
        Returns:
            Array of derivatives (dy/dt) for each species
        """

        yprime = np.zeros(len(y))

        for reaction_class in reactions:
            yprime += reaction_class.reaction(y, species_names, parameter_values)

        return yprime

    def run(self, 
            reactions: list['Reaction'], 
            species_names: list[str], 
            species_values: list[float], 
            parameter_values: list[float], 
            time: np.ndarray) -> np.ndarray:
        """Solve the ODE system using scipy.integrate.odeint.
        
        Args:
            reactions: List of reaction objects
            species_names: Ordered species names
            species_values: Initial species concentrations
            parameter_values: Parameter values
            time: Time points for integration
            
        Returns:
            Solution array with shape (n_timepoints, n_species)
        """


        y0 = np.array(species_values, dtype=float)
        y = integrate.odeint(self.deriv, y0, time, 
                           args=(reactions, species_names, parameter_values), 
                           mxstep=self.mxsteps)
        return y


