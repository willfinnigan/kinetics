from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kinetics.models.reaction_class import Reaction
    import numpy as np


class ODESolver(ABC):
    """Abstract base class for ODE solvers.
    
    Defines the interface for solving ordinary differential equations
    in kinetic models. Implementations should handle the integration
    of reaction rate equations over time.
    """
    
    @abstractmethod
    def run(self, reactions: list['Reaction'], species_names: list[str], 
            species_values: list[float], parameters: list[float], 
            time: 'np.ndarray') -> 'np.ndarray':
        """Solve the ODE system for given reactions and conditions.
        
        Args:
            reactions: List of reaction objects
            species_names: Ordered list of species names
            species_values: Initial species concentrations
            parameters: Parameter values
            time: Time points for integration
            
        Returns:
            Solution array with shape (n_timepoints, n_species)
        """
        pass