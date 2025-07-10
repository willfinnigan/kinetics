from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Callable, List, Tuple


class Sampler(ABC):
    """Abstract base class for parameter sampling.
    
    Defines the interface for sampling from parameter and species distributions
    to enable uncertainty quantification and sensitivity analysis.
    """
    
    @abstractmethod
    def __init__(self, num_samples: int):
        """Initialize sampler with number of samples.
        
        Args:
            num_samples: Number of samples to generate
        """
        self.num_samples = num_samples

    @abstractmethod
    def sample(self,
               parameter_distributions: dict,
               species_distributions: dict) -> List[Tuple[dict, dict]]:
        """Generate samples from parameter and species distributions.
        
        Args:
            parameter_distributions: Dictionary of parameter distributions
            species_distributions: Dictionary of species distributions
            
        Returns:
            List of tuples (parameter_dict, species_dict) for each sample
        """
        pass








