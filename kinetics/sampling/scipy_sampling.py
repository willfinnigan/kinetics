from __future__ import annotations

from typing import List, Tuple

import numpy as np

from kinetics.sampling.sampling_interface import Sampler


class ScipyDist_Sampler(Sampler):
    """Sampler using scipy.stats distributions.
    
    Generates samples from parameter and species distributions using
    scipy.stats random variable objects. Supports rejection sampling
    to ensure non-negative values for specified parameters.
    
    Args:
        num_samples: Number of samples to generate
        negative_allowed: List of parameter names that can be negative
    """

    def __init__(self, num_samples: int, negative_allowed: List[str] = None):
        """Initialize scipy distribution sampler.
        
        Args:
            num_samples: Number of samples to generate
            negative_allowed: Parameter names that can have negative values
        """
        self.num_samples = num_samples
        self.negative_allowed = negative_allowed if negative_allowed is not None else []

    def sample(self,
               parameter_distributions: dict,
               species_distributions: dict) -> List[Tuple[dict, dict]]:
        """Generate samples from parameter and species distributions.
        
        Uses rejection sampling to ensure non-negative values where required.
        
        Args:
            parameter_distributions: Dictionary of parameter distributions
            species_distributions: Dictionary of species distributions
            
        Returns:
            List of tuples (parameter_dict, species_dict) for each sample
        """

        samples = []
        for i in range(self.num_samples):  # make samples
            parameter_dict, species_dict = {}, {}  # Initialize empty dicts for this sample

            # Sample parameters from distributions
            for name, distribution in parameter_distributions.items():
                sample = None
                while not self._check_not_neg(sample, name):
                    sample = distribution.rvs()

                if type(sample) == np.ndarray:
                    parameter_dict[name] = sample[0]
                else:
                    parameter_dict[name] = sample

            # Handle species
            for name, value in species_distributions.items():
                if hasattr(value, 'rvs'):  # It's a distribution
                    sample = None
                    while not self._check_not_neg(sample, name):
                        sample = value.rvs()
                    species_dict[name] = sample
                else:  # It's a fixed value
                    species_dict[name] = value

            samples.append((parameter_dict, species_dict))

        return samples

    def _check_not_neg(self, sample, name: str) -> bool:
        """Check if sample is non-negative (unless explicitly allowed).
        
        Args:
            sample: Sample value to check
            name: Parameter name
            
        Returns:
            True if sample is valid, False otherwise
        """
        def check_sample(sample_to_check, name_to_check: str) -> bool:
            """Check individual sample value."""
            if (sample_to_check <= 0) and (name_to_check not in self.negative_allowed):
                return False
            return True

        if type(sample) == np.ndarray:
            for s in sample:
                if check_sample(s, name) == False:
                    return False

        elif type(sample) == np.float64:
            if check_sample(sample, name) == False:
                return False

        elif sample == None:
            return False

        else:
            # Handle other numeric types (int, float, etc.)
            if check_sample(sample, name) == False:
                return False

        return True


