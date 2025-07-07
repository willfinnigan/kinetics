from typing import List, Tuple

import numpy as np

from kinetics.sampling.sampling_interface import Sampler

class ScipyDist_Sampler(Sampler):

    def __init__(self, num_samples: int, negative_allowed: List[str] = None):
        self.num_samples = num_samples
        self.negative_allowed = negative_allowed if negative_allowed is not None else []

    def sample(self,
               parameter_distributions: dict,
               species_distributions: dict) -> List[Tuple[dict, dict]]:

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

    def _check_not_neg(self, sample, name):
        def check_sample(sample_to_check, name_to_check):
            if (sample_to_check <= 0) and (name_to_check not in self.negative_allowed):
                return False

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


