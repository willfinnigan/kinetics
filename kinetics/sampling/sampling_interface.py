from typing import Callable

SamplingMethod = Callable

def random_sample(species_distributions: dict,
                  parameter_distributions: dict,
                  num: int):
    """Take random samples from scipy distributions."""

    samples = []
    for i in range(num):  # make samples
        parameter_dict, species_dict = {}, {}  # Initialize empty dicts for this sample





