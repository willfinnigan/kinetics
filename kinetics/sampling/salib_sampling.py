from __future__ import annotations

from typing import List, Tuple
import numpy as np
from SALib.sample import latin, saltelli

from kinetics.sampling.sampling_interface import Sampler


def build_salib_problem(parameter_distributions: dict, species_distributions: dict, log_parameters: List[str] = None) -> dict:
    """Build SALib problem dictionary for sensitivity analysis.
    
    Args:
        parameter_distributions: Dictionary of parameter distributions
        species_distributions: Dictionary of species distributions
        log_parameters: List of parameters to sample in log space
        
    Returns:
        SALib problem dictionary
    """
    if log_parameters is None:
        log_parameters = []
    
    names = []
    bounds = []
    
    # Add parameter bounds
    for name, distribution in parameter_distributions.items():
        names.append(name)
        bounds.append([distribution[0], distribution[1]])
    
    # Add species bounds (only for distributions, not single values)
    for name, value in species_distributions.items():
        if isinstance(value, (tuple, list)) and len(value) == 2:
            names.append(name)
            bounds.append([value[0], value[1]])
    
    problem = {
        'num_vars': len(names),
        'names': names,
        'bounds': bounds
    }
    
    # Convert to log space if needed
    return problem_to_log_space(problem, log_parameters)


def problem_to_log_space(problem: dict, log_parameters: List[str]) -> dict:
    """Convert specified parameters to log space for sampling.
    
    Args:
        problem: SALib problem dictionary
        log_parameters: Parameter names to convert to log space
        
    Returns:
        Modified problem dictionary with log space bounds
    """
    for i, name in enumerate(problem['names']):
        if name in log_parameters:
            lower = np.log(problem['bounds'][i][0])
            upper = np.log(problem['bounds'][i][1])
            problem['bounds'][i] = [lower, upper]
    return problem


def samples_to_normal_space(samples: np.ndarray, problem: dict, log_parameters: List[str]) -> np.ndarray:
    """Convert log space samples back to normal space.
    
    Args:
        samples: Sample array in log space
        problem: SALib problem dictionary
        log_parameters: Parameter names that were in log space
        
    Returns:
        Sample array in normal space
    """
    for i, name in enumerate(problem['names']):
        if name in log_parameters:
            for j in range(len(samples)):
                samples[j][i] = float(np.exp(samples[j][i]))
    return samples


def parse_samples(samples: np.ndarray, 
                 parameter_distributions: dict, 
                 species_distributions: dict) -> List[Tuple[dict, dict]]:
    """Parse samples into parameter and species dictionaries.
    
    Args:
        samples: Sample array from SALib
        parameter_distributions: Dictionary of parameter distributions
        species_distributions: Dictionary of species distributions
        
    Returns:
        List of tuples (parameter_dict, species_dict) for each sample
    """
    parsed_samples = []
    parameter_names = list(parameter_distributions.keys())
    
    # Only include species that are distributions (not single values)
    species_distribution_names = [name for name, value in species_distributions.items() 
                                  if isinstance(value, (tuple, list)) and len(value) == 2]
    
    for sample in samples:
        # Extract parameters
        parameter_dict = {}
        for i, name in enumerate(parameter_names):
            parameter_dict[name] = sample[i]
        
        # Extract species distributions
        species_dict = {}
        for i, name in enumerate(species_distribution_names):
            species_dict[name] = sample[i + len(parameter_names)]
        
        # Add single value species
        for name, value in species_distributions.items():
            if not isinstance(value, (tuple, list)) or len(value) != 2:
                species_dict[name] = value
        
        parsed_samples.append((parameter_dict, species_dict))
    
    return parsed_samples


class SalibLatinHypercubeSampler(Sampler):
    """Latin Hypercube Sampling using SALib.
    
    Generates samples using Latin Hypercube Sampling (LHS) which provides
    good space-filling properties for parameter exploration.
    
    Args:
        num_samples: Number of samples to generate
        log_parameters: List of parameters to sample in log space
    """
    
    def __init__(self, num_samples: int, log_parameters: List[str] = None):
        """Initialize Latin Hypercube sampler.
        
        Args:
            num_samples: Number of samples to generate
            log_parameters: Parameter names to sample in log space
        """
        self.num_samples = num_samples
        self.log_parameters = log_parameters if log_parameters is not None else []
    
    def sample(self,
               parameter_distributions: dict,
               species_distributions: dict) -> List[Tuple[dict, dict]]:
        
        # Build SALib problem
        problem = build_salib_problem(parameter_distributions, species_distributions, self.log_parameters)
        
        # Generate samples
        samples = latin.sample(problem, self.num_samples)
        
        # Convert log space back to normal space
        samples = samples_to_normal_space(samples, problem, self.log_parameters)
        
        # Parse samples into parameter and species dictionaries
        return parse_samples(samples, parameter_distributions, species_distributions)


class SalibSaltelliSampler(Sampler):
    """Saltelli sampling for sensitivity analysis using SALib.
    
    Generates samples using the Saltelli method which is designed for
    Sobol sensitivity analysis. Creates N*(2D+2) samples where N is
    the base sample size and D is the number of parameters.
    
    Args:
        num_samples: Base number of samples (actual samples will be larger)
        second_order: Whether to calculate second-order indices
        log_parameters: List of parameters to sample in log space
    """
    
    def __init__(self, num_samples: int, second_order: bool = False, log_parameters: List[str] = None):
        """Initialize Saltelli sampler.
        
        Args:
            num_samples: Base number of samples
            second_order: Whether to calculate second-order indices
            log_parameters: Parameter names to sample in log space
        """
        self.num_samples = num_samples
        self.second_order = second_order
        self.log_parameters = log_parameters if log_parameters is not None else []
        self.problem = None  # Store the problem for sensitivity analysis
    
    def sample(self,
               parameter_distributions: dict,
               species_distributions: dict) -> List[Tuple[dict, dict]]:
        
        # Build SALib problem
        self.problem = build_salib_problem(parameter_distributions, species_distributions, self.log_parameters)
        
        # Generate samples
        samples = saltelli.sample(self.problem, self.num_samples, calc_second_order=self.second_order)
        
        # Convert log space back to normal space
        samples = samples_to_normal_space(samples, self.problem, self.log_parameters)
        
        # Parse samples into parameter and species dictionaries
        return parse_samples(samples, parameter_distributions, species_distributions)