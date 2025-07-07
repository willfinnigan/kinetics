from typing import List, Tuple
import numpy as np
from SALib.sample import latin, saltelli

from kinetics.sampling.sampling_interface import Sampler


def build_salib_problem(parameter_distributions: dict, species_distributions: dict, log_parameters: List[str] = None) -> dict:
    """Build SALib problem dictionary"""
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
    """Convert specified parameters to log space"""
    for i, name in enumerate(problem['names']):
        if name in log_parameters:
            lower = np.log(problem['bounds'][i][0])
            upper = np.log(problem['bounds'][i][1])
            problem['bounds'][i] = [lower, upper]
    return problem


def samples_to_normal_space(samples: np.ndarray, problem: dict, log_parameters: List[str]) -> np.ndarray:
    """Convert log space samples back to normal space"""
    for i, name in enumerate(problem['names']):
        if name in log_parameters:
            for j in range(len(samples)):
                samples[j][i] = float(np.exp(samples[j][i]))
    return samples


def parse_samples(samples: np.ndarray, 
                 parameter_distributions: dict, 
                 species_distributions: dict) -> List[Tuple[dict, dict]]:
    """Parse samples into parameter and species dictionaries"""
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
    """Latin Hypercube Sampling using SALib"""
    
    def __init__(self, num_samples: int, log_parameters: List[str] = None):
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
    """Saltelli sampling for sensitivity analysis using SALib"""
    
    def __init__(self, num_samples: int, second_order: bool = False, log_parameters: List[str] = None):
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