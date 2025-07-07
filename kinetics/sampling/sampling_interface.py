from abc import ABC, abstractmethod
from typing import Callable, List, Tuple

class Sampler(ABC):
    @abstractmethod
    def __init__(self, num_samples):
        self.num_samples = num_samples

    @abstractmethod
    def sample(self,
               parameter_distributions: dict,
               species_distributions: dict) -> List[Tuple[dict, dict]]:
        pass








