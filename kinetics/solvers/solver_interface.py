from abc import ABC, abstractmethod



class ODESolver(ABC):
    @abstractmethod
    def run(self, reactions, species_names, species_values, parameters, time):
        pass