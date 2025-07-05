from abc import ABC, abstractmethod



class ODESolver(ABC):
    @abstractmethod
    def run(self, reactions, species, parameters, time):
        pass