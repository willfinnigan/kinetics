import numpy as np
import pandas as pd
from scipy import integrate
import matplotlib.pyplot as plt

from kinetics.solvers.scipy_solver import SciPy_Solver

class TimeSeries(object):

    def __init__(self, start: int, end: int, steps: int, mxsteps=10000):
        """
        This function sets the time parameters for the model.  This is how long the model will simulate

        Args:
            start (int): the start time - usually 0
            end (int): the end time (default is 100)
            steps (int): the number of timepoints for the output
        """

        self.start = start
        self.end = end
        self.steps = steps
        self.t = np.linspace(self.start, self.end, self.steps)
        self.mxsteps = mxsteps



class Model(object):
    """
    The model class is central.  It inherits from a list.  Reactions are appended to this list to build the model.
    Upon creating a new object logging can be turned off by passing in logging=False

    1.  Create model, append reactions and set time and species.
    2.  setup_model()
    3.  run_model()

    Attributes:
        species (dict): The starting species concentrations.  For example {'Substrate_1' : 100}

        species_distributions (dict): The starting species concentrations, with uncertainty using probability distributions from scipy.
                                      For example {'Substrate_1' : norm(100, 10)}

        parameters (dict): Parameters.  These are loaded from the appended reactions upon running setup_model(). For example {'param_1' : 100}
        parameter_distributions (dict):  Parameter scipy distributions.  These are loaded from the appended reactions upon running setup_model(). For example {'param_1' : norm(100, 10)}

        y (numpy array): a numpy array of 2 dimensions. Time by substrate.  Filled upon running run_model()
                         The first dimension gives a list of all the substrate concentrations at that timepoint.
                         The first dimension is the same size as self.time.
                         Each index in self.time relates to an index in the first dimension of y.

        logging (bool): True gives text feedback upon running some commands

        start (int): Model start time
        end (int): Model end time
        steps (int): The number of timpoints in the model output
        mxsteps (int): mxsteps used by scipy.integrate.odeint
        time (np.linspace(self.start, self.end, self.steps)):  The timepoints of the model

    """

    def __init__(self):

        """ Reactions - a list of reaction classes """
        self.reactions = []
        self.species = {}
        self.parameters = {}
        self.timeseries: TimeSeries = None

        """ Species and parameters used when the model is ran. These are changed each run when doing ua/sa """
        self.run_model_species = {}
        self.run_model_species_names = []
        self.run_model_species_starting_values = []
        self.run_model_parameters = {}

        self.y = []

    def set_time(self, start: int, end: int, steps: int, mxsteps=10000):
        self.timeseries = TimeSeries(start, end, steps, mxsteps)

    def set_species(self, species):
        self.species.update(species)

    def add_reaction(self, reaction):
        self.reactions.append(reaction)
        reaction.set_parameter_defaults_to_mean()

        # Load the model's parameters from the reaction's parameters
        for name in reaction.parameters:
            if name not in self.parameters:
                self.parameters[name] = reaction.parameters[name]
            else:
                print('Warning - parameter ' + name + ' already set in model, not overwriting with reaction value')

        # Load the model's species from the reaction's species
        for substrate in reaction.substrates + reaction.products + reaction.reaction_substrate_names:
            if substrate not in self.species:
                self.species[substrate] = 0

    def run_model(self):
        solver = SciPy_Solver()
        self.y = solver.run(self.reactions, self.species, self.parameters, self.timeseries.t)
        return self.y


    # Export results as dataframe and plot
    def results_dataframe(self):
        """
        Gives the results of a model run as a dataframe

        Returns:
            Pandas dataframe of results
        """
        ys_at_t = {'Time' : self.timeseries.t}
        species_names = list(self.species.keys())

        for i in range(len(species_names)):
            name = species_names[i]
            ys_at_t[name] = []

            for t in range(len(self.timeseries.t)):
                ys_at_t[name].append(self.y[t][i])

        df = pd.DataFrame(ys_at_t)

        return df
