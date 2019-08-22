import numpy as np
import pandas as pd
from scipy import integrate
import matplotlib.pyplot as plt

class Model(list):
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


    def __init__(self, logging=True):
        # Model inherits from list - reaction classes are held in this self list.
        super(Model, self).__init__()

        """ Time """
        self.start = 0
        self.end = 100
        self.steps = 100
        self.mxsteps = 10000
        self.time = np.linspace(self.start, self.end, self.steps)

        """ Species - used to reset the model, or as the bounds to run ua/sa """
        self.species = {}
        self.species_distributions = {}

        """ Parameters - used to reset the model, or as the bounds to run ua/sa.  Set by self.set_parameters_from_reactions() """
        self.parameters = {}
        self.parameter_distributions = {}

        """ Species and parameters used when the model is ran. These are changed each run when doing ua/sa """
        self.run_model_species = {}
        self.run_model_species_names = []
        self.run_model_species_starting_values = []
        self.run_model_parameters = {}

        self.y = []

        self.logging = logging

    # Time
    def set_time(self, start, end, steps):
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
        self.time = np.linspace(self.start, self.end, self.steps)

    # Setup Model
    def set_parameters_from_reactions(self):
        """
        Sets all the parameter variables from those set in the reaction classes attached to the model

        For each reaction_class, updates self.parameters and self.parameter_distributions with the dictionaries held in each reaction_class.
        This will add new keys, or overwrite existing ones.

        Where only a distribution is set, the median value of this distribution will be used for the parameter value.

        Called by self.setup_model()
        """

        self.run_model_parameters = {}

        if self.logging == True:
            print('-- Setting default parameters, using means of distributions where undefined: --')
        for reaction_class in self:
            reaction_class.set_parameter_defaults_to_mean()
            if self.logging==True:
                print(reaction_class.parameters)

            # if parameter not already set in model, load it from reaction
            for name in reaction_class.parameters:
                if name not in self.parameters:
                    self.parameters[name] = reaction_class.parameters[name]

            # if parameter_distribution not already set in model, load it from reaction
            for name in reaction_class.parameter_distributions:
                if name not in self.parameter_distributions:
                    self.parameter_distributions[name] = reaction_class.parameter_distributions[name]

            # if parameter not set in model, and hasn't been loaded from reaction, take mean of model_distribution
            for name in self.parameter_distributions:
                if name not in self.parameters:
                    self.parameters[name] = self.parameter_distributions[name].mean()
                    if self.logging == True:
                        print(str(name) + ' - ' + str(self.parameters[name]))


            self.run_model_parameters.update(self.parameters)

    def update_species(self, species_dict):
        """
        This func is used by to update starting species values used by the model
        Called by: self.setup_model() and self.reset_model_to_defaults()
        """

        self.run_model_species.update(species_dict)
        self.run_model_species_names = list(self.run_model_species.keys())
        self.run_model_species_starting_values = list(self.run_model_species.values())

    def load_species_from_reactions(self):
        """
        Loads species which are present in one of the reaction_classes but not in
        either self.species or self.species_distributions.  Loads them as self.species[name] = 0.

        Called by self.setup_model()
        """
        if self.logging == True:
            print('-- Load unspecified species as default = 0 --')
        for reaction in self:
            for substrate in reaction.substrates + reaction.products + reaction.reaction_substrate_names:
                if substrate not in self.species:
                    self.species[substrate] = 0
                    if self.logging == True:
                        print(str(substrate) + ' ', end='')
        if self.logging == True:
            print()

    def set_species_defaults_to_mean(self):
        """
        For any species defined in self.species_distributions, but not in self.species,
        set self.species[name] to the median of self.species_distributions[name]

        Called by self.setup_model()
        """

        for name in self.species_distributions:
            if name not in self.species:
                self.species[name] = self.species_distributions[name].mean()

    def setup_model(self):
        """
        Run methods to setup the model.
        1. set_species_defaults_to_median()
        2. load_species_from_reactions()
        3. update_species(self.species())
        4. set_parameters_from_reactions()
        """

        # Species
        self.set_species_defaults_to_mean()
        self.load_species_from_reactions()
        self.update_species(self.species)

        # Parameters
        self.set_parameters_from_reactions()

    # Reset the model
    def reset_reaction_indexes(self):
        """
        Called at the end of run_model() to reset the indexes of the substrates and parameters in the reaction classes.
        May not be necessary - need to look into this.
        """
        for reaction_class in self:
            reaction_class.reset_reaction()

    def reset_model_to_defaults(self):
        """
        Reset the model back to the default settings

        This uses self.species and self.parameters to set the run_model attibutes, which are used when calling run_model
        When running ua the run_model attributes are the ones that are changed.
        """

        self.update_species(self.species)
        self.run_model_parameters = self.parameters
        self.y = []

    # Run the model
    def deriv(self, y, t):
        """
        deriv function called by integrate.odeint(self.deriv, y0, self.time)

        For each step when the model is run, the rate for each reaction is calculated and changes in substrates and products calculated.
        These are returned by this function as y_prime, which are added to y which is returned by run_model

        Args:
            y (list): ordered list of substrate values at this current timepoint. Has the same order as self.run_model_species_names
            t (): time, not used in this function but required for some reason

        Returns:
            y_prime - ordered list the same as y, y_prime is the new set of y's for this timepoint.
        """

        yprime = np.zeros(len(y))

        for reaction_class in self:
            yprime += reaction_class.reaction(y, self.run_model_species_names, self.run_model_parameters)

        return yprime

    def run_model(self):
        """
        Runs the model and outputs y

        Uses self.run_model_species, run_model_species_names, self.run_model_species_starting_values and self.run_model_parameters.
        These are loaded by calling self.setup_model() before running.

        Outputs saved to self.y
        """

        y0 = np.array(self.run_model_species_starting_values)
        self.y = integrate.odeint(self.deriv, y0, self.time, mxstep=self.mxsteps)
        self.reset_reaction_indexes()

        return self.y

    # Export results as dataframe and plot
    def results_dataframe(self):
        """
        Gives the results of a model run as a dataframe

        Returns:
            Pandas dataframe of results
        """
        ys_at_t = {'Time' : self.time}

        for i in range(len(self.run_model_species_names)):
            name = self.run_model_species_names[i]
            ys_at_t[name] = []

            for t in range(len(self.time)):
                ys_at_t[name].append(self.y[t][i])

        df = pd.DataFrame(ys_at_t)

        return df

    def plot_substrate(self, substrate, plot=False):
        """
        Plot a graph of substrate concentration vs time.

        Args:
            substrate (str): Name of substrate to plot
            plot (bool): Default False.  If True calls plt.show()
        """
        ys_at_t = []
        i = self.run_model_species_names.index(substrate)
        for t in range(len(self.time)):
            ys_at_t.append(self.y[t][i])

        plt.plot(self.time, ys_at_t, label=substrate)
        plt.legend()

        if plot == True:
            plt.show()

    # Check parameters when contraining parameter space
    def check_parameter_limits(self):
        all_within_limits = True
        for reaction_class in self:
            if reaction_class.sampling_limits(self.run_model_parameters) == False:
                all_within_limits = False
        return all_within_limits