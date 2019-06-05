import numpy as np
import pandas as pd
from scipy import integrate
import copy
import time
from kinetics.Uncertainty_Sensitivity_Analysis import get_bounds_from_pc_error
import matplotlib.pyplot as plt


class Model(list):
    """
    The model class inherits from a list.
      The user made reaction classes are appended to this self list.

    When the model is run it uses only:
      self.parameters
      self.species_names
      self.species_starting_values

    When the model is run, it calls integrate.odeint(self.deriv, y0, self.time).
      self.deriv(self, y, t) runs each reaction_class.reaction(y, self.species_names, self.parameters) in turn,
        with the output added to y_prime as the relevent index (determined by self.species_names)

    Each reaction class contains parameter_defaults and parameter_bounds.
      These are used to set self.parameter_defaults and self.parameter_bounds, when self.set_parameters_from_reactions() is called.

    Species are set using the set_species_defaults and set_species_bounds functions.

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
        self.species_defaults = {}
        self.species_distributions = {}

        """ Parameters - used to reset the model, or as the bounds to run ua/sa.  Set by self.set_parameters_from_reactions() """
        self.parameter_defaults = {}
        self.parameter_distributions = {}

        """ Species and parameters used when the model is ran. These are changed each run when doing ua/sa """
        self.species = {}
        self.species_names = []
        self.species_starting_values = []
        self.parameters = {}

        self.y = []

        self.logging = logging

    # Time
    def set_time(self, start, end, steps):
        """
        This function sets all the time parameters for the model.

        :param start: integer - the start time - usually 0
        :param end: integer - the end time
        :param steps: integer - the number of timepoints for the output
        :param mxsteps: integer - the max number of steps that should be used to solve the odes
        """
        self.start = start
        self.end = end
        self.steps = steps
        self.time = np.linspace(self.start, self.end, self.steps)

    # Parameters
    def set_parameters_from_reactions(self):
        """
        Sets all the parameter variables from those set in the reaction classes

        For each reaction_class, updates self.parameters, self_parameter_defaults and self.parameter_bounds,
        with the dictionaries held in each reaction_class.
        This will add new keys, or overwrite existing ones.
        """

        self.parameters = {}
        self.parameter_defaults = {}
        self.parameter_distributions = {}

        for reaction_class in self:
            if reaction_class.parameter_defaults == {}:
                reaction_class.set_parameter_defaults_to_means()
                if self.logging==True:
                    print('Set parameters to means:')
                    print(reaction_class.parameter_defaults)

            self.parameters.update(reaction_class.parameter_defaults)
            self.parameter_defaults.update(copy.deepcopy(reaction_class.parameter_defaults))
            self.parameter_distributions.update(reaction_class.parameter_distributions)

    # Species
    def update_species(self, species_dict):
        self.species.update(species_dict)
        self.species_names = list(self.species.keys())
        self.species_starting_values = list(self.species.values())

    def load_species_from_reactions(self):
        for reaction in self:
            for substrate in reaction.substrates + reaction.products + reaction.reaction_substrate_names:
                if substrate not in self.species_defaults:
                    self.species_defaults[substrate] = 0
                    if self.logging == True:
                        print("Loaded species '" + str(substrate) + "' as 0")

    def set_species_defaults_to_mean(self):
        for name in self.species_distributions:
            self.species_defaults[name] = self.species_distributions[name].mean()

        return self.species_defaults

    def setup_model(self):
        self.set_species_defaults_to_mean()
        self.load_species_from_reactions()
        self.set_parameters_from_reactions()
        self.update_species(self.species_defaults)

    # Run the model
    def deriv(self, y, t):
        """
        deriv function called by integrate.odeint(self.deriv, y0, self.time)

        For each step when the model is run, the rate for each reaction is calculated and changes in substrates and products calculated.
        These are returned by this function as y_prime, which are added to y which is returned by run_model

        :param y: ordered list of substrate values at this current timepoint. Has the same order as self.substrate_names
        :param t: time, not used in this function but required for some reason
        :return: y_prime - ordered list the same as y, y_prime is the new set of y's for this timepoint.
        """

        yprime = np.zeros(len(y))

        for reaction_class in self:
            yprime += reaction_class.reaction(y, self.species_names, self.parameters)

        return yprime

    def run_model(self):
        """
        Runs the model and outputs y

        This will use self.species_names, self.species_starting_values and self.parameters to run the model.

        Outputs y which is a numpy array of 2 dimensions.
          The first dimension gives a list of all the substrate concentrations at that timepoint.
          The first dimension is the same size as self.time.
          Each index in self.time relates to an index in the first dimension of y.

          eg.  y[0] will return a list of all the starting substrate concentrations
               y[1] gives the substrate concentrations at the first timepoints

               y[0][2] gives the substrate concentration of the second substrate at the first timepoint.

        :return: y - a numpy array of 2 dimensions. Time by substrate.
        """
        y0 = np.array(self.species_starting_values)
        self.y = integrate.odeint(self.deriv, y0, self.time, mxstep=self.mxsteps)
        self.reset_reactions()

        return self.y

    # Reset the model back to default settings
    def reset_reactions(self):
        for reaction_class in self:
            reaction_class.reset_reaction()

    def reset_model(self):
        """
        Reset the model back to the default settings

        This uses self.species_defaults and self.parameter_defaults
          to set self.species_names, self.species_starting_values and self.parameters
          back to the original default settings.  These the variables used to run the model.
        """

        self.species = self.species_defaults
        self.species_names, self.species_starting_values = get_species_positions(self.species)
        self.parameters = self.parameter_defaults
        self.y = []

    # Export results as dataframe
    def results_dataframe(self):
        ys_at_t = {'Time' : self.time}

        for i in range(len(self.species_names)):
            name = self.species_names[i]
            ys_at_t[name] = []

            for t in range(len(self.time)):
                ys_at_t[name].append(self.y[t][i])

        df = pd.DataFrame(ys_at_t)

        return df

    # Check parameters
    def check_parameter_limits(self):
        all_within_limits = True
        for reaction_class in self:
            if reaction_class.sampling_limits(self.parameters) == False:
                all_within_limits = False
        return all_within_limits

    def plot_substrate(self, substrate):
        ys_at_t = []
        i = self.species_names.index(substrate)
        for t in range(len(self.time)):
            ys_at_t.append(self.y[t][i])

        plt.plot(self.time, ys_at_t)


"""Functions to add or substract the rate from yprime at the correct index's"""
def calculate_yprime(y, rate, substrates, products, substrate_names):
    """
    This function is used by the rate classes the user creates.

    It takes the numpy array for y_prime,
      and adds or subtracts the amount in rate to all the substrates or products listed
    Returns the new y_prime

    :param y_prime: a numpy array for the substrate values, the same order as y
    :param rate:   the rate calculated by the user made rate equation
    :param substrates: list of substrates for which rate should be subtracted
    :param products: list of products for which rate should be added
    :param substrate_names: the ordered list of substrate names in the model.  Used to get the position of each substrate or product in y_prime
    :return: y_prime: following the addition or subtraction of rate to the specificed substrates
    """

    y_prime = np.zeros(len(y))

    for name in substrates:
        y_prime[substrate_names.index(name)] -= rate

    for name in products:
        y_prime[substrate_names.index(name)] += rate

    return y_prime

def check_positive(y_prime):

    for i in range(len(y_prime)):
        if y_prime[i] < 0:
            y_prime[i] = 0

    return y_prime

