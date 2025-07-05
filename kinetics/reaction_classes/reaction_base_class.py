import numpy as np
import copy
try:
    import jax.numpy as jnp
    HAS_JAX = True
except ImportError:
    jnp = None
    HAS_JAX = False

def calculate_yprime(y, rate, substrates, products, substrate_names):
    """
    This function is used by the rate classes the user creates.

    It takes the numpy array for y_prime,
    and adds or subtracts the amount in rate to all the substrates or products listed
    Returns the new y_prime

    Args:
        y: a numpy array for the substrate values, the same order as y
        rate: the rate calculated by the user made rate equation
        substrates: list of substrates for which rate should be subtracted
        products: list of products for which rate should be added
        substrate_names: the ordered list of substrate names in the model.  Used to get the position of each substrate or product in y_prime

    Returns:
        y_prime: following the addition or subtraction of rate to the specificed substrates
    """

    # Create zeros array compatible with input array type (JAX or NumPy)
    if HAS_JAX and hasattr(y, 'shape') and str(type(y)).startswith('<class \'jax'):
        y_prime = jnp.zeros_like(y)
    else:
        y_prime = np.zeros(len(y))

    for name in substrates:
        idx = substrate_names.index(name)
        if HAS_JAX and hasattr(y, 'shape') and str(type(y)).startswith('<class \'jax'):
            y_prime = y_prime.at[idx].add(-rate)
        else:
            y_prime[idx] -= rate

    for name in products:
        idx = substrate_names.index(name)
        if HAS_JAX and hasattr(y, 'shape') and str(type(y)).startswith('<class \'jax'):
            y_prime = y_prime.at[idx].add(rate)
        else:
            y_prime[idx] += rate

    return y_prime

def check_positive(y_prime):
    """
    Chack that substrate values are not negative when they shouldnt be
    """

    for i in range(len(y_prime)):
        if y_prime[i] < 0:
            y_prime[i] = 0

    return y_prime

class Reaction():

    def __init__(self):

        # These are set by the user
        self.parameters = {}
        self.parameter_distributions = {}

        # indexes used to access values during model run
        self.substrate_indexes = []
        self.parameter_indexes = []

        # These are set when the reaction is set up
        self.reaction_substrate_names = []
        self.parameter_names = []
        self.substrates = []
        self.products = []

        # These are added as needed
        self.modifiers = []
        self.check_positive = False
        self.check_limits_functions = []

    def set_parameter_defaults_to_mean(self):
        for name in self.parameter_distributions:
            if name not in self.parameters:
                if type(self.parameter_distributions[name]) == list or type(self.parameter_distributions[name]) == tuple:
                    self.parameters[name] = (self.parameter_distributions[name][0] + self.parameter_distributions[name][1]) / 2
                else:
                    self.parameters[name] = self.parameter_distributions[name].mean()

    def setup_reaction(self, species_names, parameter_names):

        # get indexes
        self.substrate_indexes = []
        for name in self.reaction_substrate_names:
            self.substrate_indexes.append(species_names.index(name))
        self.parameter_indexes = []
        for name in self.parameter_names:
            self.parameter_indexes.append(parameter_names.index(name))

        # set up modifiers
        for modifier in self.modifiers:
            modifier.get_substrate_indexes(self.reaction_substrate_names)
            modifier.get_parameter_indexes(self.parameter_names)

    def add_modifier(self, modifier):
        for name in modifier.parameter_names:
            if name not in self.parameter_names:
                self.parameter_names.append(name)

        for name in modifier.substrate_names:
            if name not in self.reaction_substrate_names:
                self.reaction_substrate_names.append(name)

        self.modifiers.append(modifier)

    def calculate_rate(self, substrates, parameters):
        return 0

    def reaction(self, y, substrate_names, parameter_values):
        """Calculate the rate of the reaction and return the change in substrate concentrations (y_prime)."""

        # Get the substrates from y using the substrate indexes
        substrates = []
        for index in self.substrate_indexes:
            substrates.append(y[index])

        # Get the parameters using the parameter indexes
        parameters = []
        for index in self.parameter_indexes:
            parameters.append(parameter_values[index])

        # calculate the effects of any modifiers
        for modifier in self.modifiers:
            substrates, parameters = modifier.calc_modifier(substrates, parameters)

        # calculate the rate (this function is modified by the user)
        rate = self.calculate_rate(substrates, parameters)

        # calculate the change in substrate concentrations (y_prime)
        y_prime = calculate_yprime(y, rate, self.substrates, self.products, substrate_names)
        y_prime = self.modify_product(y_prime, substrate_names)

        if self.check_positive == True:
            y_prime = check_positive(y_prime)

        return y_prime

    def modify_product(self, y_prime, substrate_names):
        return y_prime

    def sampling_limits(self, parameter_dict):
        # Return true if parameters within limits, false if not
        for func in self.check_limits_functions:
            if func(parameter_dict) == False:
                return False

        return True
