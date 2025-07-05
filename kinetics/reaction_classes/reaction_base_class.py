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

        self.reaction_substrate_names = []
        self.substrate_indexes = []
        self.substrates = []
        self.products = []
        self.parameters = {}
        self.parameter_distributions = {}
        self.parameter_names = []
        self.run_model_parameters = []
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

    def get_indexes(self, substrate_names):
        self.substrate_indexes = []
        for name in self.reaction_substrate_names:
            self.substrate_indexes.append(substrate_names.index(name))

    def get_parameters(self, parameter_dict):
        """Look up the parameters in the parameter_dict and return them in the order of self.parameter_names."""
        self.run_model_parameters = []
        for name in self.parameter_names:
            self.run_model_parameters.append(parameter_dict[name])


    def reset_reaction(self):
        self.substrate_indexes = []
        self.run_model_parameters = []

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

    def reaction(self, y, substrate_names):
        """Calculate the rate of the reaction and return the change in substrate concentrations (y_prime)."""

        # Get the substrates from y using the substrate indexes
        substrates = []
        for index in self.substrate_indexes:
            substrates.append(y[index])

        parameters = self.run_model_parameters

        # calculate the effects of any modifiers
        for modifier in self.modifiers:
            substrates, parameters = modifier.calc_modifier(substrates, parameters)

        rate = self.calculate_rate(substrates, parameters)

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
