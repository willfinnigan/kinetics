import numpy as np
from scipy import integrate

from kinetics.solvers.solver_interface import ODESolver


class SciPySolver(ODESolver):

    def __init__(self, mxsteps=5000):
        self.mxsteps = mxsteps

        # Run the model
    def deriv(self, y, t, reactions, species_names, parameter_values):
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

        for reaction_class in reactions:
            yprime += reaction_class.reaction(y, species_names, parameter_values)

        return yprime

    def run(self, reactions, species_names, species_values, parameter_values, time):
        """
        Runs the model and outputs y

        Uses self.run_model_species, run_model_species_names, self.run_model_species_starting_values and self.run_model_parameters.
        These are loaded by calling self.setup_model() before running.

        Outputs saved to self.y
        """


        y0 = np.array(species_values, dtype=float)
        y = integrate.odeint(self.deriv, y0, time, args=(reactions, species_names, parameter_values), mxstep=self.mxsteps)
        return y


