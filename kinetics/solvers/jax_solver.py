from kinetics.solvers.solver_base_class import ODESolver
import numpy as np
import pandas as pd
import jax.numpy as jnp
from diffrax import diffeqsolve, ODETerm, Dopri5, SaveAt, PIDController

class JaxSolver(ODESolver):

    # Run the model
    def deriv(self, t, y, args):
        """
        deriv function called by diffrax.diffeqsolve.

        For each step when the model is run, the rate for each reaction is calculated and changes in substrates and products calculated.
        These are returned by this function as y_prime, which are added to y which is returned by run_model

        Args:
            t: current time (required by Diffrax, not always used by underlying reaction kinetics).
            y (jax.numpy.ndarray): ordered array of substrate values at this current timepoint. Has the same order as self.run_model_species_names.
            args (tuple): Contains (run_model_species_names, run_model_parameters).

        Returns:
            y_prime (jax.numpy.ndarray): ordered array the same shape as y, y_prime is the new set of dy/dt for this timepoint.
        """
        (reactions, species_names, parameters) = args

        # y is now a JAX array, and reaction_class.reaction is expected to handle it.
        # Initialize yprime as a JAX array.
        yprime = jnp.zeros_like(y)

        for reaction_class in reactions:
            yprime += reaction_class.reaction(y, species_names)

        return yprime

    def run(self, reactions, species, parameters, time):
        """
        Runs the model and outputs y

        Uses self.run_model_species, run_model_species_names, self.run_model_species_starting_values and self.run_model_parameters.
        These are loaded by calling self.setup_model() before running.

        Outputs saved to self.y
        """

        # y0 = np.array(list(species.values()), dtype=float)
        # species_names = list(species.keys())
        # y = integrate.odeint(self.deriv, y0, time, args=(reactions, species_names, parameters), mxstep=self.mxsteps)

        y0 = jnp.array(list(species.values()))
        time = jnp.asarray(time)  # Ensure time is a JAX array for compatibility with Diffrax

        # Package arguments for the derivative function
        species_names = list(species.keys())
        deriv_args = (reactions, species_names, parameters)

        term = ODETerm(self.deriv)
        solver = Dopri5()

        saveat = SaveAt(ts=time)

        # Using default rtol/atol for PIDController, can be tuned if needed.
        stepsize_controller = PIDController(rtol=1e-7, atol=1e-7)

        # Determine initial dt0; a small fraction of the first time step or a fixed small value.
        # Ensure dt0 is a concrete Python float.
        if len(time) > 1:
            # Convert JAX array elements to NumPy/Python float for this calculation
            time_np = np.asarray(time)
            dt_val = (time_np[1] - time_np[0]) / 10.0
            dt0 = float(dt_val)
        else:
            # Default dt0 if only one time point or for safety.
            dt0 = 0.1 # This is already a Python float.

        solution = diffeqsolve(term,
                               solver,
                               t0=time[0],
                               t1=time[-1],
                               dt0=dt0,
                               y0=y0,
                               args=deriv_args,
                               saveat=saveat,
                               stepsize_controller=stepsize_controller)

        # Convert JAX array solution to NumPy array for compatibility with existing methods
        y = np.asarray(solution.ys)

        return y
