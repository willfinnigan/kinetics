# from __future__ import division

import pytest
import kinetics
import numpy as np
from numpy.testing import assert_allclose

from scipy.stats import reciprocal, uniform, norm

from kinetics.sampling.salib_sampling import SalibLatinHypercubeSampler, SalibSaltelliSampler
from kinetics.sampling.scipy_sampling import ScipyDist_Sampler

solvers = {'scipy': kinetics.SciPySolver(),
           'jax': kinetics.JaxSolver()}

@pytest.mark.parametrize("solver_mode", ['jax', 'scipy'])
def test_distribution_two_enzyme_model(solver_mode):
    # Define reactions
    enzyme_1 = kinetics.Uni(kcat='enz1_kcat', kma='enz1_km', enz='enz_1', a='A',
                            substrates=['A'], products=['B'])

    enzyme_1.parameter_distributions = {'enz1_kcat': norm(100, 12),
                                        'enz1_km': uniform(2000, 6000)}

    enzyme_2 = kinetics.Uni(kcat='enz2_kcat', kma='enz2_km', enz='enz_2', a='B',
                            substrates=['B'], products=['C'])

    enzyme_2.parameter_distributions = {'enz2_kcat': norm(30, 5),
                                        'enz2_km': reciprocal(1, 10000)}

    # Set up the model
    model = kinetics.Model()
    model.set_time(0, 1000, 100)
    model.add_reaction(enzyme_1)
    model.add_reaction(enzyme_2)

    solver = solvers[solver_mode]
    sampler = ScipyDist_Sampler(num_samples=1000)

    result = model.run_multi({"A": 10000, "enz_1": 5},
                             sampler,
                             solver)


samplers = {'latin': SalibLatinHypercubeSampler,
            'satelli': SalibSaltelliSampler}

@pytest.mark.parametrize("sampler_mode", ['latin', 'satelli'])
def test_salib_two_enzyme_model(sampler_mode):
    # Define reactions
    enzyme_1 = kinetics.Uni(kcat='enz1_kcat', kma='enz1_km', enz='enz_1', a='A',
                            substrates=['A'], products=['B'])

    enzyme_1.parameter_distributions = {'enz1_kcat': (90, 110),
                                        'enz1_km': (2000, 6000)}

    enzyme_2 = kinetics.Uni(kcat='enz2_kcat', kma='enz2_km', enz='enz_2', a='B',
                            substrates=['B'], products=['C'])

    enzyme_2.parameter_distributions = {'enz2_kcat': (25, 35),
                                        'enz2_km': (1, 10000)}

    # Set up the model
    model = kinetics.Model()
    model.set_time(0, 1000, 100)
    model.add_reaction(enzyme_1)
    model.add_reaction(enzyme_2)

    solver = solvers['scipy']

    sampler = samplers[sampler_mode](num_samples=100,
                                     log_parameters=['enz2_km'])

    result = model.run_multi({"A": 10000, "enz_1": 5},
                             sampler,
                             solver)

