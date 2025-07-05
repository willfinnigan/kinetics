# from __future__ import division

import pytest
import kinetics
import numpy as np
from numpy.testing import assert_allclose

@pytest.mark.parametrize("solver_mode", ['jax', 'scipy', 'jax_gpu'])
def test_simple_model(solver_mode):
    model = kinetics.Model()
    model.set_time(0, 1000, 100)

    enzyme_1 = kinetics.Uni(kcat='enz1_kcat', kma='enz1_km', enz='enz_1', a='A',
                            substrates=['A'], products=['B'])

    enzyme_1.parameters = {'enz1_kcat': 100,
                           'enz1_km': 10000}

    model.add_reaction(enzyme_1)

    model.set_species({"A": 10000, "enz_1": 5})

    result = model.run_model(mode=solver_mode)
    df = result.results_dataframe()

    start = df['A'][0]
    end = df['B'][99]

    expected = np.array([10000.0, 10000.0])
    actual = np.array([start, end])

    assert_allclose(expected, actual, atol=1, rtol=1)





