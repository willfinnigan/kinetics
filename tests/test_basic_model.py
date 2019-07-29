# from __future__ import division

import kinetics
import numpy as np
from numpy.testing import assert_allclose

def test_simple_model():
    model = kinetics.Model()
    model.set_time(0, 1000, 100)

    enzyme_1 = kinetics.Uni(kcat='enz1_kcat', kma='enz1_km', enz='enz_1', a='A',
                            substrates=['A'], products=['B'])

    enzyme_1.parameters = {'enz1_kcat': 100,
                           'enz1_km': 10000}

    model.append(enzyme_1)

    model.species = {"A": 10000,
                     "enz_1": 5}
    model.setup_model()

    model.run_model()
    df = model.results_dataframe()

    start = df['A'][0]
    end = df['B'][99]

    expected = np.array([10000.0,10000.0])
    actual = np.array([start, end])

    assert_allclose(expected, actual, atol=1, rtol=1)





