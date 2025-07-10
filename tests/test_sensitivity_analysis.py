import pytest
import kinetics
import numpy as np
from numpy.testing import assert_allclose

from kinetics import (SalibSaltelliSampler, 
                      analyze_sensitivity_at_timepoint, 
                      analyze_sensitivity_time_to_concentration)


@pytest.fixture
def two_enzyme_model_result():
    """Create a simple 2-enzyme model and run it with Saltelli sampler for sensitivity analysis"""
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
    
    solver = kinetics.JaxSolver()
    sampler = SalibSaltelliSampler(num_samples=100, log_parameters=['enz2_km'])
    
    result = model.run_multi({"A": 10000, "enz_1": 5, "enz_2": 2},
                             sampler,
                             solver)
    
    return result, sampler


def test_analyze_sensitivity_at_timepoint(two_enzyme_model_result):
    """Test sensitivity analysis for concentration at a specific timepoint"""
    result, sampler = two_enzyme_model_result
    
    # Test the new simplified API
    sa_result = analyze_sensitivity_at_timepoint(
        result=result,
        sampler=sampler,
        species='B',
        timepoint=500,
        threshold=0.01
    )
    
    # Check that we get a result object with expected properties
    assert hasattr(sa_result, 'data')
    assert hasattr(sa_result, 'plot')
    
    # Check that data is a DataFrame with expected columns
    assert 'S1' in sa_result.data.columns
    assert 'ST' in sa_result.data.columns
    
    # Check that filtering worked (no values below threshold)
    assert (sa_result.data['ST'] >= 0.01).all()
    
    # Check that plot method works
    sa_result.plot()


def test_analyze_sensitivity_time_to_concentration(two_enzyme_model_result):
    """Test sensitivity analysis for time to reach a concentration"""
    result, sampler = two_enzyme_model_result
    
    # Test the new simplified API
    sa_result = analyze_sensitivity_time_to_concentration(
        result=result,
        sampler=sampler,
        species='B',
        concentration=5000,
        mode='>=',
        threshold=0.01
    )
    
    # Check that we get a result object with expected properties
    assert hasattr(sa_result, 'data')
    assert hasattr(sa_result, 'plot')
    
    # Check that data is a DataFrame with expected columns
    assert 'S1' in sa_result.data.columns
    assert 'ST' in sa_result.data.columns
    
    # Check that filtering worked (no values below threshold)
    assert (sa_result.data['ST'] >= 0.01).all()
    
    # Check that plot method works
    sa_result.plot()


def test_analyze_sensitivity_time_to_concentration_less_than(two_enzyme_model_result):
    """Test sensitivity analysis for time to reach a concentration with <= mode"""
    result, sampler = two_enzyme_model_result
    
    # Test with <= mode
    sa_result = analyze_sensitivity_time_to_concentration(
        result=result,
        sampler=sampler,
        species='A',  # Use A which should decrease over time
        concentration=5000,
        mode='<=',
        threshold=0.01
    )
    
    # Check that we get a result object
    assert hasattr(sa_result, 'data')
    assert hasattr(sa_result, 'plot')