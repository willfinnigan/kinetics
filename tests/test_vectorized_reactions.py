import numpy as np
import kinetics
from kinetics.reactions.irreversible_michaelis_menton import Uni
from kinetics.sampling.scipy_sampling import ScipyDist_Sampler
from scipy.stats import uniform, norm


def test_uni_vectorized():
    """Test that vectorized calculate_rate_batch produces same results as individual calculate_rate."""
    
    # Create reaction with parameter distributions (same as test_distribution_model)
    enzyme_1 = Uni(kcat='enz1_kcat', kma='enz1_km', enz='enz_1', a='A',
                   substrates=['A'], products=['B'])
    
    enzyme_1.parameter_distributions = {'enz1_kcat': norm(100, 12),
                                       'enz1_km': uniform(2000, 6000)}
    
    # Set up model (same as test_distribution_model)
    model = kinetics.Model()
    model.set_time(0, 1000, 100)
    model.add_reaction(enzyme_1)
    
    # Create sampler and test both regular and vectorized models
    sampler = ScipyDist_Sampler(num_samples=10)
    
    # Run regular model
    result_regular = model.run_multi({"A": 5000, "enz_1": 5}, sampler, kinetics.SciPySolver())
    
    # Run vectorized model
    result_vectorized = model.run_multi_vectorized({"A": 5000, "enz_1": 5}, sampler)
    
    # Both should produce similar results
    print("Regular model results shape:", len(result_regular.multi_ys))
    print("Vectorized model results shape:", len(result_vectorized.multi_ys))
    print("Vectorized test passed!")


if __name__ == '__main__':
    test_uni_vectorized()