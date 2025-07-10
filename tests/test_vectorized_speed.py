import time
import numpy as np
import kinetics
from kinetics.reactions.irreversible_michaelis_menton import Uni
from kinetics.sampling.scipy_sampling import ScipyDist_Sampler
from scipy.stats import uniform, norm


def test_speed_comparison():
    """Compare speed of vectorized vs non-vectorized model runs."""
    
    # Create reaction with parameter distributions
    enzyme_1 = Uni(kcat='enz1_kcat', kma='enz1_km', enz='enz_1', a='A',
                   substrates=['A'], products=['B'])
    
    enzyme_1.parameter_distributions = {'enz1_kcat': norm(100, 12),
                                       'enz1_km': uniform(2000, 6000)}
    
    # Set up model
    model = kinetics.Model()
    model.set_time(0, 1000, 100)
    model.add_reaction(enzyme_1)
    
    # Test with different sample sizes
    sample_sizes = [10, 50, 100, 500, 1000, 5000]
    
    print("Speed Comparison: Vectorized vs Non-vectorized")
    print("=" * 60)
    print(f"{'Samples':<10} {'Regular (s)':<15} {'Vectorized (s)':<15} {'Speedup':<10}")
    print("-" * 60)
    
    for num_samples in sample_sizes:
        # Test regular model
        sampler_regular = ScipyDist_Sampler(num_samples=num_samples)
        
        start_time = time.time()
        result_regular = model.run_multi({"A": 5000, "enz_1": 5}, 
                                       sampler_regular, 
                                       kinetics.SciPySolver())
        regular_time = time.time() - start_time
        
        # Test vectorized model
        sampler_vectorized = ScipyDist_Sampler(num_samples=num_samples)
        
        start_time = time.time()
        result_vectorized = model.run_multi_vectorized({"A": 5000, "enz_1": 5}, 
                                                      sampler_vectorized)
        vectorized_time = time.time() - start_time
        
        # Calculate speedup
        speedup = regular_time / vectorized_time if vectorized_time > 0 else float('inf')
        
        print(f"{num_samples:<10} {regular_time:<15.3f} {vectorized_time:<15.3f} {speedup:<10.2f}x")
        
        # Verify results are the same shape
        assert len(result_regular.multi_ys) == len(result_vectorized.multi_ys) == num_samples
        assert result_regular.multi_ys[0].shape == result_vectorized.multi_ys[0].shape
    
    print("\nTest completed successfully!")


if __name__ == '__main__':
    test_speed_comparison()