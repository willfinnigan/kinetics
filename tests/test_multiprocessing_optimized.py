import time
import numpy as np
import kinetics
from kinetics.reactions.irreversible_michaelis_menton import Uni
from kinetics.sampling.scipy_sampling import ScipyDist_Sampler
from kinetics.solvers.scipy_multiprocessing_solver_optimized import SciPyMultiprocessingSolverOptimized
from scipy.stats import uniform, norm


def test_multiprocessing_methods():
    """Compare different multiprocessing methods."""
    
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
    sample_sizes = [50, 200, 1000]
    mp_methods = ['map', 'imap', 'imap_unordered', 'futures']
    
    print("Multiprocessing Method Comparison")
    print("=" * 80)
    
    for num_samples in sample_sizes:
        print(f"\nTesting {num_samples} samples:")
        print("-" * 40)
        
        # Test regular model as baseline
        sampler_regular = ScipyDist_Sampler(num_samples=num_samples)
        start_time = time.time()
        result_regular = model.run_multi({"A": 5000, "enz_1": 5}, 
                                       sampler_regular, 
                                       kinetics.SciPySolver())
        regular_time = time.time() - start_time
        
        print(f"Regular:          {regular_time:.3f}s")
        
        # Test vectorized model
        sampler_vectorized = ScipyDist_Sampler(num_samples=num_samples)
        start_time = time.time()
        result_vectorized = model.run_multi_vectorized({"A": 5000, "enz_1": 5}, 
                                                      sampler_vectorized)
        vectorized_time = time.time() - start_time
        
        print(f"Vectorized:       {vectorized_time:.3f}s ({regular_time/vectorized_time:.2f}x)")
        
        # Test different multiprocessing methods
        for mp_method in mp_methods:
            try:
                sampler_mp = ScipyDist_Sampler(num_samples=num_samples)
                
                # Create solver with specific method
                solver = SciPyMultiprocessingSolverOptimized(mp_method=mp_method)
                
                # Use run_multi_multiprocessing but with our optimized solver
                start_time = time.time()
                result_mp = model.run_multi_multiprocessing({"A": 5000, "enz_1": 5}, 
                                                          sampler_mp, solver)
                mp_time = time.time() - start_time
                
                speedup = regular_time / mp_time
                print(f"MP {mp_method:12}: {mp_time:.3f}s ({speedup:.2f}x)")
                
                # Verify results
                assert len(result_mp.multi_ys) == num_samples
                
            except Exception as e:
                print(f"MP {mp_method:12}: FAILED - {e}")
    
    # Show system info
    solver = SciPyMultiprocessingSolverOptimized()
    print(f"\nSystem info:")
    print(f"CPU cores available: {solver.n_processes}")
    print("\nTest completed!")


if __name__ == '__main__':
    test_multiprocessing_methods()