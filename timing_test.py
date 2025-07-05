import time
import kinetics

def time_solver(solver_mode, runs=3):
    """Time multiple runs of the solver"""
    times = []
    
    for i in range(runs):
        start = time.time()
        
        model = kinetics.Model()
        model.set_time(0, 1000, 100)

        enzyme_1 = kinetics.Uni(kcat='enz1_kcat', kma='enz1_km', enz='enz_1', a='A',
                                substrates=['A'], products=['B'])

        enzyme_1.parameters = {'enz1_kcat': 100,
                               'enz1_km': 10000}

        model.add_reaction(enzyme_1)
        model.set_species({"A": 10000, "enz_1": 5})

        result = model.run_model(mode=solver_mode)
        
        end = time.time()
        times.append(end - start)
        print(f"Run {i+1}: {(end - start)*1000:.1f}ms")
    
    avg_time = sum(times) / len(times)
    print(f"\n{solver_mode} average: {avg_time*1000:.1f}ms")
    return times

if __name__ == "__main__":
    print("=== SciPy Solver ===")
    scipy_times = time_solver('scipy')
    
    print("\n=== JAX Solver ===")
    jax_times = time_solver('jax')
    
    print(f"\nSpeed ratio: JAX is {jax_times[0]/scipy_times[0]:.1f}x slower on first run")
    if len(jax_times) > 1:
        print(f"Speed ratio: JAX is {jax_times[1]/scipy_times[1]:.1f}x slower on second run")