import kinetics
from scipy.stats import reciprocal, uniform, norm


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

solver = kinetics.SciPySolver()
sampler = kinetics.ScipyDist_Sampler(num_samples=1000)

result = model.run_multi({"A": 10000, "enz_1": 5},
                         sampler,
                         solver)