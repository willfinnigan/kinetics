import kinetics

# Define reactions
enzyme_1 = kinetics.Uni(kcat='enz1_kcat', kma='enz1_km', enz='enz_1', a='A',
                            substrates=['A'], products=['B'])

enzyme_1.parameters = {'enz1_kcat': 100,
                       'enz1_km': 10000}

enzyme_2 = kinetics.Uni(kcat='enz2_kcat', kma='enz2_km', enz='enz_2', a='B',
                        substrates=['B'], products=['C'])

enzyme_2.parameters = {'enz2_kcat': 100,
                       'enz2_km': 10000}

# Set up the model
model = kinetics.Model()
model.set_time(0, 1000, 100)
model.add_reaction(enzyme_1)
model.add_reaction(enzyme_2)

solver = kinetics.SciPySolver()

result = model.run_single({"A": 10000, "enz_1": 5},
                         solver)