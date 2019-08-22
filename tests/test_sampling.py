

import kinetics
from scipy.stats import multivariate_normal, norm



# Define a multivariate_normal distribution, sample from it,  be able to run a model.


def test_multivariate():
    enzyme_1 = kinetics.Uni(kcat='enz1_kcat', kma='enz1_km', enz='enz_1', a='A',
                            substrates=['A'], products=['B'])

    enzyme_1.parameter_distributions = {'enz1_kcat': multivariate_normal([98, 114], [[46, 9],[9.19, 3.4]]),
                                        'enz1_km': ('enz1_kcat', 1)}

    enzyme_1.parameters = {'enz1_kcat': 98,
                           'enz1_km': 114}

    # Set up the model
    model = kinetics.Model(logging=False)
    model.append(enzyme_1)
    model.set_time(0, 120, 1000)

    # Set starting concentrations
    model.species = {"A": 10000}
    model.species_distributions = {"enz_1": norm(4, 4 * 0.05),
                                   "enz_2": norm(10, 10 * 0.05)}
    model.setup_model()

    # Run the model 1000 times, sampling from distributions
    samples = kinetics.sample_distributions(model, num_samples=50)
    outputs = kinetics.run_all_models(model, samples, logging=True)