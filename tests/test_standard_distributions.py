

import kinetics
import kinetics.Uncertainty as ua
from scipy.stats import reciprocal, uniform, norm

def test_distributions_model():
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
    model = kinetics.Model(logging=False)
    model.append(enzyme_1)
    model.append(enzyme_2)
    model.set_time(0, 120, 1000)

    # Set starting concentrations
    model.species = {"A": 10000}
    model.species_distributions = {"enz_1": norm(4, 4 * 0.05),
                                   "enz_2": norm(10, 10 * 0.05)}
    model.setup_model()

    # Run the model 1000 times, sampling from distributions
    samples = ua.sample_distributions(model, num_samples=1000)
    outputs = ua.run_all_models(model, samples, logging=True)

    model.run_model()
    model.plot_substrate('A')
    model.plot_substrate('B')
    model.plot_substrate('C')

    # Plot model runs at 95% CI
    ci_dataframes = ua.dataframes_quartiles(model, outputs)
    ua.plot_ci_intervals(['A', 'B', 'C'], ci_dataframes, colours=['blue', 'darkorange', 'green'])

    # Plot all model runs
    all_runs_dataframes = ua.dataframes_all_runs(model, outputs)
    ua.plot_substrate('A', all_runs_dataframes, colour='blue', alpha=0.01, linewidth=5)
    ua.plot_substrate('B', all_runs_dataframes, colour='darkorange', alpha=0.01, linewidth=5)
    ua.plot_substrate('C', all_runs_dataframes, colour='green', alpha=0.01, linewidth=5)