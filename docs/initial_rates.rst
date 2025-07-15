Initial Rates Analysis
======================

The initial rates analysis module provides tools for calculating initial reaction rates from kinetic models, both for single parameter values and uncertainty analysis with parameter distributions.

Overview
--------

Initial rates analysis is a fundamental technique in enzyme kinetics where reaction rates are measured at the beginning of the reaction when substrate concentrations are still high and product inhibition is minimal. This module provides functions to:

* Calculate initial rates for single parameter values
* Perform uncertainty analysis with parameter distributions
* Generate substrate concentration series around Km values
* Fit Michaelis-Menten equations to experimental data
* Visualize results with confidence intervals

Key Functions
-------------

Single Parameter Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: kinetics.analysis.initial_rates.calc_initial_rates_single

This function calculates initial rates for a series of substrate concentrations using fixed parameter values.

**Example:**

.. code-block:: python

    import kinetics
    from kinetics.analysis.initial_rates import calc_initial_rates_single

    # Create a simple enzyme model
    model = kinetics.Model()
    enzyme = kinetics.Uni(kcat='kcat', kma='km', enz='enzyme', a='substrate',
                          substrates=['substrate'], products=['product'])
    enzyme.parameters = {'kcat': 100.0, 'km': 1000.0}
    model.add_reaction(enzyme)

    # Calculate initial rates
    substrate_concs = [100, 500, 1000, 2000, 5000]
    starting_concentrations = {'enzyme': 1.0}
    
    rates = calc_initial_rates_single(
        model=model,
        substrate_name='substrate',
        enzyme_name='enzyme',
        substrate_concs=substrate_concs,
        starting_concentrations=starting_concentrations,
        time=1.0
    )

Uncertainty Analysis
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: kinetics.analysis.initial_rates.calc_initial_rates_multi

This function performs Monte Carlo simulation to estimate initial rates and their uncertainty.

**Example:**

.. code-block:: python

    from kinetics.analysis.initial_rates import calc_initial_rates_multi
    from kinetics.sampling.scipy_sampling import ScipyDist_Sampler
    from scipy.stats import norm

    # Add parameter distributions to the model
    enzyme.parameter_distributions = {
        'kcat': norm(loc=100.0, scale=10.0),
        'km': norm(loc=1000.0, scale=100.0)
    }

    # Create sampler
    sampler = ScipyDist_Sampler(num_samples=1000)

    # Calculate initial rates with uncertainty
    rate_quartiles, rate_all = calc_initial_rates_multi(
        model=model,
        substrate_name='substrate',
        enzyme_name='enzyme',
        substrate_concs=substrate_concs,
        starting_concentrations=starting_concentrations,
        sampler=sampler,
        quartile_range=95
    )

Utility Functions
-----------------

Unit Conversion
~~~~~~~~~~~~~~~

.. autofunction:: kinetics.analysis.initial_rates.kcat_to_umolminmg

Converts kcat units from μM/min/μM_enzyme to μmol/min/mg_enzyme.

Concentration Series Generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: kinetics.analysis.initial_rates.concentrations_around_km

Generates substrate concentrations as multiples of the Km value for systematic kinetic analysis.

**Example:**

.. code-block:: python

    from kinetics.analysis.initial_rates import concentrations_around_km

    # Generate concentrations around Km
    concentrations = concentrations_around_km(
        reaction=enzyme,
        km_param_name='km',
        datapoints=(0, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)
    )

Curve Fitting
~~~~~~~~~~~~~

.. autofunction:: kinetics.analysis.initial_rates.fit_mm

Fits Michaelis-Menten equation to experimental data.

.. autofunction:: kinetics.analysis.initial_rates.standard_mm_equation

Standard Michaelis-Menten equation implementation.

**Example:**

.. code-block:: python

    from kinetics.analysis.initial_rates import fit_mm
    import numpy as np

    # Fit Michaelis-Menten curve to data
    x_data = np.array([100, 500, 1000, 2000, 5000])
    y_data = np.array([8.5, 27.5, 40.0, 60.0, 75.0])
    
    fit_result = fit_mm(x_data, y_data, verbose=True)
    print(f"Km = {fit_result['Km'][0]} ± {fit_result['Km'][1]}")
    print(f"Kcat = {fit_result['Kcat'][0]} ± {fit_result['Kcat'][1]}")

Visualization
-------------

.. autofunction:: kinetics.analysis.initial_rates.plot_scatter_all_runs

Plots scatter points for all simulation runs from uncertainty analysis.

.. autofunction:: kinetics.analysis.initial_rates.plot_fit_ua

Plots fitted curves with confidence intervals for uncertainty analysis results.

**Example:**

.. code-block:: python

    from kinetics.analysis.initial_rates import plot_scatter_all_runs, plot_fit_ua
    import matplotlib.pyplot as plt

    # Plot uncertainty analysis results
    plt.figure(figsize=(10, 6))
    
    # Plot all individual runs as scatter points
    plot_scatter_all_runs((rate_quartiles, rate_all), colour='gray', alpha=0.3)
    
    # Plot fitted curves with confidence intervals
    plot_fit_ua((rate_quartiles, rate_all), substrate_concs, colour='blue')
    
    plt.xlabel('Substrate Concentration (μM)')
    plt.ylabel('Initial Rate (μM/min/μM_enzyme)')
    plt.title('Initial Rates Analysis with Uncertainty')
    plt.show()

Best Practices
--------------

1. **Time Point Selection**: Choose a time point that is early enough to avoid product inhibition but late enough to provide measurable substrate depletion.

2. **Substrate Concentration Range**: Use a wide range of substrate concentrations, including some well below and above the expected Km value.

3. **Parameter Distributions**: When defining parameter distributions for uncertainty analysis, use realistic estimates based on experimental data or literature values.

4. **Sample Size**: Use sufficient Monte Carlo samples (typically 500-1000) for stable uncertainty estimates.

5. **Enzyme Concentration**: Ensure enzyme concentrations are much lower than substrate concentrations to maintain pseudo-steady-state conditions.

Error Handling
--------------

The module includes comprehensive error handling for common issues:

* Missing enzyme in starting concentrations
* Invalid parameter names in reactions
* Curve fitting failures
* Missing parameter distributions

All functions provide informative error messages to help diagnose and fix issues.

Integration with v2 Architecture
---------------------------------

The initial rates analysis module is fully compatible with the v2 kinetics architecture:

* Uses the new :class:`Model` class with :meth:`run_single` and :meth:`run_multi` methods
* Works with the new :class:`SingleModelResult` and :class:`MultiModelResult` classes
* Integrates with the sampling framework via :class:`Sampler` interfaces
* Follows the new reaction architecture with parameter and parameter_distributions attributes

This ensures seamless integration with other v2 components while providing backward compatibility through clear API design.