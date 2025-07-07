===================
Sensitivity Analysis
===================

This tutorial demonstrates how to perform sensitivity analysis on enzyme kinetic models using the kinetics package. Sensitivity analysis helps identify which parameters have the most influence on your model outputs.

Overview
--------

The sensitivity analysis uses **Sobol indices** to quantify parameter importance:

* **S1**: First-order sensitivity (direct effect of each parameter)
* **ST**: Total sensitivity (includes interaction effects with other parameters)

Basic Two-Enzyme Pathway Example
---------------------------------

Let's analyze a simple two-enzyme pathway: A → B → C

.. code-block:: python

    import kinetics
    import numpy as np
    from kinetics import (SalibSaltelliSampler, 
                          analyze_sensitivity_at_timepoint, 
                          analyze_sensitivity_time_to_concentration)

    # Step 1: Define the enzymatic reactions
    enzyme_1 = kinetics.Uni(
        kcat='enz1_kcat', 
        kma='enz1_km', 
        enz='enz_1', 
        a='A',
        substrates=['A'], 
        products=['B']
    )

    enzyme_2 = kinetics.Uni(
        kcat='enz2_kcat', 
        kma='enz2_km', 
        enz='enz_2', 
        a='B',
        substrates=['B'], 
        products=['C']
    )

    # Step 2: Define parameter distributions for uncertainty analysis
    # Format: (min_value, max_value)
    enzyme_1.parameter_distributions = {
        'enz1_kcat': (90, 110),      # kcat varies from 90-110 s⁻¹
        'enz1_km': (2000, 6000)      # Km varies from 2000-6000 µM
    }

    enzyme_2.parameter_distributions = {
        'enz2_kcat': (25, 35),       # kcat varies from 25-35 s⁻¹
        'enz2_km': (1, 10000)        # Km varies from 1-10000 µM (log scale)
    }

    # Step 3: Set up the model
    model = kinetics.Model()
    model.set_time(0, 1000, 100)  # 0-1000 seconds, 100 time points
    model.add_reaction(enzyme_1)
    model.add_reaction(enzyme_2)

    # Step 4: Configure the sampler
    # SalibSaltelliSampler uses Saltelli sampling for Sobol analysis
    sampler = SalibSaltelliSampler(
        num_samples=1000,              # Number of samples (more = better statistics)
        log_parameters=['enz2_km']     # Parameters to sample on log scale
    )

    # Step 5: Run the multi-parameter simulation
    solver = kinetics.SciPySolver()
    initial_conditions = {
        "A": 10000,    # 10 mM substrate A
        "enz_1": 5,    # 5 µM enzyme 1
        "enz_2": 2     # 2 µM enzyme 2
    }

    result = model.run_multi(initial_conditions, sampler, solver)

    print(f"Completed {len(result.results)} simulations")

Analysis 1: Sensitivity at Specific Timepoint
----------------------------------------------

Analyze which parameters most influence the concentration of species B at t=500 seconds:

.. code-block:: python

    # Analyze sensitivity for [B] at t=500 seconds
    sa_result = analyze_sensitivity_at_timepoint(
        result=result,
        sampler=sampler,
        species='B',
        timepoint=500,
        threshold=0.01  # Only show parameters with >1% sensitivity
    )

    # Display results
    print("Sensitivity Analysis for [B] at t=500s:")
    print(sa_result.data)

    # Visualize results
    sa_result.plot()

Expected output::

    Sensitivity Analysis for [B] at t=500s:
              S1    ST    S1_conf    ST_conf
    enz1_kcat  0.45  0.48      0.03      0.04
    enz1_km    0.32  0.35      0.02      0.03
    enz2_kcat  0.15  0.18      0.02      0.02

Analysis 2: Time to Reach Target Concentration
-----------------------------------------------

Analyze which parameters most influence the time to reach 5000 µM of species B:

.. code-block:: python

    # Analyze sensitivity for time to reach [B] = 5000 µM
    sa_result = analyze_sensitivity_time_to_concentration(
        result=result,
        sampler=sampler,
        species='B',
        concentration=5000,
        mode='>=',       # Time to reach or exceed concentration
        threshold=0.01
    )

    print("Sensitivity Analysis for time to reach [B] = 5000 µM:")
    print(sa_result.data)

    # Visualize results
    sa_result.plot()

Interpreting Results
--------------------

Sobol Indices
~~~~~~~~~~~~~

* **S1 (First-order)**: Direct effect of the parameter
* **ST (Total)**: Total effect including interactions
* **Difference (ST - S1)**: Interaction effects with other parameters

Key Insights
~~~~~~~~~~~~

* Higher values = greater parameter importance
* S1 ≈ ST suggests minimal parameter interactions
* ST >> S1 indicates strong parameter interactions
* Parameters below threshold are filtered out

Confidence Intervals
~~~~~~~~~~~~~~~~~~~~

* ``S1_conf`` and ``ST_conf`` show uncertainty in sensitivity estimates
* Increase ``num_samples`` for tighter confidence intervals

Advanced Usage
--------------

Using SALib for Sensitivity Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The sensitivity analysis functionality is powered by `SALib (Sensitivity Analysis Library) <https://salib.readthedocs.io/>`_, 
a Python library for performing global sensitivity analysis. SALib provides robust methods for:

* **Sobol Analysis**: Variance-based sensitivity analysis
* **Saltelli Sampling**: Efficient sampling scheme for Sobol indices
* **Parameter Space Exploration**: Comprehensive parameter uncertainty analysis

For more information on sensitivity analysis methods and theory, visit the `SALib documentation <https://salib.readthedocs.io/>`_.


Performance Tips
----------------

1. **Start with fewer samples** (100-500) for initial exploration
2. **Use log-scale sampling** for parameters spanning multiple orders of magnitude
3. **Increase samples** (1000-5000) for publication-quality results
4. **Use thresholds** to focus on important parameters
5. **Consider computational cost** - sensitivity analysis requires many model runs

Troubleshooting
---------------

**Error: "Not enough samples"**
    * Increase ``num_samples`` in the sampler
    * SALib requires sufficient samples for robust Sobol analysis

**Error: "Parameter not found"**
    * Check parameter names match those in reaction definitions
    * Verify parameter_distributions keys are correct

**Warning: "High confidence intervals"**
    * Increase ``num_samples`` for better statistical precision
    * Check if parameter ranges are appropriate