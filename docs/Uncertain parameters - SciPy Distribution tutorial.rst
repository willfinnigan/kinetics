==================================================
Uncertain parameters - SciPy Distribution tutorial
==================================================

This tutorial uses the same example as the simple tutorial, but demonstrates the use of probability distributions rather than single parameter values.

Again I recommend running this using a jupyter notebook or google colab.

Here is a direct link to this example in google colab.
Or here is a direct link to this workbook in google colab.
https://colab.research.google.com/github/willfinnigan/kinetics/blob/master/examples/Advanced_example.ipynb



.. image:: images/example_system.png
   :scale: 20
   :alt: graphical abstract

Define reactions with uncertain parameters
------------------------------------------
We can describe the uncertainty we have for a parameter using a probability distribution.

Where parameters have been characterised, resulting in a standard error, we can use this to describe a normal distribution.

Where we roughly know where a parameter is, but don't want to make any suggestion as to more or less likley values, we can use a uniform distribution.

Where we have absolutely no idea what value a parameter takes, we can use a log-uniform distribution (reciprocal) to describe it was being equally as likely to be within a number of orders of magnitude.

Log-normal distributions can also be used.  A recent paper describes how log-normal distributions can be created taking into account numerous literature values.

To use these probability distributions in our modelling, we can employ the probability distributions available through `scipy
<https://docs.scipy.org/doc/scipy/reference/stats.html>`_.

Define our reactions as before, but this time we specify reaction.parameter_distributions.  Make sure to import these from scipy.

.. code:: python

    import kinetics
    from scipy.stats import reciprocal, uniform, norm

    # Define reactions
    enzyme_1 = kinetics.Uni(kcat='enz1_kcat', kma='enz1_km', enz='enz_1', a='A',
                            substrates=['A'], products=['B'])

    enzyme_1.parameter_distributions = {'enz1_kcat' : norm(100,12),
                                        'enz1_km' : uniform(2000, 6000)}

    enzyme_2 = kinetics.Uni(kcat='enz2_kcat', kma='enz2_km', enz='enz_2', a='B',
                            substrates=['B'], products=['C'])

    enzyme_2.parameter_distributions = {'enz2_kcat' : norm(30, 5),
                                        'enz2_km' : reciprocal(1,10000)}

Next we define our model as before.

.. code:: python

    # Set up the model
    model = kinetics.Model()
    model.add_reaction(enzyme_1)
    model.add_reaction(enzyme_2)
    model.set_time(0, 120, 1000)

We can include uncertainty in some or all of the starting species concentrations.
Here we have specified uncertainty in the enzyme concentraions using a normal distribution with a standard deviation of 5% of the starting value.
We have specified no uncertainty in the starting concentration of A.

.. code:: python

    # Set starting concentrations (mix of fixed values and distributions)
    species_dict = {"A": 10000,
                    "enz_1": norm(4, 4*0.05),
                    "enz_2": norm(10, 10*0.05)}

Running the model with a single set of parameter values
-------------------------------------------------------
We can run the model exactly as in the simple example, and we will get a single prediction for each substrate.
Running the model this way will use the mean of each probability distribution specified, unless a different value is specified.

.. code:: python

    result = model.run_single(species_dict)
    result.plot('A')
    result.plot('B')
    result.plot('C')
    plt.show()

.. image:: images/simple_example1.png
   :scale: 25
   :alt: example plot

Running the model by sampling within the probability distributions
------------------------------------------------------------------
However we would like to run lots of models, sampling within our probability distributions.

In the new API, we use a sampler class to generate samples from distributions and then run multiple models using ``model.run_multi()``.

.. code:: python

    # Create a sampler and run the model 1000 times, sampling from distributions
    sampler = kinetics.ScipyDist_Sampler(num_samples=1000)
    result = model.run_multi({"A": 10000, "enz_1": 4, "enz_2": 10}, sampler)

Plotting the data
-----------------
The result object from ``model.run_multi()`` provides methods to access and plot the data.

``result.dataframe()`` returns a dictionary containing dataframes for each species showing all runs.

``result.dataframe_quartiles(quartile=95)`` returns dataframes with confidence intervals (High, Low, Mean).

These dataframes can be exported for further use, or the result object provides built-in plotting methods.

Plotting graphs with confidence intervals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Plotting the 95% confidence intervals can look neater, but we lose some information on the outliers by doing this.

.. code:: python

    # Plot model runs with 95% confidence intervals
    result.plot('A', quartile=95)
    result.plot('B', quartile=95)  
    result.plot('C', quartile=95)
    plt.show()

.. image:: images/advanced_example1.png
   :scale: 25
   :alt: example plot

Plotting graphs showing all runs (spagetti plots)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Alternatively we can plot every single run.  With 1000 runs this can look a bit chaotic, and it may be clearer to plot each substrate on its own graph.

.. code:: python

    # Plot all individual model runs
    result.plot_all('A')
    result.plot_all('B')
    result.plot_all('C')
    plt.show()

.. image:: images/advanced_example2.png
   :scale: 25
   :alt: example plot

Of course the dataframes are also available to be used as the output, possibly to create your own graphs or for other analysis.

Complete code
----------------------------------------

.. code:: python

    import kinetics
    import matplotlib.pyplot as plt
    from scipy.stats import reciprocal, uniform, norm
    %config InlineBackend.figure_format ='retina'

    # Define reactions
    enzyme_1 = kinetics.Uni(kcat='enz1_kcat', kma='enz1_km', enz='enz_1', a='A',
                            substrates=['A'], products=['B'])

    enzyme_1.parameter_distributions = {'enz1_kcat' : norm(100,12),
                                        'enz1_km' : uniform(2000, 6000)}

    enzyme_2 = kinetics.Uni(kcat='enz2_kcat', kma='enz2_km', enz='enz_2', a='B',
                            substrates=['B'], products=['C'])

    enzyme_2.parameter_distributions = {'enz2_kcat' : norm(30, 5),
                                        'enz2_km' : reciprocal(1,10000)}

    # Set up the model
    model = kinetics.Model()
    model.add_reaction(enzyme_1)
    model.add_reaction(enzyme_2)
    model.set_time(0, 120, 1000)

    # Set starting concentrations (fixed values and distributions)
    species_dict = {"A": 10000,
                    "enz_1": norm(4, 4*0.05),
                    "enz_2": norm(10, 10*0.05)}

    # Run a single model first with mean values
    single_result = model.run_single(species_dict)
    single_result.plot('A')
    single_result.plot('B')
    single_result.plot('C')

    # Run the model 1000 times, sampling from distributions
    sampler = kinetics.ScipyDist_Sampler(num_samples=1000)
    multi_result = model.run_multi(species_dict, sampler)

    # Plot model runs with 95% confidence intervals
    multi_result.plot('A', quartile=95)
    multi_result.plot('B', quartile=95)
    multi_result.plot('C', quartile=95)
    plt.show()

    # Plot all individual model runs
    multi_result.plot_all('A')
    multi_result.plot_all('B')
    multi_result.plot_all('C')
    plt.show()


