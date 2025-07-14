============
Installation
============

It is highly recommended to use a distribution such as anaconda, which already contains many of the required packages.

Installing kinetics
-------------------

Installation
~~~~~~~~~~~~

To install kinetics with all functionality:

::

    pip install kinetics


Development Installation
~~~~~~~~~~~~~~~~~~~~~~~~

To install for development:

::

    git clone https://github.com/willfinnigan/kinetics.git
    cd kinetics
    pip install -e .


Running in google colab
-----------------------
An easy way to get started quickly is to use a `google colab <https://colab.research.google.com/>`_. notebook.

In the first cell of the notebook, run  ``!pip install kinetics`` to install the kinetics package.

Try this block of code in a `google colab <https://colab.research.google.com/>`_. notebook to get started quickly..

.. code:: python

    # In Jupyter notebook, run: !pip install kinetics
    import kinetics

    # Define reactions
    enzyme_1 = kinetics.Uni(kcat='enz1_kcat', kma='enz1_km', enz='enz_1', a='A',
                            substrates=['A'], products=['B'])

    enzyme_1.parameters = {'enz1_kcat' : 100,
                           'enz1_km' : 8000}

    # Set up the model
    model = kinetics.Model()
    model.add_reaction(enzyme_1)
    model.set_time(0, 120, 1000) # 120 mins, 1000 timepoints.

    # Set starting concentrations
    starting_concentrations = {"A": 10000, "enz_1": 4}

    # Run the model
    result = model.run_single(starting_concentrations)
    result.plot('A')
    result.plot('B')


Dependencies
------------

All Dependencies Included
~~~~~~~~~~~~~~~~~~~~~~~~~

kinetics installs all these packages automatically:

- `NumPy <http://www.numpy.org/>`_ - numerical computing
- `SciPy <http://www.scipy.org/>`_ - scientific computing  
- `matplotlib <http://matplotlib.org/>`_ - plotting
- `pandas <http://pandas.pydata.org>`_ - data structures
- `tqdm <https://tqdm.github.io>`_ - progress bars
- `SALib <https://salib.readthedocs.io>`_ - sensitivity analysis
- `seaborn <http://seaborn.pydata.org>`_ - statistical plotting
- `pytest <https://docs.pytest.org/>`_ - testing framework
- `JAX <https://jax.readthedocs.io/>`_ - high-performance numerical computing
- `diffrax <https://docs.kidger.site/diffrax/>`_ - JAX-based differential equation solver
- `deap <https://deap.readthedocs.io/en/master/>`_ - genetic algorithms and optimization

Modern Package Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

kinetics now uses modern Python packaging with ``pyproject.toml`` configuration. 
The project follows PEP 517/518 standards for build systems and dependency management.

