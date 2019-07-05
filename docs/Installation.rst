============
Installation
============

It is highly recommended to use a distribution such as anaconda, which already contains many of the required packages.

Installing kinetics
-------------------

To install the latest stable version of kinetics using pip, together with all the
dependencies, run the following command:

::

    pip install kinetics


Running in google colab
-----------------------
An easy way to get started quickly is to use a `google colab <https://colab.research.google.com/>`_. notebook.

In the first cell of the notebook, run  ``!pip install kinetics`` to install the kinetics package.

Try this block of code in a `google colab <https://colab.research.google.com/>`_. notebook to get started quickly..

.. code:: python

    !pip install kinetics
    import kinetics

    # Define reactions
    enzyme_1 = kinetics.Uni(kcat='enz1_kcat', kma='enz1_km', enz='enz_1', a='A',
                            substrates=['A'], products=['B'])

    enzyme_1.parameters = {'enz1_kcat' : 100,
                           'enz1_km' : 8000}

    # Set up the model
    model = kinetics.Model(logging=False)
    model.append(enzyme_1)
    model.set_time(0, 120, 1000) # 120 mins, 1000 timepoints.

    # Set starting concentrations
    model.species = {"A" : 10000,
                     "enz_1" : 4,}
    model.setup_model()

    # Run the model
    model.run_model()
    model.plot_substrate('A')
    model.plot_substrate('B', plot=True)


Prerequisite Software
---------------------
You shouldn't need to worry about this if using the **anaconda** python distribution, and ``pip install kinetics``.


kinetics requires `NumPy <http://www.numpy.org/>`_, `SciPy <http://www.scipy.org/>`_,
`matplotlib <http://matplotlib.org/>`_, `tqdm <https://tqdm.github.io>`_, `pandas <http://pandas.pydata.org>`_,
`SALib <https://salib.readthedocs.io>`_, `seaborn <http://seaborn.pydata.org>`_, and `deap <https://deap.readthedocs.io/en/master/>`_,installed on your computer.
Using `pip <https://pip.pypa.io/en/stable/installing/>`_, these libraries can be installed with the following command:

::

    pip install numpy
    pip install scipy
    pip install matplotlib
    pip install tqdm
    pip install pandas
    pip install salib
    pip install deap
    pip install seaborn

The packages are normally included with most Python bundles, such as Anaconda.
In any case, they are installed automatically when using pip or setuptools to install
kinetics.

