===============
Getting Started
===============


Installing kinetics
-------------------

To install the latest stable version of kinetics using pip, together with all the
dependencies, run the following command:

::

    pip install kinetics


To install the latest development version of kinetics, run the following
commands.  Note that the development version may be unstable and include bugs.

::

    git clone -branch dev https://github.com/willfinnigan/kinetics.git
    cd kinetics
    python setup.py develop

Installing Prerequisite Software
--------------------------------

kinetics requires `NumPy <http://www.numpy.org/>`_, `SciPy <http://www.scipy.org/>`_,
`matplotlib <http://matplotlib.org/>`_, `tqdm <https://tqdm.github.io>`_, `pandas <http://pandas.pydata.org>`_,
`SALib <https://salib.readthedocs.io>`_, and `deap <https://deap.readthedocs.io/en/master/>`_,installed on your computer.
Using `pip <https://pip.pypa.io/en/stable/installing/>`_, these libraries can be installed with the following command:

::

    pip install numpy
    pip install scipy
    pip install matplotlib
    pip install tqdm
    pip install pandas
    pip install salib
    pip install deap

The packages are normally included with most Python bundles, such as Anaconda and Canopy.
In any case, they are installed automatically when using pip or setuptools to install
kinetics.