kinetics
--------

kinetics is a package for modelling reactions using ordinary differential equations.
It's primarily aimed at modelling enzyme reactions, although can be used for other purposes.

kinetics uses scipy.integrate.odeint to solve ordinary differential equations,
but extends upon this to allow the use of parameter distributions rather than single parameter values.
This allows error to be incorporated into the modelling.

kinetics uses `scipy's probability distributions <https://docs.scipy.org/doc/scipy/reference/stats.html/>`_, with a large number of distributions to choose from.
Typically uniform , normal, log-uniform or log-normal distributions are used.

Python implementations of commonly used sensitivity analysis methods.
Useful in systems modeling to calculate the effects of model inputs or
exogenous factors on outputs of interest.


**Documentation:** `ReadTheDocs <http://salib.readthedocs.org>`__

**Requirements:** `NumPy <http://www.numpy.org/>`__,
`SciPy <http://www.scipy.org/>`__,
`matplotlib <http://matplotlib.org/>`__,
`pandas <http://https://pandas.pydata.org/>`__,
Python 3 (from SALib v1.2 onwards SALib does not officially support Python 2)

**Installation:** ``pip install SALib`` or ``python setup.py install`` or ``conda install SALib``