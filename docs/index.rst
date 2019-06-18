kinetics
========

kinetics is a module for modelling reactions using ordinary differential equations.
It's primarily aimed at modelling enzyme reactions, although can be used for other purposes.

kinetics uses scipy.integrate.odeint to solve ordinary differential equations,
but extends upon this to allow the use of parameter distributions rather than single parameter values.
Currently only uniform distributions are used, which will hopefully be extended in the near future.

.. toctree::
   :maxdepth: 2

   getting-started


Features
--------



Installation
------------

Install kinetics by running:

    pip install kinetics

Contribute
----------

- Source Code: http://github.com/willfinnigan/kinetics

Support
-------

wjafinnigan@gmail.com

License
-------

The project is licensed under the MIT license.

.. automodule:: io
   :kinetics.Model:
