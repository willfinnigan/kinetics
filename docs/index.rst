kinetics
========

kinetics is a package for modelling reactions using ordinary differential equations.
It's primarily aimed at modelling enzyme reactions, although can be used for other purposes.

kinetics uses scipy.integrate.odeint to solve ordinary differential equations,
but extends upon this to allow the use of parameter distributions rather than single parameter values.

Features
Construct systems of ODEs simply by selecting suitable rate equations and naming parameters and species.
Use either simple parameter values or probability distributions.
Run sensitivity analysis using SALib
Easily plot model runs using predefined plotting functions







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
