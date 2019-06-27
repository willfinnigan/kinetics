kinetics
========

kinetics is a package for modelling reactions using ordinary differential equations.
It's primarily aimed at modelling enzyme reactions, although can be used for other purposes.

kinetics uses scipy.integrate.odeint to solve ordinary differential equations,
but extends upon this to allow the use of parameter distributions rather than single parameter values.

.. image:: images/simple_example1.png
   :scale: 50
   :alt: example plot

.. toctree::
   :maxdepth: 2

   getting-started



Features
--------
Construct systems of ODEs simply by selecting suitable rate equations and naming parameters and species.
Use either simple parameter values or probability distributions.
Run sensitivity analysis using SALib
Easily plot model runs using predefined plotting functions



Support
-------

wjafinnigan@gmail.com

License
-------

The project is licensed under the MIT license.
