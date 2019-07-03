==========
Why Model?
==========

A model is simply a mathematical description for what an experimenter ‘thinks’ will happen in a reaction.
For single enzymes, relatively simple mental or back-of-the-envelope models are often sufficient.
For example, it is fairly obvious in the majority of cases that more enzyme results in a faster reaction.
However, as we begin to build more complicated multi-enzyme reactions,
the increasing complexity of these systems requires a more methodical approach.
Constructing a kinetic model allows the dynamics of such a system to be investigated in silico,
and to ask questions such as, more of which enzyme gives a faster reaction?

What is a model?
----------------
A model is a mathematical description of how we think something is.

In the cases of the deterministic kinetic models this package deals with,
a model is a set of ordinary differential equations (ODEs),
which describe changes in substrate concentrations over time.

To define the rate at which substrates concentrations change over time, we define rate laws.

The most well known of these is the Michaelis-Menton equation.  A--enz-->B

Our rate law would be:

.. math::
    rate = \frac{c_{enz}\cdot k_{cat}\cdot c_{A}}{c_{A}+K_{M}^{A}}

And our Ordinary Differential Equations (ODEs) are:

.. math::
    \frac{dA}{dt} = -rate

.. math::
    \frac{dB}{dt} = +rate



