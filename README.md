# Kinetics
The core of this module uses the scipy.integrate.odeint to solve ordinary differential equations to model enzyme reactions

The purpose of this module is to make it easier to deal with large numbers of species and reactions.

It does this by taking species and parameters as dictionaries for which the order does not matter.  The species dictionary is converted into two sets of ordered lists, one for the names and one for the intial values.

The rate of each enzyme reaction can then be determined in turn, and the change in each substrate concentration altered using the index position of the substrate name for yprime.

This module also uses SALib to carry out uncertainty and sensitivity analysis.

The outputs from the modelling runs are saved in txt documents which can easily be copied into software such as excel or graphpad.



