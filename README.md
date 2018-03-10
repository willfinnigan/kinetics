# Kinetics
The core of this module uses the scipy.integrate.odeint to solve ordinary differential equations to model enzyme reactions

The purpose of this module is to make it easier to deal with large numbers of species and reactions.

It does this by taking species and parameters as dictionaries for which the order does not matter.  The species dictionary is converted into two sets of ordered lists, one for the names and one for the intial values.

The rate of each enzyme reaction can then be determined in turn, and the change in each substrate concentration altered using the index position of the substrate name for yprime.

This also makes it much easier to carry out uncertainty and sensitivity analysis

This module uses SALib to carry out uncertainty and sensitivity analysis, with a class provided for each.

The outputs are saved as pandas dataframes.   

This module works well using a jupyter notebook or jupyter lab.  Please see 'esterase_only.ipynb' in the example folder





