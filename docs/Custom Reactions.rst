================
Custom Reactions
================

It's not possible to pre-define every rate equation anyone would ever need.
However many of the most common rate equations are already set up (see section on Reactions).

There are two options for custom rate equations.

1.  Use the Generic Reaction Class
----------------------------------
This is a reaction class which lets you specify your own rate equation.

.. autoclass:: kinetics.Generic
    :noindex:

An example which models a UniUni enzyme

.. code:: python

    import kinetics
    import matplotlib.pyplot as plt
    from scipy.stats import norm

    step1 = kinetics.Generic(params=['k1', 'k_1'], species=['a','e','ea'],
                             rate_equation='(k1*a*e)-(k_1*ea)',
                             substrates=['a', 'e'], products=['ea'])

    step2 = kinetics.Generic(params=['k2','k_2'], species=['ea','e','p'],
                             rate_equation='(k2*ea)-(k_2*ea*p)',
                             substrates=['ea'], products=['e', 'p'])

    step1.parameters = {'k1' : 0.1,
                        'k_1' : 0.001}

    step2.parameters = {'k2' : 100,
                        'k_2': 0.1}

    step1.parameter_distributions = {'k1' : norm(0.1, 0.01),
                                     'k_1' : norm(0.001, 0.0001)}

    step2.parameter_distributions = {'k2' : norm(100, 10),
                                     'k_2' : norm(0.1, 0.01)}

    model = kinetics.Model()
    model.add_reaction(step1)
    model.add_reaction(step2)
    model.set_time(0, 100, 1000)
    
    starting_concentrations = {'e': 1, 'a': 100}

    # Run single model
    result = model.run_single(starting_concentrations)
    result.plot('a')
    result.plot('e')
    result.plot('ea')
    result.plot('p')
    plt.show()

    # For uncertainty analysis, use the sampling functionality
    sampler = kinetics.ScipyDist_Sampler(num_samples=1000)
    multi_result = model.run_multi(starting_concentrations, sampler)
    
    # Plot results with uncertainty bands
    multi_result.plot('a', quartile=95)
    multi_result.plot('e', quartile=95)
    multi_result.plot('ea', quartile=95)
    multi_result.plot('p', quartile=95)
    plt.show()

2.  Make your own reaction class.
---------------------------------
**This might be useful if its going to be re-used alot**

To make a reaction class for a custom rate equation we need to define a new class which inherits from ``kinetics.Reaction``

The new class needs two funcions.  an __init__() function and a calculate_rate() function.

The following code example provides an example for doing this:

.. code:: python

    class My_New_Reaction(kinetics.Reaction):

        def __init__(self,
                     param1='', param2='', species1='', species2='',
                     substrates=[], products=[]):

        # This is required to inherit from kinetics.Reaction
        super().__init__()

        # Set parameter and substrates names from the arguments passed in.  The order is important here.
        self.parameter_names=[param1, param2]
        self.reaction_substrate_names = [species1, species2]

        # Set the substrates and products from the arguments passed in.
        # Substrates are used up in the reaction, while produces are generated.
        self.substrates = substrates
        self.products = products

        def calculate_rate(self, substrates, parameters):

            # This function is used to calculate the rate at each time step in the model
            # It takes substrates and parameters as arguments, which are lists with the same order as we defined in __init__.

            # Substrates
            species1 = substrates[0]
            species2 = substrates[1]

            # Parameters
            param1 = parameters[0]
            param2 = parameters[1]

            # This is where the rate equation goes.  An example is shown.
            rate = param1*species1 + param2*species2

            return rate