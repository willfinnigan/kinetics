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
    import kinetics.Uncertainty as ua
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

    model.append(step1)
    model.append(step2)
    model.species = {'e' : 1,
                     'a' : 100}

    model.setup_model()

    model.run_model()
    model.plot_substrate('a')
    model.plot_substrate('e')
    model.plot_substrate('ea')
    model.plot_substrate('p')
    plt.show()

    samples = ua.make_samples_from_distributions(model, num_samples=1000)
    outputs = ua.run_all_models(model, samples, logging=True)
    all_runs_dataframes = ua.dataframes_all_runs(model, outputs)
    ua.plot_substrate('a', all_runs_dataframes, colour='blue', alpha=0.01, linewidth=5)
    ua.plot_substrate('e', all_runs_dataframes, colour='darkorange', alpha=0.01, linewidth=5)
    ua.plot_substrate('ea', all_runs_dataframes, colour='green', alpha=0.01, linewidth=5)
    ua.plot_substrate('p', all_runs_dataframes, colour='purple', alpha=0.01, linewidth=5)
    plt.show()

2.  Make your own reaction class.
--------------------------------
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