================
Custom Reactions
================

It's not possible to pre-define every rate equation anyone would ever need.  However many of the most common rate equations are already set up (see section on Reactions).

To make a reaction class for a custom rate equations we need to define a new class which inherits from ``kinetics.Reaction``

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