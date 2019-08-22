from kinetics.reaction_classes.reaction_base_class import Reaction

class Generic(Reaction):
    """
    This Reaction class allows you to specify your own rate equation.
    Enter the parameter names in params, and the substrate names used in the reaction in species.
    Type the rate equation as a string in rate_equation, using these same names.
    Enter the substrates used up, and the products made in the reaction as normal.
    """

    def __init__(self,
                 params=[], species=[],
                 rate_equation='',
                 substrates=[], products=[]):

        super().__init__()

        self.reaction_substrate_names = species
        self.parameter_names=params
        self.rate_equation = rate_equation
        self.substrates = substrates
        self.products = products


    def calculate_rate(self, substrates, parameters):

        for i, name in enumerate(self.reaction_substrate_names):
            locals().update({name: substrates[i]})

        for i, name in enumerate(self.parameter_names):
            locals().update({name: parameters[i]})

        rate = eval(self.rate_equation, locals(), globals())

        return rate