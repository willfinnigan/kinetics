from kinetics.reaction_classes.reaction_base_class import Reaction


class Uni_mass_action_eq(kinetics.Reaction):

    def __init__(self,
                 keq=None, kf=None,
                 a='', p='',
                 substrates=[], products=[]):

        super().__init__()

        self.parameter_names=[keq, kf]
        self.reaction_substrate_names = [a, p]
        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):

        # Substrates
        a = substrates[0]
        p = substrates[1]

        # Parameters
        keq = parameters[0]
        kf = parameters[1]

        fwd = kf*a
        rev = (kf/keq)*p

        rate = fwd-rev

        return rate