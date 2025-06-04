from kinetics.reaction_classes.reaction_base_class import Reaction

class Bi_Uni_sep_eq(Reaction):

    def __init__(self,
                 kcat=None,
                 kma=None, kmb=None,
                 kmp=None,
                 keq=None,
                 a='', b='', p='', enz='',
                 substrates=[], products=[]):

        super().__init__()

        self.parameter_names = [kcat, kma, kmb, kmp, keq]
        self.reaction_substrate_names = [a, b, p, enz]
        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):

        # Substrates
        a = substrates[0]
        b = substrates[1]
        p = substrates[2]
        enz = substrates[3]

        # Parameters
        kcat = parameters[0]
        kma = parameters[1]
        kmb = parameters[2]
        kmp = parameters[3]
        keq = parameters[4]

        catalytic_capacity = enz * kcat

        thermodynamic_driving_force = 1 - ( p/(a*b) / keq )

        subs = (a/kma)*(b/kmb)
        prods = (p/kmp)
        substrate_saturation = subs / (subs + prods)

        return catalytic_capacity * substrate_saturation * thermodynamic_driving_force


class Bi_Bi_sep_eq(Reaction):

    def __init__(self,
                 kcat=None,
                 kma=None, kmb=None,
                 kmp=None, kmq=None,
                 keq=None,
                 a='', b='', p='', q='', enz='',
                 substrates=[], products=[]):

        super().__init__()

        self.parameter_names = [kcat, kma, kmb, kmp, kmq, keq]
        self.reaction_substrate_names = [a, b, p, q, enz]
        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):

        # Substrates
        a = substrates[0]
        b = substrates[1]
        p = substrates[2]
        q = substrates[3]
        enz = substrates[4]

        # Parameters
        kcat = parameters[0]
        kma = parameters[1]
        kmb = parameters[2]
        kmp = parameters[3]
        kmq = parameters[4]
        keq = parameters[5]

        catalytic_capacity = enz * kcat

        thermodynamic_driving_force = 1 - ( (p*q)/(a*b) / keq )

        subs = (a/kma)*(b/kmb)
        prods = (p/kmp)*(q/kmq)
        substrate_saturation = subs / (subs + prods)

        return catalytic_capacity * substrate_saturation * thermodynamic_driving_force

class Tri_Tri_seq_eq(Reaction):

    def __init__(self,
                 kcat=None,
                 kma=None, kmb=None, kmc=None,
                 kmp=None, kmq=None, kmr=None,
                 keq=None,
                 a='', b='', c='', p='', q='', r='', enz='',
                 substrates=[], products=[]):

        super().__init__()

        self.parameter_names = [kcat, kma, kmb, kmc, kmp, kmq, kmr, keq]
        self.reaction_substrate_names = [a, b, c, p, q, r, enz]
        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):

        # Substrates
        a = substrates[0]
        b = substrates[1]
        c = substrates[2]
        p = substrates[3]
        q = substrates[4]
        r = substrates[5]
        enz = substrates[6]

        # Parameters
        kcat = parameters[0]
        kma = parameters[1]
        kmb = parameters[2]
        kmc = parameters[3]
        kmp = parameters[4]
        kmq = parameters[5]
        kmr = parameters[6]
        keq = parameters[7]

        catalytic_capacity = enz * kcat

        #thermodynamic_driving_force = 1 - np.exp(kjmol / (R * T))
        thermodynamic_driving_force = 1 - ( (p*q*r)/(a*b*c) / keq )     # 1-((p/s)/keq)

        subs = (a/kma)*(b/kmb)*(c/kmc)
        prods = (p/kmp)*(q/kmq)*(r/kmr)
        substrate_saturation = subs / (subs + prods)

        return catalytic_capacity * substrate_saturation * thermodynamic_driving_force