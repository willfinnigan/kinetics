from kinetics.models.reaction_class import Reaction
import numpy as np

class Uni(Reaction):
    """Irreversible unimolecular Michaelis-Menten reaction.
    
    Implements the standard Michaelis-Menten kinetics for a single substrate:
    E + S ⇌ ES → E + P
    
    Rate equation: v = (kcat * E * S) / (Km + S)
    
    Where:
    - kcat: Catalytic rate constant (turnover number)
    - Km: Michaelis constant (substrate concentration at half-maximal rate)
    - E: Enzyme concentration
    - S: Substrate concentration
    
    Examples:
        >>> # Create enzyme reaction A -> B
        >>> enzyme = kinetics.Uni(kcat='k_cat', kma='K_m', enz='enzyme', a='A',
        ...                       substrates=['A'], products=['B'])
        >>> enzyme.parameters = {'k_cat': 100, 'K_m': 1000}
        >>> model.add_reaction(enzyme)
    """

    def __init__(self,
                 kcat=None, kma=None, a=None, enz=None,
                 substrates=[], products=[]):
        """Initialize unimolecular Michaelis-Menten reaction.
        
        Args:
            kcat: Parameter name for catalytic rate constant
            kma: Parameter name for Michaelis constant for substrate A
            a: Species name for substrate A
            enz: Species name for enzyme
            substrates: List of substrate species names
            products: List of product species names
        """

        super().__init__()

        self.reaction_substrate_names = [a, enz]
        self.parameter_names = [kcat, kma]

        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):
        """Calculate reaction rate using Michaelis-Menten equation.
        
        Args:
            substrates: List containing [substrate_concentration, enzyme_concentration]
            parameters: List containing [kcat, Km]
            
        Returns:
            Reaction rate (concentration/time)
        """
        # Substrates
        a = substrates[0]
        enz = substrates[1]

        # Parameters
        kcat = parameters[0]
        kma = parameters[1]

        rate = kcat * enz * (a / (kma + a))

        return rate


class Bi(Reaction):
    
    def __init__(self,
                 kcat=None, kma=None, kmb=None,
                 a=None, b=None, enz=None,
                 substrates=[], products=[]):

        super().__init__()

        self.reaction_substrate_names = [a,b,enz]
        self.parameter_names=[kcat, kma, kmb]

        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):
        # Substrates
        a = substrates[0]
        b = substrates[1]
        enz = substrates[2]

        # Parameters
        kcat = parameters[0]
        kma = parameters[1]
        kmb = parameters[2]

        # Rate equation
        rate =  (kcat * enz) * (a / (kma + a)) * (b / (kmb + b))

        return rate

class Bi_ternary_complex(Reaction):

    def __init__(self,
                 kcat=None, kma=None, kmb=None, kia=None,
                 a=None, b=None, enz=None,
                 substrates=[], products=[]):

        super().__init__()

        self.reaction_substrate_names = [a,b,enz]
        self.parameter_names=[kcat, kma, kmb, kia]

        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):
        # Substrates
        a = substrates[0]
        b = substrates[1]
        enz = substrates[2]

        # Parameters
        kcat = parameters[0]
        kma = parameters[1]
        kmb = parameters[2]
        kia = parameters[3]

        num = kcat * enz * a * b
        den = ((kia * kmb) + (kmb * a) + (kma * b) + (a * b))

        rate = num / den

        return rate

class Bi_ping_pong(Reaction):

    def __init__(self,
                 kcat=None, kma=None, kmb=None, a=None, b=None, enz=None,
                 substrates=[], products=[]):

        super().__init__()

        self.reaction_substrate_names = [a, b, enz]
        self.parameter_names=[kcat, kma, kmb]

        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):
        # Substrates
        a = substrates[0]
        b = substrates[1]
        enz = substrates[2]

        # Parameters
        kcat = parameters[0]
        kma = parameters[1]
        kmb = parameters[2]

        rate = (kcat * enz * a * b) / ((kmb * a) + (kma * b) + (a * b))

        return rate

class Ter_seq_redam(Reaction):

    def __init__(self,
                 kcat=None, kma=None, kmb=None, kmc=None, kia=None, kib=None,
                 enz=None, a=None, b=None, c=None,
                 substrates=[], products=[]):

        super().__init__()

        self.reaction_substrate_names = [a, b, c, enz]
        self.parameter_names=[kcat, kma, kmb, kmc, kia, kib]

        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):
        # Substrates
        a = substrates[0]
        b = substrates[1]
        c = substrates[2]
        enz = substrates[3]

        # Parameters
        kcat = parameters[0]
        kma = parameters[1]
        kmb = parameters[2]
        kmc = parameters[3]
        kia = parameters[4]
        kib = parameters[5]

        numerator = kcat * enz * a * b * c

        denominator = (kia * kib * kmc) + (kib * kmc * a) + (kia * kmb * c) + (kmc * a * b) + (kmb * a * c) + (kma * b * c) + (a * b * c)

        rate = numerator / denominator

        return rate

class Ter_seq_car(Reaction):

    def __init__(self,
                 kcat=None,
                 kma=None, kmb=None, kmc=None,
                 kia=None,
                 enz=None, a=None, b=None, c=None,
                 substrates=[], products=[]):

        super().__init__()

        self.reaction_substrate_names = [a, b, c, enz]
        self.parameter_names=[kcat, kma, kmb, kmc, kia]

        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):
        # Substrates
        a = substrates[0]
        b = substrates[1]
        c = substrates[2]
        enz = substrates[3]

        # Parameters
        kcat = parameters[0]
        kma = parameters[1]
        kmb = parameters[2]
        kmc = parameters[3]
        kia = parameters[4]

        rate = (kcat * enz * a * b * c) / ((kia * c) + (kmc * a * b) + (kmb * a * c) + (kma * b * c) + (a * b * c))

        return rate

class Bi_ternary_complex_small_kma(Reaction):

    def __init__(self,
                 kcat=None, kmb=None, kia=None,
                 a=None, b=None, enz=None,
                 substrates=[], products=[]):

        super().__init__()

        self.reaction_substrate_names = [a,b,enz]
        self.parameter_names=[kcat, kmb, kia]

        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):
        # Substrates
        a = substrates[0]
        b = substrates[1]
        enz = substrates[2]

        # Parameters
        kcat = parameters[0]
        kmb = parameters[1]
        kia = parameters[2]

        rate = (kcat * enz * a * b) / ((kia * kmb) + (kmb * a) + (a * b))


        return rate