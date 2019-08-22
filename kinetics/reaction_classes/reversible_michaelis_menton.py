from kinetics.reaction_classes.reaction_base_class import Reaction

""" Michaelis-Menten reversible equations """
class UniUni_rev(Reaction):
    r"""
    """

    def __init__(self,
                 kcatf=None, kcatr=None, kma=None, kmp=None, a=None, p=None, enz=None,
                 substrates=[], products=[]):

        super().__init__()

        self.reaction_substrate_names = [a, p, enz]
        self.parameter_names=[kcatf, kcatr, kma, kmp]

        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):
        # Substrates
        a = substrates[0]
        p = substrates[1]
        enz = substrates[2]

        # Parameters
        kcatf = parameters[0]
        kcatr = parameters[1]
        kma = parameters[2]
        kmp = parameters[3]

        rate = ((kcatf * enz * a) - (kcatr * enz * p)) / (1 + (a / kma) + (p / kmp))

        return rate

class BiBi_Ordered_rev(Reaction):

    def __init__(self,
                 kcatf=None, kcatr=None,
                 kmb=None, kia=None, kib=None, kmp=None, kip=None, kiq=None,
                 enz=None, a=None, b=None, p=None, q=None,
                 substrates=[], products=[]):

        super().__init__()

        self.reaction_substrate_names = [a,b,p,q,enz]
        self.parameter_names=[kcatf, kcatr, kmb, kia, kib, kmp, kip, kiq]

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
        kcatf = parameters[0]
        kcatr = parameters[1]
        kmb = parameters[2]
        kia = parameters[3]
        kib = parameters[4]
        kmp = parameters[5]
        kip = parameters[6]
        kiq = parameters[7]

        # Rate equation
        numerator = ((enz * kcatf * a * b) / (kia * kmb)) - ((enz * kcatr * p * q) / (kmp * kiq))
        denominator = 1 + (a / kia) + (b / kib) + (q / kiq) + (p / kip) + ((a * b) / (kia * kmb)) + ((p * q) / (kmp * kiq))

        return (numerator / denominator)

class BiBi_Random_rev(Reaction):

    def __init__(self,
                 kcatf=None, kcatr=None, kmb=None, kia=None, kib=None, kmp=None, kip=None, kiq=None,
                 a=None, b=None, p=None, q=None, enz=None,
                 substrates=[], products=[]):

        super().__init__()

        self.reaction_substrate_names = [a, b, p, q, enz]
        self.parameter_names=[kcatf, kcatr, kmb, kia, kib, kmp, kip, kiq]

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
        kcatf = parameters[0]
        kcatr = parameters[1]
        kmb = parameters[2]
        kia = parameters[3]
        kib = parameters[4]
        kmp = parameters[5]
        kip = parameters[6]
        kiq = parameters[7]

        num = ((kcatf * enz * a * b) / (kia * kmb)) - ((kcatr * enz * p * q) / (kmp * kiq))

        dom = 1 + (a / kia) + (b / kib) + (p / kip) + (q / kiq) + ((a * b) / (kia * kmb)) + ((p * q) / (kmp * kiq))

        rate = num / dom

        return rate

class BiBi_Pingpong_rev(Reaction):

    def __init__(self,
                 kcatf=None, kma=None, kmb=None, kia=None,
                 kcatr=None, kmp=None, kmq=None, kip=None, kiq=None,
                 enz=None, a=None, b=None, p=None, q=None,
                 substrates=[], products=[]):

        super().__init__()

        self.reaction_substrate_names = [a, b, p, q, enz]
        self.parameter_names=[kcatf, kcatr, kma, kmb, kia, kmp, kmq, kip, kiq]

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
        kcatf = parameters[0]
        kcatr = parameters[1]
        kma = parameters[2]
        kmb = parameters[3]
        kia = parameters[4]
        kmp = parameters[5]
        kmq = parameters[6]
        kip = parameters[7]
        kiq = parameters[8]

        num = ((kcatf * enz * a * b) / (kia * kmb)) - ((kcatr * enz * p * q) / kip * kmq)

        den = (a / kia) + ((kma * b) / (kia * kmb)) + (p / kip) + ((kmp * q) / (kip * kmq)) + (
        (a * b) / (kia * kmb)) + ((a * p) / (kia * kip)) + ((kma * b * q) / (kia * kmb * kiq)) + ((p * q) / (kip * kmq))

        return num / den
