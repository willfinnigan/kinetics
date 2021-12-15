from kinetics.reaction_classes.reaction_base_class import Reaction

class BiBi_Ordered_rev_eq(kinetics.Reaction):

    def __init__(self,
                 keq=None,
                 kcatf=None, kcatr=None,
                 kma=None, kmb=None, kmp=None, kmq=None,
                 kib=None, kip=None, kia=None,
                 a='', b='', p='', q='', enz='',
                 substrates=[], products=[]):

        super().__init__()

        self.parameter_names = [keq, kcatf, kcatr, kma, kmb, kmp, kmq, kia, kib, kip]
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
        keq = parameters[0]
        kcatf = parameters[1]
        kcatr = parameters[2]
        kma = parameters[3]
        kmb = parameters[4]
        kmp = parameters[5]
        kmq = parameters[6]
        kia = parameters[7]
        kib = parameters[8]
        kip = parameters[9]


        Vf = kcatf*enz
        Vr = kcatr*enz

        num = Vf*(a*b-((p*q)/keq))

        den1 = a*b*(1+p/kip)
        den2 = (kma*b) + (kmb*(a*kia))
        den3 = (Vf/(Vr*keq)) * ( (kmq*p*(1+(a/kia))) + (q* (kmp*( 1+((kma*b)/(kia*kmb)) + p*(1+(b/kib)))) ))

        rate = num / (den1 + den2 + den3)
        return rate

class UniUni_rev_eq(kinetics.Reaction):

    def __init__(self,
                 keq=None, kcatf=None, kma=None, kmp=None,
                 a='', p='', enz='',
                 substrates=[], products=[]):

        super().__init__()

        self.parameter_names=[keq, kcatf, kma, kmp]
        self.reaction_substrate_names = [a, p, enz]
        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):

        # Substrates
        a = substrates[0]
        p = substrates[1]
        enz = substrates[2]

        # Parameters
        keq = parameters[0]
        kcatf = parameters[1]
        kma = parameters[2]
        kmp = parameters[3]

        numerator = enz*kcatf*(a-(p/keq))
        denominator = a+(kma*(1+(p/kmp)))
        rate = numerator / denominator

        return rate