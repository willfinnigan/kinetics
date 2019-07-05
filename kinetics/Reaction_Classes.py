import numpy as np
import copy

def calculate_yprime(y, rate, substrates, products, substrate_names):
    """
    This function is used by the rate classes the user creates.

    It takes the numpy array for y_prime,
    and adds or subtracts the amount in rate to all the substrates or products listed
    Returns the new y_prime

    Args:
        y: a numpy array for the substrate values, the same order as y
        rate: the rate calculated by the user made rate equation
        substrates: list of substrates for which rate should be subtracted
        products: list of products for which rate should be added
        substrate_names: the ordered list of substrate names in the model.  Used to get the position of each substrate or product in y_prime

    Returns:
        y_prime: following the addition or subtraction of rate to the specificed substrates
    """

    y_prime = np.zeros(len(y))

    for name in substrates:
        y_prime[substrate_names.index(name)] -= rate

    for name in products:
        y_prime[substrate_names.index(name)] += rate

    return y_prime

def check_positive(y_prime):
    """
    Chack that substrate values are not negative when they shouldnt be
    """

    for i in range(len(y_prime)):
        if y_prime[i] < 0:
            y_prime[i] = 0

    return y_prime

class Reaction():

    def __init__(self, substrates=[], products=[]):

        self.reaction_substrate_names = []
        self.substrate_indexes = []
        self.substrates = substrates
        self.products = products

        self.parameters = {}
        self.parameter_distributions = {}

        self.parameter_names = []
        self.run_model_parameters = []

        self.modifiers = []

        self.check_positive = False

        self.check_limits_functions = []

    def set_parameter_defaults_to_median(self):
        for name in self.parameter_distributions:
            if name not in self.parameters:
                self.parameters[name] = self.parameter_distributions[name].median()

    def get_indexes(self, substrate_names):
        self.substrate_indexes = []
        for name in self.reaction_substrate_names:
            self.substrate_indexes.append(substrate_names.index(name))

    def get_substrates(self, y):
        substrates = []
        for index in self.substrate_indexes:
            substrates.append(y[index])

        return substrates

    def get_parameters(self, parameter_dict):
        parameters = []
        for name in self.parameter_names:
            parameters.append(parameter_dict[name])

        return parameters

    def reset_reaction(self):
        self.substrate_indexes = []
        self.run_model_parameters = []

    def add_modifier(self, modifier):
        for name in modifier.parameter_names:
            if name not in self.parameter_names:
                self.parameter_names.append(name)

        for name in modifier.substrate_names:
            if name not in self.reaction_substrate_names:
                self.reaction_substrate_names.append(name)

        self.modifiers.append(modifier)

    def calculate_modifiers(self, substrates, parameters):
        for modifier in self.modifiers:
            substrates, parameters = modifier.calc_modifier(substrates, parameters)

        return substrates, parameters

    def calculate_rate(self, substrates, parameters):
        return 0

    def reaction(self, y, substrate_names, parameter_dict):
        if self.substrate_indexes == []:
            self.get_indexes(substrate_names) # need to move this to the model

        if self.run_model_parameters == []:
            self.run_model_parameters = self.get_parameters(parameter_dict)

        for modifier in self.modifiers:
            if modifier.substrate_indexes == []:
                modifier.get_substrate_indexes(self.reaction_substrate_names)
            if modifier.parameter_indexes == []:
                modifier.get_parameter_indexes(self.parameter_names)

        substrates = self.get_substrates(y)

        substrates, parameters = self.calculate_modifiers(substrates, copy.copy(self.run_model_parameters))

        rate = self.calculate_rate(substrates, parameters)

        y_prime = calculate_yprime(y, rate, self.substrates, self.products, substrate_names)
        y_prime = self.modify_product(y_prime, substrate_names)

        if self.check_positive == True:
            y_prime = check_positive(y_prime)

        return y_prime

    def modify_product(self, y_prime, substrate_names):
        return y_prime

    def sampling_limits(self, parameter_dict):
        # Return true if parameters within limits, false if not
        for func in self.check_limits_functions:
            if func(parameter_dict) == False:
                return False

        return True

""" Michaelis-Menten irreversible equations """
class Uni(Reaction):
    r"""
    Uni
    """

    def __init__(self,
                 kcat=None, kma=None, a=None, enz=None,
                 substrates=[], products=[]):

        super().__init__()

        self.reaction_substrate_names = [a, enz]
        self.parameter_names=[kcat, kma]

        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):
        # Substrates
        a = substrates[0]
        enz = substrates[1]

        # Parameters
        kcat = parameters[0]
        kma = parameters[1]

        rate = kcat * enz * (a / (kma + a))

        return rate

class Bi(Reaction):
    r"""
    """

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
    r"""
    """

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
    r"""

    """

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
    r"""
    """

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


""" Other rate equations """
class FirstOrderRate(Reaction):

    def __init__(self,
                 k=None, a=None,
                 substrates=[], products=[]):

        super().__init__()

        self.reaction_substrate_names = [a]
        self.parameter_names=[k]

        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):
        # Substrates
        a = substrates[0]

        # Parameters
        k = parameters[0]

        return k*a

class Binding(Reaction):

    def __init__(self, k1=None, kminus1=None,
                 a=None, b=None, c=None,
                 substrates=[], products=[]):

        super().__init__()

        self.reaction_substrate_names = [a, b, c]
        self.parameter_names=[k1, kminus1]

        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):
        # Substrates
        a = substrates[0]
        b = substrates[1]
        c = substrates[2]

        # Parameters
        k1 = parameters[0]
        kminus1 = parameters[1]

        rate = (k1*a*b) - (kminus1*c)

        return rate

class BiSecondOrderRate(Reaction):

    def __init__(self,
                 k=None, a=None, b=None,
                 substrates=[], products=[]):

        super().__init__()

        self.reaction_substrate_names = [a, b]
        self.parameter_names=[k]

        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):
        # Substrates
        a = substrates[0]
        b = substrates[1]

        # Parameters
        k = parameters[0]

        return k*a*b

class Binding_kd(Reaction):

    def __init__(self, kd=None, k1=None,
                 a=None, b=None, c=None,
                 substrates=[], products=[]):

        super().__init__()

        self.reaction_substrate_names = [a, b, c]
        self.parameter_names=[kd, k1]

        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):
        # Substrates
        a = substrates[0]
        b = substrates[1]
        c = substrates[2]

        # Parameters
        kd = parameters[0]
        k1 = parameters[1]

        kminus1 = kd*k1

        rate = (k1*a*b) - (kminus1*c)

        return rate

class DiffusionEquilibrium(Reaction):
    def __init__(self, kd=None, k1=None, org_c=None, aq_c=None):
        super().__init__()

        self.reaction_substrate_names = [org_c, aq_c]
        self.parameter_names = [kd, k1]

        self.substrates = [org_c]
        self.products = [aq_c]

    def calculate_rate(self, substrates, parameters):
        # Substrates
        org_c = substrates[0]
        aq_c = substrates[1]

        # Parameters
        kd = parameters[0]
        k1 = parameters[1]

        kminus1 = kd * k1

        rate = (k1 * a * b) - (kminus1 * c)

        return rate

class OxygenDiffusion(Reaction):

    def __init__(self,
                 kl=None, area=None, o2sat=None,
                 o2aq=None,
                 substrates=[], products=[]):

        super().__init__()

        self.reaction_substrate_names = [o2aq]
        self.parameter_names=[kl, area, o2sat]

        self.substrates = substrates
        self.products = products

    def calculate_rate(self, substrates, parameters):
        # Substrates
        o2aq = substrates[0]

        # Parameters
        kl = parameters[0]
        area = parameters[1]
        o2sat = parameters[2]

        rate = -kl * area * (o2aq - o2sat)

        return rate

class Flow(Reaction):
    def __init__(self,
                 flow_rate=None, column_volume=None,
                 input_substrates=[], substrates=[],
                 compartment_name=''):

        super().__init__()

        self.reaction_substrate_names = substrates
        self.parameter_names = [flow_rate, column_volume]

        self.substrates = substrates
        self.input_substrates = input_substrates
        self.input_substrates_indexes = []

        self.compartment_name = compartment_name

    def get_input_indexes(self, substrate_names):
        self.input_substrates_indexes = []
        for name in self.input_substrates:
            self.input_substrates_indexes.append(substrate_names.index(name))

    def reset_reaction(self):
        self.substrate_indexes = []
        self.input_substrates_indexes = []
        self.parameters = []

    def reaction(self, y, substrate_names, parameter_dict):
        if self.substrate_indexes == []:
            self.get_indexes(substrate_names)  # need to move this to the model

        if self.input_substrates_indexes == []:
            self.get_input_indexes(substrate_names)

        if self.run_model_parameters == []:
            self.run_model_parameters = self.get_parameters(parameter_dict)

        fr_over_cv = self.run_model_parameters[0] / self.run_model_parameters[1]

        y_prime = np.zeros(len(y))

        for index, input_index in zip(self.substrate_indexes, self.input_substrates_indexes):
            uM_current = y[index]
            uM_input = y[input_index]
            rate = fr_over_cv*(uM_input-uM_current)

            y_prime[index] += rate

        return y_prime


""" Modifiers (eg inhibtion) """
class Modifier():

    def __init__(self):
        self.substrate_names = []
        self.substrate_indexes = []

        self.parameter_names = []
        self.parameter_indexes = []

    def get_substrate_indexes(self, substrate_names):
        self.substrate_indexes = []
        for name in self.substrate_names:
            self.substrate_indexes.append(substrate_names.index(name))

    def get_parameter_indexes(self, parameter_names):
        self.parameter_indexes = []
        for name in self.parameter_names:
            self.parameter_indexes.append(parameter_names.index(name))

    def calc_modifier(self, substrates, parameters):
        # the substrate indexes will be stored in self.substrate_indexes,
        # in the order that they are named in self.substrate_names
        # same for parameters
        # use these indexes to write the equation here.

        return substrates, parameters

class SubstrateInhibition(Modifier):

    def __init__(self, ki=None, a=None):
        super().__init__()
        self.substrate_names = [a]
        self.parameter_names = [ki]

    def calc_modifier(self, substrates, parameters):
        ki = parameters[self.parameter_indexes[0]]
        a = substrates[self.substrate_indexes[0]]

        substrates[self.substrate_indexes[0]] = a * (1 + a / ki)

        return substrates, parameters

class CompetitiveInhibition(Modifier):

    def __init__(self, km=None, ki=None, i=None):
        super().__init__()
        self.substrate_names = [i]
        self.parameter_names = [km, ki]

    def calc_modifier(self, substrates, parameters):
        km = parameters[self.parameter_indexes[0]]
        ki = parameters[self.parameter_indexes[1]]
        i = substrates[self.substrate_indexes[0]]

        parameters[self.parameter_indexes[0]] = km * (1 + i/ki)

        return substrates, parameters

class MixedInhibition(Modifier):

    def __init__(self, kcat=None, km=None, ki=None, alpha=None, i=None):
        super().__init__()
        self.substrate_names = [i]
        self.parameter_names = [kcat, km, ki, alpha]

    def calc_modifier(self, substrates, parameters):
        kcat = parameters[self.parameter_indexes[0]]
        km = parameters[self.parameter_indexes[1]]
        ki = parameters[self.parameter_indexes[2]]
        alpha = parameters[self.parameter_indexes[3]]

        i = substrates[self.substrate_indexes[0]]

        parameters[self.parameter_indexes[0]] = kcat / (1 + i / (alpha * ki))
        parameters[self.parameter_indexes[1]] = km * (1 + i / ki) / (1 + i / (alpha * ki))

        return substrates, parameters

class FirstOrder_Modifier(Modifier):

    def __init__(self, kcat=None, k=None, s=None):
        super().__init__()
        self.substrate_names = [s]
        self.parameter_names = [kcat, k]

    def calc_modifier(self, substrates, parameters):
        kcat = parameters[self.parameter_indexes[0]]
        k = parameters[self.parameter_indexes[1]]
        s = substrates[self.substrate_indexes[0]]

        parameters[self.parameter_indexes[0]] = s*k*kcat

        return substrates, parameters






