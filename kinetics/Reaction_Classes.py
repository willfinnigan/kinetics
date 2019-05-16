from kinetics import Equations
from kinetics.Model import calculate_yprime, check_positive
import numpy as np
import copy

class Reaction():

    def __init__(self):

        self.reaction_substrate_names = []
        self.substrate_indexes = []
        self.substrates = []
        self.products = []

        self.parameter_names = []
        self.parameter_defaults = {}
        self.parameter_bounds = {}

        self.parameters = []

        self.modifiers = []

        self.check_positive = False

    def set_substrates_and_products(self, substrates, products):
        self.substrates = substrates
        self.products = products

    def set_parameters(self, parameter_defaults={}, parameter_bounds={}):
        self.parameter_defaults = parameter_defaults
        self.parameter_bounds = parameter_bounds

    def set_parameter_defaults_to_mean_of_bounds(self):
        for name in self.parameter_bounds:
            lower = self.parameter_bounds[name][0]
            upper = self.parameter_bounds[name][1]
            mean_value = (lower + upper) / 2
            self.parameter_defaults[name] = mean_value

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
        self.parameters = []

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

        if self.parameters == []:
            self.parameters = self.get_parameters(parameter_dict)

        for modifier in self.modifiers:
            if modifier.substrate_indexes == []:
                modifier.get_substrate_indexes(self.reaction_substrate_names)
            if modifier.parameter_indexes == []:
                modifier.get_parameter_indexes(self.parameter_names)

        substrates = self.get_substrates(y)

        substrates, parameters = self.calculate_modifiers(substrates, copy.copy(self.parameters))

        rate = self.calculate_rate(substrates, parameters)

        y_prime = calculate_yprime(y, rate, self.substrates, self.products, substrate_names)
        y_prime = self.modify_product(y_prime, substrate_names)

        if self.check_positive == True:
            y_prime = check_positive(y_prime)

        return y_prime

    def modify_product(self, y_prime, substrate_names):
        return y_prime

""" Michaelis-Menten irreversible equations """
class One_irr(Reaction):

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

class Two_irr(Reaction):

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

class Two_ternary_complex_irr(Reaction):

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

class Two_ping_pong_irr(Reaction):

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

        rate = Equations.two_substrate_pingpong_irr(kcat=kcat, kma=kma, kmb=kmb, enz=enz, a=a, b=b)

        return rate

class Three_seq_irr_redam(Reaction):
    # This is the mechanism RedAms use

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

        rate = Equations.three_substrate_irreversible_sequential(kcat=kcat, kma=kma, kmb=kmb, kmc=kmc,
                                                                 kia=kia, kib=kib,
                                                                 enz=enz, a=a, b=b, c=c)
        return rate

class Three_seq_irr_car(Reaction):
    # This is the mechanism CARs use

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

        rate = Equations.three_substrate_irreversible_ter_ordered(kcat=kcat,
                                                                  kma=kma, kmb=kmb, kia=kia, kmc=kmc,
                                                                  enz=enz, a=a, b=b, c=c)
        return rate

class Two_ternary_complex_small_kma(Reaction):

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
class Two_Ordered_rev(Reaction):

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

        rate = Equations.o2_diffusion(kl=kl, area=area, o2sat=o2sat, o2aq=o2aq)
        return rate

class Flow(Reaction):

    def __init__(self,
                 flow_rate=None, column_volume=None,
                 input_substrates=[], substrates=[],
                 compartment_name=''):

        super().__init__()

        self.reaction_substrate_names = substrates
        self.parameter_names=[flow_rate, column_volume]

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
            self.get_indexes(substrate_names) # need to move this to the model

        if self.input_substrates_indexes == []:
            self.get_input_indexes(substrate_names)

        if self.parameters == []:
            self.parameters = self.get_parameters(parameter_dict)

        # Parameters (convert to L from ml)
        flow_rate = self.parameters[0] / 1000
        column_volume = self.parameters[1] / 1000

        def calculate_uM_per_min(uM_initial, uM_input, column_volume, flow_rate):
            umols_initial = uM_initial * column_volume
            umols_in_flow_leaving = uM_initial * flow_rate
            umols_in_flow_entering = uM_input * flow_rate

            umols_now_in_column = umols_initial - umols_in_flow_leaving + umols_in_flow_entering
            uM_now_in_column = umols_now_in_column / column_volume

            uM_change_per_min = uM_now_in_column - uM_initial

            return uM_change_per_min

        y_prime = np.zeros(len(y))

        for index, input_index in zip(self.substrate_indexes, self.input_substrates_indexes):
            uM_initial = y[index]
            uM_input = y[input_index]

            rate = calculate_uM_per_min(uM_initial, uM_input, column_volume, flow_rate)
            y_prime[index] += rate

        return y_prime



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
