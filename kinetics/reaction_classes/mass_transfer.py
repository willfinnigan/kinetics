from kinetics.reaction_classes.reaction_base_class import Reaction
import numpy as np

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

        rate = (k1*org_c) - (kminus1*aq_c)

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

        self.parameters = {}
        self.parameter_distributions = {}

    def get_input_indexes(self, substrate_names):
        self.input_substrates_indexes = []
        for name in self.input_substrates:
            self.input_substrates_indexes.append(substrate_names.index(name))

    def reset_reaction(self):
        self.substrate_indexes = []
        self.input_substrates_indexes = []
        self.run_model_parameters = []


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