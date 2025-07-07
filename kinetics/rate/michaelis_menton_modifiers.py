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

class MixedInhibition2(Modifier):

    def __init__(self, kcat=None, km=None, kic=None, kiu=None, i=None):
        super().__init__()
        self.substrate_names = [i]
        self.parameter_names = [kcat, km, kic, kiu]

    def calc_modifier(self, substrates, parameters):
        kcat = parameters[self.parameter_indexes[0]]
        km = parameters[self.parameter_indexes[1]]
        kic = parameters[self.parameter_indexes[2]]
        kiu = parameters[self.parameter_indexes[3]]

        i = substrates[self.substrate_indexes[0]]

        parameters[self.parameter_indexes[0]] = kcat / (1 + i / kiu)
        parameters[self.parameter_indexes[1]] = km * (1 + i / kic) / (1 + i / (kiu))

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
