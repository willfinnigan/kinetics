import numpy as np
from scipy import integrate
import datetime
import copy

"""
The model class, which derives from a list.  
The list is a list of the reaction functions.  
When the model is run each reaction is modelled in turn and the output added to yprime"""


class Model(list):
    def __init__(self):
        super(Model, self).__init__()

        self.start = 0
        self.end = 100
        self.steps = 100
        self.mxsteps = 10000

        self.time = np.linspace(self.start, self.end, self.steps)

        self.species_defaults = {}
        self.species = {}
        self.species_names = []
        self.species_starting_values = []

        self.parameter_defaults = {}
        self.parameters = {}

    def set_time(self, start, end, steps, mxsteps=10000):
        self.start = start
        self.end = end
        self.steps = steps
        self.mxsteps = mxsteps

        self.time = np.linspace(self.start, self.end, self.steps)

    def set_parameter_defaults(self, parameters):
        self.parameters = parameters
        self.parameter_defaults = copy.deepcopy(parameters)


    def update_parameters(self, parameters):
        self.parameters.update(parameters)

    def set_species_defaults(self, species_defaults):
        self.species_defaults = copy.deepcopy(species_defaults)
        self.species = species_defaults
        self.species_names, self.species_starting_values = get_species_positions(self.species)


    def update_species(self, species_dict):
        self.species.update(species_dict)
        self.species_names, self.species_starting_values = get_species_positions(self.species)

    def deriv(self, y, t):
        yprime = 0

        for reaction in self:
            yprime += reaction(y, self.species_names, self.parameters)

        return yprime

    def run_model(self):
        y0 = np.array(self.species_starting_values)
        y = integrate.odeint(self.deriv, y0, self.time, mxstep=self.mxsteps)

        return y

    def reset_model(self):
        self.species = self.species_defaults
        self.species_names, self.species_starting_values = get_species_positions(self.species)
        self.parameters = self.parameter_defaults


"""Functions for formatting species and parameters dicts to the correct format"""


def get_species_positions(species):
    species_names = []
    species_starting_values = []

    for name in species:
        species_names.append(name)
        if type(species[name]) == list:
            species_starting_values.append(species[name][0])
        else:
            species_starting_values.append(species[name])

    return species_names, species_starting_values


def set_parameter_defaults(parameters_with_error):
    parameters = {}
    for name in parameters_with_error:
        parameters[name] = parameters_with_error[name][0]

    return parameters


def set_species_defaults(species_with_error):
    species = {}
    for name in species_with_error:
        if type(species_with_error[name]) == list or type(species_with_error[name]) == tuple:
            species[name] = species_with_error[name][0]
        else:
            species[name] = species_with_error[name]

    return species


"""Functions to add or substract the rate from yprime at the correct index's"""


def yprime_plus(y_prime, rate, substrates, s_names):
    for name in substrates:
        y_prime[s_names.index(name)] += rate

    return y_prime

def yprime_minus(y_prime, rate, substrates, s_names):
    for name in substrates:
        y_prime[s_names.index(name)] -= rate

    return y_prime

def calculate_yprime(y, rate, substrates, products, substrate_names):
    y_prime = np.zeros(len(y))

    for name in substrates:
        y_prime[substrate_names.index(name)] -= rate

    for name in products:
        y_prime[substrate_names.index(name)] += rate

    return y_prime

"""Functions to output basic model output (with only 1 set of y values)"""

def print_model_output(y, time, species_names, names_to_output):
    print("Time" + ", ", end="")

    for name in names_to_output:
        print(str(name) + ", ", end="")
    print()

    for i in range(len(time)):
        print(str(time[i]) + ", ", end="")

        for name in names_to_output:
            print(str(y[i][species_names.index(name)]) + ", ", end="")
        print()

def save_model_ouput(y, model, names_to_output, filename=''):
    import os
    directory = os.path.dirname('output/')

    try:
        os.stat(directory)
    except:
        os.mkdir(directory)

    if filename == '':
        filename = os.path.basename(__file__)
        filename = filename[0:-3] + '.txt'

    filename = directory + "/" + filename

    file = open(filename, "w")

    file.write("Date: " + str(datetime.datetime.now()) + str("\n"))
    file.write("\n")

    file.write("---Species concentrations--- \n")
    for i in range(len(model.species_names)):
        file.write(str(model.species_names[i]) + " : " + str(model.species_starting_values[i]) + "\n")

    file.write("\n")

    file.write("--- Parameters --- \n")
    for name in model.parameters:
        file.write(name + " : " + str(model.parameters[name]) + "\n")

    file.write("\n")

    file.write("Time" + ", ")

    for name in names_to_output:
        file.write(str(name) + ", ")
    file.write("\n")

    for i in range(len(model.time)):
        file.write(str(model.time[i]) + ", ")

        for name in names_to_output:
            file.write(str(y[i][model.species_names.index(name)]) + ", ")
        file.write("\n")

    file.close()
