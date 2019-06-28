from deap import creator, base, tools, algorithms
import random
from tqdm import tqdm
import matplotlib.pyplot as plt

class GA_Base_Class(object):

    def __init__(self, model=None, metrics=None, weights=(1,), bounds={}):

        self.names_list = []
        self.bounds_list = []
        self.bounds_dict = bounds

        self.indpb_mate = 0.5
        self.mu=0
        self.sigma=0.4
        self.indpb_mutate = 0.5

        self.cxpb = 0.5
        self.mutpb = 0.5

        self.toolbox = base.Toolbox()

        self.model = model

        self.metrics=metrics

        self.weights = weights

        self.all_pops = []

        self.initial_pop_size = 100
        self.generations = 10
        self.num_to_select = 75
        self.num_children = 75

        self.logging=True

        self.flow = False

    def set_ga_settings(self, indpb_mate=0.5, mu=0, sigma=0.4, indpb_mutate=0.5):
        self.indpb_mate = indpb_mate
        self.mu=mu
        self.sigma=sigma
        self.indpb_mutate = indpb_mutate

    def setup(self):
        self.names_list = list(self.bounds_dict.keys())
        self.bounds_list = list(self.bounds_dict.values())

        creator.create("FitnessMax", base.Fitness, weights=self.weights)
        creator.create("Individual", list, fitness=creator.FitnessMax)

        self.toolbox.register("make_ind", self.make_ind, self.bounds_list)
        self.toolbox.register("individual", tools.initIterate, creator.Individual, self.toolbox.make_ind)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)
        self.toolbox.register("mate", tools.cxUniform, indpb=self.indpb_mate)

        self.toolbox.register("select", tools.selNSGA2)
        self.toolbox.register("mutate", tools.mutGaussian, mu=self.mu, sigma=self.sigma, indpb=self.indpb_mutate)
        self.toolbox.register("evaluate", self.evaluate)

    def make_ind(self,bounds):
        # bounds = [(0,100), (3, 4), (12, 15)...]

        ind = []
        for tuple in bounds:
            selection = random.uniform(tuple[0], tuple[1])
            ind.append(selection)

        return ind

    def check_neg_substrate(self, ind):
        for substrate in ind:
            if substrate <= 0:
                return True

        return False

    def low_fitness(self):
        fit = []
        for weight in self.weights:
            fit.append(-99999 * weight)

        return fit

    def update_model_for_evaluation(self, ind):
        # Update model species with the concentrations in ind
        for i in range(len(self.names_list)):
            name = self.names_list[i]
            if name == 'Time':
                self.model.set_end_time(ind[i])
            elif 'Parameter_' in name:
                parameter_name = name[10:]
                self.model.parameters[parameter_name] = ind[i]
                self.model.run_model_parameters[parameter_name] = ind[i]
            else:
                self.model.species[name] = ind[i]

        self.metrics.refresh_metrics(model=self.model, flow_rate=self.flow)

    def evaluate(self, ind):

        # Check that the GA hasn't evolved towards negative substrate
        if self.check_neg_substrate(ind) == True:
                return self.low_fitness()

        self.update_model_for_evaluation(ind)

        # Calculate fitness
        fitness = self.fitness()

        return fitness

    def fitness(self):
        return 1

    def run_ga(self, initial_pop=False, plot=False):
        if initial_pop != False:
            population = initial_pop
        else:
            population = self.toolbox.population(n=self.initial_pop_size)

        fitnesses = list(map(self.toolbox.evaluate, population))
        for ind, fit in zip(population, fitnesses):
            ind.fitness.values = fit

        self.all_pops.append(population)

        for n in tqdm(range(self.generations)):
            self.all_pops.append(population)
            offspring = algorithms.varOr(population, self.toolbox, cxpb=self.cxpb, mutpb=self.mutpb, lambda_=self.num_children)

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = map(self.toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit

            new_population = population + offspring
            population = self.toolbox.select(new_population, k=self.num_to_select)

            if plot==True:
                self.model.plot_substrate(self.metrics.substrate)
                self.model.plot_substrate(self.metrics.product)
                plt.show()


