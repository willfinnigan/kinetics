from SALib.sample import latin, saltelli

def parse_samples(samples, parameter_names, species_names):
    """
    Takes a list of samples and converts this to a list of tuples of dictionaries which can update the Model class.

    Samples are created using the relevent function (eg LHC or Saltelli)
    Samples will be a list.
    Each entry in the list, is a list of values for the species and parameters,
    in the order they were defined in the 'problem' dictionary.
    This will be the parameters first followed by the species.

    For example, Samples = [ [value1, value2, value3..],  [value1, value2, value3], ...]

    The ordered lists of parameter and species names are used to unpack these lists back into dictionaries.
    These dictionaries can be used to update the model.
    The dictionaries are saved in a new list called parsed samples

    Parsed samples = [ (parameter_dict1, species_dict1), (parameter_dict2, species_dict2) ..]


    :param samples:  a list of lists containing the samples values, in the order they were defined in problem dict.
    :param parameter_names: ordered list of parameter names
    :param species_names: ordered list of species names
    :return: Parsed samples - [ (parameter_dict1, species_dict1), (parameter_dict2, species_dict2) ..]
    """

    parsed_samples = []

    def unpack_sampled_params_and_species(sample_list, parameter_names_list, species_names_list):
        """
        Returns a dictionaries for parameters and species

        :param sample_list: a list of sample values,
                            parameters first then species, of the same order as the names lists
        :param parameter_names_list: an ordered list of parameter_names which matches the first part of sample_list
        :param species_names_list:  an ordered list of species_names which matches the first part of sample_list
        :return: (parameters_dictionary, species_dictionary)   - dictionary format is {'name' : value, ..}
        """

        parameters_dict = {}
        count = 0
        for i in range(len(parameter_names_list)):
            name = parameter_names_list[i]
            parameters_dict[name] = sample_list[i]
            count += 1

        species_dict = {}
        for i in range(count, len(species_names_list) + count):
            name = species_names_list[i - count]
            species_dict[name] = sample_list[i]

        return parameters_dict, species_dict

    for sampled in samples:
        parameters, species = unpack_sampled_params_and_species(sampled, parameter_names, species_names)
        parsed_samples.append([parameters, species])

    # returns a list of tuples containing [(parameter_dict, species_dict), ] for each sample
    return parsed_samples

def check_not_neg(sample, name, negative_allowed):
    if sample == None:
        return False

    if (sample < 0) and (name not in negative_allowed):
        return False

    return True


""" -- Generate Samples --"""
def sample_distributions(model, num_samples=1000, negative_allowed=[]):
    """
    Makes a set of samples from the species and parameter distributions in the model.

    Args:
        model (kinetics.Model): A model object
        num_samples (int): Number of samples to make (default 1000)
        negative_allowed (list): A list of any distributions that can be negative.

    Returns:
        A list of samples.  Each entry in the list is a tuple containing (parameter_dict, species_dict) for the samples.
    """

    samples = []
    for i in range(num_samples):
        parameter_dict = {}
        species_dict = {}

        for name, distribution in model.parameter_distributions.items():
            sample = None
            while not check_not_neg(sample, name, negative_allowed):
                sample = distribution.rvs()
            parameter_dict[name] = sample

        for name, distribution in model.species_distributions.items():
            sample=None
            while not check_not_neg(sample, name, negative_allowed):
                sample = distribution.rvs()
            species_dict[name] = sample

        samples.append([parameter_dict, species_dict])

    return samples # samples will = [ (parameter_dict1, species_dict1), (parameter_dict2, species_dict2) ..]

def sample_uniforms(model, num_samples=1000):

    bounds = []
    names = []
    for name, tuple in model.parameter_distributions.items():
        bounds.append(tuple)
        names.append(names)

    for name, tuple in  model.species_distributions.items():
        bounds.append(tuple)
        names.append(names)

    problem = {'num_vars': len(names),
               'names': names,
               'bounds': bounds}

    lhc_samples = latin.sample(problem, num_samples)

    samples = parse_samples(lhc_samples, list(model.parameter_distributions.keys()), list(model.species_distributions.keys()))

    return samples

def salib_problem_with_bounds(model, negative_allowed=[], ppf=(0,1)):
    """
    Make a salib problem by specifying bounds using ppf of scipy distributions

    Args:
        model (Model): The model object
        negative_allowed (list): Any distributions which are allowed negative samples
        ppf (Tuple): Percent point functions to take.  Default is 0 and 1 which will give the absolute upper and lower bounds of a distribution.

    Returns:
        An SALib problem
    """

    names = list(model.parameter_distributions.keys()) + list(model.species_distributions.keys())

    bounds = []

    for name, distribution in model.parameter_distributions.items():
        upper = distribution.ppf(ppf[1])
        lower = None
        lower_ppf = ppf[0]
        while not check_not_neg(lower, name, negative_allowed):
            lower = distribution.ppf(lower_ppf)
            lower_ppf += 0.01
        bounds.append([lower, upper])

    for name, distribution in model.species_distributions.items():
        upper = distribution.ppf(ppf[1])
        lower = None
        lower_ppf = ppf[0]
        while not check_not_neg(lower, name, negative_allowed):
            lower = distribution.ppf(lower_ppf)
            lower_ppf += 0.01
        bounds.append([lower, upper])

    problem = {'num_vars': len(names),
               'names': names,
               'bounds': bounds}

    return problem

def make_saltelli_samples(model, salib_problem, num_samples, second_order=False):
    """
    Use SALib to make saltelli samples

    Args:
        model (Model): The model object
        salib_problem (dict): An SALib Problem
        num_samples (int): number of samples to take
        second_order (bool): look at second order interactions

    Returns:
        Samples from salib which have been parsed into a set of samples which run_all_models can take.

    """

    saltelli_samples = saltelli.sample(salib_problem, num_samples, calc_second_order=second_order)

    samples = parse_samples(saltelli_samples, list(model.parameter_distributions.keys()), list(model.species_distributions.keys()))

    return samples