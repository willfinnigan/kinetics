class Reaction {
    constructor() {
        // These are set by the user
        this.parameters = {};
        this.parameterDistributions = {};

        // indexes used to access values during model run
        this.substrateIndexes = [];
        this.parameterIndexes = [];

        // These are set when the reaction is set up
        this.reactionSubstrateNames = [];
        this.parameterNames = [];
        this.substrates = [];
        this.products = [];

        // These are added as needed
        this.modifiers = [];
        this.checkPositive = false;
        this.checkLimitsFunctions = [];
    }

    setParameterDefaultsToMean() {
        for (const name in this.parameterDistributions) {
            if (!this.parameters.hasOwnProperty(name)) {
                if (Array.isArray(this.parameterDistributions[name])) {
                    this.parameters[name] = (this.parameterDistributions[name][0] + this.parameterDistributions[name][1]) / 2;
                } else {
                    // Assuming the distribution object has a mean() method
                    this.parameters[name] = this.parameterDistributions[name].mean();
                }
            }
        }
    }

    setupReaction(speciesNames, parameterNames) {
        // get indexes for rate calculation
        this.substrateIndexes = this.reactionSubstrateNames.map(name => speciesNames.indexOf(name));
        this.parameterIndexes = this.parameterNames.map(name => parameterNames.indexOf(name));

        // set up modifiers
        for (const modifier of this.modifiers) {
            modifier.getSubstrateIndexes(this.reactionSubstrateNames);
            modifier.getParameterIndexes(this.parameterNames);
        }
    }

    addModifier(modifier) {
        for (const name of modifier.parameterNames) {
            if (!this.parameterNames.includes(name)) {
                this.parameterNames.push(name);
            }
        }

        for (const name of modifier.substrateNames) {
            if (!this.reactionSubstrateNames.includes(name)) {
                this.reactionSubstrateNames.push(name);
            }
        }

        this.modifiers.push(modifier);
    }

    calculateRate(substrates, parameters) {
        return 0;
    }

    reaction(y, substrateNames, parameterValues) {
        // Get the substrates from y using the substrate indexes
        const substrates = this.substrateIndexes.map(index => y.get(index));

        // Get the parameters using the parameter indexes
        let parameters = this.parameterIndexes.map(index => parameterValues[index]);

        // calculate the effects of any modifiers
        for (const modifier of this.modifiers) {
            [substrates, parameters] = modifier.calcModifier(substrates, parameters);
        }

        // calculate the rate (this function is modified by the user)
        const rate = this.calculateRate(substrates, parameters);

        // calculate the change in substrate concentrations (y_prime)
        const substrateIndices = this.substrates.map(name => substrateNames.indexOf(name));
        const productIndices = this.products.map(name => substrateNames.indexOf(name));

        const yPrime = new Array(y.shape[0]).fill(0);

        for (const index of substrateIndices) {
            yPrime[index] -= rate;
        }

        for (const index of productIndices) {
            yPrime[index] += rate;
        }


        let finalYPrime = this.modifyProduct(yPrime, substrateNames);

        if (this.checkPositive) {
            finalYPrime = finalYPrime.map(value => Math.max(0, value));
        }

        return finalYPrime;
    }

    modifyProduct(yPrime, substrateNames) {
        return yPrime;
    }

    samplingLimits(parameterDict) {
        for (const func of this.checkLimitsFunctions) {
            if (!func(parameterDict)) {
                return false;
            }
        }
        return true;
    }
}

class Model {
    constructor() {
        this._reactions = [];
        this.timeseries = null;
        this.y = [];
        this.ts = null;
    }

    set_time(start, end, steps, mode = 'linear') {
        if (mode === 'linear') {
            const stepSize = (end - start) / (steps - 1);
            this.ts = Array.from({
                length: steps
            }, (_, i) => start + i * stepSize);
        } else if (mode === 'log') {
            const startLog = Math.log10(start);
            const endLog = Math.log10(end);
            const stepSize = (endLog - startLog) / (steps - 1);
            this.ts = Array.from({
                length: steps
            }, (_, i) => Math.pow(10, startLog + i * stepSize));
        } else {
            throw new Error(`Unknown time series mode: ${mode}. Use 'linear' or 'log'.`);
        }
    }

    add_reaction(reaction) {
        this._reactions.push(reaction);
    }

    _parameters_and_species_from_reactions() {
        const species = {};
        const parameters = {};

        for (const reaction of this._reactions) {
            reaction.setParameterDefaultsToMean();

            for (const name in reaction.parameters) {
                if (!parameters.hasOwnProperty(name)) {
                    parameters[name] = reaction.parameters[name];
                } else {
                    throw new Error(`Parameter ${name} already set in model, can't overwrite`);
                }
            }

            for (const substrate of [...reaction.substrates, ...reaction.products, ...reaction.reaction_substrate_names]) {
                if (!species.hasOwnProperty(substrate)) {
                    species[substrate] = 0;
                }
            }
        }

        return [species, parameters];
    }

    get_parameter_distributions() {
        const parameterDistributions = {};
        for (const reaction of this._reactions) {
            for (const name in reaction.parameterDistributions) {
                if (parameterDistributions.hasOwnProperty(name)) {
                    throw new Error(`Parameter ${name} already set in model, can't overwrite`);
                }
                parameterDistributions[name] = reaction.parameterDistributions[name];
            }
        }
        return parameterDistributions;
    }

    _setup_model(species_names, parameter_names) {
        for (const reaction of this._reactions) {
            reaction.setupReaction(species_names, parameter_names);
        }
    }

    _set_default_species(species) {
        for (const name in species) {
            if (species[name] && typeof species[name].rvs === 'function') {
                species[name] = species[name].mean();
            }
            if (Array.isArray(species[name]) && species[name].length === 2) {
                species[name] = (species[name][0] + species[name][1]) / 2;
            }
        }
        return species;
    }

    run_single(starting_concentrations, solver = new Solver()) {
        let [species, parameters] = this._parameters_and_species_from_reactions();

        starting_concentrations = this._set_default_species(starting_concentrations);
        Object.assign(species, starting_concentrations);

        const species_names = Object.keys(species);
        const species_values = Object.values(species);
        const parameter_names = Object.keys(parameters);
        const parameter_values = Object.values(parameters);

        this._setup_model(species_names, parameter_names);

        const y = solver.run(this._reactions, species_names, species_values, parameter_values, this.ts);

        return {
            y: y,
            species_names: species_names
        };
    }

    run_multi(starting_concentrations, sampler, solver = new Solver()) {
        if (!sampler) {
            throw new Error("Sampler must be provided for multi-parameter simulations.");
        }

        let [species, parameters] = this._parameters_and_species_from_reactions();

        this._setup_model(Object.keys(species), Object.keys(parameters));

        Object.assign(species, starting_concentrations);

        const parameter_distributions = this.get_parameter_distributions();
        const samples = sampler.sample(parameter_distributions, species);

        const results = [];
        for (const [parameter_dict, species_dict] of samples) {
            const current_species = { ...species,
                ...species_dict
            };
            const current_parameters = { ...parameters,
                ...parameter_dict
            };

            const species_names = Object.keys(current_species);
            const species_values = Object.values(current_species);
            const parameter_names = Object.keys(current_parameters);
            const parameter_values = Object.values(current_parameters);

            const y = solver.run(this._reactions, species_names, species_values, parameter_values, this.ts);
            results.push(y);
        }

        return {
            results: results,
            species_names: Object.keys(species)
        };
    }

    plot(species_name, result, options = {}) {
        const plotly = require('plotly.js');

        const species_index = result.species_names.indexOf(species_name);
        if (species_index === -1) {
            throw new Error(`Species ${species_name} not found in model`);
        }

        const x = this.ts;
        const y = result.y.map(row => row[species_index]);

        const trace = {
            x: x,
            y: y,
            mode: 'lines',
            name: species_name
        };

        const layout = {
            title: options.title || `Concentration of ${species_name} over time`,
            xaxis: {
                title: options.xlabel || 'Time'
            },
            yaxis: {
                title: options.ylabel || 'Concentration'
            }
        };

        const figure = {
            data: [trace],
            layout: layout
        };

        // This will create a div with the plot, which can be then saved to a file or displayed in a browser.
        // Since I cannot run the code, I will just return the figure object.
        return figure;
    }
}

class Modifier {
    constructor() {
        this.substrate_names = [];
        this.substrate_indexes = [];

        this.parameter_names = [];
        this.parameter_indexes = [];
    }

    getSubstrateIndexes(substrate_names) {
        this.substrate_indexes = [];
        for (const name of this.substrate_names) {
            this.substrate_indexes.push(substrate_names.indexOf(name));
        }
    }

    getParameterIndexes(parameter_names) {
        this.parameter_indexes = [];
        for (const name of this.parameter_names) {
            this.parameter_indexes.push(parameter_names.indexOf(name));
        }
    }

    calc_modifier(substrates, parameters) {
        return [substrates, parameters];
    }
}

class Uni extends Reaction {
    constructor(kcat, kma, a, enz, substrates = [], products = []) {
        super();
        this.reaction_substrate_names = [a, enz];
        this.parameter_names = [kcat, kma];
        this.substrates = substrates;
        this.products = products;
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const a = substrates[0];
        const enz = substrates[1];

        // Parameters
        const kcat = parameters[0];
        const kma = parameters[1];

        const rate = kcat * enz * (a / (kma + a));

        return rate;
    }
}

class Solver {
    constructor() {}

    run(reactions, species_names, species_values, parameter_values, time) {
        const ndarray = require('ndarray');
        const rk4 = require('ode-rk4');
        const y0 = ndarray(new Float64Array(species_values));

        const deriv = (dydt, y, t) => {
            for (let i = 0; i < dydt.shape[0]; i++) {
                dydt.set(i, 0);
            }
            for (const reaction of reactions) {
                const y_prime = reaction.reaction(y, species_names, parameter_values);
                for (let i = 0; i < y_prime.length; i++) {
                    dydt.set(i, dydt.get(i) + y_prime[i]);
                }
            }
        };

        const integrator = rk4(y0, deriv, 0, time[1] - time[0]);

        const y = [y0.data.slice()];
        for (let i = 1; i < time.length; i++) {
            integrator.step();
            y.push(integrator.y.data.slice());
        }

        return y;
    }
}

module.exports = {
    Model,
    Reaction,
    Modifier,
    Uni,
    Solver
};
