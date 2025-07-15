const Solver = require('../solvers/Solver');

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
}

module.exports = Model;
