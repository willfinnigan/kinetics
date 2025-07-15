const Reaction = require('../models/Reaction');

class Generic extends Reaction {
    constructor(params = [], species = [], rate_equation = '', substrates = [], products = []) {
        super();
        this.reaction_substrate_names = species;
        this.parameter_names = params;
        this.rate_equation = rate_equation;
        this.substrates = substrates;
        this.products = products;
        this.rate_function = this._create_rate_function();
    }

    _create_rate_function() {
        const param_names = this.parameter_names;
        const substrate_names = this.reaction_substrate_names;
        const all_names = [...param_names, ...substrate_names];
        return new Function(...all_names, `return ${this.rate_equation}`);
    }


    calculateRate(substrates, parameters) {
        const all_values = [...parameters, ...substrates];
        return this.rate_function(...all_values);
    }
}

module.exports = {
    Generic
};
