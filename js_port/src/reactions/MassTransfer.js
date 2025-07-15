const Reaction = require('../models/Reaction');

class FirstOrderRate extends Reaction {
    constructor(k, a, substrates = [], products = []) {
        super();
        this.reaction_substrate_names = [a];
        this.parameter_names = [k];
        this.substrates = substrates;
        this.products = products;
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const a = substrates[0];

        // Parameters
        const k = parameters[0];

        return k * a;
    }
}

class Binding extends Reaction {
    constructor(k1, kminus1, a, b, c, substrates = [], products = []) {
        super();
        this.reaction_substrate_names = [a, b, c];
        this.parameter_names = [k1, kminus1];
        this.substrates = substrates;
        this.products = products;
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const a = substrates[0];
        const b = substrates[1];
        const c = substrates[2];

        // Parameters
        const k1 = parameters[0];
        const kminus1 = parameters[1];

        const rate = (k1 * a * b) - (kminus1 * c);

        return rate;
    }
}

class BiSecondOrderRate extends Reaction {
    constructor(k, a, b, substrates = [], products = []) {
        super();
        this.reaction_substrate_names = [a, b];
        this.parameter_names = [k];
        this.substrates = substrates;
        this.products = products;
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const a = substrates[0];
        const b = substrates[1];

        // Parameters
        const k = parameters[0];

        return k * a * b;
    }
}

class Binding_kd extends Reaction {
    constructor(kd, k1, a, b, c, substrates = [], products = []) {
        super();
        this.reaction_substrate_names = [a, b, c];
        this.parameter_names = [kd, k1];
        this.substrates = substrates;
        this.products = products;
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const a = substrates[0];
        const b = substrates[1];
        const c = substrates[2];

        // Parameters
        const kd = parameters[0];
        const k1 = parameters[1];

        const kminus1 = kd * k1;

        const rate = (k1 * a * b) - (kminus1 * c);

        return rate;
    }
}

class DiffusionEquilibrium extends Reaction {
    constructor(kd, k1, org_c, aq_c) {
        super();
        this.reaction_substrate_names = [org_c, aq_c];
        this.parameter_names = [kd, k1];
        this.substrates = [org_c];
        this.products = [aq_c];
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const org_c = substrates[0];
        const aq_c = substrates[1];

        // Parameters
        const kd = parameters[0];
        const k1 = parameters[1];

        const kminus1 = kd * k1;

        const rate = (k1 * org_c) - (kminus1 * aq_c);

        return rate;
    }
}

class OxygenDiffusion extends Reaction {
    constructor(kl, area, o2sat, o2aq, substrates = [], products = []) {
        super();
        this.reaction_substrate_names = [o2aq];
        this.parameter_names = [kl, area, o2sat];
        this.substrates = substrates;
        this.products = products;
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const o2aq = substrates[0];

        // Parameters
        const kl = parameters[0];
        const area = parameters[1];
        const o2sat = parameters[2];

        const rate = -kl * area * (o2aq - o2sat);

        return rate;
    }
}

class Flow extends Reaction {
    constructor(flow_rate, column_volume, input_substrates = [], substrates = [], compartment_name = '') {
        super();
        this.reaction_substrate_names = substrates;
        this.parameter_names = [flow_rate, column_volume];
        this.substrates = substrates;
        this.input_substrates = input_substrates;
        this.input_substrates_indexes = [];
        this.compartment_name = compartment_name;
        this.parameters = {};
        this.parameter_distributions = {};
    }

    getInputIndexes(substrate_names) {
        this.input_substrates_indexes = this.input_substrates.map(name => substrate_names.indexOf(name));
    }

    reaction(y, substrate_names, parameter_values) {
        if (this.substrate_indexes.length === 0) {
            this.setupReaction(substrate_names, this.parameter_names);
        }
        if (this.input_substrates_indexes.length === 0) {
            this.getInputIndexes(substrate_names);
        }

        const fr_over_cv = parameter_values[this.parameter_indexes[0]] / parameter_values[this.parameter_indexes[1]];

        const y_prime = new Array(y.shape[0]).fill(0);

        for (let i = 0; i < this.substrate_indexes.length; i++) {
            const index = this.substrate_indexes[i];
            const input_index = this.input_substrates_indexes[i];
            const uM_current = y.get(index);
            const uM_input = y.get(input_index);
            const rate = fr_over_cv * (uM_input - uM_current);
            y_prime[index] += rate;
        }

        return y_prime;
    }
}


module.exports = {
    FirstOrderRate,
    Binding,
    BiSecondOrderRate,
    Binding_kd,
    DiffusionEquilibrium,
    OxygenDiffusion,
    Flow
};
