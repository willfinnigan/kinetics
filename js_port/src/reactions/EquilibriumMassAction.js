const Reaction = require('../models/Reaction');

class Uni_mass_action_eq extends Reaction {
    constructor(keq, kf, a, p, substrates = [], products = []) {
        super();
        this.parameter_names = [keq, kf];
        this.reaction_substrate_names = [a, p];
        this.substrates = substrates;
        this.products = products;
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const a = substrates[0];
        const p = substrates[1];

        // Parameters
        const keq = parameters[0];
        const kf = parameters[1];

        const catalytic_capacity = kf;

        if (a === 0 || p === 0) {
            return 0.0;
        }

        const thermodynamic_driving_force = 1 - (p / a / keq);

        return catalytic_capacity * thermodynamic_driving_force;
    }
}

module.exports = {
    Uni_mass_action_eq
};
