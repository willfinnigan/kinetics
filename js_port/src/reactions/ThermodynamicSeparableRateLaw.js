const Reaction = require('../models/Reaction');

class Bi_Uni_sep_eq extends Reaction {
    constructor(kcat, kma, kmb, kmp, keq, a, b, p, enz, substrates = [], products = []) {
        super();
        this.parameter_names = [kcat, kma, kmb, kmp, keq];
        this.reaction_substrate_names = [a, b, p, enz];
        this.substrates = substrates;
        this.products = products;
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const a = substrates[0];
        const b = substrates[1];
        const p = substrates[2];
        const enz = substrates[3];

        // Parameters
        const kcat = parameters[0];
        const kma = parameters[1];
        const kmb = parameters[2];
        const kmp = parameters[3];
        const keq = parameters[4];

        const catalytic_capacity = enz * kcat;
        const thermodynamic_driving_force = 1 - (p / (a * b) / keq);
        const subs = (a / kma) * (b / kmb);
        const prods = (p / kmp);
        const substrate_saturation = subs / (subs + prods);

        return catalytic_capacity * substrate_saturation * thermodynamic_driving_force;
    }
}

class Bi_Bi_sep_eq extends Reaction {
    constructor(kcat, kma, kmb, kmp, kmq, keq, a, b, p, q, enz, substrates = [], products = []) {
        super();
        this.parameter_names = [kcat, kma, kmb, kmp, kmq, keq];
        this.reaction_substrate_names = [a, b, p, q, enz];
        this.substrates = substrates;
        this.products = products;
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const a = substrates[0];
        const b = substrates[1];
        const p = substrates[2];
        const q = substrates[3];
        const enz = substrates[4];

        // Parameters
        const kcat = parameters[0];
        const kma = parameters[1];
        const kmb = parameters[2];
        const kmp = parameters[3];
        const kmq = parameters[4];
        const keq = parameters[5];

        const catalytic_capacity = enz * kcat;
        const thermodynamic_driving_force = 1 - ((p * q) / (a * b) / keq);
        const subs = (a / kma) * (b / kmb);
        const prods = (p / kmp) * (q / kmq);
        const substrate_saturation = subs / (subs + prods);

        return catalytic_capacity * substrate_saturation * thermodynamic_driving_force;
    }
}

class Tri_Tri_seq_eq extends Reaction {
    constructor(kcat, kma, kmb, kmc, kmp, kmq, kmr, keq, a, b, c, p, q, r, enz, substrates = [], products = []) {
        super();
        this.parameter_names = [kcat, kma, kmb, kmc, kmp, kmq, kmr, keq];
        this.reaction_substrate_names = [a, b, c, p, q, r, enz];
        this.substrates = substrates;
        this.products = products;
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const a = substrates[0];
        const b = substrates[1];
        const c = substrates[2];
        const p = substrates[3];
        const q = substrates[4];
        const r = substrates[5];
        const enz = substrates[6];

        // Parameters
        const kcat = parameters[0];
        const kma = parameters[1];
        const kmb = parameters[2];
        const kmc = parameters[3];
        const kmp = parameters[4];
        const kmq = parameters[5];
        const kmr = parameters[6];
        const keq = parameters[7];

        const catalytic_capacity = enz * kcat;
        const thermodynamic_driving_force = 1 - ((p * q * r) / (a * b * c) / keq);
        const subs = (a / kma) * (b / kmb) * (c / kmc);
        const prods = (p / kmp) * (q / kmq) * (r / kmr);
        const substrate_saturation = subs / (subs + prods);

        return catalytic_capacity * substrate_saturation * thermodynamic_driving_force;
    }
}

module.exports = {
    Bi_Uni_sep_eq,
    Bi_Bi_sep_eq,
    Tri_Tri_seq_eq
};
