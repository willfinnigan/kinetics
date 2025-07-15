const Reaction = require('../models/Reaction');

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

class Bi extends Reaction {
    constructor(kcat, kma, kmb, a, b, enz, substrates = [], products = []) {
        super();
        this.reaction_substrate_names = [a, b, enz];
        this.parameter_names = [kcat, kma, kmb];
        this.substrates = substrates;
        this.products = products;
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const a = substrates[0];
        const b = substrates[1];
        const enz = substrates[2];

        // Parameters
        const kcat = parameters[0];
        const kma = parameters[1];
        const kmb = parameters[2];

        // Rate equation
        const rate = (kcat * enz) * (a / (kma + a)) * (b / (kmb + b));

        return rate;
    }
}

class Bi_ternary_complex extends Reaction {
    constructor(kcat, kma, kmb, kia, a, b, enz, substrates = [], products = []) {
        super();
        this.reaction_substrate_names = [a, b, enz];
        this.parameter_names = [kcat, kma, kmb, kia];
        this.substrates = substrates;
        this.products = products;
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const a = substrates[0];
        const b = substrates[1];
        const enz = substrates[2];

        // Parameters
        const kcat = parameters[0];
        const kma = parameters[1];
        const kmb = parameters[2];
        const kia = parameters[3];

        const num = kcat * enz * a * b;
        const den = ((kia * kmb) + (kmb * a) + (kma * b) + (a * b));

        const rate = num / den;

        return rate;
    }
}

class Bi_ping_pong extends Reaction {
    constructor(kcat, kma, kmb, a, b, enz, substrates = [], products = []) {
        super();
        this.reaction_substrate_names = [a, b, enz];
        this.parameter_names = [kcat, kma, kmb];
        this.substrates = substrates;
        this.products = products;
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const a = substrates[0];
        const b = substrates[1];
        const enz = substrates[2];

        // Parameters
        const kcat = parameters[0];
        const kma = parameters[1];
        const kmb = parameters[2];

        const rate = (kcat * enz * a * b) / ((kmb * a) + (kma * b) + (a * b));

        return rate;
    }
}

class Ter_seq_redam extends Reaction {
    constructor(kcat, kma, kmb, kmc, kia, kib, enz, a, b, c, substrates = [], products = []) {
        super();
        this.reaction_substrate_names = [a, b, c, enz];
        this.parameter_names = [kcat, kma, kmb, kmc, kia, kib];
        this.substrates = substrates;
        this.products = products;
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const a = substrates[0];
        const b = substrates[1];
        const c = substrates[2];
        const enz = substrates[3];

        // Parameters
        const kcat = parameters[0];
        const kma = parameters[1];
        const kmb = parameters[2];
        const kmc = parameters[3];
        const kia = parameters[4];
        const kib = parameters[5];

        const numerator = kcat * enz * a * b * c;
        const denominator = (kia * kib * kmc) + (kib * kmc * a) + (kia * kmb * c) + (kmc * a * b) + (kmb * a * c) + (kma * b * c) + (a * b * c);
        const rate = numerator / denominator;

        return rate;
    }
}

class Ter_seq_car extends Reaction {
    constructor(kcat, kma, kmb, kmc, kia, enz, a, b, c, substrates = [], products = []) {
        super();
        this.reaction_substrate_names = [a, b, c, enz];
        this.parameter_names = [kcat, kma, kmb, kmc, kia];
        this.substrates = substrates;
        this.products = products;
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const a = substrates[0];
        const b = substrates[1];
        const c = substrates[2];
        const enz = substrates[3];

        // Parameters
        const kcat = parameters[0];
        const kma = parameters[1];
        const kmb = parameters[2];
        const kmc = parameters[3];
        const kia = parameters[4];

        const rate = (kcat * enz * a * b * c) / ((kia * c) + (kmc * a * b) + (kmb * a * c) + (kma * b * c) + (a * b * c));

        return rate;
    }
}

class Bi_ternary_complex_small_kma extends Reaction {
    constructor(kcat, kmb, kia, a, b, enz, substrates = [], products = []) {
        super();
        this.reaction_substrate_names = [a, b, enz];
        this.parameter_names = [kcat, kmb, kia];
        this.substrates = substrates;
        this.products = products;
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const a = substrates[0];
        const b = substrates[1];
        const enz = substrates[2];

        // Parameters
        const kcat = parameters[0];
        const kmb = parameters[1];
        const kia = parameters[2];

        const rate = (kcat * enz * a * b) / ((kia * kmb) + (kmb * a) + (a * b));

        return rate;
    }
}


module.exports = {
    Uni,
    Bi,
    Bi_ternary_complex,
    Bi_ping_pong,
    Ter_seq_redam,
    Ter_seq_car,
    Bi_ternary_complex_small_kma
};
