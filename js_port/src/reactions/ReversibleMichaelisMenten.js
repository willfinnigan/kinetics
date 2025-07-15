const Reaction = require('../models/Reaction');

class UniUni_rev extends Reaction {
    constructor(kcatf, kcatr, kma, kmp, a, p, enz, substrates = [], products = []) {
        super();
        this.reaction_substrate_names = [a, p, enz];
        this.parameter_names = [kcatf, kcatr, kma, kmp];
        this.substrates = substrates;
        this.products = products;
    }

    calculateRate(substrates, parameters) {
        // Substrates
        const a = substrates[0];
        const p = substrates[1];
        const enz = substrates[2];

        // Parameters
        const kcatf = parameters[0];
        const kcatr = parameters[1];
        const kma = parameters[2];
        const kmp = parameters[3];

        const rate = (((kcatf / kma) * enz * a) - ((kcatr / kmp) * enz * p)) / (1 + (a / kma) + (p / kmp));

        return rate;
    }
}

class BiBi_Ordered_rev extends Reaction {
    constructor(kcatf, kcatr, kmb, kia, kib, kmp, kip, kiq, enz, a, b, p, q, substrates = [], products = []) {
        super();
        this.reaction_substrate_names = [a, b, p, q, enz];
        this.parameter_names = [kcatf, kcatr, kmb, kia, kib, kmp, kip, kiq];
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
        const kcatf = parameters[0];
        const kcatr = parameters[1];
        const kmb = parameters[2];
        const kia = parameters[3];
        const kib = parameters[4];
        const kmp = parameters[5];
        const kip = parameters[6];
        const kiq = parameters[7];

        // Rate equation
        const numerator = ((enz * kcatf * a * b) / (kia * kmb)) - ((enz * kcatr * p * q) / (kmp * kiq));
        const denominator = 1 + (a / kia) + (b / kib) + (q / kiq) + (p / kip) + ((a * b) / (kia * kmb)) + ((p * q) / (kmp * kiq));

        return (numerator / denominator);
    }
}

class BiBi_Random_rev extends Reaction {
    constructor(kcatf, kcatr, kmb, kia, kib, kmp, kip, kiq, a, b, p, q, enz, substrates = [], products = []) {
        super();
        this.reaction_substrate_names = [a, b, p, q, enz];
        this.parameter_names = [kcatf, kcatr, kmb, kia, kib, kmp, kip, kiq];
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
        const kcatf = parameters[0];
        const kcatr = parameters[1];
        const kmb = parameters[2];
        const kia = parameters[3];
        const kib = parameters[4];
        const kmp = parameters[5];
        const kip = parameters[6];
        const kiq = parameters[7];

        const num = ((kcatf * enz * a * b) / (kia * kmb)) - ((kcatr * enz * p * q) / (kmp * kiq));
        const dom = 1 + (a / kia) + (b / kib) + (p / kip) + (q / kiq) + ((a * b) / (kia * kmb)) + ((p * q) / (kmp * kiq));
        const rate = num / dom;

        return rate;
    }
}

class BiBi_Pingpong_rev extends Reaction {
    constructor(kcatf, kma, kmb, kia, kcatr, kmp, kmq, kip, kiq, enz, a, b, p, q, substrates = [], products = []) {
        super();
        this.reaction_substrate_names = [a, b, p, q, enz];
        this.parameter_names = [kcatf, kcatr, kma, kmb, kia, kmp, kmq, kip, kiq];
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
        const kcatf = parameters[0];
        const kcatr = parameters[1];
        const kma = parameters[2];
        const kmb = parameters[3];
        const kia = parameters[4];
        const kmp = parameters[5];
        const kmq = parameters[6];
        const kip = parameters[7];
        const kiq = parameters[8];

        const num = ((kcatf * enz * a * b) / (kia * kmb)) - ((kcatr * enz * p * q) / kip * kmq);
        const den = (a / kia) + ((kma * b) / (kia * kmb)) + (p / kip) + ((kmp * q) / (kip * kmq)) + ((a * b) / (kia * kmb)) + ((a * p) / (kia * kip)) + ((kma * b * q) / (kia * kmb * kiq)) + ((p * q) / (kip * kmq));

        return num / den;
    }
}

module.exports = {
    UniUni_rev,
    BiBi_Ordered_rev,
    BiBi_Random_rev,
    BiBi_Pingpong_rev
};
