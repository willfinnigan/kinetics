const Modifier = require('../models/Modifier');

class SubstrateInhibition extends Modifier {
    constructor(ki, a) {
        super();
        this.substrate_names = [a];
        this.parameter_names = [ki];
    }

    calc_modifier(substrates, parameters) {
        const ki = parameters[this.parameter_indexes[0]];
        const a = substrates[this.substrate_indexes[0]];

        substrates[this.substrate_indexes[0]] = a * (1 + a / ki);

        return [substrates, parameters];
    }
}

class CompetitiveInhibition extends Modifier {
    constructor(km, ki, i) {
        super();
        this.substrate_names = [i];
        this.parameter_names = [km, ki];
    }

    calc_modifier(substrates, parameters) {
        const km = parameters[this.parameter_indexes[0]];
        const ki = parameters[this.parameter_indexes[1]];
        const i = substrates[this.substrate_indexes[0]];

        parameters[this.parameter_indexes[0]] = km * (1 + i / ki);

        return [substrates, parameters];
    }
}

class MixedInhibition extends Modifier {
    constructor(kcat, km, ki, alpha, i) {
        super();
        this.substrate_names = [i];
        this.parameter_names = [kcat, km, ki, alpha];
    }

    calc_modifier(substrates, parameters) {
        const kcat = parameters[this.parameter_indexes[0]];
        const km = parameters[this.parameter_indexes[1]];
        const ki = parameters[this.parameter_indexes[2]];
        const alpha = parameters[this.parameter_indexes[3]];
        const i = substrates[this.substrate_indexes[0]];

        parameters[this.parameter_indexes[0]] = kcat / (1 + i / (alpha * ki));
        parameters[this.parameter_indexes[1]] = km * (1 + i / ki) / (1 + i / (alpha * ki));

        return [substrates, parameters];
    }
}

class MixedInhibition2 extends Modifier {
    constructor(kcat, km, kic, kiu, i) {
        super();
        this.substrate_names = [i];
        this.parameter_names = [kcat, km, kic, kiu];
    }

    calc_modifier(substrates, parameters) {
        const kcat = parameters[this.parameter_indexes[0]];
        const km = parameters[this.parameter_indexes[1]];
        const kic = parameters[this.parameter_indexes[2]];
        const kiu = parameters[this.parameter_indexes[3]];
        const i = substrates[this.substrate_indexes[0]];

        parameters[this.parameter_indexes[0]] = kcat / (1 + i / kiu);
        parameters[this.parameter_indexes[1]] = km * (1 + i / kic) / (1 + i / (kiu));

        return [substrates, parameters];
    }
}

class FirstOrder_Modifier extends Modifier {
    constructor(kcat, k, s) {
        super();
        this.substrate_names = [s];
        this.parameter_names = [kcat, k];
    }

    calc_modifier(substrates, parameters) {
        const kcat = parameters[this.parameter_indexes[0]];
        const k = parameters[this.parameter_indexes[1]];
        const s = substrates[this.substrate_indexes[0]];

        parameters[this.parameter_indexes[0]] = s * k * kcat;

        return [substrates, parameters];
    }
}

module.exports = {
    SubstrateInhibition,
    CompetitiveInhibition,
    MixedInhibition,
    MixedInhibition2,
    FirstOrder_Modifier
};
