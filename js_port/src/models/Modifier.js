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

module.exports = Modifier;
