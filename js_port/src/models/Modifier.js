class Modifier {
    constructor() {
        this.substrateNames = [];
        this.substrateIndexes = [];

        this.parameterNames = [];
        this.parameterIndexes = [];
    }

    getSubstrateIndexes(substrate_names) {
        this.substrateIndexes = [];
        for (const name of this.substrateNames) {
            this.substrateIndexes.push(substrate_names.indexOf(name));
        }
    }

    getParameterIndexes(parameter_names) {
        this.parameterIndexes = [];
        for (const name of this.parameterNames) {
            this.parameterIndexes.push(parameter_names.indexOf(name));
        }
    }

    calc_modifier(substrates, parameters) {
        return [substrates, parameters];
    }
}

module.exports = Modifier;
