class Reaction {
    constructor() {
        // These are set by the user
        this.parameters = {};
        this.parameterDistributions = {};

        // indexes used to access values during model run
        this.substrateIndexes = [];
        this.parameterIndexes = [];

        // These are set when the reaction is set up
        this.reactionSubstrateNames = [];
        this.parameterNames = [];
        this.substrates = [];
        this.products = [];

        // These are added as needed
        this.modifiers = [];
        this.checkPositive = false;
        this.checkLimitsFunctions = [];
    }

    setParameterDefaultsToMean() {
        for (const name in this.parameterDistributions) {
            if (!this.parameters.hasOwnProperty(name)) {
                if (Array.isArray(this.parameterDistributions[name])) {
                    this.parameters[name] = (this.parameterDistributions[name][0] + this.parameterDistributions[name][1]) / 2;
                } else {
                    // Assuming the distribution object has a mean() method
                    this.parameters[name] = this.parameterDistributions[name].mean();
                }
            }
        }
    }

    setupReaction(speciesNames, parameterNames) {
        // get indexes for rate calculation
        this.substrateIndexes = this.reactionSubstrateNames.map(name => speciesNames.indexOf(name));
        this.parameterIndexes = this.parameterNames.map(name => parameterNames.indexOf(name));

        // set up modifiers
        for (const modifier of this.modifiers) {
            modifier.getSubstrateIndexes(this.reactionSubstrateNames);
            modifier.getParameterIndexes(this.parameterNames);
        }
    }

    addModifier(modifier) {
        for (const name of modifier.parameterNames) {
            if (!this.parameterNames.includes(name)) {
                this.parameterNames.push(name);
            }
        }

        for (const name of modifier.substrateNames) {
            if (!this.reactionSubstrateNames.includes(name)) {
                this.reactionSubstrateNames.push(name);
            }
        }

        this.modifiers.push(modifier);
    }

    calculateRate(substrates, parameters) {
        return 0;
    }

    reaction(y, substrateNames, parameterValues) {
        // Get the substrates from y using the substrate indexes
        const substrates = this.substrateIndexes.map(index => y.get(index));

        // Get the parameters using the parameter indexes
        let parameters = this.parameterIndexes.map(index => parameterValues[index]);

        // calculate the effects of any modifiers
        for (const modifier of this.modifiers) {
            [substrates, parameters] = modifier.calc_modifier(substrates, parameters);
        }

        // calculate the rate (this function is modified by the user)
        const rate = this.calculateRate(substrates, parameters);

        // calculate the change in substrate concentrations (y_prime)
        const substrateIndices = this.substrates.map(name => substrateNames.indexOf(name));
        const productIndices = this.products.map(name => substrateNames.indexOf(name));

        const yPrime = new Array(y.shape[0]).fill(0);

        for (const index of substrateIndices) {
            yPrime[index] -= rate;
        }

        for (const index of productIndices) {
            yPrime[index] += rate;
        }


        let finalYPrime = this.modifyProduct(yPrime, substrateNames);

        if (this.checkPositive) {
            finalYPrime = finalYPrime.map(value => Math.max(0, value));
        }

        return finalYPrime;
    }

    modifyProduct(yPrime, substrateNames) {
        return yPrime;
    }

    samplingLimits(parameterDict) {
        for (const func of this.checkLimitsFunctions) {
            if (!func(parameterDict)) {
                return false;
            }
        }
        return true;
    }
}

module.exports = Reaction;
