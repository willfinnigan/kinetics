const ndarray = require('ndarray');
const rk4 = require('ode-rk4');

class Solver {
    constructor() {}

    run(reactions, species_names, species_values, parameter_values, time) {
        // Use regular JavaScript arrays instead of ndarrays for ode-rk4
        const y0 = species_values.slice();

        const deriv = (dydt, y, t) => {
            // Initialize dydt to zeros
            for (let i = 0; i < y.length; i++) {
                dydt[i] = 0;
            }
            
            // Convert y to ndarray for reaction calculations
            const y_ndarray = ndarray(new Float64Array(y));
            
            // Sum up reaction contributions
            for (const reaction of reactions) {
                const y_prime = reaction.reaction(y_ndarray, species_names, parameter_values);
                for (let i = 0; i < Math.min(y_prime.length, y.length); i++) {
                    dydt[i] += y_prime[i];
                }
            }
        };

        const integrator = rk4(y0, deriv, 0, time[1] - time[0]);

        const y = [y0.slice()];
        for (let i = 1; i < time.length; i++) {
            integrator.step();
            y.push(integrator.y.slice());
        }

        return y;
    }
}

module.exports = Solver;
