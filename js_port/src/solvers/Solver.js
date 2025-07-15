const ndarray = require('ndarray');
const rk4 = require('ode-rk4');

class Solver {
    constructor() {}

    run(reactions, species_names, species_values, parameter_values, time) {
        const y0 = ndarray(new Float64Array(species_values));

        const deriv = (dydt, y, t) => {
            for (let i = 0; i < dydt.shape[0]; i++) {
                dydt.set(i, 0);
            }
            for (const reaction of reactions) {
                const y_prime = reaction.reaction(y, species_names, parameter_values);
                for (let i = 0; i < y_prime.length; i++) {
                    dydt.set(i, dydt.get(i) + y_prime[i]);
                }
            }
        };

        const integrator = rk4(y0, deriv, 0, time[1] - time[0]);

        const y = [y0.data.slice()];
        for (let i = 1; i < time.length; i++) {
            integrator.step();
            y.push(integrator.y.data.slice());
        }

        return y;
    }
}

module.exports = Solver;
