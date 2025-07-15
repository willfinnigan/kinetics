const {
    Model
} = require('../src/models/Model');
const {
    Uni
} = require('../src/reactions/IrreversibleMichaelisMenten');

// Define reactions
const enzyme_1 = new Uni('enz1_kcat', 'enz1_km', 'A', 'enz_1', ['A'], ['B']);
enzyme_1.parameters = {
    'enz1_kcat': 200,
    'enz1_km': 8000
};

const enzyme_2 = new Uni('enz2_kcat', 'enz2_km', 'B', 'enz_2', ['B'], ['C']);
enzyme_2.parameters = {
    'enz2_kcat': 30,
    'enz2_km': 2000
};

// Create model
const model = new Model();
model.add_reaction(enzyme_1);
model.add_reaction(enzyme_2);
model.set_time(0, 120, 100);

// Set starting concentrations and run
const starting_concs = {
    "A": 10000,
    "enz_1": 4,
    "enz_2": 10
};

const result = model.run_single(starting_concs);

// The result object contains the simulation data.
// The data can be plotted using a library like Plotly.js or Chart.js.
// For example, to plot the concentration of species 'A' over time:
console.log("Concentration of A:", result.y.map(y => y[result.species_names.indexOf('A')]));
console.log("Concentration of B:", result.y.map(y => y[result.species_names.indexOf('B')]));
console.log("Concentration of C:", result.y.map(y => y[result.species_names.indexOf('C')]));
