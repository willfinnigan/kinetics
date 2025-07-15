const Model = require('../src/models/Model');
const {
    Uni
} = require('../src/reactions/IrreversibleMichaelisMenten');
const Solver = require('../src/solvers/Solver');

// Helper function to check if two arrays are close enough
function assert_allclose(actual, expected, atol = 1e-8, rtol = 1e-5) {
    expect(actual.length).toBe(expected.length);
    for (let i = 0; i < actual.length; i++) {
        expect(actual[i]).toBeCloseTo(expected[i], 5);
    }
}


test('simple one enzyme model', () => {
    const model = new Model();
    model.set_time(0, 1000, 100);

    const enzyme_1 = new Uni('enz1_kcat', 'enz1_km', 'A', 'enz_1', ['A'], ['B']);
    enzyme_1.parameters = {
        'enz1_kcat': 100,
        'enz1_km': 10000
    };
    model.add_reaction(enzyme_1);

    const solver = new Solver();

    const result = model.run_single({
        "A": 10000,
        "enz_1": 5
    }, solver);

    const df = {};
    for (let i = 0; i < result.species_names.length; i++) {
        df[result.species_names[i]] = result.y.map(row => row[i]);
    }


    const start = df['A'][0];
    const end = df['B'][99];

    const expected = [10000.0, 10000.0];
    const actual = [start, end];

    assert_allclose(actual, expected, 1, 1);
});


test('simple two enzyme model', () => {
    const model = new Model();
    model.set_time(0, 1000, 100);

    const enzyme_1 = new Uni('enz1_kcat', 'enz1_km', 'A', 'enz_1', ['A'], ['B']);
    enzyme_1.parameters = {
        'enz1_kcat': 100,
        'enz1_km': 10000
    };
    model.add_reaction(enzyme_1);

    const enzyme_2 = new Uni('enz2_kcat', 'enz2_km', 'B', 'enz_2', ['B'], ['C']);
    enzyme_2.parameters = {
        'enz2_kcat': 100,
        'enz2_km': 10000
    };
    model.add_reaction(enzyme_2);

    const solver = new Solver();

    const result = model.run_single({
        "A": 10000,
        "enz_1": 5,
        "enz_2": 5
    }, solver);

    const df = {};
    for (let i = 0; i < result.species_names.length; i++) {
        df[result.species_names[i]] = result.y.map(row => row[i]);
    }

    const start = df['A'][0];
    const end = df['C'][99];

    const expected = [10000.0, 10000.0];
    const actual = [start, end];

    assert_allclose(actual, expected, 1, 1);
});
