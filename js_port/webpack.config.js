const path = require('path');

module.exports = {
    entry: './kinetics.js',
    output: {
        filename: 'kinetics.bundle.js',
        path: path.resolve(__dirname, 'dist'),
        library: 'kinetics',
        libraryTarget: 'umd',
    },
};
