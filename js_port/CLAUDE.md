# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a JavaScript port of the Python kinetics library for modeling enzyme reactions using ordinary differential equations. The package enables modeling of enzymatic reactions with single parameter values and provides the foundation for future probability distribution support.

## Common Commands

### Testing
```bash
npm test                    # Run all tests with Jest
npm test -- --watch        # Run tests in watch mode
```

### Development
```bash
npm install                 # Install dependencies
node examples/simple_example.js  # Run example simulation
```

### Building
```bash
npx webpack                 # Build bundle with webpack
```

## Code Architecture

### Core System Design

The library follows a modular architecture with clear separation of concerns:

1. **Model System** (`src/models/`)
   - `Model` class: Main orchestrator managing reactions, time series, and execution
   - `Reaction` base class: Abstract base for all reaction types with parameter/modifier support
   - `Modifier` class: Base class for reaction modifiers (inhibitors, activators, etc.)

2. **Reaction Types** (`src/reactions/`)
   - Various kinetic models implementing specific rate equations
   - Examples: `IrreversibleMichaelisMenten` (Uni, Bi, Ter variants), `ReversibleMichaelisMenten`, `MassTransfer`
   - Each reaction extends `Reaction` and implements `calculateRate()` method

3. **Solver System** (`src/solvers/`)
   - `Solver` class: ODE integration using Runge-Kutta 4th order (ode-rk4)
   - Takes reactions, species, parameters, and time points to produce concentration trajectories

4. **Sampling Framework** (`src/sampling/`)
   - `Sampler` base class: Interface for parameter sampling (placeholder for future Monte Carlo features)
   - `ScipyDistSampler`: Stub for future probability distribution sampling

### Key Data Flow

1. Create `Model` instance and configure time series with `set_time()`
2. Create `Reaction` objects with parameter names and set `.parameters` dict
3. Add reactions to model with `add_reaction()`
4. Call `run_single()` with starting concentrations and optional solver
5. Model extracts all species/parameters from reactions and sets up indexing
6. Solver integrates ODEs using reaction rate equations over time points
7. Returns object with trajectory data (`y`) and species names

### Architecture Patterns

- **Plugin-based Reactions**: Each reaction type implements `calculateRate()` for different kinetic models
- **Index-based Performance**: Reactions pre-compute parameter/species indices for fast array access during integration
- **Modifier System**: Reactions can have modifiers (inhibitors, activators) that alter rates
- **Dual Execution Modes**: `run_single()` for deterministic runs, `run_multi()` for future Monte Carlo support

### Dependencies

- **Core**: `ndarray` (n-dimensional arrays), `ode-rk4` (ODE solver)
- **Utility**: `jstat` (statistical functions)
- **Testing**: `jest` (test framework)
- **Build**: `webpack` (bundling)

### File Structure Notes

- `kinetics.js`: Main entry point exposing all classes
- `examples/`: Usage examples showing typical workflows
- `tests/`: Jest test files following `*.test.js` pattern
- Single-file reaction implementations for easy extension