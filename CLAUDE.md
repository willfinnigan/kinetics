# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Python package for modeling enzyme reactions using ordinary differential equations. The package allows both single parameter values and probability distributions for parameters, enabling uncertainty analysis and sensitivity studies.

## Common Commands

### Testing
```bash
pytest tests/                    # Run all tests
pytest tests/test_single_param_model.py  # Run specific test file
```

When writing code for a new feature, if possible write a failing test first.  When writing tests, please refer to the documentation and to existing tests for how code should be written.  Follow TTD - you should then make this test pass.  Once the test is passing, consider whether any refactoring would help (and ensure tests still pass once this is complete).  

### Installation and Setup
```bash
pip install -e .                # Install package in development mode
pip install -r requirements.txt # Install dependencies (if requirements.txt exists)
```

### Building/Distribution
```bash
python -m build                    # Build distribution packages
```

### Documentation
Each time you modify any documentation, please check that it builds ok without errors - fix if necessary.
```bash
cd docs/                           # Navigate to documentation directory
make html                          # Build HTML documentation with Sphinx
make clean                         # Clean build artifacts
```

**Documentation Structure:**
- `docs/` - Sphinx documentation source files (RST format)
- `docs/_build/html/` - Built HTML documentation
- `docs/index.rst` - Main documentation index
- Available tutorials: Installation, Simple Tutorial, Advanced Tutorial, Advanced Tutorial 2, Sensitivity Analysis, Reactions, Custom Reactions
- API documentation auto-generated from docstrings

## Code Architecture

### Core Components

1. **Model System** (`kinetics/models/`)
   - `Model` class: Main orchestrator that manages reactions, time series, and execution
   - `Reaction` class: Base class for different reaction types with parameter handling
   - Result classes: `SingleModelResult` and `MultiModelResult` for storing outputs

2. **Reaction Types** (`kinetics/reactions/`)
   - Various kinetic models: Michaelis-Menten variants, mass action, thermodynamic models
   - Each reaction implements rate equations and parameter management
   - Support for both NumPy and JAX arrays for performance

3. **Solvers** (`kinetics/solvers/`)
   - `ODESolver` interface for pluggable ODE solvers
   - `SciPySolver`: Uses scipy.integrate for solving ODEs
   - `JAXSolver`: JAX-based solver for performance-critical applications

4. **Sampling** (`kinetics/sampling/`)
   - `Sampler` interface for parameter sampling
   - `ScipyDist_Sampler`: Uses scipy.stats distributions
   - `SALib_Sampler`: Integration with SALib for sensitivity analysis

### Key Design Patterns

- **Plugin Architecture**: Solvers and samplers implement interfaces for easy extension
- **Dual Array Support**: Code supports both NumPy and JAX arrays with runtime detection
- **Parameter Distribution**: Parameters can be single values or probability distributions
- **Reaction Composition**: Models are built by adding multiple reaction objects

### Data Flow

1. Create `Model` instance and set time series
2. Add `Reaction` objects with parameters (values or distributions)
3. Model extracts species and parameters from reactions
4. Solver integrates ODEs over time using reaction rate equations
5. Results stored in result classes with plotting capabilities

### Recent Refactoring

The codebase is currently being refactored (branch: v2):
- `kinetics/rate/` → `kinetics/reactions/` (reaction types moved)
- `kinetics/analysis/` and `kinetics/optimisation/` → `kinetics/code_to_be_refactored/`
- Focus on cleaner separation between models, reactions, solvers, and sampling

### Dependencies

Core: numpy, scipy, matplotlib, pandas, tqdm
Advanced: jax, diffrax (for performance), SALib (sensitivity analysis), deap (optimization)