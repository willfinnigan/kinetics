==========
Change Log
==========

2.0.0 (Current)
- **BREAKING CHANGE**: Refactored package structure 
- Moved `kinetics.rate` → `kinetics.reactions` (all reaction classes)
- Moved `kinetics.analysis` and `kinetics.optimisation` → `kinetics.code_to_be_refactored` (deprecated modules)
- Improved code organization with cleaner separation between models, reactions, solvers, and sampling
- Updated imports - all reaction classes still accessible as `kinetics.ReactionName`
- **NEW**: Modern packaging with `pyproject.toml` configuration (PEP 517/518)
- **NEW**: Optional dependency groups for performance, analysis, and optimization
- **NEW**: Separated core dependencies from optional ones
- **REMOVED**: `setup.py` - now using pyproject.toml only
- Added comprehensive API documentation with autoclass directives
- Enhanced documentation with :members: for better class documentation
- Updated all tutorials and examples to reflect new structure

1.3.7
- Reorganised code, Uncertainty and Senstivity modules are now imported directly into kinetics.
- Changed docs to reflect changine in Uncertainty module import

1.3.6
- Refactored Uncertainty.make_samples_from_distributions(..) to Uncertainty.sample_distributions(..)
- Added Uncertainty.sample_uniforms
- Added test_simple_model
- Added MixedInhibition2, which takes kic and kiu
- Refactored Senstivity module into Uncertainty
- Added check to model, if no parameter is set in either reactions or model, take mean of model.parameter_distribution
