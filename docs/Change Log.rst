==========
Change Log
==========

1.3.6
- Refactored Uncertainty.make_samples_from_distributions(..) to Uncertainty.sample_distributions(..)
- Added Uncertainty.sample_uniforms
- Added test_simple_model
- Added MixedInhibition2, which takes kic and kiu
- Refactored Senstivity module into Uncertainty
- Added check to model, if no parameter is set in either reactions or model, take mean of model.parameter_distribution

1.3.7
- Reorganised code, Uncertainty and Senstivity modules are now imported directly into kinetics.
- Changed to 'from setuptools import setup' in setup.py
- Changed docs to reflect changine in Uncertainty module import
