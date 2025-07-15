===
API
===

Core Classes
============

Model
-----

.. autoclass:: kinetics.Model
   :members:

Reaction
--------

.. autoclass:: kinetics.Reaction
   :members:

Modifier
--------
.. autoclass:: kinetics.Modifier
   :members:

Result Classes
==============

SingleModelResult
------------------

.. autoclass:: kinetics.SingleModelResult
   :members:

MultiModelResult
----------------

.. autoclass:: kinetics.MultiModelResult
   :members:

Solvers
=======

ODESolver Interface
-------------------

.. autoclass:: kinetics.ODESolver
   :members:

SciPySolver
-----------

.. autoclass:: kinetics.SciPySolver
   :members:

JaxSolver
---------

.. autoclass:: kinetics.JaxSolver
   :members:

Sampling
========

Sampler Interface
-----------------

.. autoclass:: kinetics.Sampler
   :members:

ScipyDist_Sampler
-----------------

.. autoclass:: kinetics.ScipyDist_Sampler
   :members:

SalibLatinHypercubeSampler
--------------------------

.. autoclass:: kinetics.SalibLatinHypercubeSampler
   :members:

SalibSaltelliSampler
--------------------

.. autoclass:: kinetics.SalibSaltelliSampler
   :members:

Sensitivity Analysis
====================

SensitivityResult
-----------------

.. autoclass:: kinetics.sampling.analysis.sensitivity_analysis.SensitivityResult
   :members:

Analysis Functions
------------------

.. autofunction:: kinetics.analyze_sensitivity_at_timepoint

.. autofunction:: kinetics.analyze_sensitivity_time_to_concentration
