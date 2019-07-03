=========
Reactions
=========

kinetics works by specifying a ``Model()`` object to which reactions are added.

Reactions are first defined, selected from one of the reaction objects described here.  Parameters are set and the reaction added to the model.

For example:

.. code:: python

    import kinetics

    enzyme_1 = kinetics.Uni(kcat='kcat1', kma='kma1', a='a', enz='enz1',
                            substrates=['a'], productions=['b'])

    # To specify single parameter values
    enzyme_1.parameters = {'kcat1': 100,
                           'kma1': 500}

    # To specify parameter distributions
    enzyme_1.parameter_distributions = {'kcat1': norm(100,10),
                                        'kma1': uniform(25,50)}

    model = kinetics.Model()
    model.append(enzyme1)


Michaelis-Menten kinetics, irreversible.
----------------------------------------

Uni
~~~
The classic Miachelis-Menton equation for a single substrate.

.. autoclass:: kinetics.Uni

Bi
~~
Not strictly a true Miachaelis-Menton equation.
Use with caution.  Will give a reasonable prediction if one substrate is saturating, otherwise is likely wrong.

.. autoclass:: kinetics.Bi

Bi Ternary Complex
~~~~~~~~~~~~~~~~~~
For reactions with two substrates which have an sequential mechanism (either ordered or random).

.. autoclass:: kinetics.Bi_ternary_complex

Bi Ping Pong
~~~~~~~~~~~~
For reactions with two substrates which have a ping-pong mechanism

.. autoclass:: kinetics.Bi_ping_pong

Ter seq redam
~~~~~~~~~~~~~
A three substrate rate equation which can be used for Reductive Aminase enzymes.

.. autoclass:: kinetics.Ter_seq_redam

Ter seq car
~~~~~~~~~~~~~
A three substrate rate equation which can be used for Carboxylic Acid Reductase enzymes.

.. autoclass:: kinetics.Ter_seq_car

Bi ternary complex small kma
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A special case of Bi Ternary Complex where kma << kia.

.. autoclass:: kinetics.Bi_ternary_complex_small_kma



Michaelis-Menten kinetics, reversible.
--------------------------------------






