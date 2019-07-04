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

.. autoclass:: kinetics.Uni

The classic Miachelis-Menton equation for a single substrate.

.. math::
    rate = \frac{c_{enz}\cdot k_{cat}\cdot c_{A}}{c_{A}+K_{M}^{A}}

Bi
~~

.. autoclass:: kinetics.Bi

Not strictly a true Miachaelis-Menton equation.
Use with caution.  Will give a reasonable prediction if one substrate is saturating, otherwise is likely wrong.

.. math::
    rate = c_{enz}\cdot k_{cat}\cdot \frac{c_{A}}{c_{A}+K_{M}^{A}} \cdot \frac{c_{B}}{c_{B}+K_{M}^{B}}

Bi Ternary Complex
~~~~~~~~~~~~~~~~~~

.. autoclass:: kinetics.Bi_ternary_complex

For reactions with two substrates which have an sequential mechanism (either ordered or random).

.. math::
    rate = \frac{c_{enz}\cdot k_{cat}\cdot c_{A}\cdot c_{B}}{(K_{I}^{A}\cdot K_{M}^{B}) + (K_{M}^{B}\cdot c_{A}) + (K_{M}^{A}\cdot c_{B}) + (c_{A} \cdot c_{B})}

Bi Ping Pong
~~~~~~~~~~~~

.. autoclass:: kinetics.Bi_ping_pong

For reactions with two substrates which have a ping-pong mechanism

.. math::
    rate = \frac{c_{enz}\cdot k_{cat}\cdot c_{A}\cdot c_{B}}{(K_{M}^{B}\cdot c_{A}) + (K_{M}^{A}\cdot c_{B}) + (c_{A} \cdot c_{B})}


Ter seq redam
~~~~~~~~~~~~~

.. autoclass:: kinetics.Ter_seq_redam

A three substrate rate equation which can be used for Reductive Aminase enzymes.

.. math::
    rate = \frac{c_{enz}\cdot k_{cat}\cdot c_{A}\cdot c_{B}\cdot c_{C}}
                     {(K_{I}^{A}\cdot K_{I}^{B} \cdot K_{M}^{C}) +
                     (K_{I}^{B}\cdot K_{M}^{C} \cdot c_{A}) +
                     (K_{I}^{A}\cdot K_{M}^{B}\cdot c_{C}) +
                     (K_{M}^{C}\cdot c_{A}\cdot c_{B}) +
                     (K_{M}^{B}\cdot c_{A}\cdot c_{C}) +
                     (K_{M}^{A}\cdot c_{B}\cdot c_{C}) +
                     (c_{A} \cdot c_{B} \cdot c_{C})}


Ter seq car
~~~~~~~~~~~

.. autoclass:: kinetics.Ter_seq_car

A three substrate rate equation which can be used for Carboxylic Acid Reductase enzymes.

.. math::
    rate = \frac{c_{enz}\cdot k_{cat}\cdot c_{A}\cdot c_{B} \cdot c_{C}}
                {(K_{I}^{A}\cdot c_{C}) +
                (K_{M}^{C}\cdot c_{A} \cdot c_{B}) +
                (K_{M}^{B}\cdot c_{A} \cdot c_{C}) +
                (K_{M}^{A}\cdot c_{B} \cdot c_{C}) +
                (c_{A} \cdot c_{B} \cdot c_{C})}


Bi ternary complex small kma
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: kinetics.Bi_ternary_complex_small_kma

A special case of Bi Ternary Complex where kma << kia.

.. math::
    rate = \frac{c_{enz}\cdot k_{cat}\cdot c_{A}\cdot c_{B}}
                {(K_{I}^{A}\cdot K_{M}^{B}) +
                (K_{M}^{B}\cdot c_{A}) +
                (c_{A} \cdot c_{B})}



Michaelis-Menten kinetics, reversible.
--------------------------------------

UniUni Reversible
~~~~~~~~~~~~~~~~~

.. autoclass:: kinetics.UniUni_rev

.. math::
    rate = \frac{(c_{enz}\cdot k_{cat}^{fwd}\cdot c_{A}) - (c_{enz}\cdot k_{cat}^{rev}\cdot c_{P})}
                     {1 +
                     \frac{c_{A}}{K_{M}^{A}} +
                     \frac{c_{P}}{K_{M}^{P}} }

BiBi Ordered Rev
~~~~~~~~~~~~~~~~

.. autoclass:: kinetics.BiBi_Ordered_rev

BiBi Random Rev
~~~~~~~~~~~~~~~

.. autoclass:: kinetics.BiBi_Random_rev

BiBi Pingpong Rev
~~~~~~~~~~~~~~~~~

.. autoclass:: kinetics.BiBi_Pingpong_rev


Modifiers of Michaelis-Menten kinetics eg for Inhibition
--------------------------------------------------------
Modifications to rate equations for things like competitive inhibition can applied as follows:

(Remember to add new parameters to the reaction parameters)

Modifications are applied at each timestep of the model, for example calculating the apparent Km resulting from competitive inhibtion.

This feature allows the easy modification of the pre-defined rate equations.

.. code:: python

    enzyme_1.add_modifier(kinetics.CompetitiveInhibition(km='kma1', ki='ki1', i='I'))
    enzyme_1.parameters.update({'ki1': 25})


.. autoclass:: kinetics.SubstrateInhibition

.. autoclass:: kinetics.CompetitiveInhibition

.. autoclass:: kinetics.MixedInhibition

.. autoclass:: kinetics.FirstOrder_Modifier











