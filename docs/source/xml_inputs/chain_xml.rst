.. _chain_xml:

=========
Chain.xml
=========

Initially this files copying from OpenMC repo with some modifications to 
uncetainty analysis calculation. A depletion chain file has a 
``<depletion_chain>`` root element with one or more ``<nuclide>`` child elements.
The decay, reaction, and fission product data for each nuclide appears as child
elements of ``<nuclide>``.

---------------------
``<nuclide>`` Element
---------------------

The ``<nuclide>`` element contains information on the decay modes, reactions,
and fission product yields for a given nuclide in the depletion chain. This
element may have the following attributes:

  :name:
    Name of the nuclide

  :half_life:
    Half-life of the nuclide in [s]

  :d_half_life:
    Uncertainty of half-life value in [s]

  :decay_modes:
    Number of decay modes present

  :decay_energy:
    Decay energy released in [eV]

  :reactions:
    Number of reactions present

For each decay mode, a :ref:`chain_decay` appears as a child of
``<nuclide>``. For each reaction present, a :ref:`chain_reaction` appears as
a child of ``<nuclide>``. If the nuclide is fissionable, a :ref:`chain_nfy`
appears as well.

.. _chain_decay:

-------------------
``<decay>`` Element
-------------------

The ``<decay>`` element represents a single decay mode and has the following
attributes:

  :type:
    The type of the decay, e.g. 'ec/beta+'

  :target:
    The daughter nuclide produced from the decay

  :branching_ratio:
    The branching ratio for this decay mode

  :d_branching_ratio:
    The uncetainity of branching ratio values for this decay mode

.. _chain_reaction:

----------------------
``<reaction>`` Element
----------------------

The ``<reaction>`` element represents a single transmutation reaction. This
element has the following attributes:

  :type:
    The type of the reaction, e.g., '(n,gamma)'

  :Q:
    The Q value of the reaction in [eV]

  :target:
    The nuclide produced in the reaction (absent if the type is 'fission')

  :branching_ratio:
    The branching ratio for the reaction

.. _chain_nfy:

------------------------------------
``<neutron_fission_yields>`` Element
------------------------------------

The ``<neutron_fission_yields>`` element provides yields of fission products for
fissionable nuclides. Normally, it has the follow sub-elements:

  :energies:
    Energies in [eV] at which yields for products are tabulated

  :fission_yields:

    Fission product yields for a single energy point. This element itself has a
    number of attributes/sub-elements:

      :energy:
        Energy in [eV] at which yields are tabulated

      :products:
        Names of fission products

      :data:
        Independent yields for each fission product

