.. _material_xml:

=============
Materials.xml
=============

Materials.xml is a number of ``<material>`` records under the root section 
``<materials>``.

----------------------
``<material>`` Element
----------------------

The ``<material>`` element contains information of volume, mass and
external radiation power 

  :name:
    Name of the materials

  :volume:
    volume of material in [cm^3]

  :power:
    energy released in material by exposure of external radiation source [Wt]

  :mass:
    material mass, optional [g]

    *Default*: 0.0

For each nuclide in ``<nuclides>`` string array should be presented 
according value in an array with nuclear concentration ``<conc>``,
for uncertainty analysis ``<dconc>`` array can be placed too.

----------------------
``<nuclides>`` Element
----------------------

The array with names of nuclides

------------------
``<conc>`` Element
------------------

The array with nuclear concentrations in [atoms/cm^3]

-------------------
``<dconc>`` Element
-------------------

The array with uncertainty values of nuclear concentrations in [atoms/cm^3]
