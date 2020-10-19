.. _xslib_xml:

==========
Xslibs.xml
==========

Multigroup cross-section data could be presented both as an element
of reactions.xml file and independent from ``<impxslib>`` which filename and 
path are described in configure.xml. The only one common condition assumed is
presence of ``<energies>`` sections.   

----------------------
``<energies>`` Element
----------------------

The energies stored in dictionary like manner where attribute 

   :ng:
     the number of energy group

is a key and values is a real type array in text section separated by blank
space with size ``ng + 1``.

--------------------
``<xslibs>`` Element
--------------------

The ``<xslibs>`` element contains data specification and `<xslib>` elements

  :type:
    Define type of reaction description: by reaction-rates "rxs" or
    by cross-section presentation "cs"

    *Default*: "rxs"

  :ng:
    number of energy multigroups

-------------------
``<xslib>`` Element
-------------------

The ``<xslib>`` can store energy-dependent either cross-sections or reactions 
rates in form of the ``ng`` values separated by blank space.

  :name:
    Name of nuclide with data
  :type:
    Type of interaction in form "(p, r)" where p is the radiation particle, r -
    the result of interacation. The name to take account into depletion equation
    should match with reaction name in ``<reaction>`` ``chain`` element. 
  :ng:
    number of energy multigroup

