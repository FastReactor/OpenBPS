.. _reactions_xml:

=============
Reactions.xml
=============

The characterisctic of interactions between material nuclei and induced 
radiation are presented in ``<composition>`` form. The root element of the file
is ``<compositions>``.The presented cross-section data could be extentd  by 
including ``<impxslibs>`` data.

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

The ``<xslibs>`` element contains data specification and ``<xslib>`` elements

  :type:
    Define type of reaction description: by reaction-rates "rxs" or
    by cross-section presentation "cs"

    *Default*: "rxs"

  :ng:
    number of energy multigroup

-------------------------
``<composition>`` Element
-------------------------

The ``<compostion>`` store cross-section/reaction-rates data along with 
radiation flux arrays in multigroup form.

  :name:
    Name of material with cross-section data for correct work should match with
    one of the ``<materials>``

------------------
``<flux>`` Element
------------------

The particle multigroup flow distribution in the composition. Contains a real
array separated by blank space with size of ``ng`` attribute:

  :ng:
    number of energy multigroup

----------------------
``<spectrum>`` Element
----------------------

The same as ``<flux>`` but for normed distribution:

  :ng:
    number of energy multigroup

-------------------
``<dflux>`` Element
-------------------

The uncertainties values of the ``<flux>`` array:

  :ng:
    number of energy multigroup

-----------------------
``<dspectrum>`` Element
-----------------------

The uncertainties values of the  ``<spectrum>`` array:

  :ng:
    number of energy multigroup




