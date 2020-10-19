.. _overview:

--------
Overview
--------

The OpenBPS is developing to solve nuclide :ref:`transfromation problem
<decay_equation>` due to the decay process. The minimal information is
necessary to solve such problem is information about material with 
number of isotope and their nuclear concentrations and a file
with all required transition from one nuclides to others called a chain.
File format specification is .xml.The short description of inputs in the
table below:

.. table:: Table 1. Description of input OpenBPS files

    +----------------------+-----------------------------------------+
    | Filename             | Description                             |
    +======================+=========================================+
    |  :ref:`configure_xml`| Main configuration file with            |
    |                      | settings and information about placement|    
    |                      | of other files. configure.xml   should  |
    |                      | be placed itself by default in directory|
    |                      | where OpenBPS executive are called from |
    +----------------------+-----------------------------------------+
    | :ref:`material_xml`  | All about material nuclides composition |
    | (name of file can be | along with `power` [Wt] normalization   |
    | arbitrary)           | and `volume` [cm^3] attributes.         | 
    +----------------------+-----------------------------------------+ 
    | :ref:`chain_xml`     | Nuclides changing through decay process |
    | (name of file can be | information including data for          |
    | arbitrary)           | transformation by external field. For   |
    |                      | neutrons it could be capture, fission,  |
    |                      | e.t.c.                                  |
    +----------------------+-----------------------------------------+
    | :ref:`nuclide_xml`   | The individual nuclide data including   |
    | (name of file can be | mass [a.e.m], index (ZAM), and other    |
    | arbitrary)           | specific features.                      |
    | `optional`           |                                         |
    +----------------------+-----------------------------------------+
    | :ref:`reactions_xml` | For material under induced radiation is |
    | (name of file can be | provided values with reaction rates,    |
    | arbitrary)           | flux and/or cross-section structure.    |
    | `optional`           |                                         |
    +----------------------+-----------------------------------------+
    | Number of            | Additional cross-section data in multi- |
    | :ref:`xslib_xml`     | group energy discretization. The data   |
    | (names of files can  | could be used only with reactions.xml   |
    | be arbitrary)        | file. Number of cross-section files is  |
    | `optional`           | not restricted                          |
    +----------------------+-----------------------------------------+
