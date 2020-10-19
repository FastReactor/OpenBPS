.. _reactions:

------------
Compositions
------------

To calcualte nuclear concentrations are changed by  induced external particle
field is necesserary to know the value of :ref:`reactions rate <reactions_xml>`.
Reactions rate could be expressed as
external flux [particle/sec/cm^2] x cross-section[cm^2].
To describe that values  multigroup energy model are used in the OpenBPS. Every
energy-dependent distribution with dimension (``NG``) should be presented by 
energy discretization array with a size ``NG + 1``. Meantime it could be a 
number of values of flux/cross-sections/reaction-rates. For example:

.. code-block:: xml
 
    <energy ng="4">0.0252 0.215 0.465 1.0 2.15</energy>
    <energy ng="1">0.0252 2.15</energy>
    <xslibs ng="1" typex="cs">
    </xslibs>
    <flux ng="4">1.0 0.00000E+00 0.00000E+00 0.00000E+00</flux>

Data of ``<flux>`` matches to energy discretization with ``ng="4"`` meanwhile 
``<xslibs>`` to one with ``ng="1"``. To bind external field interactions data
with nuclear concentration, composition attribute ``name`` should be equal to
the same attribute in materials description.

    **Note:** For neutron induced radiation it is very important
    to have ``<flux>`` or normed ``<spectrum>`` array to calculate right
    fraction of fission product yield which are energy dependent too.
