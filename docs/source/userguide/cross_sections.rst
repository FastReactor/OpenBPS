.. _cross_sections:

--------------
Cross-sections
--------------

Total availiable amount of nuclides with cross-section data are restricted
especially for neutron induced depletion with a large amout of fission products.
But for calcualtion with a number of materials when cross-section data 
are repeated is more comfortable to store arrays in independent libriaries.
Description of externals paths are provided in :ref:`configuration` in the 
section ``<impxslibs>``. Number of external :ref:`libraries <xslib_xml>` are not
restricted:

.. code-block:: xml
 
    <impxslibs>
    <impxslib>/pathto/implib.xml</impxslib>
    </impxslibs>

The structure of ``implib.xml`` is identicall to :ref:`xslibs <xslib_xml>` 
section of composition file. The minimun information is energy discrization
and ``<xslibs>`` section with attribute ``<type>`` means cross-section data
if equal ``cs`` or reaction rates with  ``rx`` attribute value.

.. code-block:: xml

    <xslibs ng="4" type="cs">
    <energy ng="4">0.0252 0.215 0.465 1.0 2.15</energy>
    <xslib name="Cm245" ng="4" reaction="fission">376.3989 237.6379 217.4472 273.0525<xslib>
    </xslibs>


