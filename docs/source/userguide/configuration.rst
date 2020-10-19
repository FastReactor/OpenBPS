.. _configuration:

-------------
Configuration
-------------

User defined settings are stored in configuration.xml file. The file is read
by its own name in executable directory. Input directory could be switched by
flag -i (--inputs) as follows:

.. code-block:: sh

    openbps -i /path/to/configure.xml

All OpenBPS option are listed below:

-i, --inputs           Placement of configure.xml
-o, --output           Output director to store results file
-v, --verbose          Printing out execution log to the console

The minimum required information to start simple decay problem (i.m. without 
external source) for 1 second with 1 step containing in :ref:`configure_xml` 
presented below:

.. code-block:: xml
 
    <configure>
      <chain>/pathto/chain_endfb71.xml</chain>
      <inpmaterials>/pathto/materials.xml</inpmaterials>
      <outpmaterials>/pathto/materials.xml</outpmaterials>
      <numbers>1</numbers>
      <timestep>1.0</timestep>
      <decaykey>true</decaykey>
    </configure>

.. _keys:

The text of each record in example matches the value of according parametres.

----------------
Calculation keys
----------------

``decaykey`` is a calculation key which means the calculation without 
reading a reactions file. The others most imortant keys are:

    * ``decay_print`` printing output information in outlog.csv  with
      constitution of one defined by ``<filters>`` section(by default False)
    * ``uncertanties`` switched on a uncertainty calcaulation in main execution
      process. Allowed only with ``iteration`` ``<mathod>``
    *  ``outwrite`` defined whether or not printing information in 
       ``<outpmaterials>`` xml file (by default True)

.. _filters:

-------
Filters
-------

Output data written in outlog.csv could be manipulated by :ref:`filter_xml` 
block in :ref:`configure_xml`. By default in file printing 2 values
in 4 columns by time placed in row : Act - activity in [Bk] and decay heat in
[eV] with their according uncertainties along with time. One filters differs
from another by ``type`` attribute.
``header`` modificate default 4 values (by except time);
``nuclide`` including new columns with nuclear concentration for time for 
presented nuclides;
``time`` define time intevals for which data be printing;
``materials`` restrict a materials for which data be printing;
``exnuclide`` for nuclides presented in that list their fraction be excluded
from values of activity and decay heat.

