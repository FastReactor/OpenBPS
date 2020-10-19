.. _configure_xml:

=============
Configure.xml
=============

The main configuration settings are specified in configure.xml.

-------------------
``<chain>`` Element
-------------------

The path to chain .xml files

----------------------
``<nuclides>`` Element
----------------------

The path to nuclides .xml file.

--------------------------
``<inpmaterials>`` Element
--------------------------

The path to materials input .xml file.

---------------------------
``<outpmaterials>`` Element
---------------------------

The path to the materials output .xml file.

----------------------
``<reaction>`` Element
----------------------

The path to the reactions .xml file.

-----------------------
``<impxslibs>`` Element
-----------------------

The sections should contain a number of xslib .xml with cross-section
multigroup data.

--------------------------
``<impxslib>`` Sub-Element
--------------------------

Parent element is ``<impxslibs>``
The path to the impxslib .xml file.

--------------------
``<output>`` Element
--------------------

This section contains the output directory of OpenBPS

  *Default*: "."

----------------------
``<timestep>`` Element
----------------------

The total simulation time

------------------------
``<timerecord>`` Element
------------------------

The total simulation time described as follows formats with attributes:

* year - int number of years;
* month - int number of months
* day - int number of days;
* hour - int number of hours;
* minute - int number of minutes;
* second - int number of seconds;

If both ``<timerecord>`` and ``<timestep>`` are presented total simulation time
is calculated from ``<timerecord>`` information

--------------------
``<method>`` Element
--------------------

One of the OpenBPS calculation method : ``iteration``, ``chebyshev``,``baetman``

  *Default*: "iteration"

-----------------
``<epb>`` Element
-----------------

The calculation accuracy of the ``iteration`` method

  *Default*: 1.e-3

-------------------
``<order>`` Element
-------------------

The order of the ``chebyshev`` method. Possible values are
16; 48.  With increasing order value code performs calculations problem slower. 

  *Default*: 16

---------------------------
``<is_outrewrite>`` Element
---------------------------

The bool value signs whether write down the output materials xml file

  *Default*: True

-------------------------
``<decay_print>`` Element
-------------------------

The bool value indicates writing a calculation results expressed with 
activity and decay heat by default and managed by ``<filters>`` section

  *Default*: False

--------------------------
``<uncertanties>`` Element
--------------------------

The bool value if true a calculation performs with uncertainties analysis

    .. note :: Works only with ``iteration`` method

  *Default*: False

----------------------
``<decaykey>`` Element
----------------------

The bool value means a simple decay calculation without external sources and
not reading reactions.xml file

  *Default*: False

.. _filter_xml:

---------------------
``<filters>`` Element
---------------------

This section contains ``<filter>`` subelements manipulate over outlog.csv file

------------------------
``<filter>`` Sub-Element
------------------------

The parent element is ``<filters>``

The main purpose of this section is moderating content of outlog.csv file by
the next attributes:

   :type:
     Strings values is one of the list : 
     * material;
     * time;
     * nuclides;
     * header;
     * exnuclices.
   
For every type text content of the sections can differs from each other. For
``material`` type it would be name of materials for which log should be printed.
``time`` filter consist with array of real type data which means a time 
intervals to print out information. ``header`` can restrict the default headers
by pointed theese for which data be printed {``Act``,sec-1``;``Q, Mev``;
``dAct, sec-1``;``dQ, Mev``}. The ``nuclides`` contains the names of nuclides
for those nuclear concentration during simulation will be written to the file. 
The ``exnuclides`` excludes the fraction of presented  nuclides from values of 
decay rate and decay heat.

:Example:

The simple example to show how filters work.

.. code-block:: xml

    <filters>
      <filter type="materials">"test"</filter>
      <filter type="time">10.0 12.0 20.0 23.0</filter>
      <filter type="header">"Q, Mev", "dQ, Mev"</filter>
      <filter type="nuclide">Cs133</filter>
      <filter type="exnuclide">Pu238 Pu239 Pu240 Pu241 Pu242</filter>
    </filters>

With ouputs:

    dt;Q, Mev;dQ, Mev;Cs133
    test;
    11.0;8.521e+13;3.47642e+12;5.17185e-19;
    21.0;5.77999e+14;2.40786e+13;3.462e-18;
    22.0;6.37238e+14;2.66199e+13;3.81353e-18;

that means information is only for "test" ``material`` for time values which
lay in ``time`` intervals from 10.0 to 12.0 and 20.0 to 23.0 only for decay 
heat and its uncertainty values excluding from cummulative plutonium isotope 
additives with nuclear concentration of `Cs133`.
  


