.. OpenBPS documentation master file, created by
   sphinx-quickstart on Wed Oct  7 23:25:16 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=================================================
The OpenBPS Open source code of BurnuP Simulation
=================================================

The OpenBPS is open source code is developing to solve wide range problems with
nuclide concetration changing over time of modelling process including as 
simulation with induced external partical source like depletion, transmutation,
activation (e.g. by neutrons in a nuclear reactor) well as simple nucleus decay. 
The code implements three main approaches to solve Baetman's equation:

1.  **Iteration method with uncetainity analysis option** [1].
2.  **Chebyshev reational approximation method (CRAM)** [2].
3.  **Direct analytical Baetman method** [3] (this method currently under 
    development).


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install.rst
   theory/index
   userguide/index
   examples/index
   xml_inputs/index
   license.rst
   reference.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
