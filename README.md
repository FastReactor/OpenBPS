# OpenBPS
Open source code of BurnuP Simulation
## Overview
The OpenBPS is open source code is developing to solve wide range problems with nuclide concetration changed over time of modelling process including as simulation
with induced external partical source like depletion, transmutation, activation (e.g. by neutrons in a nuclear reactor) well as simple nucleus decay. The code implments three main approaches to solve Baetmans equation:
1.  **Iteration method with uncetainity analysis option** [1].
2.  **Chebyshev reational approximation method (CRAM)** [2].
3.  **Direct analytical Baetman method** [3] (this method currently under development).

The OpenPBS relies on state of the art open source, open access instruments and data sources both to work and prepare input data and to solve linear equation system.
The information about every nucleus could be obtained using python modules of the OpenMC[4] is a community-developed Monte Carlo neutron and photon transport simulation code and
the next ENDF formated nuclear data libraries:
*   ENDFB-VII.1, ENDFB-VIII.0 https://www.nndc.bnl.gov/exfor
*   JEFF-3.3 https://www.oecd-nea.org/dbdata/jeff
*   JENDL-4.0 https://wwwndc.jaea.go.jp/ftpnd/jendl/jendl40he.html
*   Talys based TENDL-2017[5] https://tendl.web.psi.ch/tendl_2017/tendl2017.html

To solve nuclear reactor depletion problem is necessary to have multigroup cross-section library file which could be obtained by:
*   Utilizing open access code PREPRO19 https://www-nds.iaea.org/public/endf/prepro and ENDF formated induced neutron data presented above with python  automatization scripts
*   Another neutron cross-section multigroup librairies (for example https://www.polymtl.ca/merlin/version5.htm)

## Platforms

Platforms:

*   Linux
*   Mac OS X
*   Windows
## Requirenments
* a C++11-standard-compliant compiler
* Python 3.5+ (optional)
## External libraries

libraries:

*   Xtl vers 0.6.13 (https://github.com/xtensor-stack/xtl)
*   Xtensor vers 0.21.4 (https://github.com/xtensor-stack/xtensor)
*   Eigen vers 3.3.7 (https://gitlab.com/libeigen/eigen)
*   Google-test vers 1.10.x (https://github.com/google/googletest)
*   Python Openmc package 0.12.0+ (optional to proceed ENDF formatted nuclear data libraries: https://github.com/openmc-dev/openmc)

## Getting started
### Installation
Build OpenBPS using CMake:
``` bash
    mkdir build && cd build
    cmake ..
    with Optimization:
    cmake -Doptimize=ON ..
```
### Config
Config default
``` xml
<?xml version="1.0"?>
<configure>

  <chain>./examples/chain_endfb71.xml</chain>
  <nuclides>./examples/nuclides.xml</nuclides>
  <reaction>./examples/reactions.xml</reaction>
  <inpmaterials>./examples/materials_test.xml</inpmaterials>
  <numbers>1</numbers>
  <timestep>720000</timestep>
  <method>chebyshev</method>
  <decaykey>true</decaykey>
  <epb>1.e-5</epb>
  <cram_order>16</cram_order>

</configure>
```
