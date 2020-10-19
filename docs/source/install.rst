.. _install:

===================
Installation
===================

----------
Platforms:
----------
*   Linux
*   Mac OS X
*   Windows

-------------
Requirenments
-------------
* a C++11-standard-compliant compiler
* Python 3.5+ (optional)

----------
Libraries:
----------
*   Xtl vers 0.6.13 (https://github.com/xtensor-stack/xtl)
*   Xtensor vers 0.21.4 (https://github.com/xtensor-stack/xtensor)
*   Eigen vers 3.3.7 (https://gitlab.com/libeigen/eigen)
*   Google-test vers 1.10.x (https://github.com/google/googletest)
*   Python Openmc package 0.12.0+ (optional to proceed ENDF formatted nuclear data libraries: https://github.com/openmc-dev/openmc)

---------------
Getting started
---------------

Build OpenBPS using CMake:

.. code-block:: sh

    mkdir build && cd build
    cmake ..
    # with Optimization:
    cmake -Doptimize=ON ..

