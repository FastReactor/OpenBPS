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

tag **chain** : path to xml file with chain info (in examples  dir `chain_endfb71.xml`)

tag **nuclides** : path to xml file with nuclides (in examples  dir `nuclides.xml`)

tag **inpmaterials** : path to xml file with material nuclide conentration composition (in examples  dir `materials_test.xml`)

tag **reaction** : (optional) path to xml file with particle induced nuclide reactions(in examples  dir `reactions.xml`)

tag **numbers** : int number of simulation steps

tag **timestep** : double time of each step simulation

tag **method** : str {'iteration', 'chebyshev', 'baetman'} 

tag **decaykey** : bool calculation with/without external source

tag **epb** : double the accuracy of iteration calculation

tag **cram_order** : int {16, 48} order of CRAM

### Run
After successfull installation user can run OpenBPS from console from directory containing config.xml.
For example:
``` bash
    ./build/bin/openbps
```
**Output** could be presented in form xml file if tag **outmaterials** is pointed out in configure.xml or 
activity and decay heat with uncertanties could be write by **decay_print** _true_ key in file outlog.csv.

## Theory
   
The system of equations nuclide kinetics is inhomogeneous linear system of equations:
####  ``` 1 ``` 
<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{d&space;y_{i}}{d&space;t}=\sum_{j=1}^{n}&space;a_{i&space;j}&space;y_{j}&plus;q_{i},&space;\quad(i=1,2,&space;\dots,&space;n)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{d&space;y_{i}}{d&space;t}=\sum_{j=1}^{n}&space;a_{i&space;j}&space;y_{j}&plus;q_{i},&space;\quad(i=1,2,&space;\dots,&space;n)" title="\frac{d y_{i}}{d t}=\sum_{j=1}^{n} a_{i j} y_{j}+q_{i}, \quad(i=1,2, \dots, n)" /></a>

with the initial condition:
####  ``` 2 ``` 
<a href="https://www.codecogs.com/eqnedit.php?latex=y_{i}(0)=y_{i&space;0},&space;\quad(i=1,2,&space;\dots,&space;n)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?y_{i}(0)=y_{i&space;0},&space;\quad(i=1,2,&space;\dots,&space;n)" title="y_{i}(0)=y_{i 0}, \quad(i=1,2, \dots, n)" /></a>

where y<sub>i</sub> is the concentration of the i th nuclide, a<sub>ij</sub>  are the coefficients characterizing the channels of the transformations of the i th nuclide from the j–th nuclide, qi , an external source, n is the number of nuclides, y<sub>i0</sub>  is the concentration of the i th nuclide at time t<sub>0</sub> .

- ## Iteration Method 
  
The system of equations (1) – (2) in matrix form:
####  ``` 3 ``` 
<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{d&space;\vec{y}}{d&space;t}=\hat{A}&space;\vec{y}&plus;\vec{q},&space;\quad\left(A&space;\equiv\left[a_{i&space;j}\right]\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{d&space;\vec{y}}{d&space;t}=\hat{A}&space;\vec{y}&plus;\vec{q},&space;\quad\left(A&space;\equiv\left[a_{i&space;j}\right]\right)" title="\frac{d \vec{y}}{d t}=\hat{A} \vec{y}+\vec{q}, \quad\left(A \equiv\left[a_{i j}\right]\right)" /></a>

####  ``` 4 ``` 

<a href="https://www.codecogs.com/eqnedit.php?latex=\vec{y}(0)=\vec{y}_{0}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\vec{y}(0)=\vec{y}_{0}" title="\vec{y}(0)=\vec{y}_{0}" /></a>

The system of equations (1) – (2) can be solved analytically as the sum of its particular solution and the General solution of the corresponding homogeneous linear system of equations. In this case, there is a need matrix of size n × n multiplication. Therefore, this method of solution is often used for a partial transition matrices, for example, triangular.
Another widely known method of solving the problem is to decompose the solution in the range of the exponential function, which leads to the need to use a recurrence relation.
In the code the solution of equations (1) – (2) is the iterative method. For the k th iteration solution is obtained:

####  ``` 5 ``` 
<a href="https://www.codecogs.com/eqnedit.php?latex=y_{k}^{i}(\tau)=y_{k-1}^{i}(\tau)&plus;\sum_{j&space;\neq&space;i}&space;d&space;y_{k-1}^{j}&space;\frac{\lambda^{j&space;\rightarrow&space;i}}{\lambda_{p}^{j}}&space;\frac{1-\exp&space;\left(-\lambda_{p}^{i}&space;\tau\right)}{\lambda_{p}^{i}&space;\tau},&space;\quad(j=1,2,&space;\ldots,&space;k)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?y_{k}^{i}(\tau)=y_{k-1}^{i}(\tau)&plus;\sum_{j&space;\neq&space;i}&space;d&space;y_{k-1}^{j}&space;\frac{\lambda^{j&space;\rightarrow&space;i}}{\lambda_{p}^{j}}&space;\frac{1-\exp&space;\left(-\lambda_{p}^{i}&space;\tau\right)}{\lambda_{p}^{i}&space;\tau},&space;\quad(j=1,2,&space;\ldots,&space;k)" title="y_{k}^{i}(\tau)=y_{k-1}^{i}(\tau)+\sum_{j \neq i} d y_{k-1}^{j} \frac{\lambda^{j \rightarrow i}}{\lambda_{p}^{j}} \frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}, \quad(j=1,2, \ldots, k)" /></a>


####  ``` 6 ``` 
<a href="https://www.codecogs.com/eqnedit.php?latex=d&space;y_{k}^{i}=\sum_{j&space;\neq&space;i}&space;d&space;y_{k-1}^{j}&space;\frac{\lambda^{j-i}}{\lambda_{p}^{j}}\left[1-\frac{1-\exp&space;\left(-\lambda_{p}^{i}&space;\tau\right)}{\lambda_{p}^{i}&space;\tau}\right],&space;\quad(j=1,2,&space;\ldots,&space;k)" target="_blank"><img src="https://latex.codecogs.com/png.latex?d&space;y_{k}^{i}=\sum_{j&space;\neq&space;i}&space;d&space;y_{k-1}^{j}&space;\frac{\lambda^{j-i}}{\lambda_{p}^{j}}\left[1-\frac{1-\exp&space;\left(-\lambda_{p}^{i}&space;\tau\right)}{\lambda_{p}^{i}&space;\tau}\right],&space;\quad(j=1,2,&space;\ldots,&space;k)" title="d y_{k}^{i}=\sum_{j \neq i} d y_{k-1}^{j} \frac{\lambda^{j-i}}{\lambda_{p}^{j}}\left[1-\frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}\right], \quad(j=1,2,\ldots,k)"/></a>

where τ is the time step; λ<sup>j→i</sup> is the rate of formation of the i th nuclide from the j th nuclide, taking into account the probability of such a process (with the possibility of branching); and the speed of withdrawal of the nuclide j and i respectively at the expense of all processes.

The final solution is:

####  ``` 7 ```
<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%26y_%7Bk%7D%5E%7Bi%7D(%5Ctau)%3Dy_%7Bk-1%7D%5E%7Bi%7D(%5Ctau)%2B%5C%5C%0A%26%5Cbegin%7Barray%7D%7Bl%7D%0A%2B%5Csum_%7Bj_%7Bk%7D%20%5Cneq%20i%7D%20%5Cfrac%7B%5Clambda%5E%7Bj_%7Bk%7D%20%5Crightarrow%20i%7D%7D%7B%5Clambda_%7Bp%7D%5E%7Bj_%7Bk%7D%7D%7D%20%5Cfrac%7B1-%5Cexp%20%5Cleft(-%5Clambda_%7Bp%7D%5E%7Bi%7D%20%5Ctau%5Cright)%7D%7B%5Clambda_%7Bp%7D%5E%7Bi%7D%20%5Ctau%7D%5Cleft(1-%5Cfrac%7B1-%5Cexp%20%5Cleft(-%5Clambda_%7Bp%7D%5E%7Bi%7D%20%5Ctau%5Cright)%7D%7B%5Clambda_%7Bp%7D%5E%7Bi%7D%20%5Ctau%7D%5Cright)%5E%7Bk-1%7D%20%5Cprod_%7Bm%3D1%7D%5E%7Bk-1%7D%5Cleft(%5Csum_%7Bj_%7Bn%7D%20%5Cneq%20i%7D%20%5Cfrac%7B%5Clambda%5E%7Bj_%7Bm%7D%20%5Crightarrow%20i%7D%7D%7B%5Clambda_%7Bp%7D%5E%7Bj_%7Bn%7D%7D%7D%5Cright)%20d%20y_%7B0%7D%5E%7Bj_%7Bk%7D%7D%20%5Capprox%20%5C%5C%0A%5Capprox%20y_%7B0%7D%5E%7Bi%7D(%5Ctau)%2B%5Csum_%7Bj%20%5Cneq%20i%7D%5Cleft%5C%7Bd%20y_%7B0%7D%5E%7Bj%7D%20%5Cfrac%7B%5Clambda%5E%7Bj%20%5Crightarrow%20i%7D%7D%7B%5Clambda_%7Bp%7D%5E%7Bj%7D%7D%20%5Cfrac%7B1-%5Cexp%20%5Cleft(-%5Clambda_%7Bp%7D%5E%7Bi%7D%20%5Ctau%5Cright)%7D%7B%5Clambda_%7Bp%7D%5E%7Bi%7D%20%5Ctau%7D%5Cleft%5B1%2B%5Cleft(1-%5Cfrac%7B1-%5Cexp%20%5Cleft(-%5Clambda_%7Bp%7D%5E%7Bi%7D%20%5Ctau%5Cright)%7D%7B%5Clambda_%7Bp%7D%5E%7Bi%7D%20%5Ctau%7D%5Cright)%20%5Csum_%7Bj%20%5Cneq%20i%7D%20%5Cfrac%7B%5Clambda%5E%7Bj%20%5Crightarrow%20i%7D%7D%7B%5Clambda_%7Bp%7D%5E%7Bj%7D%7D%2B%5Cright.%5Cright.%0A%5Cend%7Barray%7D%5C%5C%0A%26%5Cleft.%5Cleft.%2B%5Cldots%2B%5Cleft(1-%5Cfrac%7B1-%5Cexp%20%5Cleft(-%5Clambda_%7Bp%7D%5E%7Bi%7D%20%5Ctau%5Cright)%7D%7B%5Clambda_%7Bp%7D%5E%7Bi%7D%20%5Ctau%7D%5Cright)%5E%7Bk-1%7D%20%5Cprod_%7Bm%3D1%7D%5E%7Bk-1%7D%5Cleft(%5Csum_%7Bj_%7Bn%7D%20%5Cneq%20i%7D%20%5Cfrac%7B%5Clambda%5E%7Bj_%7Bn%7D-i%7D%7D%7B%5Clambda_%7Bp%7D%5E%7Bj_%7Bn%7D%7D%7D%5Cright)%5Cright%5D%5Cright%5C%7D%0A%5Cend%7Baligned%7D"/>

Equation (5) describes the concentration of a nuclide, given the possible departure newcomers nuclei of that nuclide at the end of the time step. Equation (6) takes into account the formation of new nuclei of a nuclide of the other, broken at the k-th iteration.
The algorithm takes into account two main types of channels of nuclear transformations: radioactive decay of nuclei and nuclear reactions initiated by neutrons.

The rate of nuclear reactions for the i-th nuclide is described by the expression:
####  ``` 8 ``` 
<a href="https://www.codecogs.com/eqnedit.php?latex=\lambda_{r}=\frac{1}{V}\iint_{V_{F}}\sigma_{r}(E,\vec{r})\varphi(E,\vec{r})dEdV" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lambda_{r}=\frac{1}{V}\iint_{V_{F}}\sigma_{r}(E,\vec{r})\varphi(E,\vec{r})dEdV" title="\lambda_{r}=\frac{1}{V}\iint_{V_{F}}\sigma_{r}(E,\vec{r})\varphi(E,\vec{r})dEdV" /></a>


where λ is the reaction rate, V – volume of the computational cell, σ is micro reactions caused by neutrons with energy E in the point φ is the flux density of neutrons.
For fission expression (8) takes the form:
####  ``` 9 ``` 
<a href="https://www.codecogs.com/eqnedit.php?latex=\lambda_{fi}=\frac{1}{V}\iint_{VE}\sigma_{f&space;i}(E,\vec{r})\varphi(E,\vec{r})dEdV" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lambda_{fi}=\frac{1}{V}\iint_{VE}\sigma_{f&space;i}(E,\vec{r})\varphi(E,\vec{r})dEdV" title="\lambda_{fi}=\frac{1}{V}\iint_{VE}\sigma_{f i}(E,\vec{r})\varphi(E,\vec{r})dEdV" /></a>

For speed capture process λ<sub>ci</sub>:

####  ``` 10 ``` 
<a href="https://www.codecogs.com/eqnedit.php?latex=\lambda_{c&space;i}=\frac{1}{V}&space;\iint_{V&space;E}&space;\sigma_{c&space;i}(E,\vec{r})\varphi(E,\vec{r})dEdV" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lambda_{c&space;i}=\frac{1}{V}&space;\iint_{V&space;E}&space;\sigma_{c&space;i}(E,\vec{r})\varphi(E,\vec{r})dEdV" title="\lambda_{c i}=\frac{1}{V} \iint_{V E} \sigma_{c i}(E,\vec{r})\varphi(E,\vec{r})dEdV" /></a>

When calculating the speeds of neutron reactions assumes the immutability of the absolute density of the neutron flow on the time interval, i.e. throughout step τ used constant speed processes.
The transfer speed of the i th nuclide in the j th nuclide in the radioactive decay with half-lives T <sub>½</sub> and the probability of decay channel of the ε<sub>i→j</sub> is described by the expression:

####  ``` 11 ``` 
<a href="https://www.codecogs.com/eqnedit.php?latex=\lambda^{i&space;\rightarrow&space;j}=\varepsilon_{i&space;\rightarrow&space;j}&space;\frac{\ln&space;(2)}{T_{1&space;/&space;2}^{i}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lambda^{i&space;\rightarrow&space;j}=\varepsilon_{i&space;\rightarrow&space;j}&space;\frac{\ln&space;(2)}{T_{1&space;/&space;2}^{i}}" title="\lambda^{i \rightarrow j}=\varepsilon_{i \rightarrow j} \frac{\ln (2)}{T_{1 / 2}^{i}}" /></a>

Given the properties of the expression
####  ``` 12 ``` 
<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{aligned}&space;&\left(1-\frac{1-\exp&space;\left(-\lambda_{p}^{i}&space;\tau\right)}{\lambda_{p}^{i}&space;\tau}\right)<1,&space;\quad&space;\lambda_{p}^{i}&space;\rightarrow&space;\infty\\&space;&\left(1-\frac{1-\exp&space;\left(-\lambda_{p}^{i}&space;\tau\right)}{\lambda_{p}^{i}&space;\tau}\right)&space;\rightarrow&space;0,&space;\quad&space;\lambda_{p}^{i}&space;\rightarrow&space;0&space;\end{aligned}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{aligned}&space;&\left(1-\frac{1-\exp&space;\left(-\lambda_{p}^{i}&space;\tau\right)}{\lambda_{p}^{i}&space;\tau}\right)<1,&space;\quad&space;\lambda_{p}^{i}&space;\rightarrow&space;\infty\\&space;&\left(1-\frac{1-\exp&space;\left(-\lambda_{p}^{i}&space;\tau\right)}{\lambda_{p}^{i}&space;\tau}\right)&space;\rightarrow&space;0,&space;\quad&space;\lambda_{p}^{i}&space;\rightarrow&space;0&space;\end{aligned}" title="\begin{aligned} &\left(1-\frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}\right)<1, \quad \lambda_{p}^{i} \rightarrow \infty\\ &\left(1-\frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}\right) \rightarrow 0, \quad \lambda_{p}^{i} \rightarrow 0 \end{aligned}" /></a>

you can verify that the formation of "new" nuclei in (6) at the end of the time step tend to zero, which ensures convergence of the iteration process.
Moreover, the solution of (5) – (6) is always non-negative, i.e. the solution of the system of equations (1) – (2) y<sub>i</sub> exists and it is always positive, since in the sum (5) is always at least one summand is different from zero because the initial concentration of at least one nuclide is positive (otherwise, the task loses physical meaning).
The iterative process of solving equations (5) – (6) continues until the condition is met:
####  ``` 13 ``` 
<a href="https://www.codecogs.com/eqnedit.php?latex=\left|1-\frac{y_{n}^{k}}{y_{n-1}^{k}}\right|&space;\leq&space;\delta_{\max&space;},&space;\quad&space;\forall&space;k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left|1-\frac{y_{n}^{k}}{y_{n-1}^{k}}\right|&space;\leq&space;\delta_{\max&space;},&space;\quad&space;\forall&space;k" title="\left|1-\frac{y_{n}^{k}}{y_{n-1}^{k}}\right| \leq \delta_{\max }, \quad \forall k" /></a>

where δ<sub>max</sub> is the maximum allowable value of the variation of the concentrations of nuclides at two adjacent iterations specified by the user in the calculation options. The accuracy of the iterative solution the default is 10-5, but can be changed by the user.
From equation (7) it follows that when k → ∞ the last term in the series is actually a product of two multiplicands, each of which is raised to the power of (k - 1) by the number of obviously smaller units and tend to zero, which ensures convergence. A number of (7) can be interpreted as the Neumann series, convergence has been proven.
During the iterative process operated with only non-negative values of equations (5) and (6). Thus, a solution exists, it is positive, because the equation (5) is always one summand is different from zero.
The calculation of the decay heat  is inextricably linked with the concentration of radioactive nuclides in the material and the energy released during radioactive decay of the nucleus, and is determined by the formula:
####  ``` 14 ``` 
<a href="https://www.codecogs.com/eqnedit.php?latex=W_{o&space;c&space;m}=1,60219&space;\cdot&space;10^{-13}&space;\cdot&space;\sum_{j}&space;y_{j}&space;\lambda_{j}&space;E_{j}" target="_blank"><img src="https://latex.codecogs.com/png.latex?W_{o&space;c&space;m}=1,60219&space;\cdot&space;10^{-13}&space;\cdot&space;\sum_{j}&space;y_{j}&space;\lambda_{j}&space;E_{j}" title="W_{o c m}=1,60219 \cdot 10^{-13} \cdot \sum_{j} y_{j} \lambda_{j} E_{j}" /></a>

where Ej is the heat generation by the decay of nuclide j, MeV/dis.; λj – decay constant of j-th radioactive nuclide with-1; 1,60219 ∙ 10-13 – the conversion factor from MeV to watts.
The calculation of activity of nuclide Aj(t) at time t is a product of its nuclear density yi(t) by a constant decay in the whole volume of the computational cell:
####  ``` 15 ``` 
<a href="https://www.codecogs.com/eqnedit.php?latex=A_{j}(t)=\lambda_{j}&space;\int_{0}^{V}&space;y_{j}(\vec{r})&space;d&space;v" target="_blank"><img src="https://latex.codecogs.com/png.latex?A_{j}(t)=\lambda_{j}&space;\int_{0}^{V}&space;y_{j}(\vec{r})&space;d&space;v" title="A_{j}(t)=\lambda_{j} \int_{0}^{V} y_{j}(\vec{r}) d v" /></a>
where r - the coordinate of the calculated point, V - its volume.

- ## Chebyshev Rational Approximation 
The Chebyshev rational approximation method (CRAM) is a relatively straight- forward algorithm. A rational function 𝑟^𝑘,𝑘(𝑥) is found that minimizes the max- imum error with regard to the scalar exponent along the negative real axis
The defining equation is Equation (21), where 𝜋𝑘,𝑘 is the set of all rational functions with numerators and denominators of order 𝑘. As 𝑘 increases, the accuracy of the approximation also increases.

####  ``` 16 ``` 
<a href="https://www.codecogs.com/eqnedit.php?latex=\sup&space;_{x&space;\in&space;\mathbb{R}_{-}}\left|\hat{r}_{k,&space;k}(x)-e^{x}\right|=\inf&space;_{r_{k,&space;k}&space;\in&space;\pi_{k,&space;k}}\left\{\sup&space;_{x&space;\in&space;\mathbb{R}_{-}}\left|r_{k,&space;k}(x)-e^{x}\right|\right\}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\sup&space;_{x&space;\in&space;\mathbb{R}_{-}}\left|\hat{r}_{k,&space;k}(x)-e^{x}\right|=\inf&space;_{r_{k,&space;k}&space;\in&space;\pi_{k,&space;k}}\left\{\sup&space;_{x&space;\in&space;\mathbb{R}_{-}}\left|r_{k,&space;k}(x)-e^{x}\right|\right\}" title="\sup _{x \in \mathbb{R}_{-}}\left|\hat{r}_{k, k}(x)-e^{x}\right|=\inf _{r_{k, k} \in \pi_{k, k}}\left\{\sup _{x \in \mathbb{R}_{-}}\left|r_{k, k}(x)-e^{x}\right|\right\}" /></a>

Once the function 𝑟^𝑘,𝑘(𝑥) is known, it can be rearranged to reduce costs further or to improve numerical stability. The incomplete partial fraction (IPF) form, shown in Equation (16), is a good combination of numerical stability and efficiency. The values 𝛼𝑙 and 𝜃𝑙 are tabulated and are available for a variety of values of 𝑘 up to 48. In the IPF form, only sparse matrix solves are necessary to compute the action on a vector.

####  ``` 17 ``` 
<a href="https://www.codecogs.com/eqnedit.php?latex=\hat{r}_{k,&space;k}(x)=\alpha_{0}&space;\prod_{l=1}^{k&space;/&space;2}\left(1&plus;2&space;\Re\left\{\frac{\tilde{\alpha}_{l}}{x-\theta_{l}}\right\}\right)" target="_blank"><img src="https://latex.codecogs.com/png.latex?\hat{r}_{k,&space;k}(x)=\alpha_{0}&space;\prod_{l=1}^{k&space;/&space;2}\left(1&plus;2&space;\Re\left\{\frac{\tilde{\alpha}_{l}}{x-\theta_{l}}\right\}\right)" title="\hat{r}_{k, k}(x)=\alpha_{0} \prod_{l=1}^{k / 2}\left(1+2 \Re\left\{\frac{\tilde{\alpha}_{l}}{x-\theta_{l}}\right\}\right)" /></a>

CRAM is both efficient and highly accurate over the domain in which it is derived. However, eigenvalues with extremely large imaginary components or positive real components will reduce the accuracy. As such, CRAM is not recommended for use in highly oscillatory problems or those with possible exponential growth such as reactor dynamics.
## References
1.    E.F. Seleznev, A.A. Belov, V.I. Belousov, I.S.Chernova BPSD code upgrade for solving the nuclear kinetics problem [https://doi.org/10.26583/npe.2018.4.10]
2.    M. Pusa Rational Approximations to the Matrix Exponential in Burnup Calculation, Nucl. Sci. Eng., 169, 2, p.155-167 (2011)
3.    H. Bateman: Proc.Cambridge Phil.Soc.15, 15, 423 (1910)
4.    Paul K. Romano, Nicholas E. Horelik, Bryan R. Herman, Adam G. Nelson, Benoit Forget, and Kord Smith, “OpenMC: A State-of-the-Art Monte Carlo Code for Research and Development,” Ann. Nucl. Energy, 82, 90–97 (2015).
5.    M. Fleming, J.C. Sublet, J. Kopecky, D. Rochman and A.J. Koning, "Probing experimental and systematic trends of the neutron-induced TENDL-2014 nuclear data library", CCFE report UKAEA-R(15)29, October 2015
