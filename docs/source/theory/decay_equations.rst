.. _decay_equation:

-----------------------
Nuclide decay equations
-----------------------

The system of equations nuclide kinetics is heterogeneous linear system 
of equations:

.. math::
    :label: decay

    \frac{dy_{i}}{dt}=\sum_{j=1}^{n}a_{i;j}*y_{j}+q_{i};\quad(i=1,2,\dots,n)

with the initial condition:

.. math::
    :label: decay_ic1

    y_{i}(0)=y_{i0},\quad(i=1,2,\dots,n)

where y\ :sub:`i`\ is the concentration of the `i` th nuclide, a\ :sub:`ij`\ are
the coefficients characterizing the channels of the transformations of the `i` th
nuclide from the `j`–th nuclide,q\ :sub:`i`\ an external source,
`n` is the number of nuclides, y\ :sub:`i0`\  is the concentration of the `i` th
nuclide at time t\ :sub:`0`\.

.. _decay_heat:

----------
Decay heat
----------

The calculation of the decay heat  is inextricably linked with the concentration
of radioactive nuclides in the material and the energy released during
radioactive decay of the nucleus, and is determined by the expression:

.. math::
    :label: decay_ic2

    W_{dh}=1,60219\cdot10^{-13}\cdot\sum_{j}y_{j}\lambda_{j}E_{j}

where E\ :sub:`j`\ is the heat generation by the decay of nuclide `j`, MeV/rs.;
λ\ :sub:`j`\ – decay constant of `j`-th radioactive nuclide with-1;
1,60219 ∙ 10-13 – the conversion factor from MeV to watts.
The calculation of activity of nuclide A\ :sub:`j`\(t) at time t is a product of
its nuclear density y\ :sub:`i`\(t) by a constant decay in the whole volume
of the computational cell:

.. math::
    :label: decay_ic3

    A_{j}(t)=\lambda_{j}\int_{0}^{V}y_{j}(\vec{r})dv

where `r` - the coordinate of the calculated point, `V` - its volume.

.. _iteration_method:

----------------
Iteration method
----------------

The system of equations :eq:`decay` – :eq:`decay_ic1` can be solved analytically
as the sum of its particular solution and the General solution of the 
corresponding homogeneous linear system of equations. In this case, there is a 
needed matrix of size n × n dimension. Therefore, this method of solution is
often used for a partial transition matrices, for example, triangular.
Another widely known method of solving the problem is to decompose the solution 
in the range of the exponential function, which leads to the need to use a 
recurrence relation. In the code the solution of equations :eq:`decay` – 
:eq:`decay_ic1` is the iterative method. The `k` th iteration could be written
as:

.. math::
    :label: iteration1

    y_{k}^{i}(\tau) = y_{k-1}^{i}(\tau) + 
    \sum_{j\neq{i}}dy_{k-1}^{j}\frac{\lambda^{j\rightarrow{i}}}{\lambda_{p}^{j}}\frac{1-\exp\left(-\lambda_{p}^{i}\tau\right)}{\lambda_{p}^{i}\tau},\quad(j=1,2,\ldots,k)

.. math::
    :label: iteration2

    y_{k}^{i}=\sum_{j\neq{i}}dy_{k-1}^{j}\frac{\lambda^{j-i}}{\lambda_{p}^{j}}\left[1-\frac{1-\exp\left(-\lambda_{p}^{i}\tau\right)}{\lambda_{p}^{i}\tau}\right],\quad(j=1,2,\ldots,k)

where τ is the time step; λ\ :sup:`j→i`\ is the rate of formation of the `i` th
nuclide from the `j` th nuclide, taking into account the probability of such a
process (with the possibility of branching) and the speed of disappearence of 
the nuclide `j` and `i` respectively at the expense of all processes.

The final solution is:

.. math::
    :label: iteration3

    \begin{aligned}
    &y_{k}^{i}(\tau)=y_{k-1}^{i}(\tau)+\\
    &\begin{array}{l}
    +\sum_{j_{k} \neq i} \frac{\lambda^{j_{k} \rightarrow i}}{\lambda_{p}^{j_{k}}} \frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}\left(1-\frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}\right)^{k-1} \prod_{m=1}^{k-1}\left(\sum_{j_{n} \neq i} \frac{\lambda^{j_{m} \rightarrow i}}{\lambda_{p}^{j_{n}}}\right) d y_{0}^{j_{k}} \approx \\\approx y_{0}^{i}(\tau)+\sum_{j \neq i}\left\{d y_{0}^{j} \frac{\lambda^{j \rightarrow i}}{\lambda_{p}^{j}} \frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}\left[1+\left(1-\frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}\right) \sum_{j \neq i} \frac{\lambda^{j \rightarrow i}}{\lambda_{p}^{j}}+
    \right.\right.
    \end{array}\\
    &\left.\left.+\ldots+\left(1-\frac{1-\exp \left(-\lambda_{p}^{i} \tau\right)}{\lambda_{p}^{i} \tau}\right)^{k-1} \prod_{m=1}^{k-1}\left(\sum_{j_{n} \neq i} \frac{\lambda^{j_{n}-i}}{\lambda_{p}^{j_{n}}}  \right)\right]\right\}\end{aligned}


Equation :eq:`iteration1` describes the concentration of a nuclide, given the 
possible expiration of newcomers nuclei of that nuclide at the end of the time
step. Equation :eq:`iteration2` takes into account the formation of new nuclei
from other nuclides decayed at the `k`-th iteration. The algorithm takes into
account two main types of channels of nuclear transformations: nuclide
radioactive decay and nuclear reactions induced by external particle source.

The rate of nuclear reactions for the `i` -th nuclide is described by the expression:

.. math::
    :label: iteration4
    
    \lambda_{r}=\frac{1}{V}\iint_{V_{F}}\sigma_{r}(E,\vec{r})\varphi(E,\vec{r})dEdV

where λ is the reaction rate, 
`V` – volume of the computational cell, σ is micro reactions caused by particles
with energy `E` in the point φ is the flux density of external particle source.
For neutron external source fission above expression takes the form:

.. math::
    :label: iteration5
    
    \lambda_{fi}=\frac{1}{V}\iint_{VE}\sigma_{f{i}}(E,\vec{r})\varphi(E,\vec{r})dEdV

For radiactive capture reaction rate λ\ :sub:`ci`\ :

.. math::
    :label: iteration6

    \lambda_{c{i}}=\frac{1}{V}\iint_{VE}\sigma_{c{i}}(E,\vec{r})\varphi(E,\vec{r})dEdV
    
When calculating the reaction rates assumes the immutability of the absolute 
flux density on the time interval, i.e. throughout step τ used constant 
processes rate. The transfer rate of the `i` th nuclide to the `j` th nuclide
in the radioactive decay with half-lives `T` \ :sub:`½`\ and the probability
of decay channel of the ε<\ :sub:`i→j`\ is described by the expression:

.. math::
    :label: iteration7

    \lambda^{i\rightarrow{j}}=\varepsilon_{i\rightarrow{j}}\frac{\ln(2)}{T_{1/2}^{i}}

Given the properties of the expression:

.. math::
    :label: iteration8
  
    \begin{aligned}&\left(1-\frac{1-\exp\left(-\lambda_{p}^{i}\tau\right)}{\lambda_{p}^{i}\tau}\right)<1,\quad\lambda_{p}^{i}\rightarrow\infty\\&\left(1-\frac{1-\exp\left(-\lambda_{p}^{i}\tau\right)}{\lambda_{p}^{i}\tau}\right)\rightarrow0,\quad\lambda_{p}^{i}\rightarrow0\end{aligned}

you can verify that the formation of "new" nuclei in :eq:`iteration2` at the end
of the time step tend to zero, which ensures convergence of the iteration 
process. Moreover, the solution of :eq:`iteration1` – :eq:`iteration2` is always
non-negative, i.e. the solution of the system of equations :eq:`decay` – 
:eq:`decay_ic1` y\ :sub:`i`\ exists and it is always positive, since in the sum
:eq:`iteration1` is always at least one summ and is different from zero because
the initial concentration of at least one nuclide is positive (otherwise,
the task loses physical meaning). The iterative process of solving equations
:eq:`iteration1` – :eq:`iteration2` continues until the condition is met:

.. math::
    :label: iteration9

    \left|1-\frac{y_{n}^{k}}{y_{n-1}^{k}}\right|\leq\delta_{\max},\quad\forall k

where δ\ :sub:`max`\ is the maximum acceptable value of the variation of the
nuclides concentrations at two adjacent iterations specified by the user in the
calculation options. The accuracy of the iterative solution the default is 10-3,
but can be changed by the user. From equation :eq:`iteration3` it follows that
when `k` → ∞ the last term in the series is actually a product of two 
multiplicands, each of which is raised to the power of (k - 1) by the number of
obviously smaller units and tend to zero, which ensures convergence. A number of
:eq:`iteration3` can be interpreted as the Neumann series, convergence has been
proven. During the iterative process operated with only non-negative values of
equations :eq:`iteration1` and :eq:`iteration2`. Thus, a solution exists, it is
positive, because the equation :eq:`iteration1` is always one summand is
different from zero.

.. _chebyshev_method:

--------------------------------
Chebyshev Rational Approximation
--------------------------------

The Chebyshev rational approximation method (CRAM) is a relatively straight-
forward algorithm. A rational function 𝑟^𝑘,𝑘(𝑥) is found that minimizes the 
maximum error with regard to the scalar exponent along the negative real axis
The defining equation is Equation :eq:`cram1`, where 𝜋𝑘,𝑘 is the set of all 
rational functions with numerators and denominators of order 𝑘. As 𝑘 increases,
the accuracy of the approximation also increases.

.. math::
    :label: cram1

    \sup_{x\in\mathbb{R}_{-}}\left|\hat{r}_{k,k}(x)-e^{x}\right|=\inf_{r_{k,k}\in\pi_{k,k}}\left\{\sup_{x\in\mathbb{R}_{-}}\left|r_{k,k}(x)-e^{x}\right|\right\}

Once the function 𝑟^𝑘,𝑘(𝑥) is known, it can be rearranged to reduce costs
further or to improve numerical stability. The values 𝛼𝑙 and 𝜃𝑙 are tabulated
and are available for a variety of values of 𝑘 up to 48. In the IPF form, only
sparse matrix solves are necessary to compute the action on a vector.

.. math::
    :label: cram2

    \hat{r}_{k,k}(x)=\alpha_{0}\prod_{l=1}^{k/2}\left(1 + 2\Re\left\{\frac{\tilde{\alpha}_{l}}{x-\theta_{l}}\right\}\right)

.. _baetman_method:

--------------
Baetman method
--------------

IN PROGRESS...
    






    







   
