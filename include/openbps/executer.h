#ifndef SRC_EXECUTER_H_
#define SRC_EXECUTER_H_
#include <complex>
#include "xtensor/xarray.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xcomplex.hpp"
#include "xtensor/xbuilder.hpp"
#include "configure.h"
#include "matrix.h"
#include "chain.h"

namespace openbps {

namespace executer {

//==============================================================================
// Global variables
//==============================================================================
//! CRAM constants
extern xt::xarray<std::complex<double>> theta16; //!< Pole values CRAM 8th order
extern xt::xarray<std::complex<double>> alpha16; //!< Resides values CRAM 8th
extern xt::xarray<std::complex<double>> theta48; //!< Pole values CRAM 24th ord
extern xt::xarray<std::complex<double>> alpha48; //!< Resides values CRAM 24th
extern xt::xarray<double> theta_48r;             //!< Real part value of CRAM
                                                 //!< 24th order poles
extern xt::xarray<double> theta_48i;             //!< Imag part value of CRAM
                                                 //!< 24th order poles
extern xt::xarray<double> alpha_48r;             //!< Real part value of CRAM
                                                 //!< 24th order resides
extern xt::xarray<double> alpha_48i;             //!< Imag part value of CRAM
                                                 //!< 24th order resides
extern double alpha480;                          //!< CRAM coefficient 24th ord
extern double alpha160;                          //!< CRAM coefficient 8th ord

//==============================================================================
// Non class methods
//==============================================================================

//! Run a calculation process
//!
void run_solver();

//! Initialize a calcaulation process
//!
void init_solver();

void executethr(int imat, Chain& chain, DecayMatrix& dm, DecayMatrix& ddm);
//! Iterative method implementation v.1.
//! Method based on algorythm described at paper Seleznev E.F., Belov A.A.,
//! Belousov V.I., Chernova I.S. "BPSD CODE UPGRADE FOR SOLVING
//! THE NUCLEAR KINETICS PROBLEM" Izvestiya VUZov: Yadernaya energetica N 4 2018 
//! pp 115-127. The main advantage is ability to get results uncertainty in the
//! main calculation process. This overloading method is without uncertainty 
//! analysis.
//
//!\param[in] matrix a nuclide transition matrix with zero values 
//! diagonal elements
//!\param[in] sigp vector with positive diagonal main decay matrix elements
//!\param[inout] y nuclear concentrations vector
void iterative(xt::xarray<double>& matrix, xt::xarray<double>& sigp,
               xt::xarray<double>& y);

//! Iterative method implementation v.2.
//! Method based on algorythm described at paper Seleznev E.F., Belov A.A.,
//! Belousov V.I., Chernova I.S. "BPSD CODE UPGRADE FOR SOLVING
//! THE NUCLEAR KINETICS PROBLEM" Izvestiya VUZov: Yadernaya energetica N 4 2018
//! pp 115-127. The main advantage is ability to get results uncertainty in the
//! main calculation process. This overloading method is with uncertainty 
//! analysis.
//
//!\param[in] matrix a nuclide transition matrix with zero values 
//! diagonal elements
//!\param[in] sigp vector with positive diagonal main decay matrix elements
//!\param[inout] y nuclear concentrations vector
//!\param[in] dmatrix a deviative part of nuclides transition matrix
//!\param[in] dsigp vector with deviative part of diagonal matrix elements
//!\param[inout] dy uncertainty values of nuclear concentration 
void iterative (xt::xarray<double>& matrix, xt::xarray<double>& sigp, 
                xt::xarray<double>& y,
                xt::xarray<double>& dmatrix, xt::xarray<double>& dsigp,
                xt::xarray<double>& dy);

//! Chebyshev rational appproximation method (CRAM)
//! Method are based on work published for the Chebyshev Rational Approximation
//! Method (CRAM), as described in the following paper: M. Pusa, "`Higher-Order
//! Chebyshev Rational Approximation Method and Application to Burnup Equations
//! <https://doi.org/10.13182/NSE15-26>`_," Nucl. Sci. Eng., 182:3, 297-318.and
//! its method implementation in OpenMC project on SciPy Python library:
//! https://github.com/openmc-dev/openmc/blob/develop/openmc/deplete/cram.py
//!
//!
//!\param[in] matrix a full nuclide transition matrix
//!\param[inout] y nuclear concentrations vector
//!\param[in] alpha resides constatnts depending on input value 'order'
//!\param[in] theta resides constants depending on input value 'order'
//!\param[in] order of CRAM
//!\param[in] alpha0 real constant with value depending on input value 'order'
void cram(xt::xarray<double>& matrix, xt::xarray<double>& y,
          xt::xarray<std::complex<double>>& alpha,
          xt::xarray<std::complex<double>>& theta,
          int order, double alpha0);

//! Apply filters to result
//!
//! \param[in] chainer chain object
//! \param[in] matname name of material
void apply_filters(const Chain &chainer, const std::string& matname);

} // namespace executer

} // namespace openbps

#endif /* SRC_EXECUTER_H_ */
