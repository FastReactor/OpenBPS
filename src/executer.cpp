#include "openbps/executer.h"
#include <vector>
#include <thread>
#include <cmath>
#include <iostream>
#include <fstream>
#include <complex>
#include "xtensor/xarray.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xcomplex.hpp"
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include "openbps/parse.h"
#include "openbps/timeproc.h"
#include "openbps/uncertainty.h"
#include "openbps/nuclide.h"
#include "openbps/chain.h"
#include "openbps/materials.h"
#include "openbps/configure.h"
#include "openbps/functionals.h"
#include "openbps/reactions.h"
#include "openbps/filter.h"
#include "openbps/matrix.h"

//==============================================================================
// Type definitions
//==============================================================================
typedef Eigen::SparseMatrix<double> SpMat;

using namespace std::complex_literals;

namespace openbps {

namespace executer {

//==============================================================================
// Global variables
//==============================================================================
//! CRAM constants description
xt::xarray<std::complex<double>> theta16 {
    +3.509103608414918 + 8.436198985884374i,
    +5.948152268951177 + 3.587457362018322i,
    -5.264971343442647 + 16.22022147316793i,
    +1.419375897185666 + 10.92536348449672i,
    +6.416177699099435 + 1.194122393370139i,
    +4.993174737717997 + 5.996881713603942i,
    -1.413928462488886 + 13.49772569889275i,
    -10.84391707869699 + 19.27744616718165i};

xt::xarray<std::complex<double>> alpha16{
    +5.464930576870210e+3 - 3.797983575308356e+4i,
    +9.045112476907548e+1 - 1.115537522430261e+3i,
    +2.344818070467641e+2 - 4.228020157070496e+2i,
    +9.453304067358312e+1 - 2.951294291446048e+2i,
    +7.283792954673409e+2 - 1.205646080220011e+5i,
    +3.648229059594851e+1 - 1.155509621409682e+2i,
    +2.547321630156819e+1 - 2.639500283021502e+1i,
    +2.394538338734709e+1 - 5.650522971778156e+0i};

xt::xarray<double> theta_48r{
    -4.465731934165702e+1, -5.284616241568964e+0,
    -8.867715667624458e+0, +3.493013124279215e+0,
    +1.564102508858634e+1, +1.742097597385893e+1,
    -2.834466755180654e+1, +1.661569367939544e+1,
    +8.011836167974721e+0, -2.056267541998229e+0,
    +1.449208170441839e+1, +1.853807176907916e+1,
    +9.932562704505182e+0, -2.244223871767187e+1,
    +8.590014121680897e-1, -1.286192925744479e+1,
    +1.164596909542055e+1, +1.806076684783089e+1,
    +5.870672154659249e+0, -3.542938819659747e+1,
    +1.901323489060250e+1, +1.885508331552577e+1,
    -1.734689708174982e+1, +1.316284237125190e+1};

xt::xarray<double> theta_48i{
    +6.233225190695437e+1, +4.057499381311059e+1,
    +4.325515754166724e+1, +3.281615453173585e+1,
    +1.558061616372237e+1, +1.076629305714420e+1,
    +5.492841024648724e+1, +1.316994930024688e+1,
    +2.780232111309410e+1, +3.794824788914354e+1,
    +1.799988210051809e+1, +5.974332563100539e+0,
    +2.532823409972962e+1, +5.179633600312162e+1,
    +3.536456194294350e+1, +4.600304902833652e+1,
    +2.287153304140217e+1, +8.368200580099821e+0,
    +3.029700159040121e+1, +5.834381701800013e+1,
    +1.194282058271408e+0, +3.583428564427879e+0,
    +4.883941101108207e+1, +2.042951874827759e+1};

xt::xarray<double> alpha_48r{
    +6.387380733878774e+2, +1.909896179065730e+2,
    +4.236195226571914e+2, +4.645770595258726e+2,
    +7.765163276752433e+2, +1.907115136768522e+3,
    +2.909892685603256e+3, +1.944772206620450e+2,
    +1.382799786972332e+5, +5.628442079602433e+3,
    +2.151681283794220e+2, +1.324720240514420e+3,
    +1.617548476343347e+4, +1.112729040439685e+2,
    +1.074624783191125e+2, +8.835727765158191e+1,
    +9.354078136054179e+1, +9.418142823531573e+1,
    +1.040012390717851e+2, +6.861882624343235e+1,
    +8.766654491283722e+1, +1.056007619389650e+2,
    +7.738987569039419e+1, +1.041366366475571e+2};

xt::xarray<double> alpha_48i{
    -6.743912502859256e+2, -3.973203432721332e+2,
    -2.041233768918671e+3, -1.652917287299683e+3,
    -1.783617639907328e+4, -5.887068595142284e+4,
    -9.953255345514560e+3, -1.427131226068449e+3,
    -3.256885197214938e+6, -2.924284515884309e+4,
    -1.121774011188224e+3, -6.370088443140973e+4,
    -1.008798413156542e+6, -8.837109731680418e+1,
    -1.457246116408180e+2, -6.388286188419360e+1,
    -2.195424319460237e+2, -6.719055740098035e+2,
    -1.693747595553868e+2, -1.177598523430493e+1,
    -4.596464999363902e+3, -1.738294585524067e+3,
    -4.311715386228984e+1, -2.777743732451969e+2};
xt::xarray<std::complex<double>> theta48{theta_48r+theta_48i*1.i};
xt::xarray<std::complex<double>> alpha48 {alpha_48r+alpha_48i*1.i};
double alpha480{2.258038182743983e-47};
double alpha160{2.124853710495224e-16};
//! CRAM constants description

//==============================================================================
// Non class methods implementation
//==============================================================================

//! Initialize a calculation process
void init_solver() {
    // Read a nuclides information
    read_nuclide_xml(configure::nuclide_file);
    // Read a materials
    read_materials_from_inp(configure::inmaterials_file);
    if (!configure::decay_extra_out) {
        // If mode is not decay only calculation
        read_reactions_xml();
        // Make bind between materials and compositions
        matchcompositions();
    }
}

void executethr(int imat, Chain& chain, DecayMatrix& dm, DecayMatrix& ddm) {
    //std::vector<std::unique_ptr<Materials>>::iterator mat = materials.begin() + imat;
//    auto mat = materials[imat];
    xt::xarray<double> dy;
    bool isMaterial {false};
    xt::xarray<double> y =
                make_concentration(chain, materials[imat]->namenuclides, materials[imat]->conc);
        // Prepare dump to store every step of calculation results
        if (configure::outwrite) {
            configure::dumpoutput.clear();
            configure::dumpoutput.resize(configure::numstep);
            for (int k = 0; k < configure::numstep; k++) {
                configure::dumpoutput[k].resize(y.size());
            }
        }
        switch (configure::calcmode) {
        // Analytical method
        case Mode::baetman:
            std::cout << "Baetman implementation in progress" << std::endl;
            break;
        // Iteration method
        case Mode::iteration:
        {
            IterMatrix im(dm);
            auto matrix = im.matrixreal(chain, *materials[imat]);
            auto sigp = im.sigp(chain, *materials[imat]);
            // Perform calculation with uncertanties analysis
            if (configure::uncertantie_mod) {
                dy = make_concentration(chain,materials[imat]->namenuclides,
                                        materials[imat]->conc, true);
                IterMatrix imim(ddm);
                auto dmatrix = imim.matrixdev(chain, *materials[imat]);
                auto dsigp = imim.dsigp(chain, *materials[imat]);
                iterative(matrix, sigp, y,
                          dmatrix, dsigp, dy);
            } else {
                iterative(matrix, sigp, y);
            }
        }
            break;
        // Chebyshev rational approximation method
        case Mode::chebyshev:
        {
            CramMatrix cm(dm);
            auto matrix = cm.matrixreal(chain, *materials[imat]);
            if (configure::order == 8) {
                std::vector<std::size_t> shape = {y.size()};
                xt::xarray<double> dy2  = xt::zeros<double>(shape);
                
                cram(matrix, y, alpha16, theta16,
                     configure::order, alpha160);
                if (configure::uncertantie_mod) {
                    configure::outwrite = false;
                    
                    CramMatrix devcm(ddm);
                    dy = make_concentration(chain, materials[imat]->namenuclides,
                                        materials[imat]->conc, true);
                    auto dmatrix = devcm.matrixdev(chain, *materials[imat]);
                    cram(matrix, dy, alpha16, theta16,
                     configure::order, alpha160);
                    
                    std::copy(y.begin(), y.end(), dy2.begin());
                    cram(dmatrix, dy2, alpha16, theta16,
                     configure::order, alpha160);
                    dy = dy + xt::abs(y - dy2);
                    configure::outwrite = true;
                    
                }
            } else {
                cram(matrix, y, alpha48, theta48,
                     configure::order, alpha480);
            }
        }
            break;
        }//switch case
        size_t j {0};
        for (auto& item : chain.name_idx) {
            if (configure::verbose)
                std::cout << item.first << " = " << y[j] << " err: "<<
                             dy[j]<< std::endl;
            if (configure::rewrite) {
                if (configure::uncertantie_mod) {
                    materials[imat]->add_nuclide(item.first, udouble(y[j], dy[j]));
                } else {
                    materials[imat]->add_nuclide(item.first, y[j]);
                }
            }
            j++;
        }
        // Write dump to the text *.csv file
        if (configure::outwrite) {
            // If materials fitler is present
            // then applying it
            if (materialfilter != nullptr) {
                materialfilter->apply(materials[imat]->Name(), isMaterial);
            } else {
                isMaterial = true;
            }
            if (isMaterial)
                apply_filters(chain, materials[imat]->Name());
         }
}

//! Run a calculation process
void run_solver() {
    bool isMaterial {false};
    if (configure::verbose)
        std::cout << "We start execution\n";
    // Read a chain from xml file
    Chain chain = read_chain_xml(configure::chain_file);
    // Form a calculation matrix
    DecayMatrix dm(chain.name_idx.size());
    dm.form_matrixreal(chain);
    DecayMatrix ddm(chain.name_idx.size());
    ddm.form_matrixdev(chain);
    std::vector<std::thread> threads;
    for (int i = 0; i < materials.size(); i++) {
        std::thread thr=std::thread(executethr, i, std::ref(chain), std::ref(dm), std::ref(ddm));
        threads.emplace_back(std::move(thr));
    }
    for(auto& thr : threads) {
        thr.join();
    }

    std::cout << "Done!" << std::endl;
    // Write down the getting nuclear concentration for every material
    // into *.xml file
    if (configure::rewrite)
        form_materials_xml(configure::outmaterials_file);
}

//! Iterative method implementation v.1. without uncertanties
void iterative(xt::xarray<double>& matrix, xt::xarray<double>& sigp,
               xt::xarray<double>& y){
    // Time step
    double dt {configure::timestep/configure::numstep};
    // Array shapes
    size_t N {y.size()};
    std::vector<std::size_t> shape = { y.size() };
    std::vector<std::size_t> dshape = { y.size(), y.size() };
    // Auxilary arrays declaration
    xt::xarray<double> ro = xt::zeros<double>(shape);
    xt::xarray<double> roo = xt::zeros<double>(shape);
    xt::xarray<double> rrr = xt::zeros<double>(shape);
    xt::xarray<double> arr = xt::zeros<double>(shape);
    xt::xarray<double> arrtemp = xt::zeros<double>(dshape);//?!
    xt::xarray<double> rest = xt::zeros<double>(shape);
    xt::xarray<double> disr = xt::zeros<double>(shape);
    xt::xarray<double> et = xt::ones<double>(shape);
    xt::xarray<double> er = xt::ones<double>(shape);
    double aa {0.0};
    // Implementation
    arr = sigp * dt;
    // Rest part of radioactive nuclide through decay for time dt
    rest = xt::exp(-arr);
    // Disapperead part of radioactive nuclide through decay for time dt
    disr = 1 - rest;
    // Fill auxilary arrays
    for (size_t i = 0; i < N; i++) {
        if (disr(i) < 1.e-10) disr(i) = arr(i);
        if (arr(i) > 0.0) er(i) = disr(i) / arr(i);
        for (size_t j = 0; j < N; j++) {
            if (sigp(j) > 0) matrix(i, j) = matrix(i, j) / sigp(j);
        }
    }
    et = 1.0 - er;
    roo = y;
    // Time iteration
    for (int k = 0; k < configure::numstep; k++) {
        bool proceed {true};
        int t = 0;
        rrr = y * disr;       //!< disappear at time dt
        y = y * rest;         //!< rest at time dt
        ro = y;
        // Main loop
        while (proceed) {
            t++;
            // Transition of nuclear concentration part from parent
            // to child nuclear
            for (size_t iparent=0; iparent < N; iparent++) {
                for (size_t ichild=0; ichild < N; ichild++) {
                    arrtemp(ichild, iparent) = matrix(ichild, iparent) *
                            rrr(iparent);
                }
            }
            // Additive nucleus through decay
            for (size_t ichild=0; ichild < N; ichild++) {
                aa = xt::sum(xt::row(arrtemp, ichild))(0);
                ro(ichild) += aa;
            }
            // Get additives to child nucleus in decay assumption
            arr = ro - y;
            rrr = arr * et;
            y += arr * er;
            ro = y;
            // Check convergence criteria
            for (int iparent=0; iparent < N; iparent++) {
                if (roo(iparent) > 0.0) {
                    if ( abs(1.0 - ro(iparent)/roo(iparent)) > configure::epb) {
                        proceed=true;
                        break;
                    } else {
                        proceed=false;
                    }
                }
            }
            // If number of iteration exceeds maximum constant stop process
            if (t > MAXITERATION) {

                std::cout << "Exceed :number of iteration " << std::endl;
                proceed=false;
            }
            roo = y;
        } // Iteration loop
        // Fill the output dump
        if (configure::outwrite) {
            for (int i = 0; i < y.size(); i++) {
                configure::dumpoutput[k][i][0] = y(i);
            }
        }
    }//numstep
}

//! Iterative method implementation v.2.
void iterative(xt::xarray<double>& matrix, xt::xarray<double>& sigp,
               xt::xarray<double>& y,
               xt::xarray<double>& dmatrix, xt::xarray<double>& dsigp,
               xt::xarray<double>& dy) {
    // Time step
    double dt {configure::timestep/configure::numstep};
    // Array shapes
    std::vector<std::size_t> shape = {y.size()};
    std::vector<std::size_t> dshape = {y.size(), y.size()};
    size_t N {y.size()};
    // Auxilary arrays declaration
    xt::xarray<double> ro = xt::zeros<double>(shape);
    xt::xarray<double> dro = xt::zeros<double>(shape);
    xt::xarray<double> roo = xt::zeros<double>(shape);
    xt::xarray<double> rrr = xt::zeros<double>(shape);
    xt::xarray<double> drr = xt::zeros<double>(shape);
    xt::xarray<double> arr = xt::zeros<double>(shape);
    xt::xarray<double> arrtemp = xt::zeros<double>(dshape);
    xt::xarray<double> drrtemp = xt::zeros<double>(dshape);
    xt::xarray<double> rest = xt::zeros<double>(shape);
    xt::xarray<double> disr = xt::zeros<double>(shape);
    xt::xarray<double> et = xt::ones<double>(shape);
    xt::xarray<double> er = xt::ones<double>(shape);
    xt::xarray<double> es = xt::ones<double>(shape);
    xt::xarray<double> ds = xt::ones<double>(shape);
    xt::xarray<double> dr = xt::ones<double>(shape);
    xt::xarray<double> ddr = xt::ones<double>(shape);
    double aa {0.0};
    double dd {0.0};
    // Implementation
    arr = sigp * dt;
    // Rest part of radioactive nuclide through decay for time dt
    rest = xt::exp(-arr);
    // Disapperead part of radioactive nuclide through decay for time dt
    disr = 1 - rest;
    // Fill auxilary arrays
    for (size_t i=0; i < N; i++) {
        if (disr(i) < 1.e-10) disr(i) = arr(i);
        if (arr(i) > 0.0) er(i) = disr(i) / arr(i);
        for (size_t j=0; j < N; j++) {
            if (sigp(j) > 0) matrix(i,j) = matrix(i,j) / sigp(j);
            if (sigp(j) > 0) dmatrix(i,j) = dmatrix(i,j) / sigp(j);

        }
    }
    // Uncertanty estimation
    et = 1.0 - er;
    es = (er + rest) * dsigp;
    roo = y;
    // Time iteration
    for (int k = 0; k < configure::numstep; k++) {
        bool proceed {true};
        int t = 0;
        rrr = y * disr;       //!< disappear at time dt
        // Uncertainties arrays initialization
        drr = dy * disr;
        ddr = (sigp * dt * dsigp + disr * dy);
        dr = y * rest * ddr;
        dy = dy + sigp * dt * dsigp;
        ds = y * rest * dy;
        dro = dy;
        y = y * rest;         //!< rest at time dt
        ro = y;
        // Main loop
        while (proceed) {
            t++;
            // Transition of nuclear concentration part from parent
            // to child nuclear
            for (size_t iparent=0; iparent < N; iparent++) {
                for (size_t ichild=0; ichild < N; ichild++) {
                    arrtemp(ichild, iparent) = matrix(ichild, iparent) *
                            rrr(iparent);
                    if (t == 0) {
                        drrtemp(ichild, iparent) = dmatrix(ichild, iparent) +
                                dy(iparent);
                    } else {
                        drrtemp(ichild, iparent) = dmatrix(ichild, iparent) +
                                dr(iparent);

                    }
                    drrtemp(ichild, iparent) = arrtemp(ichild, iparent) *
                            drrtemp(ichild, iparent);
                }
            }
            // Additive nucleus through decay
            for (size_t ichild=0; ichild < N; ichild++) {
                aa = xt::sum(xt::row(arrtemp, ichild))(0);
                dd = xt::sum(xt::row(drrtemp, ichild))(0);
                ro(ichild) += aa;
                dr(ichild) = aa * es(ichild);
                ds(ichild) += dd * er(ichild) + dr(ichild);
                dr(ichild) += dd * et(ichild);
            }
            // Get additives to child nucleus in decay assumption
            arr = ro - y;

            rrr = arr * et;

            y += arr * er;
            for (int iparent=0; iparent < N; iparent++) {
                if (y(iparent) < 0.0) y(iparent) = 0.0;
                if (arr(iparent) > 0.0) dr(iparent) = dr(iparent) /
                        arr(iparent);
                if (dr(iparent) < 0.0) dr(iparent) = 0.0;
                if (y(iparent) > 0.0) dy(iparent) = ds(iparent);
                if (y(iparent) < 1.e-80 || dy(iparent) < 0.) dy(iparent) = 0.0;
            }
            ro = y;
            // Check convergence criteria
            for (int iparent=0; iparent < N; iparent++) {
                if (roo(iparent) > 0.0) {
                    if ( abs(1.0 - ro(iparent)/roo(iparent)) > configure::epb) {
                        proceed=true;
                        break;
                    } else {
                        proceed=false;
                    }
                }
            }
            // If number of iteration exceeds maximum constant stop process
            if (t > MAXITERATION) {
                std::cout << "Exceed :number of iteration " << std::endl;
                proceed=false;
            }
            roo = y;
        }// Iteration loop
        // Fill the output dump
        if (configure::outwrite) {
            for (int i = 0; i < y.size(); i++) {
                configure::dumpoutput[k][i][0] = y(i);
                configure::dumpoutput[k][i][1] = dy(i);
            }
        }
    }//numstep
}

//! Chebyshev rational appproximation method (CRAM)
void cram(xt::xarray<double>& matrix, xt::xarray<double>& y,
          xt::xarray<std::complex<double>>& alpha,
          xt::xarray<std::complex<double>>& theta,
          int order, double alpha0) {
    // Time step
    double dt {configure::timestep/configure::numstep};
    // Problem dimensions
    size_t n {y.size()};
    matrix = matrix * dt;
    Eigen::VectorXd ytemp(n);
    Eigen::VectorXcd x(n), zytemp(n);
    Eigen::SparseMatrix<std::complex<double>> A(n,n);
    std::complex<double> c(1, 0);
    // Reserve space in sparse matrix
    A.reserve(n * 10);
    // Solver initialization
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>,
            Eigen::COLAMDOrdering<int>>  solver;
    // Fill A and b;
    for (size_t i=0; i < n; i++) {
        ytemp(i) = y(i);
        zytemp(i) = c * y(i);
        for (size_t j = 0; j < n; j++)
            if (matrix(i, j) != 0.) A.insert(i, j) = c * matrix(i, j);
    }
    // Prepare matrix
    A.makeCompressed();
    // Time iteration
    for (int t=0; t < configure::numstep; t++) {
        // CRAM order iteration
        for (int it = 0; it < order; it++) {
            if (configure::verbose)
                std :: cout << "Number of iterations is "<< it << std::endl;
            for (size_t j=0; j < n; j++) {
                A.coeffRef(j,j) = c * matrix(j, j) - theta(it);
            }
            // Compute the ordering permutation vector from the structural
            // pattern of A
            solver.analyzePattern(A);
            // Compute the numerical facmatrixtorization
            solver.factorize(A);
            //Use the factors to solve the linear system
            zytemp.real() = ytemp;
            x = solver.solve(zytemp);
            x = alpha(it) * x;
            ytemp += 2 * x.real();
        } // CRAM iteration

        ytemp = ytemp * alpha0;
        for (size_t i=0; i < n; i++)
            y(i) = ytemp(i);
        // Fill the dump
        if (configure::outwrite) {
            for (int i = 0; i < y.size(); i++) {
                if (y(i) < 0.0)
                    y(i) = 0.0;
                configure::dumpoutput[t][i][0] = y(i);
            }
        }
    } // time iteration
}

//! Apply filters to result
void apply_filters(const Chain &chainer, const std::string& matname) {
    std::array<double, 4> storeval {0., 0., 0., 0.};
    std::ofstream output(configure::path_output + "outlog.csv",
                         std::ofstream::app);
    if (!output.is_open()) {
        std::cout << "Warning!: " <<
                     " Output file to write reuslts is not open" << std::endl;
        return;
    }
    std::vector<int> xscale;
    std::vector<int> yscale;
    std::vector<int> mainheader;
    std::vector<int> addheader;
    std::vector<std::string> headnames;
    std::vector<std::string> nuclnames;
    xscale.reserve(configure::numstep);
    yscale.reserve(chainer.name_idx.size());
    nuclnames.reserve(chainer.name_idx.size());
    // Applying a time filter to result
    if (timefilter != nullptr) {
        timefilter->apply(configure::timestep / configure::numstep,
                          configure::numstep, xscale);
    } else {
        for (int i = 0; i < configure::numstep; i++)
            xscale.push_back(i);
    }
    // Fill headers and nuclide names array
    std::for_each(chainer.name_idx.begin(),
                  chainer.name_idx.end(),
                  [&nuclnames](std::pair<std::string, size_t> item){
        nuclnames.push_back(item.first);
    });
    std::for_each(configure::header_names.begin(),
                  configure::header_names.end(),
                  [&headnames](std::string item){
        headnames.push_back(item);
    });
    // Search for ordinary string filters
    auto searchexnuclfilter = std::find_if(filters.begin(),
                                           filters.end(),
                                           [] (Filter& f) {
            return f.type == "exnuclide";});
    auto searchnuclfilter = std::find_if(filters.begin(),
                                         filters.end(),
                                         [] (Filter& f) {
            return f.type == "nuclide";});
    auto searchheaderfilter = std::find_if(filters.begin(),
                                           filters.end(),
                                           [] (Filter& f) {
            return f.type == "header";});
    // If nuclide filter is present -> apply it
    if (searchexnuclfilter != filters.end()) {
        searchexnuclfilter->apply(nuclnames, yscale);
    }
    // If header filter is present -> apply it
    if (searchheaderfilter != filters.end()) {
        mainheader.push_back(0);
        searchheaderfilter->apply(headnames, mainheader);
    } else {
        for (int i = 0; i < configure::header_names.size(); i++)
            mainheader.push_back(i);
    }
    // Addition to header if nuclide concentration specifier is present
    if (searchnuclfilter != filters.end()) {
        searchnuclfilter->apply(nuclnames, addheader);
    }
    for (size_t k = 0; k < mainheader.size(); k++) {
        switch (mainheader[k]) {
        case 0:
            output << "dt" << ";";
            break;
        case 1:
            output << "Act, sec-1" << ";";
            break;
        case 2:
            output << "Q, Mev" << ";";
            break;
        case 3:
            output << "dAct, sec-1" << ";";
            break;
        case 4:
            output << "dQ, Mev" << ";";
            break;
        }
    }
    for (auto& k : addheader)
        output << nuclnames[k] << ";";
    output << std::endl;

    output << matname << ";" << std::endl;
    for (auto t : xscale) {
        output << configure :: timestep / configure :: numstep *
                  (t + 1) << ";";
        for (size_t i = 0; i < nuclnames.size(); i++) {
            if (std::find(yscale.begin(), yscale.end(), i) ==
                    yscale.end()) {
                size_t ijk {chainer.name_idx.at(
                                nuclnames[i])};
                udouble uhalflife {nuclides[chainer.name_idx.at(
                                nuclnames[i])]->half_life};
                udouble uenergy {nuclides[ijk]->decay_energy};
                for (size_t k = 0; k < mainheader.size(); k++) {
                    if (mainheader[k] == 1 &&
                            uhalflife.Real() > 0)//Decay-rate
                        storeval[0] +=
                                log(2.0) /
                                uhalflife.Real()
                                * configure::dumpoutput[t][i][0] *
                                1.e+24;
                    if (mainheader[k] == 2 &&
                            uhalflife.Real() > 0) //Decay heat
                        storeval[1] +=
                                log(2.0) /
                                uhalflife.Real() *
                                uenergy.Real() *
                                configure::dumpoutput[t][i][0] *
                                1.e+24 * 1.e-6;
                    if (mainheader[k] == 3 &&
                            uhalflife.Real() > 0) //d(Decay-rate)
                        storeval[2] +=
                                (log(2.0) /
                                 uhalflife).Dev() *
                                configure::dumpoutput[t][i][0] *
                                1.e+24 +
                                log(2.0) /
                                uhalflife.Real()
                                * configure::dumpoutput[t][i][1] *
                                1.e+24;
                    if (mainheader[k] == 4 &&
                            uhalflife.Real() > 0) //d(Decay heat)
                        storeval[3] +=
                                ((log(2.0) /
                                  uhalflife * uenergy
                                  ).Dev() *
                                 configure::dumpoutput[t][i][0] +
                                log(2.0) /
                                uhalflife.Real() *
                                uenergy.Real() *
                                configure::dumpoutput[t][i][1]) *
                                1.e+24 * 1.e-6;
                } // for with if
            }
        } // for nuclides
        for (size_t k = 1; k < mainheader.size(); k++) {
            output << storeval[mainheader[k] - 1] << ";";
            storeval[mainheader[k] - 1] = 0.0;
        }
        for (auto& k : addheader)
            output << configure::dumpoutput[t][k][0] << ";";
        output << std::endl;
    } // for time

    output.close();
}

}//namespace executer

}//namespace openbps
