#include "openbps/matrix.h"
#include <vector>
#include <memory>
#include <initializer_list>
#include <cmath>
#include "xtensor/xarray.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xview.hpp"
#include "openbps/uncertainty.h"
#include "openbps/functionals.h"
#include "openbps/nuclide.h"
#include "openbps/chain.h"
#include "openbps/reactions.h"
#include "openbps/materials.h"
#include "openbps/configure.h"

namespace openbps {

//==============================================================================
// Decay matrix implementation
//==============================================================================

//! Form a decay nuclide matrix
void DecayMatrix::form_matrixreal(Chain& chain) {
    size_t i {0};
    size_t inucl {0};
    size_t k;
    udouble decay_ {0.0, 0.0};
    //Run over all nuclides
    for (auto it = chain.name_idx.begin(); it != chain.name_idx.end(); it++) {
        inucl = it->second; // from nuclides
        i = chain.get_nuclide_index(it->first);
        if (nuclides[inucl]->half_life.Real() > 0) {
            decay_ = (log(2.0) / nuclides[inucl]->half_life);
            this->data_[i][i] = -decay_.Real();
            // Cycle over all decay modes
            for (auto& v : nuclides[inucl]->get_decaybr()) {
                if (v.first != std::string("Nothing")) {
                    // To nuclide
                    k = chain.get_nuclide_index(v.first);
                    this->data_[k][i] +=
                            decay_.Real() * v.second.Real();
                }
            }

        }
    }

}

//! Form a deviation decay nuclide matrix
void DecayMatrix::form_matrixdev(Chain& chain) {
    size_t i {0};
    size_t inucl {0};
    size_t k;
    udouble decay_ {0.0, 0.0};
    //Run over all nuclides
    for (auto it = chain.name_idx.begin(); it != chain.name_idx.end(); it++) {
        inucl = it->second; // from nuclide
        i = chain.get_nuclide_index(it->first);
        if (nuclides[inucl]->half_life.Real() > 0) {
            decay_ = (log(2.0) / nuclides[inucl]->half_life);
            this->data_[i][i] = -decay_.Dev();
            // Cycle over all decay modes
            for (auto& v : nuclides[inucl]->get_decaybr()) {
                if (v.first != std::string("Nothing")) {
                    // To nuclide
                    k = chain.get_nuclide_index(v.first);
                    this->data_[k][i] +=
                            (decay_ * v.second).Dev();
                }
            }

        }
    }
}

//==============================================================================
// IterativeMatrix implementation
//==============================================================================

//! Methods
//! Form a real decay nuclide matrix
xt::xarray<double> IterMatrix::matrixreal(Chain& chain,
                                            const Materials& mat) {    
    // Variables for calculation fissiop product yields by spectrum
    std::pair<std::vector<double>, std::vector<double>> pair2;
    std::vector<double> weight;
    size_t k {0};
    //
    size_t NN {chain.name_idx.size()}; //!< Nuclide number
    std::vector<std::size_t> shape = { NN, NN };
    xt::xarray<double> result(shape, 0.0);
    for (size_t i = 0; i < NN; i++) {
        std::copy(&this->data_[i][0], &this->data_[i][NN],
                  (result.begin() + i * NN));
        result(i, i) = 0.0;
    }
    int icompos {mat.numcomposition};
    if (icompos > -1) {
        // Get neutron flux - energy value descretization
        pair2 = compositions[icompos]->get_fluxenergy();
        for (auto it = chain.name_idx.begin();
             it != chain.name_idx.end(); it++) {
            size_t inucl = it->second; // from nuclide
            size_t i = chain.get_nuclide_index(it->first);
            // For xslib
            for (auto& obj : compositions[icompos]->xslib) {
            // If cross section presented for nuclides
            if (obj.xsname == it->first) {
                double rr {0.0};
                // Sum reaction-rates over energy group
                for (auto& r: obj.rxs)
                    rr += r.Real();
                // Iterate over chain nuclides transition
                for (auto& r : nuclides[inucl]->get_reactions()) {
                    // If match
                    if (r.first == obj.xstype) {
                        if (obj.xstype != "fission") {
                            // And non-fission then add
                            size_t k = chain.get_nuclide_index(r.second);
                            result(k, i) += rr * PWD * mat.normpower;
                        } else {
                            // Considering fission reaction by energy separately
                            std::vector<double> energies =
                                    nuclides[inucl]->get_nfy_energies();
                            // Fission yields by product
                            for (auto& item: nuclides[inucl]->
                                     get_yield_product()) {
                                k = chain.get_nuclide_index(item.first);
                                double br {0.0};
                                double norm {0.0};
                                if (weight.empty())
                                    weight = transition(pair2.first,
                                                        pair2.second,
                                                        energies);
                                    for (int l = 0; l < weight.size(); l++) {
                                         br += weight[l] * item.second[l];
                                         norm += weight[l];
                                    } // for weight

                                    result(k, i) += br / norm * rr * PWD *
                                                   mat.normpower;
                                    norm = 1.0;
                            } // for product
                            weight.clear();
                        } // fission yields
                    } // if reaction = chain.reaction
                } // run over reaction
            } // if nuclide is in crossection data and chain
        } // for xslib in composition
    } // if composition is presented
    } // for all nuclides

    return result;
}

//! Form a deviation decay nuclide matrix for unceratanties analysis
xt::xarray<double> IterMatrix::matrixdev(Chain& chain,
                                           const Materials& mat) {
    // Variables for calculation fissiop product yields by spectrum
    std::pair<std::vector<double>, std::vector<double>> pair2;
    std::vector<double> weight;
    size_t k {0};
    //
    size_t NN {chain.name_idx.size()}; //!< Nuclide number
    std::vector<std::size_t> shape = { NN, NN };
    xt::xarray<double> result(shape, 0.0);
    for (size_t i = 0; i < NN; i++) {
        std::copy(&this->data_[i][0], &this->data_[i][NN],
                  (result.begin() + i * NN));
        result(i, i) = 0.0;
    }
    int icompos {mat.numcomposition};
    if (icompos > -1) {
        // Get neutron flux - energy value descretization
        pair2 = compositions[icompos]->get_fluxenergy();
        for (auto it = chain.name_idx.begin();
             it != chain.name_idx.end(); it++) {
            size_t inucl = it->second; // from nuclide
            size_t i = chain.get_nuclide_index(it->first);
            // For xslib
            for (auto& obj : compositions[icompos]->xslib) {
            // If cross section presented for nuclides
            if (obj.xsname == it->first) {
                double rr {0.0};
                // Sum reaction-rates over energy group
                for (auto& r: obj.rxs)
                    rr += r.Dev();
                // Iterate over chain nuclides transition
                for (auto& r : nuclides[inucl]->get_reactions()) {
                    // If match
                    if (r.first == obj.xstype) {
                        if (obj.xstype != "fission") {
                            // And non-fission then add
                            size_t k = chain.get_nuclide_index(r.second);
                            result(k, i) += rr * PWD * mat.normpower;
                        } else {
                            // Considering fission reaction by energy separately
                            std::vector<double> energies =
                                    nuclides[inucl]->get_nfy_energies();
                            // Fission yields by product
                            for (auto& item: nuclides[inucl]->
                                     get_yield_product()) {
                                k = chain.get_nuclide_index(item.first);
                                double br {0.0};
                                double norm {0.0};
                                if (weight.empty())
                                    weight = transition(pair2.first,
                                                        pair2.second,
                                                        energies);
                                    for (int l = 0; l < weight.size(); l++) {
                                         br += weight[l] * item.second[l];
                                         norm += weight[l];
                                    } // for weight

                                    result(k, i) += br / norm * rr * PWD *
                                                   mat.normpower;
                                    norm = 1.0;
                            } // for product
                            weight.clear();
                        } // fission yields
                    } // if reaction = chain.reaction
                } // run over reaction
            } // if nuclide is in crossection data and chain
        } // for xslib in composition
    } // if composition is presented
    } // for all nuclides

    return result;
}

//! Form a dev decay nuclide vector for all transition from every nuclide
xt::xarray<double> IterMatrix::sigp(Chain& chain,
                                      const Materials& mat) {
    std::vector<size_t> shape {chain.name_idx.size()};
    xt::xarray<double> result(shape, 0.0);
    udouble decay_ {0.0, 0.0};
    //Run over all nuclides
    if (configure::verbose)
        std::cout << "SIGP: "<<std::endl;
    for (auto it = chain.name_idx.begin(); it != chain.name_idx.end(); it++) {
        size_t inucl = it->second; // from nuclide

        size_t i = chain.get_nuclide_index(it->first);
        if (nuclides[inucl]->half_life.Real() > 0) {
            decay_ = (log(2.0) / nuclides[inucl]->half_life);
            result(i) = decay_.Real();
        }
        int icompos {mat.numcomposition};
        if (icompos > -1) {
            // For xslib
            for (auto& obj : compositions[icompos]->xslib) {
                // If cross section presented for nuclides
                if (obj.xsname == it->first) {
                    double rr {0.0};
                    // Sum reaction-rates over energy group
                    for (auto& r: obj.rxs)
                        rr += r.Real();
                    result(i) += rr * PWD * mat.normpower;
                }
           }
        }
        if (configure::verbose)
            std::cout << "sigp: " << it->first << " " << result(i)<<std::endl;
    }

    return result;
}

//! Form a dev decay nuclide vector for all transition from every nuclide
xt::xarray<double> IterMatrix::dsigp(Chain& chain,
                                     const Materials& mat) {
    std::vector<size_t> shape {chain.name_idx.size()};
    xt::xarray<double> result(shape, 0.0);
    udouble decay_ {0.0, 0.0};
    //Run over all nuclides
    for (auto it = chain.name_idx.begin(); it != chain.name_idx.end(); it++) {
        size_t inucl = it->second; // from nuclide

        size_t i = chain.get_nuclide_index(it->first);
        if (nuclides[inucl]->half_life.Real() > 0) {
            decay_ = (log(2.0) / nuclides[inucl]->half_life);
            result(i) = decay_.Dev();
        }
        int icompos {mat.numcomposition};
        if (icompos > -1) {
            // For xslib
            for (auto& obj : compositions[icompos]->xslib) {
                // If cross section presented for nuclides
                if (obj.xsname == it->first) {
                    double rr {0.0};
                    // Sum reaction-rates over energy group
                    for (auto& r: obj.rxs)
                        rr += r.Dev();
                    result(i) += rr * PWD * mat.normpower;
                }
           }
        }
    }
    return result;
}

//==============================================================================
// ChebyshevMatrix implementation
//==============================================================================

//! Form a real decay nuclide matrix
xt::xarray<double> CramMatrix::matrixreal(Chain& chain,
                                            const Materials& mat) {
    // Variables for calculation fissiop product yields by spectrum
    std::pair<std::vector<double>, std::vector<double>> pair2;
    std::vector<double> weight;
    size_t k {0};
    //
    size_t NN {chain.name_idx.size()}; //!< Nuclide number
    std::vector<std::size_t> shape = { NN, NN };
    xt::xarray<double> result(shape, 0.0);
    for (size_t i = 0; i < NN; i++)
        std::copy(&this->data_[i][0], &this->data_[i][NN],
                  (result.begin() + i * NN));
    int icompos {mat.numcomposition};
    if (icompos > -1) {
        // Get neutron flux - energy value descretization
        pair2 = compositions[icompos]->get_fluxenergy();
        for (auto it = chain.name_idx.begin();
             it != chain.name_idx.end(); it++) {
            size_t inucl = it->second; // from nuclide
            size_t i = chain.get_nuclide_index(it->first);
            // For xslib
            for (auto& obj : compositions[icompos]->xslib) {
            // If cross section presented for nuclides
            if (obj.xsname == it->first) {
                double rr {0.0};
                // Sum reaction-rates over energy group
                for (auto& r: obj.rxs)
                    rr += r.Real();
                // Iterate over chain nuclides transition
                for (auto& r : nuclides[inucl]->get_reactions()) {
                    // If match
                    if (r.first == obj.xstype) {
                        if (obj.xstype != "fission") {
                            // And non-fission then add
                            size_t k = chain.get_nuclide_index(r.second);
                            result(k, i) += rr * PWD * mat.normpower;
                        } else {
                            // Considering fission reaction by energy separately
                            std::vector<double> energies =
                                    nuclides[inucl]->get_nfy_energies();
                            // Fission yields by product
                            for (auto& item: nuclides[inucl]->
                                     get_yield_product()) {
                                k = chain.get_nuclide_index(item.first);
                                double br {0.0};
                                double norm {0.0};
                                if (weight.empty())
                                    weight = transition(pair2.first,
                                                        pair2.second,
                                                        energies);
                                    for (int l = 0; l < weight.size(); l++) {
                                         br += weight[l] * item.second[l];
                                         norm += weight[l];
                                    } // for weight

                                    result(k, i) += br / norm * rr * PWD *
                                                   mat.normpower;
                                    norm = 1.0;
                            } // for product
                            weight.clear();
                        } // fission yields
                    } // if reaction = chain.reaction
                } // run over reaction
                result(i, i) -= rr * PWD * mat.normpower;
            } // if nuclide is in crossection data and chain
        } // for xslib in composition
    } // if composition is presented
    } // for all nuclides

    return result;
}

//==============================================================================
// Non class methods
//==============================================================================

//! Get a concentration vector according to chain nuclide list
xt::xarray<double> make_concentration(Chain& chainer, 
                                      const std::vector<std::string>& nameconc,
                                      const std::vector<udouble>& ro,
                                      bool isDev) {
    xt::xarray<size_t>::shape_type shape = {chainer.name_idx.size()};
    xt::xarray<double> result = xt::zeros<double>(shape);
    for (size_t i = 0; i < ro.size(); i++) {
        if (std::find_if(chainer.name_idx.begin(), chainer.name_idx.end(),
                         [&] (std::pair<std::string, size_t> item){
                         return nameconc[i] == item.first;}) !=
                         chainer.name_idx.end()) {
            if (isDev) {
                result[chainer.get_nuclide_index(nameconc[i])] = ro[i].Dev();
            } else {
                result[chainer.get_nuclide_index(nameconc[i])] = ro[i].Real(); 
            }
        }
    }
    return result;
}

//! Find out power normalization coefficient to reaction-rate
void power_normalization(Materials& mat) {
     double R {0.0};
     for (size_t i = 0; i < mat.conc.size(); i++) {
         auto search = std::find_if(nuclides.begin(),
                                    nuclides.end(),
                                    [&mat, 
                                     i](std::unique_ptr<ChainNuclide>& item) {
                                    return item->name_ == mat.namenuclides[i];});
         if (search != nuclides.end()) {
             if (mat.numcomposition > -1) {
                 for (auto& obj : compositions[mat.numcomposition]->xslib) {
                      if (obj.xsname == mat.namenuclides[i]) {
                           size_t index = std::distance(nuclides.begin(), search);
                           auto releases = nuclides[index]->get_qvalue();
                           if (releases.find(obj.xstype) != releases.end()) {
                               double rr {0.0};
                               // Sum reaction-rates over energy group
                               for (auto& r: obj.rxs)
                                   rr += r.Real();
                               R += rr * releases[obj.xstype] *
                                       mat.conc[i].Real();//eV or Mev
                           }
                      }
                 }
             }
         }
     }
     if (R > 0 && mat.Volume() > 0.0) {
         mat.normpower = mat.Power() * PWRC / (R * mat.Volume());
     }
}

} //namespace openbps
