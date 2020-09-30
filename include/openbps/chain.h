#ifndef CHAIN_H
#define CHAIN_H

#include <iostream>
#include <string.h>
#include <vector>
#include <map>
#include <memory>
#include <tuple>
#include <sstream>
#include "../extern/pugiData/pugixml.h"
#include "parse.h"
#include "nuclide.h"

namespace openbps {

//==============================================================================
// Class to handle with nuclide chain structure
//==============================================================================

class Chain {
public:
    //--------------------------------------------------------------------------
    //! Constructors, destructors, factory functions
    Chain(){}
    ~Chain() = default;
    Chain(pugi::xml_node node);
    //--------------------------------------------------------------------------
    //! Methods
    //! Get the information about chain nuclide names
    //!
    //! \return vector with pairs all nuclides in chain and its number
    std::vector<std::pair<int, std::string>> form_idx_name();
    //! Get the information about chain nuclide decay constant [sec-1]
    //!
    //! \return vector with pair nuclide number in global array-decay constant
    std::vector<std::pair<int, double>> form_idx_lambda();
    //! Get the information about chain nuclide decay branching ratio
    //!
    //! \return vector with triple nuclide and target nuclide number in global
    //!          array and branching ration
    std::vector<std::tuple<int, int, udouble> > form_idx_decay();
    //! Get the information about chain from nuclide to nuclide by neutron
    //! reaction
    //! \return map with key energy name and vector of pair from nuclide to
    //! target
    std::map<std::string, std::vector<std::pair<int, int>>> form_reaction();
    //! Get the information about chain yield fissions
    //!
    //! \return map with key energy value and vector of pair from nuclide to
    //! target
    //!         with fraction
    std::map<double, std::vector<std::tuple<int, int, double>>>
    form_yield_map();
    //! Get the information about chain yield fissions
    //!
    //! \return a pair of vector energies and fractions
    std::pair<std::vector<double>, std::vector<double>>
    get_yield_map(size_t father, const std::string& daughter);
    //! Get an index in chain nuclides array
    //!
    //! \param[in] name of nuclide
    //! \return index in chain order
    int get_nuclide_index(const std::string& name);
    //--------------------------------------------------------------------------
    //! Attributes
    std::map<std::string, size_t> name_idx; //!< Dictionary between \keyword
                                            //!< nuclide name and an index
                                            //!< in nuclides vector;
};

//==============================================================================
// Non class methods
//==============================================================================
//! Get a chain from xml file
//!
//! \param[in] filepath xml file containing chain
Chain read_chain_xml(const std::string& filename);
/*xt::xarray<double> form_matrix(Chain& chainer, Materials& mat);
    xt::xarray<double> form_sigp(Chain& chainer, Materials& mat);
    xt::xarray<double> form_dmatrix(Chain& chainer, Materials& mat);
    xt::xarray<double> form_dsigp(Chain& chainer, Materials& mat);
    xt::xarray<double> make_concentration(Chain& chainer, std::vector<std::string>& nameconc,
                                          std::vector<double>& ro);*/

} // namespace openbps

#endif // CHAIN_H
