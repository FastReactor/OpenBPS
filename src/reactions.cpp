#include "openbps/reactions.h"
#include <vector>
#include <map>
#include <list>
#include <iostream>
#include <memory>
#include <string.h>
#include "openbps/configure.h"
#include "openbps/functionals.h"
#include "openbps/parse.h"
#include "../extern/pugiData/pugixml.h"

namespace openbps {

//==============================================================================
// Global variables
//==============================================================================

std::vector<std::unique_ptr<Composition>> compositions;
std::map<std::string, int> composmap;
int indexall {-1};
std::vector<Xslibs> externxslibs;

//==============================================================================
// Class cross-section implementations
//==============================================================================

Sxs::Sxs(pugi::xml_node node,const std::string& rxs, const std::string& redex) {
    this->xstype = node.attribute("reaction").value();
    this->xsname = node.attribute("name").value();
    if (rxs == "cs") {
        this->xs_ = get_node_array<udouble>(node, redex.c_str());
    } else {
        this->rxs = get_node_array<udouble>(node, redex.c_str());
    }
}

Xslibs::Xslibs(pugi::xml_node node) {
    // Get an energy multigroup structure
    if (check_for_node(node, "energy")) {
        for (pugi::xml_node tool : node.children("energy"))
            this->energies_=get_node_array<double>(tool,
                                      "energy");
        this->numgroup = this->energies_.size();
    }
    if (check_for_node(node, "xslibs")) {
        parse_xml_xslibs_(node, this->xsdata);
    }
}

//==============================================================================
// Composition class implementation
//==============================================================================

Composition::Composition(pugi::xml_node node) {
    // Composition name
    if (check_for_node(node, "name")) {
        this->name_ = node.attribute("name").value();
    }
    // Read energy discretization
    if (check_for_node(node, "energy")) {
        for (pugi::xml_node tool : node.children("energy")) {
            size_t cng = std::stoi(tool.attribute("ng").value());
            this->energies_.insert({cng, get_node_array<double>(tool,
                                    "energy")});
            this->energy_number_++;
        }
    }
    // Read a flux
    if (check_for_node(node, "flux")) {
        for (pugi::xml_node tool : node.children("flux")) {
            this->flux_= get_node_array<udouble>(tool, "flux");
        }
    }
    if (check_for_node(node, "dflux")) {
        for (pugi::xml_node tool : node.children("dflux")) {
            std::vector<double> dflux {get_node_array<double>(tool, "dflux")};
            for (size_t j = 0; j < dflux.size() && j < flux_.size(); j++)
                    flux_[j].Adddeviation(dflux[j]);

        }
    }
    // Read a spectrum (if a flux not presented)
    if (check_for_node(node, "spectrum")) {
        for (pugi::xml_node tool : node.children("spectrum")) {
            size_t cng = std::stoi(tool.attribute("ng").value());
            this->spectrum_= get_node_array<udouble>(tool, "spectrum");
        }
    }

    if (check_for_node(node, "dspectrum")) {
        for (pugi::xml_node tool : node.children("dspectrum")) {
            std::vector<double> dspectrum {get_node_array<double>(tool,
                                                                  "dspectrum")};
            for (size_t j = 0; j < dspectrum.size() && j< spectrum_.size(); j++)
                    spectrum_[j].Adddeviation(dspectrum[j]);

        }
    }
    // Read a cross-section data
    if (check_for_node(node, "xslibs"))
        parse_xml_xslibs_(node, this->xslib);

    this->spectrum_.size() > 0 ?
            this->energy_number_ = this->spectrum_.size() :
            this->energy_number_ = this->flux_.size();
}

//! Auxilary function to copy data from xslib
Sxs parse_xs_xml_
(pugi::xml_node node,const std::string& rxs, const std::string& redex) {
    Sxs result(node, rxs, redex);
    return result;
}

//! Parse xslibs
void parse_xml_xslibs_(pugi::xml_node node, std::vector<Sxs>& xssource) {
    std::string rxs {node.child("xslibs").attribute("typex").value()};
    for (pugi::xml_node tool : node.child("xslibs").children("xslib")) {
        xssource.push_back(parse_xs_xml_(tool, rxs, "xslib"));
    }
    for (pugi::xml_node tool : node.child("xslibs").children("dxslib")) {
        // Read deviation part of cross section
        Sxs deriv_rxs(tool, rxs, "dxslib");
        for (auto ixs = xssource.begin();
             ixs != xssource.end(); ixs++) {
            if (ixs->xsname == deriv_rxs.xsname &&
                    ixs->xstype == deriv_rxs.xstype) {
                for (size_t j = 0; j < deriv_rxs.rxs.size() &&
                     j< ixs->rxs.size(); j++)
                        ixs->rxs[j].Adddeviation(deriv_rxs.rxs[j]);
                for (size_t j = 0; j < deriv_rxs.xs_.size() &&
                     j< ixs->xs_.size(); j++)
                        ixs->xs_[j].Adddeviation(deriv_rxs.xs_[j]);
            }
        }
    }
}

//! Auxilary function to copy data from xslib
void Composition::depcopymap_(std::map<size_t, std::vector<double>>& fmap,
                              std::map<size_t, std::vector<double>>& smap) {

    std::map<size_t, std::vector<double>>::iterator it;
    if (!smap.empty()) {
        for (it = smap.begin(); it != smap.end(); ++it) {
            if (fmap.find(it->first) == fmap.end()){
                fmap[it->first] = it->second;
            }
        }
    }

}

//! Copy data from composition marked name "all" in xml file
void Composition::deploy_all(Composition& externcompos) {

    if (externcompos.name_ == this->name_) {
        return;
    } else {

        this->depcopymap_(this->energies_, externcompos.energies_);
        if (!externcompos.spectrum_.empty() && this->spectrum_.empty()) {
            this->spectrum_.resize(externcompos.spectrum_.size());
            std::copy(externcompos.spectrum_.begin(),
                      externcompos.spectrum_.end(),
                      this->spectrum_.begin());
        }
        if (!externcompos.flux_.empty() && this->flux_.empty()) {
            this->flux_.resize(externcompos.flux_.size());
            std::copy(externcompos.flux_.begin(),externcompos.flux_.end(),
                      this->flux_.begin());
        }
        if (!externcompos.xslib.empty() && this->xslib.empty()){
            this->xslib.resize(externcompos.xslib.size());
            std::copy(externcompos.xslib.begin(),externcompos.xslib.end(),
                      this->xslib.begin());

        }
    }
}

//! Calculate reaction-rate from cross-section data and flux
void Composition::calculate_rr_(Sxs& ixs,const std::vector<double>& extenergy,
                   udouble& rxs) {
    size_t ng = ixs.xs_.size();
    std::vector<udouble> cuflux(ng, 1.0);
    if (this->flux_.size() > 0) {
        if (this->flux_.size() == ng) {
            cuflux = this->flux_;
        } else {
            cuflux = collapsing<udouble>(this->energies_[
                                         this->flux_.size()],
                                         this->flux_,
                                         extenergy);
        }
    }
    for (int i = 0; i != ng; i++)
        rxs = rxs + cuflux[i] * ixs.xs_[i];
}

//! Calculate reaction rate for all reactions in xslib
void  Composition::get_reaction() {
    for (auto ixs = this->xslib.begin(); ixs != this->xslib.end(); ixs++) {
        if (ixs->rxs.empty()) {
            ixs->rxs.resize(1);
            ixs->rxs[0] = 0.0;
            if (!ixs->xs_.empty()) {
                calculate_rr_(*ixs, this->energies_[ixs->xs_.size()],
                        ixs->rxs[0]);
            } // if xs_ not empty
        } // if rxs is empty
    } //for
}

//! Import data from external cross section source
void Composition::import_xsdata(Xslibs& implibs) {
    for (auto& xs : implibs.xsdata) {
        if (std::find_if(xslib.begin(), xslib.end(), [&xs](Sxs& lib) {
                         return xs.xsname == lib.xsname &&
                         xs.xstype == lib.xstype;
        }) == xslib.end()) {
            xslib.push_back(Sxs());
            xslib[xslib.size() - 1].xsname = xs.xsname;
            xslib[xslib.size() - 1].xstype = xs.xstype;
            xslib[xslib.size() - 1].rxs.resize(1);
            calculate_rr_(xs, implibs.get_egroups(),
                          xslib[xslib.size() - 1].rxs[0]);
        }
    }
}

//! Comparator function
bool compfe (std::pair <double, double> a,std::pair <double, double> b) {
    return a.first < b.first;
}

std::pair<std::vector<double>, std::vector<double>>
Composition::get_fluxenergy() {
    std::pair<std::vector<double>, std::vector<double>> result;
    size_t ng;
    double fluxnorm {0.0};
    if (spectrum_.empty()) {
        std::for_each(flux_.begin(), flux_.end(), [&] (udouble n) {
            fluxnorm += n.Real();
        });
        for (auto& f: flux_) {
            spectrum_.push_back(f / fluxnorm);
        }
    }
    ng = spectrum_.size();

    result = std::make_pair(energies_[ng], usplit<double>(spectrum_).first);
    return result;
}

//==============================================================================
// Non class methods implementation
//==============================================================================
//! Read compositions from xml file
void read_reactions_xml() {
    using namespace configure;
    pugi::xml_document doc;
    auto result = doc.load_file(reaction_file.c_str());
    if (!result) {
        std::cout << "Warning: file reactions.xml not found!" << std::endl;
        // If reactions.xml not found run in decay only mode
        configure::decay_extra_out = true;
        return;
    }
    pugi::xml_node root_node = doc.child("compositions");
    if (configure::verbose)
        std::cout << "I' m in reactions.xml parser" << std::endl;

    for (pugi::xml_node tool : root_node.children("composit")) {
        compositions.push_back(std::unique_ptr<Composition>(new
                                                            Composition(tool)));
        int index = compositions.size() - 1;
        if (compositions[index]->Name() == "all") {
            indexall = index;
        }
        composmap.insert({compositions[index]->Name(), index});
    }

    // Read imported cross-section libraries
    if (libs.size() > 0)
        read_importedlib_xml();

    for (size_t i = 0; i < compositions.size(); i++)
        if (i != indexall) {
            if (indexall > -1)
                compositions[i]->deploy_all(*compositions[indexall]);
            compositions[i]->get_reaction();
            if (externxslibs.size() > 0)
                for (auto& v: externxslibs)
                    compositions[i]->import_xsdata(v);
    }
}

//! Read an imported cross section data library from *.xml files
void read_importedlib_xml() {
    for (auto& v : configure::libs) {
        pugi::xml_document doc;
        auto result = doc.load_file(v.c_str());
        if (!result) {
            std::cout << "Warning: file " << v << " not found!" << std::endl;
            return;
        }
        pugi::xml_node root_node = doc.child("importlib");
        externxslibs.emplace_back(root_node);
    }
}

} //namespace openbps
