#include "openbps/nuclide.h"

#include <iostream>
#include <string.h>
#include <vector>
#include <map>
#include <sstream>
#include <memory>
#include "../extern/pugiData/pugixml.h"
#include "openbps/configure.h"

namespace openbps {

//==============================================================================
// Global variables
//==============================================================================

bool isNuclidePresent {false};
// Model nuclides container
std::vector<std::unique_ptr<ChainNuclide>> nuclides;

//==============================================================================
// Nuclide class implementation
//==============================================================================
ChainNuclide::s_decay ChainNuclide::parse_decay_(pugi::xml_node node) {
    ChainNuclide::s_decay temp;

    temp.type = node.attribute("type").value();
    temp.target = node.attribute("target").value();
    temp.branching_ratio = udouble(atof(node.
                                        attribute("branching_ratio").value()),
                                   atof(node.
                                        attribute("d_branching_ratio").value()));
    temp.q = udouble(atof(node.
                          attribute("q").value()),
                     atof(node.
                          attribute("dq").value()));
    return temp;
}

ChainNuclide::s_reaction ChainNuclide::parse_reaction_(pugi::xml_node node) {
    ChainNuclide::s_reaction temp;

    temp.type = node.attribute("type").value();
    temp.q = atof(node.attribute("Q").value());
    temp.target = node.attribute("target").value();
    return temp;
}

ChainNuclide::s_yield ChainNuclide::parse_yield_(pugi::xml_node node) {
    ChainNuclide::s_yield temp;
    std::vector<std::string> lnuclides {
        get_node_array<std::string>(node, "products")};
    std::vector<double> numbers {
        get_node_array<double>(node, "data")};

    temp.energy = atof(node.attribute("energy").value());
    for (std::size_t i = 0; i < lnuclides.size(); i++) {
        temp.product_data.insert({lnuclides[i], numbers[i]});
    }
    return temp;
}

ChainNuclide::s_nfy ChainNuclide::parse_nfy_(pugi::xml_node node) {
    ChainNuclide::s_nfy temp;

    temp.energies = get_node_array<double>(node, "energies");
    for (pugi::xml_node tool : node.children("fission_yields")) {
        temp.yield_arr.push_back(parse_yield_(tool));
    }
    return temp;
}

//! Fill the yieldproduct map
void ChainNuclide::make_product_yields_() {

    for (size_t i = 0; i < this->nfy.energies.size(); i++) {
        for (auto& items : this->nfy.yield_arr[i].product_data)
            if (this->nfy.yieldproduct.find(items.first) ==
                    this->nfy.yieldproduct.end()) {
                this->nfy.yieldproduct[items.first] = std::vector<double>();
                this->nfy.yieldproduct[items.first].push_back(items.second);
            } else {
                this->nfy.yieldproduct[items.first].push_back(items.second);
            }
    }

}

// Read info from chain.xml file
void ChainNuclide::read_from_chain_xml(pugi::xml_node nuclide_node) {

    this->name_ = nuclide_node.attribute("name").value();
    this->decay_modes = atoi(nuclide_node.attribute("decay_modes").value());
    this->reactions = atoi(nuclide_node.attribute("reactions").value());
    this->half_life = udouble(atof(nuclide_node.
                                   attribute("half_life").value()),
                              atof(nuclide_node.
                                   attribute("d_half_life").value()));
    // Decay modes procssing
    if (this->decay_modes > 0) {
        this->decay_energy = udouble(0.0, 0.0);
                //?atof(nuclide_node.attribute("decay_energy").value());
        for (pugi::xml_node tool : nuclide_node.children("decay")) {
            this->decay_arr.push_back(parse_decay_(tool));
        }
        this->decay_energy =
                atof(nuclide_node.attribute("decay_energy").value());
    } else {
        this->decay_energy = 0;
    }

    // Reactions processing
    if (this->reactions > 0) {
        for (pugi::xml_node tool : nuclide_node.children("reaction")) {
            this->reaction_arr.push_back(parse_reaction_(tool));
        }
    }

    // Neutron fission yields processing
    if (nuclide_node.child("neutron_fission_yields")) {
        this->nfy = parse_nfy_(nuclide_node.child("neutron_fission_yields"));
        this->make_product_yields_();
    }
    if (configure::verbose)
        std::cout << "Parse finish for " << this->name_ << std::endl;
}

ChainNuclide::ChainNuclide(pugi::xml_node nuclide_node) {

    this->name_ = nuclide_node.attribute("name").value();
    this->Z_ = atoi(nuclide_node.attribute("z").value());
    this->A_ = atoi(nuclide_node.attribute("a").value());
    this->m_ = atoi(nuclide_node.attribute("m").value());
    this->awr_ = atof(nuclide_node.attribute("awr").value());
}

//! Get a decay transition information
std::vector<std::pair<std::string,udouble>> ChainNuclide::get_decaybr() {
    std::vector<std::pair<std::string,udouble>> result;
    for (int i = 0; i < this->decay_arr.size(); i++) {
        result.push_back(std::make_pair(this->decay_arr[i].target,
                                        this->decay_arr[i].branching_ratio));
    }
    return result;
}

//! Get a reactions for nuclide
std::vector<std::pair<std::string, std::string>> ChainNuclide::get_reactions() {
    std::vector<std::pair<std::string, std::string>> result;
    for (int i = 0; i < this->reaction_arr.size(); i++) {
        result.push_back(std::make_pair(this->reaction_arr[i].type,
                         this->reaction_arr[i].target));
    }
    return result;
}

//! Get an energy release per reaction for nuclide
std::map<std::string, double> ChainNuclide::get_qvalue() {
    std::map<std::string, double> result;
    for (int i = 0; i < this->reaction_arr.size(); i++) 
        if (this->reaction_arr[i].q > 0)
            result[this->reaction_arr[i].type] = this->reaction_arr[i].q;
    return result;

}

//! Get a neutron fission energies
std::vector<double> ChainNuclide::get_nfy_energies() {
    std::sort(this->nfy.energies.begin(), this->nfy.energies.end(),
              [](double& felem, double& selem) {return felem < selem;});
    return this->nfy.energies;
}

//! Get a neutron fission data by energy index
std::map<std::string, double> ChainNuclide::get_product_data(int index) {
    return this->nfy.yield_arr[index].product_data;
}

//! Get a fission yields for group energy discretization
std::map<std::string,
std::vector<double>> ChainNuclide::get_yield_product() {

    return this->nfy.yieldproduct;
}

//==============================================================================
// Non class methods
//==============================================================================
//! Get a nuclides from xml file
void read_nuclide_xml(const std::string& filepath) {
    // Parse configure.xml file
    pugi::xml_document doc;
    auto result = doc.load_file(filepath.c_str());
    if (!result) {
        std::cerr << "Error while processing nuclide.xml file" ;
        return;
    }
    // Get root element
    pugi::xml_node root = doc.document_element();
    // Proceed all nuclides
    for (pugi::xml_node tool : root.children("Nuclide")) {
        nuclides.push_back(std::unique_ptr<ChainNuclide>(new ChainNuclide(tool))); //!TODO make_unique
    }
    if (nuclides.size() > 0) isNuclidePresent = true;

}

//! Get an index in nucludes array by name
int get_nuclidarray_index(const std::string& name) {
    int result {-1};
    auto it = std::find_if(nuclides.begin(), nuclides.end(),
                           [name](std::unique_ptr<ChainNuclide>& v) {
        return v->name_ == name;});
    if (it != nuclides.end())
        result = std::distance(nuclides.begin(), it);

    return result;
}

} //namespace openbps

