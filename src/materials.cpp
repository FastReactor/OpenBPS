#include "openbps/materials.h"
#include "openbps/configure.h"
#include "openbps/parse.h"
#include "openbps/nuclide.h"
#include <memory>
#include <algorithm>
#include "../extern/pugiData/pugixml.h"
#include "openbps/reactions.h"

namespace openbps {

std::vector<std::unique_ptr<Materials>> materials; //!< vector consisting
                                                   //!< all materials

//==============================================================================
// Materials class implementation
//==============================================================================
Materials::Materials(){}

void Materials::xml_add_material(pugi::xml_node node) {
    auto material = node.append_child("material");
    material.append_attribute("name") =  this->Name().c_str();
    material.append_attribute("volume") =this->Volume();
    material.append_attribute("mass") =  this->Mass();
    material.append_attribute("power") = this->Power();

    auto nameofn = material.append_child("namenuclides");
    nameofn.append_child(pugi::node_pcdata).set_value(join(namenuclides,
                                                           " ").c_str());
    auto concpair = usplit(this->conc);
    auto concnod = material.append_child("conc");
    concnod.append_child(pugi::node_pcdata).set_value(joinDouble(concpair.first,
                                                                 " ").c_str());
    auto dconcnod = material.append_child("dconc");
    dconcnod.append_child(pugi::node_pcdata).set_value(joinDouble(concpair.second,
                                                                  " ").c_str());
}

void Materials::add_nuclide(const std::string& extname, udouble extconc) {
    auto it = std::find(this->namenuclides.begin(),
                        this->namenuclides.end(), extname);
    if (it == this->namenuclides.end()) {
        this->namenuclides.push_back(extname);
        this->conc.push_back(extconc);
        this->indexnuclides.push_back(get_nuclidarray_index(extname));
    } else {
        auto index = std::distance(this->namenuclides.begin(), it);
        this->conc[index] = extconc;
    }
}

/*void Materials::bindcomposition(Composition& extcompos){
  this->compos = std::make_shared<Composition>(extcompos);
}*/

//==============================================================================
// Non class method implementation
//==============================================================================

//! Write down materials vector information into xml file
void form_materials_xml(std::string xml_path) {
    pugi::xml_document doc;
    auto declarationNode = doc.append_child(pugi::node_declaration);
    declarationNode.append_attribute("version") = "1.0";
    declarationNode.append_attribute("encoding") = "UTF-8";

    auto materialsnod = doc.append_child("materials");
    for(int i = 0; i < materials.size(); i++) {
        materials[i]->xml_add_material(materialsnod);
        if (configure::verbose)
            std :: cout << "material was included" << std ::endl;
    }
    bool saveSucceeded = doc.save_file(xml_path.c_str(), PUGIXML_TEXT("  "));
}

//! Read materials from xml file to a vector
void read_materials_from_inp(std::string inp_path
                             ) {
    pugi::xml_document doc;
    // Load a file
	auto result = doc.load_file(inp_path.c_str());
	if (!result) {
	    std::cerr << "Error: file not found!" << std::endl;
	}
	pugi::xml_node root_node = doc.child("materials");
    // Iterate over materials records
	for (pugi::xml_node tool : root_node.children("material")) {
        std::string mname {tool.attribute("name").value()};
        auto m = new Materials(mname,
                    atof(tool.attribute("volume").value()),
                    atof(tool.attribute("volume").value()),
                    atof(tool.attribute("power").value()));
        m->namenuclides = get_node_array<std::string>(tool, "namenuclides");
        std::for_each(m->namenuclides.begin(), m->namenuclides.end(),
                      [m](std::string& nucname) {
                        m->indexnuclides.push_back(get_nuclidarray_index(
                                                      nucname));
                      }
        );
        // Read a concentration
        m->conc = get_node_array<udouble>(tool, "conc");
        if (check_for_node(tool, "dconc")) {
            std::vector<double> dconc = get_node_array<double>(tool, "dconc");
            int j = 0;
            if (dconc.size() <= m->conc.size())
                std::for_each(dconc.begin(), dconc.end(), [&j, m](double& dc) {
                    m->conc[j].Adddeviation(dc);
                    j++; }
                );
        }
        materials.push_back(std::unique_ptr<Materials>(m));
    } // for material
    return;
}

//! Look up all compositions and find a match with material name
//!
void matchcompositions() {
    for (auto&  mat: materials) {
        for  (size_t j = 0; j < compositions.size(); j++ ) {
            if (mat->Name()==compositions[j]->Name()) {
                mat->numcomposition = j;
                break;
            }
        }
    }
}

} //namespace openbps
