#include "openbps/chain.h"

#include "../extern/pugiData/pugixml.h"
#include "openbps/parse.h"
#include "openbps/nuclide.h"
#include <iostream>
#include <tuple>
#include <list>
#include <algorithm>

namespace openbps {

//==============================================================================
// Chain implementation
//==============================================================================

Chain::Chain(pugi::xml_node node) {
    std::string current_name;
    size_t num;
    size_t col {0};
    for (pugi::xml_node tool : node.children("nuclide")) {
        current_name = get_node_value(tool, "name");
        if (current_name == "Nothing")
            break;
        // search nuclide in nuclides
        auto search = [current_name](std::unique_ptr<ChainNuclide>& item) {
            return item->name_ == current_name;
        };
        auto it = std::find_if(nuclides.begin(), nuclides.end(),
                               search);
        if (it != nuclides.end()) {
            // If nuclides have been read define a number
            num = std::distance(nuclides.begin(), it);
            this->name_idx.insert({current_name, num});
            col++;
        } else {
            // If nuclides.xml weren't read
            // we create new nuclide and add one to nuclids
            this->name_idx.insert({current_name, nuclides.size()});
            nuclides.push_back(std::unique_ptr<ChainNuclide>(new ChainNuclide(tool))); //!TODO make_unique
            num = nuclides.size() - 1;
            isNuclidePresent = false;
            col++;
        }
        nuclides[num]->read_from_chain_xml(tool);
    }

}

//! Get the information about chain nuclide names
std::vector<std::pair<int, std::string>> Chain::form_idx_name() {
    std::vector<std::pair<int, std::string>> out;

    for(auto& item : this->name_idx) {
        out.push_back(std::pair<int, std::string>(item.second, item.first));
    }
    return out;
}

//! Get the information about chain nuclide decay constant [sec^-1]
std::vector<std::pair<int, double>> Chain::form_idx_lambda() {
    std::vector<std::pair<int, double>> out;

    for(auto& item : this->name_idx) {
        out.push_back(std::make_pair(item.second,
                                   nuclides[item.second]->half_life.Real()));
    }
    return out;
}

//! Get the information about chain nuclide decay branching ratio
std::vector<std::tuple<int, int, udouble>> Chain::form_idx_decay() {
    std::vector<std::tuple<int, int, udouble>> out;

    for(auto& item : this->name_idx) {
        for (auto& elem : nuclides[item.second]->get_decaybr()) {
            out.push_back(std::make_tuple(item.second,
                                          this->name_idx[elem.first],
                                          elem.second));
    }
    }
    return out;
}

//! Get the information about chain from nuclide to nuclide by neutron reaction
/*
* to test the following parser use:
*
*
  auto react_map = chain.form_reaction();
  for(int i = 0; i < react_map["(n,gamma)"].size(); i++)
  {
      std::cout << react_map["(n,gamma)"][i].first << "  " << react_map["(n,gamma)"][i].second << "  "<< 1 << "  " <<std::endl;
  }
  output:
0  1  1
1  2  1
7  8  1
17  18  1
18  19  1
28  29  1
40  41  1
41  42  1
...
*/
std::map<std::string, std::vector<std::pair<int, int>>> Chain::form_reaction() {
    std::list<std::string> reactionlist;
    std::map<std::string, std::vector<std::pair<int, int>>> out;
    for (auto& item : this->name_idx) {
         for (auto& elem : nuclides[item.second]->get_reactions()) {
              auto search = std::find(reactionlist.begin(),
                                      reactionlist.end(), elem.first);
              auto searchelem = std::find_if(name_idx.begin(), name_idx.end(),
                                             [elem](std::pair<std::string,
                                             size_t> index) {
                                             return index.first == elem.second;
                                             });
              if (searchelem != name_idx.end()) {
              if (search == reactionlist.end()) {
                  std::vector<std::pair<int, int>> vec;
                  vec.push_back(std::make_pair(item.second,
                                               this->name_idx[elem.second]));
                  reactionlist.push_back(elem.first);
                  out[elem.first] = vec;
              } else {
                  out[elem.first].push_back(std::make_pair(
                                                 item.second,
                                                 this->name_idx[elem.second]));
              }
              }
         }
    }
    return out;
}

/*
* To test this parser (only Pu239 have 2000000.0):
  auto map = chain.form_yield_map();
  for(int i = 0; i < map[2000000.0].size(); i++)
  {
      std::cout << map[2000000.0][i][0] << "  " << map[2000000.0][i][1] << "  "<< map[2000000.0][i][2] << "  " <<std::endl;
  }
output
...
3557  1346  0
3557  1350  0
3557  1351  0
3557  1352  6.12916e-14
3557  1353  1.65977e-13
3557  1354  2.47967e-12
and so on
*/
//! Get the information about chain yield fissions
std::map<double, std::vector<
std::tuple<int, int, double>>> Chain::form_yield_map() {

    std::map<double,
            std::vector<
            std::tuple<int, int, double>>> out;
    std::vector<std::tuple<int, int, double>> insertion;
    std::tuple<int, int, double> elem;

    for(int i = 0; i < nuclides.size(); i++) {
        for(auto& item : nuclides[i]->get_nfy_energies()) {
            if(out.count(item) == 0)
                out.insert({item, insertion});
        }
    }
    for(int i = 0; i < nuclides.size(); i++) {
        for(int j = 0 ; j < nuclides[i]->get_nfy_energies().size(); j++) {
            auto n_map = nuclides[i]->get_product_data(j);
            for(auto& item : n_map) {
                elem = std::make_tuple(i, name_idx[item.first], item.second);
                out[nuclides[i]->get_nfy_energies()[j]].push_back(elem);
            }
        }
    }
    return out;
}
//! Local comparer
bool comp (std::pair <double, double> a,std::pair <double, double> b) {
    return a.first < b.first;
}

//! Convert fission yields to map for decay matrix
std::pair<std::vector<double>, std::vector<double>>
Chain::get_yield_map(size_t father, const std::string& daughter) {
    std::pair<std::vector<double>, std::vector<double>> out;
    std::vector<std::pair<double, double>> fracenergy_;
    std::vector<double> energies;
    std::vector<double> fractions;
    std::vector<double> nuclidenergies {nuclides[father]->get_nfy_energies()};
    for(int j = 0 ; j < nuclidenergies.size(); j++) {
        auto n_map = nuclides[father]->get_product_data(j);
        auto it = n_map.find(daughter);
        if (it != n_map.end()) {
            fracenergy_.push_back(
                        std::make_pair(nuclidenergies[j],
                                       it->second));
        }
    }
    std::sort(fracenergy_.begin(), fracenergy_.end(), comp);
    for (auto it=fracenergy_.begin(); it != fracenergy_.end(); it++){
        energies.push_back(it->first);
        fractions.push_back(it->second);
    }
    out = std::make_pair(energies, fractions);

    return out;
}

//! Get an index in chain nuclides array
int Chain::get_nuclide_index(const std::string& name) {
    int result {0};
    for (auto& item : name_idx)
        if (item.first == name) {
            return result;
        } else {
            result++;
        }
    return result;
}

// Read a chain from xml
Chain read_chain_xml(const std::string& filename) {
    pugi::xml_document doc;
    auto result = doc.load_file(filename.c_str());
    if (!result) {
        std::cerr << "Error: file not found!" << std::endl;
    }

    pugi::xml_node chain_node = doc.child("depletion_chain");

    Chain chainer(chain_node);
    
    return chainer;
}

} // namespace openbps



