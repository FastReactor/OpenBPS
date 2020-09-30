#ifndef PARSE_H
#define PARSE_H

#include <iostream>
#include <string.h>
#include <vector>
#include <map>
#include <sstream>
#include "../extern/pugiData/pugixml.h"
#include <chrono>

namespace openbps {

//==============================================================================
// Methods to work with xml inputs
//==============================================================================

std::string get_node_value(pugi::xml_node node, const char *name);
bool get_node_value_bool(pugi::xml_node node, const char *name);
bool check_for_node(pugi::xml_node node, const char *name);

std::vector<double> splitAtof(const std::string &s, char delimiter);
std::vector<std::string> split(const std::string &s, char delimiter);
std::string joinDouble(std::vector<double> initList, const std::string& separator);
std::string join( std::vector<std::string> initList, const std::string& separator);
template <typename T>
std::vector<T> get_node_array(pugi::xml_node node, const char* name,
                              bool lowercase=false)
{
    // Get value of node attribute/child
    std::string s {get_node_value(node, name)};
    // Read values one by one into vector
    std::stringstream iss {s};
    T value;
    std::vector<T> values;
    while (iss >> value)
        values.push_back(value);

    return values;
}

} // namespace openbps

#endif // PARSE_H
