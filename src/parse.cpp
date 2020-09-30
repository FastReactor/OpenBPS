#include "openbps/parse.h"
#include "../extern/pugiData/pugixml.h"

extern pugi::xml_document docx;

namespace openbps {

//==============================================================================
// Parse xml implementation
//==============================================================================

std::string get_node_value(pugi::xml_node node, const char *name) {
    // Search for either an attribute or child tag and get the data as a char*.
    const pugi::char_t *value_char;
    std::string nnode {node.name()};
    if (node.attribute(name)) {
        value_char = node.attribute(name).value();
    } else if (node.child(name)) {
        value_char = node.child_value(name);
    } else if (nnode == name) {
        value_char = node.child_value();
    } else {

        std::cerr << "Node \"" << name << "\" is not a member of the \""
                  << node.name() << "\" XML node";
    }
    std::string value{value_char};
    return value;
}

// Boolean value for function check_for_node
bool
check_for_node(pugi::xml_node node, const char *name) {
    return node.attribute(name) || node.child(name);
}

std::string join( std::vector<std::string> initList, const std::string& separator)
{
    std::string s;
    for(const auto& i : initList)
    {
        if(s.empty())
        {
            s = i;
        }
        else
        {
            s += separator + i;
        }
    }
    return s;
}
std::string joinDouble(std::vector<double> initList, const std::string& separator)
    {
        std::string s;
        std::ostringstream streamObj;
        for(const auto& i : initList)
        {
            if(s.empty())
            {
                streamObj << i;
                s = streamObj.str();
                streamObj.str("");
            }
            else
            {
                streamObj << i;
                s += separator + streamObj.str();
                streamObj.str("");
            }
        }
        return s;
    }

bool get_node_value_bool(pugi::xml_node node, const char *name) {
    if (node.attribute(name)) {
        return node.attribute(name).as_bool();
    } else if (node.child(name)) {
        return node.child(name).text().as_bool();
    } else {

        std::cerr << "Node \"" << name << "\" is not a member of the \""
                  << node.name() << "\" XML node";
    }
    return false;
}
// split a text string by a separator
std::vector<std::string> split(const std::string &s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}
// split a text string by a separator
std::vector<double> splitAtof(const std::string &s, char delimiter) {
    std::vector<double> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(atof(token.c_str()));
    }
    return tokens;
}


} //namespace openbps
