#include "openbps/configure.h"
#include <string>
#include <algorithm>
#include <sstream>
#include "../extern/pugiData/pugixml.h"
#include "openbps/parse.h"
#include "openbps/timeproc.h"
#include "openbps/filter.h"
#include "openbps/uncertainty.h"
#include "openbps/chain.h"
#include "openbps/nuclide.h"
#include "openbps/materials.h"

namespace openbps {

//==============================================================================
// Global variable initialization
//==============================================================================

namespace configure {

std::string path_input;           //!< Directory where main .xml files resides
std::string path_output {
    ""};                          //!< Directory where output files are written
std::string chain_file;           //!< Chain-filename.xml
std::string nuclide_file;         //!< Nuclides database in *.xml
std::string reaction_file;        //!< Reaction-filename.xml
std::string inmaterials_file;     //!< Input materials-filename.xml
std::string outmaterials_file;    //!< Output materials-filename.xml
std::vector<std::string> libs;    //!< External libs with cross-sections
double timestep;                  //!< Time of simulation
int numstep;                      //!< Number of time step
double epb {1.e-03};              //!< accuracy of calculation
pugi::xml_document docx;          //!< Xml document
Mode calcmode {Mode::iteration};  //!< Type of solver time depended
                                  //!< exponental equation:
                                  //!< 1 - analogue Baetman method
                                  //!< 2 - iteration method by
                                  //!< E.F. Seleznev and I.V Chernova 2018
                                  //!< 3 - Chebyshev rational approximation by
                                  //!< Pussa, Josey, etc.
int order {8};                    //!< CRAM order in {8, 24} by default CRAM16
bool rewrite {true};              //!< Whether to rewrite a concentration
                                  //!< data by including nuclid from chain
bool outwrite{true};              //!< Write calculation result in file
std::vector<std::vector<std::array<double, 2>>>
dumpoutput;                       //!< Ouput dump
bool uncertantie_mod{false};      //!< Calculation mode with taking account 
                                  //!< uncertanties
bool decay_extra_out{false};      //!< Print out more information about
                                  //!< energy decay
bool verbose{false};              //!< Print information in out pipe
std::array<std::string,5>
header_names {"dt",
              "heat",
              "decay-rate",
              "dheat",
              "ddr"};             //!< Name of header in the ouput file

//==============================================================================
// Non class methods implementation
//==============================================================================


//! Parse init line
int parse_command_line(int argc, char* argv[])
{
    for (int i=1; i < argc; ++i) {
        std::string arg {argv[i]};
        if (arg[0] == '-') {
            if (arg == "-i" || arg == "--input") {
                path_input = "";
            }
            if (arg == "-o" || arg == "--output") {
                configure::outwrite = true;
            }
        }

    }

    return 0;
}

//! Read configure from XML file
void read_conigure_xml()
{
    using namespace configure;

    std::string configfile;
    if (path_input.length() > 1) {
        configfile = path_input + "configure.xml";
    } else {
        configfile = "configure.xml";
    }
    // Parse configure.xml file
    pugi::xml_document doc;
    auto result = doc.load_file(configfile.c_str());
    if (!result) {
        std::cerr << "Error while processing configure.xml file" ;
    }
    // Get root element
    pugi::xml_node root = doc.document_element();
    // Read a name of chain input *.xml file
    if (check_for_node(root, "chain")) {
        chain_file = get_node_value(root, "chain");
    }
    // Read a name of reaction input *.xml file (if presented)
    if (check_for_node(root, "reaction")) {
        reaction_file = get_node_value(root, "reaction");
    }
    // Read a name of nuclides database *.xml file
    if (check_for_node(root, "nuclides")) {
        nuclide_file = get_node_value(root, "nuclides");
    }
    // Read a name of materials input *.xml file
    if (check_for_node(root, "inpmaterials")) {
        inmaterials_file = get_node_value(root, "inpmaterials");
    }
    // Read a name of materials output *.xml file
    if (check_for_node(root, "outmaterials")) {
        outmaterials_file = get_node_value(root, "outmaterials");
    }
    // Read an external cross-sections *.xml files
    if (check_for_node(root, "impxslibs")) {
        for (pugi::xml_node tool :
             root.child("impxslibs").children("impxslib")) {
            libs.push_back(get_node_value(tool, "impxslib"));
        }
    }
    // Read an output directory name
    if (check_for_node(root, "output")) {
        path_output = get_node_value(root, "output");
    }
    // Read a number of steps in calculation
    if (check_for_node(root, "numbers")) {
        numstep = std::stoi(get_node_value(root, "numbers"));
    }
    // Read an overal simulation time
    //! TODO: improve time description
    if (check_for_node(root, "timestep")) {
        timestep = std::stod(get_node_value(root, "timestep"));
    }
    // Read a time duration
    //! TODO: improve time description
    if (check_for_node(root, "timerecord")) {
        Timexec duration(root.child("timerecord"));
        timestep = duration.get_seconds();
    }
    // Calculation accuracy (for iteration method)
    if (check_for_node(root, "epb")) {
        epb = std::stod(get_node_value(root, "epb"));
    }
    // Method selection
    if (check_for_node(root, "method")) {
        std::string method = get_node_value(root, "method");
        if (method == "baetman")  calcmode = Mode::baetman;
        if (method == "chebyshev") calcmode = Mode::chebyshev;
        if (method == "iteration") calcmode = Mode::iteration;
    }
    // CRAM order
    if (check_for_node(root, "cram_order")) {
        int ival = std::stoi(get_node_value(root, "cram_order"));
        switch (ival) {
        case 16:
            order = 8;
            break;
        case 48:
            order = 24;
        }
    }
    // Filters
    if (check_for_node(root, "filters")) {
        openbps::read_fitlers_from_xml(root.child("filters"));
    }
    // Simulation keys:
    if (check_for_node(root, "is_outrewrite")) {
        rewrite = get_node_value_bool(root, "is_outrewrite");
    }
    if (check_for_node(root, "decay_print")) {
        outwrite = get_node_value_bool(root, "decay_print");
    }
    if (check_for_node(root, "uncertanties")) {
        uncertantie_mod = get_node_value_bool(root, "uncertanties");
    }
    if (check_for_node(root, "decaykey")) {
        decay_extra_out = get_node_value_bool(root, "decaykey");
    }
}

} // namespace configure

//! Read input XML files
void read_input_xml()
{
    configure::read_conigure_xml();
}

} // namespace openbps





