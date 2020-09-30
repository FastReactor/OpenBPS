#ifndef SRC_CONFIGURE_H_
#define SRC_CONFIGURE_H_

#include "../extern/pugiData/pugixml.h"
#include <sstream>
#include <vector>
#include <array>


namespace openbps {
//==============================================================================
// Datatype definition
//==============================================================================
enum class Mode {
    baetman,  //!> Analogue direct Baetman method solution /* Not implemented in v.1.0*/
    iteration,//!> Iterative alogrithm
    chebyshev //!> Chebyshev rational approximation method
};
//==============================================================================
// Global variable declarations
//==============================================================================
namespace configure {

extern std::string path_input;            //!< Directory where main .xml
                                          //!< files resides
extern std::string path_output;           //!< Directory where output files
                                          //!< are written
extern std::string chain_file;            //!< Chain-filename.xml
extern std::string nuclide_file;          //!< Nuclides database in file *.xml
extern std::string reaction_file;         //!< Reaction-filename.xml
extern std::string inmaterials_file;      //!< Input materials-filename.xml
extern std::string outmaterials_file;     //!< Output materials-filename.xml
extern std::vector<std::string> libs;     //!< External libs with cross-sections
extern int numstep;                       //!< Number of substep per
                                          //!< one time step
extern double timestep ;                  //!< Length of time interval
                                          //!< for decesion
extern double epb;                        //!< Accuracy of calculation
extern pugi::xml_document docx;           //!< file handle for chain.xml
extern Mode calcmode;                     //!< Mode of calculation
extern int order;                         //!< CRAM order in {8, 24}
extern bool rewrite;                      //!< Whether to rewrite a
                                          //!< concentration data by including
                                          //!< nuclides from chain
extern bool outwrite;                     //!< Write calculation result in file
extern std::vector<std::vector<std::array<double, 2>>>
dumpoutput;                               //!< Ouput dump
extern bool uncertantie_mod;              //!< Calculation mode with
                                          //!< uncertanties taking account
extern bool decay_extra_out;              //!< Print out more information
                                          //!< about energy decay
extern bool verbose;                      //!< Print information in out pipe
extern std::array<std::string,5>
header_names;                             //!< Name of header in the ouput file


//==============================================================================
// Non class methods
//==============================================================================

//! Parse init line
int parse_command_line(int argc, char* argv[]);

//! Read configure from XML file
void read_conigure_xml();
}

//! Read input XML files
void read_input_xml();

} //namespace openbps

#endif /* SRC_CONFIGURE_H_ */
