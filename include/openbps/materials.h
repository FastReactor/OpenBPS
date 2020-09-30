#ifndef SRC_MATERIALS_H_
#define SRC_MATERIALS_H_
#include <vector>
#include <map>
#include <iostream>
#include <memory>
#include "../extern/pugiData/pugixml.h"
#include "uncertainty.h"

namespace openbps {

//==============================================================================
// Global variables
//==============================================================================

class BaseMaterials;
class Materials;

extern std::vector<std::unique_ptr<Materials>> materials; //!< vector consisting
                                                          //!< all materials

//==============================================================================
// Materials class description
//==============================================================================

class BaseMaterials {
public:
    //--------------------------------------------------------------------------
    // Constructors, destructors, factory functions
    BaseMaterials() {}
    BaseMaterials(std::string name, double volume, double power,
                  double mass): name_{name}, volume_{volume}, power_{power},
                                mass_{mass} {}
    virtual ~BaseMaterials() = default;
    //--------------------------------------------------------------------------
    // Methods
    //! Get name
    //! \return Material name
    std::string Name() {return name_;}

    //! Get the radiation power
    //! \return power [Wt]
    double Power()     {return power_;}

    //! Get material volume
    //! \return material volume [cm^3]
    double Volume()    {return volume_;}

    //! Get material mass
    //! \return material mass [g]
    double Mass()      {return mass_;}
private:
    //--------------------------------------------------------------------------
    // Attributes
    std::string name_; //!< Material name
    double mass_;      //!< Material mass [g]
    double volume_;    //!< Material volume [cm^3]
    double power_;     //!< Radiation power of material [Wt]
};

class Materials : public BaseMaterials {
public:
    //--------------------------------------------------------------------------
    // Constructors, destructors, factory functions
    Materials();
    Materials(std::string name, double volume, double power,
              double mass) : BaseMaterials(name, volume, power, mass) {}
    //--------------------------------------------------------------------------
    // Methods
    //! Add the material to xml node
    //!
    //! \param[inout] node pugixml node to write in
    void xml_add_material(pugi::xml_node node);
    //! Check a node corrections
    //!
    //! \param[in] node pugixml node to check it
    void check_node(pugi::xml_node node);
    //! Add nuclide to the material
    //!
    //! \param[in] extname name of new nuclide
    //! \param[in] exconc concentration of nuclide [atom/b-cm]
    void add_nuclide(const std::string& extname, udouble extconc);

    //--------------------------------------------------------------------------
    // Attributes
    std::vector<int> indexnuclides;       //!< Indexes of nuclides ownded
                                          //!< by the material in global nuclids
    std::vector<std::string> namenuclides;//!< Name of material nuclides
    std::vector<udouble> conc;            //!< Nuclide concentrations [atom/b-cm]
    int numcomposition {-1};              //!< Number of compoistion in container
    double normpower   {1.0};             //!< Value of power normalization

};

//==============================================================================
// Non class methods
//==============================================================================

//! Write down global materials vector information into xml file
//!
//! \param[in] xml_path filename of a directory to write
void form_materials_xml(std::string xml_path);

//! Read materials from xml file into a global materials vector
//!
//! \param[in]  inp_path input file layout
void read_materials_from_inp(std::string inp_path);
std::vector<Materials> parse_xml_materials();
//! Look up all compositions and find a match with material name
//!
void matchcompositions();
}//namespace openbps

#endif /* SRC_MATERIALS_H_ */
