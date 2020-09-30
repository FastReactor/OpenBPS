#ifndef SRC_REACTIONS_H_
#define SRC_REACTIONS_H_
#include <vector>
#include <list>
#include <map>
#include <iostream>
#include <algorithm>
#include <memory>
#include "../extern/pugiData/pugixml.h"
#include "uncertainty.h"

namespace openbps {

//==============================================================================
// Global variables
//==============================================================================

class BaseCompostion;
class Composition;
class Sxs;
class Xslibs;
extern std::vector<
std::unique_ptr<Composition>> compositions;   //!< Vector with all compositions
extern int indexall;                          //!< Index of "all" record with
                                              //!< common shared data
extern std::map<std::string, int> composmap;  //!< Composition map with name:int
                                              //!< values
extern std::vector<Xslibs> externxslibs;      //!< External cross-section data
                                              //!< source

//==============================================================================
// Class with cross-section description
//==============================================================================
class Sxs {
public:
    //--------------------------------------------------------------------------
    // Constructors, destructors, factory functions
    Sxs() {}
    Sxs(pugi::xml_node node,const std::string& rxs, const std::string& redex);
    //--------------------------------------------------------------------------
    // Attributes
    std::string xsname;      //!< Name of cross-section/reaction
    std::string xstype;      //!< Cross-section/reaction type
    std::vector<udouble> rxs;//!< Reactions vector by energy group
    std::vector<udouble> xs_;//!< Cross-sections values by energy group

};

class Xslibs {
public:
    //--------------------------------------------------------------------------
    // Constructors, destructors, factory functions
    Xslibs() {}
    Xslibs(pugi::xml_node node);
    //--------------------------------------------------------------------------
    // Attributes
    std::vector<Sxs> xsdata; //!< Data source with cross-sections
    size_t numgroup;         //!< Energy group number for external xslib
    //--------------------------------------------------------------------------
    // Methods
    std::vector<double> get_egroups() { //!< Get an energy group structure for
    return energies_;}                  //!< external cross-section source
private:
    std::vector<double> energies_; //!< Energies for external cross-section lib
};

//==============================================================================
// Compositions classes descriptions
//==============================================================================
class BasicComposition {
public:
    //--------------------------------------------------------------------------
    // Constructors, destructors, factory functions
    BasicComposition() {}
    BasicComposition(std::string name, size_t nuclidnumber, size_t energynumber):
                     name_{name}, nuclide_number_{nuclidnumber},
                     energy_number_{energynumber} {}
    virtual ~BasicComposition() = default;
    //--------------------------------------------------------------------------
    // Methods
    //! Get name
    //!
    //! \return Composition name
    std::string Name() {return name_;}
    //! Get nuclid number
    //!
    //! \return number of nuclides
    size_t NuclidNumber() {return nuclide_number_;}
    //! Get a number of energy discretezation interval
    //!
    //! \return number of energy points
    size_t EnergyNumber() {return energy_number_;}

protected:
    //--------------------------------------------------------------------------
    // Attributes
    size_t nuclide_number_; //!< Nuclide number
    size_t energy_number_;  //!< Energy points number
    std::string name_;      //!< Composition Name
};

class Composition : public BasicComposition {
public:
    //--------------------------------------------------------------------------
    // Constructors, destructors, factory functions
    Composition(std::string name, size_t nuclidnumber, size_t energynumber):
                     BasicComposition(name, nuclidnumber, energynumber) {}

    Composition(pugi::xml_node node);
    //--------------------------------------------------------------------------
    // Methods
    //! Copy data from composition marked name "all" in xml file
    //!
    //! \param[in] externcompos external composition
    void deploy_all(Composition &externcompos);
    //! Calculate reaction rate for all reactions in xslib
    void get_reaction();
    //! Import data from external cross section source
    //!
    //! \param[in] implibs cross section library to calculate reaction rates
    void import_xsdata(Xslibs& implibs);
    //! Get spectrum energy distribution
    //!
    //! \return pair with energies and spectrum distr
    std::pair<std::vector<double>, std::vector<double>> get_fluxenergy();
    //--------------------------------------------------------------------------
    // Attributes
    std::vector<Sxs> xslib; //!< Cross-section/reactions data

private:
    //--------------------------------------------------------------------------
    // Methods
    //! Auxilary function to copy data from xslib
    //!
    //! \param[in] fmap external composition data
    //! \param[out] smap the copied data
    void depcopymap_(std::map<size_t, std::vector<double>>& fmap,
                     std::map<size_t, std::vector<double>>& smap);

    //! Calculate reaction-rate from cross-section data and flux
    //!
    //! \param[in] ixs cross-section data
    //! \param[in] extenergy energy range of cross section data
    //! \param[out]rxs calculated reaction rate
    void calculate_rr_(Sxs& ixs, const std::vector<double> &extenergy,
                       udouble& rxs);

    //--------------------------------------------------------------------------
    // Attributes
    std::map<size_t, std::vector<double>> energies_; //!< Energies
                                                     //!< discretization
    std::vector<udouble> spectrum_;                  //!< Energy spectrum
    std::vector<udouble> flux_;                      //!< Energy flux

};

//==============================================================================
// Non class methods
//==============================================================================
//! Auxilary function to copy data from xslib
//!
//! \param[in] node xml node element with xslib name
//! \param[in] rxs the sign of reactions (rxs)/ cross-section(xs) data
Sxs parse_xs_xml_
(pugi::xml_node node, const std::string& rxs,const std::string& redex);

//! Parse xslibs
//!
//! \param[in] node xml node with xslib data structure
//! \param[inout] xssource data source to store xslib in
void parse_xml_xslibs_(pugi::xml_node node, std::vector<Sxs>& xssource);

//! Read compositions from xml file
void read_reactions_xml();

//! Read an imported cross section data library from *.xml files
void read_importedlib_xml();

}

#endif /* SRC_REACTIONS_H_ */

