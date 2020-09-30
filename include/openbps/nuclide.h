#ifndef NUCLIDE_H
#define NUCLIDE_H

#include <iostream>
#include <string.h>
#include <vector>
#include <map>
#include <memory>
#include <sstream>
#include "../extern/pugiData/pugixml.h"
#include "parse.h"
#include "uncertainty.h"

namespace openbps {

//==============================================================================
// Global variables
//==============================================================================
extern bool isNuclidePresent;

class BaseNuclide;
class ChainNuclide;

extern std::vector<std::unique_ptr<ChainNuclide>> nuclides;

//==============================================================================
// Nuclide class description
//==============================================================================

class BaseNuclide {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors, factory functions
  BaseNuclide() {}
  BaseNuclide(std::string name) : name_{name} {}
  BaseNuclide(std::string name, int Z, int A, int M, double awr) : name_{name},
                                                               Z_ {Z},
                                                               A_ {A},
                                                               m_ {M}, 
                                                               awr_{awr} {}
  virtual ~BaseNuclide() = default;

  //----------------------------------------------------------------------------
  // Methods

  //----------------------------------------------------------------------------
  // Attributes

  std::string name_; //!< Name of nuclide, e.g. "U235"
  int Z_;            //!< Atomic number
  int A_;            //!< Mass number
  int m_;            //!< Metastable state
  double awr_;       //!< Atomic weight ratio
  int i_nuclide_;    //!< Index in the nuclides array
};

class ChainNuclide : public BaseNuclide {
public:
  // Types, aliases
  struct s_yield {
    double energy;                              //!< Energy of flying neutron                              
    std::map<std::string, double> product_data; //!< Fission product with yields
  };

  struct s_nfy {
    std::vector<double> energies;     //!< Energies of flying neutron
    std::vector<s_yield> yield_arr;   //!< Yield data
    std::map<std::string,
    std::vector<double>> yieldproduct;//!< Number of yields for \keyword nuclide

  };

  struct s_decay {
    std::string type;        //!< Decay type
    std::string target;      //!< Name of target nucleus
    udouble branching_ratio; //!< Branching ratio
    udouble q;               //!< Energy release per decay [eV]
  };

  struct s_reaction {
    std::string type;  //!< Type of nucleus-particle iteraction    
    std::string target;//!< Name of resulting nucleus
    double q;          //!< Energy releasing/absorbing at the reaction
  };
  //----------------------------------------------------------------------------
  // Constructors, destructors, factory functions

  ChainNuclide(){}
  ChainNuclide(std::string name, int Z, int A, int M, double awr) :
      BaseNuclide(name, Z, A, M, awr) {}
  ChainNuclide(pugi::xml_node nuclide_node);

  //----------------------------------------------------------------------------
  // Methods
  //! Read a nuclides from xml chain file
  //
  //! \param[in] nuclide_node pugixml node witch Chain headers
  void read_from_chain_xml(pugi::xml_node nuclide_node);

  //! Get a decay transition information
  //!
  //! \return vector with pair nuclide name and branching ratio
  std::vector<std::pair<std::string, udouble> > get_decaybr();

  //! Get a reactions for nuclide
  //!
  //! \return vector consisting pairs of reactions and reaction target
  std::vector<std::pair<std::string, std::string>> get_reactions();

  //! Get an energy release per reaction for nuclide
  //!
  //! \return map with pair reaction type name -> qvalue [eV]
  std::map<std::string, double> get_qvalue();

  //! Get a neutron fission energies
  //!
  //! \return list of neutron fission yields available enegies [eV]
  std::vector<double> get_nfy_energies();

  //! Get a neutron fission data by energy index
  //!
  //! \param[in] an index of energies array
  //! \return map with target nuclide name and its fraction
  std::map<std::string, double> get_product_data(int index);

  //! Get a fission yields for group energy discretization
  //!
  //! \return map with \keyword a target nuclide name and \value yields by
  //! neutron energy
  std::map<std::string,
  std::vector<double>> get_yield_product();
  //----------------------------------------------------------------------------
  // Attributes
  udouble half_life;
  udouble decay_energy;

private:
  //----------------------------------------------------------------------------
  // Attributes
  size_t decay_modes;    //!< Number of decay modes for nuclide
  size_t reactions;      //!< Number of nuclide-particle reactions
  std::vector<s_decay>
  decay_arr;             //!< Decay data for nuclide
  std::vector<s_reaction>
  reaction_arr;          //!< Data for particle-nuclide reactions (if presented)
  s_nfy nfy;             //!< Data for nuclear fission yields (if presented)
  //----------------------------------------------------------------------------
  // Methods
  //! Fill the yieldproduct map
  void make_product_yields_();
  //! Get a decay data from xml node
  s_decay parse_decay_(pugi::xml_node node);
  //! Get a reaction data from xml node
  s_reaction parse_reaction_(pugi::xml_node node);
  //! Get an yields from xml node
  s_yield parse_yield_(pugi::xml_node node);
  //! Get overall fission yield data from xml node
  s_nfy parse_nfy_(pugi::xml_node node);

};

//==============================================================================
// Non class methods
//==============================================================================
//! Get a nuclides from xml file
//!
//! \param[in] filepath xml file containing all nuclides
void read_nuclide_xml(const std::string& filepath);

//! Get an index in nucludes array by name
//!
//! \param[in] name of nuclide
//! \return result index of mentioned nuclide position
int get_nuclidarray_index(const std::string& name);

} // namespace openbps
#endif // CHAIN_H
