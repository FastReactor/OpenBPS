#ifndef FILTER_H
#define FILTER_H
#include "../extern/pugiData/pugixml.h"
#include "parse.h"
#include <memory>
#include <iostream>

namespace openbps {

class BaseFilter;
class Filter;
class MaterialFilter;
class TimeFilter;
//==============================================================================
// Global variables
//==============================================================================
extern std::vector<Filter> filters;           //!<Random filters
extern std::unique_ptr<MaterialFilter>
materialfilter;                               //!<Material output filter
extern std::unique_ptr<TimeFilter> timefilter;//!<Time filter for all problems

//==============================================================================
// Class to handle with calculation results
//==============================================================================

class BaseFilter
{
public:
    //--------------------------------------------------------------------------
    //! Constructors, destructors, factory functions
    BaseFilter() {}
    virtual ~BaseFilter() = default;
    //--------------------------------------------------------------------------
    //! Methods
    //! Check filter
    //virtual bool isapplying() = 0;
    //! Attributes
    //!
    std::string type;//!< Type of using filter

};

class Filter : public BaseFilter {
public:
    //--------------------------------------------------------------------------
    //! Constructors, destructors, factory functions
    Filter() {}
    Filter(pugi::xml_node node);
    //--------------------------------------------------------------------------
    //! Methods
    //!
    //! Apply filter to results
    //!
    //! \param[in] input vector with filtering names
    //! \param[out] indices indices to iterate for
    void apply(const std::vector<std::string>&input,
               std::vector<int> &indices);
    //! Check filter
    //bool isapplying() override {return false;}
    //----------------------------------------------------------------------------
    //! Attributes
    //!
    std::vector<std::string> bins_;//!< Filter conndition

};

//==============================================================================
// Special filter for materials
//==============================================================================

class MaterialFilter : public BaseFilter {
public:
    //----------------------------------------------------------------------------
    //! Constructors, destructors, factory functions
    MaterialFilter() {}
    MaterialFilter(pugi::xml_node node);
    //----------------------------------------------------------------------------
    //! Methods
    //!
    //! Apply filter to results
    //!
    //! \param[in] matname material name for applying
    void apply(const std::string& matname,
               bool &isValid);
    //! Check filter
    //virtual bool isapplying() override {return false;}
    //----------------------------------------------------------------------------
    //! Attributes
    //!
    std::vector<std::string> bins_;//!< Name of materials
};

//==============================================================================
// Special filter for time output manipulation
//==============================================================================

class TimeFilter : public BaseFilter {
public:
    //--------------------------------------------------------------------------
    //! Constructors, destructors, factory functions
    TimeFilter() {}
    TimeFilter(pugi::xml_node node);
    //--------------------------------------------------------------------------
    //! Methods
    //!
    //! Apply filter to results
    //!
    //! \param[in] dt time of calculation step [sec]
    //! \param[in] numstep number of calculation step
    //! \param[out] indices numerical values
    void apply(double dt, int numstep, std::vector<int>& indices);
    //! Check filter
    //virtual bool isapplying() override {return false;}
    //--------------------------------------------------------------------------
    //! Attributes
    //!
    std::vector<double> bins_;//!< Time intervals for which information be printed
};

//==============================================================================
// Non - class methods description
//==============================================================================
//! Reading data from configure.xml
//!
//! \param[in] root_node node containing filters information
void read_fitlers_from_xml(pugi::xml_node root_node);


} // namespace openbps

#endif // FILTER_H
