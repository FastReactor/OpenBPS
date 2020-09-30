#include "openbps/filter.h"
#include <sstream>
#include <iostream>
#include <memory>
#include <algorithm>
#include "../extern/pugiData/pugixml.h"
#include "openbps/parse.h"


namespace openbps {

//==============================================================================
// Global variables
//==============================================================================
std::vector<Filter> filters;
std::unique_ptr<MaterialFilter> materialfilter;
std::unique_ptr<TimeFilter> timefilter;

//==============================================================================
// Filters class implementation
//==============================================================================

Filter::Filter(pugi::xml_node node)
{
   type = node.attribute("type").value();
   bins_ = get_node_array<std::string>(node, "filter");
}

//! Apply filter to results
void Filter::apply(const std::vector<std::string>&input,
           std::vector<int>& indices) {
    for (int i = 0; i < input.size(); i++) {
        auto search =
                std::find_if(bins_.begin(), bins_.end(),
                             [&input, i](const std::string& bname) {
            return bname == input[i];
        });
        if (search != bins_.end()) {
            indices.push_back(i);
        }
    }
}

MaterialFilter::MaterialFilter(pugi::xml_node node)
{
   type = "material";
   bins_ = get_node_array<std::string>(node, "filter");
}

//! Apply filter to results
void MaterialFilter::apply(const std::string &matname, bool& isValid) {
    auto search =
            std::find_if(bins_.begin(), bins_.end(),
                         [&matname](const std::string& bname) {
                             return bname == matname;
    });
    isValid = (search != bins_.end());

}

TimeFilter::TimeFilter(pugi::xml_node node)
{
   type = "time";
   bins_ = get_node_array<double>(node, "filter");
   if (bins_.size() % 2 == 1) {
       std::cout << "Intervals number should be even" << std::endl;
       bins_.erase(bins_.begin() + bins_.size() - 1, bins_.end());
   }

}

//! Apply filter to results
void TimeFilter::apply(double dt, int numstep, std::vector<int>& indices) {
    size_t j = 0;
    for (size_t k = 0; k < bins_.size() / 2; k++)
        while(j < numstep) {
            if ((j + 1) * dt <= bins_[2 * k + 1] &&
                    (j + 1) * dt > bins_[2 * k]) {
                 indices.push_back(j);

            }
            if ((j + 1) * dt > bins_[2 * k + 1])
                break;
             j++;
       }
}
//==============================================================================
// Non - class methods implementation
//==============================================================================
//! Reading data from configure.xml
void read_fitlers_from_xml(pugi::xml_node root_node) {
    // Proceed all filters
    for (pugi::xml_node tool : root_node.children("filter")) {
        std::string current_type;
        current_type = tool.attribute("type").value();
        if (current_type == "time") {
            timefilter = std::unique_ptr<TimeFilter>(new TimeFilter(tool));
        } else if (current_type == "material") {
            materialfilter = std::unique_ptr<MaterialFilter>(new MaterialFilter(tool));
        } else {
            Filter f(tool);
            filters.push_back(f);
        }
    }
}

} //namespace openbps
