#ifndef TIMEPROC_H
#define TIMEPROC_H
#include "../extern/pugiData/pugixml.h"
#include "parse.h"
#include <chrono>
#include <string>

namespace openbps {

//==============================================================================
// Class to execution time code measure
//==============================================================================

class Timer {

    using clock_t = std::chrono::high_resolution_clock;
    using microseconds = std::chrono::microseconds;

public:
    //----------------------------------------------------------------------------
    //! Constructors, destructors, factory functions
    Timer() : start_(clock_t::now()) {}

    ~Timer() {
        const auto finish = clock_t::now();
        const auto us =
                std::chrono::duration_cast<microseconds>(finish - start_).count();
        std::cout << us << " us" << std::endl;
    }

private:
    //----------------------------------------------------------------------------
    //! Attributes
    const clock_t::time_point start_;
};

//==============================================================================
// Class to store time of problem execution
//==============================================================================

class Timexec {
public:
    //----------------------------------------------------------------------------
    //! Constructors, destructors, factory functions
    Timexec(int year, int day, int hour,
            int minute, int second) : year_(year), day_(day),
                                      hour_(hour), minute_(minute),
                                      second_(second) {}
    Timexec(pugi::xml_node node) {
        year_ = atoi(node.attribute("year").value());
        day_ = atoi(node.attribute("day").value());
        hour_ = atoi(node.attribute("hour").value());
        minute_ = atoi(node.attribute("minute").value());
        second_ = atoi(node.attribute("second").value());
    }
    //----------------------------------------------------------------------------
    //! Methods
    //! Get an overall time duration in seconds
    double get_seconds() {
        double result = year_ * 365.0 * 24.0 * 60.0 * 60.0 +
                        day_  * 24.0 * 60.0 * 60.0 +
                        hour_ * 60.0 * 60.0 +
                        minute_ * 60.0 +
                        second_;
        return result;
    }
private:
    //----------------------------------------------------------------------------
    //! Attributes
    int year_  {0}; //!< Time duration in years
    int day_   {0}; //!< Time duration in days
    int hour_  {0}; //!< Time duration in hours
    int minute_{0}; //!< Time duration in minutes
    int second_{0}; //!< Time duration in seconds
};

} // namespace openbps
#endif // TIMEPROC_H
