#ifndef SRC_FUNCTIONALS_H_
#define SRC_FUNCTIONALS_H_
#include <vector>
#include "uncertainty.h"
//
namespace openbps {
//==============================================================================
// Global variables and constants
//==============================================================================
constexpr double PWD {1.E-24};
constexpr int MAXITERATION {25};
constexpr double PWRC {6.2415093E+18};
//==============================================================================
// Non class methods implementation
//==============================================================================
//! Perform collapsing energy distribution from
//!
//! \param[in] yval vector of values according to
//! \param[in] xval vector of arguments to new
//! \param[in] xtarget vector of arguments
//! \return a vector with collapsing values
template <typename T>
std::vector<T> collapsing(const std::vector<double>& xval,
                          const std::vector<T>& yval,
                          const std::vector<double>& xtarget) {
    // Check a containers sizes
    if ((xval.size() < 2) || (xtarget.size() < 2)) {
        std::cerr << "Incorrect vector sizes to collapse it\n";
    }

    if (yval.size() != xval.size() - 1){
        std::cerr << "Incorrect vector sizes to collapse it\n";
    }
    // Output vector
    std::vector<T> ytarget(xtarget.size() - 1); // +1??

    if (xtarget[xtarget.size()-1] < xval[0]) {
        for (int i = 0; i < xtarget.size(); i++){
            ytarget[i] = yval[0];
        }
        return ytarget;
    }

    if (xtarget[0] > xval[xval.size() - 1]) {
        for (int i = 0; i < xtarget.size(); i++){
            ytarget[i] = yval[0];
        }
        return ytarget;
    }

    int jc {1};
    T val ;//{yval[0]};
    T normx ;//{0.0, 0.0};
    val = yval[0];
    normx = 0.0;
    double bx {xtarget[0]};
    // The main cycle
    for (int i = 1; i < xval.size(); i++) {

        while ((xtarget[jc] <= xval[i]) && (jc < xtarget.size())){
            normx = normx + val * (xtarget[jc] - bx);
            ytarget[jc-1] = normx / (xtarget[jc] - xtarget[jc-1]);
            bx = xtarget[jc];
            normx = 0;
            jc += 1;
        } // while

        if (jc >= xtarget.size()){
            return ytarget;
        } else {
            if (bx < xval[i])
                normx = normx + val * (xval[i] - bx);
            bx = xval[i];
            if (i < xval.size() - 1)  val = yval[i];
        }

    } // for
    // Finall checkup
    if (jc < xtarget.size()){
        for (int j = jc; j < xtarget.size(); j++) {
            normx = normx + val * (xtarget[j] - bx);
            ytarget[j-1] = normx / (xtarget[j] - xtarget[j-1]);
            bx = xtarget[j];
            normx = 0;
        } // for
    } // if
    return ytarget;

}

template <typename T>
std::vector<T> func(std::vector<double>& d1, std::vector<T>& tt) {
    std::vector<T> y(d1.size());
    for (size_t i=0; i<d1.size(); i++) {
        y[i] = d1[i] / tt[i];
    }
    return y;
}
// Future function to data translation
std::vector<std::pair<int, double>>
translating(const std::vector<double>& xval, const std::vector<double>& xtarget);


//! Form a new discreatization based on
//!
//! \param[in] xval arguments with
//! \param[in] yval values
//! \param[in] xtarget new arguments
//! \return a vector with new discretization
std::vector<double> transition(const std::vector<double>& xval,
                               const std::vector<double>& yval,
                               const std::vector<double>& xtarget);

} //namespace openbps
#endif /* SRC_FUNCTIONALS_H_ */
