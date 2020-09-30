#include "openbps/functionals.h"
#include <iostream>
#include <vector>
#include "openbps/chain.h"
#include "openbps/uncertainty.h"

namespace openbps {

//==============================================================================
// Non class methods implementation
//==============================================================================
// Future function to data translation
std::vector<std::pair<int, double>>
translating(const std::vector<double>& xval,
            const std::vector<double>& xtarget) {

    std::vector<std::pair<int, double>> ytranslater;

    return ytranslater;
}

// Form a new discreatization
std::vector<double> transition(const std::vector<double>& xval,
                               const std::vector<double>& yval,
                               const std::vector<double>& xtarget) {
    std::vector<double> result;
    result.resize(xtarget.size());
    // Input container sizes checkup
    if ((xtarget[0] > xval[xval.size()-1]) ||
            (xtarget[xtarget.size()-1] < xval[0])) {
        return result;
    }

    std::vector<double> xhalfv;
    // Define a middle average values
    for (int i=0; i<xval.size()-1; i++)
        xhalfv.push_back((xval[i]+ xval[i+1])/2);

    int icurr{0};
    double right {0.};
    double left {0.};
    // The main cycle
    for (int i=0; i < xhalfv.size(); i++) {

        if (xhalfv[i] < xtarget[0]) {
            result[0]=result[0] + yval[i];
        } else if (xhalfv[i] >= xtarget[xtarget.size()-1]) {
            result[xtarget.size()-1]=result[xtarget.size()-1] + yval[i];
        } else {
            while ((icurr < xtarget.size() - 1) &&
                   (xhalfv[i] > xtarget[icurr+1]))
                icurr++;
            left = (xhalfv[i] - xtarget[icurr]) /
                   (xtarget[icurr+1] - xtarget[icurr]);
            right = (-xhalfv[i] + xtarget[icurr+1]) /
                    ( xtarget[icurr+1] - xtarget[icurr]);
            result[icurr] = result[icurr] + yval[i] * (1 - left);
            result[icurr+1] = result[icurr+1] + yval[i] * (1 - right);
        }
    } // for

    return result;
}

} // namespace openbps




