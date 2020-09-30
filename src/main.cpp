#include <iostream>
#include <vector>
#include "openbps/nuclide.h"
#include "openbps/chain.h"
#include "openbps/uncertainty.h"
#include "openbps/configure.h"
#include "openbps/materials.h"
#include "openbps/reactions.h"
#include "openbps/functionals.h"
#include "openbps/matrix.h"
#include "openbps/timeproc.h"
#include "openbps/executer.h"
#include "xtensor/xarray.hpp"
#include <initializer_list>
using namespace std;

int main(int argc, char* argv[])
{
    using namespace openbps;
    Timer t;
    std::cout <<"Openbps start!\n";
    // Read input data
    configure::parse_command_line(argc, argv);
    // Read configure.xml
    read_input_xml();
    // Initialize solver
    executer::init_solver();
    // Run execution
    executer::run_solver();
    std::cout <<"Openbps finish!\n";
    return 0;
}
