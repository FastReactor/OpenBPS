#include <iostream>
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
