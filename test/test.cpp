#include "gtest/gtest.h"
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
#include "xtensor/xbuilder.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xview.hpp"

using namespace std::complex_literals;
using namespace openbps;

TEST(ConfigureTest, test4) {
  using namespace openbps;
  pugi::xml_document doc;
  const std::string input_filepath = "./test/test_xml/configure.xml";
  auto result = doc.load_file(input_filepath.c_str());
  // Get root element
  pugi::xml_node root = doc.document_element();
  int numstep = std::stoi(get_node_value(root, "numbers"));
  double timestep = std::stod(get_node_value(root, "timestep"));
  EXPECT_EQ(numstep, 30);
  EXPECT_EQ(numstep, 30);
  
}

TEST(XtensorTEST, test)
{
    using namespace xt;
    xarray<double> a = {{ 2, 1, 1},
                        {-1, 1,-1},
                        { 1, 2, 3}};

    xarray<double> vec = {2, 3, -10};
    xarray<double> expected = {3, 1, -5};
    auto rows {xt::row(a, 1)};
    auto columns {xt::col(a, 2)};
    auto d1=1.0 - columns;
    auto d2=xt::exp2(d1);
    double res1 {-1.0};
    double res2 {3.0};
    EXPECT_EQ(d2(0), 1.0);
    EXPECT_EQ(xt::sum(xt::col(a, 2))(0), res2);
}

TEST(Uncerttest, test1)
{
    Uncertainty<double> u1(1.0, 1.0);
    Uncertainty<double> u2(2.0, 0.5);
    udouble u3 {u1 + 2.0};
    EXPECT_EQ(u2.Real(), 2.0);
    EXPECT_EQ(u2.Dev(), 0.5);
    EXPECT_EQ(u3.Real(), 3.0);
    EXPECT_EQ(u3.Dev(), 1.0);
}

TEST(Uncerttest, test2)
{
    std::vector<udouble> u4 {{1.,1.},{2.,2.},
                             {3.,3.},{4.,4.}};
    EXPECT_EQ(u4[0].Real(), 1.0);
    EXPECT_EQ(u4[0].Dev(), 1.0);
    EXPECT_EQ(u4[2].Real(), 3.0);
    EXPECT_EQ(u4[2].Dev(), 3.0);
}

TEST(Matrixtest, test3)
{
    BaseMatrix<double> m1({{1., 0.}, {1., 1.}});
    BaseMatrix<double> m2({{4., -2.}, {2., 5.}});
    m1[1][1] = -3.0;
    double a = m1[1][1];
    m1 = m1 + m2;
    BaseMatrix<double> m3(2,2);
    BaseMatrix<double> m4 = m1 + m2;
    m4 = m1 + m3;
    m4 = m1 | m2;
    BaseMatrix<double> m5(4,4);
    m5 = std::move(m2);
    std::vector<BaseMatrix<double>> vm;
    vm.push_back(BaseMatrix<double>(4,2));
    vm.push_back(std::move(m5));
    vm.push_back(m3);
    vm.push_back(m4);
    for (auto& v : vm[0]) {
        v = 2;
    }
    EXPECT_EQ(a, -3.0);
    EXPECT_EQ(m1[1][1], 2.0);
    EXPECT_EQ(m4[1][3], 5.0);
    EXPECT_EQ(vm[0][0][0], 2);
}
