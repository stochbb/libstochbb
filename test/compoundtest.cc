#include "compoundtest.hh"
#include "lib/api.hh"

using namespace stochbb;
using namespace stochbb::UnitTest;


CompoundTest::CompoundTest()
{
  // pass...
}

void
CompoundTest::testNormalCompound() {
  // parameter distr.
  Var mu = normal(10, 10);
  // compound distr.
  Var X  = normal(mu, 10);
  // analytic distribution.
  Var Y  = normal(10, std::sqrt(200));

  size_t N = 200;
  Eigen::VectorXd dX(N), dY(N);
  X.density().eval(-20, 40, dX);
  Y.density().eval(-20, 40, dY);
  for (size_t i=0; i<N; i++) {
    UT_ASSERT_NEAR_EPS(dX(i), dY(i), 1e-6);
  }

  X.density().evalCDF(-20,40, dX);
  Y.density().evalCDF(-20,40, dY);
  for (size_t i=0; i<N; i++) {
    UT_ASSERT_NEAR_EPS(dX(i), dY(i), 1e-5);
  }
}

void
CompoundTest::testGammaCompound() {
  Var k = uniform(0,4);
  Var X = gamma(5*k+5, 10);

  size_t N = 100;
  Eigen::VectorXd dX(N);
  X.density().eval(0, 1200, dX);
}

TestSuite *
CompoundTest::suite() {
  TestSuite *suite = new TestSuite("Compound");
  suite->addTest(new TestCaller<CompoundTest>("normalCompound", &CompoundTest::testNormalCompound));
  suite->addTest(new TestCaller<CompoundTest>("gammaCompound", &CompoundTest::testGammaCompound));
  return suite;
}
