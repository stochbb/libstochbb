#include "compoundtest.hh"
#include "src/api.hh"
#include "src/density.hh"
#include <iostream>

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
  Var X = normal(mu, 10);
  // analytic distribution.
  Var Y  = normal(10, std::sqrt(200));

  size_t N = 200;
  Eigen::VectorXd dX(N), dY(N);
  X.density().eval(-20, 40, dX);
  Y.density().eval(-20, 40, dY);
  for (size_t i=0; i<N; i++) {
    UT_ASSERT_NEAR_EPS(dX(i), dY(i), 1e-5);
  }

  X.density().evalCDF(-20,40, dX);
  Y.density().evalCDF(-20,40, dY);
  for (size_t i=0; i<N; i++) {
    UT_ASSERT_NEAR_EPS(dX(i), dY(i), 1e-5);
  }
}

void
CompoundTest::testNormalCompoundReduction() {
  // parameter distr.
  Var mu = normal(10, 10);
  // compound distr.
  Var X = normal(mu, 10);
  UT_ASSERT(X.density().is<AtomicDensity>());
}

void
CompoundTest::testNormalGammaCompound() {
  size_t N = 100;
  Var mu  = normal(0, 100);
  Var sig = invgamma(1,1);
  Var X = normal(mu, sig);
  Var Y = normal(0, 102);

  Eigen::VectorXd dX(N), dY(N);
  X.density().eval(-600, 600, dX);
  Y.density().eval(-600, 600, dY);
  std::cerr << "Got " << dX.transpose() << std::endl;
  std::cerr << "Exp " << dY.transpose() << std::endl;
  for (size_t i=0; i<N; i++) {
    UT_ASSERT_NEAR_EPS(dX(i), dY(i), 1e-4);
  }

  X.density().evalCDF(-600, 600, dX);
  Y.density().evalCDF(-600, 600, dY);
  for (size_t i=0; i<N; i++) {
    UT_ASSERT_NEAR_EPS(dX(i), dY(i), 1e-4);
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
  suite->addTest(new TestCaller<CompoundTest>("normal compound", &CompoundTest::testNormalCompound));
  suite->addTest(new TestCaller<CompoundTest>("normal compound reduction", &CompoundTest::testNormalCompoundReduction));
  suite->addTest(new TestCaller<CompoundTest>("normal gamma compound", &CompoundTest::testNormalGammaCompound));
  suite->addTest(new TestCaller<CompoundTest>("gamma compound", &CompoundTest::testGammaCompound));
  return suite;
}

