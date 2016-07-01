#include "independencetest.hh"
#include "src/api.hh"
#include "src/randomvariable.hh"
#include <set>
#include <iostream>


using namespace stochbb;
using namespace stochbb::UnitTest;


IndependenceTest::IndependenceTest()
{
  // pass...
}

void
IndependenceTest::testMinMax() {
  Var X = gamma(1, 100);
  Var Y = gamma(1, 100);
  Var minXY = minimum(X,Y);

  UT_ASSERT(independent(X,Y));
  UT_ASSERT(! independent(X, minXY));
  UT_ASSERT(! independent(Y, minXY));
}

UnitTest::TestSuite *
IndependenceTest::suite() {
  TestSuite *suite = new TestSuite("Independence");
  suite->addTest(new TestCaller<IndependenceTest>("minimum/maximum", &IndependenceTest::testMinMax));
  return suite;
}
