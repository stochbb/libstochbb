#ifndef INDEPENDENCETEST_HH
#define INDEPENDENCETEST_HH

#include "src/unittest.hh"


namespace stochbb {

class IndependenceTest: public UnitTest::TestCase
{
public:
  IndependenceTest();

  void testMinMax();

public:
  static UnitTest::TestSuite *suite();
};

}

#endif // INDEPENDENCETEST_HH
