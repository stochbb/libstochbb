#ifndef COMPOUNDTEST_HH
#define COMPOUNDTEST_HH

#include "src/unittest.hh"

namespace stochbb {

class CompoundTest: public UnitTest::TestCase
{
public:
  CompoundTest();

  void testNormalCompound();
  void testNormalCompoundReduction();
  void testNormalGammaCompound();
  void testNormalInvGammaCompound();
  void testGammaCompound();

public:
  static UnitTest::TestSuite *suite();
};

}

#endif // COMPOUNDTEST_HH
