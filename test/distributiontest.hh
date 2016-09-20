#ifndef DISTRIBUTIONTEST_HH
#define DISTRIBUTIONTEST_HH

#include "src/unittest.hh"

namespace stochbb {

class DistributionTest: public UnitTest::TestCase
{
public:
  DistributionTest();

  void testUniformQuantiles();
  void testNormalQuantiles();
  void testGammaPDF();
  void testGammaCDF();
  void testGammaQuantiles();
  void testInvGammaQuantiles();

public:
  static UnitTest::TestSuite *suite();
};

}

#endif // DISTRIBUTIONTEST_HH
