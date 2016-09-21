#ifndef DISTRIBUTIONTEST_HH
#define DISTRIBUTIONTEST_HH

#include "src/unittest.hh"

namespace stochbb {

class DistributionTest: public UnitTest::TestCase
{
public:
  DistributionTest();

  void testUniformQuantiles();
  void testNormalPDF();
  void testNormalCDF();
  void testNormalQuantiles();
  void testGammaPDF();
  void testGammaCDF();
  void testGammaQuantiles();
  void testInvGammaPDF();
  void testInvGammaCDF();
  void testInvGammaQuantiles();
  void testWeibullPDF();
  void testWeibullCDF();
  void testWeibullQuantiles();
  void testStudtPDF();
  void testStudtCDF();
  void testStudtQuantiles();

public:
  static UnitTest::TestSuite *suite();
};

}

#endif // DISTRIBUTIONTEST_HH
