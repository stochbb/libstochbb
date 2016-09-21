#include "distributiontest.hh"
#include "src/api.hh"
#include "src/distribution.hh"

using namespace stochbb;
using namespace stochbb::UnitTest;


DistributionTest::DistributionTest()
  : UnitTest::TestCase()
{
  // pass...
}

void
DistributionTest::testUniformQuantiles() {
  Distribution dist(new UniformDistributionObj());
  Eigen::VectorXd param(2); param << 0, 1;
  double lower, upper;

  dist.quantile(lower, upper, 0.05, param);
  UT_ASSERT_NEAR_EPS(lower, 0.05, 1e-6);
  UT_ASSERT_NEAR_EPS(upper, 0.95, 1e-6);

  dist.quantile(lower, upper, 0.01, param);
  UT_ASSERT_NEAR_EPS(lower, 0.01, 1e-6);
  UT_ASSERT_NEAR_EPS(upper, 0.99, 1e-6);
}

void
DistributionTest::testNormalPDF() {
  Distribution dist(new NormalDistributionObj());
  Eigen::VectorXd param(2); param << 0, 1;
  Eigen::VectorXd pdf1(5), pdf2(5);
  pdf2 << 0.05399097, 0.24197072, 0.39894228, 0.24197072, 0.05399097;
  dist.pdf(-2, 3, pdf1, param);
  for (size_t i=0; i<5; i++) {
    UT_ASSERT_NEAR_EPS(pdf1(i), pdf2(i), 1e-8);
  }
}

void
DistributionTest::testNormalCDF() {
  Distribution dist(new NormalDistributionObj());
  Eigen::VectorXd param(2); param << 0, 1;
  Eigen::VectorXd cdf1(5), cdf2(5);
  cdf2 << 0.02275013, 0.15865525, 0.50000000, 0.84134475, 0.97724987;
  dist.cdf(-2, 3, cdf1, param);
  for (size_t i=0; i<5; i++) {
    UT_ASSERT_NEAR_EPS(cdf1(i), cdf2(i), 1e-8);
  }
}

void
DistributionTest::testNormalQuantiles() {
  Distribution dist(new NormalDistributionObj());
  Eigen::VectorXd param(2); param << 0, 1;
  double lower, upper;

  dist.quantile(lower, upper, 0.05, param);
  UT_ASSERT_NEAR_EPS(lower, -1.644854, 1e-6);
  UT_ASSERT_NEAR_EPS(upper,  1.644854, 1e-6);

  dist.quantile(lower, upper, 0.01, param);
  UT_ASSERT_NEAR_EPS(lower, -2.326348, 1e-6);
  UT_ASSERT_NEAR_EPS(upper,  2.326348, 1e-6);
}

void
DistributionTest::testGammaPDF() {
  Distribution dist(new GammaDistributionObj());
  Eigen::VectorXd param(3); param << 1, 1, 0;
  Eigen::VectorXd pdf1(5), pdf2(5);
  pdf2 << 1.00000000, 0.36787944, 0.13533528, 0.04978707, 0.01831564;
  dist.pdf(0, 5, pdf1, param);
  for (size_t i=0; i<5; i++) {
    UT_ASSERT_NEAR_EPS(pdf1(i), pdf2(i), 1e-8);
  }
}

void
DistributionTest::testGammaCDF() {
  Distribution dist(new GammaDistributionObj());
  Eigen::VectorXd param(3); param << 1, 1, 0;
  Eigen::VectorXd cdf1(5), cdf2(5);
  cdf2 << 0.0000000, 0.6321206, 0.8646647, 0.9502129, 0.9816844;
  dist.cdf(0, 5, cdf1, param);
  for (size_t i=0; i<5; i++) {
    UT_ASSERT_NEAR_EPS(cdf1(i), cdf2(i), 1e-7);
  }
}

void
DistributionTest::testGammaQuantiles() {
  Distribution dist(new GammaDistributionObj());
  Eigen::VectorXd param(3); param << 1, 1, 0;
  double lower, upper;

  dist.quantile(lower, upper, 0.05, param);
  UT_ASSERT_NEAR_EPS(lower, 0.05129329, 1e-6);
  UT_ASSERT_NEAR_EPS(upper, 2.99573227, 1e-6);

  dist.quantile(lower, upper, 0.01, param);
  UT_ASSERT_NEAR_EPS(lower, 0.01005034, 1e-6);
  UT_ASSERT_NEAR_EPS(upper, 4.60517019, 1e-6);
}

void
DistributionTest::testInvGammaPDF() {
  Distribution dist(new InvGammaDistributionObj());
  Eigen::VectorXd param(3); param << 1, 1, 0;
  Eigen::VectorXd pdf1(5), pdf2(5);
  pdf2 << 0.00000000, 0.36787944, 0.15163266, 0.07961459, 0.04867505;
  dist.pdf(0, 5, pdf1, param);
  for (size_t i=0; i<5; i++) {
    UT_ASSERT_NEAR_EPS(pdf1(i), pdf2(i), 1e-8);
  }
}

void
DistributionTest::testInvGammaCDF() {
  Distribution dist(new InvGammaDistributionObj());
  Eigen::VectorXd param(3); param << 1, 1, 0;
  Eigen::VectorXd cdf1(5), cdf2(5);
  cdf2 << 0.0000000, 0.3678794, 0.6065307, 0.7165313, 0.7788008;
  dist.cdf(0, 5, cdf1, param);
  for (size_t i=0; i<5; i++) {
    UT_ASSERT_NEAR_EPS(cdf1(i), cdf2(i), 1e-7);
  }
}

void
DistributionTest::testInvGammaQuantiles() {
  Distribution dist(new InvGammaDistributionObj());
  Eigen::VectorXd param(3); param << 1, 1, 0;
  double lower, upper;

  dist.quantile(lower, upper, 0.05, param);
  UT_ASSERT_NEAR_EPS(lower,  0.3338082, 1e-6);
  UT_ASSERT_NEAR_EPS(upper, 19.4957257, 1e-6);

  dist.quantile(lower, upper, 0.01, param);
  UT_ASSERT_NEAR_EPS(lower,  0.2171472, 1e-6);
  UT_ASSERT_NEAR_EPS(upper, 99.4991625, 1e-6);
}

void
DistributionTest::testWeibullPDF() {
  Distribution dist(new WeibullDistributionObj());
  Eigen::VectorXd param(3); param << 1, 1, 0;
  Eigen::VectorXd pdf1(5), pdf2(5);
  pdf2 << 1.00000000, 0.36787944, 0.13533528, 0.04978707, 0.01831564;
  dist.pdf(0, 5, pdf1, param);
  for (size_t i=0; i<5; i++) {
    UT_ASSERT_NEAR_EPS(pdf1(i), pdf2(i), 1e-8);
  }
}

void
DistributionTest::testWeibullCDF() {
  Distribution dist(new WeibullDistributionObj());
  Eigen::VectorXd param(3); param << 1, 1, 0;
  Eigen::VectorXd cdf1(5), cdf2(5);
  cdf2 << 0.0000000, 0.6321206, 0.8646647, 0.9502129, 0.9816844;
  dist.cdf(0, 5, cdf1, param);
  for (size_t i=0; i<5; i++) {
    UT_ASSERT_NEAR_EPS(cdf1(i), cdf2(i), 1e-7);
  }
}

void
DistributionTest::testWeibullQuantiles() {
  Distribution dist(new WeibullDistributionObj());
  Eigen::VectorXd param(3); param << 1, 1, 0;
  double lower, upper;

  dist.quantile(lower, upper, 0.05, param);
  UT_ASSERT_NEAR_EPS(lower, 0.05129329, 1e-6);
  UT_ASSERT_NEAR_EPS(upper, 2.99573227, 1e-6);

  dist.quantile(lower, upper, 0.01, param);
  UT_ASSERT_NEAR_EPS(lower, 0.01005034, 1e-6);
  UT_ASSERT_NEAR_EPS(upper, 4.60517019, 1e-6);
}

void
DistributionTest::testStudtPDF() {
  Distribution dist(new StudtDistributionObj());
  Eigen::VectorXd param(3); param << 1, 1, 0;
  Eigen::VectorXd pdf1(5), pdf2(5);
  pdf2 << 0.31830989, 0.15915494, 0.06366198, 0.03183099, 0.01872411;
  dist.pdf(0, 5, pdf1, param);
  for (size_t i=0; i<5; i++) {
    UT_ASSERT_NEAR_EPS(pdf1(i), pdf2(i), 1e-8);
  }
}

void
DistributionTest::testStudtCDF() {
  Distribution dist(new StudtDistributionObj());
  Eigen::VectorXd param(3); param << 1, 1, 0;
  Eigen::VectorXd cdf1(5), cdf2(5);
  cdf2 << 0.5000000, 0.7500000, 0.8524164, 0.8975836, 0.9220209;
  dist.cdf(0, 5, cdf1, param);
  for (size_t i=0; i<5; i++) {
    UT_ASSERT_NEAR_EPS(cdf1(i), cdf2(i), 1e-7);
  }
}

void
DistributionTest::testStudtQuantiles() {
  Distribution dist(new StudtDistributionObj());
  Eigen::VectorXd param(3); param << 1, 1, 0;
  double lower, upper;

  dist.quantile(lower, upper, 0.05, param);
  UT_ASSERT_NEAR_EPS(lower, -6.313752, 1e-6);
  UT_ASSERT_NEAR_EPS(upper,  6.313752, 1e-6);

  dist.quantile(lower, upper, 0.01, param);
  UT_ASSERT_NEAR_EPS(lower, -31.820516, 1e-6);
  UT_ASSERT_NEAR_EPS(upper,  31.820516, 1e-6);
}

TestSuite *
DistributionTest::suite() {
  TestSuite *suite = new TestSuite("Distributions");
  suite->addTest(new TestCaller<DistributionTest>("uniform quantiles", &DistributionTest::testUniformQuantiles));
  suite->addTest(new TestCaller<DistributionTest>("normal pdf", &DistributionTest::testNormalPDF));
  suite->addTest(new TestCaller<DistributionTest>("normal cdf", &DistributionTest::testNormalCDF));
  suite->addTest(new TestCaller<DistributionTest>("normal quantiles", &DistributionTest::testNormalQuantiles));
  suite->addTest(new TestCaller<DistributionTest>("gamma pdf", &DistributionTest::testGammaPDF));
  suite->addTest(new TestCaller<DistributionTest>("gamma cdf", &DistributionTest::testGammaCDF));
  suite->addTest(new TestCaller<DistributionTest>("gamma quantiles", &DistributionTest::testGammaQuantiles));
  suite->addTest(new TestCaller<DistributionTest>("inv-gamma pdf", &DistributionTest::testInvGammaPDF));
  suite->addTest(new TestCaller<DistributionTest>("inv-gamma cdf", &DistributionTest::testInvGammaCDF));
  suite->addTest(new TestCaller<DistributionTest>("inv-gamma quantiles", &DistributionTest::testInvGammaQuantiles));
  suite->addTest(new TestCaller<DistributionTest>("weibull pdf", &DistributionTest::testWeibullPDF));
  suite->addTest(new TestCaller<DistributionTest>("weibull cdf", &DistributionTest::testWeibullCDF));
  suite->addTest(new TestCaller<DistributionTest>("weibull quantiles", &DistributionTest::testWeibullQuantiles));
  suite->addTest(new TestCaller<DistributionTest>("student-t pdf", &DistributionTest::testStudtPDF));
  suite->addTest(new TestCaller<DistributionTest>("student-t cdf", &DistributionTest::testStudtCDF));
  suite->addTest(new TestCaller<DistributionTest>("student-t quantiles", &DistributionTest::testStudtQuantiles));
  return suite;
}

