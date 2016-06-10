#include "memorytest.hh"
#include "src/api.hh"

using namespace stochbb;
using namespace stochbb::UnitTest;

MemoryTest::MemoryTest()
  : TestCase()
{
  // pass...
}

void
MemoryTest::testNetReassemble() {
  Var A = normal(0,1), B1 = gamma(1,1), B2 = gamma(2,2), C = gamma(3,3);
  Var X = A + minimum(B1, B2) + C;
  for (size_t i=0; i<100; i++) {
    B1 = gamma(4,4);
    B2 = gamma(5,5);
    // X = A + minimum(B1, B2) + C;
  }
}


UnitTest::TestSuite *
MemoryTest::suite() {
  TestSuite *suite = new TestSuite("MemoryTest");
  suite->addTest(new TestCaller<MemoryTest>("network reassemble", &MemoryTest::testNetReassemble));
  return suite;
}
