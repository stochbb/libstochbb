#ifndef MEMORYTEST_HH
#define MEMORYTEST_HH

#include "src/unittest.hh"

namespace stochbb {

class MemoryTest : public UnitTest::TestCase
{
public:
  MemoryTest();
  void testNetReassemble();

public:
  static UnitTest::TestSuite *suite();
};

}

#endif // MEMORYTEST_HH
