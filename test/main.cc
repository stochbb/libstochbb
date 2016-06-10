#include "src/unittest.hh"
#include <iostream>

#include "memorytest.hh"
#include "convolutiontest.hh"
#include "compoundtest.hh"
#include "exactsamplertest.hh"
#include "src/logger.hh"

using namespace stochbb;
using namespace stochbb::UnitTest;


int main(int argc, char *argv[]) {

  //Logger::addHandler(IOLogHandler());

  // Construct test-runner
  TestRunner runner(std::cout);

  // Add suites
  runner.addSuite(MemoryTest::suite());
  runner.addSuite(ConvolutionTest::suite());
  runner.addSuite(CompoundTest::suite());
  runner.addSuite(ExactSamplerTest::suite());

  // Exec tests:
  runner();

  return 0;
}
