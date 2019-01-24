#define CATCH_CONFIG_MAIN
#include "third-party/Catch/single_include/catch.hpp"

#include <iostream>

#include "base/Ptr.h"
#include "base/vector.h"
#include "tools/Random.h"

#include "ProgSynthBenchmarkSuite.h"

TEST_CASE("NumberIOTest", "[prog_synth_benchmarks]") {

  constexpr int SEED = 2;
  emp::Random random(SEED);

  // Make empty number io test
  NumberIOTest testzeros;
  std::cout << "All zeros test: "; 
  testzeros.Print(); std::cout << std::endl;
  REQUIRE(testzeros.GetInt() == 0);
  REQUIRE(testzeros.GetFloat() == 0);

  NumberIOTest rndtest(random, -50,50,-1000.0,100.0);
  REQUIRE(rndtest.Validate( -50,50,-1000.0,100.0));

  REQUIRE(rndtest.Evaluate(rndtest.GetInt() + rndtest.GetFloat()));
  REQUIRE(!rndtest.Evaluate(rndtest.GetInt() + rndtest.GetFloat() + 1));
}

TEST_CASE("SmallOrLargeTest", "[prog_synth_benchmarks]") {

  constexpr int SEED = 2;
  emp::Random random(SEED);

  // Make empty number io test
  SmallOrLargeTest emptytest;
  REQUIRE(emptytest.GetN() == 0);
  REQUIRE(emptytest.Evaluate("small"));
  REQUIRE(emptytest.Evaluate(SmallOrLargeTest::SMALL));

  SmallOrLargeTest randtest(random, -100, 100);
  std::cout << "Random test: ";
  randtest.Print(); std::cout << std::endl;
  REQUIRE(randtest.Validate(random, -100, 100));

}

// TODO
TEST_CASE("ForLoopIndexTest", "[prog_synth_benchmarks]") {
  constexpr int SEED = 2;
  emp::Random random(SEED);

  ForLoopIndexTest emptytest;

}

// TODO
TEST_CASE("CollatzNumbersTest", "[prog_synth_benchmarks]") {
  constexpr int SEED = 2;
  emp::Random random(SEED);

  CollatzNumbersTest emptytest;

}

// TODO
TEST_CASE("EvenSquaresTest", "[prog_synth_benchmarks]") {
  constexpr int SEED = 2;
  emp::Random random(SEED);

  EvenSquaresTest emptytest;

}

// TODO
TEST_CASE("WallisPiTest", "[prog_synth_benchmarks]") {
  constexpr int SEED = 2;
  emp::Random random(SEED);

  WallisPiTest emptytest;

}

// TODO
TEST_CASE("LastIndexOfZeroTest", "[prog_synth_benchmarks]") {
  constexpr int SEED = 2;
  emp::Random random(SEED);

  LastIndexOfZeroTest emptytest;

}

// TODO
TEST_CASE("VectorAverageTest", "[prog_synth_benchmarks]") {
  constexpr int SEED = 2;
  emp::Random random(SEED);

  VectorAverageTest emptytest;

}

// TODO
TEST_CASE("CountOddsTest", "[prog_synth_benchmarks]") {
  constexpr int SEED = 2;
  emp::Random random(SEED);

  CountOddsTest emptytest;

}

// TODO
TEST_CASE("MirrorImageTest", "[prog_synth_benchmarks]") {
  constexpr int SEED = 2;
  emp::Random random(SEED);

  MirrorImageTest emptytest;

}

// TODO
TEST_CASE("SumOfSquaresTest", "[prog_synth_benchmarks]") {
  constexpr int SEED = 2;
  emp::Random random(SEED);

  SumOfSquaresTest emptytest;

}

// TODO
TEST_CASE("VectorsSummedTest", "[prog_synth_benchmarks]") {
  constexpr int SEED = 2;
  emp::Random random(SEED);

  VectorsSummedTest emptytest;

}

// TODO
TEST_CASE("NegativeToZeroTest", "[prog_synth_benchmarks]") {
  constexpr int SEED = 2;
  emp::Random random(SEED);

  NegativeToZeroTest emptytest;

}

// TODO
TEST_CASE("DigitsTest", "[prog_synth_benchmarks]") {
  constexpr int SEED = 2;
  emp::Random random(SEED);

  DigitsTest emptytest;
}

// TODO
TEST_CASE("GradeTest", "[prog_synth_benchmarks]") {
  constexpr int SEED = 2;
  emp::Random random(SEED);

  GradeTest emptytest;

}

// TODO
TEST_CASE("MedianTest", "[prog_synth_benchmarks]") {
  constexpr int SEED = 2;
  emp::Random random(SEED);

  MedianTest emptytest;

}

// TODO
TEST_CASE("SmallestTest", "[prog_synth_benchmarks]") {
  constexpr int SEED = 2;
  emp::Random random(SEED);

  SmallestTest emptytest;

}