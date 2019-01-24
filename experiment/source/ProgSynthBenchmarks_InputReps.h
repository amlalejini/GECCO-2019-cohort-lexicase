#ifndef PROG_SYNTH_BENCHMARK_INPUT_REPRESENTATIONS_H
#define PROG_SYNTH_BENCHMARK_INPUT_REPRESENTATIONS_H

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <utility>
#include <unordered_set>

#include "base/Ptr.h"
#include "base/vector.h"
#include "tools/Random.h"
#include "tools/random_utils.h"
#include "tools/math.h"
#include "tools/vector_utils.h"
#include "tools/sequence_utils.h"
#include "tools/string_utils.h"
#include "tools/stats.h"

#include "parser.hpp"

#include "TestCaseSet.h"
#include "Utilities.h"

/*

Genetic representations for different input types used across the programming
synthesis benchmarks.

Ideas:
- Validate take restriction lambdas(?) 

Input/Output types:
- NumberIO: Pair<integer, float>
- SmallOrLarge: Integer
- ForLoopIndex: Array<Integer, 3>
- Compare String Lengths: Array<String, 3>
- Double Letters: String
- Collatz Numbers: Integer
- Replace Space with Newline: String
- String Differences: Array<String, 2>
- Even Squares: Integer
- Wallis Pi: Integer
- String Lengths Backwards: Vector<String>
- Last Index of Zero: Vector<Integer>
- Vector Average: Vector<Float>
- Count Odds: Vector<Integer>
- Mirror Image: Array<Vector<Integer>, 2>
- Super Anagrams: Array<String, 2>
- Sum of Squares: Integer
- Vectors Summed: Array<Vector<Integer>, 2>
- X-Word Lines: Pair<Integer, String>
- Pig Latin: String
- Negative To Zero: Vector<Integer>
- Scrabble Score: String
- Word Stats: File
- Checksum: String
- Digits: Integer
- Grade: Array<Integer, 5>
- Median: Array<Integer, 3>
- Smallest: Array<Integer, 4>
- Syllables: String

Aggregated types:
- Pair<TYPE1, TYPE2>
- Integer
- String
- Array<TYPE, SIZE>
- Vector<TYPE>
*/

// Problem-specific, static variables

// -- PROBLEM: Small or Large --
constexpr int SmallOrLarge__SMALL_THRESH = 1000;
constexpr int SmallOrLarge__LARGE_THRESH = 2000;
const std::string SmallOrLarge__SMALL_STR = "small";
const std::string SmallOrLarge__LARGE_STR = "large";
const std::string SmallOrLarge__NONE_STR = "";

const std::string Grade__A_STR = "A";
const std::string Grade__B_STR = "B";
const std::string Grade__C_STR = "C";
const std::string Grade__D_STR = "D";
const std::string Grade__F_STR = "F";

// ================== Problem I/O Type Aliases ==================
using Problem_NumberIO_input_t = std::pair<int, double>;
using Problem_NumberIO_output_t = double;

using Problem_SmallOrLarge_input_t = int;
using Problem_SmallOrLarge_output_t = std::string;

using Problem_ForLoopIndex_input_t = std::array<int, 3>;
using Problem_ForLoopIndex_output_t = emp::vector<int>;

using Problem_CompareStringLengths_input_t = std::array<std::string, 3>;
using Problem_CompareStringLengths_output_t = bool;

using Problem_DoubleLetters_input_t = std::string;
using Problem_DoubleLetters_output_t = std::string;

using Problem_CollatzNumbers_input_t = int;
using Problem_CollatzNumbers_output_t = int;

using Problem_ReplaceSpaceWithNewline_input_t = std::string;
using Problem_ReplaceSpaceWithNewline_output_t = std::pair<std::string, int>;

using Problem_StringDifferences_input_t = std::array<std::string, 2>;
using Problem_StringDifferences_output_t = std::string;

using Problem_EvenSquares_input_t = int;
using Problem_EvenSquares_output_t = int;

using Problem_WallisPi_input_t = int;
using Problem_WallisPi_output_t = double;

using Problem_StringLengthsBackwards_input_t = emp::vector<std::string>;
using Problem_StringLengthsBackwards_output_t = emp::vector<size_t>;

using Problem_LastIndexOfZero_input_t = emp::vector<int>;
using Problem_LastIndexOfZero_output_t = int;

using Problem_VectorAverage_input_t = emp::vector<double>;
using Problem_VectorAverage_output_t = double;

using Problem_CountOdds_input_t = emp::vector<int>;
using Problem_CountOdds_output_t = int;

using Problem_MirrorImage_input_t = std::array<emp::vector<int>, 2>;
using Problem_MirrorImage_output_t = bool;

using Problem_SuperAnagrams_input_t = std::array<std::string, 2>;
using Problem_SuperAnagrams_output_t = bool;

using Problem_SumOfSquares_input_t = int;
using Problem_SumOfSquares_output_t = int;

using Problem_VectorsSummed_input_t = std::array<emp::vector<int>, 2>;
using Problem_VectorsSummed_output_t = emp::vector<int>;

using Problem_XWordLines_input_t = std::pair<int, std::string>;
using Problem_XWordLines_output_t = std::string;

using Problem_PigLatin_input_t = std::string;
using Problem_PigLatin_output_t = std::string;

using Problem_NegativeToZero_input_t = emp::vector<int>;
using Problem_NegativeToZero_output_t = emp::vector<int>;

using Problem_ScrabbleScore_input_t = std::string;
using Problem_ScrabbleScore_output_t = int;

using Problem_Checksum_input_t = std::string;
using Problem_Checksum_output_t = std::string;

using Problem_Digits_input_t = int;
using Problem_Digits_output_t = emp::vector<int>;

using Problem_Grade_input_t = std::array<int, 5>;
using Problem_Grade_output_t = std::string;

using Problem_Median_input_t = std::array<int, 3>;
using Problem_Median_output_t = int;

using Problem_Smallest_input_t = std::array<int, 4>;
using Problem_Smallest_output_t = int;

using Problem_Syllables_input_t = std::string;
using Problem_Syllables_output_t = std::string;


// ================ Problem mutators ================

// Problem_NumberIO_input_t = std::pair<int, double>;
struct Pair_IntDouble_Mutator {
  int MIN_INT;
  int MAX_INT;

  double MIN_DOUBLE;
  double MAX_DOUBLE;

  double PER_INT_RATE;
  double PER_DOUBLE_RATE;

  size_t Mutate(emp::Random & rnd, std::pair<int, double> & mut_pair) {
    size_t muts = 0;
    if (rnd.P(PER_INT_RATE)) {
      mut_pair.first = rnd.GetInt(MIN_INT, MAX_INT+1);
      ++muts;
    }
    if (rnd.P(PER_DOUBLE_RATE)) {
      mut_pair.second = rnd.GetDouble(MIN_DOUBLE, MAX_DOUBLE);
      ++muts;
    }
    return muts;
  }
};

struct Int_Mutator {
  int MIN_INT;
  int MAX_INT;

  double PER_INT_RATE;

  size_t Mutate(emp::Random & rnd, int & mut_int) {
    if (rnd.P(PER_INT_RATE)) {
      mut_int = rnd.GetInt(MIN_INT, MAX_INT+1);
      return 1;
    }
    return 0;
  }
};

// ================ Problem Organism Classes ================

// Test org base class
class TestOrg_Base {
  public:

    struct Phenotype {
      emp::vector<double> test_scores;
      double total_score;

      emp::vector<bool> test_passes;
      size_t num_passes;
      size_t num_fails;

      void Reset(size_t s=0) {
        test_scores.clear();
        test_scores.resize(s, 0);
        test_passes.clear();
        test_passes.resize(s, false);
        num_passes = 0;
        num_fails = 0;
      }

      void RecordScore(size_t id, double val) {
        emp_assert(id < test_scores.size());
        total_score += val;
        test_scores[id] = val;
      }

      void RecordPass(size_t id, bool pass) {
        emp_assert(id < test_passes.size());
        if (pass) ++num_passes;
        else ++num_fails;
        test_passes[id] = pass;
      }
    };

    using phenotype_t = Phenotype;

  protected:
    phenotype_t phenotype;

  public:

    phenotype_t & GetPhenotype() { return phenotype; }
    const phenotype_t & GetPhenotype() const { return phenotype; }

    virtual ~TestOrg_Base() = default;

    virtual void CalcOut() = 0;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: NumberIO
// - Input type: [double, integer]
// - Output type: double 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Problem_NumberIO_input_t GenRandomTestInput_NumberIO(emp::Random & rand, const std::pair<int, int> & int_range, const std::pair<double, double> & double_range) {
  emp_assert(double_range.first < double_range.second);
  emp_assert(int_range.first < int_range.second);
  return Problem_NumberIO_input_t{rand.GetInt(int_range.first, int_range.second), rand.GetDouble(double_range.first, double_range.second)};
}
  
Problem_NumberIO_output_t GenCorrectOut_NumberIO(const Problem_NumberIO_input_t & input) {
  return input.first + input.second;
}

void SetCorrectOut_NumberIO(const Problem_NumberIO_input_t & input, Problem_NumberIO_output_t & output) {
  output = input.first + input.second;
}

/// Calculate pass/fail score on NumberIO problem.
std::pair<double, bool> CalcScorePassFail_NumberIO(const Problem_NumberIO_output_t & correct_test_output, double sub) {
  const bool pass = (double)sub == correct_test_output;
  return {(double)pass, pass};
}

std::pair<double, bool> CalcScoreGradient_NumberIO(const Problem_NumberIO_output_t & correct_test_output, double sub, double MAX_ERROR) {
  // If output is correct, return a score of 1.0 and mark that submission passes.
  if (correct_test_output == sub) {
    return {1.0, true}; 
  } else { // Otherwise, return {score=[0:1], false}
    double error = emp::Abs(correct_test_output - sub);
    emp_assert(error != 0, "Error shouldn't equal zero here");
    double score = (error <= MAX_ERROR) ? 1 - (error/MAX_ERROR) : 0;
    return {score, false};
  }  
}

/// ProblemOrg: NumberIO
/// NumberIO: Pair<integer, float>
class TestOrg_NumberIO : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = Problem_NumberIO_input_t;
    using out_t = Problem_NumberIO_output_t;

  protected:
    genome_t genome;
    out_t out;
  
  public:
    TestOrg_NumberIO(const genome_t & _g) : genome(_g), out() { ; }

    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    out_t & GetCorrectOut() { return out; }
    const out_t & GetCorrectOut() const { return out; }

    // std::pair<int, double>;
    int GetTestIntegerInput() const { return genome.first; }
    double GetTestDoubleInput() const { return genome.second; }

    void SetOut(const out_t & _out) { out = _out; }

    void CalcOut() { SetCorrectOut_NumberIO(genome, out); }

    void Print(std::ostream & os=std::cout) const {
      os << genome.first << "," << genome.second;
    }
};

struct ProblemUtilities_NumberIO {
  using input_t = Problem_NumberIO_input_t;
  using output_t = Problem_NumberIO_output_t;
  
  using testcase_set_t = TestCaseSet<input_t,output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  emp::vector<emp::Ptr<TestOrg_NumberIO>> testingset_pop;

  emp::vector<emp::vector<output_t>> population_validation_outputs;

  // A few useful things for use within a test evaluation
  emp::Ptr<TestOrg_NumberIO> cur_eval_test_org;
  bool submitted;
  double submitted_val; // if going to do string thing, we can have a submission_str.

  Pair_IntDouble_Mutator mutator;

  emp::vector<std::function<double(TestOrg_NumberIO &)>> lexicase_fit_set;

  ProblemUtilities_NumberIO() 
    : testing_set(ProblemUtilities_NumberIO::LoadTestCaseFromLine),
      training_set(ProblemUtilities_NumberIO::LoadTestCaseFromLine),
      submitted(false), submitted_val(0.0)
  { ; }

  ~ProblemUtilities_NumberIO() {
    for (size_t i = 0; i < testingset_pop.size(); ++i) testingset_pop[i].Delete();
  }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  void ResetTestEval() {
    submitted = false;
    submitted_val = 0.0;
  }

  void Submit(double val) {
    submitted = true;
    submitted_val = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const std::string & line) {
    emp::vector<std::string> split_line = emp::slice(line, ',');
    input_t input;
    output_t output;
    input.second = std::atof(split_line[0].c_str());
    input.first = std::atof(split_line[1].c_str());
    output = std::atof(split_line[2].c_str());
    return {input, output};
  }

  void GenerateTestingSetPop() {
    for (size_t i = 0; i < testing_set.GetSize(); ++i) {
      testingset_pop.emplace_back(emp::NewPtr<TestOrg_NumberIO>(testing_set.GetInput(i)));
      testingset_pop[i]->CalcOut();
    }
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << in.first << "," << in.second;
    os << "\"";
  }

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: Small or large
// - Input: integer
// - Output: string {'small', 'large', ''}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Generate random test input
Problem_SmallOrLarge_input_t GenRandomTestInput_SmallOrLarge(emp::Random & rand, const std::pair<int,int> & int_range) {
  emp_assert(int_range.first < int_range.second);
  return rand.GetInt(int_range.first, int_range.second+1);
}

Problem_SmallOrLarge_output_t GenCorrectOut_SmallOrLarge(const Problem_SmallOrLarge_input_t & input) {
  if (input < SmallOrLarge__SMALL_THRESH) return SmallOrLarge__SMALL_STR;
  else if (input >= SmallOrLarge__LARGE_THRESH) return SmallOrLarge__LARGE_STR;
  else return SmallOrLarge__NONE_STR;
}

/// SmallOrLarge: Integer
class TestOrg_SmallOrLarge : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = Problem_SmallOrLarge_input_t;
    using out_t = Problem_SmallOrLarge_output_t;

  protected:
    genome_t genome;
    out_t out;

  public:
    TestOrg_SmallOrLarge(const genome_t & _g) : genome(_g), out() { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    out_t & GetCorrectOut() { return out; }
    const out_t & GetCorrectOut() const { return out; }

    void CalcOut() { out = GenCorrectOut_SmallOrLarge(genome); }

    void Print(std::ostream & os=std::cout) const {
      os << genome;
    }
};

struct ProblemUtilities_SmallOrLarge { 
  using this_t = ProblemUtilities_SmallOrLarge;
  using problem_org_t = TestOrg_SmallOrLarge;
  using input_t = Problem_SmallOrLarge_input_t;
  using output_t = Problem_SmallOrLarge_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  emp::vector<emp::Ptr<problem_org_t>> testingset_pop;

  emp::vector<emp::vector<output_t>> population_validation_outputs;

  // --- Useful during a test evaluation ---
  emp::Ptr<problem_org_t> cur_eval_test_org;
  bool submitted;
  std::string submitted_str;

  // Mutation
  Int_Mutator mutator;

  // Selection
  emp::vector<std::function<double(problem_org_t &)>> lexicase_fit_set;

  ProblemUtilities_SmallOrLarge()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_str("")
  { ; }

  ~ProblemUtilities_SmallOrLarge() {
    for (size_t i = 0; i < testingset_pop.size(); ++i) testingset_pop[i].Delete();
  }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  void ResetTestEval() {
    submitted = false;
    submitted_str = "";
  }

  void Submit(const std::string & val) {
    submitted = true;
    submitted_str = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const std::string & line) {
    emp::vector<std::string> split_line = emp::slice(line + " ", ',');
    input_t input;    // int
    output_t output;  // std::string
    // std::cout << "LINE=" << line << std::endl;
    input = std::atof(split_line[0].c_str());
    output = split_line[1];
    if (!(output == " " || output == "small " || output == "large ")) {
      std::cout << "ERROR! Bad output ("<<output<<") from line: " << line << std::endl;
      exit(-1);
    }
    if (output == " ") output = "";
    else if (output == "small ") output = "small";
    else if (output == "large ") output = "large";
    return {input, output};
  }

  void GenerateTestingSetPop() {
    for (size_t i = 0; i < testing_set.GetSize(); ++i) {
      testingset_pop.emplace_back(emp::NewPtr<problem_org_t>(testing_set.GetInput(i)));
      testingset_pop[i]->CalcOut();
    }
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << in;
    os << "\"";
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: ForLoopIndex
// - Input: std::array<int, 3>;
// - Output: emp::vector<int>;
// - Description: Given 3 integer inputs (start, end, step), print the integers in the sequence:
//              - n0 = start
//              - ni = ni-1 + step
//              - for each ni < end
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Generate random test input
Problem_ForLoopIndex_input_t GenRandomTestInput_ForLoopIndex(emp::Random & rand,
                                                         const std::pair<int,int> & start_end_range,
                                                         const std::pair<int,int> & step_range,
                                                         bool promise_multistep=false) {
  // Guarantee that start comes before end.
  int start, end, step;
  start = rand.GetInt(start_end_range.first, start_end_range.second+1);
  end = rand.GetInt(start_end_range.first, start_end_range.second+1);
  step = rand.GetInt(step_range.first, step_range.second+1);

  // Need to follow the following rules:
  // (1) start < end               // -> Start should come before end
  // (2) start + (20xstep)+1 > end // -> Enumeration should take no more than 20 steps
  if (promise_multistep) {
    while (true) {
      if (end < start) std::swap(end, start);
      if ( ( start + (20*step) + 1 > end ) && ( (start + step) < end ) ) break; // Guarantee that output will at minimum be [start, start+step, ...]
      start = rand.GetInt(start_end_range.first, start_end_range.second+1);
      end = rand.GetInt(start_end_range.first, start_end_range.second+1);
      step = rand.GetInt(step_range.first, step_range.second+1);
    }
  } else {
    while (true) {
      if (end < start) std::swap(end, start);
      if ( (start + (20*step) + 1 > end) ) break;
      start = rand.GetInt(start_end_range.first, start_end_range.second+1);
      end = rand.GetInt(start_end_range.first, start_end_range.second+1);
      step = rand.GetInt(step_range.first, step_range.second+1);
    }
  }

  return {start, end, step};
}

/// Generate correct output
Problem_ForLoopIndex_output_t GenCorrectOut_ForLoopIndex(const Problem_ForLoopIndex_input_t & input) {
  emp::vector<int> out;
  for (int i = input[0]; i < input[1]; i+= input[2]) {
    out.emplace_back(i);
  }
  return out;
}

class TestOrg_ForLoopIndex : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = std::array<int, 3>;
    using out_t = Problem_ForLoopIndex_output_t;

  protected:
    genome_t genome;
    out_t out;

  public:
    TestOrg_ForLoopIndex(const genome_t & _g) : genome(_g), out() { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { out = GenCorrectOut_ForLoopIndex(genome); }

    out_t & GetCorrectOut() { return out; }
    const out_t & GetCorrectOut() const { return out; }   

    void Print(std::ostream & os=std::cout) {
      os << genome[0] << "," << genome[1] << "," << genome[2];
    }
};

struct ProblemUtilities_ForLoopIndex { 
  using this_t = ProblemUtilities_ForLoopIndex;
  using problem_org_t = TestOrg_ForLoopIndex;
  using input_t = Problem_ForLoopIndex_input_t;
  using output_t = Problem_ForLoopIndex_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  emp::vector<emp::Ptr<problem_org_t>> testingset_pop;

  emp::vector<emp::vector<output_t>> population_validation_outputs;

  // --- Useful during a test evaluation ---
  emp::Ptr<problem_org_t> cur_eval_test_org;
  bool submitted;
  emp::vector<int> submitted_vec;

  // Mutation - Handle here...
  int MIN_START_END;
  int MAX_START_END;
  int MIN_STEP;
  int MAX_STEP;
  double MUT_RATE;  

  size_t Mutate(emp::Random & rnd, Problem_ForLoopIndex_input_t & mut_input) {
    size_t muts = 0;

    int start = mut_input[0];
    int end = mut_input[1];
    int step = mut_input[2];

    // Mutate start
    if (rnd.P(MUT_RATE)) { // Mutate start?
      start = rnd.GetInt(MIN_START_END, MAX_START_END+1);
      muts++;
    }
    // Mutate end
    if (rnd.P(MUT_RATE)) {
      end = rnd.GetInt(MIN_START_END, MAX_START_END+1);
      muts++;
    }
    // Mutate step
    if (rnd.P(MUT_RATE)) {
      step = rnd.GetInt(MIN_STEP, MAX_STEP+1);
      muts++;
    }
    
    // Ensure that start still comes before end.  
    if (end == start) end++;
    if (end < start) std::swap(start, end);
    if (!(start + (20*step) + 1 > end)) {
      // Move end closer to start
      end = rnd.GetInt(start, start + (20*step)+1);
    }

    mut_input[0] = start;
    mut_input[1] = end;
    mut_input[2] = step;
    
    return muts;
  }

  // Selection
  emp::vector<std::function<double(problem_org_t &)>> lexicase_fit_set;

  ProblemUtilities_ForLoopIndex()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_vec()
  { ; }

  ~ProblemUtilities_ForLoopIndex() {
    for (size_t i = 0; i < testingset_pop.size(); ++i) testingset_pop[i].Delete();
  }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  void ResetTestEval() {
    submitted = false;
    submitted_vec.clear();
  }

  void Submit(int val) {
    submitted = true;
    submitted_vec.emplace_back(val);
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const std::string & line) {
    emp::vector<std::string> split_line = emp::slice(line, ',');
    input_t input;    // int
    output_t output;  // std::string
    // Start = line[0]
    input[0] = std::atof(split_line[0].c_str());
    // End = line[1]
    input[1] = std::atof(split_line[1].c_str());
    // Step = line[2]
    input[2] = std::atof(split_line[2].c_str());
    output = GenCorrectOut_ForLoopIndex(input);  
    return {input, output};
  }

  void GenerateTestingSetPop() {
    for (size_t i = 0; i < testing_set.GetSize(); ++i) {
      testingset_pop.emplace_back(emp::NewPtr<problem_org_t>(testing_set.GetInput(i)));
      testingset_pop[i]->CalcOut();
    }
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  std::pair<double, bool> CalcScoreGradient(const output_t & correct_test_output, const output_t & sub) {
    const double max_dist = emp::Max(correct_test_output.size(), sub.size());
    double dist = emp::calc_edit_distance(correct_test_output, sub);
    if (dist == 0) {
      return {1.0, true};
    } else {
      return {(max_dist - dist)/max_dist, false};
    }
  } // todo - test this

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << in[0] << "," << in[1] << "," << in[2];
    os << "\"";
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: CompareStringLengths
// - Input: std::array<std::string, 3>;
// - Output: bool;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// /// Generate random test input
Problem_CompareStringLengths_input_t GenRandomTestInput_CompareStringLengths(emp::Random & rand,
                                                                             std::pair<size_t, size_t> str_size_range) {
  // Valid characters: \n, \t, [32, 127)
  emp::vector<char> valid_chars = {'\n', '\t'};
  for (size_t i = 32; i < 127; ++i) valid_chars.emplace_back((char)i);

  // String 1
  size_t str_size = rand.GetUInt(str_size_range.first, str_size_range.second+1);
  std::string str1;
  for (size_t i = 0; i < str_size; ++i) str1.push_back(valid_chars[rand.GetUInt(valid_chars.size())]);

  // String 2
  str_size = rand.GetUInt(str_size_range.first, str_size_range.second+1);
  std::string str2;
  for (size_t i = 0; i < str_size; ++i) str2.push_back(valid_chars[rand.GetUInt(valid_chars.size())]);

  // String 3
  str_size = rand.GetUInt(str_size_range.first, str_size_range.second+1);
  std::string str3;
  for (size_t i = 0; i < str_size; ++i) str3.push_back(valid_chars[rand.GetUInt(valid_chars.size())]);
  // std::cout << "---- RANDOM STRING ----" << std::endl;
  // std::cout << "String 1: " << str1 << std::endl;
  // std::cout << "String 2: " << str2 << std::endl;
  // std::cout << "String 3: " << str3 << std::endl;
  return {str1, str2, str3};
}

/// Generate correct output
Problem_CompareStringLengths_output_t GenCorrectOut_CompareStringLengths(const Problem_CompareStringLengths_input_t & input) {
  if (input[0].size() < input[1].size() && input[1].size() < input[2].size()) return true;
  return false;
}

/// Compare String Lengths: Array<String, 3>
class TestOrg_CompareStringLengths: public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = std::array<std::string, 3>;
    using out_t = Problem_CompareStringLengths_output_t;

  protected:
    genome_t genome;
    out_t out;

  public:
    TestOrg_CompareStringLengths(const genome_t & _g) : genome(_g), out() { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { out = GenCorrectOut_CompareStringLengths(genome); }

    out_t & GetCorrectOut() { return out; }
    const out_t & GetCorrectOut() const { return out; }   

    void Print(std::ostream & os=std::cout) {
      os << "[{BEGIN-STR}" << genome[0] << "{END-STR},{BEGIN-STR}" << genome[1] << "{END-STR},{BEGIN-STR}" << genome[2] << "{END-STR}]";
    }

};

struct ProblemUtilities_CompareStringLengths { 
  using this_t = ProblemUtilities_CompareStringLengths;
  using problem_org_t = TestOrg_CompareStringLengths;
  using input_t = Problem_CompareStringLengths_input_t;
  using output_t = Problem_CompareStringLengths_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  emp::vector<emp::Ptr<problem_org_t>> testingset_pop;
  emp::vector<emp::vector<output_t>> population_validation_outputs;

  // --- Useful during a test evaluation ---
  emp::Ptr<problem_org_t> cur_eval_test_org;
  bool submitted;
  bool submitted_val;

  // Mutation - Handle here...
  size_t MIN_STR_LEN;
  size_t MAX_STR_LEN;
  double PER_SITE_INS_RATE;
  double PER_SITE_DEL_RATE;
  double PER_SITE_SUB_RATE; 
  double PER_STR_SWAP_RATE; 
  emp::vector<char> valid_chars;
  
  size_t Mutate(emp::Random & rnd, input_t & mut_input) {
    size_t muts = 0;

    for (size_t i = 0; i < mut_input.size(); ++i) {
      std::string & str = mut_input[i];
      int expected_size = (int)str.size();
      std::string new_string = "";

      for (size_t s = 0; s < str.size(); ++s) {
        
        // Per-site insertions
        if (rnd.P(PER_SITE_INS_RATE) && (expected_size+1 <= (int)MAX_STR_LEN)) {
          // Insert a random character.
          new_string.push_back(valid_chars[rnd.GetUInt(valid_chars.size())]);
          ++expected_size;
          ++muts;
        }
        // Per-site deletions
        if (rnd.P(PER_SITE_DEL_RATE) && (expected_size-1 >= (int)MIN_STR_LEN)) {
          --expected_size;
          ++muts;
          continue;
        }
        size_t whead = new_string.size();
        new_string.push_back(str[s]); 
        // Per-site substitutions
        if (rnd.P(PER_SITE_SUB_RATE)) {
          new_string[whead] = valid_chars[rnd.GetUInt(valid_chars.size())];
          ++muts;
        }
      }
      emp_assert(new_string.size() >= (size_t)MIN_STR_LEN);
      emp_assert(new_string.size() <= (size_t)MAX_STR_LEN);
      mut_input[i] = new_string;
    }

    for (size_t i = 0; i < mut_input.size(); ++i) {
      if (rnd.P(PER_STR_SWAP_RATE) && mut_input.size() > 1) {
        ++muts;
        // Swap position i with a random different position
        size_t other = rnd.GetUInt(mut_input.size());
        while (other == i) other = rnd.GetUInt(mut_input.size());
        std::swap(mut_input[i], mut_input[other]);
      }
    }

    return muts;
  } // todo - test

  // Selection
  emp::vector<std::function<double(problem_org_t &)>> lexicase_fit_set;

  ProblemUtilities_CompareStringLengths()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val(false)
  { 
    valid_chars = {'\n', '\t'};
    for (size_t i = 32; i < 127; ++i) valid_chars.emplace_back((char)i); 
  }

  ~ProblemUtilities_CompareStringLengths() {
    for (size_t i = 0; i < testingset_pop.size(); ++i) testingset_pop[i].Delete();
  }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  void ResetTestEval() {
    submitted = false;
    submitted_val = false;
  }

  void Submit(bool val) {
    submitted = true;
    submitted_val = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    // Load input.
    input[0] = line[0];
    input[1] = line[1];
    input[2] = line[2];
    // Load output.
    if (line[3] == "false") output = false;
    else if (line[3] == "true") output = true;
    else {
      std::cout << "ERROR when loading test case output (" << line[3] << "). Exiting." << std::endl;
      exit(-1);
    }
    return {input, output};
  }

  void GenerateTestingSetPop() {
    for (size_t i = 0; i < testing_set.GetSize(); ++i) {
      testingset_pop.emplace_back(emp::NewPtr<problem_org_t>(testing_set.GetInput(i)));
      testingset_pop[i]->CalcOut();
    }
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"[{BEGIN-STR}" << in[0] << "{END-STR},{BEGIN-STR}" << in[1] << "{END-STR},{BEGIN-STR}" << in[2] << "{END-STR}]\"";
    // os << "\"" << in[0] << "\"," << "\"" << in[1] << "\"," << "\"" << in[2] << "\"";
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem:
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Double Letters: String
class TestOrg_DoubleLetters : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = std::string;
  protected:
    genome_t genome;

  public:
    TestOrg_DoubleLetters(const genome_t & _g) : genome(_g) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { ; }
};

struct ProblemUtilities_DoubleLetters { emp::vector<std::function<double(TestOrg_DoubleLetters &)>> lexicase_fit_set; };

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: Collatz Numbers
// - Input: int
// - Output: int
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Generate random test input
Problem_CollatzNumbers_input_t GenRandomTestInput_CollatzNumbers(emp::Random & rand,
                                                                 const std::pair<int, int> & num_range) {
  return rand.GetInt(num_range.first, num_range.second+1);
}

Problem_CollatzNumbers_output_t GenCorrectOut_CollatzNumbers(const Problem_CollatzNumbers_input_t & input) {
  emp_assert(input > 0);
  // std::cout << "Generating collatz seq - input = " << input << std::endl;
  emp::vector<int> seq;
  seq.emplace_back(input);
  while (seq.back() != 1) {
    if (seq.back() % 2 == 0) seq.emplace_back(seq.back() / 2);
    else seq.emplace_back((3*seq.back())+1);
  }
  // std::cout << "Generating collatz seq - output = " << seq.size() << std::endl;
  return seq.size();
}

/// Collatz Numbers: Integer
class TestOrg_CollatzNumbers : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = Problem_CollatzNumbers_input_t;
    using out_t = Problem_CollatzNumbers_output_t;

  protected:
    genome_t genome;
    out_t out;

    emp::Ptr<std::unordered_map<int, int>> out_cache;

    void CalcOutCache() {
      // std::cout << "Calc out w/cache!" << std::endl;
      emp_assert(out_cache != nullptr);
      if (emp::Has(*out_cache, genome)) {
        out = (*out_cache)[genome];
      } else {
        out = GenCorrectOut_CollatzNumbers(genome);
        (*out_cache)[genome] = out;
      }
    }

    void CalcOutNoCache() {
      // std::cout << "Calc out w/ no cache!" << std::endl;
      out = GenCorrectOut_CollatzNumbers(genome);
    }
    
  public:
    TestOrg_CollatzNumbers(const genome_t & _g) : genome(_g), out(), out_cache(nullptr) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void SetCache(emp::Ptr<std::unordered_map<int, int>> cache_ptr) {
      out_cache = cache_ptr;
    }

    void CalcOut() { 
      if (out_cache) {
        CalcOutCache();
      } else {
        CalcOutNoCache();
      }
    }

    out_t & GetCorrectOut() { return out; }
    const out_t & GetCorrectOut() const { return out; }  

    void Print(std::ostream & os=std::cout) {
      os << genome;
    }

};

struct ProblemUtilities_CollatzNumbers {  
  using this_t = ProblemUtilities_CollatzNumbers;
  using problem_org_t = TestOrg_CollatzNumbers;
  using input_t = Problem_CollatzNumbers_input_t;
  using output_t = Problem_CollatzNumbers_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  emp::Ptr<std::unordered_map<int, int>> out_cache;

  emp::vector<emp::Ptr<problem_org_t>> testingset_pop;

  emp::vector<emp::vector<output_t>> population_validation_outputs;

  // // --- Useful during a test evaluation ---
  emp::Ptr<problem_org_t> cur_eval_test_org;
  bool submitted;
  int submitted_val;
  
  // Error
  int MAX_ERROR;

  // // Mutation - Handle here...
  int MIN_NUM;
  int MAX_NUM;
  double NUM_SUB_RATE;
  
  size_t Mutate(emp::Random & rnd, input_t & mut_input) {
    size_t muts = 0;

    if (rnd.P(NUM_SUB_RATE)) {
      ++muts;
      mut_input = rnd.GetUInt(MIN_NUM, MAX_NUM);
    }

    return muts;
  } // todo - test

  // Selection
  emp::vector<std::function<double(problem_org_t &)>> lexicase_fit_set;

  ProblemUtilities_CollatzNumbers()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val(0)
  { 
    out_cache = emp::NewPtr<std::unordered_map<int, int>>(); 
  }

  ~ProblemUtilities_CollatzNumbers() {
    out_cache.Delete();
    for (size_t i = 0; i < testingset_pop.size(); ++i) testingset_pop[i].Delete();
  }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  void ResetTestEval() {
    submitted = false;
    submitted_val = 0;
  }

  void Submit(int val) {
    submitted = true;
    submitted_val = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    // Load input.
    input = (int)std::atof(line[0].c_str());
    // Load output.
    output = (int)std::atof(line[1].c_str());
    emp_assert(output == GenCorrectOut_CollatzNumbers(input));
    return {input, output};
  }

  void GenerateTestingSetPop() {
    for (size_t i = 0; i < testing_set.GetSize(); ++i) {
      testingset_pop.emplace_back(emp::NewPtr<problem_org_t>(testing_set.GetInput(i)));
      testingset_pop.back()->SetCache(out_cache);
      testingset_pop[i]->CalcOut();
    }
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  std::pair<double, bool> CalcScoreGradient(const output_t & correct_test_output, const output_t & sub) {
    if (correct_test_output == sub) {
      return {1.0, true};
    } else {
      double error = (double)emp::Abs(correct_test_output - sub);
      emp_assert(error != 0, "Error shouldn't be zero here!");
      double score = (error <= MAX_ERROR) ? 1 - (error/(double)MAX_ERROR) : 0.0;
      return {score, false};
    }
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << in;
    os << "\"";
  }

};



/// Replace Space with Newline: String
class TestOrg_ReplaceSpaceWithNewline : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = std::string;
  protected:
    genome_t genome;

  public:
    TestOrg_ReplaceSpaceWithNewline(const genome_t & _g) : genome(_g) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { ; }
};

struct ProblemUtilities_ReplaceSpaceWithNewline { emp::vector<std::function<double(TestOrg_ReplaceSpaceWithNewline &)>> lexicase_fit_set; };



/// String Differences: Array<String, 2>
class TestOrg_StringDifferences : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = std::array<std::string, 2>;
  protected:
    genome_t genome;

  public:
    TestOrg_StringDifferences(const genome_t & _g) : genome(_g) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { ; }
};

struct ProblemUtilities_StringDifferences { emp::vector<std::function<double(TestOrg_StringDifferences &)>> lexicase_fit_set; };



/// Even Squares: Integer
class TestOrg_EvenSquares : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = int;
  protected:
    genome_t genome;

  public:
    TestOrg_EvenSquares(const genome_t & _g) : genome(_g) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { ; }
};

struct ProblemUtilities_EvenSquares { emp::vector<std::function<double(TestOrg_EvenSquares &)>> lexicase_fit_set; };




/// Wallis Pi: Integer
class TestOrg_WallisPi : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = int;
  protected:
    genome_t genome;

  public:
    TestOrg_WallisPi(const genome_t & _g) : genome(_g) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { ; }
};

struct ProblemUtilities_WallisPi { emp::vector<std::function<double(TestOrg_WallisPi &)>> lexicase_fit_set; };


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: StringLengthsBackwards
// - Input: emp::vector<std::string> 
// - Output: emp::vector<size_t>
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Generate random test input
Problem_StringLengthsBackwards_input_t GenRandomTestInput_StringLengthsBackwards(emp::Random & rand,
                                                                                 const std::pair<size_t, size_t> & str_cnt_range,
                                                                                 const std::pair<size_t, size_t> & str_size_range) {

  emp::vector<char> valid_chars = {'\n', '\t'};
  for (size_t i = 32; i < 127; ++i) valid_chars.emplace_back((char)i);

  emp::vector<std::string> rand_input;
  // How many strings should we generate?
  const size_t str_cnt = rand.GetUInt(str_cnt_range.first, str_cnt_range.second);
  // Generate each string randomly.
  for (size_t i = 0; i < str_cnt; ++i) {
    const size_t str_size = rand.GetUInt(str_size_range.first, str_size_range.second);
    rand_input.emplace_back("");
    for (size_t i = 0; i < str_size; ++i) rand_input.back().push_back(valid_chars[rand.GetUInt(valid_chars.size())]);
  }
  return rand_input;
}

/// Generate correct output
Problem_StringLengthsBackwards_output_t GenCorrectOut_StringLengthsBackwards(const Problem_StringLengthsBackwards_input_t & input) {
  emp::vector<size_t> output;
  for (int i = input.size()-1; i >= 0; --i) {
    output.emplace_back(input[i].size());
  }
  return output;
}

/// String Lengths Backwards: Vector<String>
class TestOrg_StringLengthsBackwards : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = emp::vector<std::string>;
    using out_t = Problem_StringLengthsBackwards_output_t;

  protected:
    genome_t genome;
    out_t out;

  public:
    TestOrg_StringLengthsBackwards(const genome_t & _g) : genome(_g), out() { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { out = GenCorrectOut_StringLengthsBackwards(genome); }

    out_t & GetCorrectOut() { return out; }
    const out_t & GetCorrectOut() const { return out; }   

    void Print(std::ostream & os=std::cout) {
      os << "[";
      for (size_t i = 0; i < genome.size(); ++i) {
        if (i) os << ",";
        os << "{BEGIN-STR}" << genome[i] << "{END-STR}";
      }
      os << "]"; 
    }
};

struct ProblemUtilities_StringLengthsBackwards { 
  using this_t = ProblemUtilities_StringLengthsBackwards;
  using problem_org_t = TestOrg_StringLengthsBackwards;
  using input_t = Problem_StringLengthsBackwards_input_t;
  using output_t = Problem_StringLengthsBackwards_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  emp::vector<emp::Ptr<problem_org_t>> testingset_pop;

  emp::vector<emp::vector<output_t>> population_validation_outputs;

  // --- Useful during a test evaluation ---
  emp::Ptr<problem_org_t> cur_eval_test_org;
  bool submitted;
  emp::vector<size_t> submitted_vec;

  // Mutation
  size_t MIN_STR_LEN;
  size_t MAX_STR_LEN;
  size_t MIN_STR_CNT;
  size_t MAX_STR_CNT;
  double PER_STR_SWAP_RATE;
  double PER_STR_DUP_RATE;
  double PER_STR_DEL_RATE;
  double PER_CHAR_INS_RATE;
  double PER_CHAR_DEL_RATE;
  double PER_CHAR_SUB_RATE;
  emp::vector<char> valid_chars;

  size_t Mutate(emp::Random & rnd, input_t & mut_input) {
    size_t muts = 0;
    
    // Mutate each string.
    for (size_t i = 0; i < mut_input.size(); ++i) {
      std::string & str = mut_input[i];
      int expected_size = (int)str.size();
      std::string new_string = "";

      for (size_t s = 0; s < str.size(); ++s) {  
        // Per-site insertions
        if (rnd.P(PER_CHAR_INS_RATE) && (expected_size+1 <= (int)MAX_STR_LEN)) {
          // Insert a random character.
          new_string.push_back(valid_chars[rnd.GetUInt(valid_chars.size())]);
          ++expected_size;
          ++muts;
        }
        // Per-site deletions
        if (rnd.P(PER_CHAR_DEL_RATE) && (expected_size-1 >= (int)MIN_STR_LEN)) {
          --expected_size;
          ++muts;
          continue;
        }
        size_t whead = new_string.size();
        new_string.push_back(str[s]); 
        // Per-site substitutions
        if (rnd.P(PER_CHAR_SUB_RATE)) {
          new_string[whead] = valid_chars[rnd.GetUInt(valid_chars.size())];
          ++muts;
        }
      }
      emp_assert(new_string.size() >= MIN_STR_LEN);
      emp_assert(new_string.size() <= MAX_STR_LEN);
      mut_input[i] = new_string;
    }    

    // String-level mutations.
    // Swaps
    for (size_t i = 0; i < mut_input.size(); ++i) {
      if (rnd.P(PER_STR_SWAP_RATE) && mut_input.size() > 1) {
        ++muts;
        // Swap position i with a random different position
        size_t other = rnd.GetUInt(mut_input.size());
        while (other == i) other = rnd.GetUInt(mut_input.size());
        std::swap(mut_input[i], mut_input[other]);
      }
    }
    // Dups
    const size_t osize = mut_input.size();
    for (size_t i = 0; i < osize; ++i) {
      if (rnd.P(PER_STR_DUP_RATE) && mut_input.size() < MAX_STR_CNT) {
        ++muts;
        // Duplicate position i
        mut_input.emplace_back(mut_input[i]);
      }
    }
    // Deletions
    emp::vector<std::string> new_input;
    for (size_t i = 0; i < mut_input.size(); ++i) {
      if (rnd.P(PER_STR_DUP_RATE) && mut_input.size() > MIN_STR_CNT) {
        ++muts;
        continue;
      }
      new_input.emplace_back(mut_input[i]);
    }
    
    mut_input = new_input;

    emp_assert(mut_input.size() <= MAX_STR_CNT);
    emp_assert(mut_input.size() >= MIN_STR_CNT);

    return muts;
  }
  
  emp::vector<std::function<double(problem_org_t &)>> lexicase_fit_set; 

  ProblemUtilities_StringLengthsBackwards()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_vec()
  { 
    valid_chars = {'\n', '\t'};
    for (size_t i = 32; i < 127; ++i) valid_chars.emplace_back((char)i); 
  }

  ~ProblemUtilities_StringLengthsBackwards() {
    for (size_t i = 0; i < testingset_pop.size(); ++i) testingset_pop[i].Delete();
  }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  void ResetTestEval() {
    submitted = false;
    submitted_vec.clear();
  }

  void Submit(size_t val) {
    submitted = true;
    submitted_vec.emplace_back(val);
  }

  void Submit(const emp::vector<size_t> & vec) {
    submitted = true;
    submitted_vec = vec;
  }
  
  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    
    // Load input.
    // std::cout << "==== LOAD TEST CASE FROM LINE ====" << std::endl;
    // Parse line[0]
    std::string input_str = line[0];
    if (input_str.front() == '[') { input_str.erase(0, 1); }
    if (input_str.back() == ']') { input_str.pop_back(); }
    input_str += "\n";
    
    // std::cout << "Input str =" << input_str << std::endl;
    // std::cout << "Line[0] = " << line[0] << std::endl;
    // std::cout << "(1) Replace commas" << std::endl;
    input_str = StrReplace(input_str, ",", "{COMMA}");                  // Get rid of commas to avoid confusing CSV parser
    // std::cout << "(1.25) Replace escaped backslashes" << std::endl;
    input_str = StrReplace(input_str, "\\\\", "{BSLASH}");              // Get rid of backslashes so we can get rid of escaped quotes to avoid confusing the CSV parser
    // std::cout << "(1.5) Replace escaped uotes" << std::endl; 
    input_str = StrReplace(input_str, "\\\"", "{QUOTE}");               // Get rid of quotes to avoid confusing the parser
    // std::cout << "  Input str = " << input_str << std::endl; 
    
    // We use a csv parser to tackle the input line (after a bit of processing things that confuse the parser...)
    std::istringstream instr(input_str);
    aria::csv::CsvParser parser(instr);
    parser.delimiter(' ');
    parser.quote('"');
    parser.terminator('\n');

    size_t cnt = 0;
    for (auto & row : parser) {
      ++cnt;
      for (auto & field : row) {
        std::string fstr = StrReplace(field, "{COMMA}", ",");
        fstr = StrReplace(fstr, "\\t", "\t");
        fstr = StrReplace(fstr, "\\n", "\n");
        fstr = StrReplace(fstr, "{QUOTE}", "\"");
        fstr = StrReplace(fstr, "{BSLASH}", "\\");
        // std::cout << "- Instr Field=" << fstr << std::endl;
        input.emplace_back(fstr);
      }
    }

    emp::vector<size_t> correct_out;
    emp::vector<size_t> read_out;

    // What do we *expect* the correct output to be?
    for (int i = input.size()-1; i >= 0; --i) { correct_out.emplace_back(input[i].size()); }
    
    // Handle output
    if (line.size() > 1) {
      // What do we read the correct output as?
      emp::vector<std::string> sliced_line = emp::slice(line[1], '\n');
      for (size_t i = 0; i < sliced_line.size(); ++i) {
        read_out.emplace_back(std::atoi(sliced_line[i].c_str()));
      }
      // If read out and correct out are not the same size, throw an error.
      if (read_out.size() != correct_out.size()) {
        std::cout << "ERROR! Read len(" << read_out.size() << ") != gen len (" << correct_out.size() << "). Exiting" << std::endl;
        exit(-1);
      }
      // If read out and correct out are not identical, throw an error.
      if (read_out != correct_out) {
        std::cout << "ERROR! Read output different from calculated output!" << std::endl;
        
        std::cout << "Read out: ";
        for (size_t i = 0; i < read_out.size(); ++i) {
          if (i) std::cout << ",";
          std::cout << read_out[i];
        } std::cout << std::endl;

        std::cout << "Calculated out: ";
        for (size_t i = 0; i < correct_out.size(); ++i) {
          if (i) std::cout << ",";
          std::cout << correct_out[i];
        } std::cout << std::endl;
        
        std::cout << "Exiting." << std::endl;
        exit(-1);
      }

    }

    output = correct_out;

    // std::cout << "output = [";
    // for (size_t i = 0; i < output.size(); ++i) {
      // if (i) std::cout << ",";
      // std::cout << output[i];
    // } std::cout << "]" << std::endl;

    return {input, output};
  }

  void GenerateTestingSetPop() {
    for (size_t i = 0; i < testing_set.GetSize(); ++i) {
      testingset_pop.emplace_back(emp::NewPtr<problem_org_t>(testing_set.GetInput(i)));
      testingset_pop[i]->CalcOut();
    }
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  std::pair<double, bool> CalcScoreGradient(const output_t & correct_test_output, const output_t & sub) {
    const double max_dist = emp::Max(correct_test_output.size(), sub.size());
    double dist = emp::calc_edit_distance(correct_test_output, sub);
    if (dist == 0) {
      return {1.0, true};
    } else {
      return {(max_dist - dist)/max_dist, false};
    }
  } // todo - test this

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << "[";
    for (size_t i = 0; i < in.size(); ++i) {
      if (i) os << ",";
      os << "{BEGIN-STR}" << in[i] << "{END-STR}";
    }
    os << "]"; 
    os << "\"";
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: LastIndexOfZero
// - Input: emp::vector<int>; (at least one value *must* be zero)
// - Output: int;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Generate random test input
Problem_LastIndexOfZero_input_t GenRandomTestInput_LastIndexOfZero(emp::Random & rand,
                                                                          const std::pair<size_t, size_t> & vec_size_range,
                                                                          const std::pair<int, int> & vec_val_range) {
  emp::vector<int> input;
  size_t zero_cnt = 0;
  // Randomize size.
  const size_t vec_size = rand.GetUInt(vec_size_range.first, vec_size_range.second+1);
  for (size_t i = 0; i < vec_size; ++i) {
    int val = rand.GetInt(vec_val_range.first, vec_val_range.second+1);
    input.emplace_back(val);
    if (val == 0) zero_cnt++;
  }
  // ensure there's at least one zero
  if (zero_cnt == 0) {
    input[rand.GetUInt(0, input.size())] = 0;
  }
  emp_assert(emp::Has(input, 0));
  return input;
}

Problem_LastIndexOfZero_output_t GenCorrectOut_LastIndexOfZero(const Problem_LastIndexOfZero_input_t & input) {
  emp_assert(emp::Has(input, 0));
  for (int i = input.size()-1; i >= 0; --i) {
    if (input[i] == 0) return i;
  }
  std::cout << "ERROR! Failed to find 0 in last index of zero output!" << std::endl;
  exit(-1);
  return 0; // Should never get here!
}

/// Last Index of Zero: Vector<Integer>
class TestOrg_LastIndexOfZero : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = emp::vector<int>;
    using out_t = Problem_LastIndexOfZero_output_t;

  protected:
    genome_t genome;
    out_t out;

  public:
    TestOrg_LastIndexOfZero(const genome_t & _g) : genome(_g), out() { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { out = GenCorrectOut_LastIndexOfZero(genome); }

    out_t & GetCorrectOut() { return out; }
    const out_t & GetCorrectOut() const { return out; }   

    void Print(std::ostream & os=std::cout) {
      os << "[";
      for (size_t i = 0; i < genome.size(); ++i) {
        if (i) os << ",";
        os << genome[i];
      }
      os << "]"; 
    }
};

struct ProblemUtilities_LastIndexOfZero { 
  using this_t = ProblemUtilities_LastIndexOfZero;
  using problem_org_t = TestOrg_LastIndexOfZero;
  using input_t = Problem_LastIndexOfZero_input_t;
  using output_t = Problem_LastIndexOfZero_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  emp::vector<emp::Ptr<problem_org_t>> testingset_pop;

  emp::vector<emp::vector<output_t>> population_validation_outputs;

  // --- Useful during a test evaluation ---
  emp::Ptr<problem_org_t> cur_eval_test_org;
  bool submitted;
  int submitted_val;

  int MAX_ERROR;

  // Mutation
  size_t MIN_VEC_LEN;
  size_t MAX_VEC_LEN;
  int MIN_NUM;
  int MAX_NUM;
  double PER_NUM_SWAP_RATE;
  double PER_NUM_DEL_RATE;
  double PER_NUM_INS_RATE;
  double PER_NUM_SUB_RATE;

  size_t Mutate(emp::Random & rnd, input_t & mut_input) {
    size_t muts = 0; 
    emp::vector<int> new_input;
    int expected_size = mut_input.size();
    for (size_t i = 0; i < mut_input.size(); ++i) {
      // Deletion?
      if (rnd.P(PER_NUM_DEL_RATE) && (expected_size-1)>=(int)MIN_VEC_LEN) {
        ++muts;
        --expected_size;
        continue;
      } 
      const size_t whead = new_input.size();
      new_input.emplace_back(mut_input[i]);
      // Substitution?
      if (rnd.P(PER_NUM_SUB_RATE)) {
        ++muts;
        new_input[whead] = rnd.GetInt(MIN_NUM, MAX_NUM+1);
      }
      // Insertion?
      if (rnd.P(PER_NUM_INS_RATE) && (expected_size+1)<=(int)MAX_VEC_LEN) {
        ++muts;
        ++expected_size;
        new_input.emplace_back(rnd.GetInt(MIN_NUM, MAX_NUM+1));
      }
    }
    mut_input = new_input; // update mut_input
    // Swaps?
    for (size_t i = 0; i < mut_input.size(); ++i) {
      if (rnd.P(PER_NUM_SWAP_RATE) && mut_input.size() > 1) {
        ++muts;
        // Swap position i with a random different position
        size_t other = rnd.GetUInt(mut_input.size());
        while (other == i) other = rnd.GetUInt(mut_input.size());
        std::swap(mut_input[i], mut_input[other]);
      }
    }
    // Guarantee that vector still has at least one zero!
    if (!emp::Has(mut_input, 0)) {
      mut_input[rnd.GetUInt(mut_input.size())] = 0;
    }
    emp_assert(emp::Has(mut_input, 0));
    return muts;
  }
  
  emp::vector<std::function<double(problem_org_t &)>> lexicase_fit_set; 

  ProblemUtilities_LastIndexOfZero()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val(0)
  { ; }

  ~ProblemUtilities_LastIndexOfZero() {
    for (size_t i = 0; i < testingset_pop.size(); ++i) testingset_pop[i].Delete();
  }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  void ResetTestEval() {
    submitted = false;
    submitted_val = 0;
  }

  void Submit(int val) {
    submitted = true;
    submitted_val = val;
  }


  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    
    std::string input_str = line[0];
    if (input_str.front() == '[') { input_str.erase(0, 1); }
    if (input_str.back() == ']') { input_str.pop_back(); }
    emp::vector<std::string> sliced_input_str = emp::slice(input_str, ' ');
    
    for (size_t i = 0; i < sliced_input_str.size(); ++i) {
      input.emplace_back(std::atoi(sliced_input_str[i].c_str()));
    }
    
    // Calculate correct output given loaded input
    int calc_out = GenCorrectOut_LastIndexOfZero(input);
    // Get output from file
    output = std::atoi(line[1].c_str());
    // make sure generated output and read output match
    if (calc_out != output) {
      std::cout << "ERROR! Generated last index of zero output does not match read output! Exiting." << std::endl;
      exit(-1);
    }

    return {input, output};
  }

  void GenerateTestingSetPop() {
    for (size_t i = 0; i < testing_set.GetSize(); ++i) {
      testingset_pop.emplace_back(emp::NewPtr<problem_org_t>(testing_set.GetInput(i)));
      testingset_pop[i]->CalcOut();
    }
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  std::pair<double, bool> CalcScoreGradient(const output_t & correct_test_output, const output_t & sub) {
    if (correct_test_output == sub) {
      return {1.0, true};
    } else {
      double error = (double)emp::Abs(correct_test_output - sub);
      emp_assert(error != 0, "Error shouldn't be zero here!");
      double score = (error <= MAX_ERROR) ? 1 - (error/(double)MAX_ERROR) : 0.0;
      return {score, false};
    }
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << "[";
    for (size_t i = 0; i < in.size(); ++i) {
      if (i) os << ",";
      os << in[i];
    }
    os << "]"; 
    os << "\"";
  }

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: VectorAverage
// - Input: emp::vector<double>
// - Output: double (must be within epsilon)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Generate random test input
Problem_VectorAverage_input_t GenRandomTestInput_VectorAverage(emp::Random & rand,
                                                              const std::pair<size_t, size_t> & vec_size_range,
                                                              const std::pair<double, double> & vec_val_range) {
  emp::vector<double> input;
  const size_t vec_size = rand.GetUInt(vec_size_range.first, vec_size_range.second);
  for (size_t i = 0; i < vec_size; ++i) {
    double val = rand.GetDouble(vec_val_range.first, vec_val_range.second+1);
    if (val > vec_val_range.second) val = vec_val_range.second;
    input.emplace_back(val);
  }
  return input;
}

/// Generate correct test output
Problem_VectorAverage_output_t GenCorrectOut_VectorAverage(const Problem_VectorAverage_input_t & input) {
  emp_assert(input.size() > 0);
  if (input.size()) return (double)emp::Sum(input) / (double)input.size();
  else return 0;
}

/// Vector Average: Vector<Float>
class TestOrg_VectorAverage : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = emp::vector<double>;
    using out_t = Problem_VectorAverage_output_t;

  protected:
    genome_t genome;
    out_t out;

  public:
    TestOrg_VectorAverage(const genome_t & _g) : genome(_g), out() { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { out = GenCorrectOut_VectorAverage(genome); }

    out_t & GetCorrectOut() { return out; }
    const out_t & GetCorrectOut() const { return out; }   

    void Print(std::ostream & os=std::cout) {
      os << "[";
      for (size_t i = 0; i < genome.size(); ++i) {
        if (i) os << ",";
        os << genome[i];
      }
      os << "]"; 
    }
};

struct ProblemUtilities_VectorAverage {
  using this_t = ProblemUtilities_VectorAverage;
  using problem_org_t = TestOrg_VectorAverage;
  using input_t = Problem_VectorAverage_input_t;
  using output_t = Problem_VectorAverage_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  emp::vector<emp::Ptr<problem_org_t>> testingset_pop;

  emp::vector<emp::vector<output_t>> population_validation_outputs;

  // --- Useful during a test evaluation ---
  emp::Ptr<problem_org_t> cur_eval_test_org;
  bool submitted;
  double submitted_val;
  double MAX_ERROR;

  double EPSILON; // How much error do we allow submitted values to have to still be correct?

  // Mutation
  size_t MIN_VEC_LEN;
  size_t MAX_VEC_LEN;
  double MIN_NUM;
  double MAX_NUM;
  double INS_RATE;
  double DEL_RATE;
  double SUB_RATE;
  
  size_t Mutate(emp::Random & rnd, input_t & mut_input) {
    size_t muts = 0; 
    emp::vector<double> new_input;
    int expected_size = mut_input.size();
    for (size_t i = 0; i < mut_input.size(); ++i) {
      // Deletion?
      if (rnd.P(DEL_RATE) && (expected_size-1)>=(int)MIN_VEC_LEN) {
        ++muts;
        --expected_size;
        continue;
      } 
      const size_t whead = new_input.size();
      new_input.emplace_back(mut_input[i]);
      // Substitution?
      if (rnd.P(SUB_RATE)) {
        ++muts;
        double val = rnd.GetDouble(MIN_NUM, MAX_NUM+1);
        if (val > MAX_NUM) val = MAX_NUM;
        new_input[whead] = val;
      }
      // Insertion?
      if (rnd.P(INS_RATE) && (expected_size+1)<=(int)MAX_VEC_LEN) {
        ++muts;
        ++expected_size;
        double val = rnd.GetDouble(MIN_NUM, MAX_NUM+1);
        if (val > MAX_NUM) val = MAX_NUM;
        new_input.emplace_back(val);
      }
    }
    mut_input = new_input; // update mut_input
    return muts;
  }
  
  emp::vector<std::function<double(problem_org_t &)>> lexicase_fit_set; 

  ProblemUtilities_VectorAverage()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val(0.0)
  { ; }

  ~ProblemUtilities_VectorAverage() {
    for (size_t i = 0; i < testingset_pop.size(); ++i) testingset_pop[i].Delete();
  }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  void ResetTestEval() {
    submitted = false;
    submitted_val = 0.0;
  }

  void Submit(double val) {
    submitted = true;
    submitted_val = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    
    std::string input_str = line[0];
    if (input_str.front() == '[') { input_str.erase(0, 1); }
    if (input_str.back() == ']') { input_str.pop_back(); }
    emp::vector<std::string> sliced_input_str = emp::slice(input_str, ' ');
    
    for (size_t i = 0; i < sliced_input_str.size(); ++i) {
      input.emplace_back(std::atof(sliced_input_str[i].c_str()));
    }
    
    // Calculate correct output given loaded input
    double calc_out = GenCorrectOut_VectorAverage(input);
    // Get output from file
    output = std::atof(line[1].c_str());
    // make sure generated output and read output match
    if (! ((output <= (calc_out + 0.00005)) && (output >= (calc_out - 0.00005))) ) {
      std::cout << "ERROR! Generated output does not match read output! Exiting." << std::endl;
      exit(-1);
    }

    return {input, output};
  }

  bool EpsilonEqu(double target, double val) {
    if (val <= target + EPSILON && val >= target - EPSILON) return true;
    else return false;
  }

  void GenerateTestingSetPop() {
    for (size_t i = 0; i < testing_set.GetSize(); ++i) {
      testingset_pop.emplace_back(emp::NewPtr<problem_org_t>(testing_set.GetInput(i)));
      testingset_pop[i]->CalcOut();
    }
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = EpsilonEqu(correct_test_output, sub);
    return {(double)pass, pass};
  }

  std::pair<double, bool> CalcScoreGradient(const output_t & correct_test_output, const output_t & sub) {
    if (EpsilonEqu(correct_test_output, sub)) {
      return {1.0, true};
    } else {
      double error = (double)emp::Abs(correct_test_output - sub);
      emp_assert(error != 0, "Error shouldn't be zero here!");
      double score = (error <= MAX_ERROR) ? 1 - (error/(double)MAX_ERROR) : 0.0;
      return {score, false};
    }
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << "[";
    for (size_t i = 0; i < in.size(); ++i) {
      if (i) os << ",";
      os << in[i];
    }
    os << "]"; 
    os << "\"";
  }

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: CountOdds
// - Input: emp::vector<int>
// - Output: int
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Generate random test input
Problem_CountOdds_input_t GenRandomTestInput_CountOdds(emp::Random & rand,
                                                       const std::pair<size_t, size_t> & vec_size_range,
                                                       const std::pair<int, int> & vec_val_range) {
  emp::vector<int> input;
  const size_t vec_size = rand.GetUInt(vec_size_range.first, vec_size_range.second+1);
  for (size_t i = 0; i < vec_size; ++i) {
    int val = rand.GetInt(vec_val_range.first, vec_val_range.second+1);
    input.emplace_back(val);
  }
  return input;
}

/// Generate correct test output
Problem_CountOdds_output_t GenCorrectOut_CountOdds(const Problem_CountOdds_input_t & input) {
  size_t odd_cnt = 0;
  for (size_t i = 0; i < input.size(); ++i) {
    if (input[i] & 1) ++odd_cnt;
  }
  return odd_cnt;
}

class TestOrg_CountOdds : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = emp::vector<int>;
    using out_t = Problem_CountOdds_output_t;

  protected:
    genome_t genome;
    out_t out;

  public:
    TestOrg_CountOdds(const genome_t & _g) : genome(_g), out() { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { out = GenCorrectOut_CountOdds(genome); }

    out_t & GetCorrectOut() { return out; }
    const out_t & GetCorrectOut() const { return out; }   

    void Print(std::ostream & os=std::cout) {
      os << "[";
      for (size_t i = 0; i < genome.size(); ++i) {
        if (i) os << ",";
        os << genome[i];
      }
      os << "]"; 
    }
};

struct ProblemUtilities_CountOdds { 
  using this_t = ProblemUtilities_CountOdds;
  using problem_org_t = TestOrg_CountOdds;
  using input_t = Problem_CountOdds_input_t;
  using output_t = Problem_CountOdds_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  emp::vector<emp::Ptr<problem_org_t>> testingset_pop;

  emp::vector<emp::vector<output_t>> population_validation_outputs;

  // --- Useful during a test evaluation ---
  emp::Ptr<problem_org_t> cur_eval_test_org;
  bool submitted;
  int submitted_val;
  int MAX_ERROR;

  // Mutation
  size_t MIN_VEC_LEN;
  size_t MAX_VEC_LEN;
  int MIN_NUM;
  int MAX_NUM;
  double PER_NUM_SWAP_RATE;
  double PER_NUM_DEL_RATE;
  double PER_NUM_INS_RATE;
  double PER_NUM_SUB_RATE;

  size_t Mutate(emp::Random & rnd, input_t & mut_input) {
    size_t muts = 0; 
    emp::vector<int> new_input;
    int expected_size = mut_input.size();
    for (size_t i = 0; i < mut_input.size(); ++i) {
      // Deletion?
      if (rnd.P(PER_NUM_DEL_RATE) && (expected_size-1)>=(int)MIN_VEC_LEN) {
        ++muts;
        --expected_size;
        continue;
      } 
      const size_t whead = new_input.size();
      new_input.emplace_back(mut_input[i]);
      // Substitution?
      if (rnd.P(PER_NUM_SUB_RATE)) {
        ++muts;
        new_input[whead] = rnd.GetInt(MIN_NUM, MAX_NUM+1);
      }
      // Insertion?
      if (rnd.P(PER_NUM_INS_RATE) && (expected_size+1)<=(int)MAX_VEC_LEN) {
        ++muts;
        ++expected_size;
        new_input.emplace_back(rnd.GetInt(MIN_NUM, MAX_NUM+1));
      }
    }
    mut_input = new_input; // update mut_input
    // Swaps?
    for (size_t i = 0; i < mut_input.size(); ++i) {
      if (rnd.P(PER_NUM_SWAP_RATE) && mut_input.size() > 1) {
        ++muts;
        // Swap position i with a random different position
        size_t other = rnd.GetUInt(mut_input.size());
        while (other == i) other = rnd.GetUInt(mut_input.size());
        std::swap(mut_input[i], mut_input[other]);
      }
    }
    return muts;
  }
  
  emp::vector<std::function<double(problem_org_t &)>> lexicase_fit_set; 

  ProblemUtilities_CountOdds()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val(0)
  { ; }

  ~ProblemUtilities_CountOdds() {
    for (size_t i = 0; i < testingset_pop.size(); ++i) testingset_pop[i].Delete();
  }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  void ResetTestEval() {
    submitted = false;
    submitted_val = 0;
  }

  void Submit(int val) {
    submitted = true;
    submitted_val = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    
    std::string input_str = line[0];
    if (input_str.front() == '[') { input_str.erase(0, 1); }
    if (input_str.back() == ']') { input_str.pop_back(); }
    emp::vector<std::string> sliced_input_str = emp::slice(input_str, ' ');
    
    for (size_t i = 0; i < sliced_input_str.size(); ++i) {
      input.emplace_back(std::atoi(sliced_input_str[i].c_str()));
    }
    
    // Calculate correct output given loaded input
    int calc_out = GenCorrectOut_CountOdds(input);
    // Get output from file
    output = std::atoi(line[1].c_str());
    // make sure generated output and read output match
    if (calc_out != output) {
      std::cout << "ERROR! Generated count odds output does not match read output! Exiting." << std::endl;
      exit(-1);
    }

    return {input, output};
  }

  void GenerateTestingSetPop() {
    for (size_t i = 0; i < testing_set.GetSize(); ++i) {
      testingset_pop.emplace_back(emp::NewPtr<problem_org_t>(testing_set.GetInput(i)));
      testingset_pop[i]->CalcOut();
    }
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  std::pair<double, bool> CalcScoreGradient(const output_t & correct_test_output, const output_t & sub) {
    if (correct_test_output == sub) {
      return {1.0, true};
    } else {
      double error = (double)emp::Abs(correct_test_output - sub);
      emp_assert(error != 0, "Error shouldn't be zero here!");
      double score = (error <= MAX_ERROR) ? 1 - (error/(double)MAX_ERROR) : 0.0;
      return {score, false};
    }
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << "[";
    for (size_t i = 0; i < in.size(); ++i) {
      if (i) os << ",";
      os << in[i];
    }
    os << "]"; 
    os << "\"";
  }

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: MirrorImage
// - Input: std::array<emp::vector<int>, 2>;
// - Output: bool;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Generate random input
/// Generate random test input
Problem_MirrorImage_input_t GenRandomTestInput_MirrorImage(emp::Random & rand,
                                                           const std::pair<size_t, size_t> & vec_size_range,
                                                           const std::pair<int, int> & vec_val_range) {
  std::array<emp::vector<int>, 2> input;
  // 4 random cases: equal, random, mirrored, nearly mirrored
  const size_t vec_size = rand.GetUInt(vec_size_range.first, vec_size_range.second+1);
  const size_t rand_case = rand.GetUInt(0, 4);
  switch (rand_case) {
    case 0: { // Equal
      for (size_t k = 0; k < vec_size; ++k) {
        input[0].emplace_back(rand.GetInt(vec_val_range.first, vec_val_range.second+1));
      }
      input[1] = input[0];
      break; 
    }
    case 1: { // Random
      for (size_t i = 0; i < input.size(); ++i) {
        for (size_t k = 0; k < vec_size; ++k) {
          input[i].emplace_back(rand.GetInt(vec_val_range.first, vec_val_range.second+1));
        }
      } 
      break; 
    }
    case 2: { // Mirrored
      for (size_t k = 0; k < vec_size; ++k) {
        input[0].emplace_back(rand.GetInt(vec_val_range.first, vec_val_range.second+1));
      }
      input[1] = input[0];
      std::reverse(std::begin(input[1]), std::end(input[1]));
      break; 
    }
    case 3: { // Nearly mirrored
      for (size_t k = 0; k < vec_size; ++k) {
        input[0].emplace_back(rand.GetInt(vec_val_range.first, vec_val_range.second+1));
      }
      input[1] = input[0];
      std::reverse(std::begin(input[1]), std::end(input[1]));
      const size_t num_randos = rand.GetUInt(input[0].size());
      for (size_t i = 0; i < num_randos; ++i) {
        input[0][rand.GetUInt(input[0].size())] = rand.GetInt(vec_val_range.first, vec_val_range.second+1);
      }
      break; 
    }
  }
  return input;
}

/// Generate correct output given input
Problem_MirrorImage_output_t GenCorrectOut_MirrorImage(const Problem_MirrorImage_input_t & input) {
  emp_assert(input[0].size() == input[1].size());
  if (input[0].size() != input[1].size()) return false;
  const size_t vec_size = input[0].size();
  for (size_t i = 0; i < vec_size; ++i) {
    const size_t ri = (vec_size - 1) - i;
    if (input[0][i] != input[1][ri]) return false;
  }
  return true;
}

/// Mirror Image: Array<Vector<Integer>, 2>
class TestOrg_MirrorImage : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = std::array<emp::vector<int>, 2>;
    using out_t = Problem_MirrorImage_output_t;

  protected:
    genome_t genome;
    out_t out;

  public:
    TestOrg_MirrorImage(const genome_t & _g) : genome(_g) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { out = GenCorrectOut_MirrorImage(genome); }

    out_t & GetCorrectOut() { return out; }
    const out_t & GetCorrectOut() const { return out; }   

    void Print(std::ostream & os=std::cout) {
      os << "[[";
      for (size_t i = 0; i < genome[0].size(); ++i) {
        if (i) os << ",";
        os << genome[0][i];
      }
      os << "],["; 
      for (size_t i = 0; i < genome[1].size(); ++i) {
        if (i) os << ",";
        os << genome[1][i];
      }
      os << "]]"; 
    }
};

struct ProblemUtilities_MirrorImage {
  using this_t = ProblemUtilities_MirrorImage;
  using problem_org_t = TestOrg_MirrorImage;
  using input_t = Problem_MirrorImage_input_t;
  using output_t = Problem_MirrorImage_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  emp::vector<emp::Ptr<problem_org_t>> testingset_pop;

  emp::vector<emp::vector<output_t>> population_validation_outputs;

  // --- Useful during a test evaluation ---
  emp::Ptr<problem_org_t> cur_eval_test_org;
  bool submitted;
  bool submitted_val;

  // Mutation
  size_t MIN_VEC_LEN;
  size_t MAX_VEC_LEN;
  int MIN_NUM;
  int MAX_NUM;
  double PER_VEC_RANDOMIZE_VAL_RATE;  // How often do we randomize a single value in a vector?
  double PER_VEC_MIRROR_RATE;              // How often do we mirror a vector?
  double COPY_RATE;                   // How often do we copy one vector over the other?
  double INS_RATE;                    // How often do we insert a random value on the end of both vectors?
  double DEL_RATE;                    // How often do we delete
  double PER_VEC_SHUFFLE_RATE;

  size_t Mutate(emp::Random & rnd, input_t & mut_input) {
    size_t muts = 0; 
    emp_assert(mut_input[0].size() == mut_input[1].size());
    if (mut_input[0].size() != mut_input[1].size()) {
      std::cout << "ERROR! MirrorImage vectors do not have same length! Exiting." << std::endl;
      exit(-1);
    }
    size_t vec_size = mut_input[0].size();

    // Should we randomize any values?
    for (size_t i = 0; i < mut_input.size(); ++i) {
      if (rnd.P(PER_VEC_RANDOMIZE_VAL_RATE) && mut_input[i].size() > 0) {
        ++muts;
        mut_input[i][rnd.GetUInt(mut_input[i].size())] = rnd.GetInt(MIN_NUM, MAX_NUM+1);
      }
    }
    // Should we copy one vector over the other?
    if (rnd.P(COPY_RATE)) {
      ++muts;
      // Which should we copy?
      if (rnd.P(0.5)) {
        mut_input[1] = mut_input[0];
      } else {
        mut_input[0] = mut_input[1];
      }
    }
    
    // Should we insert anything into each vector? (has weird side effects, but don't really care in this case)
    if (rnd.P(INS_RATE) && vec_size < MAX_VEC_LEN) {
      ++muts;
      // Where?
      const size_t loc = rnd.GetUInt(vec_size);
      const int val = rnd.GetInt(MIN_NUM, MAX_NUM+1);
      mut_input[0].emplace_back(val); std::swap(mut_input[0][loc], mut_input[0][vec_size]);
      mut_input[1].emplace_back(val); std::swap(mut_input[1][loc], mut_input[1][vec_size]);
      vec_size = mut_input[0].size();
    }

    // Should we delete anything in each vector? (has weird side effects, but don't really care in this case)
    if (rnd.P(DEL_RATE) && vec_size > MIN_VEC_LEN) {
      ++muts;
      const size_t loc = rnd.GetUInt(vec_size);
      --vec_size;
      std::swap(mut_input[0][loc], mut_input[0][vec_size]); mut_input[0].resize(vec_size);
      std::swap(mut_input[1][loc], mut_input[1][vec_size]); mut_input[1].resize(vec_size);
    } 

    // Should we mirror any of the vectors?
    for (size_t i = 0; i < mut_input.size(); ++i) {
      if (rnd.P(PER_VEC_MIRROR_RATE)) {
        ++muts;
        std::reverse(std::begin(mut_input[i]), std::end(mut_input[i]));
      }
    }

    // Should we randomize a vector?
    for (size_t i = 0; i < mut_input.size(); ++i) {
      if (rnd.P(PER_VEC_SHUFFLE_RATE)) {
        ++muts;
        emp::Shuffle(rnd, mut_input[i]);
      }
    }

    emp_assert(mut_input[0].size() == mut_input[1].size());
    emp_assert(mut_input[0].size() <= MAX_VEC_LEN);
    emp_assert(mut_input[0].size() >= MIN_VEC_LEN);
    return muts;
  }
  
  emp::vector<std::function<double(problem_org_t &)>> lexicase_fit_set; 

  ProblemUtilities_MirrorImage()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val(0)
  { ; }

  ~ProblemUtilities_MirrorImage() {
    for (size_t i = 0; i < testingset_pop.size(); ++i) testingset_pop[i].Delete();
  }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  void ResetTestEval() {
    submitted = false;
    submitted_val = false;
  }

  void Submit(bool val) {
    submitted = true;
    submitted_val = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    
    // Vector 1
    std::string input_str = line[0];
    if (input_str.front() == '[') { input_str.erase(0, 1); }
    if (input_str.back() == ']') { input_str.pop_back(); }
    emp::vector<std::string> sliced_input_str = emp::slice(input_str, ' ');
    for (size_t i = 0; i < sliced_input_str.size(); ++i) {
      input[0].emplace_back(std::atoi(sliced_input_str[i].c_str()));
    }

    // Vector 2
    input_str = line[1];
    if (input_str.front() == '[') { input_str.erase(0, 1); }
    if (input_str.back() == ']') { input_str.pop_back(); }  
    sliced_input_str = emp::slice(input_str, ' ');
    for (size_t i = 0; i < sliced_input_str.size(); ++i) {
      input[1].emplace_back(std::atoi(sliced_input_str[i].c_str()));
    }   

    // Calculate correct output given loaded input
    bool calc_out = GenCorrectOut_MirrorImage(input);
    // Get output from file
    if (line[2] == "false") output = false;
    else if (line[2] == "true") output = true;
    else {
      std::cout << "Unrecognized output value for mirror image examples! Exiting." << std::endl;
      exit(-1);
    }

    // make sure generated output and read output match
    if (calc_out != output) {
      std::cout << "ERROR! Generated output does not match read output! Exiting." << std::endl;
      exit(-1);
    }

    return {input, output};
  }

  void GenerateTestingSetPop() {
    for (size_t i = 0; i < testing_set.GetSize(); ++i) {
      testingset_pop.emplace_back(emp::NewPtr<problem_org_t>(testing_set.GetInput(i)));
      testingset_pop[i]->CalcOut();
    }
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << "[[";
    for (size_t i = 0; i < in[0].size(); ++i) {
      if (i) os << ",";
      os << in[0][i];
    }
    os << "],["; 
    for (size_t i = 0; i < in[1].size(); ++i) {
      if (i) os << ",";
      os << in[1][i];
    }
    os << "]]"; 
    os << "\"";
  }

};



/// Super Anagrams: Array<String, 2>
class TestOrg_SuperAnagrams : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = std::array<std::string, 2>;
  protected:
    genome_t genome;

  public:
    TestOrg_SuperAnagrams(const genome_t & _g) : genome(_g) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { ; }
};

struct ProblemUtilities_SuperAnagrams { emp::vector<std::function<double(TestOrg_SuperAnagrams &)>> lexicase_fit_set; };

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: SumOfSquares (NOTE: only for: 1 >= n <= 100)
// - Input: int
// - Output: int
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::array<int, 101> SumOfSquaresLookup = {0, 1, 5, 14, 30, 55, 91, 140, 204, 285, 385, 506, 650, 819, 1015, 1240, 1496, 1785, 2109, 2470, 2870, 3311, 3795, 4324, 4900, 5525, 6201, 6930, 7714, 8555, 9455, 10416, 11440, 12529, 13685, 14910, 16206, 17575, 19019, 20540, 22140, 23821, 25585, 27434, 29370, 31395, 33511, 35720, 38024, 40425, 42925, 45526, 48230, 51039, 53955, 56980, 60116, 63365, 66729, 70210, 73810, 77531, 81375, 85344, 89440, 93665, 98021, 102510, 107134, 111895, 116795, 121836, 127020, 132349, 137825, 143450, 149226, 155155, 161239, 167480, 173880, 180441, 187165, 194054, 201110, 208335, 215731, 223300, 231044, 238965, 247065, 255346, 263810, 272459, 281295, 290320, 299536, 308945, 318549, 328350, 338350};

Problem_SumOfSquares_input_t GenRandomTestInput_SumOfSquares(emp::Random & rnd, const std::pair<int, int> & num_range) {
  return rnd.GetInt(num_range.first, num_range.second+1);
}


Problem_SumOfSquares_output_t GenCorrectOut_SumOfSquares(const Problem_SumOfSquares_input_t & input) {
  if (input < (int)SumOfSquaresLookup.size()) {
    return SumOfSquaresLookup[(size_t)input];
  } else {
    return ((input) * (input + 1) * ((2*input)+1))/6;
  }
}

/// Sum of Squares: Integer
class TestOrg_SumOfSquares : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = int;
    using out_t = Problem_SumOfSquares_output_t; 

  protected:
    genome_t genome;
    out_t out;

  public:
    TestOrg_SumOfSquares(const genome_t & _g) : genome(_g) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { out = GenCorrectOut_SumOfSquares(genome); }

    out_t & GetCorrectOut() { return out; }
    const out_t & GetCorrectOut() const { return out; }  

    void Print(std::ostream & os=std::cout) {
      os << genome;
    }
};

struct ProblemUtilities_SumOfSquares { 
  using this_t = ProblemUtilities_SumOfSquares;
  using problem_org_t = TestOrg_SumOfSquares;
  using input_t = Problem_SumOfSquares_input_t;
  using output_t = Problem_SumOfSquares_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  emp::vector<emp::Ptr<problem_org_t>> testingset_pop;
  emp::vector<emp::vector<output_t>> population_validation_outputs;

  // // --- Useful during a test evaluation ---
  emp::Ptr<problem_org_t> cur_eval_test_org;
  bool submitted;
  int submitted_val;
  int MAX_ERROR;

  // // Mutation - Handle here...
  int MIN_NUM;
  int MAX_NUM;
  double NUM_MUT_RATE;
  
  size_t Mutate(emp::Random & rnd, input_t & mut_input) {
    size_t muts = 0;
    emp_assert(MIN_NUM >= 0);

    if (rnd.P(NUM_MUT_RATE)) {
      ++muts;
      mut_input = rnd.GetInt(MIN_NUM, MAX_NUM+1);
    }

    return muts;
  }

  // Selection
  emp::vector<std::function<double(problem_org_t &)>> lexicase_fit_set;

  ProblemUtilities_SumOfSquares()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val(0)
  { ; }

  ~ProblemUtilities_SumOfSquares() {
    for (size_t i = 0; i < testingset_pop.size(); ++i) testingset_pop[i].Delete();
  }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  void ResetTestEval() {
    submitted = false;
    submitted_val = 0;
  }

  void Submit(int val) {
    submitted = true;
    submitted_val = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    // Load input.
    input = std::atof(line[0].c_str());
    // Load output.
    output = std::atof(line[1].c_str());
    emp_assert(output == GenCorrectOut_SumOfSquares(input));
    return {input, output};
  }

  void GenerateTestingSetPop() {
    for (size_t i = 0; i < testing_set.GetSize(); ++i) {
      testingset_pop.emplace_back(emp::NewPtr<problem_org_t>(testing_set.GetInput(i)));
      testingset_pop[i]->CalcOut();
    }
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  std::pair<double, bool> CalcScoreGradient(const output_t & correct_test_output, const output_t & sub) {
    if (correct_test_output == sub) {
      return {1.0, true};
    } else {
      double error = (double)emp::Abs(correct_test_output - sub);
      emp_assert(error != 0, "Error shouldn't be zero here!");
      double score = (error <= MAX_ERROR) ? 1 - (error/(double)MAX_ERROR) : 0.0;
      return {score, false};
    }
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << in;
    os << "\"";
  }

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: VectorsSummed
// - Input: std::array<emp::vector<int>, 2>;
// - Output: emp::vector<int>;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Generate random test input
Problem_VectorsSummed_input_t GenRandomTestInput_VectorsSummed(emp::Random & rand,
                                                               const std::pair<size_t, size_t> & vec_size_range,
                                                                const std::pair<int, int> & vec_val_range) {
  std::array<emp::vector<int>, 2> input; // Input
  const size_t vec_size = rand.GetUInt(vec_size_range.first, vec_size_range.second+1);
  for (size_t i = 0; i < vec_size; ++i) {
    input[0].emplace_back(rand.GetInt(vec_val_range.first, vec_val_range.second+1));
    input[1].emplace_back(rand.GetInt(vec_val_range.first, vec_val_range.second+1));
  }
  return input;
}

/// Generate correct output given input
Problem_VectorsSummed_output_t GenCorrectOut_VectorsSummed(const Problem_VectorsSummed_input_t & input) {
  emp::vector<int> output;
  emp_assert(input[0].size() == input[1].size());
  const size_t vec_size = input[0].size();
  for (size_t i = 0; i < vec_size; ++i) {
    const int res = input[0][i] + input[1][i];
    output.emplace_back(res);
  }
  emp_assert(output.size() == input[0].size());
  return output;
}

/// Vectors Summed: Array<Vector<Integer>, 2>
class TestOrg_VectorsSummed: public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = std::array<emp::vector<int>, 2>;
    using out_t = Problem_VectorsSummed_output_t;

  protected:
    genome_t genome;
    out_t out;

  public:
    TestOrg_VectorsSummed(const genome_t & _g) : genome(_g) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { out = GenCorrectOut_VectorsSummed(genome); }

    out_t & GetCorrectOut() { return out; }
    const out_t & GetCorrectOut() const { return out; }   

    void Print(std::ostream & os=std::cout) {
      os << "[[";
      for (size_t i = 0; i < genome[0].size(); ++i) {
        if (i) os << ",";
        os << genome[0][i];
      }
      os << "],["; 
      for (size_t i = 0; i < genome[1].size(); ++i) {
        if (i) os << ",";
        os << genome[1][i];
      }
      os << "]]"; 
    }
};

struct ProblemUtilities_VectorsSummed { 
  using this_t = ProblemUtilities_VectorsSummed;
  using problem_org_t = TestOrg_VectorsSummed;
  using input_t = Problem_VectorsSummed_input_t;
  using output_t = Problem_VectorsSummed_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  emp::vector<emp::Ptr<problem_org_t>> testingset_pop;
  emp::vector<emp::vector<output_t>> population_validation_outputs;

  // --- Useful during a test evaluation ---
  emp::Ptr<problem_org_t> cur_eval_test_org;
  bool submitted;
  emp::vector<int> submitted_vec;
  int MAX_ERROR;

  // Mutation
  size_t MIN_VEC_LEN;
  size_t MAX_VEC_LEN;
  int MIN_NUM;
  int MAX_NUM;
  double PER_NUM_SUB_RATE;
  double COPY_RATE;                   // How often do we copy one vector over the other?
  double INS_RATE;                    // How often do we insert a random values on the end of both vectors?
  double DEL_RATE;                    // How often do we delete things from each vector
  
  size_t Mutate(emp::Random & rnd, input_t & mut_input) {
    size_t muts = 0; 
    emp_assert(mut_input[0].size() == mut_input[1].size());
    size_t vec_size = mut_input[0].size();
    
    // Should we copy one vector over the other?
    if (rnd.P(COPY_RATE)) {
      ++muts;
      if (rnd.P(0.5)) {
        mut_input[0] = mut_input[1];
      } else {
        mut_input[1] = mut_input[0];
      }
    }

    // Should we do an insertion?
    if (rnd.P(INS_RATE) && vec_size < MAX_VEC_LEN) {
      ++muts;
      mut_input[0].emplace_back(rnd.GetInt(MIN_NUM, MAX_NUM+1));
      mut_input[1].emplace_back(rnd.GetInt(MIN_NUM, MAX_NUM+1));
      ++vec_size;
    }

    // Should we do a deletion?
    if (rnd.P(DEL_RATE) && vec_size > MIN_VEC_LEN) {
      ++muts;
      --vec_size;
      mut_input[0].resize(vec_size);
      mut_input[1].resize(vec_size);
    }
    
    emp_assert(vec_size == mut_input[0].size());
    
    // Per-number substitutions
    for (size_t i = 0; i < vec_size; ++i) {
      if (rnd.P(PER_NUM_SUB_RATE)) {
        ++muts;
        mut_input[0][i] = rnd.GetInt(MIN_NUM, MAX_NUM+1);
      }
      if (rnd.P(PER_NUM_SUB_RATE)) {
        ++muts;
        mut_input[1][i] = rnd.GetInt(MIN_NUM, MAX_NUM+1);
      }
    }
    
    emp_assert(mut_input[0].size() == mut_input[1].size());
    emp_assert(mut_input[0].size() <= MAX_VEC_LEN);
    emp_assert(mut_input[0].size() >= MIN_VEC_LEN);
    return muts;
  }
  
  emp::vector<std::function<double(problem_org_t &)>> lexicase_fit_set; 

  ProblemUtilities_VectorsSummed()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_vec(0)
  { ; }

  ~ProblemUtilities_VectorsSummed() {
    for (size_t i = 0; i < testingset_pop.size(); ++i) testingset_pop[i].Delete();
  }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  void ResetTestEval() {
    submitted = false;
    submitted_vec.clear();
  }

  void Submit(const emp::vector<int> & vec) {
    submitted = true;
    submitted_vec = vec;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    
    // Vector 1
    std::string input_str = line[0];
    if (input_str.front() == '[') { input_str.erase(0, 1); }
    if (input_str.back() == ']') { input_str.pop_back(); }
    emp::vector<std::string> sliced_input_str = emp::slice(input_str, ' ');
    for (size_t i = 0; i < sliced_input_str.size(); ++i) {
      input[0].emplace_back(std::atoi(sliced_input_str[i].c_str()));
    }

    // Vector 2
    input_str = line[1];
    if (input_str.front() == '[') { input_str.erase(0, 1); }
    if (input_str.back() == ']') { input_str.pop_back(); }  
    sliced_input_str = emp::slice(input_str, ' ');
    for (size_t i = 0; i < sliced_input_str.size(); ++i) {
      input[1].emplace_back(std::atoi(sliced_input_str[i].c_str()));
    }   

    // Output vector
    input_str = line[2];
    if (input_str.front() == '[') { input_str.erase(0, 1); }
    if (input_str.back() == ']') { input_str.pop_back(); }  
    sliced_input_str = emp::slice(input_str, ' ');
    for (size_t i = 0; i < sliced_input_str.size(); ++i) {
      output.emplace_back(std::atoi(sliced_input_str[i].c_str()));
    }

    // Calculate correct output given loaded input
    emp::vector<int> calc_out = GenCorrectOut_VectorsSummed(input);
    
    // make sure generated output and read output match
    if (calc_out != output) {
      std::cout << "ERROR! Generated output does not match read output! Exiting." << std::endl;
      exit(-1);
    }

    return {input, output};
  }

  void GenerateTestingSetPop() {
    for (size_t i = 0; i < testing_set.GetSize(); ++i) {
      testingset_pop.emplace_back(emp::NewPtr<problem_org_t>(testing_set.GetInput(i)));
      testingset_pop[i]->CalcOut();
    }
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  std::pair<double, bool> CalcScoreGradient(const output_t & correct_test_output, const output_t & sub) {
    if (correct_test_output == sub) {
      return {1.0, true};
    } else {
      // double error = (double)emp::Abs(correct_test_output - sub);
      double error = 0;
      for (size_t i = 0; i < correct_test_output.size(); ++i) {
        if (i < sub.size()) {
          // Add error.
          error += emp::Abs((double)((double)correct_test_output[i] - (double)sub[i]));
        } else {
          // Add max error.
          error += (2*MAX_NUM);
        }
      }
      // Add error for each extra thing in sub
      if (sub.size() > correct_test_output.size()) {
        error += (2*MAX_NUM) * (sub.size() - correct_test_output.size());
      }
      // Make sure programs that overflow error to be negative (WTF?) don't get really high scores.
      double score = (error <= MAX_ERROR && error >= 0) ? 1 - (error/(double)MAX_ERROR) : 0.0;

      // if (score > 1000) {
      //   std::cout << "Score is greater than 1000!" << std::endl;
      //   std::cout << "  Score = " << score << std::endl;
      //   std::cout << "  Max number = " << MAX_NUM << std::endl;
      //   std::cout << "  Max error = " << MAX_ERROR << std::endl;
      //   std::cout << "  Calculated error = " << error << std::endl;
      //   std::cout << "  Correct output: [";
      //   for (size_t i = 0; i < correct_test_output.size(); ++i) {
      //     if (i) std::cout << ",";
      //     std::cout << correct_test_output[i];
      //   } std::cout << "]" << std::endl;
      //   std::cout << "  Sub output: [";
      //   for (size_t i = 0; i < sub.size(); ++i) {
      //     if (i) std::cout << ",";
      //     std::cout << sub[i];
      //   } std::cout << "]" << std::endl;
      // }

      return {score, false};
    }
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << "[[";
    for (size_t i = 0; i < in[0].size(); ++i) {
      if (i) os << ",";
      os << in[0][i];
    }
    os << "],["; 
    for (size_t i = 0; i < in[1].size(); ++i) {
      if (i) os << ",";
      os << in[1][i];
    }
    os << "]]"; 
    os << "\"";
  }
};




/// X-Word Lines: Pair<Integer, String>
class TestOrg_XWordLines : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = std::pair<int, std::string>;
  protected:
    genome_t genome;

  public:
    TestOrg_XWordLines(const genome_t & _g) : genome(_g) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { ; }
};

struct ProblemUtilities_XWordLines { emp::vector<std::function<double(TestOrg_XWordLines &)>> lexicase_fit_set; };



/// Pig Latin: String
class TestOrg_PigLatin : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = std::string;
  protected:
    genome_t genome;

  public:
    TestOrg_PigLatin(const genome_t & _g) : genome(_g) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { ; }
};

struct ProblemUtilities_PigLatin { emp::vector<std::function<double(TestOrg_PigLatin &)>> lexicase_fit_set; };



/// Negative To Zero: Vector<Integer>
class TestOrg_NegativeToZero : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = emp::vector<int>;
  protected:
    genome_t genome;

  public:
    TestOrg_NegativeToZero(const genome_t & _g) : genome(_g) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { ; }
};

struct ProblemUtilities_NegativeToZero { emp::vector<std::function<double(TestOrg_NegativeToZero &)>> lexicase_fit_set; };



/// Scrabble Score: String
class TestOrg_ScrabbleScore : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = std::string;
  protected:
    genome_t genome;

  public:
    TestOrg_ScrabbleScore(const genome_t & _g) : genome(_g) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { ; }
};

struct ProblemUtilities_ScrabbleScore { emp::vector<std::function<double(TestOrg_ScrabbleScore &)>> lexicase_fit_set; };



/// Word Stats: File
// class TestOrg_WordStats : public TestOrg_Base {
//   public:
//     using parent_t = TestOrg_Base;
//     using parent_t::phenotype;

//     using genome_t = ;
//   protected:
//     genome_t genome;

//   public:
//     TestOrg_WordStats (const genome_t & _g) : genome(_g) { ; }
    
//     genome_t & GetGenome() { return genome; }
//     const genome_t & GetGenome() const { return genome; }
// };


/// Checksum: String
class TestOrg_Checksum : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = std::string;
  protected:
    genome_t genome;

  public:
    TestOrg_Checksum(const genome_t & _g) : genome(_g) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { ; }
};

struct ProblemUtilities_Checksum { emp::vector<std::function<double(TestOrg_Checksum &)>> lexicase_fit_set; };



/// Digits: Integer
class TestOrg_Digits : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = int;
  protected:
    genome_t genome;

  public:
    TestOrg_Digits(const genome_t & _g) : genome(_g) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { ; }

};

struct ProblemUtilities_Digits { emp::vector<std::function<double(TestOrg_Digits &)>> lexicase_fit_set; };


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: Grade
// - Input: Array<Integer, 5>
// - Output: String
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Generate random test input
Problem_Grade_input_t GenRandomTestInput_Grade(emp::Random & rand, const std::pair<int, int> & num_range) {
  emp_assert(num_range.first < num_range.second);
  emp_assert(num_range.first > -1);
  Problem_Grade_input_t input;
  input[0] = -1;
  input[1] = -1;
  input[2] = -1;
  input[3] = -1;
  input[4] = -1;

  for (size_t i = 0; i < 4; ++i) {
    int val = rand.GetInt(num_range.first+1, num_range.second);
    // while (emp::Has(input, val)) val = rand.GetInt(num_range.first+1, num_range.second);
    while (input[0] == val || input[1] == val || input[2] == val || input[3] == val) val = rand.GetInt(num_range.first+1, num_range.second);
    input[i] = val;
  }
  // Sort input
  std::sort(input.begin(), input.end());
  std::reverse(std::begin(input), std::end(input));
  input[4] = rand.GetInt(num_range.first, num_range.second+1); // Grade?

  // Test!
  // std::cout << "A thresh: " << input[0] << std::endl;
  // std::cout << "B thresh: " << input[1] << std::endl;
  // std::cout << "C thresh: " << input[2] << std::endl;
  // std::cout << "D thresh: " << input[3] << std::endl;
  // std::cout << "Grade: " << input[4] << std::endl;
  emp_assert(100 >= input[0], input[0], num_range.second);
  emp_assert(num_range.second >= input[0]);
  emp_assert(input[0] > input[1]);
  emp_assert(input[1] > input[2]);
  emp_assert(input[2] > input[3]);
  emp_assert(input[3] >= num_range.first);
  emp_assert(input[3] >= 0);
  return input;
}

// Generate correct output for a given test
Problem_Grade_output_t GenCorrectOut_Grade(const Problem_Grade_input_t & input) {
  if (input[4] >= input[0]) { return Grade__A_STR; }
  else if (input[4] >= input[1]) { return Grade__B_STR; }
  else if (input[4] >= input[2]) { return Grade__C_STR; }
  else if (input[4] >= input[3]) { return Grade__D_STR; }
  else { return Grade__F_STR; }
}

/// Grade: Array<Integer, 5>
class TestOrg_Grade : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = Problem_Grade_input_t;
    using out_t = Problem_Grade_output_t;

  protected:
    genome_t genome;
    out_t out;

  public:
    TestOrg_Grade(const genome_t & _g) : genome(_g), out() { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    out_t & GetCorrectOut() { return out; }
    const out_t & GetCorrectOut() const { return out; }

    void CalcOut() { out = GenCorrectOut_Grade(genome); }    

    void Print(std::ostream & os=std::cout) {
      for (size_t i = 0; i < genome.size(); ++i) {
        if (i) os << ",";
        os << genome[i];
      }
    }
  };

struct ProblemUtilities_Grade { 
  using this_t = ProblemUtilities_Grade;
  using problem_org_t = TestOrg_Grade;
  using input_t = Problem_Grade_input_t;
  using output_t = Problem_Grade_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  emp::vector<emp::Ptr<problem_org_t>> testingset_pop;
  emp::vector<emp::vector<output_t>> population_validation_outputs;

  // // --- Useful during a test evaluation ---
  emp::Ptr<problem_org_t> cur_eval_test_org;
  bool submitted;
  std::string submitted_str;

  // // Mutation - Handle here...
  int MIN_NUM;
  int MAX_NUM;
  double PER_NUM_ADJUST_RATE;
  double PER_NUM_RANDOMIZE_RATE;
  
  size_t Mutate(emp::Random & rnd, input_t & mut_input) {
    // std::cout << "Mutate!" << std::endl;
    // std::cout << "A thresh: " << mut_input[0] << std::endl;
    // std::cout << "B thresh: " << mut_input[1] << std::endl;
    // std::cout << "C thresh: " << mut_input[2] << std::endl;
    // std::cout << "D thresh: " << mut_input[3] << std::endl;
    // std::cout << "Grade: " << mut_input[4] << std::endl;
    // Are grade thresholds valid?
    emp_assert(MAX_NUM  >= mut_input[0], MAX_NUM, mut_input[0]);
    emp_assert(mut_input[0] > mut_input[1]);
    emp_assert(mut_input[1] > mut_input[2]);
    emp_assert(mut_input[2] > mut_input[3]);
    emp_assert(mut_input[3] >= MIN_NUM);
    // Is grade valid?
    emp_assert(mut_input[4] >= MIN_NUM);
    emp_assert(mut_input[4] <= MAX_NUM);

    size_t muts = 0;
    const size_t grade_id = 4;

    // Adjust mutations
    for (size_t i = 0; i < grade_id; ++i) {
      int upper_cap = MAX_NUM;
      if (i) upper_cap = mut_input[i-1];

      int lower_cap = MIN_NUM;
      if (i < 3) lower_cap = mut_input[i+1]+1;
      
      if (rnd.P(PER_NUM_ADJUST_RATE)) { 
        mut_input[i] = rnd.GetInt(lower_cap, upper_cap);
        ++muts;
      }
    }

    int orig_grade = mut_input[grade_id];
    mut_input[grade_id] = -1;

    // Randomize mutations
    for (size_t i = 0; i < grade_id; ++i) {
      if (rnd.P(PER_NUM_RANDOMIZE_RATE)) {
        mut_input[i] = -1;
        int val = rnd.GetInt(MIN_NUM+1, MAX_NUM);
        while (mut_input[0] == val || mut_input[1] == val || mut_input[2] == val || mut_input[3] == val) val = rnd.GetInt(MIN_NUM+1, MAX_NUM);
        mut_input[i] = val;
        ++muts;
      }
    }

    // Guarantee distinctness
    for (size_t i = 0; i < grade_id; ++i) {
      int val = mut_input[i];
      mut_input[i] = -1;
      while (mut_input[0] == val || mut_input[1] == val || mut_input[2] == val || mut_input[3] == val) val = rnd.GetInt(MIN_NUM+1, MAX_NUM);
      mut_input[i] = val;
    }

    // Repair ordering
    std::sort(mut_input.begin(), mut_input.end());
    std::reverse(std::begin(mut_input), std::end(mut_input));
    
    // Randomize grade?
    if (rnd.P(PER_NUM_RANDOMIZE_RATE)) {
      mut_input[4] = rnd.GetInt(MIN_NUM, MAX_NUM+1); 
      ++muts;
    } else {
      mut_input[4] = orig_grade;
    }

    // if (!(mut_input[2] > mut_input[3])) {
    //   std::cout << "vvvv" << std::endl;
    //   std::cout << "A thresh: " << mut_input[0] << std::endl;
    //   std::cout << "B thresh: " << mut_input[1] << std::endl;
    //   std::cout << "C thresh: " << mut_input[2] << std::endl;
    //   std::cout << "D thresh: " << mut_input[3] << std::endl;
    //   std::cout << "Grade: " << mut_input[4] << std::endl;
    // }
    // std::cout << "---------------" << std::endl;

    
    
    // Are grade thresholds valid?
    emp_assert(MAX_NUM  >= mut_input[0], MAX_NUM, mut_input[0]);
    emp_assert(mut_input[0] > mut_input[1], mut_input[0], mut_input[1]);
    emp_assert(mut_input[1] > mut_input[2], mut_input[1], mut_input[2]);
    emp_assert(mut_input[2] > mut_input[3], mut_input[2], mut_input[3]);
    emp_assert(mut_input[3] >= MIN_NUM, mut_input[3], MIN_NUM);
    // Is grade valid?
    emp_assert(mut_input[4] >= MIN_NUM, mut_input[4], MIN_NUM);
    emp_assert(mut_input[4] <= MAX_NUM, mut_input[4], MAX_NUM);

    return muts;
  }

  // Selection
  emp::vector<std::function<double(problem_org_t &)>> lexicase_fit_set;

  ProblemUtilities_Grade()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_str("")
  { ; }

  ~ProblemUtilities_Grade() {
    for (size_t i = 0; i < testingset_pop.size(); ++i) testingset_pop[i].Delete();
  }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  void ResetTestEval() {
    submitted = false;
    submitted_str = "";
  }

  void Submit(const std::string & val) {
    submitted = true;
    submitted_str = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;
    output_t output;

    // Load input
    input[0] = std::atof(line[0].c_str());
    input[1] = std::atof(line[1].c_str());
    input[2] = std::atof(line[2].c_str());
    input[3] = std::atof(line[3].c_str());
    input[4] = std::atof(line[4].c_str());

    // std::cout << "A thresh: " << input[0] << std::endl;
    // std::cout << "B thresh: " << input[1] << std::endl;
    // std::cout << "C thresh: " << input[2] << std::endl;
    // std::cout << "D thresh: " << input[3] << std::endl;
    // std::cout << "Grade: " << input[4] << std::endl;

    if (line[5] == "Student has a A grade.") {
      output = Grade__A_STR;
    } else if (line[5] == "Student has a B grade.") {
      output = Grade__B_STR;
    } else if (line[5] == "Student has a C grade.") {
      output = Grade__C_STR;
    } else if (line[5] == "Student has a D grade.") {
      output = Grade__D_STR;
    } else if (line[5] == "Student has a F grade.") {
      output = Grade__F_STR;
    } else {
      std::cout << "ERROR ERROR! OH NO! INVALID OUTPUT FROM GRADE EXAMPLES!" << std::endl;
      exit(-1);
    }

    output_t gen_out = GenCorrectOut_Grade(input);
    emp_assert(gen_out == output);

    emp_assert(100 >= input[0]);
    emp_assert(input[0] > input[1]);
    emp_assert(input[1] > input[2]);
    emp_assert(input[2] > input[3]);
    emp_assert(input[3] >= 0);

    return {input, output};

  }

  void GenerateTestingSetPop() {
    for (size_t i = 0; i < testing_set.GetSize(); ++i) {
      testingset_pop.emplace_back(emp::NewPtr<problem_org_t>(testing_set.GetInput(i)));
      testingset_pop[i]->CalcOut();
    }
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    for (size_t i = 0; i < in.size(); ++i) {
      if (i) os << ",";
      os << in[i];
    }
    os << "\"";
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: Median
// - Input: std::array<int, 3>;
// - Output: int;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Generate random test input
Problem_Median_input_t GenRandomTestInput_Median(emp::Random & rand,
                                                 const std::pair<int, int> & num_range) {
  Problem_Median_input_t input;
  for (size_t i = 0; i < input.size(); ++i) {
    input[i] = rand.GetInt(num_range.first, num_range.second+1);
  }
  return input;
}

/// Generate correct output
Problem_Median_output_t GenCorrectOut_Median(const Problem_Median_input_t & input) {
  const int min_val = emp::Min(input[0], input[1], input[2]);
  const int max_val = emp::Max(input[0], input[1], input[2]);
  // std::cout << "Input: " << input[0] << ", " << input[1] << ", " << input[2] << "; Output: " << (input[0] + input[1] + input[2]) - min_val - max_val << std::endl;
  return (input[0] + input[1] + input[2]) - min_val - max_val;
}

/// Median: Array<Integer, 3>
class TestOrg_Median : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = std::array<int, 3>;
    using out_t = Problem_Median_output_t;

  protected:
    genome_t genome;
    out_t out;

  public:
    TestOrg_Median(const genome_t & _g) : genome(_g), out() { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { out = GenCorrectOut_Median(genome); }

    out_t & GetCorrectOut() { return out; }
    const out_t & GetCorrectOut() const { return out; }  

    void Print(std::ostream & os=std::cout) {
      for (size_t i = 0; i < genome.size(); ++i) {
        if (i) os << ",";
        os << genome[i];
      }
    }
};

struct ProblemUtilities_Median { 
  using this_t = ProblemUtilities_Median;
  using problem_org_t = TestOrg_Median;
  using input_t = Problem_Median_input_t;
  using output_t = Problem_Median_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  emp::vector<emp::Ptr<problem_org_t>> testingset_pop;
  emp::vector<emp::vector<output_t>> population_validation_outputs;

  // // --- Useful during a test evaluation ---
  emp::Ptr<problem_org_t> cur_eval_test_org;
  bool submitted;
  int submitted_val;

  // // Mutation - Handle here...
  int MIN_NUM;
  int MAX_NUM;
  double PER_NUM_COPY_RATE;
  double PER_NUM_SUB_RATE;
  double PER_NUM_SWAP_RATE;
  
  size_t Mutate(emp::Random & rnd, input_t & mut_input) {
    size_t muts = 0;

    for (size_t i = 0; i < mut_input.size(); ++i) {
      if (rnd.P(PER_NUM_SUB_RATE)) {
        ++muts;
        mut_input[i] = rnd.GetInt(MIN_NUM, MAX_NUM+1);
      }
    }

    for (size_t i = 0; i < mut_input.size(); ++i) {
      if (rnd.P(PER_NUM_SWAP_RATE)) {
        ++muts;
        size_t other_pos = rnd.GetUInt(mut_input.size());
        while (other_pos == i) other_pos = rnd.GetUInt(mut_input.size());
        std::swap(mut_input[i], mut_input[other_pos]);
      }
    }

    for (size_t i = 0; i < mut_input.size(); ++i) {
      if (rnd.P(PER_NUM_COPY_RATE)) {
        ++muts;
        size_t other_pos = rnd.GetUInt(mut_input.size());
        while (other_pos == i) other_pos = rnd.GetUInt(mut_input.size());
        mut_input[other_pos] = mut_input[i];
      }
    }

    return muts;
  } // todo - test

  // Selection
  emp::vector<std::function<double(problem_org_t &)>> lexicase_fit_set;

  ProblemUtilities_Median()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val(0)
  { ; }

  ~ProblemUtilities_Median() {
    for (size_t i = 0; i < testingset_pop.size(); ++i) testingset_pop[i].Delete();
  }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  void ResetTestEval() {
    submitted = false;
    submitted_val = 0;
  }

  void Submit(int val) {
    submitted = true;
    submitted_val = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    // Load input.
    input[0] = std::atof(line[0].c_str());
    input[1] = std::atof(line[1].c_str());
    input[2] = std::atof(line[2].c_str());
    // Load output.
    output = std::atof(line[3].c_str());
    emp_assert(output == GenCorrectOut_Median(input));
    return {input, output};
  }

  void GenerateTestingSetPop() {
    for (size_t i = 0; i < testing_set.GetSize(); ++i) {
      testingset_pop.emplace_back(emp::NewPtr<problem_org_t>(testing_set.GetInput(i)));
      testingset_pop[i]->CalcOut();
    }
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    for (size_t i = 0; i < in.size(); ++i) {
      if (i) os << ",";
      os << in[i];
    }
    os << "\"";
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: Smallest
// - Input: Array<Integer, 4>
// - Output: int
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Generate random test input
Problem_Smallest_input_t GenRandomTestInput_Smallest(emp::Random & rand,
                                                     const std::pair<int, int> & num_range) {
  Problem_Smallest_input_t input;
  for (size_t i = 0; i < input.size(); ++i) {
    input[i] = rand.GetInt(num_range.first, num_range.second+1);
  }
  return input;
}

/// Generate correct output
Problem_Smallest_output_t GenCorrectOut_Smallest(const Problem_Smallest_input_t & input) {
  int smallest = input[0];
  for (size_t i = 1; i < input.size(); ++i) {
    if (input[i] < smallest) smallest = input[i];
  }
  return smallest;
}

/// Smallest: Array<Integer, 4>
class TestOrg_Smallest : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = std::array<int, 4>;
    using out_t = Problem_Smallest_output_t;

  protected:
    genome_t genome;
    out_t out;

  public:
    TestOrg_Smallest(const genome_t & _g) : genome(_g) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { out = GenCorrectOut_Smallest(genome); }

    out_t & GetCorrectOut() { return out; }
    const out_t & GetCorrectOut() const { return out; }  

    void Print(std::ostream & os=std::cout) {
      for (size_t i = 0; i < genome.size(); ++i) {
        if (i) os << ",";
        os << genome[i];
      }
    }
};

struct ProblemUtilities_Smallest { 
  using this_t = ProblemUtilities_Smallest;
  using problem_org_t = TestOrg_Smallest;
  using input_t = Problem_Smallest_input_t;
  using output_t = Problem_Smallest_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  emp::vector<emp::Ptr<problem_org_t>> testingset_pop;
  emp::vector<emp::vector<output_t>> population_validation_outputs;

  // // --- Useful during a test evaluation ---
  emp::Ptr<problem_org_t> cur_eval_test_org;
  bool submitted;
  int submitted_val;

  // // Mutation - Handle here...
  int MIN_NUM;
  int MAX_NUM;
  double PER_NUM_SUB_RATE;
  double PER_NUM_SWAP_RATE;
  
  size_t Mutate(emp::Random & rnd, input_t & mut_input) {
    size_t muts = 0;

    for (size_t i = 0; i < mut_input.size(); ++i) {
      if (rnd.P(PER_NUM_SUB_RATE)) {
        ++muts;
        mut_input[i] = rnd.GetInt(MIN_NUM, MAX_NUM+1);
      }
    }

    for (size_t i = 0; i < mut_input.size(); ++i) {
      if (rnd.P(PER_NUM_SWAP_RATE)) {
        ++muts;
        size_t other_pos = rnd.GetUInt(mut_input.size());
        while (other_pos == i) other_pos = rnd.GetUInt(mut_input.size());
        std::swap(mut_input[i], mut_input[other_pos]);
      }
    }

    return muts;
  } // todo - test

  // Selection
  emp::vector<std::function<double(problem_org_t &)>> lexicase_fit_set;

  ProblemUtilities_Smallest()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val(false)
  { ; }

  ~ProblemUtilities_Smallest() {
    for (size_t i = 0; i < testingset_pop.size(); ++i) testingset_pop[i].Delete();
  }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  void ResetTestEval() {
    submitted = false;
    submitted_val = 0;
  }

  void Submit(int val) {
    submitted = true;
    submitted_val = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    // Load input.
    input[0] = std::atof(line[0].c_str());
    input[1] = std::atof(line[1].c_str());
    input[2] = std::atof(line[2].c_str());
    input[3] = std::atof(line[3].c_str());
    // Load output.
    output = std::atof(line[4].c_str());
    emp_assert(output == GenCorrectOut_Smallest(input));
    return {input, output};
  }

  void GenerateTestingSetPop() {
    for (size_t i = 0; i < testing_set.GetSize(); ++i) {
      testingset_pop.emplace_back(emp::NewPtr<problem_org_t>(testing_set.GetInput(i)));
      testingset_pop[i]->CalcOut();
    }
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    for (size_t i = 0; i < in.size(); ++i) {
      if (i) os << ",";
      os << in[i];
    }
    os << "\"";
  }

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: Syllables
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Syllables: String
class TestOrg_Syllables : public TestOrg_Base {
  public:
    using parent_t = TestOrg_Base;
    using parent_t::phenotype;

    using genome_t = std::string;
  protected:
    genome_t genome;

  public:
    TestOrg_Syllables(const genome_t & _g) : genome(_g) { ; }
    
    genome_t & GetGenome() { return genome; }
    const genome_t & GetGenome() const { return genome; }

    void CalcOut() { ; }
};


struct ProblemUtilities_Syllables { emp::vector<std::function<double(TestOrg_Syllables &)>> lexicase_fit_set; };

#endif