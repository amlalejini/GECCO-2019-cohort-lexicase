#ifndef PROGRAMMING_SYNTHESIS_EXP_H
#define PROGRAMMING_SYNTHESIS_EXP_H

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <utility>
#include <unordered_set>
#include <tuple>
// #include <variant>
#include <cassert>

#include "base/Ptr.h"
#include "base/vector.h"
#include "control/Signal.h"
#include "Evolve/World.h"
#include "Evolve/World_select.h"
#include "tools/Random.h"
#include "tools/random_utils.h"
#include "tools/math.h"
#include "tools/string_utils.h"
#include "tools/stats.h"
#include "tools/tuple_utils.h"

#include "TagLinearGP.h"
#include "Selection.h"
#include "Mutators.h"

#include "ProgOrg.h"
#include "ProgSynthConfig.h"
#include "ProgSynthBenchmarks_InputReps.h"
#include "TestCaseSet.h"
#include "TagLinearGP.h"
#include "TagLinearGP_InstLib.h"
#include "TagLinearGP_Utilities.h"

//////////////////////////////////////////
// --- Notes/Todos ---
// - May need to generate more training examples for problems(?)
// - INSTRUCTION SET
//  - [ ] Add LoadAllSetInputs
// - SELECTION
//  - [ ] Pools (as in Cliff's implementation of lexicase) (?)
// - SCORING
//  - [x] Assume pass/fail only at first, next add gradient (NOTE - will need to update how we screen/add more things to phenotypes).
//    - [x] Simplest thing to do would be to add a pass_vector + score_vector to test/program phenotypes
//  - Add submission test case(?) to encourage submitting *something* ==> For now, not doing this.
//  - [x] Break apart score w/passes
//  - [x] Add selection pressure for small sizes
// - DIAGNOSTICS
//  - [ ] Clean up printing format
// - Cleaning
//  - [ ] Flesh out config comments
//  - [ ] Round out instruction descriptions
// - ISSUES
//  - Memory indexing seems to be the biggest impact on evaluation speed. 
//    Every instruction argument requires a linear scan of memory for best match.
// --------------
// Scratch
// -----
// How do we want to evaluate on the testing set?
// - Evaluate
//////////////////////////////////////////

constexpr size_t TAG_WIDTH = 16;
constexpr size_t MEM_SIZE = TAG_WIDTH;

// How do training examples change over time?
// - coevolution - training examples co-evolve with programs
// - static - training examples are static
// - random - training examples randomly change over time
// - static gen - use loaded training examples, but generate more random examples
enum TRAINING_EXAMPLE_MODE_TYPE { COEVOLUTION=0, STATIC, RANDOM, STATIC_GEN, STATIC_COEVO };
enum EVALUATION_TYPE { COHORT=0, FULL=1, PROG_ONLY_COHORT=2, TEST_DOWNSAMPLING };
enum SELECTION_TYPE { LEXICASE=0, COHORT_LEXICASE, TOURNAMENT, DRIFT, PROG_ONLY_COHORT_LEXICASE, TEST_DOWNSAMPLING_LEXICASE };

enum PROBLEM_ID { NumberIO=0,
                  SmallOrLarge,
                  ForLoopIndex,
                  CompareStringLengths,
                  DoubleLetters,
                  CollatzNumbers,
                  ReplaceSpaceWithNewline,
                  StringDifferences,
                  EvenSquares,
                  WallisPi,
                  StringLengthsBackwards,
                  LastIndexOfZero,
                  VectorAverage,
                  CountOdds,
                  MirrorImage,
                  SuperAnagrams,
                  SumOfSquares,
                  VectorsSummed,
                  XWordLines,
                  PigLatin,
                  NegativeToZero,
                  ScrabbleScore,
                  Checksum,
                  Digits,
                  Grade,
                  Median,
                  Smallest,
                  Syllables
};

struct ProblemInfo {
  PROBLEM_ID id; 
  std::string training_fname;
  std::string testing_fname;
  
  ProblemInfo(PROBLEM_ID _id, const std::string & _training_fname, const std::string & _testing_fname) 
    : id(_id), training_fname(_training_fname), testing_fname(_testing_fname)
  { ; }
  
  ProblemInfo(const ProblemInfo &) = default;
  ProblemInfo(ProblemInfo &&) = default;

  ProblemInfo & operator=(const ProblemInfo &) = default;
  ProblemInfo & operator=(ProblemInfo &&) = default;

  const std::string & GetTestingSetFilename() const { return testing_fname; }
  const std::string & GetTrainingSetFilename() const { return training_fname; }

};

std::unordered_map<std::string, ProblemInfo> problems = {
  {"number-io", {PROBLEM_ID::NumberIO, "training-examples-number-io.csv", "testing-examples-number-io.csv"}},
  {"small-or-large", {PROBLEM_ID::SmallOrLarge, "training-examples-small-or-large.csv", "testing-examples-small-or-large.csv"}},
  {"for-loop-index", {PROBLEM_ID::ForLoopIndex, "training-examples-for-loop-index.csv", "testing-examples-for-loop-index.csv"}},
  {"compare-string-lengths", {PROBLEM_ID::CompareStringLengths, "training-examples-compare-string-lengths.csv", "testing-examples-compare-string-lengths.csv"}},
  {"collatz-numbers", {PROBLEM_ID::CollatzNumbers, "training-examples-collatz-numbers.csv", "testing-examples-collatz-numbers.csv"}},
  {"string-lengths-backwards", {PROBLEM_ID::StringLengthsBackwards, "training-examples-string-lengths-backwards.csv", "testing-examples-string-lengths-backwards.csv"}},
  {"last-index-of-zero", {PROBLEM_ID::LastIndexOfZero, "training-examples-last-index-of-zero.csv", "testing-examples-last-index-of-zero.csv"}},
  {"count-odds", {PROBLEM_ID::CountOdds, "training-examples-count-odds.csv", "testing-examples-count-odds.csv"}},
  {"mirror-image", {PROBLEM_ID::MirrorImage, "training-examples-mirror-image.csv", "testing-examples-mirror-image.csv"}},
  {"vectors-summed", {PROBLEM_ID::VectorsSummed, "training-examples-vectors-summed.csv", "testing-examples-vectors-summed.csv"}},
  {"sum-of-squares", {PROBLEM_ID::SumOfSquares, "training-examples-sum-of-squares.csv", "testing-examples-sum-of-squares.csv"}},
  {"vector-average", {PROBLEM_ID::VectorAverage, "training-examples-vector-average.csv", "testing-examples-vector-average.csv"}},
  {"median", {PROBLEM_ID::Median, "training-examples-median.csv", "testing-examples-median.csv"}},
  {"smallest", {PROBLEM_ID::Smallest, "training-examples-smallest.csv", "testing-examples-smallest.csv"}},
  {"grade", {PROBLEM_ID::Grade, "training-examples-grade.csv", "testing-examples-grade.csv"}}
};

class ProgramSynthesisExperiment {
public:
  using hardware_t = typename TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using inst_lib_t = typename TagLGP::InstLib<hardware_t>;
  using inst_t = typename hardware_t::inst_t;

  using prog_org_t = ProgOrg<TAG_WIDTH>;
  using prog_org_phen_t = typename prog_org_t::Phenotype;
  using prog_org_gen_t = typename prog_org_t::genome_t;

  using prog_world_t = emp::World<prog_org_t>;

  using test_org_phen_t = TestOrg_Base::Phenotype;

  using prog_taxon_t = typename emp::Systematics<prog_org_t, prog_org_gen_t>::taxon_t;

  // test world aliases
  using prob_NumberIO_world_t = emp::World<TestOrg_NumberIO>;
  using prob_SmallOrLarge_world_t = emp::World<TestOrg_SmallOrLarge>;
  using prob_ForLoopIndex_world_t = emp::World<TestOrg_ForLoopIndex>;
  using prob_CompareStringLengths_world_t = emp::World<TestOrg_CompareStringLengths>;
  using prob_DoubleLetters_world_t = emp::World<TestOrg_DoubleLetters>;
  using prob_CollatzNumbers_world_t = emp::World<TestOrg_CollatzNumbers>;
  using prob_ReplaceSpaceWithNewline_world_t = emp::World<TestOrg_ReplaceSpaceWithNewline>;
  using prob_StringDifferences_world_t = emp::World<TestOrg_StringDifferences>;
  using prob_EvenSquares_world_t = emp::World<TestOrg_EvenSquares>;
  using prob_WallisPi_world_t = emp::World<TestOrg_WallisPi>;
  using prob_StringLengthsBackwards_world_t = emp::World<TestOrg_StringLengthsBackwards>;
  using prob_LastIndexOfZero_world_t = emp::World<TestOrg_LastIndexOfZero>;
  using prob_VectorAverage_world_t = emp::World<TestOrg_VectorAverage>;
  using prob_CountOdds_world_t = emp::World<TestOrg_CountOdds>;
  using prob_MirrorImage_world_t = emp::World<TestOrg_MirrorImage>;
  using prob_SuperAnagrams_world_t = emp::World<TestOrg_SuperAnagrams>;
  using prob_SumOfSquares_world_t = emp::World<TestOrg_SumOfSquares>;
  using prob_VectorsSummed_world_t = emp::World<TestOrg_VectorsSummed>;
  using prob_XWordLines_world_t = emp::World<TestOrg_XWordLines>;
  using prob_PigLatin_world_t = emp::World<TestOrg_PigLatin>;
  using prob_NegativeToZero_world_t = emp::World<TestOrg_NegativeToZero>;
  using prob_ScrabbleScore_world_t = emp::World<TestOrg_ScrabbleScore>;
  using prob_Checksum_world_t = emp::World<TestOrg_Checksum>;
  using prob_Digits_world_t = emp::World<TestOrg_Digits>;
  using prob_Grade_world_t = emp::World<TestOrg_Grade>;
  using prob_Median_world_t = emp::World<TestOrg_Median>;
  using prob_Smallest_world_t = emp::World<TestOrg_Smallest>;
  using prob_Syllables_world_t = emp::World<TestOrg_Syllables>;

protected:

  // Useful experiment structs
  struct Cohorts {
    protected:
      emp::vector<size_t> population_ids;
      emp::vector<emp::vector<size_t>> cohorts;

      size_t cohort_size;

      bool init;

      /// Used internally to access population_ids
      size_t GetOrgPopID(size_t cohortID, size_t memberID) const {
        emp_assert(init);
        emp_assert(cohortID < cohorts.size());
        emp_assert(memberID < cohorts[cohortID].size());
        return (cohortID * cohort_size) + memberID;
      }

    public:
      Cohorts() : population_ids(), cohorts(), cohort_size(0), init(false) { ; }

      /// Setup cohorts
      void Setup(size_t _pop_size, size_t _cohort_size) {
        emp_assert(_cohort_size > 0);
        population_ids.clear();
        cohorts.clear();
        cohort_size = _cohort_size;
        const size_t num_cohorts = _pop_size / _cohort_size;
        // Initialize population IDS
        for (size_t i = 0; i < _pop_size; ++i) population_ids.emplace_back(i);
        init = true;
        // Initialize cohort vectors
        for (size_t cID = 0; cID < num_cohorts; ++cID) {
          cohorts.emplace_back(cohort_size);
          for (size_t i = 0; i < cohort_size; ++i) cohorts[cID][i] = population_ids[GetOrgPopID(cID, i)];
        }
      }

      size_t GetCohortCnt() const { return cohorts.size(); }
      size_t GetCohortSize() const { return cohort_size; }

      const emp::vector<size_t> & GetCohort(size_t cohortID) const {
        emp_assert(cohortID < cohorts.size());
        return cohorts[cohortID];
      }
      
      /// Used to get world ID of an organism given by its cohort ID and member ID.
      size_t GetWorldID(size_t cohortID, size_t memberID) const {
        emp_assert(init);
        emp_assert(cohortID < cohorts.size());
        emp_assert(memberID < cohorts[cohortID].size());
        return cohorts[cohortID][memberID];
      }

      /// Randomize cohort assignment.
      void Randomize(emp::Random & rnd) {
        // Shuffle population ids
        emp::Shuffle(rnd, population_ids);
        // Reassign cohorts.
        for (size_t cID = 0; cID < cohorts.size(); ++cID) {
          for (size_t i = 0; i < cohort_size; ++i) cohorts[cID][i] = population_ids[GetOrgPopID(cID, i)];
        }
      }

      /// Print cohorts.
      void Print(std::ostream & out=std::cout) const {
        for (size_t cID = 0; cID < cohorts.size(); ++cID) {
          if (cID) out << "\n";
          out << "Cohort[" << cID << "]: [";
          for (size_t i = 0; i < cohorts[cID].size(); ++i) {
            if (i) out << ", ";
            out << cohorts[cID][i];
          }
          out << "]";
        }
      }
  };

  // --- Localized experiment parameters ---
  int SEED;
  size_t GENERATIONS;
  size_t PROG_POP_SIZE;
  size_t TEST_POP_SIZE;
  size_t EVALUATION_MODE;
  size_t PROG_COHORT_SIZE;
  size_t TEST_COHORT_SIZE;
  size_t TRAINING_EXAMPLE_MODE;
  std::string PROBLEM;
  std::string BENCHMARK_DATA_DIR;

  size_t PROG_SELECTION_MODE;
  size_t TEST_SELECTION_MODE;
  size_t PROG_LEXICASE_MAX_FUNS;
  size_t PROG_COHORTLEXICASE_MAX_FUNS;
  size_t TEST_LEXICASE_MAX_FUNS;
  size_t TEST_COHORTLEXICASE_MAX_FUNS;
  size_t PROG_TOURNAMENT_SIZE;
  size_t TEST_TOURNAMENT_SIZE;
  bool DISCRIMINATORY_LEXICASE_TESTS;

  size_t MIN_PROG_SIZE;
  size_t MAX_PROG_SIZE;
  size_t PROG_EVAL_TIME;
  double PROG_MUT__PER_BIT_FLIP;
  double PROG_MUT__PER_INST_SUB;
  double PROG_MUT__PER_INST_INS;
  double PROG_MUT__PER_INST_DEL;
  double PROG_MUT__PER_PROG_SLIP;
  double PROG_MUT__PER_MOD_DUP;
  double PROG_MUT__PER_MOD_DEL;

  double MIN_TAG_SPECIFICITY;
  size_t MAX_CALL_DEPTH;

  double PROB_NUMBER_IO__DOUBLE_MIN;
  double PROB_NUMBER_IO__DOUBLE_MAX;
  int PROB_NUMBER_IO__INT_MIN;
  int PROB_NUMBER_IO__INT_MAX;
  double PROB_NUMBER_IO__MUTATION__PER_INT_RATE;
  double PROB_NUMBER_IO__MUTATION__PER_DOUBLE_RATE; 

  int PROB_SMALL_OR_LARGE__INT_MIN;
  int PROB_SMALL_OR_LARGE__INT_MAX;
  double PROB_SMALL_OR_LARGE__MUTATION__PER_INT_RATE;

  int PROB_FOR_LOOP_INDEX__START_END_MIN;
  int PROB_FOR_LOOP_INDEX__START_END_MAX;
  int PROB_FOR_LOOP_INDEX__STEP_MIN;
  int PROB_FOR_LOOP_INDEX__STEP_MAX;
  double PROB_FOR_LOOP_INDEX__MUTATION__MUT_RATE;
  bool PROB_FOR_LOOP_INDEX__PROMISE_MULTISTEP_TESTCASES;

  size_t PROB_COMPARE_STRING_LENGTHS__MIN_STR_LEN;
  size_t PROB_COMPARE_STRING_LENGTHS__MAX_STR_LEN;
  double PROB_COMPARE_STRING_LENGTHS__PER_SITE_INS_RATE;
  double PROB_COMPARE_STRING_LENGTHS__PER_SITE_DEL_RATE;
  double PROB_COMPARE_STRING_LENGTHS__PER_SITE_SUB_RATE;
  double PROB_COMPARE_STRING_LENGTHS__PER_STR_SWAP_RATE;

  int PROB_COLLATZ_NUMBERS__MIN_NUM;
  int PROB_COLLATZ_NUMBERS__MAX_NUM;
  double PROB_COLLATZ_NUMBERS__MUTATION__PER_NUM_SUB_RATE;

  size_t PROB_STRING_LENGTHS_BACKWARDS__MIN_STR_LEN;
  size_t PROB_STRING_LENGTHS_BACKWARDS__MAX_STR_LEN;
  size_t PROB_STRING_LENGTHS_BACKWARDS__MIN_STR_CNT;
  size_t PROB_STRING_LENGTHS_BACKWARDS__MAX_STR_CNT;
  double PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_CHAR_INS_RATE;
  double PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_CHAR_DEL_RATE;
  double PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_CHAR_SUB_RATE;
  double PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_STR_SWAP_RATE;
  double PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_STR_DUP_RATE;
  double PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_STR_DEL_RATE;

  size_t PROB_LAST_INDEX_OF_ZERO__MIN_VEC_LEN;
  size_t PROB_LAST_INDEX_OF_ZERO__MAX_VEC_LEN;
  int PROB_LAST_INDEX_OF_ZERO__MIN_NUM;
  int PROB_LAST_INDEX_OF_ZERO__MAX_NUM;
  double PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_SWAP_RATE;
  double PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_DEL_RATE;
  double PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_INS_RATE;
  double PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_SUB_RATE;

  size_t PROB_COUNT_ODDS__MIN_VEC_LEN; 
  size_t PROB_COUNT_ODDS__MAX_VEC_LEN; 
  int PROB_COUNT_ODDS__MIN_NUM; 
  int PROB_COUNT_ODDS__MAX_NUM; 
  double PROB_COUNT_ODDS__MUTATION__PER_NUM_SWAP_RATE; 
  double PROB_COUNT_ODDS__MUTATION__PER_NUM_DEL_RATE; 
  double PROB_COUNT_ODDS__MUTATION__PER_NUM_INS_RATE; 
  double PROB_COUNT_ODDS__MUTATION__PER_NUM_SUB_RATE; 

  size_t PROB_MIRROR_IMAGE__MIN_VEC_LEN; 
  size_t PROB_MIRROR_IMAGE__MAX_VEC_LEN; 
  int PROB_MIRROR_IMAGE__MIN_NUM; 
  int PROB_MIRROR_IMAGE__MAX_NUM; 
  double PROB_MIRROR_IMAGE__MUTATION__PER_VEC_RANDOMIZE_VAL_RATE; 
  double PROB_MIRROR_IMAGE__MUTATION__PER_VEC_MIRROR_RATE; 
  double PROB_MIRROR_IMAGE__MUTATION__COPY_RATE; 
  double PROB_MIRROR_IMAGE__MUTATION__INS_RATE; 
  double PROB_MIRROR_IMAGE__MUTATION__DEL_RATE; 
  double PROB_MIRROR_IMAGE__MUTATION__PER_VEC_SHUFFLE_RATE; 

  size_t PROB_VECTORS_SUMMED__MIN_VEC_LEN; 
  size_t PROB_VECTORS_SUMMED__MAX_VEC_LEN; 
  int PROB_VECTORS_SUMMED__MIN_NUM; 
  int PROB_VECTORS_SUMMED__MAX_NUM; 
  double PROB_VECTORS_SUMMED__MUTATION__PER_NUM_SUB_RATE; 
  double PROB_VECTORS_SUMMED__MUTATION__COPY_RATE; 
  double PROB_VECTORS_SUMMED__MUTATION__INS_RATE; 
  double PROB_VECTORS_SUMMED__MUTATION__DEL_RATE; 

  int PROB_SUM_OF_SQUARES__MIN_NUM;
  int PROB_SUM_OF_SQUARES__MAX_NUM;
  double PROB_SUM_OF_SQUARES__MUTATION__NUM_MUT_RATE;

  double PROB_VECTOR_AVERAGE__EPSILON;
  size_t PROB_VECTOR_AVERAGE__MIN_VEC_LEN;
  size_t PROB_VECTOR_AVERAGE__MAX_VEC_LEN;
  double PROB_VECTOR_AVERAGE__MIN_NUM;
  double PROB_VECTOR_AVERAGE__MAX_NUM;
  double PROB_VECTOR_AVERAGE__MUTATION__INS_RATE;
  double PROB_VECTOR_AVERAGE__MUTATION__DEL_RATE;
  double PROB_VECTOR_AVERAGE__MUTATION__SUB_RATE;

  int PROB_MEDIAN__MIN_NUM;
  int PROB_MEDIAN__MAX_NUM;
  double PROB_MEDIAN__MUTATION__PER_NUM_COPY_RATE;
  double PROB_MEDIAN__MUTATION__PER_NUM_SUB_RATE;
  double PROB_MEDIAN__MUTATION__PER_NUM_SWAP_RATE;

  int PROB_SMALLEST__MIN_NUM;
  int PROB_SMALLEST__MAX_NUM;
  double PROB_SMALLEST__MUTATION__PER_NUM_SUB_RATE;
  double PROB_SMALLEST__MUTATION__PER_NUM_SWAP_RATE;

  int PROB_GRADE__MIN_NUM;
  int PROB_GRADE__MAX_NUM;
  double PROB_GRADE__MUTATION__PER_NUM_RANDOMIZE_RATE;
  double PROB_GRADE__MUTATION__PER_NUM_ADJUST_RATE;

  // - Data collection group -
  std::string DATA_DIRECTORY;
  size_t SUMMARY_STATS_INTERVAL;
  size_t SNAPSHOT_INTERVAL;
  size_t SOLUTION_SCREEN_INTERVAL;

  // Experiment variables
  bool setup;
  size_t update;

  size_t dominant_prog_id;
  size_t dominant_test_id;

  size_t PROGRAM_MAX_PASSES;
  size_t PROGRAM_EVALUATION_TESTCASE_CNT; ///< How many test cases are organisms evaluated on during evaluation?
  size_t NUM_COHORTS;

  size_t STATIC_COEVO__NUM_STATIC_TESTCASES;

  emp::Ptr<emp::Random> random;

  emp::Ptr<inst_lib_t> inst_lib;
  emp::Ptr<hardware_t> eval_hardware;
  size_t eval_time;

  size_t smallest_prog_sol_size;
  bool solution_found;
  size_t update_first_solution_found;

  emp::BitSet<TAG_WIDTH> call_tag;
  
  emp::Ptr<prog_world_t> prog_world;
  emp::Ptr<emp::Systematics<prog_org_t, prog_org_gen_t>> prog_genotypic_systematics;
  emp::Ptr<prog_taxon_t> mrca_taxa_ptr;
  size_t mrca_changes;
  
  emp::vector<std::function<double(prog_org_t &)>> lexicase_prog_fit_set;

  TagLGPMutator<TAG_WIDTH> prog_mutator;

  // Test worlds
  emp::Ptr<prob_NumberIO_world_t> prob_NumberIO_world;
  emp::Ptr<prob_SmallOrLarge_world_t> prob_SmallOrLarge_world;
  emp::Ptr<prob_ForLoopIndex_world_t> prob_ForLoopIndex_world;
  emp::Ptr<prob_CompareStringLengths_world_t> prob_CompareStringLengths_world;
  emp::Ptr<prob_DoubleLetters_world_t> prob_DoubleLetters_world;
  emp::Ptr<prob_CollatzNumbers_world_t> prob_CollatzNumbers_world;
  emp::Ptr<prob_ReplaceSpaceWithNewline_world_t> prob_ReplaceSpaceWithNewline_world;
  emp::Ptr<prob_StringDifferences_world_t> prob_StringDifferences_world;
  emp::Ptr<prob_EvenSquares_world_t> prob_EvenSquares_world;
  emp::Ptr<prob_WallisPi_world_t> prob_WallisPi_world;
  emp::Ptr<prob_StringLengthsBackwards_world_t> prob_StringLengthsBackwards_world;
  emp::Ptr<prob_LastIndexOfZero_world_t> prob_LastIndexOfZero_world;
  emp::Ptr<prob_VectorAverage_world_t> prob_VectorAverage_world;
  emp::Ptr<prob_CountOdds_world_t> prob_CountOdds_world;
  emp::Ptr<prob_MirrorImage_world_t> prob_MirrorImage_world;
  emp::Ptr<prob_SuperAnagrams_world_t> prob_SuperAnagrams_world;
  emp::Ptr<prob_SumOfSquares_world_t> prob_SumOfSquares_world;
  emp::Ptr<prob_VectorsSummed_world_t> prob_VectorsSummed_world;
  emp::Ptr<prob_XWordLines_world_t> prob_XWordLines_world;
  emp::Ptr<prob_PigLatin_world_t> prob_PigLatin_world;
  emp::Ptr<prob_NegativeToZero_world_t> prob_NegativeToZero_world;
  emp::Ptr<prob_ScrabbleScore_world_t> prob_ScrabbleScore_world;
  emp::Ptr<prob_Checksum_world_t> prob_Checksum_world;
  emp::Ptr<prob_Digits_world_t> prob_Digits_world;
  emp::Ptr<prob_Grade_world_t> prob_Grade_world;
  emp::Ptr<prob_Median_world_t> prob_Median_world;
  emp::Ptr<prob_Smallest_world_t> prob_Smallest_world;
  emp::Ptr<prob_Syllables_world_t> prob_Syllables_world;
  
  // Problem utilities
  ProblemUtilities_NumberIO prob_utils_NumberIO;
  ProblemUtilities_SmallOrLarge prob_utils_SmallOrLarge;
  ProblemUtilities_ForLoopIndex prob_utils_ForLoopIndex;
  ProblemUtilities_CompareStringLengths prob_utils_CompareStringLengths;
  ProblemUtilities_DoubleLetters prob_utils_DoubleLetters;
  ProblemUtilities_CollatzNumbers prob_utils_CollatzNumbers;
  ProblemUtilities_ReplaceSpaceWithNewline prob_utils_ReplaceSpaceWithNewline;
  ProblemUtilities_StringDifferences prob_utils_StringDifferences;
  ProblemUtilities_EvenSquares prob_utils_EvenSquares;
  ProblemUtilities_WallisPi prob_utils_WallisPi;
  ProblemUtilities_StringLengthsBackwards prob_utils_StringLengthsBackwards;
  ProblemUtilities_LastIndexOfZero prob_utils_LastIndexOfZero;
  ProblemUtilities_VectorAverage prob_utils_VectorAverage;
  ProblemUtilities_CountOdds prob_utils_CountOdds;
  ProblemUtilities_MirrorImage prob_utils_MirrorImage;
  ProblemUtilities_SuperAnagrams prob_utils_SuperAnagrams;
  ProblemUtilities_SumOfSquares prob_utils_SumOfSquares;
  ProblemUtilities_VectorsSummed prob_utils_VectorsSummed;
  ProblemUtilities_XWordLines prob_utils_XWordLines;
  ProblemUtilities_PigLatin prob_utils_PigLatin;
  ProblemUtilities_NegativeToZero prob_utils_NegativeToZero;
  ProblemUtilities_ScrabbleScore prob_utils_ScrabbleScore;
  ProblemUtilities_Checksum prob_utils_Checksum;
  ProblemUtilities_Digits prob_utils_Digits;
  ProblemUtilities_Grade prob_utils_Grade;
  ProblemUtilities_Median prob_utils_Median;
  ProblemUtilities_Smallest prob_utils_Smallest;
  ProblemUtilities_Syllables prob_utils_Syllables;

  Cohorts prog_cohorts;
  Cohorts test_cohorts;

  emp::Ptr<emp::DataFile> solution_file;
  emp::Ptr<emp::DataFile> prog_phen_diversity_file;

  struct TestResult {
    double score;
    bool pass;
    bool sub;
    TestResult(double sc=0, bool p=false, bool sb=false) : score(sc), pass(p), sub(sb) { ; }
  };

  /// StatsUtil is useful for managing target program/test during snapshots (gets
  /// captured in lambda).
  struct StatsUtil {
    size_t cur_progID;
    size_t cur_testID;

    emp::vector<TestResult> current_program__validation__test_results;
    double current_program__validation__total_score;
    size_t current_program__validation__total_passes;
    bool current_program__validation__is_solution;

    StatsUtil(size_t pID=0, size_t tID=0) : cur_progID(pID), cur_testID(tID) { ; }
  } stats_util;

  struct ProgramStats {
    // Generic stuff
    std::function<size_t(void)> get_id;
    std::function<double(void)> get_fitness;

    // Fitness evaluation stats
    std::function<double(void)> get_fitness_eval__total_score;          // - get_fitness_eval__total_score
    std::function<size_t(void)> get_fitness_eval__num_passes;           // - get_fitness_eval__num_passes
    std::function<size_t(void)> get_fitness_eval__num_fails;            // - get_fitness_eval__num_fails
    std::function<size_t(void)> get_fitness_eval__num_tests;            // - get_fitness_eval__num_tests
    std::function<std::string(void)> get_fitness_eval__passes_by_test;  // - get_fitness_eval__passes_by_test

    // program validation stats
    std::function<size_t(void)> get_validation_eval__num_passes;          // - get_validation_eval__num_passes;
    std::function<size_t(void)> get_validation_eval__num_tests;           // - get_validation_eval__num_tests
    std::function<std::string(void)> get_validation_eval__passes_by_test; // - get_validation_eval__passes_by_test
    std::function<double(void)> get_prog_behavioral_diversity;
    std::function<double(void)> get_prog_unique_behavioral_phenotypes;
    
    // program 'morphology' stats
    std::function<size_t(void)> get_program_len;
    std::function<std::string(void)> get_program;

  } program_stats;

  // Experiment signals
  emp::Signal<void(void)> do_evaluation_sig;
  emp::Signal<void(void)> do_selection_sig;
  emp::Signal<void(void)> do_update_sig;

  emp::Signal<void(void)> do_pop_snapshot_sig;
  
  emp::Signal<void(void)> end_setup_sig;
  emp::Signal<void(void)> on_destruction; ///< Triggered on experiment destruction

  // Program evaluation signals.
  emp::Signal<void(prog_org_t &)> begin_program_eval;
  emp::Signal<void(prog_org_t &)> end_program_eval;

  emp::Signal<void(prog_org_t &, emp::Ptr<TestOrg_Base>)> begin_program_test;
  emp::Signal<void(prog_org_t &, emp::Ptr<TestOrg_Base>)> do_program_test;
  emp::Signal<void(prog_org_t &, emp::Ptr<TestOrg_Base>)> end_program_test;

  emp::Signal<void(prog_org_t &)> do_program_advance;

  // Functions to be setup depending on experiment configuration (e.g., what problem we're solving)
  std::function<void(void)> UpdateTestCaseWorld;
  
  std::function<TestResult(prog_org_t &, size_t)> EvaluateWorldTest;                ///< Evaluate given program org on world test (specified by given ID). Return the test result.
  std::function<TestResult(prog_org_t &, TestOrg_Base &)> CalcProgramResultOnTest;  ///< Calculate the test result for a given program on a given test organism.
  
  std::function<void(prog_org_t &)> DoTestingSetValidation;     ///< Run program on full validation testing set.
  std::function<bool(prog_org_t &)> ScreenForSolution;          ///< Run program on validation testing set. Return true if program is a solution; false otherwise.
  
  std::function<test_org_phen_t&(size_t)> GetTestPhenotype;     ///< Utility function used to get test phenotype of given test (test type agnostic).
  std::function<void(void)> SetupTestMutation;                  ///< Test world configuration utility. To be defined by test setup.
  std::function<void(void)> SetupTestFitFun;                    ///< Test world configuration utility. To be defined by test setup.
  std::function<void(void)> SnapshotTests;                      ///< Snapshot test population.
  

  // Internal function signatures.
  void InitConfigs(const ProgramSynthesisConfig & config);

  void InitProgPop_Random();    ///< Randomly initialize the program population.
  
  void SetupHardware();         ///< Setup virtual hardware.
  void SetupEvaluation();       ///< Setup evaluation
  void SetupSelection();        ///< Setup selection (?)
  void SetupMutation();         ///< Setup mutation (?)
  void SetupFitFuns();
  void SetupDataCollection();   ///< Setup data collection

  void SetupProgramSelection(); ///< Setup program selection scheme
  void SetupProgramMutation();  ///< Setup program mutations
  void SetupProgramFitFun();
  void SetupProgramStats();

  void AddDefaultInstructions(const std::unordered_set<std::string> & includes);

  void SnapshotPrograms();

  void SetupProblem();
  void SetupProblem_NumberIO();
  void SetupProblem_SmallOrLarge();
  void SetupProblem_ForLoopIndex();
  void SetupProblem_CompareStringLengths();
  void SetupProblem_DoubleLetters();
  void SetupProblem_CollatzNumbers();
  void SetupProblem_ReplaceSpaceWithNewline();
  void SetupProblem_StringDifferences();
  void SetupProblem_EvenSquares();
  void SetupProblem_WallisPi();
  void SetupProblem_StringLengthsBackwards();
  void SetupProblem_LastIndexOfZero();
  void SetupProblem_VectorAverage();
  void SetupProblem_CountOdds();
  void SetupProblem_MirrorImage();
  void SetupProblem_SuperAnagrams();
  void SetupProblem_SumOfSquares();
  void SetupProblem_VectorsSummed();
  void SetupProblem_XWordLines();
  void SetupProblem_PigLatin();
  void SetupProblem_NegativeToZero();
  void SetupProblem_ScrabbleScore();
  void SetupProblem_Checksum();
  void SetupProblem_Digits();
  void SetupProblem_Grade();
  void SetupProblem_Median();
  void SetupProblem_Smallest();
  void SetupProblem_Syllables();

  // ---- Some useful world-type generic setup functions ----

  // Create a new world and do initial configuration.
  template<typename WORLD_ORG_TYPE>
  void NewTestCaseWorld(emp::Ptr<emp::World<WORLD_ORG_TYPE>> & w, emp::Random & rnd, const std::string & wname) {
    if (w != nullptr) { w.Delete(); }
    w.New(rnd, wname);
    w->SetPopStruct_Mixed(true);
  }

  template<typename WORLD_ORG_TYPE>
  void SetupTestCaseWorldUpdate(emp::Ptr<emp::World<WORLD_ORG_TYPE>> w) {
    // Tell experiment how to update the world (i.e., which world to update).
    switch (TRAINING_EXAMPLE_MODE) {
      case (size_t)TRAINING_EXAMPLE_MODE_TYPE::COEVOLUTION: {
        std::cout << "COEVOLUTION training example mode detected, configuring test world to update." << std::endl;
        UpdateTestCaseWorld = [w]() {
          // std::cout << "-COEVO world update" << std::endl;
          w->Update();
          w->ClearCache();
        };
        break;
      }
      case (size_t)TRAINING_EXAMPLE_MODE_TYPE::STATIC_GEN:
        std::cout << "STATIC-GEN training example mode detected, configuring test world to NOT update (reusing STATIC code)." << std::endl;
      case (size_t)TRAINING_EXAMPLE_MODE_TYPE::STATIC: {
        std::cout << "STATIC training example mode detected, configuring test world to NOT update." << std::endl;
        UpdateTestCaseWorld = [w]() {
          // std::cout << "-STATIC world update" << std::endl;
          w->ClearCache();
          // Reset test phenotypes
          for (size_t i = 0; i < w->GetSize(); ++i) {
            emp_assert(w->IsOccupied(i));
            TestOrg_Base & org = static_cast<TestOrg_Base&>(w->GetOrg(i));
            org.GetPhenotype().Reset(org.GetPhenotype().test_passes.size());
          }
        };
        break;
      }
      case (size_t)TRAINING_EXAMPLE_MODE_TYPE::STATIC_COEVO: {
        std::cout << "STATIC_COEVO training example mode detected, configuring test world to update." << std::endl;
        UpdateTestCaseWorld = [w]() {
          w->Update();
          w->ClearCache();
        };
        break;
      }
      case (size_t)TRAINING_EXAMPLE_MODE_TYPE::RANDOM: {
        std::cout << "RANDOM training example mode detected, configuring test world to NOT update. Instead, RANDOMIZE population." << std::endl;
        UpdateTestCaseWorld = [w]() {
          // std::cout << "-RANDOM world update" << std::endl;
          w->ClearCache();
          // Randomize population
          w->DoMutations();
          // Reset test phenotypes
          for (size_t i = 0; i < w->GetSize(); ++i) {
            emp_assert(w->IsOccupied(i));
            TestOrg_Base & org = static_cast<TestOrg_Base&>(w->GetOrg(i));
            org.GetPhenotype().Reset(org.GetPhenotype().test_passes.size());
            org.CalcOut();
          }
        };
        break;
      }
      default: {
        std::cout << "Unknown TRAINING_EXAMPLE_MODE (" << TRAINING_EXAMPLE_MODE << "). Exiting." << std::endl;
        exit(-1);
      }
    };
  }

  template<typename WORLD_ORG_TYPE>
  void SetupTestSelection(emp::Ptr<emp::World<WORLD_ORG_TYPE>> w, emp::vector<std::function<double(WORLD_ORG_TYPE &)>> & lexicase_fit_set) {
    std::cout << "Setting up test selection." << std::endl; 

    if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::STATIC_COEVO) {
      std::cout << "STATIC_COEVO mode, need to copy first STATIC_NUM (" << STATIC_COEVO__NUM_STATIC_TESTCASES << ") organisms." << std::endl;
      do_selection_sig.AddAction([this, w]() mutable {
        for (size_t i = 0; i < STATIC_COEVO__NUM_STATIC_TESTCASES; ++i) {
          w->DoBirth(w->GetGenomeAt(i), i);
          // note - may need to reset phenotype
        }
      });
    }

    switch (TEST_SELECTION_MODE) {
      case (size_t)SELECTION_TYPE::LEXICASE: {
        emp_assert(EVALUATION_MODE == (size_t)EVALUATION_TYPE::FULL);
        std::cout << "  Setting up test LEXICASE selection." << std::endl;
        // Setup lexicase selection.
        // - 1 lexicase function for every program organism.
        if (DISCRIMINATORY_LEXICASE_TESTS) {
          std::cout << "    Tests configured to be DISCRIMINATORY." << std::endl;
          for (size_t i = 0; i < PROG_POP_SIZE; ++i) {
            lexicase_fit_set.push_back([i](WORLD_ORG_TYPE & test_org) {
              TestOrg_Base & org = static_cast<TestOrg_Base&>(test_org);
              if (org.GetPhenotype().test_passes[i]) {
                return (double)org.GetPhenotype().num_fails;
              } else if (org.GetPhenotype().num_passes == 0) {
                return 0.5;
              } else {
                return 0.0;
              }        
            });
          }
        } else {
          for (size_t i = 0; i < PROG_POP_SIZE; ++i) {
            lexicase_fit_set.push_back([i](WORLD_ORG_TYPE & test_org) {
              TestOrg_Base & org = static_cast<TestOrg_Base&>(test_org);
              return (size_t)(!org.GetPhenotype().test_passes[i]);        // NOTE - test case fitness functions do not use gradient; use pass/fail.
            });
          }
        }
        // Add selection action.
        do_selection_sig.AddAction([this, w, lexicase_fit_set]() mutable { // todo - check that capture is working as expected
          emp::LexicaseSelect_NAIVE(*w, lexicase_fit_set, TEST_POP_SIZE, TEST_LEXICASE_MAX_FUNS);
        });
        break;
      }
      case (size_t)SELECTION_TYPE::COHORT_LEXICASE: {
        emp_assert(EVALUATION_MODE == (size_t)EVALUATION_TYPE::COHORT);
        std::cout << "  Setting up test COHORT LEXICASE selection." << std::endl;
        // Setup cohort lexicase.
        // - 1 lexicase function for every program cohort member.
        if (DISCRIMINATORY_LEXICASE_TESTS) {
          std::cout << "    Tests configured to be DISCRIMINATORY." << std::endl;
          for (size_t i = 0; i < PROG_COHORT_SIZE; ++i) {
            lexicase_fit_set.push_back([i](WORLD_ORG_TYPE & test_org) {
              TestOrg_Base & org = static_cast<TestOrg_Base&>(test_org);
              if (org.GetPhenotype().test_passes[i]) {
                return (double)org.GetPhenotype().num_fails;
              } else if (org.GetPhenotype().num_passes == 0) {
                return 0.5;
              } else {
                return 0.0;
              }
            });
          }
        } else {
          for (size_t i = 0; i < PROG_COHORT_SIZE; ++i) {
            lexicase_fit_set.push_back([i](WORLD_ORG_TYPE & test_org) {
              TestOrg_Base & org = static_cast<TestOrg_Base&>(test_org);
              return (size_t)(!org.GetPhenotype().test_passes[i]);
            });
          }
        }
        // Add selection action.
        emp_assert(TEST_COHORT_SIZE * test_cohorts.GetCohortCnt() == TEST_POP_SIZE);
        do_selection_sig.AddAction([this, w, lexicase_fit_set]() mutable { // todo - check that capture is working as expected
          // For each cohort, run selection.
          for (size_t cID = 0; cID < test_cohorts.GetCohortCnt(); ++cID) {
            emp::CohortLexicaseSelect_NAIVE(*w,
                                            lexicase_fit_set,
                                            test_cohorts.GetCohort(cID),
                                            TEST_COHORT_SIZE,
                                            TEST_COHORTLEXICASE_MAX_FUNS);
          }
        });
        break;
      }
      case (size_t)SELECTION_TYPE::TOURNAMENT: {
        std::cout << "  Setting up test TOURNAMENT selection." << std::endl;
        do_selection_sig.AddAction([this, w]() mutable {
          emp::TournamentSelect(*w, TEST_TOURNAMENT_SIZE, TEST_POP_SIZE);
        });
        break;
      }
      case (size_t)SELECTION_TYPE::DRIFT: {
        std::cout << "  Setting up test DRIFT selection." << std::endl;
        do_selection_sig.AddAction([this, w]() mutable {
          emp::RandomSelect(*w, TEST_POP_SIZE);
        });
        break;
      }
      default: {
        std::cout << "Unknown TEST_SELECTION_MODE (" << TEST_SELECTION_MODE << "). Exiting." << std::endl;
        exit(-1);
      }
    }

    if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::STATIC_COEVO) {
      std::cout << "STATIC_COEVO mode, need to fix population size post-selection." << std::endl;
      do_selection_sig.AddAction([this, w]() mutable {
        // std::cout << "Population size pre-resize? " << w->GetSize() << std::endl;
        // Shuffle the non-static test cases to be fair to late cohorts
        for (size_t i = STATIC_COEVO__NUM_STATIC_TESTCASES; i < w->GetSize(); ++i) {
          w->Swap(i, w->GetRandom().GetUInt(STATIC_COEVO__NUM_STATIC_TESTCASES, w->GetSize()));
        }
        w->Resize(TEST_POP_SIZE);
        // std::cout << "Population size? " << w->GetSize() << std::endl;
      });
    }

  }

  template<typename WORLD_ORG_TYPE>
  void SetupTestSystematics(emp::Ptr<emp::World<WORLD_ORG_TYPE>> w, const std::function<void(std::ostream &, const typename emp::World<WORLD_ORG_TYPE>::genome_t &)> & print_test) {
    using w_t = emp::World<WORLD_ORG_TYPE>;
    std::function<typename w_t::genome_t(const WORLD_ORG_TYPE & org)> sys_get_genome = [](const WORLD_ORG_TYPE & org) { return org.GetGenome(); };
    auto sys = w->AddSystematics(sys_get_genome, true, true, false, true, "test_genotype");
    sys->AddEvolutionaryDistinctivenessDataNode();
    sys->AddPairwiseDistanceDataNode();
    sys->AddPhylogeneticDiversityDataNode();
    auto & sys_file = w->SetupSystematicsFile("test_genotype", DATA_DIRECTORY + "/test_gen_sys.csv", false);
    sys_file.SetTimingRepeat(SUMMARY_STATS_INTERVAL); 
    sys_file.AddStats(*sys->GetDataNode("evolutionary_distinctiveness") , "evolutionary_distinctiveness", "evolutionary distinctiveness for a single update", true, true);
    sys_file.AddStats(*sys->GetDataNode("pairwise_distances"), "pairwise_distance", "pairwise distance for a single update", true, true);
    sys_file.AddCurrent(*sys->GetDataNode("phylogenetic_diversity"), "current_phylogenetic_diversity", "current phylogenetic_diversity", true, true);
    sys_file.template AddFun<size_t>([sys]() { return sys->GetTreeSize(); }, "tree_size", "Phylogenetic tree size");
    sys_file.PrintHeaderKeys();
    using to_taxon_t = typename emp::Systematics<WORLD_ORG_TYPE, typename w_t::genome_t>::taxon_t;
    sys->AddSnapshotFun([print_test](const to_taxon_t & t) {
      std::ostringstream stream;
      print_test(stream, t.GetInfo());
      return stream.str();
    }, "test", "Test (input)");
    do_pop_snapshot_sig.AddAction([this, sys]() mutable {
      sys->Snapshot(DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate()) + "/test_phylogeny_" + emp::to_string(prog_world->GetUpdate()) + ".csv");
    });
  }

  void OnPlacement_ActiveTestCaseWorld(const std::function<void(size_t)> & fun) {
      if (prob_NumberIO_world != nullptr) prob_NumberIO_world->OnPlacement(fun);
      else if (prob_SmallOrLarge_world != nullptr) prob_SmallOrLarge_world->OnPlacement(fun);
      else if (prob_ForLoopIndex_world != nullptr) prob_ForLoopIndex_world->OnPlacement(fun);
      else if (prob_CompareStringLengths_world != nullptr) prob_CompareStringLengths_world->OnPlacement(fun);
      else if (prob_DoubleLetters_world != nullptr) prob_DoubleLetters_world->OnPlacement(fun);
      else if (prob_CollatzNumbers_world != nullptr) prob_CollatzNumbers_world->OnPlacement(fun);
      else if (prob_ReplaceSpaceWithNewline_world != nullptr) prob_ReplaceSpaceWithNewline_world->OnPlacement(fun);
      else if (prob_StringDifferences_world != nullptr) prob_StringDifferences_world->OnPlacement(fun);
      else if (prob_EvenSquares_world != nullptr) prob_EvenSquares_world->OnPlacement(fun);
      else if (prob_WallisPi_world != nullptr) prob_WallisPi_world->OnPlacement(fun);
      else if (prob_StringLengthsBackwards_world != nullptr) prob_StringLengthsBackwards_world->OnPlacement(fun);
      else if (prob_LastIndexOfZero_world != nullptr) prob_LastIndexOfZero_world->OnPlacement(fun);
      else if (prob_VectorAverage_world != nullptr) prob_VectorAverage_world->OnPlacement(fun);
      else if (prob_CountOdds_world != nullptr) prob_CountOdds_world->OnPlacement(fun);
      else if (prob_MirrorImage_world != nullptr) prob_MirrorImage_world->OnPlacement(fun);
      else if (prob_SuperAnagrams_world != nullptr) prob_SuperAnagrams_world->OnPlacement(fun);
      else if (prob_SumOfSquares_world != nullptr) prob_SumOfSquares_world->OnPlacement(fun);
      else if (prob_VectorsSummed_world != nullptr) prob_VectorsSummed_world->OnPlacement(fun);
      else if (prob_XWordLines_world != nullptr) prob_XWordLines_world->OnPlacement(fun);
      else if (prob_PigLatin_world != nullptr) prob_PigLatin_world->OnPlacement(fun);
      else if (prob_NegativeToZero_world != nullptr) prob_NegativeToZero_world->OnPlacement(fun);
      else if (prob_ScrabbleScore_world != nullptr) prob_ScrabbleScore_world->OnPlacement(fun);
      else if (prob_Checksum_world != nullptr) prob_Checksum_world->OnPlacement(fun);
      else if (prob_Digits_world != nullptr) prob_Digits_world->OnPlacement(fun);
      else if (prob_Grade_world != nullptr) prob_Grade_world->OnPlacement(fun);
      else if (prob_Median_world != nullptr) prob_Median_world->OnPlacement(fun);
      else if (prob_Smallest_world != nullptr) prob_Smallest_world->OnPlacement(fun);
      else if (prob_Syllables_world != nullptr) prob_Syllables_world->OnPlacement(fun);
      else { std::cout << "AHH! More than one test case world has been created. Exiting." << std::endl; exit(-1); }
  }
  
  /// Test set is test set to use if doing static initialization
  /// rand gen fun is used to generate random genome
  template<typename WORLD_ORG_TYPE, typename TEST_IN_TYPE, typename TEST_OUT_TYPE>
  void SetupTestCasePop_Init(emp::Ptr<emp::World<WORLD_ORG_TYPE>> w, TestCaseSet<TEST_IN_TYPE, TEST_OUT_TYPE> & test_set,
                             const std::function<typename emp::World<WORLD_ORG_TYPE>::genome_t(void)> & gen_rand_test) {
    // Configure how population should be initialized -- TODO - maybe move this into functor(?)
    
    if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::STATIC) {
      TEST_POP_SIZE = test_set.GetSize();
      std::cout << "In STATIC training example mode, adjusting TEST_POP_SIZE to: " << TEST_POP_SIZE << std::endl;
      end_setup_sig.AddAction([this, w, test_set]() {                     // TODO - test that this is actually working!
        InitTestCasePop_TrainingSet(w, test_set);
      });
    } else if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::STATIC_GEN) {
      std::cout << "In STATIC_GEN training example mode bolstering training examples to match TEST_POP_SIZE (" << TEST_POP_SIZE << ")" << std::endl;
      end_setup_sig.AddAction([this, w, test_set, gen_rand_test]() {       // TODO - test that this is actually working!
        if (PROBLEM == "sum-of-squares") {
          InitTestCasePop_TrainingSetBolstered(w, test_set, gen_rand_test, false);
        } else {
          InitTestCasePop_TrainingSetBolstered(w, test_set, gen_rand_test, true);
        }
      });

    } else if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::STATIC_COEVO) {
      std::cout << "In STATIC_COEVO training example mode, bolstering training examples to match pop size." << std::endl;
      STATIC_COEVO__NUM_STATIC_TESTCASES = test_set.GetSize(); // Number of test cases to never mutate
      end_setup_sig.AddAction([this, w, test_set, gen_rand_test]() {       // TODO - test that this is actually working!
        if (PROBLEM == "sum-of-squares") {
          InitTestCasePop_TrainingSetBolstered(w, test_set, gen_rand_test, false);
        } else {
          InitTestCasePop_TrainingSetBolstered(w, test_set, gen_rand_test, true);
        }
      });
    } else {
      std::cout << "Generating training example population randomly." << std::endl;
      end_setup_sig.AddAction([this, w, test_set, gen_rand_test]() {      // TODO - test that this is actually working!
        InitTestCasePop_Random(w, gen_rand_test);
      });
    }
  }

  // Initialize given world's population with training examples in given test case set.
  template<typename WORLD_ORG_TYPE, typename TEST_IN_TYPE, typename TEST_OUT_TYPE>
  void InitTestCasePop_TrainingSet(emp::Ptr<emp::World<WORLD_ORG_TYPE>> w, const TestCaseSet<TEST_IN_TYPE, TEST_OUT_TYPE> & test_set) {
    std::cout << "Initializing test case population from a training set." << std::endl;
    for (size_t i = 0; i < test_set.GetSize(); ++i) {
      w->Inject(test_set.GetInput(i), 1);
    }
  }

  // Initialize given world's population with training examples in given test case set.
  template<typename WORLD_ORG_TYPE, typename TEST_IN_TYPE, typename TEST_OUT_TYPE>
  void InitTestCasePop_TrainingSetBolstered(emp::Ptr<emp::World<WORLD_ORG_TYPE>> w, 
                                            const TestCaseSet<TEST_IN_TYPE, TEST_OUT_TYPE> & test_set,
                                            const std::function<typename emp::World<WORLD_ORG_TYPE>::genome_t(void)> & gen_rand,
                                            bool guarantee_unique=true) {
    std::cout << "Initializing test case population from a training set then bolstering to TEST_POP_SIZE." << std::endl;
    std::cout << "  Test set size = " << test_set.GetSize() << std::endl;

    emp::vector<size_t> tids;
    for (size_t i = 0; i < test_set.GetSize(); ++i) { tids.emplace_back(i); }
    emp::Shuffle(*random, tids);
    std::cout << "Test IDS: ";
    for (size_t i = 0; i < test_set.GetSize(); ++i) { std::cout << tids[i] << ";"; }
    std::cout << std::endl;

    while (w->GetSize() < TEST_POP_SIZE) {
      if (w->GetSize() < test_set.GetSize()) {
        // std::cout << "Injecting test id " << w->GetSize() << std::endl;
        w->Inject(test_set.GetInput(tids[w->GetSize()]), 1);
      } else {
        // std::cout << "Injecting randomly generated input " << std::endl;
        bool dup = true;  // Assume a dup, prove it's not.
        typename emp::World<WORLD_ORG_TYPE>::genome_t rand_genome = gen_rand();
        while (dup && guarantee_unique) {
          dup = false; // We have no reason to believe random is actually a duplicate.
          for (size_t i = 0; i < w->GetSize(); ++i) {
            if (w->GetGenomeAt(i) == rand_genome) {
              dup = true;
              rand_genome = gen_rand();
              break;
            }
          }
        }
        w->Inject(rand_genome, 1);
      }
    }
    std::cout << "  World size = " << w->GetSize() << std::endl;
    
    // ensure uniqueness
    bool unique = true;
    for (size_t i = 0; i < w->GetSize() && unique; ++i) {
      for (size_t k = i+1; k < w->GetSize() && unique; ++k) {
        // std::cout << "Comparing " << i << " " << k << std::endl;
        if (w->GetGenomeAt(i) == w->GetGenomeAt(k)) {
          unique = false;
        }
      }
    }

    if (guarantee_unique) { emp_assert(unique); }
    std::cout << "  Unique pop? " << unique << std::endl;

  }

  // Initialize given world's population randomly.
  template<typename WORLD_ORG_TYPE>
  void InitTestCasePop_Random(emp::Ptr<emp::World<WORLD_ORG_TYPE>> w, const std::function<typename emp::World<WORLD_ORG_TYPE>::genome_t(void)> & fun) {
    std::cout << "Initializing test case population randomly." << std::endl;
    for (size_t i = 0; i < TEST_POP_SIZE; ++i) {
      w->Inject(fun(), 1);
    }
  }

  // ---- Problem-specific instruction signatures ----
  // -- NumberIO --
  void Inst_LoadInt_NumberIO(hardware_t & hw, const inst_t & inst);
  void Inst_LoadDouble_NumberIO(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_NumberIO(hardware_t & hw, const inst_t & inst); 
  // -- SmallOrLarge --
  void Inst_LoadInt_SmallOrLarge(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitSmall_SmallOrLarge(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitLarge_SmallOrLarge(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNone_SmallOrLarge(hardware_t & hw, const inst_t & inst);
  // -- ForLoopIndex --
  void Inst_LoadStart_ForLoopIndex(hardware_t & hw, const inst_t & inst);
  void Inst_LoadEnd_ForLoopIndex(hardware_t & hw, const inst_t & inst);
  void Inst_LoadStep_ForLoopIndex(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_ForLoopIndex(hardware_t & hw, const inst_t & inst);
  // -- CompareStringLengths --
  void Inst_LoadStr1_CompareStringLengths(hardware_t & hw, const inst_t & inst);
  void Inst_LoadStr2_CompareStringLengths(hardware_t & hw, const inst_t & inst);
  void Inst_LoadStr3_CompareStringLengths(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitTrue_CompareStringLengths(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitFalse_CompareStringLengths(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitVal_CompareStringLengths(hardware_t & hw, const inst_t & inst);
  // -- CollatzNumbers ---
  void Inst_LoadNum_CollatzNumbers(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_CollatzNumbers(hardware_t & hw, const inst_t & inst);
  // -- StringLengthsBackwards --
  void Inst_LoadStrVec_StringLengthsBackwards(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitVal_StringLengthsBackwards(hardware_t & hw, const inst_t & inst);
  // -- LastIndexOfZero --
  void Inst_LoadVec_LastIndexOfZero(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_LastIndexOfZero(hardware_t & hw, const inst_t & inst);
  // -- CountOdds --
  void Inst_LoadVec_CountOdds(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_CountOdds(hardware_t & hw, const inst_t & inst);
  // -- MirrorImage --
  void Inst_LoadVec1_MirrorImage(hardware_t & hw, const inst_t & inst);
  void Inst_LoadVec2_MirrorImage(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitVal_MirrorImage(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitTrue_MirrorImage(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitFalse_MirrorImage(hardware_t & hw, const inst_t & inst);
  // -- VectorsSummed --
  void Inst_LoadVec1_VectorsSummed(hardware_t & hw, const inst_t & inst);
  void Inst_LoadVec2_VectorsSummed(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitVec_VectorsSummed(hardware_t & hw, const inst_t & inst);
  // -- Sum of Squares --
  void Inst_LoadNum_SumOfSquares(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_SumOfSquares(hardware_t & hw, const inst_t & inst);
  // -- VectorAverage --
  void Inst_LoadVec_VectorAverage(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_VectorAverage(hardware_t & hw, const inst_t & inst);
  // -- Median --
  void Inst_LoadNum1_Median(hardware_t & hw, const inst_t & inst);
  void Inst_LoadNum2_Median(hardware_t & hw, const inst_t & inst);
  void Inst_LoadNum3_Median(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_Median(hardware_t & hw, const inst_t & inst);
  // -- Smallest --
  void Inst_LoadNum1_Smallest(hardware_t & hw, const inst_t & inst);
  void Inst_LoadNum2_Smallest(hardware_t & hw, const inst_t & inst);
  void Inst_LoadNum3_Smallest(hardware_t & hw, const inst_t & inst);
  void Inst_LoadNum4_Smallest(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_Smallest(hardware_t & hw, const inst_t & inst);
  // -- Grade --
  void Inst_LoadThreshA_Grade(hardware_t & hw, const inst_t & inst);
  void Inst_LoadThreshB_Grade(hardware_t & hw, const inst_t & inst);
  void Inst_LoadThreshC_Grade(hardware_t & hw, const inst_t & inst);
  void Inst_LoadThreshD_Grade(hardware_t & hw, const inst_t & inst);
  void Inst_LoadGrade_Grade(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitA_Grade(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitB_Grade(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitC_Grade(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitD_Grade(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitF_Grade(hardware_t & hw, const inst_t & inst);

public:
  ProgramSynthesisExperiment() 
    : setup(false), update(0), solution_found(false)
  {
    std::cout << "Problem info:" << std::endl;
    for (const auto & info : problems) {
      std::cout << "  - Problem name: " << info.first << std::endl;
      std::cout << "    - Training examples file: " << info.second.GetTrainingSetFilename() << std::endl;
      std::cout << "    - Testing examples file: " << info.second.GetTestingSetFilename() << std::endl;
    }
  }

  ~ProgramSynthesisExperiment() {
    if (setup) {
      solution_file.Delete();
      prog_phen_diversity_file.Delete();
      eval_hardware.Delete();
      inst_lib.Delete();
      prog_world.Delete();
      // prog_genotypic_systematics.Delete();

      if (prob_NumberIO_world != nullptr) prob_NumberIO_world.Delete();
      if (prob_SmallOrLarge_world != nullptr) prob_SmallOrLarge_world.Delete();
      if (prob_ForLoopIndex_world != nullptr) prob_ForLoopIndex_world.Delete();
      if (prob_CompareStringLengths_world != nullptr) prob_CompareStringLengths_world.Delete();
      if (prob_DoubleLetters_world != nullptr) prob_DoubleLetters_world.Delete();
      if (prob_CollatzNumbers_world != nullptr) prob_CollatzNumbers_world.Delete();
      if (prob_ReplaceSpaceWithNewline_world != nullptr) prob_ReplaceSpaceWithNewline_world.Delete();
      if (prob_StringDifferences_world != nullptr) prob_StringDifferences_world.Delete();
      if (prob_EvenSquares_world != nullptr) prob_EvenSquares_world.Delete();
      if (prob_WallisPi_world != nullptr) prob_WallisPi_world.Delete();
      if (prob_StringLengthsBackwards_world != nullptr) prob_StringLengthsBackwards_world.Delete();
      if (prob_LastIndexOfZero_world != nullptr) prob_LastIndexOfZero_world.Delete();
      if (prob_VectorAverage_world != nullptr) prob_VectorAverage_world.Delete();
      if (prob_CountOdds_world != nullptr) prob_CountOdds_world.Delete();
      if (prob_MirrorImage_world != nullptr) prob_MirrorImage_world.Delete();
      if (prob_SuperAnagrams_world != nullptr) prob_SuperAnagrams_world.Delete();
      if (prob_SumOfSquares_world != nullptr) prob_SumOfSquares_world.Delete();
      if (prob_VectorsSummed_world != nullptr) prob_VectorsSummed_world.Delete();
      if (prob_XWordLines_world != nullptr) prob_XWordLines_world.Delete();
      if (prob_PigLatin_world != nullptr) prob_PigLatin_world.Delete();
      if (prob_NegativeToZero_world != nullptr) prob_NegativeToZero_world.Delete();
      if (prob_ScrabbleScore_world != nullptr) prob_ScrabbleScore_world.Delete();
      if (prob_Checksum_world != nullptr) prob_Checksum_world.Delete();
      if (prob_Digits_world != nullptr) prob_Digits_world.Delete();
      if (prob_Grade_world != nullptr) prob_Grade_world.Delete();
      if (prob_Median_world != nullptr) prob_Median_world.Delete();
      if (prob_Smallest_world != nullptr) prob_Smallest_world.Delete();
      if (prob_Syllables_world != nullptr) prob_Syllables_world.Delete();

      random.Delete();
    }
  }

  /// Configure the experiment.
  void Setup(const ProgramSynthesisConfig & config);

  /// Run the experiment start->finish.
  void Run();

  /// Progress experiment by a single time step (generation):
  /// - (1) evaluate population(s)
  /// - (2) select individuals for reproduction
  /// - (3) update world(s)
  void RunStep();

  /// -- Accessors --
  emp::Random & GetRandom() { return *random; }
  prog_world_t & GetProgramWorld() { return *prog_world; }

};

/// ================ Public facing implementations ================

/// Configure the experiment.
void ProgramSynthesisExperiment::Setup(const ProgramSynthesisConfig & config) {
  std::cout << "Running ProgramSynthesisExperiment setup." << std::endl;
  emp_assert(setup == false, "Can only run setup once because I'm lazy.");
  // Localize experiment configuration.
  InitConfigs(config);
  // Initialization/cleanup depending on whether or not this is the first call to setup.
  if (!setup) { // First call to Setup.
    // Create a random number generator.
    random = emp::NewPtr<emp::Random>(SEED);
    // Create the program world.
    prog_world = emp::NewPtr<prog_world_t>(*random, "program world");
  } else {
    // Fail!
    std::cout << "Setup only allowed once per experiment! (because I'm too lazy to handle multiple setups) Exiting." << std::endl;
    exit(-1);
  }

  smallest_prog_sol_size = MAX_PROG_SIZE + 1;
  solution_found = false;
  update_first_solution_found = GENERATIONS + 1;
  mrca_taxa_ptr = nullptr;
  mrca_changes = 0;

  // Configure the program world.
  prog_world->SetPopStruct_Mixed(true);
  
  // Configure how programs should be initialized.
  end_setup_sig.AddAction([this]() {
    // Initialize program population.
    InitProgPop_Random();
    std::cout << "Program population size=" << prog_world->GetSize() << std::endl;
  });

  // Configure On Update signal.
  do_update_sig.AddAction([this]() {
    std::cout << "Update: " << update << "; ";
    std::cout << "best program score: " << prog_world->CalcFitnessID(dominant_prog_id) << "; ";
    std::cout << "solution found? " << solution_found << "; ";
    std::cout << "smallest solution? " << smallest_prog_sol_size << std::endl;

    if (update % SNAPSHOT_INTERVAL == 0 || update_first_solution_found == update || update == GENERATIONS) do_pop_snapshot_sig.Trigger();

    if (update_first_solution_found == update && update % SUMMARY_STATS_INTERVAL != 0) {
      prog_world->GetFile(DATA_DIRECTORY + "/prog_gen_sys.csv").Update(); // Update the program systematics files
    } 

    prog_world->Update();
    prog_world->ClearCache();

    UpdateTestCaseWorld();
  });

  // Setup the virtual hardware used to evaluate programs.
  std::cout << "==== EXPERIMENT SETUP => evaluation hardware ====" << std::endl;
  SetupHardware(); 
  
  // Setup problem that we're evolving programs to solve. The particular problem
  // we setup depends on experiment configuration.
  std::cout << "==== EXPERIMENT SETUP => problem ====" << std::endl;
  SetupProblem();  //...many many todos embedded in this one...

  // Setup program (& test) evaluation.
  std::cout << "==== EXPERIMENT SETUP => evaluation ====" << std::endl;
  SetupEvaluation(); 
  
  // Setup program (& test) selection.
  std::cout << "==== EXPERIMENT SETUP => selection ====" << std::endl;
  SetupSelection();

  // Setup program (& test) mutation.
  std::cout << "==== EXPERIMENT SETUP => mutation ====" << std::endl;
  SetupMutation();

  // Setup program (& test) fitness calculations.
  std::cout << "==== EXPERIMENT SETUP => fitness functions ====" << std::endl;
  SetupFitFuns();

  #ifndef EMSCRIPTEN
  // If we're not compiling to javascript, setup data collection for both programs
  // and tests.
  std::cout << "==== EXPERIMENT SETUP => data collection ====" << std::endl;
  SetupDataCollection();
  #endif

  // Signal that setup is done. This will trigger any post-setup things that need
  // to run before doing evolution (e.g., population initialization).
  std::cout << "==== EXPERIMENT SETUP => triggering end setup signal ====" << std::endl;
  end_setup_sig.Trigger();

  // Flag setup as done.
  setup = true;

  // TODO - assert that only one test world ptr is not nullptr

  std::cout << "==== EXPERIMENT SETUP => DONE! ====" << std::endl;
}

/// Run the experiment start->finish [update=0 : update=config.GENERATIONS].
void ProgramSynthesisExperiment::Run() {
  // For each generation, advance 'time' by one step.
  for (update = 0; update <= GENERATIONS; ++update) {
    RunStep();
  }
}

/// Run a single step of the experiment
void ProgramSynthesisExperiment::RunStep() {
  // std::cout << "-- Doing Evaluation --" << std::endl;
  do_evaluation_sig.Trigger();  // (1) Evaluate all members of program (& test) population(s).
  // std::cout << "-- Doing Selection --" << std::endl;
  do_selection_sig.Trigger();   // (2) Select who gets to reproduce!
  // std::cout << "-- Doing Update --" << std::endl;
  do_update_sig.Trigger();      // (3) Run update on relevant worlds (population turnover, etc).
}

// ================ Internal function implementations ================
/// Localize configs.
void ProgramSynthesisExperiment::InitConfigs(const ProgramSynthesisConfig & config) {
  SEED = config.SEED();
  GENERATIONS = config.GENERATIONS();
  PROG_POP_SIZE = config.PROG_POP_SIZE();
  TEST_POP_SIZE = config.TEST_POP_SIZE();
  EVALUATION_MODE = config.EVALUATION_MODE();
  PROG_COHORT_SIZE = config.PROG_COHORT_SIZE();
  TEST_COHORT_SIZE = config.TEST_COHORT_SIZE();
  TRAINING_EXAMPLE_MODE = config.TRAINING_EXAMPLE_MODE();
  PROBLEM = config.PROBLEM();
  BENCHMARK_DATA_DIR = config.BENCHMARK_DATA_DIR();

  // -- Selection settings --
  PROG_SELECTION_MODE = config.PROG_SELECTION_MODE();
  TEST_SELECTION_MODE = config.TEST_SELECTION_MODE();
  PROG_LEXICASE_MAX_FUNS = config.PROG_LEXICASE_MAX_FUNS();
  PROG_COHORTLEXICASE_MAX_FUNS = config.PROG_COHORTLEXICASE_MAX_FUNS();
  TEST_LEXICASE_MAX_FUNS = config.TEST_LEXICASE_MAX_FUNS();
  TEST_COHORTLEXICASE_MAX_FUNS = config.TEST_COHORTLEXICASE_MAX_FUNS();
  PROG_TOURNAMENT_SIZE = config.PROG_TOURNAMENT_SIZE();
  TEST_TOURNAMENT_SIZE = config.TEST_TOURNAMENT_SIZE();
  DISCRIMINATORY_LEXICASE_TESTS = config.DISCRIMINATORY_LEXICASE_TESTS();

  // -- Hardware settings --
  MIN_TAG_SPECIFICITY = config.MIN_TAG_SPECIFICITY();
  MAX_CALL_DEPTH = config.MAX_CALL_DEPTH();

  // -- Program settings --
  MIN_PROG_SIZE = config.MIN_PROG_SIZE();
  MAX_PROG_SIZE = config.MAX_PROG_SIZE();
  PROG_EVAL_TIME  = config.PROG_EVAL_TIME();
  PROG_MUT__PER_BIT_FLIP = config.PROG_MUT__PER_BIT_FLIP();
  PROG_MUT__PER_INST_SUB = config.PROG_MUT__PER_INST_SUB();
  PROG_MUT__PER_INST_INS = config.PROG_MUT__PER_INST_INS();
  PROG_MUT__PER_INST_DEL = config.PROG_MUT__PER_INST_DEL();
  PROG_MUT__PER_PROG_SLIP = config.PROG_MUT__PER_PROG_SLIP();
  PROG_MUT__PER_MOD_DUP = config.PROG_MUT__PER_MOD_DUP();
  PROG_MUT__PER_MOD_DEL = config.PROG_MUT__PER_MOD_DEL();

  // -- Number IO settings --
  PROB_NUMBER_IO__DOUBLE_MIN = config.PROB_NUMBER_IO__DOUBLE_MIN();
  PROB_NUMBER_IO__DOUBLE_MAX = config.PROB_NUMBER_IO__DOUBLE_MAX();
  PROB_NUMBER_IO__INT_MIN = config.PROB_NUMBER_IO__INT_MIN();
  PROB_NUMBER_IO__INT_MAX = config.PROB_NUMBER_IO__INT_MAX();
  PROB_NUMBER_IO__MUTATION__PER_INT_RATE = config.PROB_NUMBER_IO__MUTATION__PER_INT_RATE();
  PROB_NUMBER_IO__MUTATION__PER_DOUBLE_RATE = config.PROB_NUMBER_IO__MUTATION__PER_DOUBLE_RATE();

  PROB_SMALL_OR_LARGE__INT_MIN = config.PROB_SMALL_OR_LARGE__INT_MIN();
  PROB_SMALL_OR_LARGE__INT_MAX = config.PROB_SMALL_OR_LARGE__INT_MAX();
  PROB_SMALL_OR_LARGE__MUTATION__PER_INT_RATE = config.PROB_SMALL_OR_LARGE__MUTATION__PER_INT_RATE();

  PROB_FOR_LOOP_INDEX__START_END_MIN = config.PROB_FOR_LOOP_INDEX__START_END_MIN();
  PROB_FOR_LOOP_INDEX__START_END_MAX = config.PROB_FOR_LOOP_INDEX__START_END_MAX();
  PROB_FOR_LOOP_INDEX__STEP_MIN = config.PROB_FOR_LOOP_INDEX__STEP_MIN();
  PROB_FOR_LOOP_INDEX__STEP_MAX = config.PROB_FOR_LOOP_INDEX__STEP_MAX();
  PROB_FOR_LOOP_INDEX__MUTATION__MUT_RATE = config.PROB_FOR_LOOP_INDEX__MUTATION__MUT_RATE();
  PROB_FOR_LOOP_INDEX__PROMISE_MULTISTEP_TESTCASES = config.PROB_FOR_LOOP_INDEX__PROMISE_MULTISTEP_TESTCASES();

  PROB_COMPARE_STRING_LENGTHS__MIN_STR_LEN = config.PROB_COMPARE_STRING_LENGTHS__MIN_STR_LEN();
  PROB_COMPARE_STRING_LENGTHS__MAX_STR_LEN = config.PROB_COMPARE_STRING_LENGTHS__MAX_STR_LEN();
  PROB_COMPARE_STRING_LENGTHS__PER_SITE_INS_RATE = config.PROB_COMPARE_STRING_LENGTHS__PER_SITE_INS_RATE();
  PROB_COMPARE_STRING_LENGTHS__PER_SITE_DEL_RATE = config.PROB_COMPARE_STRING_LENGTHS__PER_SITE_DEL_RATE();
  PROB_COMPARE_STRING_LENGTHS__PER_SITE_SUB_RATE = config.PROB_COMPARE_STRING_LENGTHS__PER_SITE_SUB_RATE();
  PROB_COMPARE_STRING_LENGTHS__PER_STR_SWAP_RATE = config.PROB_COMPARE_STRING_LENGTHS__PER_STR_SWAP_RATE();

  PROB_COLLATZ_NUMBERS__MIN_NUM = config.PROB_COLLATZ_NUMBERS__MIN_NUM();
  PROB_COLLATZ_NUMBERS__MAX_NUM = config.PROB_COLLATZ_NUMBERS__MAX_NUM();
  PROB_COLLATZ_NUMBERS__MUTATION__PER_NUM_SUB_RATE = config.PROB_COLLATZ_NUMBERS__MUTATION__PER_NUM_SUB_RATE();

  PROB_STRING_LENGTHS_BACKWARDS__MIN_STR_LEN = config.PROB_STRING_LENGTHS_BACKWARDS__MIN_STR_LEN();
  PROB_STRING_LENGTHS_BACKWARDS__MAX_STR_LEN = config.PROB_STRING_LENGTHS_BACKWARDS__MAX_STR_LEN();
  PROB_STRING_LENGTHS_BACKWARDS__MIN_STR_CNT = config.PROB_STRING_LENGTHS_BACKWARDS__MIN_STR_CNT();
  PROB_STRING_LENGTHS_BACKWARDS__MAX_STR_CNT = config.PROB_STRING_LENGTHS_BACKWARDS__MAX_STR_CNT();
  PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_CHAR_INS_RATE = config.PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_CHAR_INS_RATE();
  PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_CHAR_DEL_RATE = config.PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_CHAR_DEL_RATE();
  PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_CHAR_SUB_RATE = config.PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_CHAR_SUB_RATE();
  PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_STR_SWAP_RATE = config.PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_STR_SWAP_RATE();
  PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_STR_DUP_RATE = config.PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_STR_DUP_RATE();
  PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_STR_DEL_RATE = config.PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_STR_DEL_RATE();

  PROB_LAST_INDEX_OF_ZERO__MIN_VEC_LEN = config.PROB_LAST_INDEX_OF_ZERO__MIN_VEC_LEN();
  PROB_LAST_INDEX_OF_ZERO__MAX_VEC_LEN = config.PROB_LAST_INDEX_OF_ZERO__MAX_VEC_LEN();
  PROB_LAST_INDEX_OF_ZERO__MIN_NUM = config.PROB_LAST_INDEX_OF_ZERO__MIN_NUM();
  PROB_LAST_INDEX_OF_ZERO__MAX_NUM = config.PROB_LAST_INDEX_OF_ZERO__MAX_NUM();
  PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_SWAP_RATE = config.PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_SWAP_RATE();
  PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_DEL_RATE = config.PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_DEL_RATE();
  PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_INS_RATE = config.PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_INS_RATE();
  PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_SUB_RATE = config.PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_SUB_RATE();

  PROB_COUNT_ODDS__MIN_VEC_LEN = config.PROB_COUNT_ODDS__MIN_VEC_LEN();
  PROB_COUNT_ODDS__MAX_VEC_LEN = config.PROB_COUNT_ODDS__MAX_VEC_LEN();
  PROB_COUNT_ODDS__MIN_NUM = config.PROB_COUNT_ODDS__MIN_NUM();
  PROB_COUNT_ODDS__MAX_NUM = config.PROB_COUNT_ODDS__MAX_NUM();
  PROB_COUNT_ODDS__MUTATION__PER_NUM_SWAP_RATE = config.PROB_COUNT_ODDS__MUTATION__PER_NUM_SWAP_RATE();
  PROB_COUNT_ODDS__MUTATION__PER_NUM_DEL_RATE = config.PROB_COUNT_ODDS__MUTATION__PER_NUM_DEL_RATE();
  PROB_COUNT_ODDS__MUTATION__PER_NUM_INS_RATE = config.PROB_COUNT_ODDS__MUTATION__PER_NUM_INS_RATE();
  PROB_COUNT_ODDS__MUTATION__PER_NUM_SUB_RATE = config.PROB_COUNT_ODDS__MUTATION__PER_NUM_SUB_RATE();

  PROB_MIRROR_IMAGE__MIN_VEC_LEN = config.PROB_MIRROR_IMAGE__MIN_VEC_LEN();
  PROB_MIRROR_IMAGE__MAX_VEC_LEN = config.PROB_MIRROR_IMAGE__MAX_VEC_LEN();
  PROB_MIRROR_IMAGE__MIN_NUM = config.PROB_MIRROR_IMAGE__MIN_NUM();
  PROB_MIRROR_IMAGE__MAX_NUM = config.PROB_MIRROR_IMAGE__MAX_NUM();
  PROB_MIRROR_IMAGE__MUTATION__PER_VEC_RANDOMIZE_VAL_RATE = config.PROB_MIRROR_IMAGE__MUTATION__PER_VEC_RANDOMIZE_VAL_RATE();
  PROB_MIRROR_IMAGE__MUTATION__PER_VEC_MIRROR_RATE = config.PROB_MIRROR_IMAGE__MUTATION__PER_VEC_MIRROR_RATE();
  PROB_MIRROR_IMAGE__MUTATION__COPY_RATE = config.PROB_MIRROR_IMAGE__MUTATION__COPY_RATE();
  PROB_MIRROR_IMAGE__MUTATION__INS_RATE = config.PROB_MIRROR_IMAGE__MUTATION__INS_RATE();
  PROB_MIRROR_IMAGE__MUTATION__DEL_RATE = config.PROB_MIRROR_IMAGE__MUTATION__DEL_RATE();
  PROB_MIRROR_IMAGE__MUTATION__PER_VEC_SHUFFLE_RATE = config.PROB_MIRROR_IMAGE__MUTATION__PER_VEC_SHUFFLE_RATE();

  PROB_VECTORS_SUMMED__MIN_VEC_LEN = config.PROB_VECTORS_SUMMED__MIN_VEC_LEN();
  PROB_VECTORS_SUMMED__MAX_VEC_LEN = config.PROB_VECTORS_SUMMED__MAX_VEC_LEN();
  PROB_VECTORS_SUMMED__MIN_NUM = config.PROB_VECTORS_SUMMED__MIN_NUM();
  PROB_VECTORS_SUMMED__MAX_NUM = config.PROB_VECTORS_SUMMED__MAX_NUM();
  PROB_VECTORS_SUMMED__MUTATION__PER_NUM_SUB_RATE = config.PROB_VECTORS_SUMMED__MUTATION__PER_NUM_SUB_RATE();
  PROB_VECTORS_SUMMED__MUTATION__COPY_RATE = config.PROB_VECTORS_SUMMED__MUTATION__COPY_RATE();
  PROB_VECTORS_SUMMED__MUTATION__INS_RATE = config.PROB_VECTORS_SUMMED__MUTATION__INS_RATE();
  PROB_VECTORS_SUMMED__MUTATION__DEL_RATE = config.PROB_VECTORS_SUMMED__MUTATION__DEL_RATE();

  PROB_SUM_OF_SQUARES__MIN_NUM = config.PROB_SUM_OF_SQUARES__MIN_NUM();
  PROB_SUM_OF_SQUARES__MAX_NUM = config.PROB_SUM_OF_SQUARES__MAX_NUM();
  PROB_SUM_OF_SQUARES__MUTATION__NUM_MUT_RATE = config.PROB_SUM_OF_SQUARES__MUTATION__NUM_MUT_RATE();

  PROB_VECTOR_AVERAGE__EPSILON = config.PROB_VECTOR_AVERAGE__EPSILON();
  PROB_VECTOR_AVERAGE__MIN_VEC_LEN = config.PROB_VECTOR_AVERAGE__MIN_VEC_LEN();
  PROB_VECTOR_AVERAGE__MAX_VEC_LEN = config.PROB_VECTOR_AVERAGE__MAX_VEC_LEN();
  PROB_VECTOR_AVERAGE__MIN_NUM = config.PROB_VECTOR_AVERAGE__MIN_NUM();
  PROB_VECTOR_AVERAGE__MAX_NUM = config.PROB_VECTOR_AVERAGE__MAX_NUM();
  PROB_VECTOR_AVERAGE__MUTATION__INS_RATE = config.PROB_VECTOR_AVERAGE__MUTATION__INS_RATE();
  PROB_VECTOR_AVERAGE__MUTATION__DEL_RATE = config.PROB_VECTOR_AVERAGE__MUTATION__DEL_RATE();
  PROB_VECTOR_AVERAGE__MUTATION__SUB_RATE = config.PROB_VECTOR_AVERAGE__MUTATION__SUB_RATE();

  PROB_MEDIAN__MIN_NUM = config.PROB_MEDIAN__MIN_NUM();
  PROB_MEDIAN__MAX_NUM = config.PROB_MEDIAN__MAX_NUM();
  PROB_MEDIAN__MUTATION__PER_NUM_COPY_RATE = config.PROB_MEDIAN__MUTATION__PER_NUM_COPY_RATE();
  PROB_MEDIAN__MUTATION__PER_NUM_SUB_RATE = config.PROB_MEDIAN__MUTATION__PER_NUM_SUB_RATE();
  PROB_MEDIAN__MUTATION__PER_NUM_SWAP_RATE = config.PROB_MEDIAN__MUTATION__PER_NUM_SWAP_RATE();
  
  PROB_SMALLEST__MIN_NUM = config.PROB_SMALLEST__MIN_NUM();
  PROB_SMALLEST__MAX_NUM = config.PROB_SMALLEST__MAX_NUM();
  PROB_SMALLEST__MUTATION__PER_NUM_SUB_RATE = config.PROB_SMALLEST__MUTATION__PER_NUM_SUB_RATE();
  PROB_SMALLEST__MUTATION__PER_NUM_SWAP_RATE = config.PROB_SMALLEST__MUTATION__PER_NUM_SWAP_RATE();

  PROB_GRADE__MIN_NUM = config.PROB_GRADE__MIN_NUM();
  PROB_GRADE__MAX_NUM = config.PROB_GRADE__MAX_NUM();
  PROB_GRADE__MUTATION__PER_NUM_RANDOMIZE_RATE = config.PROB_GRADE__MUTATION__PER_NUM_RANDOMIZE_RATE();
  PROB_GRADE__MUTATION__PER_NUM_ADJUST_RATE = config.PROB_GRADE__MUTATION__PER_NUM_ADJUST_RATE();

  DATA_DIRECTORY = config.DATA_DIRECTORY();
  SUMMARY_STATS_INTERVAL = config.SUMMARY_STATS_INTERVAL();
  SNAPSHOT_INTERVAL = config.SNAPSHOT_INTERVAL();
  SOLUTION_SCREEN_INTERVAL = config.SOLUTION_SCREEN_INTERVAL();

}

/// Setup the appropriate problem by calling its problem-specific setup function.
/// Which problem is setup will depend on configuration.
void ProgramSynthesisExperiment::SetupProblem() {
  emp_assert(emp::Has(problems, PROBLEM), "Unknown problem!", PROBLEM);
  // Big ol' switch statement to select appropriate problem to setup.
  switch (problems.at(PROBLEM).id) {
    case PROBLEM_ID::NumberIO: { SetupProblem_NumberIO(); break; }
    case PROBLEM_ID::SmallOrLarge: { SetupProblem_SmallOrLarge(); break; }
    case PROBLEM_ID::ForLoopIndex: { SetupProblem_ForLoopIndex(); break; }
    case PROBLEM_ID::CompareStringLengths: { SetupProblem_CompareStringLengths(); break; }
    case PROBLEM_ID::DoubleLetters: { SetupProblem_DoubleLetters(); break; }
    case PROBLEM_ID::CollatzNumbers: { SetupProblem_CollatzNumbers(); break; }
    case PROBLEM_ID::ReplaceSpaceWithNewline: { SetupProblem_ReplaceSpaceWithNewline(); break; }
    case PROBLEM_ID::StringDifferences: { SetupProblem_StringDifferences(); break; }
    case PROBLEM_ID::EvenSquares: { SetupProblem_EvenSquares(); break; }
    case PROBLEM_ID::WallisPi: { SetupProblem_WallisPi(); break; }
    case PROBLEM_ID::StringLengthsBackwards: { SetupProblem_StringLengthsBackwards(); break; }
    case PROBLEM_ID::LastIndexOfZero: { SetupProblem_LastIndexOfZero(); break; }
    case PROBLEM_ID::VectorAverage: { SetupProblem_VectorAverage(); break; }
    case PROBLEM_ID::CountOdds: { SetupProblem_CountOdds(); break; }
    case PROBLEM_ID::MirrorImage: { SetupProblem_MirrorImage(); break; }
    case PROBLEM_ID::SuperAnagrams: { SetupProblem_SuperAnagrams(); break; }
    case PROBLEM_ID::SumOfSquares: { SetupProblem_SumOfSquares(); break; }
    case PROBLEM_ID::VectorsSummed: { SetupProblem_VectorsSummed(); break; }
    case PROBLEM_ID::XWordLines: { SetupProblem_XWordLines(); break; }
    case PROBLEM_ID::PigLatin: { SetupProblem_PigLatin(); break; }
    case PROBLEM_ID::NegativeToZero: { SetupProblem_NegativeToZero(); break; }
    case PROBLEM_ID::ScrabbleScore: { SetupProblem_ScrabbleScore(); break; }
    case PROBLEM_ID::Checksum: { SetupProblem_Checksum(); break; }
    case PROBLEM_ID::Digits: { SetupProblem_Digits(); break; }
    case PROBLEM_ID::Grade: { SetupProblem_Grade(); break; }
    case PROBLEM_ID::Median: { SetupProblem_Median(); break; }
    case PROBLEM_ID::Smallest: { SetupProblem_Smallest(); break; }
    case PROBLEM_ID::Syllables: { SetupProblem_Syllables(); break; }
    default: {
      std::cout << "Unknown problem (" << PROBLEM << "). Exiting." << std::endl;
      exit(-1);
    }
  }
}

/// Setup hardware used to evaluate evolving programs.
void ProgramSynthesisExperiment::SetupHardware() {
  // Create new instruction library.
  inst_lib = emp::NewPtr<inst_lib_t>();
  // Create evaluation hardware.
  eval_hardware = emp::NewPtr<hardware_t>(inst_lib, random);
  // Configure the CPU.
  eval_hardware->SetMemSize(MEM_SIZE);                       // Configure size of memory.
  eval_hardware->SetMinTagSpecificity(MIN_TAG_SPECIFICITY);  // Configure minimum tag specificity required for tag-based referencing.
  eval_hardware->SetMaxCallDepth(MAX_CALL_DEPTH);            // Configure maximum depth of call stack (recursion limit).
  eval_hardware->SetMemTags(GenHadamardMatrix<TAG_WIDTH>()); // Configure memory location tags. Use Hadamard matrix for given TAG_WIDTH.

  // Configure call tag (tag used to call initial module during test evaluation).
  call_tag.Clear(); // Set initial call tag to all 0s.

  // What do at beginning of program evaluation (about to be run on potentially many tests)?
  // - Reset virtual hardware (reset hardware, clear module definitions, clear program).
  // - Set program to given organism's 'genome'.
  begin_program_eval.AddAction([this](prog_org_t & prog_org) {
    eval_hardware->Reset();
    eval_hardware->SetProgram(prog_org.GetGenome());
  });

  // What do at end of program evaluation (after being run on some number of tests)?
  // - Currently, nothing.
  
  // What to do before running program on a single test?
  // - Reset virtual hardware (reset global memory, reset call stack).
  begin_program_test.AddAction([this](prog_org_t & prog_org, emp::Ptr<TestOrg_Base> test_org_ptr) {
    eval_hardware->ResetHardware();
    eval_hardware->CallModule(call_tag, MIN_TAG_SPECIFICITY, true, false);
  });

  // Specify how we 'do' a program test.
  // - For specified evaluation time, advance evaluation hardware. If at any point
  //   the program's call stack is empty, automatically finish the evaluation.
  do_program_test.AddAction([this](prog_org_t & prog_org, emp::Ptr<TestOrg_Base> test_org_ptr) {
    // std::cout << "--- DO PROGRAM TEST ---" << std::endl;
    // std::cout << "==== Initial hardware state ====" << std::endl;
    // eval_hardware->PrintHardwareState();
    for (eval_time = 0; eval_time < PROG_EVAL_TIME; ++eval_time) {
      // std::cout << "==== Time = " << eval_time << "==== " << std::endl;
      do_program_advance.Trigger(prog_org);
      // eval_hardware->PrintHardwareState();
      if (eval_hardware->GetCallStackSize() == 0) break; // If call stack is ever completely empty, program is done early.
    }
    // exit(-1);
  });
  
  // How do we advance the evaluation hardware?
  do_program_advance.AddAction([this](prog_org_t &) {
    eval_hardware->SingleProcess();
  });
}

// Setup evaluation.
void ProgramSynthesisExperiment::SetupEvaluation() {
  switch (EVALUATION_MODE) {
    // In cohort evaluation, programs and tests are evaluated in 'cohorts'. Populations
    // are divided into cohorts (the number of cohorts for tests and programs must
    // be equal -- size of cohorts may be different).
    case (size_t)EVALUATION_TYPE::COHORT: {
      emp_assert(PROG_POP_SIZE % PROG_COHORT_SIZE == 0, "Program population size must be evenly divisible by program cohort size.");
      emp_assert(TEST_POP_SIZE % TEST_COHORT_SIZE == 0, "Test population size must be evenly divisible by test cohort size.");
      std::cout << "Setting up cohorts." << std::endl;
      test_cohorts.Setup(TEST_POP_SIZE, TEST_COHORT_SIZE);
      prog_cohorts.Setup(PROG_POP_SIZE, PROG_COHORT_SIZE);
      std::cout << "  # test cohorts = " << test_cohorts.GetCohortCnt() << std::endl;
      std::cout << "  # program cohorts = " << prog_cohorts.GetCohortCnt() << std::endl;
      if (test_cohorts.GetCohortCnt() != prog_cohorts.GetCohortCnt()) {
        std::cout << "ERROR: Test cohort count must the same as program cohort count in COHORT mode. Exiting." << std::endl;
        exit(-1);
      }
      NUM_COHORTS = prog_cohorts.GetCohortCnt();
      PROGRAM_MAX_PASSES = TEST_COHORT_SIZE;

      // Setup program world on placement response.
      prog_world->OnPlacement([this](size_t pos) {
        // On placement, reset organism phenotype.
        prog_world->GetOrg(pos).GetPhenotype().Reset(TEST_COHORT_SIZE);
      });
      
      // Setup test case world on placement response.
      OnPlacement_ActiveTestCaseWorld([this](size_t pos) { 
        // On placement, reset organism phenotype.
        GetTestPhenotype(pos).Reset(PROG_COHORT_SIZE); 
      });

      // What should happen on evaluation?
      do_evaluation_sig.AddAction([this]() {
        // Randomize the cohorts.
        prog_cohorts.Randomize(*random);
        test_cohorts.Randomize(*random);
        // For each cohort, evaluate all programs against all tests in corresponding cohort.
        for (size_t cID = 0; cID < prog_cohorts.GetCohortCnt(); ++cID) {
          for (size_t pID = 0; pID < PROG_COHORT_SIZE; ++pID) {
            // Get program organism specified by pID.
            prog_org_t & prog_org = prog_world->GetOrg(prog_cohorts.GetWorldID(cID, pID));
            begin_program_eval.Trigger(prog_org);
            for (size_t tID = 0; tID < TEST_COHORT_SIZE; ++tID) {
              const size_t test_world_id = test_cohorts.GetWorldID(cID, tID);
              
              // Evaluate program on this test.
              TestResult result = EvaluateWorldTest(prog_org, test_world_id);
              
              // Grab references to relevant phenotypes.
              test_org_phen_t & test_phen = GetTestPhenotype(test_world_id);
              prog_org_phen_t & prog_phen = prog_org.GetPhenotype();
              
              // Test result contents:
              // - double score;
              // - bool pass;
              // - bool sub;
        
              // Update program phenotype.
              prog_phen.RecordScore(tID, result.score);
              prog_phen.RecordPass(tID, result.pass);
              prog_phen.RecordSubmission(result.sub);
              
              // Update test phenotype
              test_phen.RecordScore(pID, result.score);
              test_phen.RecordPass(pID, result.pass);
            }
            end_program_eval.Trigger(prog_org);
          }
        }
      });
      break;
    }
    // In full evaluation mode, evaluate all programs on all tests.
    case (size_t)EVALUATION_TYPE::FULL: {
      std::cout << "Setting up full evaluation. No cohorts here." << std::endl;
      
      // Setup program world on placement signal response.
      prog_world->OnPlacement([this](size_t pos) {
        // Reset program phenotype on placement.
        prog_world->GetOrg(pos).GetPhenotype().Reset(TEST_POP_SIZE); 
      });

      // Setup test world on placement signal response.
      OnPlacement_ActiveTestCaseWorld([this](size_t pos) {
        // Reset test phenotype on placement.
        GetTestPhenotype(pos).Reset(PROG_POP_SIZE);   
      });
      
      PROGRAM_MAX_PASSES = TEST_POP_SIZE;
      NUM_COHORTS = 0;

      // What should happen on evaluation?
      do_evaluation_sig.AddAction([this]() {
        for (size_t pID = 0; pID < PROG_POP_SIZE; ++pID) {
          emp_assert(prog_world->IsOccupied(pID));
          prog_org_t & prog_org = prog_world->GetOrg(pID);
          begin_program_eval.Trigger(prog_org);
          for (size_t tID = 0; tID < TEST_POP_SIZE; ++tID) {

            TestResult result = EvaluateWorldTest(prog_org, tID);
            
            // Grab references to test and program phenotypes.
            test_org_phen_t & test_phen = GetTestPhenotype(tID);
            prog_org_phen_t & prog_phen = prog_org.GetPhenotype();

            // Test result contents:
            // - double score;
            // - bool pass;
            // - bool sub;

            // Update program phenotype.
            prog_phen.RecordScore(tID, result.score);
            prog_phen.RecordPass(tID, result.pass);
            prog_phen.RecordSubmission(result.sub);
            
            // Update phenotypes.
            test_phen.RecordScore(pID, result.score);
            test_phen.RecordPass(pID, result.pass);
          }
          end_program_eval.Trigger(prog_org);
        }
      });

      break;
    }
    case (size_t)EVALUATION_TYPE::PROG_ONLY_COHORT: {
      std::cout << "Setting up full evaluation with PROGRAM ONLY cohorts." << std::endl;

      emp_assert(PROG_POP_SIZE % PROG_COHORT_SIZE == 0, "Program population size must be evenly divisible by program cohort size.");
      std::cout << "Setting up PROGRAM cohorts." << std::endl;
      prog_cohorts.Setup(PROG_POP_SIZE, PROG_COHORT_SIZE);
      std::cout << "  # program cohorts = " << prog_cohorts.GetCohortCnt() << std::endl;
      NUM_COHORTS = prog_cohorts.GetCohortCnt();
      PROGRAM_MAX_PASSES = TEST_POP_SIZE;  
      
      // Setup program world on placement signal response.
      prog_world->OnPlacement([this](size_t pos) {
        // Reset program phenotype on placement.
        prog_world->GetOrg(pos).GetPhenotype().Reset(TEST_POP_SIZE); 
      });

      // Setup test world on placement signal response.
      OnPlacement_ActiveTestCaseWorld([this](size_t pos) {
        // Reset test phenotype on placement.
        GetTestPhenotype(pos).Reset(PROG_POP_SIZE);   
      });
      

      // What should happen on evaluation?
      do_evaluation_sig.AddAction([this]() {

        prog_cohorts.Randomize(*random);

        for (size_t pID = 0; pID < PROG_POP_SIZE; ++pID) {
          emp_assert(prog_world->IsOccupied(pID));
          prog_org_t & prog_org = prog_world->GetOrg(pID);
          begin_program_eval.Trigger(prog_org);
          for (size_t tID = 0; tID < TEST_POP_SIZE; ++tID) {

            TestResult result = EvaluateWorldTest(prog_org, tID);
            
            // Grab references to test and program phenotypes.
            test_org_phen_t & test_phen = GetTestPhenotype(tID);
            prog_org_phen_t & prog_phen = prog_org.GetPhenotype();

            // Test result contents:
            // - double score;
            // - bool pass;
            // - bool sub;

            // Update program phenotype.
            prog_phen.RecordScore(tID, result.score);
            prog_phen.RecordPass(tID, result.pass);
            prog_phen.RecordSubmission(result.sub);
            
            // Update phenotypes.
            test_phen.RecordScore(pID, result.score);
            test_phen.RecordPass(pID, result.pass);
          }
          end_program_eval.Trigger(prog_org);
        }
      });
      break;
    }
    case (size_t)EVALUATION_TYPE::TEST_DOWNSAMPLING: {
      std::cout << "Setting up downsampled evaluation with NO program cohorts." << std::endl;
      test_cohorts.Setup(TEST_POP_SIZE, TEST_COHORT_SIZE);
      std::cout << "  # test cohorts = " << test_cohorts.GetCohortCnt() << std::endl;
      NUM_COHORTS = prog_cohorts.GetCohortCnt();
      PROGRAM_MAX_PASSES = TEST_COHORT_SIZE;  
      
      // Setup program world on placement signal response.
      prog_world->OnPlacement([this](size_t pos) {
        // Reset program phenotype on placement.
        prog_world->GetOrg(pos).GetPhenotype().Reset(TEST_COHORT_SIZE); 
      });

      // Setup test world on placement signal response.
      OnPlacement_ActiveTestCaseWorld([this](size_t pos) {
        // Reset test phenotype on placement.
        GetTestPhenotype(pos).Reset(PROG_POP_SIZE);   
      });
      

      // What should happen on evaluation?
      do_evaluation_sig.AddAction([this]() {

        test_cohorts.Randomize(*random);

        for (size_t pID = 0; pID < PROG_POP_SIZE; ++pID) {
          emp_assert(prog_world->IsOccupied(pID));
          prog_org_t & prog_org = prog_world->GetOrg(pID);
          begin_program_eval.Trigger(prog_org);
          for (size_t tID = 0; tID < TEST_COHORT_SIZE; ++tID) {
            const size_t test_world_id = test_cohorts.GetWorldID(0, tID);
            TestResult result = EvaluateWorldTest(prog_org, test_world_id);
            
            // Grab references to test and program phenotypes.
            test_org_phen_t & test_phen = GetTestPhenotype(test_world_id);
            prog_org_phen_t & prog_phen = prog_org.GetPhenotype();

            // Test result contents:
            // - double score;
            // - bool pass;
            // - bool sub;

            // Update program phenotype.
            prog_phen.RecordScore(tID, result.score);
            prog_phen.RecordPass(tID, result.pass);
            prog_phen.RecordSubmission(result.sub);
            
            // Update phenotypes.
            test_phen.RecordScore(pID, result.score);
            test_phen.RecordPass(pID, result.pass);
          }
          end_program_eval.Trigger(prog_org);
        }
      });
      break;
    }
    default: {
      std::cout << "Unknown EVALUATION_MODE (" << EVALUATION_MODE << "). Exiting." << std::endl;
      exit(-1);
    } 
  }

  // Add generic evaluation action - calculate num_passes/fails for tests and programs.
  do_evaluation_sig.AddAction([this]() {
    // Sum pass totals for programs, find 'dominant' program.
    double cur_best_score = 0;
    for (size_t pID = 0; pID < prog_world->GetSize(); ++pID) {
      emp_assert(prog_world->IsOccupied(pID));

      prog_org_t & prog_org = prog_world->GetOrg(pID);
      const size_t pass_total = prog_org.GetPhenotype().num_passes;
      const double total_score = prog_org.GetPhenotype().total_score;

      // Is this the highest pass total program this generation?
      if (total_score > cur_best_score || pID == 0) {
        dominant_prog_id = pID;
        cur_best_score = pass_total;
      }
      // At this point, program has been evaluated against all tests. .
      if (pass_total == PROGRAM_MAX_PASSES && prog_org.GetGenome().GetSize() < smallest_prog_sol_size) {
        stats_util.cur_progID = pID;
        if (ScreenForSolution(prog_org)) {
          if (!solution_found) { update_first_solution_found = prog_world->GetUpdate(); }
          solution_found = true;
          smallest_prog_sol_size = prog_org.GetGenome().GetSize();
          // Add to solutions file.
          solution_file->Update();
        }
      }
    }

    // Sum pass/fail totals for tests.
    for (size_t tID = 0; tID < TEST_POP_SIZE; ++tID) {
      test_org_phen_t & test_phen = GetTestPhenotype(tID);
      if (test_phen.num_fails > cur_best_score || tID == 0) { // NOTE - will need to be updated if switch to gradient fitness functions (will be problem specific whether or not use a gradient).
        dominant_test_id = tID;
        cur_best_score = test_phen.num_fails;
      }
    }
  });
}

/// Setup selection for programs and tests.
void ProgramSynthesisExperiment::SetupSelection() {
  // (1) Setup program selection.
  SetupProgramSelection();
  // (2) Setup test selection.
  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::COEVOLUTION || TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::STATIC_COEVO) {
    std::cout << "COEVOLUTION training example mode detected, setting up test case selection." << std::endl;
    if (prob_NumberIO_world != nullptr) { SetupTestSelection(prob_NumberIO_world, prob_utils_NumberIO.lexicase_fit_set); }
    else if (prob_SmallOrLarge_world != nullptr) { SetupTestSelection(prob_SmallOrLarge_world, prob_utils_SmallOrLarge.lexicase_fit_set); }
    else if (prob_ForLoopIndex_world != nullptr) { SetupTestSelection(prob_ForLoopIndex_world, prob_utils_ForLoopIndex.lexicase_fit_set); }
    else if (prob_CompareStringLengths_world != nullptr) { SetupTestSelection(prob_CompareStringLengths_world, prob_utils_CompareStringLengths.lexicase_fit_set); }
    else if (prob_DoubleLetters_world != nullptr) { SetupTestSelection(prob_DoubleLetters_world, prob_utils_DoubleLetters.lexicase_fit_set); }
    else if (prob_CollatzNumbers_world != nullptr) { SetupTestSelection(prob_CollatzNumbers_world, prob_utils_CollatzNumbers.lexicase_fit_set); }
    else if (prob_ReplaceSpaceWithNewline_world != nullptr) { SetupTestSelection(prob_ReplaceSpaceWithNewline_world, prob_utils_ReplaceSpaceWithNewline.lexicase_fit_set); }
    else if (prob_StringDifferences_world != nullptr) { SetupTestSelection(prob_StringDifferences_world, prob_utils_StringDifferences.lexicase_fit_set); }
    else if (prob_EvenSquares_world != nullptr) { SetupTestSelection(prob_EvenSquares_world, prob_utils_EvenSquares.lexicase_fit_set); }
    else if (prob_WallisPi_world != nullptr) { SetupTestSelection(prob_WallisPi_world, prob_utils_WallisPi.lexicase_fit_set); }
    else if (prob_StringLengthsBackwards_world != nullptr) { SetupTestSelection(prob_StringLengthsBackwards_world, prob_utils_StringLengthsBackwards.lexicase_fit_set); }
    else if (prob_LastIndexOfZero_world != nullptr) { SetupTestSelection(prob_LastIndexOfZero_world, prob_utils_LastIndexOfZero.lexicase_fit_set); }
    else if (prob_VectorAverage_world != nullptr) { SetupTestSelection(prob_VectorAverage_world, prob_utils_VectorAverage.lexicase_fit_set); }
    else if (prob_CountOdds_world != nullptr) { SetupTestSelection(prob_CountOdds_world, prob_utils_CountOdds.lexicase_fit_set); }
    else if (prob_MirrorImage_world != nullptr) { SetupTestSelection(prob_MirrorImage_world, prob_utils_MirrorImage.lexicase_fit_set); }
    else if (prob_SuperAnagrams_world != nullptr) { SetupTestSelection(prob_SuperAnagrams_world, prob_utils_SuperAnagrams.lexicase_fit_set); }
    else if (prob_SumOfSquares_world != nullptr) { SetupTestSelection(prob_SumOfSquares_world, prob_utils_SumOfSquares.lexicase_fit_set); }
    else if (prob_VectorsSummed_world != nullptr) { SetupTestSelection(prob_VectorsSummed_world, prob_utils_VectorsSummed.lexicase_fit_set); }
    else if (prob_XWordLines_world != nullptr) { SetupTestSelection(prob_XWordLines_world, prob_utils_XWordLines.lexicase_fit_set); }
    else if (prob_PigLatin_world != nullptr) { SetupTestSelection(prob_PigLatin_world, prob_utils_PigLatin.lexicase_fit_set); }
    else if (prob_NegativeToZero_world != nullptr) { SetupTestSelection(prob_NegativeToZero_world, prob_utils_NegativeToZero.lexicase_fit_set); }
    else if (prob_ScrabbleScore_world != nullptr) { SetupTestSelection(prob_ScrabbleScore_world, prob_utils_ScrabbleScore.lexicase_fit_set); }
    else if (prob_Checksum_world != nullptr) { SetupTestSelection(prob_Checksum_world, prob_utils_Checksum.lexicase_fit_set); }
    else if (prob_Digits_world != nullptr) { SetupTestSelection(prob_Digits_world, prob_utils_Digits.lexicase_fit_set); }
    else if (prob_Grade_world != nullptr) { SetupTestSelection(prob_Grade_world, prob_utils_Grade.lexicase_fit_set); }
    else if (prob_Median_world != nullptr) { SetupTestSelection(prob_Median_world, prob_utils_Median.lexicase_fit_set); }
    else if (prob_Smallest_world != nullptr) { SetupTestSelection(prob_Smallest_world, prob_utils_Smallest.lexicase_fit_set); }
    else if (prob_Syllables_world != nullptr) { SetupTestSelection(prob_Syllables_world, prob_utils_Syllables.lexicase_fit_set); }
    else { std::cout << "AHH! More than one test case world has been created. Exiting." << std::endl; exit(-1); }
  }
}

/// Setup program selection.
void ProgramSynthesisExperiment::SetupProgramSelection() {
  switch (PROG_SELECTION_MODE) {
    case (size_t)SELECTION_TYPE::LEXICASE: {
      std::cout << "Setting up program LEXICASE selection." << std::endl;
      // Lexicase selection requires full evaluation mode?
      emp_assert(EVALUATION_MODE == (size_t)EVALUATION_TYPE::FULL, "Lexicase selection requires FULL evaluation mode.");
      // Setup program fitness functions.
      // - 1 function for every test score.
      for (size_t i = 0; i < TEST_POP_SIZE; ++i) {
        lexicase_prog_fit_set.push_back([i](prog_org_t & prog_org) {
          emp_assert(i < prog_org.GetPhenotype().test_scores.size(), i, prog_org.GetPhenotype().test_scores.size());
          double score = prog_org.GetPhenotype().test_scores[i];
          return score;
        });
      }
      // Add pressure for small size
      lexicase_prog_fit_set.push_back([this](prog_org_t & prog_org) {
        if (prog_org.GetPhenotype().num_passes == PROGRAM_MAX_PASSES) {
          return (double)(MAX_PROG_SIZE - prog_org.GetGenome().GetSize());
        } 
        return 0.0;
      });
      // Add selection action
      do_selection_sig.AddAction([this]() {
        emp::LexicaseSelect_NAIVE(*prog_world,
                                  lexicase_prog_fit_set,
                                  PROG_POP_SIZE,
                                  PROG_LEXICASE_MAX_FUNS); // TODO - track lexicase fit fun stats
      });
      break;
    }
    case (size_t)SELECTION_TYPE::COHORT_LEXICASE: {
      std::cout << "Setting up program COHORT LEXICASE selection." << std::endl;
      emp_assert(EVALUATION_MODE == (size_t)EVALUATION_TYPE::COHORT, "Cohort lexicase selection requires COHORT evaluation mode.");
      emp_assert(PROG_COHORT_SIZE * prog_cohorts.GetCohortCnt() == PROG_POP_SIZE);

      // Setup program fitness functions.
      // - 1 function for every test cohort member.
      for (size_t i = 0; i < TEST_COHORT_SIZE; ++i) {
        lexicase_prog_fit_set.push_back([i](prog_org_t & prog_org) {
          emp_assert(i < prog_org.GetPhenotype().test_scores.size(), i, prog_org.GetPhenotype().test_scores.size());
          double score = prog_org.GetPhenotype().test_scores[i];
          return score;
        });
      }
      // Selection pressure for small program size.
      lexicase_prog_fit_set.push_back([this](prog_org_t & prog_org) {
        if (prog_org.GetPhenotype().num_passes == PROGRAM_MAX_PASSES) {
          return (double)(MAX_PROG_SIZE - prog_org.GetGenome().GetSize());
        } 
        return 0.0;
      });
      // Add selection action
      do_selection_sig.AddAction([this]() {
        for (size_t cID = 0; cID < prog_cohorts.GetCohortCnt(); ++cID) {
          emp::CohortLexicaseSelect_NAIVE(*prog_world,
                                          lexicase_prog_fit_set,
                                          prog_cohorts.GetCohort(cID),
                                          PROG_COHORT_SIZE,
                                          PROG_COHORTLEXICASE_MAX_FUNS);
        }
      });
      break;
    }
    case (size_t)SELECTION_TYPE::PROG_ONLY_COHORT_LEXICASE: {
      std::cout << "Setting up program PROGRAM ONLY COHORT LEXICASE selection." << std::endl;
      emp_assert(EVALUATION_MODE == (size_t)EVALUATION_TYPE::PROG_ONLY_COHORT, "Program only Cohort lexicase selection requires program only COHORT evaluation mode.");
      emp_assert(PROG_COHORT_SIZE * prog_cohorts.GetCohortCnt() == PROG_POP_SIZE);

      // Setup program fitness functions.
      // - 1 function for every test cohort member.
      for (size_t i = 0; i < TEST_POP_SIZE; ++i) {
        lexicase_prog_fit_set.push_back([i](prog_org_t & prog_org) {
          emp_assert(i < prog_org.GetPhenotype().test_scores.size(), i, prog_org.GetPhenotype().test_scores.size());
          double score = prog_org.GetPhenotype().test_scores[i];
          return score;
        });
      }
      // Selection pressure for small program size.
      lexicase_prog_fit_set.push_back([this](prog_org_t & prog_org) {
        if (prog_org.GetPhenotype().num_passes == PROGRAM_MAX_PASSES) {
          return (double)(MAX_PROG_SIZE - prog_org.GetGenome().GetSize());
        } 
        return 0.0;
      });
      // Add selection action
      do_selection_sig.AddAction([this]() {
        for (size_t cID = 0; cID < prog_cohorts.GetCohortCnt(); ++cID) {
          emp::CohortLexicaseSelect_NAIVE(*prog_world,
                                          lexicase_prog_fit_set,
                                          prog_cohorts.GetCohort(cID),
                                          PROG_COHORT_SIZE,
                                          PROG_COHORTLEXICASE_MAX_FUNS);
        }
      });

      break;
    }
    case (size_t)SELECTION_TYPE::TEST_DOWNSAMPLING_LEXICASE: {
      std::cout << "Setting up program selection - TEST_DOWNSAMPLING_LEXICASE experiment." << std::endl;
      emp_assert(EVALUATION_MODE == (size_t)EVALUATION_TYPE::TEST_DOWNSAMPLING, "Program only Cohort lexicase selection requires program only COHORT evaluation mode.");

      // Setup program fitness functions.
      // - 1 function for every test cohort member.
      for (size_t i = 0; i < TEST_COHORT_SIZE; ++i) {
        lexicase_prog_fit_set.push_back([i](prog_org_t & prog_org) {
          emp_assert(i < prog_org.GetPhenotype().test_scores.size(), i, prog_org.GetPhenotype().test_scores.size());
          double score = prog_org.GetPhenotype().test_scores[i];
          return score;
        });
      }
      // Selection pressure for small program size.
      lexicase_prog_fit_set.push_back([this](prog_org_t & prog_org) {
        if (prog_org.GetPhenotype().num_passes == PROGRAM_MAX_PASSES) {
          return (double)(MAX_PROG_SIZE - prog_org.GetGenome().GetSize());
        } 
        return 0.0;
      });
      // Add selection action
      do_selection_sig.AddAction([this]() {
          emp::LexicaseSelect_NAIVE(*prog_world,
                                     lexicase_prog_fit_set,
                                     PROG_POP_SIZE,
                                     PROG_LEXICASE_MAX_FUNS); // TODO - track lexicase fit fun stats                          
      });
      break;
    }
    case (size_t)SELECTION_TYPE::TOURNAMENT: {
      std::cout << "Setting up program TOURNAMENT selection." << std::endl;
      do_selection_sig.AddAction([this]() {
        emp::TournamentSelect(*prog_world, PROG_TOURNAMENT_SIZE, PROG_POP_SIZE);
      });
      break;
    }
    case (size_t)SELECTION_TYPE::DRIFT: {
      std::cout << "Setting up program DRIFT selection." << std::endl;
      do_selection_sig.AddAction([this]() {
        emp::RandomSelect(*prog_world, PROG_POP_SIZE);
      });
      break;
    }
    default: {
      std::cout << "Unknown PROG_SELECTION_MODE (" << PROG_SELECTION_MODE << "). Exiting." << std::endl;
      exit(-1);
    }
  }
}

/// Setup program/test mutation.
void ProgramSynthesisExperiment::SetupMutation() {
  // (1) Setup program mutations
  SetupProgramMutation();
  // (2) Setup test mutations
  SetupTestMutation();
  // (3) Setup world(s) to auto mutate.
  end_setup_sig.AddAction([this]() {
    prog_world->SetAutoMutate();      // After we've initialized populations, turn auto mutate on.
  });

  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::COEVOLUTION) {
    std::cout << "COEVOLUTION training mode detected. Setting test world to AUTO-MUTATE." << std::endl;
    if (prob_NumberIO_world != nullptr) { end_setup_sig.AddAction([this]() { prob_NumberIO_world->SetAutoMutate(); }); }
    else if (prob_SmallOrLarge_world != nullptr) { end_setup_sig.AddAction([this]() { prob_SmallOrLarge_world->SetAutoMutate(); }); }
    else if (prob_ForLoopIndex_world != nullptr) { end_setup_sig.AddAction([this]() { prob_ForLoopIndex_world->SetAutoMutate(); }); }
    else if (prob_CompareStringLengths_world != nullptr) { end_setup_sig.AddAction([this]() { prob_CompareStringLengths_world->SetAutoMutate(); }); }
    else if (prob_DoubleLetters_world != nullptr) { end_setup_sig.AddAction([this]() { prob_DoubleLetters_world->SetAutoMutate(); }); }
    else if (prob_CollatzNumbers_world != nullptr) { end_setup_sig.AddAction([this]() { prob_CollatzNumbers_world->SetAutoMutate(); }); }
    else if (prob_ReplaceSpaceWithNewline_world != nullptr) { end_setup_sig.AddAction([this]() { prob_ReplaceSpaceWithNewline_world->SetAutoMutate(); }); }
    else if (prob_StringDifferences_world != nullptr) { end_setup_sig.AddAction([this]() { prob_StringDifferences_world->SetAutoMutate(); }); }
    else if (prob_EvenSquares_world != nullptr) { end_setup_sig.AddAction([this]() { prob_EvenSquares_world->SetAutoMutate(); }); }
    else if (prob_WallisPi_world != nullptr) { end_setup_sig.AddAction([this]() { prob_WallisPi_world->SetAutoMutate(); }); }
    else if (prob_StringLengthsBackwards_world != nullptr) { end_setup_sig.AddAction([this]() { prob_StringLengthsBackwards_world->SetAutoMutate(); }); }
    else if (prob_LastIndexOfZero_world != nullptr) { end_setup_sig.AddAction([this]() { prob_LastIndexOfZero_world->SetAutoMutate(); }); }
    else if (prob_VectorAverage_world != nullptr) { end_setup_sig.AddAction([this]() { prob_VectorAverage_world->SetAutoMutate(); }); }
    else if (prob_CountOdds_world != nullptr) { end_setup_sig.AddAction([this]() { prob_CountOdds_world->SetAutoMutate(); }); }
    else if (prob_MirrorImage_world != nullptr) { end_setup_sig.AddAction([this]() { prob_MirrorImage_world->SetAutoMutate(); }); }
    else if (prob_SuperAnagrams_world != nullptr) { end_setup_sig.AddAction([this]() { prob_SuperAnagrams_world->SetAutoMutate(); }); }
    else if (prob_SumOfSquares_world != nullptr) { end_setup_sig.AddAction([this]() { prob_SumOfSquares_world->SetAutoMutate(); }); }
    else if (prob_VectorsSummed_world != nullptr) { end_setup_sig.AddAction([this]() { prob_VectorsSummed_world->SetAutoMutate(); }); }
    else if (prob_XWordLines_world != nullptr) { end_setup_sig.AddAction([this]() { prob_XWordLines_world->SetAutoMutate(); }); }
    else if (prob_PigLatin_world != nullptr) { end_setup_sig.AddAction([this]() { prob_PigLatin_world->SetAutoMutate(); }); }
    else if (prob_NegativeToZero_world != nullptr) { end_setup_sig.AddAction([this]() { prob_NegativeToZero_world->SetAutoMutate(); }); }
    else if (prob_ScrabbleScore_world != nullptr) { end_setup_sig.AddAction([this]() { prob_ScrabbleScore_world->SetAutoMutate(); }); }
    else if (prob_Checksum_world != nullptr) { end_setup_sig.AddAction([this]() { prob_Checksum_world->SetAutoMutate(); }); }
    else if (prob_Digits_world != nullptr) { end_setup_sig.AddAction([this]() { prob_Digits_world->SetAutoMutate(); }); }
    else if (prob_Grade_world != nullptr) { end_setup_sig.AddAction([this]() { prob_Grade_world->SetAutoMutate(); }); }
    else if (prob_Median_world != nullptr) { end_setup_sig.AddAction([this]() { prob_Median_world->SetAutoMutate(); }); }
    else if (prob_Smallest_world != nullptr) { end_setup_sig.AddAction([this]() { prob_Smallest_world->SetAutoMutate(); }); }
    else if (prob_Syllables_world != nullptr) { end_setup_sig.AddAction([this]() { prob_Syllables_world->SetAutoMutate(); }); }
    else { std::cout << "AHH! None of the worlds have been initialized. Exiting." << std::endl; exit(-1); }
  } else if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::STATIC_COEVO) {
    std::cout << "STATIC_COEVO training mode detected. Setting test world to AUTO-MUTATE (only if mut pos >= " << STATIC_COEVO__NUM_STATIC_TESTCASES << ")." << std::endl;
    // std::function<bool(size_t pos)> test_fun = [this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; };
    if (prob_NumberIO_world != nullptr) { end_setup_sig.AddAction([this]() { prob_NumberIO_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_SmallOrLarge_world != nullptr) { end_setup_sig.AddAction([this]() { prob_SmallOrLarge_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_ForLoopIndex_world != nullptr) { end_setup_sig.AddAction([this]() { prob_ForLoopIndex_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_CompareStringLengths_world != nullptr) { end_setup_sig.AddAction([this]() { prob_CompareStringLengths_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_DoubleLetters_world != nullptr) { end_setup_sig.AddAction([this]() { prob_DoubleLetters_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_CollatzNumbers_world != nullptr) { end_setup_sig.AddAction([this]() { prob_CollatzNumbers_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_ReplaceSpaceWithNewline_world != nullptr) { end_setup_sig.AddAction([this]() { prob_ReplaceSpaceWithNewline_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_StringDifferences_world != nullptr) { end_setup_sig.AddAction([this]() { prob_StringDifferences_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_EvenSquares_world != nullptr) { end_setup_sig.AddAction([this]() { prob_EvenSquares_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_WallisPi_world != nullptr) { end_setup_sig.AddAction([this]() { prob_WallisPi_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_StringLengthsBackwards_world != nullptr) { end_setup_sig.AddAction([this]() { prob_StringLengthsBackwards_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_LastIndexOfZero_world != nullptr) { end_setup_sig.AddAction([this]() { prob_LastIndexOfZero_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_VectorAverage_world != nullptr) { end_setup_sig.AddAction([this]() { prob_VectorAverage_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_CountOdds_world != nullptr) { end_setup_sig.AddAction([this]() { prob_CountOdds_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_MirrorImage_world != nullptr) { end_setup_sig.AddAction([this]() { prob_MirrorImage_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_SuperAnagrams_world != nullptr) { end_setup_sig.AddAction([this]() { prob_SuperAnagrams_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_SumOfSquares_world != nullptr) { end_setup_sig.AddAction([this]() { prob_SumOfSquares_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_VectorsSummed_world != nullptr) { end_setup_sig.AddAction([this]() { prob_VectorsSummed_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_XWordLines_world != nullptr) { end_setup_sig.AddAction([this]() { prob_XWordLines_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_PigLatin_world != nullptr) { end_setup_sig.AddAction([this]() { prob_PigLatin_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_NegativeToZero_world != nullptr) { end_setup_sig.AddAction([this]() { prob_NegativeToZero_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_ScrabbleScore_world != nullptr) { end_setup_sig.AddAction([this]() { prob_ScrabbleScore_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_Checksum_world != nullptr) { end_setup_sig.AddAction([this]() { prob_Checksum_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_Digits_world != nullptr) { end_setup_sig.AddAction([this]() { prob_Digits_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_Grade_world != nullptr) { end_setup_sig.AddAction([this]() { prob_Grade_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_Median_world != nullptr) { end_setup_sig.AddAction([this]() { prob_Median_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_Smallest_world != nullptr) { end_setup_sig.AddAction([this]() { prob_Smallest_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else if (prob_Syllables_world != nullptr) { end_setup_sig.AddAction([this]() { prob_Syllables_world->SetAutoMutate([this](size_t p) { return p >= STATIC_COEVO__NUM_STATIC_TESTCASES; }); }); }
    else { std::cout << "AHH! None of the worlds have been initialized. Exiting." << std::endl; exit(-1); }
  }
}

/// Setup program mutation
void ProgramSynthesisExperiment::SetupProgramMutation() {
  // Configure program mutator.
  prog_mutator.MAX_PROGRAM_LEN = MAX_PROG_SIZE;
  prog_mutator.MIN_PROGRAM_LEN = MIN_PROG_SIZE;
  prog_mutator.PER_BIT_FLIP = PROG_MUT__PER_BIT_FLIP;
  prog_mutator.PER_INST_SUB = PROG_MUT__PER_INST_SUB;
  prog_mutator.PER_INST_INS = PROG_MUT__PER_INST_INS;
  prog_mutator.PER_INST_DEL = PROG_MUT__PER_INST_DEL;
  prog_mutator.PER_PROG_SLIP = PROG_MUT__PER_PROG_SLIP;
  prog_mutator.PER_MOD_DUP = PROG_MUT__PER_MOD_DUP;
  prog_mutator.PER_MOD_DEL = PROG_MUT__PER_MOD_DEL;

  // Configure world mutation function.
  prog_world->SetMutFun([this](prog_org_t & prog_org, emp::Random & rnd) {
    return prog_mutator.Mutate(rnd, prog_org.GetGenome());
  });
}

/// Setup fitness functions for tests/programs.
void ProgramSynthesisExperiment::SetupFitFuns() {
  // (1) Setup program mutations
  SetupProgramFitFun();
  // (2) Setup test mutations
  SetupTestFitFun();
}

/// Setup program fitness function.
void ProgramSynthesisExperiment::SetupProgramFitFun() {
  prog_world->SetFitFun([this](prog_org_t & prog_org) {
    double fitness = prog_org.GetPhenotype().total_score;
    if (prog_org.GetPhenotype().num_passes == PROGRAM_MAX_PASSES) { // Add 'smallness' bonus.
      fitness += ((double)(MAX_PROG_SIZE - prog_org.GetGenome().GetSize()))/(double)MAX_PROG_SIZE;
    }
    return fitness;
  });
}

/// Setup data collection.
void ProgramSynthesisExperiment::SetupDataCollection() {
  std::cout << "Setting up data collection." << std::endl;
  // Make a data directory
  mkdir(DATA_DIRECTORY.c_str(), ACCESSPERMS);
  if (DATA_DIRECTORY.back() != '/') DATA_DIRECTORY += '/';

  // Setup snapshot program stats.
  SetupProgramStats();

  std::function<size_t(void)> get_update = [this]() { return prog_world->GetUpdate(); };
  std::function<double(void)> get_evaluations = [this]() {
    if (EVALUATION_MODE == (size_t)EVALUATION_TYPE::COHORT || EVALUATION_MODE == (size_t)EVALUATION_TYPE::TEST_DOWNSAMPLING) {
      // update * cohort size * cohort size * num cohorts
      return prog_world->GetUpdate() * PROG_COHORT_SIZE * TEST_COHORT_SIZE * NUM_COHORTS;
    } else {
      return prog_world->GetUpdate() * PROG_POP_SIZE * TEST_POP_SIZE;
    }
  };

  // Setup program systematics
  prog_genotypic_systematics = emp::NewPtr<emp::Systematics<prog_org_t, prog_org_gen_t>>([](const prog_org_t & o) { return o.GetGenome(); });
  prog_genotypic_systematics->AddEvolutionaryDistinctivenessDataNode();
  prog_genotypic_systematics->AddPairwiseDistanceDataNode();
  prog_genotypic_systematics->AddPhylogeneticDiversityDataNode();
  prog_world->AddSystematics(prog_genotypic_systematics, "prog_genotype");
  auto & prog_gen_sys_file = prog_world->SetupSystematicsFile("prog_genotype", DATA_DIRECTORY + "/prog_gen_sys.csv", false);
  prog_gen_sys_file.SetTimingRepeat(SUMMARY_STATS_INTERVAL);
  // Default systematics functions:
  // - GetNumActive (taxa)
  // - GetTotalOrgs (total orgs tracked)
  // - GetAveDepth (average phylo depth of organisms)
  // - GetNumRoots
  // - GetMRCADepth
  // - CalcDiversity (entropy of taxa in population)
  // Functions to add:
  prog_gen_sys_file.template AddFun<size_t>([this]() { return mrca_changes; }, "mrca_changes", "MRCA changes");
  prog_gen_sys_file.AddStats(*prog_genotypic_systematics->GetDataNode("evolutionary_distinctiveness") , "evolutionary_distinctiveness", "evolutionary distinctiveness for a single update", true, true);
  prog_gen_sys_file.AddStats(*prog_genotypic_systematics->GetDataNode("pairwise_distances"), "pairwise_distance", "pairwise distance for a single update", true, true);
  // - GetPhylogeneticDiversity
  prog_gen_sys_file.AddCurrent(*prog_genotypic_systematics->GetDataNode("phylogenetic_diversity"), "current_phylogenetic_diversity", "current phylogenetic_diversity", true, true);
  // - GetTreeSize
  prog_gen_sys_file.template AddFun<size_t>([this]() { return prog_genotypic_systematics->GetTreeSize(); }, "tree_size", "Phylogenetic tree size");
  // - NumSparseTaxa
  prog_gen_sys_file.template AddFun<size_t>([this]() { return prog_genotypic_systematics->GetNumSparseTaxa(); }, "num_sparse_taxa", "Number sparse taxa");
  // - mean_sparse_pairwise_distances
  prog_gen_sys_file.template AddFun<size_t>([this]() { return prog_genotypic_systematics->GetMeanPairwiseDistance(true); }, "mean_sparse_pairwise_distances", "Number sparse taxa");
  // - sum_sparse_pairwise_distances
  prog_gen_sys_file.template AddFun<size_t>([this]() { return prog_genotypic_systematics->GetSumPairwiseDistance(true); }, "sum_sparse_pairwise_distances", "Number sparse taxa");
  // - variance_sparse_pairwise_distances
  prog_gen_sys_file.template AddFun<size_t>([this]() { return prog_genotypic_systematics->GetVariancePairwiseDistance(true); }, "variance_sparse_pairwise_distances", "Number sparse taxa");

  prog_gen_sys_file.AddFun(get_evaluations, "evaluations");

  prog_gen_sys_file.PrintHeaderKeys();

  // Add function(s) to program systematics snapshot
  prog_genotypic_systematics->AddSnapshotFun([](const prog_taxon_t & t) {
    std::ostringstream stream;
    t.GetInfo().PrintCSVEntry(stream);
    return stream.str();
  }, "program", "Program");

  // Track mrca in prog gen sys
  do_update_sig.AddAction([this]() {
    emp::Ptr<prog_taxon_t> cur_taxa = prog_genotypic_systematics->GetMRCA();
    if (cur_taxa != mrca_taxa_ptr) {
      ++mrca_changes;
      mrca_taxa_ptr = cur_taxa;
    }
  });

  // Setup prog_phen_diversity_file --> Gets updated during a snapshot, so we can assume that 
  prog_phen_diversity_file = emp::NewPtr<emp::DataFile>(DATA_DIRECTORY + "/prog_phenotype_diversity.csv");
  prog_phen_diversity_file->AddFun(get_update, "update");
  prog_phen_diversity_file->AddFun(get_evaluations, "evaluations");
  // - behavioral diversity
  prog_phen_diversity_file->AddFun(program_stats.get_prog_behavioral_diversity, "behavioral_diversity", "Shannon entropy of program behaviors");
  // - unique behavioral phenotypes
  prog_phen_diversity_file->AddFun(program_stats.get_prog_unique_behavioral_phenotypes, "unique_behavioral_phenotypes", "Unique behavioral profiles in program population");
  prog_phen_diversity_file->PrintHeaderKeys();

  // Setup solution file.
  solution_file = emp::NewPtr<emp::DataFile>(DATA_DIRECTORY + "/solutions.csv");

  solution_file->AddFun(get_update, "update");
  solution_file->AddFun(get_evaluations, "evaluations");
  solution_file->AddFun(program_stats.get_id, "program_id");
  solution_file->AddFun(program_stats.get_program_len, "program_len");
  solution_file->AddFun(program_stats.get_program, "program");
  solution_file->PrintHeaderKeys();

  do_pop_snapshot_sig.AddAction([this]() {
    SnapshotPrograms();
    SnapshotTests();
  });

  // Setup program/problem[X] fitness file.
  prog_world->SetupFitnessFile(DATA_DIRECTORY + "program_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL);
  // Setup test world fitness file (just don't look...)
  if (prob_NumberIO_world != nullptr) { prob_NumberIO_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_SmallOrLarge_world != nullptr) { prob_SmallOrLarge_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_ForLoopIndex_world != nullptr) { prob_ForLoopIndex_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_CompareStringLengths_world != nullptr) { prob_CompareStringLengths_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_DoubleLetters_world != nullptr) { prob_DoubleLetters_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_CollatzNumbers_world != nullptr) { prob_CollatzNumbers_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_ReplaceSpaceWithNewline_world != nullptr) { prob_ReplaceSpaceWithNewline_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_StringDifferences_world != nullptr) { prob_StringDifferences_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_EvenSquares_world != nullptr) { prob_EvenSquares_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_WallisPi_world != nullptr) { prob_WallisPi_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_StringLengthsBackwards_world != nullptr) { prob_StringLengthsBackwards_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_LastIndexOfZero_world != nullptr) { prob_LastIndexOfZero_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_VectorAverage_world != nullptr) { prob_VectorAverage_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_CountOdds_world != nullptr) { prob_CountOdds_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_MirrorImage_world != nullptr) { prob_MirrorImage_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_SuperAnagrams_world != nullptr) { prob_SuperAnagrams_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_SumOfSquares_world != nullptr) { prob_SumOfSquares_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_VectorsSummed_world != nullptr) { prob_VectorsSummed_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_XWordLines_world != nullptr) { prob_XWordLines_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_PigLatin_world != nullptr) { prob_PigLatin_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_NegativeToZero_world != nullptr) { prob_NegativeToZero_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_ScrabbleScore_world != nullptr) { prob_ScrabbleScore_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_Checksum_world != nullptr) { prob_Checksum_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_Digits_world != nullptr) { prob_Digits_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_Grade_world != nullptr) { prob_Grade_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_Median_world != nullptr) { prob_Median_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_Smallest_world != nullptr) { prob_Smallest_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else if (prob_Syllables_world != nullptr) { prob_Syllables_world->SetupFitnessFile(DATA_DIRECTORY + "test_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL); }
  else { std::cout << "AHH! None of the worlds have been initialized. Exiting." << std::endl; exit(-1); }

  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::COEVOLUTION || TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::STATIC_COEVO) {
    // Setup test world systematics
    if (prob_NumberIO_world != nullptr) { SetupTestSystematics(prob_NumberIO_world, [this](std::ostream & out, const prob_NumberIO_world_t::genome_t & genome) { prob_utils_NumberIO.PrintTestCSV(out, genome); } ); }
    else if (prob_SmallOrLarge_world != nullptr) { SetupTestSystematics(prob_SmallOrLarge_world, [this](std::ostream & out, const prob_SmallOrLarge_world_t::genome_t & genome) { prob_utils_SmallOrLarge.PrintTestCSV(out, genome); } ); }
    else if (prob_ForLoopIndex_world != nullptr) { SetupTestSystematics(prob_ForLoopIndex_world, [this](std::ostream & out, const prob_ForLoopIndex_world_t::genome_t & genome) { prob_utils_ForLoopIndex.PrintTestCSV(out, genome); } ); }
    else if (prob_CompareStringLengths_world != nullptr) { SetupTestSystematics(prob_CompareStringLengths_world, [this](std::ostream & out, const prob_CompareStringLengths_world_t::genome_t & genome) { prob_utils_CompareStringLengths.PrintTestCSV(out, genome); } ); }
    // else if (prob_DoubleLetters_world != nullptr) { SetupTestSystematics(prob_DoubleLetters_world); }
    else if (prob_CollatzNumbers_world != nullptr) { SetupTestSystematics(prob_CollatzNumbers_world, [this](std::ostream & out, const prob_CollatzNumbers_world_t::genome_t & genome) { prob_utils_CollatzNumbers.PrintTestCSV(out, genome); } ); }
    // else if (prob_ReplaceSpaceWithNewline_world != nullptr) { SetupTestSystematics(prob_ReplaceSpaceWithNewline_world); }
    // else if (prob_StringDifferences_world != nullptr) { SetupTestSystematics(prob_StringDifferences_world); }
    // else if (prob_EvenSquares_world != nullptr) { SetupTestSystematics(prob_EvenSquares_world); }
    // else if (prob_WallisPi_world != nullptr) { SetupTestSystematics(prob_WallisPi_world); }
    else if (prob_StringLengthsBackwards_world != nullptr) { SetupTestSystematics(prob_StringLengthsBackwards_world, [this](std::ostream & out, const prob_StringLengthsBackwards_world_t::genome_t & genome) { prob_utils_StringLengthsBackwards.PrintTestCSV(out, genome); } ); }
    else if (prob_LastIndexOfZero_world != nullptr) { SetupTestSystematics(prob_LastIndexOfZero_world, [this](std::ostream & out, const prob_LastIndexOfZero_world_t::genome_t & genome) { prob_utils_LastIndexOfZero.PrintTestCSV(out, genome); } ); }
    else if (prob_VectorAverage_world != nullptr) { SetupTestSystematics(prob_VectorAverage_world, [this](std::ostream & out, const prob_VectorAverage_world_t::genome_t & genome) { prob_utils_VectorAverage.PrintTestCSV(out, genome); } ); }
    else if (prob_CountOdds_world != nullptr) { SetupTestSystematics(prob_CountOdds_world, [this](std::ostream & out, const prob_CountOdds_world_t::genome_t & genome) { prob_utils_CountOdds.PrintTestCSV(out, genome); } ); }
    else if (prob_MirrorImage_world != nullptr) { SetupTestSystematics(prob_MirrorImage_world, [this](std::ostream & out, const prob_MirrorImage_world_t::genome_t & genome) { prob_utils_MirrorImage.PrintTestCSV(out, genome); } ); }
    // else if (prob_SuperAnagrams_world != nullptr) { SetupTestSystematics(prob_SuperAnagrams_world); }
    else if (prob_SumOfSquares_world != nullptr) { SetupTestSystematics(prob_SumOfSquares_world, [this](std::ostream & out, const prob_SumOfSquares_world_t::genome_t & genome) { prob_utils_SumOfSquares.PrintTestCSV(out, genome); } ); }
    else if (prob_VectorsSummed_world != nullptr) { SetupTestSystematics(prob_VectorsSummed_world, [this](std::ostream & out, const prob_VectorsSummed_world_t::genome_t & genome) { prob_utils_VectorsSummed.PrintTestCSV(out, genome); } ); }
    // else if (prob_XWordLines_world != nullptr) { SetupTestSystematics(prob_XWordLines_world); }
    // else if (prob_PigLatin_world != nullptr) { SetupTestSystematics(prob_PigLatin_world); }
    // else if (prob_NegativeToZero_world != nullptr) { SetupTestSystematics(prob_NegativeToZero_world); }
    // else if (prob_ScrabbleScore_world != nullptr) { SetupTestSystematics(prob_ScrabbleScore_world); }
    // else if (prob_Checksum_world != nullptr) { SetupTestSystematics(prob_Checksum_world); }
    // else if (prob_Digits_world != nullptr) { SetupTestSystematics(prob_Digits_world); }
    else if (prob_Grade_world != nullptr) { SetupTestSystematics(prob_Grade_world, [this](std::ostream & out, const prob_Grade_world_t::genome_t & genome) { prob_utils_Grade.PrintTestCSV(out, genome); } ); }
    else if (prob_Median_world != nullptr) { SetupTestSystematics(prob_Median_world, [this](std::ostream & out, const prob_Median_world_t::genome_t & genome) { prob_utils_Median.PrintTestCSV(out, genome); } ); }
    else if (prob_Smallest_world != nullptr) { SetupTestSystematics(prob_Smallest_world, [this](std::ostream & out, const prob_Smallest_world_t::genome_t & genome) { prob_utils_Smallest.PrintTestCSV(out, genome); } ); }
    // else if (prob_Syllables_world != nullptr) { SetupTestSystematics(prob_Syllables_world); }
    else { std::cout << "AHH! None of the worlds have been initialized. Exiting." << std::endl; exit(-1); }

  }

}

/// Setup program stats functions (used by data collection stuff).
void ProgramSynthesisExperiment::SetupProgramStats() {
  // Setup program stats functions.

  // program_stats.get_id
  program_stats.get_id = [this]() { return stats_util.cur_progID; };

  // program_stats.get_fitness
  program_stats.get_fitness = [this]() { return prog_world->CalcFitnessID(stats_util.cur_progID); };

  // program_stats.get_fitness_eval__total_score
  program_stats.get_fitness_eval__total_score = [this]() { return prog_world->GetOrg(stats_util.cur_progID).GetPhenotype().total_score; };
  
  // program_stats.get_fitness_eval__num_passes
  program_stats.get_fitness_eval__num_passes = [this]() { return prog_world->GetOrg(stats_util.cur_progID).GetPhenotype().num_passes; };

  // program_stats.get_fitness_eval__num_fails
  program_stats.get_fitness_eval__num_fails = [this]() { return prog_world->GetOrg(stats_util.cur_progID).GetPhenotype().num_fails; };

  // program_stats.get_fitness_eval__num_tests
  program_stats.get_fitness_eval__num_tests = [this]() { return prog_world->GetOrg(stats_util.cur_progID).GetPhenotype().test_passes.size(); };

  // program_stats.get_fitness_eval__passes_by_test
  program_stats.get_fitness_eval__passes_by_test = [this]() { 
    prog_org_t & prog = prog_world->GetOrg(stats_util.cur_progID);
    std::string scores = "\"[";
    for (size_t i = 0; i < prog.GetPhenotype().test_passes.size(); ++i) {
      if (i) scores += ",";
      scores += emp::to_string((size_t)prog.GetPhenotype().test_passes[i]);
    }
    scores += "]\"";
    return scores;
  };

  // program_stats.get_validation_eval__num_passes
  program_stats.get_validation_eval__num_passes = [this]() {
    return stats_util.current_program__validation__total_passes;
  };

  // program_stats.get_validation_eval__num_tests
  program_stats.get_validation_eval__num_tests = [this]() {
    return stats_util.current_program__validation__test_results.size();
  };

  // program_stats.get_validation_eval__passes_by_test
  program_stats.get_validation_eval__passes_by_test = [this]() {
    std::string scores = "\"[";
    for (size_t i = 0; i < stats_util.current_program__validation__test_results.size(); ++i) {
      if (i) scores += ",";
      scores += emp::to_string((size_t)stats_util.current_program__validation__test_results[i].pass);
    }
    scores += "]\"";
    return scores; 
  };

  // program_stats.get_program_len
  program_stats.get_program_len = [this]() {
    return prog_world->GetOrg(stats_util.cur_progID).GetGenome().GetSize();
  };

  // program_stats.get_program
  program_stats.get_program = [this]() {
    std::ostringstream stream;
    prog_world->GetOrg(stats_util.cur_progID).GetGenome().PrintCSVEntry(stream);
    return stream.str();
  };
}

/// Take a snapshot of the program population.
void ProgramSynthesisExperiment::SnapshotPrograms() {
  std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate());
  mkdir(snapshot_dir.c_str(), ACCESSPERMS);

  emp::DataFile file(snapshot_dir + "/program_pop_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");

  // Add functions to data file.
  file.AddFun(program_stats.get_id, "program_id", "Program ID");
  
  file.AddFun(program_stats.get_fitness, "fitness");
  file.AddFun(program_stats.get_fitness_eval__total_score, "total_score__fitness_eval");
  file.AddFun(program_stats.get_fitness_eval__num_passes, "num_passes__fitness_eval");
  file.AddFun(program_stats.get_fitness_eval__num_fails, "num_fails__fitness_eval");
  file.AddFun(program_stats.get_fitness_eval__num_tests, "num_tests__fitness_eval");
  file.AddFun(program_stats.get_fitness_eval__passes_by_test, "passes_by_test__fitness_eval");

  file.AddFun(program_stats.get_validation_eval__num_passes, "num_passes__validation_eval");
  file.AddFun(program_stats.get_validation_eval__num_tests, "num_tests__validation_eval");
  file.AddFun(program_stats.get_validation_eval__passes_by_test, "passes_by_test__validation_eval");

  file.AddFun(program_stats.get_program_len, "program_len");
  file.AddFun(program_stats.get_program, "program");

  file.PrintHeaderKeys();

  // For each program in the population, dump the program and anything we want to know about it.
  for (stats_util.cur_progID = 0; stats_util.cur_progID < prog_world->GetSize(); ++stats_util.cur_progID) {
    if (!prog_world->IsOccupied(stats_util.cur_progID)) continue;
    DoTestingSetValidation(prog_world->GetOrg(stats_util.cur_progID)); // Do validation for program.
    // Update snapshot file
    file.Update();
  }
  // Take diversity snapshot
  prog_phen_diversity_file->Update();
  // Snapshot phylogeny
  prog_genotypic_systematics->Snapshot(snapshot_dir + "/program_phylogeny_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");
}

// ================= PROGRAM-RELATED FUNCTIONS ===========
void ProgramSynthesisExperiment::InitProgPop_Random() {
  std::cout << "Randomly initializing program population." << std::endl;
  for (size_t i = 0; i < PROG_POP_SIZE; ++i) {
    prog_world->Inject(TagLGP::GenRandTagGPProgram(*random, inst_lib, MIN_PROG_SIZE, MAX_PROG_SIZE), 1);
  }
  // emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  // hardware_t::Program sol(inst_lib);

  // sol.PushInst("LoadThreshA",        {matrix[0], matrix[8], matrix[8]});
  // sol.PushInst("LoadThreshB",        {matrix[1], matrix[8], matrix[8]});
  // sol.PushInst("LoadThreshC",        {matrix[2], matrix[8], matrix[8]});
  // sol.PushInst("LoadThreshD",        {matrix[3], matrix[8], matrix[8]});
  // sol.PushInst("LoadGrade",          {matrix[4], matrix[8], matrix[8]});

  // sol.PushInst("TestNumGreaterTEqu", {matrix[4], matrix[0], matrix[5]});
  // sol.PushInst("If",                 {matrix[5], matrix[8], matrix[8]});
  //   sol.PushInst("SubmitA",          {matrix[8], matrix[8], matrix[8]});
  //   sol.PushInst("Return",           {matrix[8], matrix[8], matrix[8]});
  // sol.PushInst("Close",              {matrix[8], matrix[8], matrix[8]});

  // sol.PushInst("TestNumGreaterTEqu", {matrix[4], matrix[1], matrix[5]});
  // sol.PushInst("If",                 {matrix[5], matrix[8], matrix[8]});
  //   sol.PushInst("SubmitB",          {matrix[8], matrix[8], matrix[8]});
  //   sol.PushInst("Return",           {matrix[8], matrix[8], matrix[8]});
  // sol.PushInst("Close",              {matrix[8], matrix[8], matrix[8]});

  // sol.PushInst("TestNumGreaterTEqu", {matrix[4], matrix[2], matrix[5]});
  // sol.PushInst("If",                 {matrix[5], matrix[8], matrix[8]});
  //   sol.PushInst("SubmitC",          {matrix[8], matrix[8], matrix[8]});
  //   sol.PushInst("Return",           {matrix[8], matrix[8], matrix[8]});
  // sol.PushInst("Close",              {matrix[8], matrix[8], matrix[8]});

  // sol.PushInst("TestNumGreaterTEqu", {matrix[4], matrix[3], matrix[5]});
  // sol.PushInst("If",                 {matrix[5], matrix[8], matrix[8]});
  //   sol.PushInst("SubmitD",          {matrix[8], matrix[8], matrix[8]});
  //   sol.PushInst("Return",           {matrix[8], matrix[8], matrix[8]});
  // sol.PushInst("Close",              {matrix[8], matrix[8], matrix[8]});

  // sol.PushInst("SubmitF",            {matrix[8], matrix[8], matrix[8]});
  
  // prog_world->Inject(sol, PROG_POP_SIZE);
}

void ProgramSynthesisExperiment::AddDefaultInstructions(const std::unordered_set<std::string> & includes={"Add","Sub","Mult","Div","Mod",
                                                                                                  "TestNumEqu","TestNumNEqu","TestNumLess","TestNumLessTEqu","TestNumGreater","TestNumGreaterTEqu",
                                                                                                  "Floor","Not","Inc","Dec",
                                                                                                  "CopyMem","SwapMem","Input","Output","CommitGlobal","PullGlobal","TestMemEqu","TestMemNEqu",
                                                                                                  "MakeVector","VecGet","VecSet","VecLen","VecAppend","VecPop","VecRemove","VecReplaceAll","VecIndexOf","VecOccurrencesOf","VecReverse","VecSwapIfLess","VecGetFront","VecGetBack",
                                                                                                  "StrLength","StrConcat",
                                                                                                  "IsNum","IsStr","IsVec",
                                                                                                  "If","IfNot","While","Countdown","Foreach","Close","Break","Call","Routine","Return",
                                                                                                  "ModuleDef"}) 
{
  // Configure instructions
  // Math
  if (emp::Has(includes, "Add")) inst_lib->AddInst("Add", hardware_t::Inst_Add, 3, "wmemANY[C] = wmemNUM[A] + wmemNUM[B]");
  if (emp::Has(includes, "Sub")) inst_lib->AddInst("Sub", hardware_t::Inst_Sub, 3, "wmemANY[C] = wmemNUM[A] - wmemNUM[B]");
  if (emp::Has(includes, "Mult")) inst_lib->AddInst("Mult", hardware_t::Inst_Mult, 3, "wmemANY[C] = wmemNUM[A] * wmemNUM[B]");
  if (emp::Has(includes, "Div")) inst_lib->AddInst("Div", hardware_t::Inst_Div, 3, "if (wmemNUM[B] != 0) wmemANY[C] = wmemNUM[A] / wmemNUM[B]; else NOP");
  if (emp::Has(includes, "Mod")) inst_lib->AddInst("Mod", hardware_t::Inst_Mod, 3, "if (wmemNUM[B] != 0) wmemANY[C] = int(wmemNUM[A]) % int(wmemNUM[B]); else NOP");
  if (emp::Has(includes, "TestNumEqu")) inst_lib->AddInst("TestNumEqu", hardware_t::Inst_TestNumEqu, 3, "wmemANY[C] = wmemNUM[A] == wmemNUM[B]");
  if (emp::Has(includes, "TestNumNEqu")) inst_lib->AddInst("TestNumNEqu", hardware_t::Inst_TestNumNEqu, 3, "wmemANY[C] = wmemNUM[A] != wmemNUM[B]");
  if (emp::Has(includes, "TestNumLess")) inst_lib->AddInst("TestNumLess", hardware_t::Inst_TestNumLess, 3, "wmemANY[C] = wmemNUM[A] < wmemNUM[B]");
  if (emp::Has(includes, "TestNumLessTEqu")) inst_lib->AddInst("TestNumLessTEqu", hardware_t::Inst_TestNumLessTEqu, 3, "wmemANY[C] = wmemNUM[A] <= wmemNUM[B]");
  if (emp::Has(includes, "TestNumGreater")) inst_lib->AddInst("TestNumGreater", hardware_t::Inst_TestNumGreater, 3, "wmemANY[C] = wmemNUM[A] > wmemNUM[B]");
  if (emp::Has(includes, "TestNumGreaterTEqu")) inst_lib->AddInst("TestNumGreaterTEqu", hardware_t::Inst_TestNumGreaterTEqu, 3, "wmemANY[C] = wmemNUM[A] >= wmemNUM[B]");
  if (emp::Has(includes, "Floor")) inst_lib->AddInst("Floor", hardware_t::Inst_Floor, 1, "wmemNUM[A] = floor(wmemNUM[A])");
  if (emp::Has(includes, "Not")) inst_lib->AddInst("Not", hardware_t::Inst_Not, 1, "wmemNUM[A] = !wmemNUM[A]"); 
  if (emp::Has(includes, "Inc")) inst_lib->AddInst("Inc", hardware_t::Inst_Inc, 1, "wmemNUM[A] = wmemNUM[A] + 1");
  if (emp::Has(includes, "Dec")) inst_lib->AddInst("Dec", hardware_t::Inst_Dec, 1, "wmemNUM[A] = wmemNUM[A] - 1");
  
  // Memory manipulation
  if (emp::Has(includes, "CopyMem")) inst_lib->AddInst("CopyMem", hardware_t::Inst_CopyMem, 2, "wmemANY[B] = wmemANY[A] // Copy mem[A] to mem[B]");
  if (emp::Has(includes, "SwapMem")) inst_lib->AddInst("SwapMem", hardware_t::Inst_SwapMem, 2, "swap(wmemANY[A], wmemANY[B])");
  if (emp::Has(includes, "Input")) inst_lib->AddInst("Input", hardware_t::Inst_Input, 2, "wmemANY[B] = imemANY[A]");
  if (emp::Has(includes, "Output")) inst_lib->AddInst("Output", hardware_t::Inst_Output, 2, "omemANY[B] = wmemANY[A]");
  if (emp::Has(includes, "CommitGlobal")) inst_lib->AddInst("CommitGlobal", hardware_t::Inst_CommitGlobal, 2, "gmemANY[B] = wmemANY[A]");
  if (emp::Has(includes, "PullGlobal")) inst_lib->AddInst("PullGlobal", hardware_t::Inst_PullGlobal, 2, "wmemANY[B] = gmemANY[A]");
  if (emp::Has(includes, "TestMemEqu")) inst_lib->AddInst("TestMemEqu", hardware_t::Inst_TestMemEqu, 3, "wmemANY[C] = wmemANY[A] == wmemANY[B]");
  if (emp::Has(includes, "TestMemNEqu")) inst_lib->AddInst("TestMemNEqu", hardware_t::Inst_TestMemNEqu, 3, "wmemANY[C] = wmemANY[A] != wmemANY[B]");

  // Vector-related instructions
  if (emp::Has(includes, "MakeVector")) inst_lib->AddInst("MakeVector", hardware_t::Inst_MakeVector, 3, "wmemANY[C]=Vector([wmemANY[min(A,B),max(A,B)])");  // TODO - more descriptions
  if (emp::Has(includes, "VecGet")) inst_lib->AddInst("VecGet", hardware_t::Inst_VecGet, 3, "wmemANY[C]=wmemVEC[A][wmemNUM[B]]");
  if (emp::Has(includes, "VecSet")) inst_lib->AddInst("VecSet", hardware_t::Inst_VecSet, 3, "wmemVEC[A][wmemNUM[B]]=wmemNUM/STR[C]");
  if (emp::Has(includes, "VecLen")) inst_lib->AddInst("VecLen", hardware_t::Inst_VecLen, 2, "wmemANY[B]=wmemVEC[A]");
  if (emp::Has(includes, "VecAppend")) inst_lib->AddInst("VecAppend", hardware_t::Inst_VecAppend, 2, "wmemVEC[A].Append(wmemNUM/STR[B])");
  if (emp::Has(includes, "VecPop")) inst_lib->AddInst("VecPop", hardware_t::Inst_VecPop, 2, "wmemANY[B]=wmemVEC[A].pop()");
  if (emp::Has(includes, "VecRemove")) inst_lib->AddInst("VecRemove", hardware_t::Inst_VecRemove, 2, "wmemVEC[A].Remove(wmemNUM[B])");
  if (emp::Has(includes, "VecReplaceAll")) inst_lib->AddInst("VecReplaceAll", hardware_t::Inst_VecReplaceAll, 3, "Replace all values (wmemNUM/STR[B]) in wmemVEC[A] with wmemNUM/STR[C]");
  if (emp::Has(includes, "VecIndexOf")) inst_lib->AddInst("VecIndexOf", hardware_t::Inst_VecIndexOf, 3, "wmemANY[C] = index of wmemNUM/STR[B] in wmemVEC[A]");
  if (emp::Has(includes, "VecOccurrencesOf")) inst_lib->AddInst("VecOccurrencesOf", hardware_t::Inst_VecOccurrencesOf, 3, "wmemANY[C]= occurrances of wmemNUM/STR[B] in wmemVEC[A]");
  if (emp::Has(includes, "VecReverse")) inst_lib->AddInst("VecReverse", hardware_t::Inst_VecReverse, 1, "wmemVEC[A] = Reverse(wmemVEC[A])");
  if (emp::Has(includes, "VecSwapIfLess")) inst_lib->AddInst("VecSwapIfLess", hardware_t::Inst_VecSwapIfLess, 3, "Swap two indices in wmemVEC[A] if vec[wmemNUM[A]] < vec[wmemNUM[B]].");
  if (emp::Has(includes, "VecGetFront")) inst_lib->AddInst("VecGetFront", hardware_t::Inst_VecGetFront, 2, "wmemANY[B] = front of wmemVEC[A]");
  if (emp::Has(includes, "VecGetBack")) inst_lib->AddInst("VecGetBack", hardware_t::Inst_VecGetBack, 2, "wmemANY[B] = back of wmemVEC[A]");
  
  // String
  if (emp::Has(includes, "StrLength")) inst_lib->AddInst("StrLength", hardware_t::Inst_StrLength, 2, "todo");
  if (emp::Has(includes, "StrConcat")) inst_lib->AddInst("StrConcat", hardware_t::Inst_StrConcat, 3, "todo");

  // Memory-type
  if (emp::Has(includes, "IsNum")) inst_lib->AddInst("IsNum", hardware_t::Inst_IsNum, 2, "wmemANY[B] = IsNum(wmemANY[A])");
  if (emp::Has(includes, "IsStr")) inst_lib->AddInst("IsStr", hardware_t::Inst_IsStr, 2, "wmemANY[B] = IsStr(wmemANY[A])");
  if (emp::Has(includes, "IsVec")) inst_lib->AddInst("IsVec", hardware_t::Inst_IsVec, 2, "wmemANY[B] = IsVec(wmemANY[A])");

  // Flow control
  if (emp::Has(includes, "If")) inst_lib->AddInst("If", hardware_t::Inst_If, 1, "Execute next flow if(wmemANY[A]) // To be true, mem loc must be non-zero number", {inst_lib_t::InstProperty::BEGIN_FLOW});
  if (emp::Has(includes, "IfNot")) inst_lib->AddInst("IfNot", hardware_t::Inst_IfNot, 1, "Execute next flow if(!wmemANY[A])", {inst_lib_t::InstProperty::BEGIN_FLOW});
  if (emp::Has(includes, "While")) inst_lib->AddInst("While", hardware_t::Inst_While, 1, "While loop over wmemANY[A]", {inst_lib_t::InstProperty::BEGIN_FLOW});
  if (emp::Has(includes, "Countdown")) inst_lib->AddInst("Countdown", hardware_t::Inst_Countdown, 1, "Countdown loop with wmemANY as index.", {inst_lib_t::InstProperty::BEGIN_FLOW});
  if (emp::Has(includes, "Foreach")) inst_lib->AddInst("Foreach", hardware_t::Inst_Foreach, 2, "For each thing in wmemVEC[B]", {inst_lib_t::InstProperty::BEGIN_FLOW});
  if (emp::Has(includes, "Close")) inst_lib->AddInst("Close", hardware_t::Inst_Close, 0, "Close flow", {inst_lib_t::InstProperty::END_FLOW});
  if (emp::Has(includes, "Break")) inst_lib->AddInst("Break", hardware_t::Inst_Break, 0, "Break current flow");
  if (emp::Has(includes, "Call")) inst_lib->AddInst("Call", hardware_t::Inst_Call, 1, "Call module using A for tag-based reference");
  if (emp::Has(includes, "Routine")) inst_lib->AddInst("Routine", hardware_t::Inst_Routine, 1, "Call module as a routine (don't use call stack)");
  if (emp::Has(includes, "Return")) inst_lib->AddInst("Return", hardware_t::Inst_Return, 0, "Return from current routine/call");

  // Module
  if (emp::Has(includes, "ModuleDef")) inst_lib->AddInst("ModuleDef", hardware_t::Inst_Nop, 1, "Define module with tag A", {inst_lib_t::InstProperty::MODULE});

  // Misc
  // inst_lib->AddInst("Nop", hardware_t::Inst_Nop, 3, "Do nothing");
}

// ================= PROBLEM SETUPS ======================

void ProgramSynthesisExperiment::SetupProblem_NumberIO() { 
  std::cout << "NumberIO problem setup" << std::endl; 

  using test_org_t = TestOrg_NumberIO;
  
  // Load testing examples from file (used to evaluate 'true' performance of programs).
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  prob_utils_NumberIO.GetTrainingSet().LoadTestCases(training_examples_fpath);
  prob_utils_NumberIO.GetTestingSet().LoadTestCases(testing_examples_fpath);
  prob_utils_NumberIO.GenerateTestingSetPop();
  std::cout << "Loaded training example set size = " << prob_utils_NumberIO.GetTrainingSet().GetSize() << std::endl;
  std::cout << "Loaded testing example set size = " << prob_utils_NumberIO.GetTestingSet().GetSize() << std::endl;
  std::cout << "Testing set (non-training examples used to evaluate program accuracy) size = " << prob_utils_NumberIO.testingset_pop.size() << std::endl;

  // Setup world.
  NewTestCaseWorld(prob_NumberIO_world, *random, "NumberIO Test Case World");
  
  // Configure how population should be initialized
  SetupTestCasePop_Init(prob_NumberIO_world, 
                          prob_utils_NumberIO.training_set,
                          [this]() { return GenRandomTestInput_NumberIO(*random, {PROB_NUMBER_IO__INT_MIN, PROB_NUMBER_IO__INT_MAX}, {PROB_NUMBER_IO__DOUBLE_MIN, PROB_NUMBER_IO__DOUBLE_MAX}); } );
  end_setup_sig.AddAction([this]() { std::cout << "TestCase world size = " << prob_NumberIO_world->GetSize() << std::endl; });
  
  // Tell world to calculate correct test output (given input) on placement.
  prob_NumberIO_world->OnPlacement([this](size_t pos) { prob_NumberIO_world->GetOrg(pos).CalcOut(); } );

  EvaluateWorldTest = [this](prog_org_t & prog_org, size_t testID) {
    emp::Ptr<test_org_t> test_org_ptr = prob_NumberIO_world->GetOrgPtr(testID);
    begin_program_test.Trigger(prog_org, test_org_ptr);
    do_program_test.Trigger(prog_org, test_org_ptr);
    end_program_test.Trigger(prog_org, test_org_ptr);
    return CalcProgramResultOnTest(prog_org, *test_org_ptr);
  };

  prob_utils_NumberIO.population_validation_outputs.resize(PROG_POP_SIZE);
  DoTestingSetValidation = [this](prog_org_t & prog_org) {
    // evaluate program on full testing set; update stats utils with results
    begin_program_eval.Trigger(prog_org);
    stats_util.current_program__validation__test_results.resize(prob_utils_NumberIO.testingset_pop.size());
    stats_util.current_program__validation__total_score = 0;
    stats_util.current_program__validation__total_passes = 0;
    stats_util.current_program__validation__is_solution = false;
    prob_utils_NumberIO.population_validation_outputs[stats_util.cur_progID].resize(prob_utils_NumberIO.testingset_pop.size());
    // For each test in validation set, evaluate program.
    for (size_t testID = 0; testID < prob_utils_NumberIO.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_NumberIO.testingset_pop[testID];
      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      stats_util.current_program__validation__test_results[testID] = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      stats_util.current_program__validation__total_score += stats_util.current_program__validation__test_results[testID].score;
      stats_util.current_program__validation__total_passes += (size_t)stats_util.current_program__validation__test_results[testID].pass;
      prob_utils_NumberIO.population_validation_outputs[stats_util.cur_progID][testID] = prob_utils_NumberIO.submitted_val;
    }
    stats_util.current_program__validation__is_solution = stats_util.current_program__validation__total_passes == prob_utils_NumberIO.testingset_pop.size();
    end_program_eval.Trigger(prog_org);
  }; // todo - test this out
  program_stats.get_prog_behavioral_diversity = [this]() { return emp::ShannonEntropy(prob_utils_NumberIO.population_validation_outputs); };
  program_stats.get_prog_unique_behavioral_phenotypes = [this]() { return emp::UniqueCount(prob_utils_NumberIO.population_validation_outputs); };

  ScreenForSolution = [this](prog_org_t & prog_org) {
    begin_program_eval.Trigger(prog_org);
    for (size_t testID = 0; testID < prob_utils_NumberIO.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_NumberIO.testingset_pop[testID];

      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      
      TestResult result = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      if (!result.pass) {
        end_program_eval.Trigger(prog_org);
        return false;
      }
    }
    end_program_eval.Trigger(prog_org);
    return true;
  };

  // Tell experiment how to get test phenotypes.
  GetTestPhenotype = [this](size_t testID) -> test_org_phen_t & {
    emp_assert(prob_NumberIO_world->IsOccupied(testID));
    return prob_NumberIO_world->GetOrg(testID).GetPhenotype();
  };

  // Setup how test case world updates.  
  SetupTestCaseWorldUpdate(prob_NumberIO_world);
  
  // Setup how test cases mutate.
  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::RANDOM) {
    std::cout << "RANDOM training mode detected, configuring mutation function to RANDOMIZE organisms." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_NumberIO.mutator.MIN_INT = PROB_NUMBER_IO__INT_MIN;
      prob_utils_NumberIO.mutator.MAX_INT = PROB_NUMBER_IO__INT_MAX;
      prob_utils_NumberIO.mutator.MIN_DOUBLE = PROB_NUMBER_IO__DOUBLE_MIN;
      prob_utils_NumberIO.mutator.MAX_DOUBLE = PROB_NUMBER_IO__DOUBLE_MAX;
      prob_utils_NumberIO.mutator.PER_INT_RATE = 1.0;
      prob_utils_NumberIO.mutator.PER_DOUBLE_RATE = 1.0;
      // (2) Hook mutator up to world.
      prob_NumberIO_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_NumberIO.mutator.Mutate(rnd, test_org.GetGenome());
      });
    };
  } else {
    std::cout << "Non-RANDOM training mode detected, configuring mutation function normally." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_NumberIO.mutator.MIN_INT = PROB_NUMBER_IO__INT_MIN;
      prob_utils_NumberIO.mutator.MAX_INT = PROB_NUMBER_IO__INT_MAX;
      prob_utils_NumberIO.mutator.MIN_DOUBLE = PROB_NUMBER_IO__DOUBLE_MIN;
      prob_utils_NumberIO.mutator.MAX_DOUBLE = PROB_NUMBER_IO__DOUBLE_MAX;
      prob_utils_NumberIO.mutator.PER_INT_RATE = PROB_NUMBER_IO__MUTATION__PER_INT_RATE;
      prob_utils_NumberIO.mutator.PER_DOUBLE_RATE = PROB_NUMBER_IO__MUTATION__PER_DOUBLE_RATE;
      // (2) Hook mutator up to world.
      prob_NumberIO_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_NumberIO.mutator.Mutate(rnd, test_org.GetGenome());
      });
    };
  }

  // Setup test case fitness function.
  SetupTestFitFun = [this]() {
    prob_NumberIO_world->SetFitFun([](test_org_t & test_org) {
      return (double)test_org.GetPhenotype().num_fails;
    });
  };
  
  // Tell experiment how to configure hardware inputs when running program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org, emp::Ptr<TestOrg_Base> test_org_base_ptr) {
    // Reset eval stuff
    // Set current test org.
    // prob_utils_NumberIO.cur_eval_test_org = prob_NumberIO_world->GetOrgPtr(testID); // currently only place need testID for this?
    prob_utils_NumberIO.cur_eval_test_org = test_org_base_ptr.Cast<test_org_t>(); // currently only place need testID for this?
    prob_utils_NumberIO.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 2);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      // Grab some useful references.
      Problem_NumberIO_input_t & input = prob_utils_NumberIO.cur_eval_test_org->GetGenome(); // std::pair<int, double>
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      // Set hardware input.
      wmem.Set(0, input.first);
      wmem.Set(1, input.second);
    } // TODO - confirm that this is working.
  });

  // Tell experiment how to calculate program score.
  CalcProgramResultOnTest = [this](prog_org_t & prog_org, TestOrg_Base & test_org_base) {
    // std::cout << "Calc score on test!" << std::endl;
    test_org_t & test_org = static_cast<test_org_t&>(test_org_base);
    TestResult result;
    if (!prob_utils_NumberIO.submitted) {
      result.score = 0;
      result.pass = false;
      result.sub = false;
    } else {
      double max_error = emp::Abs(PROB_NUMBER_IO__DOUBLE_MAX) * 2;
      std::pair<double, bool> r(CalcScoreGradient_NumberIO(test_org.GetCorrectOut(), prob_utils_NumberIO.submitted_val, max_error));
      result.score = r.first;
      result.pass = r.second;
      result.sub = true;
    }
    return result;
  };

  SnapshotTests = [this]() {
    std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate());
    mkdir(snapshot_dir.c_str(), ACCESSPERMS);
    
    emp::DataFile file(snapshot_dir + "/test_pop_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");
    // Test file contents:
    // - test id
    std::function<size_t(void)> get_test_id = [this]() { return stats_util.cur_testID; };
    file.AddFun(get_test_id, "test_id");

    // - test fitness
    std::function<double(void)> get_test_fitness = [this]() { return prob_NumberIO_world->CalcFitnessID(stats_util.cur_testID); };
    file.AddFun(get_test_fitness, "fitness");

    // - num passes
    std::function<size_t(void)> get_test_num_passes = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_passes; };
    file.AddFun(get_test_num_passes, "num_passes");

    // - num fails
    std::function<size_t(void)> get_test_num_fails = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_fails; };
    file.AddFun(get_test_num_fails, "num_fails");

    std::function<size_t(void)> get_num_tested = [this]() { return GetTestPhenotype(stats_util.cur_testID).test_passes.size(); };
    file.AddFun(get_num_tested, "num_programs_tested_against");

    // - test scores by program
    std::function<std::string(void)> get_passes_by_program = [this]() {
      std::string scores = "\"[";
      test_org_phen_t & phen = GetTestPhenotype(stats_util.cur_testID);
      for (size_t i = 0; i < phen.test_passes.size(); ++i) {
        if (i) scores += ",";
        scores += emp::to_string(phen.test_passes[i]);
      }
      scores += "]\"";
      return scores;
    };
    file.AddFun(get_passes_by_program, "passes_by_program");

    // - test
    std::function<std::string(void)> get_test = [this]() {
      std::ostringstream stream;
      stream << "\"";
      prob_NumberIO_world->GetOrg(stats_util.cur_testID).Print(stream);
      stream << "\"";
      return stream.str();
    };
    file.AddFun(get_test, "test", "");

    file.PrintHeaderKeys();

    // Loop over tests, snapshotting each.
    for (stats_util.cur_testID = 0; stats_util.cur_testID < prob_NumberIO_world->GetSize(); ++stats_util.cur_testID) {
      if (!prob_NumberIO_world->IsOccupied(stats_util.cur_testID)) continue;
      file.Update();
    }

    //->Snapshot(snapshot_dir + "/test_phylogeny_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");

  };
  
  AddDefaultInstructions({"Add",
                          "Sub",
                          "Mult",
                          "Div",
                          "Mod",
                          "TestNumEqu",
                          "TestNumNEqu",
                          "TestNumLess",
                          "TestNumLessTEqu",
                          "TestNumGreater",
                          "TestNumGreaterTEqu",
                          "Floor",
                          "Not",
                          "Inc",
                          "Dec",
                          "CopyMem",
                          "SwapMem",
                          "Input",
                          "Output",
                          "CommitGlobal",
                          "PullGlobal",
                          "TestMemEqu",
                          "TestMemNEqu",
                          "If",
                          "IfNot",
                          "While",
                          "Countdown",
                          "Foreach",
                          "Close",
                          "Break",
                          "Call",
                          "Routine",
                          "Return",
                          "ModuleDef"
  });

  // Add Terminals [0:16] -- TODO - may want these to be slightly more configurable.
  for (size_t i = 0; i <= 16; ++i) {
    inst_lib->AddInst("Set-" + emp::to_string(i),
      [i](hardware_t & hw, const inst_t & inst) {
        hardware_t::CallState & state = hw.GetCurCallState();
        hardware_t::Memory & wmem = state.GetWorkingMem();
        size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
        if (!hw.IsValidMemPos(posA)) return; // Do nothing
        wmem.Set(posA, (double)i);
      });
  }

  // Add custom instructions
  // - LoadInteger
  inst_lib->AddInst("LoadInt", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadInt_NumberIO(hw, inst);
  }, 1);
  // - LoadDouble
  inst_lib->AddInst("LoadDouble", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadDouble_NumberIO(hw, inst);
  }, 1);
  // - SubmitNum
  inst_lib->AddInst("SubmitNum", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitNum_NumberIO(hw, inst);
  }, 1);
}

void ProgramSynthesisExperiment::SetupProblem_SmallOrLarge() { 
  std::cout << "Setup problem SmallOrLarge." << std::endl;

  // A few useful aliases.
  using test_org_t = TestOrg_SmallOrLarge;
  // using problem_utils_t = ProblemUtilities_SmallOrLarge;
  // using problem_world_t = prob_SmallOrLarge_world_t;

  // A few useful references.
  // problem_utils_t & problem_utils = prob_utils_SmallOrLarge;
  // emp::Ptr<problem_world_t> problem_world = prob_SmallOrLarge_world;

  // Load testing examples from file.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  prob_utils_SmallOrLarge.GetTrainingSet().LoadTestCases(training_examples_fpath);
  prob_utils_SmallOrLarge.GetTestingSet().LoadTestCases(testing_examples_fpath);
  prob_utils_SmallOrLarge.GenerateTestingSetPop();
  std::cout << "Loaded training example set size = " << prob_utils_SmallOrLarge.GetTrainingSet().GetSize() << std::endl;
  std::cout << "Loaded testing example set size = " << prob_utils_SmallOrLarge.GetTestingSet().GetSize() << std::endl;
  std::cout << "Testing set (non-training examples used to evaluate program accuracy) size = " << prob_utils_SmallOrLarge.testingset_pop.size() << std::endl;

  // Setup the world.
  NewTestCaseWorld(prob_SmallOrLarge_world, *random, "SmallOrLarge world");

  // Configure how the population should be initialized.
  SetupTestCasePop_Init(prob_SmallOrLarge_world,
                        prob_utils_SmallOrLarge.training_set,
                        [this]() { return GenRandomTestInput_SmallOrLarge(*random, {PROB_SMALL_OR_LARGE__INT_MIN, PROB_SMALL_OR_LARGE__INT_MAX}); });

  end_setup_sig.AddAction([this]() { std::cout << "TestCase world size= " << prob_SmallOrLarge_world->GetSize() << std::endl; });

  // Tell the world to calculate the correct test output (given input) on placement.
  prob_SmallOrLarge_world->OnPlacement([this](size_t pos) { prob_SmallOrLarge_world->GetOrg(pos).CalcOut(); });

  // How are program results calculated on a test?
  CalcProgramResultOnTest = [this](prog_org_t & prog_org, TestOrg_Base & test_org_base) {
    test_org_t & test_org = static_cast<test_org_t&>(test_org_base);
    TestResult result;
    if (!prob_utils_SmallOrLarge.submitted) {
      result.score = 0;
      result.pass = false;
      result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_SmallOrLarge.CalcScorePassFail(test_org.GetCorrectOut(), prob_utils_SmallOrLarge.submitted_str));
      result.score = r.first;
      result.pass = r.second;
      result.sub = true;
    }
    return result;
  };
  
  // Setup how evaluation on world test should work.
  EvaluateWorldTest = [this](prog_org_t & prog_org, size_t testID) {
    emp::Ptr<test_org_t> test_org_ptr = prob_SmallOrLarge_world->GetOrgPtr(testID);
    begin_program_test.Trigger(prog_org, test_org_ptr);
    do_program_test.Trigger(prog_org, test_org_ptr);
    end_program_test.Trigger(prog_org, test_org_ptr);
    return CalcProgramResultOnTest(prog_org, *test_org_ptr);
  };

  // How should we validate programs on testing set?
  prob_utils_SmallOrLarge.population_validation_outputs.resize(PROG_POP_SIZE);
  DoTestingSetValidation = [this](prog_org_t & prog_org) { 
    // evaluate program on full testing set; update stats utils with results
    begin_program_eval.Trigger(prog_org);
    stats_util.current_program__validation__test_results.resize(prob_utils_SmallOrLarge.testingset_pop.size());
    stats_util.current_program__validation__total_score = 0;
    stats_util.current_program__validation__total_passes = 0;
    stats_util.current_program__validation__is_solution = false;
    prob_utils_SmallOrLarge.population_validation_outputs[stats_util.cur_progID].resize(prob_utils_SmallOrLarge.testingset_pop.size());
    // For each test in validation set, evaluate program.
    for (size_t testID = 0; testID < prob_utils_SmallOrLarge.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_SmallOrLarge.testingset_pop[testID];
      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      stats_util.current_program__validation__test_results[testID] = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      stats_util.current_program__validation__total_score += stats_util.current_program__validation__test_results[testID].score;
      stats_util.current_program__validation__total_passes += (size_t)stats_util.current_program__validation__test_results[testID].pass;
      prob_utils_SmallOrLarge.population_validation_outputs[stats_util.cur_progID][testID] = prob_utils_SmallOrLarge.submitted_str;
    }
    stats_util.current_program__validation__is_solution = stats_util.current_program__validation__total_passes == prob_utils_SmallOrLarge.testingset_pop.size();
    end_program_eval.Trigger(prog_org);
  };

  program_stats.get_prog_behavioral_diversity = [this]() { return emp::ShannonEntropy(prob_utils_SmallOrLarge.population_validation_outputs); };
  program_stats.get_prog_unique_behavioral_phenotypes = [this]() { return emp::UniqueCount(prob_utils_SmallOrLarge.population_validation_outputs); };

  // How should we screen for a solution?
  ScreenForSolution = [this](prog_org_t & prog_org) {
    begin_program_eval.Trigger(prog_org);
    for (size_t testID = 0; testID < prob_utils_SmallOrLarge.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_SmallOrLarge.testingset_pop[testID];

      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      
      TestResult result = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      if (!result.pass) {
        end_program_eval.Trigger(prog_org);
        return false;
      }
    }
    end_program_eval.Trigger(prog_org);
    return true;
  };

  // Tell the experiment how to get test phenotypes.
  GetTestPhenotype = [this](size_t testID) -> test_org_phen_t & {
    emp_assert(prob_SmallOrLarge_world->IsOccupied(testID));
    return prob_SmallOrLarge_world->GetOrg(testID).GetPhenotype();
  };

  // Setup how test world updates.
  SetupTestCaseWorldUpdate(prob_SmallOrLarge_world);

  // Setup how test cases mutate.
  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::RANDOM) {
    std::cout << "RANDOM training mode detected, configuring mutation function to RANDOMIZE organisms." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_SmallOrLarge.mutator.MIN_INT = PROB_SMALL_OR_LARGE__INT_MIN;
      prob_utils_SmallOrLarge.mutator.MAX_INT = PROB_SMALL_OR_LARGE__INT_MAX;
      prob_utils_SmallOrLarge.mutator.PER_INT_RATE = 1.0;
      // (2) Hook mutator up to world.
      prob_SmallOrLarge_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_SmallOrLarge.mutator.Mutate(rnd, test_org.GetGenome());
      });
    };
  } else {
    std::cout << "Non-RANDOM training mode detected, configuring mutation function normally." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_SmallOrLarge.mutator.MIN_INT = PROB_SMALL_OR_LARGE__INT_MIN;
      prob_utils_SmallOrLarge.mutator.MAX_INT = PROB_SMALL_OR_LARGE__INT_MAX;
      prob_utils_SmallOrLarge.mutator.PER_INT_RATE = PROB_SMALL_OR_LARGE__MUTATION__PER_INT_RATE;
      // (2) Hook mutator up to world.
      prob_SmallOrLarge_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_SmallOrLarge.mutator.Mutate(rnd, test_org.GetGenome());
      });
    };
  }

  // Setup test case fitness function.
  SetupTestFitFun = [this]() {
    prob_SmallOrLarge_world->SetFitFun([](test_org_t & test_org) {
      return (double)test_org.GetPhenotype().num_fails;
    });
  };

  // Tell experiment how to configure hardware inputs when running program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org, emp::Ptr<TestOrg_Base> test_org_base_ptr) {
    // Reset eval stuff
    // Set current test org.
    prob_utils_SmallOrLarge.cur_eval_test_org = test_org_base_ptr.Cast<test_org_t>(); // currently only place need testID for this?
    prob_utils_SmallOrLarge.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 1);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      // Grab some useful references.
      Problem_SmallOrLarge_input_t & input = prob_utils_SmallOrLarge.cur_eval_test_org->GetGenome(); // std::pair<int, double>
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      // Set hardware input.
      wmem.Set(0, input);
    }
  });

  // Tell experiment how to snapshot test population.
  SnapshotTests = [this]() {
    std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate());
    mkdir(snapshot_dir.c_str(), ACCESSPERMS);
    
    emp::DataFile file(snapshot_dir + "/test_pop_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");
    // Test file contents:
    // - test id
    std::function<size_t(void)> get_test_id = [this]() { return stats_util.cur_testID; };
    file.AddFun(get_test_id, "test_id");

    // - test fitness
    std::function<double(void)> get_test_fitness = [this]() { return prob_SmallOrLarge_world->CalcFitnessID(stats_util.cur_testID); };
    file.AddFun(get_test_fitness, "fitness");

    // - num passes
    std::function<size_t(void)> get_test_num_passes = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_passes; };
    file.AddFun(get_test_num_passes, "num_passes");

    // - num fails
    std::function<size_t(void)> get_test_num_fails = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_fails; };
    file.AddFun(get_test_num_fails, "num_fails");

    std::function<size_t(void)> get_num_tested = [this]() { return GetTestPhenotype(stats_util.cur_testID).test_passes.size(); };
    file.AddFun(get_num_tested, "num_programs_tested_against");

    // - test scores by program
    std::function<std::string(void)> get_passes_by_program = [this]() {
      std::string scores = "\"[";
      test_org_phen_t & phen = GetTestPhenotype(stats_util.cur_testID);
      for (size_t i = 0; i < phen.test_passes.size(); ++i) {
        if (i) scores += ",";
        scores += emp::to_string(phen.test_passes[i]);
      }
      scores += "]\"";
      return scores;
    };
    file.AddFun(get_passes_by_program, "passes_by_program");

    // - test
    std::function<std::string(void)> get_test = [this]() {
      std::ostringstream stream;
      stream << "\"";
      prob_SmallOrLarge_world->GetOrg(stats_util.cur_testID).Print(stream);
      stream << "\"";
      return stream.str();
    };
    file.AddFun(get_test, "test", "");

    file.PrintHeaderKeys();

    // Loop over tests, snapshotting each.
    for (stats_util.cur_testID = 0; stats_util.cur_testID < prob_SmallOrLarge_world->GetSize(); ++stats_util.cur_testID) {
      if (!prob_SmallOrLarge_world->IsOccupied(stats_util.cur_testID)) continue;
      file.Update();
    }
  };

  // Add default instructions to instruction set.
  AddDefaultInstructions({"Add",
                          "Sub",
                          "Mult",
                          "Div",
                          "Mod",
                          "TestNumEqu",
                          "TestNumNEqu",
                          "TestNumLess",
                          "TestNumLessTEqu",
                          "TestNumGreater",
                          "TestNumGreaterTEqu",
                          "Floor",
                          "Not",
                          "Inc",
                          "Dec",
                          "CopyMem",
                          "SwapMem",
                          "Input",
                          "Output",
                          "CommitGlobal",
                          "PullGlobal",
                          "TestMemEqu",
                          "TestMemNEqu",
                          "If",
                          "IfNot",
                          "While",
                          "Countdown",
                          "Foreach",
                          "Close",
                          "Break",
                          "Call",
                          "Routine",
                          "Return",
                          "ModuleDef"
  });

  // -- Custom instructions --
  // Add Terminals [0:16] -- TODO - may want these to be slightly more configurable.
  for (size_t i = 0; i <= 10; ++i) {
    inst_lib->AddInst("Set-" + emp::to_string(i),
      [i](hardware_t & hw, const inst_t & inst) {
        hardware_t::CallState & state = hw.GetCurCallState();
        hardware_t::Memory & wmem = state.GetWorkingMem();
        size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
        if (!hw.IsValidMemPos(posA)) return; // Do nothing
        wmem.Set(posA, (double)i);
      });
  }

  // - LoadInteger
  inst_lib->AddInst("LoadInt", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadInt_SmallOrLarge(hw, inst);
  }, 1);

  inst_lib->AddInst("SubmitSmall", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitSmall_SmallOrLarge(hw, inst);
  }, 0);

  inst_lib->AddInst("SubmitLarge", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitLarge_SmallOrLarge(hw, inst);
  }, 0);

  inst_lib->AddInst("SubmitNone", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitNone_SmallOrLarge(hw, inst);
  }, 0);

}

void ProgramSynthesisExperiment::SetupProblem_ForLoopIndex() { 
  std::cout << "Setup problem - ForLoopIndex" << std::endl;

  // A few useful aliases.
  using test_org_t = TestOrg_ForLoopIndex;

  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  std::cout << "Loading training examples." << std::endl;
  prob_utils_ForLoopIndex.GetTrainingSet().LoadTestCases(training_examples_fpath);
  std::cout << "Loading testing examples." << std::endl;
  prob_utils_ForLoopIndex.GetTestingSet().LoadTestCases(testing_examples_fpath);
  std::cout << "Generating testing set population." << std::endl;
  prob_utils_ForLoopIndex.GenerateTestingSetPop();
  std::cout << "Loaded training example set size = " << prob_utils_ForLoopIndex.GetTrainingSet().GetSize() << std::endl;
  std::cout << "Loaded testing example set size = " << prob_utils_ForLoopIndex.GetTestingSet().GetSize() << std::endl;
  std::cout << "Testing set (non-training examples used to evaluate program accuracy) size = " << prob_utils_ForLoopIndex.testingset_pop.size() << std::endl;

  // Setup the world
  NewTestCaseWorld(prob_ForLoopIndex_world, *random, "ForLoopIndex world");

  // Configure how the population should be initialized
  SetupTestCasePop_Init(prob_ForLoopIndex_world,
                        prob_utils_ForLoopIndex.training_set,
                        [this]() { return GenRandomTestInput_ForLoopIndex(*random, 
                                                                          {PROB_FOR_LOOP_INDEX__START_END_MIN, PROB_FOR_LOOP_INDEX__START_END_MAX},
                                                                          {PROB_FOR_LOOP_INDEX__STEP_MIN, PROB_FOR_LOOP_INDEX__STEP_MAX},
                                                                          PROB_FOR_LOOP_INDEX__PROMISE_MULTISTEP_TESTCASES
                                                                         ); 
                                  }
                        );
  end_setup_sig.AddAction([this]() { std::cout << "TestCase world size= " << prob_ForLoopIndex_world->GetSize() << std::endl; });

  // Tell the world to calculate the correct test output (given input) on placement.
  prob_ForLoopIndex_world->OnPlacement([this](size_t pos) { prob_ForLoopIndex_world->GetOrg(pos).CalcOut(); });

  // How are program results calculated on a test?
  CalcProgramResultOnTest = [this](prog_org_t & prog_org, TestOrg_Base & test_org_base) {
    test_org_t & test_org = static_cast<test_org_t&>(test_org_base);
    TestResult result;
    if (!prob_utils_ForLoopIndex.submitted) {
      result.score = 0;
      result.pass = false;
      result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_ForLoopIndex.CalcScoreGradient(test_org.GetCorrectOut(), prob_utils_ForLoopIndex.submitted_vec));
      result.score = r.first;
      result.pass = r.second;
      result.sub = true;
    }
    return result;
  };

  // Setup how evaluation on world test should work.
  EvaluateWorldTest = [this](prog_org_t & prog_org, size_t testID) {
    emp::Ptr<test_org_t> test_org_ptr = prob_ForLoopIndex_world->GetOrgPtr(testID);
    begin_program_test.Trigger(prog_org, test_org_ptr);
    do_program_test.Trigger(prog_org, test_org_ptr);
    end_program_test.Trigger(prog_org, test_org_ptr);
    return CalcProgramResultOnTest(prog_org, *test_org_ptr);
  };

  // How should we validate programs on testing set?
  prob_utils_ForLoopIndex.population_validation_outputs.resize(PROG_POP_SIZE);
  DoTestingSetValidation = [this](prog_org_t & prog_org) { 
    // evaluate program on full testing set; update stats utils with results
    begin_program_eval.Trigger(prog_org);
    stats_util.current_program__validation__test_results.resize(prob_utils_ForLoopIndex.testingset_pop.size());
    stats_util.current_program__validation__total_score = 0;
    stats_util.current_program__validation__total_passes = 0;
    stats_util.current_program__validation__is_solution = false;
    prob_utils_ForLoopIndex.population_validation_outputs[stats_util.cur_progID].resize(prob_utils_ForLoopIndex.testingset_pop.size());
    // For each test in validation set, evaluate program.
    for (size_t testID = 0; testID < prob_utils_ForLoopIndex.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_ForLoopIndex.testingset_pop[testID];
      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      stats_util.current_program__validation__test_results[testID] = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      stats_util.current_program__validation__total_score += stats_util.current_program__validation__test_results[testID].score;
      stats_util.current_program__validation__total_passes += (size_t)stats_util.current_program__validation__test_results[testID].pass;
      prob_utils_ForLoopIndex.population_validation_outputs[stats_util.cur_progID][testID] = prob_utils_ForLoopIndex.submitted_vec;
    }
    stats_util.current_program__validation__is_solution = stats_util.current_program__validation__total_passes == prob_utils_ForLoopIndex.testingset_pop.size();
    end_program_eval.Trigger(prog_org);
  };
  program_stats.get_prog_behavioral_diversity = [this]() { return emp::ShannonEntropy(prob_utils_ForLoopIndex.population_validation_outputs); };
  program_stats.get_prog_unique_behavioral_phenotypes = [this]() { return emp::UniqueCount(prob_utils_ForLoopIndex.population_validation_outputs); };

  // How should we screen for a solution?
  ScreenForSolution = [this](prog_org_t & prog_org) {
    begin_program_eval.Trigger(prog_org);
    for (size_t testID = 0; testID < prob_utils_ForLoopIndex.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_ForLoopIndex.testingset_pop[testID];

      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      
      TestResult result = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      if (!result.pass) {
        end_program_eval.Trigger(prog_org);
        return false;
      }
    }
    end_program_eval.Trigger(prog_org);
    return true;
  };

  // Tell the experiment how to get test phenotypes.
  GetTestPhenotype = [this](size_t testID) -> test_org_phen_t & {
    emp_assert(prob_ForLoopIndex_world->IsOccupied(testID));
    return prob_ForLoopIndex_world->GetOrg(testID).GetPhenotype();
  };

  // Setup how test world updates.
  SetupTestCaseWorldUpdate(prob_ForLoopIndex_world);

  // Setup how test cases mutate.
  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::RANDOM) {
    std::cout << "RANDOM training mode detected, configuring mutation function to RANDOMIZE organisms." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_ForLoopIndex.MIN_START_END = PROB_FOR_LOOP_INDEX__START_END_MIN;
      prob_utils_ForLoopIndex.MAX_START_END = PROB_FOR_LOOP_INDEX__START_END_MAX;
      prob_utils_ForLoopIndex.MIN_STEP = PROB_FOR_LOOP_INDEX__STEP_MIN;
      prob_utils_ForLoopIndex.MAX_STEP = PROB_FOR_LOOP_INDEX__STEP_MAX;
      prob_utils_ForLoopIndex.MUT_RATE = 1;
      // (2) Hook mutator up to world.
      prob_ForLoopIndex_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_ForLoopIndex.Mutate(rnd, test_org.GetGenome());
      });
    };
  } else {
    std::cout << "Non-RANDOM training mode detected, configuring mutation function normally." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_ForLoopIndex.MIN_START_END = PROB_FOR_LOOP_INDEX__START_END_MIN;
      prob_utils_ForLoopIndex.MAX_START_END = PROB_FOR_LOOP_INDEX__START_END_MAX;
      prob_utils_ForLoopIndex.MIN_STEP = PROB_FOR_LOOP_INDEX__STEP_MIN;
      prob_utils_ForLoopIndex.MAX_STEP = PROB_FOR_LOOP_INDEX__STEP_MAX;
      prob_utils_ForLoopIndex.MUT_RATE = PROB_FOR_LOOP_INDEX__MUTATION__MUT_RATE;
      // (2) Hook mutator up to world.
      prob_ForLoopIndex_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_ForLoopIndex.Mutate(rnd, test_org.GetGenome());
      });
    };
  }

  // Setup test case fitness function.
  SetupTestFitFun = [this]() {
    prob_ForLoopIndex_world->SetFitFun([](test_org_t & test_org) {
      return (double)test_org.GetPhenotype().num_fails;
    });
  };

  // Tell experiment how to configure hardware inputs when running program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org, emp::Ptr<TestOrg_Base> test_org_base_ptr) {
    // Reset eval stuff
    // Set current test org.
    prob_utils_ForLoopIndex.cur_eval_test_org = test_org_base_ptr.Cast<test_org_t>(); // currently only place need testID for this?
    prob_utils_ForLoopIndex.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 1);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      // Grab some useful references.
      Problem_ForLoopIndex_input_t & input = prob_utils_ForLoopIndex.cur_eval_test_org->GetGenome(); // std::pair<int, double>
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      // Set hardware input.
      wmem.Set(0, input[0]);
      wmem.Set(1, input[1]);
      wmem.Set(2, input[2]);
    }
  });

  // Tell experiment how to snapshot test population.
  SnapshotTests = [this]() {
    std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate());
    mkdir(snapshot_dir.c_str(), ACCESSPERMS);
    
    emp::DataFile file(snapshot_dir + "/test_pop_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");
    // Test file contents:
    // - test id
    std::function<size_t(void)> get_test_id = [this]() { return stats_util.cur_testID; };
    file.AddFun(get_test_id, "test_id");

    // - test fitness
    std::function<double(void)> get_test_fitness = [this]() { return prob_ForLoopIndex_world->CalcFitnessID(stats_util.cur_testID); };
    file.AddFun(get_test_fitness, "fitness");

    // - num passes
    std::function<size_t(void)> get_test_num_passes = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_passes; };
    file.AddFun(get_test_num_passes, "num_passes");

    // - num fails
    std::function<size_t(void)> get_test_num_fails = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_fails; };
    file.AddFun(get_test_num_fails, "num_fails");

    std::function<size_t(void)> get_num_tested = [this]() { return GetTestPhenotype(stats_util.cur_testID).test_passes.size(); };
    file.AddFun(get_num_tested, "num_programs_tested_against");

    // - test scores by program
    std::function<std::string(void)> get_passes_by_program = [this]() {
      std::string scores = "\"[";
      test_org_phen_t & phen = GetTestPhenotype(stats_util.cur_testID);
      for (size_t i = 0; i < phen.test_passes.size(); ++i) {
        if (i) scores += ",";
        scores += emp::to_string(phen.test_passes[i]);
      }
      scores += "]\"";
      return scores;
    };
    file.AddFun(get_passes_by_program, "passes_by_program");

    // - test
    std::function<std::string(void)> get_test = [this]() {
      std::ostringstream stream;
      stream << "\"";
      prob_ForLoopIndex_world->GetOrg(stats_util.cur_testID).Print(stream);
      stream << "\"";
      return stream.str();
    };
    file.AddFun(get_test, "test", "");

    file.PrintHeaderKeys();

    // Loop over tests, snapshotting each.
    for (stats_util.cur_testID = 0; stats_util.cur_testID < prob_ForLoopIndex_world->GetSize(); ++stats_util.cur_testID) {
      if (!prob_ForLoopIndex_world->IsOccupied(stats_util.cur_testID)) continue;
      file.Update();
    }
  };

  // Add default instructions to instruction set.
  AddDefaultInstructions({"Add",
                          "Sub",
                          "Mult",
                          "Div",
                          "Mod",
                          "TestNumEqu",
                          "TestNumNEqu",
                          "TestNumLess",
                          "TestNumLessTEqu",
                          "TestNumGreater",
                          "TestNumGreaterTEqu",
                          "Floor",
                          "Not",
                          "Inc",
                          "Dec",
                          "CopyMem",
                          "SwapMem",
                          "Input",
                          "Output",
                          "CommitGlobal",
                          "PullGlobal",
                          "TestMemEqu",
                          "TestMemNEqu",
                          "If",
                          "IfNot",
                          "While",
                          "Countdown",
                          "Foreach",
                          "Close",
                          "Break",
                          "Call",
                          "Routine",
                          "Return",
                          "ModuleDef"
  });

  // -- Custom instructions --
  inst_lib->AddInst("LoadStart", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadStart_ForLoopIndex(hw, inst);
  }, 1);

  inst_lib->AddInst("LoadStep", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadStep_ForLoopIndex(hw, inst);
  }, 1);

  inst_lib->AddInst("LoadEnd", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadEnd_ForLoopIndex(hw, inst);
  }, 1);

  inst_lib->AddInst("SubmitNum", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitNum_ForLoopIndex(hw, inst);
  }, 0);
}

void ProgramSynthesisExperiment::SetupProblem_CompareStringLengths() { 
  std::cout << "Setting up problem - CompareStringLengths." << std::endl;

  // A few useful aliases.
  using test_org_t = TestOrg_CompareStringLengths;

  // Load benchmark data for problem.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  std::cout << "Loading training examples." << std::endl;
  prob_utils_CompareStringLengths.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  std::cout << "Loading testing examples." << std::endl;
  prob_utils_CompareStringLengths.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  std::cout << "Generating testing set population." << std::endl;
  prob_utils_CompareStringLengths.GenerateTestingSetPop();
  std::cout << "Loaded training example set size = " << prob_utils_CompareStringLengths.GetTrainingSet().GetSize() << std::endl;
  std::cout << "Loaded testing example set size = " << prob_utils_CompareStringLengths.GetTestingSet().GetSize() << std::endl;
  std::cout << "Testing set (non-training examples used to evaluate program accuracy) size = " << prob_utils_CompareStringLengths.testingset_pop.size() << std::endl;

  // Setup the world
  NewTestCaseWorld(prob_CompareStringLengths_world, *random, "CompareStringLengths world");

    // Configure how the population should be initialized
  SetupTestCasePop_Init(prob_CompareStringLengths_world,
                        prob_utils_CompareStringLengths.training_set,
                        [this]() { return GenRandomTestInput_CompareStringLengths(*random, 
                                                                                  {PROB_COMPARE_STRING_LENGTHS__MIN_STR_LEN, PROB_COMPARE_STRING_LENGTHS__MAX_STR_LEN}); 
                                  }
                        );
  end_setup_sig.AddAction([this]() { std::cout << "TestCase world size= " << prob_CompareStringLengths_world->GetSize() << std::endl; });

  // Tell the world to calculate the correct test output (given input) on placement.
  prob_CompareStringLengths_world->OnPlacement([this](size_t pos) { prob_CompareStringLengths_world->GetOrg(pos).CalcOut(); });

  // How are program results calculated on a test?
  CalcProgramResultOnTest = [this](prog_org_t & prog_org, TestOrg_Base & test_org_base) {
    test_org_t & test_org = static_cast<test_org_t&>(test_org_base);
    TestResult result;
    if (!prob_utils_CompareStringLengths.submitted) {
      result.score = 0;
      result.pass = false;
      result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_CompareStringLengths.CalcScorePassFail(test_org.GetCorrectOut(), prob_utils_CompareStringLengths.submitted_val));
      result.score = r.first;
      result.pass = r.second;
      result.sub = true;
    }
    return result;
  };

  // Setup how evaluation on world test should work.
  EvaluateWorldTest = [this](prog_org_t & prog_org, size_t testID) {
    emp::Ptr<test_org_t> test_org_ptr = prob_CompareStringLengths_world->GetOrgPtr(testID);
    begin_program_test.Trigger(prog_org, test_org_ptr);
    do_program_test.Trigger(prog_org, test_org_ptr);
    end_program_test.Trigger(prog_org, test_org_ptr);
    return CalcProgramResultOnTest(prog_org, *test_org_ptr);
  };

  // How should we validate programs on testing set?
  prob_utils_CompareStringLengths.population_validation_outputs.resize(PROG_POP_SIZE);
  DoTestingSetValidation = [this](prog_org_t & prog_org) { 
    // evaluate program on full testing set; update stats utils with results
    begin_program_eval.Trigger(prog_org);
    stats_util.current_program__validation__test_results.resize(prob_utils_CompareStringLengths.testingset_pop.size());
    stats_util.current_program__validation__total_score = 0;
    stats_util.current_program__validation__total_passes = 0;
    stats_util.current_program__validation__is_solution = false;
    prob_utils_CompareStringLengths.population_validation_outputs[stats_util.cur_progID].resize(prob_utils_CompareStringLengths.testingset_pop.size());
    // For each test in validation set, evaluate program.
    for (size_t testID = 0; testID < prob_utils_CompareStringLengths.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_CompareStringLengths.testingset_pop[testID];
      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      stats_util.current_program__validation__test_results[testID] = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      stats_util.current_program__validation__total_score += stats_util.current_program__validation__test_results[testID].score;
      stats_util.current_program__validation__total_passes += (size_t)stats_util.current_program__validation__test_results[testID].pass;
      prob_utils_CompareStringLengths.population_validation_outputs[stats_util.cur_progID][testID] = prob_utils_CompareStringLengths.submitted_val;
    }
    stats_util.current_program__validation__is_solution = stats_util.current_program__validation__total_passes == prob_utils_CompareStringLengths.testingset_pop.size();
    end_program_eval.Trigger(prog_org);
  };
  program_stats.get_prog_behavioral_diversity = [this]() { return emp::ShannonEntropy(prob_utils_CompareStringLengths.population_validation_outputs); };
  program_stats.get_prog_unique_behavioral_phenotypes = [this]() { return emp::UniqueCount(prob_utils_CompareStringLengths.population_validation_outputs); };

  // How should we screen for a solution?
  ScreenForSolution = [this](prog_org_t & prog_org) {
    begin_program_eval.Trigger(prog_org);
    for (size_t testID = 0; testID < prob_utils_CompareStringLengths.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_CompareStringLengths.testingset_pop[testID];

      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      
      TestResult result = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      if (!result.pass) {
        end_program_eval.Trigger(prog_org);
        return false;
      }
    }
    end_program_eval.Trigger(prog_org);
    return true;
  };

  // Tell the experiment how to get test phenotypes.
  GetTestPhenotype = [this](size_t testID) -> test_org_phen_t & {
    emp_assert(prob_CompareStringLengths_world->IsOccupied(testID));
    return prob_CompareStringLengths_world->GetOrg(testID).GetPhenotype();
  };

  // Setup how test world updates.
  SetupTestCaseWorldUpdate(prob_CompareStringLengths_world);

  // Setup how test cases mutate.
  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::RANDOM) {
    std::cout << "RANDOM training mode detected, configuring mutation function to RANDOMIZE organisms." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Randomize organism genome on mutate.
      prob_CompareStringLengths_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        test_org.GetGenome() = GenRandomTestInput_CompareStringLengths(*random, {PROB_COMPARE_STRING_LENGTHS__MIN_STR_LEN, PROB_COMPARE_STRING_LENGTHS__MAX_STR_LEN}); 
        return 1;
      });
    };
  } else {
    std::cout << "Non-RANDOM training mode detected, configuring mutation function normally." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_CompareStringLengths.MIN_STR_LEN = PROB_COMPARE_STRING_LENGTHS__MIN_STR_LEN;
      prob_utils_CompareStringLengths.MAX_STR_LEN = PROB_COMPARE_STRING_LENGTHS__MAX_STR_LEN;
      prob_utils_CompareStringLengths.PER_SITE_INS_RATE = PROB_COMPARE_STRING_LENGTHS__PER_SITE_INS_RATE;
      prob_utils_CompareStringLengths.PER_SITE_DEL_RATE = PROB_COMPARE_STRING_LENGTHS__PER_SITE_DEL_RATE;
      prob_utils_CompareStringLengths.PER_SITE_SUB_RATE = PROB_COMPARE_STRING_LENGTHS__PER_SITE_SUB_RATE;
      prob_utils_CompareStringLengths.PER_STR_SWAP_RATE = PROB_COMPARE_STRING_LENGTHS__PER_STR_SWAP_RATE;
      // (2) Hook mutator up to world.
      prob_CompareStringLengths_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_CompareStringLengths.Mutate(rnd, test_org.GetGenome());
      });
    };
  } // todo - test mutation function

  // Setup test case fitness function.
  SetupTestFitFun = [this]() {
    prob_CompareStringLengths_world->SetFitFun([](test_org_t & test_org) {
      return (double)test_org.GetPhenotype().num_fails;
    });
  };

  // Tell experiment how to configure hardware inputs when running program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org, emp::Ptr<TestOrg_Base> test_org_base_ptr) {
    // Reset eval stuff
    // Set current test org.
    prob_utils_CompareStringLengths.cur_eval_test_org = test_org_base_ptr.Cast<test_org_t>(); // currently only place need testID for this?
    prob_utils_CompareStringLengths.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 1);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      // Grab some useful references.
      Problem_CompareStringLengths_input_t & input = prob_utils_CompareStringLengths.cur_eval_test_org->GetGenome(); // std::pair<int, double>
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      // Set hardware input.
      wmem.Set(0, input[0]);
      wmem.Set(1, input[1]);
      wmem.Set(2, input[2]);
    }
  });

  // Tell experiment how to snapshot test population.
  SnapshotTests = [this]() {
    std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate());
    mkdir(snapshot_dir.c_str(), ACCESSPERMS);
    
    emp::DataFile file(snapshot_dir + "/test_pop_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");
    // Test file contents:
    // - test id
    std::function<size_t(void)> get_test_id = [this]() { return stats_util.cur_testID; };
    file.AddFun(get_test_id, "test_id");

    // - test fitness
    std::function<double(void)> get_test_fitness = [this]() { return prob_CompareStringLengths_world->CalcFitnessID(stats_util.cur_testID); };
    file.AddFun(get_test_fitness, "fitness");

    // - num passes
    std::function<size_t(void)> get_test_num_passes = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_passes; };
    file.AddFun(get_test_num_passes, "num_passes");

    // - num fails
    std::function<size_t(void)> get_test_num_fails = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_fails; };
    file.AddFun(get_test_num_fails, "num_fails");

    std::function<size_t(void)> get_num_tested = [this]() { return GetTestPhenotype(stats_util.cur_testID).test_passes.size(); };
    file.AddFun(get_num_tested, "num_programs_tested_against");

    // - test scores by program
    std::function<std::string(void)> get_passes_by_program = [this]() {
      std::string scores = "\"[";
      test_org_phen_t & phen = GetTestPhenotype(stats_util.cur_testID);
      for (size_t i = 0; i < phen.test_passes.size(); ++i) {
        if (i) scores += ",";
        scores += emp::to_string(phen.test_passes[i]);
      }
      scores += "]\"";
      return scores;
    };
    file.AddFun(get_passes_by_program, "passes_by_program");

    // - test
    std::function<std::string(void)> get_test = [this]() {
      std::ostringstream stream;
      stream << "\"";
      prob_CompareStringLengths_world->GetOrg(stats_util.cur_testID).Print(stream);
      stream << "\"";
      return stream.str();
    };
    file.AddFun(get_test, "test", "");

    file.PrintHeaderKeys();

    // Loop over tests, snapshotting each.
    for (stats_util.cur_testID = 0; stats_util.cur_testID < prob_CompareStringLengths_world->GetSize(); ++stats_util.cur_testID) {
      if (!prob_CompareStringLengths_world->IsOccupied(stats_util.cur_testID)) continue;
      file.Update();
    }
  };

  // Add default instructions to instruction set.
  AddDefaultInstructions({"Add",
                          "Sub",
                          "Mult",
                          "Div",
                          "Mod",
                          "TestNumEqu",
                          "TestNumNEqu",
                          "TestNumLess",
                          "TestNumLessTEqu",
                          "TestNumGreater",
                          "TestNumGreaterTEqu",
                          "Floor",
                          "Not",
                          "Inc",
                          "Dec",
                          "CopyMem",
                          "SwapMem",
                          "Input",
                          "Output",
                          "CommitGlobal",
                          "PullGlobal",
                          "TestMemEqu",
                          "TestMemNEqu",
                          "If",
                          "IfNot",
                          "While",
                          "Countdown",
                          "Foreach",
                          "Close",
                          "Break",
                          "Call",
                          "Routine",
                          "Return",
                          "ModuleDef",
                          "IsNum",
                          "IsStr",
                          "StrLength",
                          "StrConcat"
  });

  // -- Custom instructions --
  inst_lib->AddInst("LoadStr1", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadStr1_CompareStringLengths(hw, inst);
  }, 1);

  inst_lib->AddInst("LoadStr2", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadStr2_CompareStringLengths(hw, inst);
  }, 1);

  inst_lib->AddInst("LoadStr3", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadStr3_CompareStringLengths(hw, inst);
  }, 1);

  inst_lib->AddInst("SubmitTrue", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitTrue_CompareStringLengths(hw, inst);
  }, 1);

  inst_lib->AddInst("SubmitFalse", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitFalse_CompareStringLengths(hw, inst);
  }, 1);

  inst_lib->AddInst("SubmitVal", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitVal_CompareStringLengths(hw, inst);
  }, 1);

}

void ProgramSynthesisExperiment::SetupProblem_DoubleLetters() { 
  std::cout << "Problem setup not yet implemented... Exiting." << std::endl;
  exit(-1); 
}

void ProgramSynthesisExperiment::SetupProblem_CollatzNumbers() { 
  std::cout << "Setting up problem - Collatz Numbers" << std::endl;

  // A few useful aliases
  using test_org_t = TestOrg_CollatzNumbers;

  // Load benchmark data for problem.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  std::cout << "Loading training examples." << std::endl;
  prob_utils_CollatzNumbers.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  std::cout << "Loading testing examples." << std::endl;
  prob_utils_CollatzNumbers.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  std::cout << "Generating testing set population." << std::endl;
  prob_utils_CollatzNumbers.GenerateTestingSetPop();
  std::cout << "Loaded training example set size = " << prob_utils_CollatzNumbers.GetTrainingSet().GetSize() << std::endl;
  std::cout << "Loaded testing example set size = " << prob_utils_CollatzNumbers.GetTestingSet().GetSize() << std::endl;
  std::cout << "Testing set (non-training examples used to evaluate program accuracy) size = " << prob_utils_CollatzNumbers.testingset_pop.size() << std::endl;
  
  prob_utils_CollatzNumbers.MAX_ERROR = 256; // Lazy, lazy. Problem-specific, magic number for max error...

  // Setup the world
  NewTestCaseWorld(prob_CollatzNumbers_world, *random, "Median world");

  // Configure how the population should be initialized
  SetupTestCasePop_Init(prob_CollatzNumbers_world,
                        prob_utils_CollatzNumbers.training_set,
                        [this]() { return GenRandomTestInput_CollatzNumbers(*random, 
                                                                           {PROB_COLLATZ_NUMBERS__MIN_NUM, PROB_COLLATZ_NUMBERS__MAX_NUM}); 
                                  }
                        );
  end_setup_sig.AddAction([this]() { std::cout << "TestCase world size= " << prob_CollatzNumbers_world->GetSize() << std::endl; });

  // Tell the world to calculate the correct test output (given input) on placement.
  prob_CollatzNumbers_world->OnPlacement([this](size_t pos) { 
    prob_CollatzNumbers_world->GetOrg(pos).SetCache(prob_utils_CollatzNumbers.out_cache);
    prob_CollatzNumbers_world->GetOrg(pos).CalcOut(); 
  });

  // How are program results calculated on a test?
  CalcProgramResultOnTest = [this](prog_org_t & prog_org, TestOrg_Base & test_org_base) {
    test_org_t & test_org = static_cast<test_org_t&>(test_org_base);
    TestResult result;
    if (!prob_utils_CollatzNumbers.submitted) {
      result.score = 0;
      result.pass = false;
      result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_CollatzNumbers.CalcScoreGradient(test_org.GetCorrectOut(), prob_utils_CollatzNumbers.submitted_val));
      result.score = r.first;
      result.pass = r.second;
      result.sub = true;
    }
    return result;
  };

  // Setup how evaluation on world test should work.
  EvaluateWorldTest = [this](prog_org_t & prog_org, size_t testID) {
    emp::Ptr<test_org_t> test_org_ptr = prob_CollatzNumbers_world->GetOrgPtr(testID);
    begin_program_test.Trigger(prog_org, test_org_ptr);
    do_program_test.Trigger(prog_org, test_org_ptr);
    end_program_test.Trigger(prog_org, test_org_ptr);
    return CalcProgramResultOnTest(prog_org, *test_org_ptr);
  };

  // How should we validate programs on testing set?
  prob_utils_CollatzNumbers.population_validation_outputs.resize(PROG_POP_SIZE);
  DoTestingSetValidation = [this](prog_org_t & prog_org) { 
    // evaluate program on full testing set; update stats utils with results
    begin_program_eval.Trigger(prog_org);
    stats_util.current_program__validation__test_results.resize(prob_utils_CollatzNumbers.testingset_pop.size());
    stats_util.current_program__validation__total_score = 0;
    stats_util.current_program__validation__total_passes = 0;
    stats_util.current_program__validation__is_solution = false;
    prob_utils_CollatzNumbers.population_validation_outputs[stats_util.cur_progID].resize(prob_utils_CollatzNumbers.testingset_pop.size());
    // For each test in validation set, evaluate program.
    for (size_t testID = 0; testID < prob_utils_CollatzNumbers.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_CollatzNumbers.testingset_pop[testID];
      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      stats_util.current_program__validation__test_results[testID] = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      stats_util.current_program__validation__total_score += stats_util.current_program__validation__test_results[testID].score;
      stats_util.current_program__validation__total_passes += (size_t)stats_util.current_program__validation__test_results[testID].pass;
      prob_utils_CollatzNumbers.population_validation_outputs[stats_util.cur_progID][testID] = prob_utils_CollatzNumbers.submitted_val; 
    }
    stats_util.current_program__validation__is_solution = stats_util.current_program__validation__total_passes == prob_utils_CollatzNumbers.testingset_pop.size();
    end_program_eval.Trigger(prog_org);
  };
  program_stats.get_prog_behavioral_diversity = [this]() { return emp::ShannonEntropy(prob_utils_CollatzNumbers.population_validation_outputs); };
  program_stats.get_prog_unique_behavioral_phenotypes = [this]() { return emp::UniqueCount(prob_utils_CollatzNumbers.population_validation_outputs); };
  
  // How should we screen for a solution?
  ScreenForSolution = [this](prog_org_t & prog_org) {
    begin_program_eval.Trigger(prog_org);
    for (size_t testID = 0; testID < prob_utils_CollatzNumbers.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_CollatzNumbers.testingset_pop[testID];

      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      
      TestResult result = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      if (!result.pass) {
        end_program_eval.Trigger(prog_org);
        return false;
      }
    }
    end_program_eval.Trigger(prog_org);
    return true;
  };

  // Tell the experiment how to get test phenotypes.
  GetTestPhenotype = [this](size_t testID) -> test_org_phen_t & {
    emp_assert(prob_CollatzNumbers_world->IsOccupied(testID));
    return prob_CollatzNumbers_world->GetOrg(testID).GetPhenotype();
  };

  // Setup how test world updates.
  SetupTestCaseWorldUpdate(prob_CollatzNumbers_world);

  // Setup how test cases mutate.
  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::RANDOM) {
    std::cout << "RANDOM training mode detected, configuring mutation function to RANDOMIZE organisms." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Randomize organism genome on mutate.
      prob_CollatzNumbers_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        test_org.GetGenome() = GenRandomTestInput_CollatzNumbers(*random, {PROB_COLLATZ_NUMBERS__MIN_NUM, PROB_COLLATZ_NUMBERS__MAX_NUM}); 
        return 1;
      });
    };
  } else {
    std::cout << "Non-RANDOM training mode detected, configuring mutation function normally." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_CollatzNumbers.MIN_NUM = PROB_COLLATZ_NUMBERS__MIN_NUM;
      prob_utils_CollatzNumbers.MAX_NUM = PROB_COLLATZ_NUMBERS__MAX_NUM;
      prob_utils_CollatzNumbers.NUM_SUB_RATE = PROB_COLLATZ_NUMBERS__MUTATION__PER_NUM_SUB_RATE;
      // (2) Hook mutator up to world.
      prob_CollatzNumbers_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_CollatzNumbers.Mutate(rnd, test_org.GetGenome());
      });
    };
  }

  // Setup test case fitness function.
  SetupTestFitFun = [this]() {
    prob_CollatzNumbers_world->SetFitFun([](test_org_t & test_org) {
      return (double)test_org.GetPhenotype().num_fails;
    });
  };

  // Tell experiment how to configure hardware inputs when running program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org, emp::Ptr<TestOrg_Base> test_org_base_ptr) {
    // Reset eval stuff
    // Set current test org.
    prob_utils_CollatzNumbers.cur_eval_test_org = test_org_base_ptr.Cast<test_org_t>(); // currently only place need testID for this?
    prob_utils_CollatzNumbers.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 3);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      // Grab some useful references.
      Problem_CollatzNumbers_input_t & input = prob_utils_CollatzNumbers.cur_eval_test_org->GetGenome(); // std::pair<int, double>
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      // Set hardware input.
      wmem.Set(0, input);
    }
  });

  // Tell experiment how to snapshot test population.
  SnapshotTests = [this]() {
    std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate());
    mkdir(snapshot_dir.c_str(), ACCESSPERMS);
    
    emp::DataFile file(snapshot_dir + "/test_pop_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");
    // Test file contents:
    // - test id
    std::function<size_t(void)> get_test_id = [this]() { return stats_util.cur_testID; };
    file.AddFun(get_test_id, "test_id");

    // - test fitness
    std::function<double(void)> get_test_fitness = [this]() { return prob_CollatzNumbers_world->CalcFitnessID(stats_util.cur_testID); };
    file.AddFun(get_test_fitness, "fitness");

    // - num passes
    std::function<size_t(void)> get_test_num_passes = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_passes; };
    file.AddFun(get_test_num_passes, "num_passes");

    // - num fails
    std::function<size_t(void)> get_test_num_fails = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_fails; };
    file.AddFun(get_test_num_fails, "num_fails");

    std::function<size_t(void)> get_num_tested = [this]() { return GetTestPhenotype(stats_util.cur_testID).test_passes.size(); };
    file.AddFun(get_num_tested, "num_programs_tested_against");

    // - test scores by program
    std::function<std::string(void)> get_passes_by_program = [this]() {
      std::string scores = "\"[";
      test_org_phen_t & phen = GetTestPhenotype(stats_util.cur_testID);
      for (size_t i = 0; i < phen.test_passes.size(); ++i) {
        if (i) scores += ",";
        scores += emp::to_string(phen.test_passes[i]);
      }
      scores += "]\"";
      return scores;
    };
    file.AddFun(get_passes_by_program, "passes_by_program");

    // - test
    std::function<std::string(void)> get_test = [this]() {
      std::ostringstream stream;
      stream << "\"";
      prob_CollatzNumbers_world->GetOrg(stats_util.cur_testID).Print(stream);
      stream << "\"";
      return stream.str();
    };
    file.AddFun(get_test, "test");

    file.PrintHeaderKeys();

    // Loop over tests, snapshotting each.
    for (stats_util.cur_testID = 0; stats_util.cur_testID < prob_CollatzNumbers_world->GetSize(); ++stats_util.cur_testID) {
      if (!prob_CollatzNumbers_world->IsOccupied(stats_util.cur_testID)) continue;
      file.Update();
    }
  };

  // Add default instructions to instruction set.
  AddDefaultInstructions({"Add",
                          "Sub",
                          "Mult",
                          "Div",
                          "Mod",
                          "TestNumEqu",
                          "TestNumNEqu",
                          "TestNumLess",
                          "TestNumLessTEqu",
                          "TestNumGreater",
                          "TestNumGreaterTEqu",
                          "Floor",
                          "Not",
                          "Inc",
                          "Dec",
                          "CopyMem",
                          "SwapMem",
                          "Input",
                          "Output",
                          "CommitGlobal",
                          "PullGlobal",
                          "TestMemEqu",
                          "TestMemNEqu",
                          "If",
                          "IfNot",
                          "While",
                          "Countdown",
                          "Foreach",
                          "Close",
                          "Break",
                          "Call",
                          "Routine",
                          "Return",
                          "ModuleDef",
                          "MakeVector",
                          "VecGet",
                          "VecSet",
                          "VecLen",
                          "VecAppend",
                          "VecPop",
                          "VecRemove",
                          "VecReplaceAll",
                          "VecIndexOf",
                          "VecOccurrencesOf",
                          "VecReverse",
                          "VecSwapIfLess",
                          "VecGetFront",
                          "VecGetBack",
                          "IsNum",
                          "IsVec"
  });

  // -- Custom Instructions --

  for (size_t i = 0; i <= 10; ++i) {
    inst_lib->AddInst("Set-" + emp::to_string(i),
      [i](hardware_t & hw, const inst_t & inst) {
        hardware_t::CallState & state = hw.GetCurCallState();
        hardware_t::Memory & wmem = state.GetWorkingMem();
        size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
        if (!hw.IsValidMemPos(posA)) return; // Do nothing
        wmem.Set(posA, (double)i);
      });
  }

  inst_lib->AddInst("LoadNum", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadNum_CollatzNumbers(hw, inst);
  }, 1);

  inst_lib->AddInst("SubmitNum", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitNum_CollatzNumbers(hw, inst);
  }, 1);

}

void ProgramSynthesisExperiment::SetupProblem_ReplaceSpaceWithNewline() { 
  std::cout << "Problem setup not yet implemented... Exiting." << std::endl;
  exit(-1); 
}

void ProgramSynthesisExperiment::SetupProblem_StringDifferences() { 
  std::cout << "Problem setup not yet implemented... Exiting." << std::endl;
  exit(-1); 
}

void ProgramSynthesisExperiment::SetupProblem_EvenSquares() { 
  std::cout << "Problem setup not yet implemented... Exiting." << std::endl;
  exit(-1); 
}

void ProgramSynthesisExperiment::SetupProblem_WallisPi() { 
  std::cout << "Problem setup not yet implemented... Exiting." << std::endl;
  exit(-1); 
}

void ProgramSynthesisExperiment::SetupProblem_StringLengthsBackwards() { 
  std::cout << "Setting up problem - StringLengthsBackwards" << std::endl;

  // A few useful aliases.
  using test_org_t = TestOrg_StringLengthsBackwards;

  // Load benchmark data for problem.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  std::cout << "Loading training examples." << std::endl;
  prob_utils_StringLengthsBackwards.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  std::cout << "Loading testing examples." << std::endl;
  prob_utils_StringLengthsBackwards.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  std::cout << "Generating testing set population." << std::endl;
  prob_utils_StringLengthsBackwards.GenerateTestingSetPop();
  std::cout << "Loaded training example set size = " << prob_utils_StringLengthsBackwards.GetTrainingSet().GetSize() << std::endl;
  std::cout << "Loaded testing example set size = " << prob_utils_StringLengthsBackwards.GetTestingSet().GetSize() << std::endl;
  std::cout << "Testing set (non-training examples used to evaluate program accuracy) size = " << prob_utils_StringLengthsBackwards.testingset_pop.size() << std::endl;

  // Setup the world
  NewTestCaseWorld(prob_StringLengthsBackwards_world, *random, "StringLengthsBackwards world");

  // Configure how the population should be initialized
  SetupTestCasePop_Init(prob_StringLengthsBackwards_world,
                        prob_utils_StringLengthsBackwards.training_set,
                        [this]() { return GenRandomTestInput_StringLengthsBackwards(*random, 
                                                                                   {PROB_STRING_LENGTHS_BACKWARDS__MIN_STR_CNT, PROB_STRING_LENGTHS_BACKWARDS__MAX_STR_CNT},
                                                                                   {PROB_STRING_LENGTHS_BACKWARDS__MIN_STR_LEN, PROB_STRING_LENGTHS_BACKWARDS__MAX_STR_LEN}); 
                                  }
                        );
  end_setup_sig.AddAction([this]() { std::cout << "TestCase world size= " << prob_StringLengthsBackwards_world->GetSize() << std::endl; });

  // Tell the world to calculate the correct test output (given input) on placement.
  prob_StringLengthsBackwards_world->OnPlacement([this](size_t pos) { prob_StringLengthsBackwards_world->GetOrg(pos).CalcOut(); });

  // How are program results calculated on a test?
  CalcProgramResultOnTest = [this](prog_org_t & prog_org, TestOrg_Base & test_org_base) {
    test_org_t & test_org = static_cast<test_org_t&>(test_org_base);
    TestResult result;
    // if (!prob_utils_StringLengthsBackwards.submitted) {
      // result.score = 0;
      // result.pass = false;
      // result.sub = false;
    // } else {
    std::pair<double, bool> r(prob_utils_StringLengthsBackwards.CalcScoreGradient(test_org.GetCorrectOut(), prob_utils_StringLengthsBackwards.submitted_vec));
    result.score = r.first;
    result.pass = r.second;
    result.sub = true;
    // }
    return result;
  };

  // Setup how evaluation on world test should work.
  EvaluateWorldTest = [this](prog_org_t & prog_org, size_t testID) {
    emp::Ptr<test_org_t> test_org_ptr = prob_StringLengthsBackwards_world->GetOrgPtr(testID);
    begin_program_test.Trigger(prog_org, test_org_ptr);
    do_program_test.Trigger(prog_org, test_org_ptr);
    end_program_test.Trigger(prog_org, test_org_ptr);
    return CalcProgramResultOnTest(prog_org, *test_org_ptr);
  };

  // How should we validate programs on testing set?
  prob_utils_StringLengthsBackwards.population_validation_outputs.resize(PROG_POP_SIZE);
  DoTestingSetValidation = [this](prog_org_t & prog_org) { 
    // evaluate program on full testing set; update stats utils with results
    begin_program_eval.Trigger(prog_org);
    stats_util.current_program__validation__test_results.resize(prob_utils_StringLengthsBackwards.testingset_pop.size());
    stats_util.current_program__validation__total_score = 0;
    stats_util.current_program__validation__total_passes = 0;
    stats_util.current_program__validation__is_solution = false;
    prob_utils_StringLengthsBackwards.population_validation_outputs[stats_util.cur_progID].resize(prob_utils_StringLengthsBackwards.testingset_pop.size());
    // For each test in validation set, evaluate program.
    for (size_t testID = 0; testID < prob_utils_StringLengthsBackwards.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_StringLengthsBackwards.testingset_pop[testID];
      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      stats_util.current_program__validation__test_results[testID] = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      stats_util.current_program__validation__total_score += stats_util.current_program__validation__test_results[testID].score;
      stats_util.current_program__validation__total_passes += (size_t)stats_util.current_program__validation__test_results[testID].pass;
      prob_utils_StringLengthsBackwards.population_validation_outputs[stats_util.cur_progID][testID] = prob_utils_StringLengthsBackwards.submitted_vec;
    }
    stats_util.current_program__validation__is_solution = stats_util.current_program__validation__total_passes == prob_utils_StringLengthsBackwards.testingset_pop.size();
    end_program_eval.Trigger(prog_org);
  }; 
  program_stats.get_prog_behavioral_diversity = [this]() { return emp::ShannonEntropy(prob_utils_StringLengthsBackwards.population_validation_outputs); };
  program_stats.get_prog_unique_behavioral_phenotypes = [this]() { return emp::UniqueCount(prob_utils_StringLengthsBackwards.population_validation_outputs); };

  // How should we screen for a solution?
  ScreenForSolution = [this](prog_org_t & prog_org) {
    begin_program_eval.Trigger(prog_org);
    for (size_t testID = 0; testID < prob_utils_StringLengthsBackwards.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_StringLengthsBackwards.testingset_pop[testID];

      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      
      TestResult result = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      if (!result.pass) {
        end_program_eval.Trigger(prog_org);
        return false;
      }
    }
    end_program_eval.Trigger(prog_org);
    return true;
  };

  // Tell the experiment how to get test phenotypes.
  GetTestPhenotype = [this](size_t testID) -> test_org_phen_t & {
    emp_assert(prob_StringLengthsBackwards_world->IsOccupied(testID));
    return prob_StringLengthsBackwards_world->GetOrg(testID).GetPhenotype();
  };

  // Setup how test world updates.
  SetupTestCaseWorldUpdate(prob_StringLengthsBackwards_world);

  // Setup how test cases mutate.
  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::RANDOM) {
    std::cout << "RANDOM training mode detected, configuring mutation function to RANDOMIZE organisms." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Randomize organism genome on mutate.
      prob_StringLengthsBackwards_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        test_org.GetGenome() = GenRandomTestInput_StringLengthsBackwards(*random, {PROB_STRING_LENGTHS_BACKWARDS__MIN_STR_CNT, PROB_STRING_LENGTHS_BACKWARDS__MAX_STR_CNT},
                                                                                  {PROB_STRING_LENGTHS_BACKWARDS__MIN_STR_LEN, PROB_STRING_LENGTHS_BACKWARDS__MAX_STR_LEN}); 
        return 1;
      });
    };
  } else {
    std::cout << "Non-RANDOM training mode detected, configuring mutation function normally." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_StringLengthsBackwards.MIN_STR_LEN = PROB_STRING_LENGTHS_BACKWARDS__MIN_STR_LEN;
      prob_utils_StringLengthsBackwards.MAX_STR_LEN = PROB_STRING_LENGTHS_BACKWARDS__MAX_STR_LEN;
      prob_utils_StringLengthsBackwards.MIN_STR_CNT = PROB_STRING_LENGTHS_BACKWARDS__MIN_STR_CNT;
      prob_utils_StringLengthsBackwards.MAX_STR_CNT = PROB_STRING_LENGTHS_BACKWARDS__MAX_STR_CNT;
      prob_utils_StringLengthsBackwards.PER_CHAR_INS_RATE = PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_CHAR_INS_RATE;
      prob_utils_StringLengthsBackwards.PER_CHAR_DEL_RATE = PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_CHAR_DEL_RATE;
      prob_utils_StringLengthsBackwards.PER_CHAR_SUB_RATE = PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_CHAR_SUB_RATE;
      prob_utils_StringLengthsBackwards.PER_STR_SWAP_RATE = PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_STR_SWAP_RATE;
      prob_utils_StringLengthsBackwards.PER_STR_DUP_RATE = PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_STR_DUP_RATE;
      prob_utils_StringLengthsBackwards.PER_STR_DEL_RATE = PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_STR_DEL_RATE;
      // (2) Hook mutator up to world.
      prob_StringLengthsBackwards_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_StringLengthsBackwards.Mutate(rnd, test_org.GetGenome());
      });
    };
  } // todo - test mutation function

  // Setup test case fitness function.
  SetupTestFitFun = [this]() {
    prob_StringLengthsBackwards_world->SetFitFun([](test_org_t & test_org) {
      return (double)test_org.GetPhenotype().num_fails;
    });
  };

  // Tell experiment how to configure hardware inputs when running program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org, emp::Ptr<TestOrg_Base> test_org_base_ptr) {
    // Reset eval stuff
    // Set current test org.
    prob_utils_StringLengthsBackwards.cur_eval_test_org = test_org_base_ptr.Cast<test_org_t>(); // currently only place need testID for this?
    prob_utils_StringLengthsBackwards.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 1);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      // Grab some useful references.
      const Problem_StringLengthsBackwards_input_t & input = prob_utils_StringLengthsBackwards.cur_eval_test_org->GetGenome(); // std::pair<int, double>
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      // Set hardware input.
      wmem.Set(0, input);
    }
  });  

  // Tell experiment how to snapshot test population.
  SnapshotTests = [this]() {
    std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate());
    mkdir(snapshot_dir.c_str(), ACCESSPERMS);
    
    emp::DataFile file(snapshot_dir + "/test_pop_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");
    // Test file contents:
    // - test id
    std::function<size_t(void)> get_test_id = [this]() { return stats_util.cur_testID; };
    file.AddFun(get_test_id, "test_id");

    // - test fitness
    std::function<double(void)> get_test_fitness = [this]() { return prob_StringLengthsBackwards_world->CalcFitnessID(stats_util.cur_testID); };
    file.AddFun(get_test_fitness, "fitness");

    // - num passes
    std::function<size_t(void)> get_test_num_passes = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_passes; };
    file.AddFun(get_test_num_passes, "num_passes");

    // - num fails
    std::function<size_t(void)> get_test_num_fails = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_fails; };
    file.AddFun(get_test_num_fails, "num_fails");

    std::function<size_t(void)> get_num_tested = [this]() { return GetTestPhenotype(stats_util.cur_testID).test_passes.size(); };
    file.AddFun(get_num_tested, "num_programs_tested_against");

    // - test scores by program
    std::function<std::string(void)> get_passes_by_program = [this]() {
      std::string scores = "\"[";
      test_org_phen_t & phen = GetTestPhenotype(stats_util.cur_testID);
      for (size_t i = 0; i < phen.test_passes.size(); ++i) {
        if (i) scores += ",";
        scores += emp::to_string(phen.test_passes[i]);
      }
      scores += "]\"";
      return scores;
    };
    file.AddFun(get_passes_by_program, "passes_by_program");

    // - test
    std::function<std::string(void)> get_test = [this]() {
      std::ostringstream stream;
      stream << "\"";
      prob_StringLengthsBackwards_world->GetOrg(stats_util.cur_testID).Print(stream);
      stream << "\"";
      return stream.str();
    };
    file.AddFun(get_test, "test", "");

    file.PrintHeaderKeys();

    // Loop over tests, snapshotting each.
    for (stats_util.cur_testID = 0; stats_util.cur_testID < prob_StringLengthsBackwards_world->GetSize(); ++stats_util.cur_testID) {
      if (!prob_StringLengthsBackwards_world->IsOccupied(stats_util.cur_testID)) continue;
      file.Update();
    }
  };

  // Add default instructions to instruction set.
  AddDefaultInstructions(); // Add them all!

  // -- Custom instructions --
  inst_lib->AddInst("LoadStrVec", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadStrVec_StringLengthsBackwards(hw, inst);
  }, 1);

  inst_lib->AddInst("SubmitVal", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitVal_StringLengthsBackwards(hw, inst);
  }, 1);

  // inst_lib->AddInst("SubmitVec", [this](hardware_t & hw, const inst_t & inst) {
  //   this->Inst_SubmitVec_StringLengthsBackwards(hw, inst);
  // }, 1);

}

void ProgramSynthesisExperiment::SetupProblem_LastIndexOfZero() { 
  std::cout << "Setting up problem - LastIndexOfZero" << std::endl;

  // A few useful aliases
  using test_org_t = TestOrg_LastIndexOfZero;

  // Load benchmark data for problem.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  std::cout << "Loading training examples." << std::endl;
  prob_utils_LastIndexOfZero.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  std::cout << "Loading testing examples." << std::endl;
  prob_utils_LastIndexOfZero.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  std::cout << "Generating testing set population." << std::endl;
  prob_utils_LastIndexOfZero.GenerateTestingSetPop();
  std::cout << "Loaded training example set size = " << prob_utils_LastIndexOfZero.GetTrainingSet().GetSize() << std::endl;
  std::cout << "Loaded testing example set size = " << prob_utils_LastIndexOfZero.GetTestingSet().GetSize() << std::endl;
  std::cout << "Testing set (non-training examples used to evaluate program accuracy) size = " << prob_utils_LastIndexOfZero.testingset_pop.size() << std::endl;

  prob_utils_LastIndexOfZero.MAX_ERROR = PROB_LAST_INDEX_OF_ZERO__MAX_VEC_LEN;

  // Setup the world
  NewTestCaseWorld(prob_LastIndexOfZero_world, *random, "LastIndexOfZero world");

  // Configure how the population should be initialized
  SetupTestCasePop_Init(prob_LastIndexOfZero_world,
                        prob_utils_LastIndexOfZero.training_set,
                        [this]() { return GenRandomTestInput_LastIndexOfZero(*random, 
                                                                             {PROB_LAST_INDEX_OF_ZERO__MIN_VEC_LEN, PROB_LAST_INDEX_OF_ZERO__MAX_VEC_LEN},
                                                                             {PROB_LAST_INDEX_OF_ZERO__MIN_NUM, PROB_LAST_INDEX_OF_ZERO__MAX_NUM}); 
                                  }
                        );
  end_setup_sig.AddAction([this]() { std::cout << "TestCase world size= " << prob_LastIndexOfZero_world->GetSize() << std::endl; });

  // Tell the world to calculate the correct test output (given input) on placement.
  prob_LastIndexOfZero_world->OnPlacement([this](size_t pos) { prob_LastIndexOfZero_world->GetOrg(pos).CalcOut(); });

  // How are program results calculated on a test?
  CalcProgramResultOnTest = [this](prog_org_t & prog_org, TestOrg_Base & test_org_base) {
    test_org_t & test_org = static_cast<test_org_t&>(test_org_base);
    TestResult result;
    if (!prob_utils_LastIndexOfZero.submitted) {
      result.score = 0;
      result.pass = false;
      result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_LastIndexOfZero.CalcScoreGradient(test_org.GetCorrectOut(), prob_utils_LastIndexOfZero.submitted_val));
      result.score = r.first;
      result.pass = r.second;
      result.sub = true;
    }
    return result;
  };

  // Setup how evaluation on world test should work.
  EvaluateWorldTest = [this](prog_org_t & prog_org, size_t testID) {
    emp::Ptr<test_org_t> test_org_ptr = prob_LastIndexOfZero_world->GetOrgPtr(testID);
    begin_program_test.Trigger(prog_org, test_org_ptr);
    do_program_test.Trigger(prog_org, test_org_ptr);
    end_program_test.Trigger(prog_org, test_org_ptr);
    return CalcProgramResultOnTest(prog_org, *test_org_ptr);
  };

  // How should we validate programs on testing set?
  prob_utils_LastIndexOfZero.population_validation_outputs.resize(PROG_POP_SIZE);
  DoTestingSetValidation = [this](prog_org_t & prog_org) { 
    // evaluate program on full testing set; update stats utils with results
    begin_program_eval.Trigger(prog_org);
    stats_util.current_program__validation__test_results.resize(prob_utils_LastIndexOfZero.testingset_pop.size());
    stats_util.current_program__validation__total_score = 0;
    stats_util.current_program__validation__total_passes = 0;
    stats_util.current_program__validation__is_solution = false;
    prob_utils_LastIndexOfZero.population_validation_outputs[stats_util.cur_progID].resize(prob_utils_LastIndexOfZero.testingset_pop.size());
    // For each test in validation set, evaluate program.
    for (size_t testID = 0; testID < prob_utils_LastIndexOfZero.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_LastIndexOfZero.testingset_pop[testID];
      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      stats_util.current_program__validation__test_results[testID] = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      stats_util.current_program__validation__total_score += stats_util.current_program__validation__test_results[testID].score;
      stats_util.current_program__validation__total_passes += (size_t)stats_util.current_program__validation__test_results[testID].pass;
      prob_utils_LastIndexOfZero.population_validation_outputs[stats_util.cur_progID][testID] = prob_utils_LastIndexOfZero.submitted_val;
    }
    stats_util.current_program__validation__is_solution = stats_util.current_program__validation__total_passes == prob_utils_LastIndexOfZero.testingset_pop.size();
    end_program_eval.Trigger(prog_org);
  };
  program_stats.get_prog_behavioral_diversity = [this]() { return emp::ShannonEntropy(prob_utils_LastIndexOfZero.population_validation_outputs); };
  program_stats.get_prog_unique_behavioral_phenotypes = [this]() { return emp::UniqueCount(prob_utils_LastIndexOfZero.population_validation_outputs); };

  // How should we screen for a solution?
  ScreenForSolution = [this](prog_org_t & prog_org) {
    begin_program_eval.Trigger(prog_org);
    for (size_t testID = 0; testID < prob_utils_LastIndexOfZero.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_LastIndexOfZero.testingset_pop[testID];

      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      
      TestResult result = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      if (!result.pass) {
        end_program_eval.Trigger(prog_org);
        return false;
      }
    }
    end_program_eval.Trigger(prog_org);
    return true;
  }; 

  // Tell the experiment how to get test phenotypes.
  GetTestPhenotype = [this](size_t testID) -> test_org_phen_t & {
    emp_assert(prob_LastIndexOfZero_world->IsOccupied(testID));
    return prob_LastIndexOfZero_world->GetOrg(testID).GetPhenotype();
  };

  // Setup how test world updates.
  SetupTestCaseWorldUpdate(prob_LastIndexOfZero_world);

  // Setup how test cases mutate.
  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::RANDOM) {
    std::cout << "RANDOM training mode detected, configuring mutation function to RANDOMIZE organisms." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Randomize organism genome on mutate.
      prob_LastIndexOfZero_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        test_org.GetGenome() = GenRandomTestInput_LastIndexOfZero(*random, 
                                                                  {PROB_LAST_INDEX_OF_ZERO__MIN_VEC_LEN, PROB_LAST_INDEX_OF_ZERO__MAX_VEC_LEN},
                                                                  {PROB_LAST_INDEX_OF_ZERO__MIN_NUM, PROB_LAST_INDEX_OF_ZERO__MAX_NUM}); 
        return 1;
      });
    };
  } else {
    std::cout << "Non-RANDOM training mode detected, configuring mutation function normally." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_LastIndexOfZero.MIN_VEC_LEN = PROB_LAST_INDEX_OF_ZERO__MIN_VEC_LEN;
      prob_utils_LastIndexOfZero.MAX_VEC_LEN = PROB_LAST_INDEX_OF_ZERO__MAX_VEC_LEN;
      prob_utils_LastIndexOfZero.MIN_NUM = PROB_LAST_INDEX_OF_ZERO__MIN_NUM;
      prob_utils_LastIndexOfZero.MAX_NUM = PROB_LAST_INDEX_OF_ZERO__MAX_NUM;
      prob_utils_LastIndexOfZero.PER_NUM_SWAP_RATE = PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_SWAP_RATE;
      prob_utils_LastIndexOfZero.PER_NUM_DEL_RATE = PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_DEL_RATE;
      prob_utils_LastIndexOfZero.PER_NUM_INS_RATE = PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_INS_RATE;
      prob_utils_LastIndexOfZero.PER_NUM_SUB_RATE = PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_SUB_RATE;
      // (2) Hook mutator up to world.
      prob_LastIndexOfZero_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_LastIndexOfZero.Mutate(rnd, test_org.GetGenome());
      });
    };
  }

  // Setup test case fitness function.
  SetupTestFitFun = [this]() {
    prob_LastIndexOfZero_world->SetFitFun([](test_org_t & test_org) {
      return (double)test_org.GetPhenotype().num_fails;
    });
  };

  // Tell experiment how to configure hardware inputs when running program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org, emp::Ptr<TestOrg_Base> test_org_base_ptr) {
    // Reset eval stuff
    // Set current test org.
    prob_utils_LastIndexOfZero.cur_eval_test_org = test_org_base_ptr.Cast<test_org_t>(); // currently only place need testID for this?
    prob_utils_LastIndexOfZero.ResetTestEval();
    prob_utils_LastIndexOfZero.MAX_ERROR = prob_utils_LastIndexOfZero.cur_eval_test_org->GetGenome().size();
    emp_assert(eval_hardware->GetMemSize() >= 3);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      // Grab some useful references.
      Problem_LastIndexOfZero_input_t & input = prob_utils_LastIndexOfZero.cur_eval_test_org->GetGenome();
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      // Set hardware input.
      wmem.Set(0, input);
    }
  });

  // Tell experiment how to snapshot test population.
  SnapshotTests = [this]() {
    std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate());
    mkdir(snapshot_dir.c_str(), ACCESSPERMS);
    
    emp::DataFile file(snapshot_dir + "/test_pop_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");
    // Test file contents:
    // - test id
    std::function<size_t(void)> get_test_id = [this]() { return stats_util.cur_testID; };
    file.AddFun(get_test_id, "test_id");

    // - test fitness
    std::function<double(void)> get_test_fitness = [this]() { return prob_LastIndexOfZero_world->CalcFitnessID(stats_util.cur_testID); };
    file.AddFun(get_test_fitness, "fitness");

    // - num passes
    std::function<size_t(void)> get_test_num_passes = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_passes; };
    file.AddFun(get_test_num_passes, "num_passes");

    // - num fails
    std::function<size_t(void)> get_test_num_fails = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_fails; };
    file.AddFun(get_test_num_fails, "num_fails");

    std::function<size_t(void)> get_num_tested = [this]() { return GetTestPhenotype(stats_util.cur_testID).test_passes.size(); };
    file.AddFun(get_num_tested, "num_programs_tested_against");

    // - test scores by program
    std::function<std::string(void)> get_passes_by_program = [this]() {
      std::string scores = "\"[";
      test_org_phen_t & phen = GetTestPhenotype(stats_util.cur_testID);
      for (size_t i = 0; i < phen.test_passes.size(); ++i) {
        if (i) scores += ",";
        scores += emp::to_string(phen.test_passes[i]);
      }
      scores += "]\"";
      return scores;
    };
    file.AddFun(get_passes_by_program, "passes_by_program");

    // - test
    std::function<std::string(void)> get_test = [this]() {
      std::ostringstream stream;
      stream << "\"";
      prob_LastIndexOfZero_world->GetOrg(stats_util.cur_testID).Print(stream);
      stream << "\"";
      return stream.str();
    };
    file.AddFun(get_test, "test");

    file.PrintHeaderKeys();

    // Loop over tests, snapshotting each.
    for (stats_util.cur_testID = 0; stats_util.cur_testID < prob_LastIndexOfZero_world->GetSize(); ++stats_util.cur_testID) {
      if (!prob_LastIndexOfZero_world->IsOccupied(stats_util.cur_testID)) continue;
      file.Update();
    }
  };

  AddDefaultInstructions({"Add",
                          "Sub",
                          "Mult",
                          "Div",
                          "Mod",
                          "TestNumEqu",
                          "TestNumNEqu",
                          "TestNumLess",
                          "TestNumLessTEqu",
                          "TestNumGreater",
                          "TestNumGreaterTEqu",
                          "Floor",
                          "Not",
                          "Inc",
                          "Dec",
                          "CopyMem",
                          "SwapMem",
                          "Input",
                          "Output",
                          "CommitGlobal",
                          "PullGlobal",
                          "TestMemEqu",
                          "TestMemNEqu",
                          "If",
                          "IfNot",
                          "While",
                          "Countdown",
                          "Foreach",
                          "Close",
                          "Break",
                          "Call",
                          "Routine",
                          "Return",
                          "ModuleDef",
                          "MakeVector",
                          "VecGet",
                          "VecSet",
                          "VecLen",
                          "VecAppend",
                          "VecPop",
                          "VecRemove",
                          "VecReplaceAll",
                          "VecIndexOf",
                          "VecOccurrencesOf",
                          "VecReverse",
                          "VecSwapIfLess",
                          "VecGetFront",
                          "VecGetBack",
                          "IsNum",
                          "IsVec"
  });

  inst_lib->AddInst("LoadVec", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadVec_LastIndexOfZero(hw, inst);
  }, 1);

  inst_lib->AddInst("SubmitNum", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitNum_LastIndexOfZero(hw, inst);
  }, 1);

  // Add Terminals
  for (size_t i = 0; i <= 10; ++i) {
    inst_lib->AddInst("Set-" + emp::to_string(i),
      [i](hardware_t & hw, const inst_t & inst) {
        hardware_t::CallState & state = hw.GetCurCallState();
        hardware_t::Memory & wmem = state.GetWorkingMem();
        size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
        if (!hw.IsValidMemPos(posA)) return; // Do nothing
        wmem.Set(posA, (double)i);
      });
  }

}

void ProgramSynthesisExperiment::SetupProblem_VectorAverage() { 
    std::cout << "Setting up problem - VectorAverage" << std::endl;

  // A few useful aliases
  using test_org_t = TestOrg_VectorAverage;

  // Load benchmark data for problem.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  std::cout << "Loading training examples." << std::endl;
  prob_utils_VectorAverage.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  std::cout << "Loading testing examples." << std::endl;
  prob_utils_VectorAverage.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  std::cout << "Generating testing set population." << std::endl;
  prob_utils_VectorAverage.GenerateTestingSetPop();
  std::cout << "Loaded training example set size = " << prob_utils_VectorAverage.GetTrainingSet().GetSize() << std::endl;
  std::cout << "Loaded testing example set size = " << prob_utils_VectorAverage.GetTestingSet().GetSize() << std::endl;
  std::cout << "Testing set (non-training examples used to evaluate program accuracy) size = " << prob_utils_VectorAverage.testingset_pop.size() << std::endl;

  // Setup epsilon
  prob_utils_VectorAverage.EPSILON = PROB_VECTOR_AVERAGE__EPSILON;

  // Setup the world
  NewTestCaseWorld(prob_VectorAverage_world, *random, "VectorAverage world");

  // Configure how the population should be initialized
  SetupTestCasePop_Init(prob_VectorAverage_world,
                        prob_utils_VectorAverage.training_set,
                        [this]() { return GenRandomTestInput_VectorAverage(*random, 
                                                                           {PROB_VECTOR_AVERAGE__MIN_VEC_LEN, PROB_VECTOR_AVERAGE__MAX_VEC_LEN},
                                                                           {PROB_VECTOR_AVERAGE__MIN_NUM, PROB_VECTOR_AVERAGE__MAX_NUM}); 
                                  }
                        );
  end_setup_sig.AddAction([this]() { std::cout << "TestCase world size= " << prob_VectorAverage_world->GetSize() << std::endl; });

  // Tell the world to calculate the correct test output (given input) on placement.
  prob_VectorAverage_world->OnPlacement([this](size_t pos) { prob_VectorAverage_world->GetOrg(pos).CalcOut(); });

  // How are program results calculated on a test?
  CalcProgramResultOnTest = [this](prog_org_t & prog_org, TestOrg_Base & test_org_base) {
    test_org_t & test_org = static_cast<test_org_t&>(test_org_base);
    TestResult result;
    if (!prob_utils_VectorAverage.submitted) {
      result.score = 0;
      result.pass = false;
      result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_VectorAverage.CalcScoreGradient(test_org.GetCorrectOut(), prob_utils_VectorAverage.submitted_val));
      result.score = r.first;
      result.pass = r.second;
      result.sub = true;
    }
    return result;
  };

  // Setup how evaluation on world test should work.
  EvaluateWorldTest = [this](prog_org_t & prog_org, size_t testID) {
    emp::Ptr<test_org_t> test_org_ptr = prob_VectorAverage_world->GetOrgPtr(testID);
    begin_program_test.Trigger(prog_org, test_org_ptr);
    do_program_test.Trigger(prog_org, test_org_ptr);
    end_program_test.Trigger(prog_org, test_org_ptr);
    return CalcProgramResultOnTest(prog_org, *test_org_ptr);
  };

  // How should we validate programs on testing set?
  prob_utils_VectorAverage.population_validation_outputs.resize(PROG_POP_SIZE);
  DoTestingSetValidation = [this](prog_org_t & prog_org) { 
    // evaluate program on full testing set; update stats utils with results
    begin_program_eval.Trigger(prog_org);
    stats_util.current_program__validation__test_results.resize(prob_utils_VectorAverage.testingset_pop.size());
    stats_util.current_program__validation__total_score = 0;
    stats_util.current_program__validation__total_passes = 0;
    stats_util.current_program__validation__is_solution = false;
    prob_utils_VectorAverage.population_validation_outputs[stats_util.cur_progID].resize(prob_utils_VectorAverage.testingset_pop.size());
    // For each test in validation set, evaluate program.
    for (size_t testID = 0; testID < prob_utils_VectorAverage.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_VectorAverage.testingset_pop[testID];
      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      stats_util.current_program__validation__test_results[testID] = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      stats_util.current_program__validation__total_score += stats_util.current_program__validation__test_results[testID].score;
      stats_util.current_program__validation__total_passes += (size_t)stats_util.current_program__validation__test_results[testID].pass;
      prob_utils_VectorAverage.population_validation_outputs[stats_util.cur_progID][testID] = prob_utils_VectorAverage.submitted_val;
    }
    stats_util.current_program__validation__is_solution = stats_util.current_program__validation__total_passes == prob_utils_VectorAverage.testingset_pop.size();
    end_program_eval.Trigger(prog_org);
  };
  program_stats.get_prog_behavioral_diversity = [this]() { return emp::ShannonEntropy(prob_utils_VectorAverage.population_validation_outputs); };
  program_stats.get_prog_unique_behavioral_phenotypes = [this]() { return emp::UniqueCount(prob_utils_VectorAverage.population_validation_outputs); };

  // How should we screen for a solution?
  ScreenForSolution = [this](prog_org_t & prog_org) {
    begin_program_eval.Trigger(prog_org);
    for (size_t testID = 0; testID < prob_utils_VectorAverage.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_VectorAverage.testingset_pop[testID];

      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      
      TestResult result = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      if (!result.pass) {
        end_program_eval.Trigger(prog_org);
        return false;
      }
    }
    end_program_eval.Trigger(prog_org);
    return true;
  }; 

  // Tell the experiment how to get test phenotypes.
  GetTestPhenotype = [this](size_t testID) -> test_org_phen_t & {
    emp_assert(prob_VectorAverage_world->IsOccupied(testID));
    return prob_VectorAverage_world->GetOrg(testID).GetPhenotype();
  };

  // Setup how test world updates.
  SetupTestCaseWorldUpdate(prob_VectorAverage_world);

  // Setup how test cases mutate.
  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::RANDOM) {
    std::cout << "RANDOM training mode detected, configuring mutation function to RANDOMIZE organisms." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Randomize organism genome on mutate.
      prob_VectorAverage_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        test_org.GetGenome() = GenRandomTestInput_VectorAverage(*random, 
                                                                {PROB_VECTOR_AVERAGE__MIN_VEC_LEN, PROB_VECTOR_AVERAGE__MAX_VEC_LEN},
                                                                {PROB_VECTOR_AVERAGE__MIN_NUM, PROB_VECTOR_AVERAGE__MAX_NUM}); 
        return 1;
      });
    };
  } else {
    std::cout << "Non-RANDOM training mode detected, configuring mutation function normally." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_VectorAverage.MIN_VEC_LEN = PROB_VECTOR_AVERAGE__MIN_VEC_LEN;
      prob_utils_VectorAverage.MAX_VEC_LEN = PROB_VECTOR_AVERAGE__MAX_VEC_LEN;
      prob_utils_VectorAverage.MIN_NUM = PROB_VECTOR_AVERAGE__MIN_NUM;
      prob_utils_VectorAverage.MAX_NUM = PROB_VECTOR_AVERAGE__MAX_NUM;
      prob_utils_VectorAverage.DEL_RATE = PROB_VECTOR_AVERAGE__MUTATION__DEL_RATE;
      prob_utils_VectorAverage.INS_RATE = PROB_VECTOR_AVERAGE__MUTATION__INS_RATE;
      prob_utils_VectorAverage.SUB_RATE = PROB_VECTOR_AVERAGE__MUTATION__SUB_RATE;
      // (2) Hook mutator up to world.
      prob_VectorAverage_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_VectorAverage.Mutate(rnd, test_org.GetGenome());
      });
    };
  }

  // Setup test case fitness function.
  SetupTestFitFun = [this]() {
    prob_VectorAverage_world->SetFitFun([](test_org_t & test_org) {
      return (double)test_org.GetPhenotype().num_fails;
    });
  };

  // Tell experiment how to configure hardware inputs when running program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org, emp::Ptr<TestOrg_Base> test_org_base_ptr) {
    // Reset eval stuff
    // Set current test org.
    prob_utils_VectorAverage.cur_eval_test_org = test_org_base_ptr.Cast<test_org_t>(); // currently only place need testID for this?
    prob_utils_VectorAverage.ResetTestEval();
    prob_utils_VectorAverage.MAX_ERROR = prob_utils_VectorAverage.cur_eval_test_org->GetGenome().size() * PROB_VECTOR_AVERAGE__MAX_NUM;
    emp_assert(eval_hardware->GetMemSize() >= 3);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      // Grab some useful references.
      Problem_VectorAverage_input_t & input = prob_utils_VectorAverage.cur_eval_test_org->GetGenome(); // std::pair<int, double>
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      // Set hardware input.
      wmem.Set(0, input);
    }
  });

  // Tell experiment how to snapshot test population.
  SnapshotTests = [this]() {
    std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate());
    mkdir(snapshot_dir.c_str(), ACCESSPERMS);
    
    emp::DataFile file(snapshot_dir + "/test_pop_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");
    // Test file contents:
    // - test id
    std::function<size_t(void)> get_test_id = [this]() { return stats_util.cur_testID; };
    file.AddFun(get_test_id, "test_id");

    // - test fitness
    std::function<double(void)> get_test_fitness = [this]() { return prob_VectorAverage_world->CalcFitnessID(stats_util.cur_testID); };
    file.AddFun(get_test_fitness, "fitness");

    // - num passes
    std::function<size_t(void)> get_test_num_passes = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_passes; };
    file.AddFun(get_test_num_passes, "num_passes");

    // - num fails
    std::function<size_t(void)> get_test_num_fails = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_fails; };
    file.AddFun(get_test_num_fails, "num_fails");

    std::function<size_t(void)> get_num_tested = [this]() { return GetTestPhenotype(stats_util.cur_testID).test_passes.size(); };
    file.AddFun(get_num_tested, "num_programs_tested_against");

    // - test scores by program
    std::function<std::string(void)> get_passes_by_program = [this]() {
      std::string scores = "\"[";
      test_org_phen_t & phen = GetTestPhenotype(stats_util.cur_testID);
      for (size_t i = 0; i < phen.test_passes.size(); ++i) {
        if (i) scores += ",";
        scores += emp::to_string(phen.test_passes[i]);
      }
      scores += "]\"";
      return scores;
    };
    file.AddFun(get_passes_by_program, "passes_by_program");

    // - test
    std::function<std::string(void)> get_test = [this]() {
      std::ostringstream stream;
      stream << "\"";
      prob_VectorAverage_world->GetOrg(stats_util.cur_testID).Print(stream);
      stream << "\"";
      return stream.str();
    };
    file.AddFun(get_test, "test");

    file.PrintHeaderKeys();

    // Loop over tests, snapshotting each.
    for (stats_util.cur_testID = 0; stats_util.cur_testID < prob_VectorAverage_world->GetSize(); ++stats_util.cur_testID) {
      if (!prob_VectorAverage_world->IsOccupied(stats_util.cur_testID)) continue;
      file.Update();
    }
  };

  AddDefaultInstructions({"Add",
                          "Sub",
                          "Mult",
                          "Div",
                          "Mod",
                          "TestNumEqu",
                          "TestNumNEqu",
                          "TestNumLess",
                          "TestNumLessTEqu",
                          "TestNumGreater",
                          "TestNumGreaterTEqu",
                          "Floor",
                          "Not",
                          "Inc",
                          "Dec",
                          "CopyMem",
                          "SwapMem",
                          "Input",
                          "Output",
                          "CommitGlobal",
                          "PullGlobal",
                          "TestMemEqu",
                          "TestMemNEqu",
                          "If",
                          "IfNot",
                          "While",
                          "Countdown",
                          "Foreach",
                          "Close",
                          "Break",
                          "Call",
                          "Routine",
                          "Return",
                          "ModuleDef",
                          "MakeVector",
                          "VecGet",
                          "VecSet",
                          "VecLen",
                          "VecAppend",
                          "VecPop",
                          "VecRemove",
                          "VecReplaceAll",
                          "VecIndexOf",
                          "VecOccurrencesOf",
                          "VecReverse",
                          "VecSwapIfLess",
                          "VecGetFront",
                          "VecGetBack",
                          "IsNum",
                          "IsVec"
  });

  inst_lib->AddInst("LoadVec", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadVec_VectorAverage(hw, inst);
  }, 1);

  inst_lib->AddInst("SubmitNum", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitNum_VectorAverage(hw, inst);
  }, 1);

  // Add Terminals
  for (size_t i = 0; i <= 10; ++i) {
    inst_lib->AddInst("Set-" + emp::to_string(i),
      [i](hardware_t & hw, const inst_t & inst) {
        hardware_t::CallState & state = hw.GetCurCallState();
        hardware_t::Memory & wmem = state.GetWorkingMem();
        size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
        if (!hw.IsValidMemPos(posA)) return; // Do nothing
        wmem.Set(posA, (double)i);
      });
  }
}

void ProgramSynthesisExperiment::SetupProblem_CountOdds() { 
  std::cout << "Setting up problem - CountOdds" << std::endl;

  // A few useful aliases
  using test_org_t = TestOrg_CountOdds;

  // Load benchmark data for problem.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  std::cout << "Loading training examples." << std::endl;
  prob_utils_CountOdds.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  std::cout << "Loading testing examples." << std::endl;
  prob_utils_CountOdds.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  std::cout << "Generating testing set population." << std::endl;
  prob_utils_CountOdds.GenerateTestingSetPop();
  std::cout << "Loaded training example set size = " << prob_utils_CountOdds.GetTrainingSet().GetSize() << std::endl;
  std::cout << "Loaded testing example set size = " << prob_utils_CountOdds.GetTestingSet().GetSize() << std::endl;
  std::cout << "Testing set (non-training examples used to evaluate program accuracy) size = " << prob_utils_CountOdds.testingset_pop.size() << std::endl;

  // Setup the world
  NewTestCaseWorld(prob_CountOdds_world, *random, "CountOdds world");

  // Configure how the population should be initialized
  SetupTestCasePop_Init(prob_CountOdds_world,
                        prob_utils_CountOdds.training_set,
                        [this]() { return GenRandomTestInput_CountOdds(*random, 
                                                                             {PROB_COUNT_ODDS__MIN_VEC_LEN, PROB_COUNT_ODDS__MAX_VEC_LEN},
                                                                             {PROB_COUNT_ODDS__MIN_NUM, PROB_COUNT_ODDS__MAX_NUM}); 
                                  }
                        );
  end_setup_sig.AddAction([this]() { std::cout << "TestCase world size= " << prob_CountOdds_world->GetSize() << std::endl; });

  // Tell the world to calculate the correct test output (given input) on placement.
  prob_CountOdds_world->OnPlacement([this](size_t pos) { prob_CountOdds_world->GetOrg(pos).CalcOut(); });

  // How are program results calculated on a test?
  CalcProgramResultOnTest = [this](prog_org_t & prog_org, TestOrg_Base & test_org_base) {
    test_org_t & test_org = static_cast<test_org_t&>(test_org_base);
    TestResult result;
    if (!prob_utils_CountOdds.submitted) {
      result.score = 0;
      result.pass = false;
      result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_CountOdds.CalcScoreGradient(test_org.GetCorrectOut(), prob_utils_CountOdds.submitted_val));
      result.score = r.first;
      result.pass = r.second;
      result.sub = true;
    }
    return result;
  };

  // Setup how evaluation on world test should work.
  EvaluateWorldTest = [this](prog_org_t & prog_org, size_t testID) {
    emp::Ptr<test_org_t> test_org_ptr = prob_CountOdds_world->GetOrgPtr(testID);
    begin_program_test.Trigger(prog_org, test_org_ptr);
    do_program_test.Trigger(prog_org, test_org_ptr);
    end_program_test.Trigger(prog_org, test_org_ptr);
    return CalcProgramResultOnTest(prog_org, *test_org_ptr);
  };

  // How should we validate programs on testing set?
  prob_utils_CountOdds.population_validation_outputs.resize(PROG_POP_SIZE);
  DoTestingSetValidation = [this](prog_org_t & prog_org) { 
    // evaluate program on full testing set; update stats utils with results
    begin_program_eval.Trigger(prog_org);
    stats_util.current_program__validation__test_results.resize(prob_utils_CountOdds.testingset_pop.size());
    stats_util.current_program__validation__total_score = 0;
    stats_util.current_program__validation__total_passes = 0;
    stats_util.current_program__validation__is_solution = false;
    prob_utils_CountOdds.population_validation_outputs[stats_util.cur_progID].resize(prob_utils_CountOdds.testingset_pop.size());
    // For each test in validation set, evaluate program.
    for (size_t testID = 0; testID < prob_utils_CountOdds.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_CountOdds.testingset_pop[testID];
      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      stats_util.current_program__validation__test_results[testID] = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      stats_util.current_program__validation__total_score += stats_util.current_program__validation__test_results[testID].score;
      stats_util.current_program__validation__total_passes += (size_t)stats_util.current_program__validation__test_results[testID].pass;
      prob_utils_CountOdds.population_validation_outputs[stats_util.cur_progID][testID] = prob_utils_CountOdds.submitted_val;
    }
    stats_util.current_program__validation__is_solution = stats_util.current_program__validation__total_passes == prob_utils_CountOdds.testingset_pop.size();
    end_program_eval.Trigger(prog_org);
  };
  program_stats.get_prog_behavioral_diversity = [this]() { return emp::ShannonEntropy(prob_utils_CountOdds.population_validation_outputs); };
  program_stats.get_prog_unique_behavioral_phenotypes = [this]() { return emp::UniqueCount(prob_utils_CountOdds.population_validation_outputs); };

  // How should we screen for a solution?
  ScreenForSolution = [this](prog_org_t & prog_org) {
    begin_program_eval.Trigger(prog_org);
    for (size_t testID = 0; testID < prob_utils_CountOdds.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_CountOdds.testingset_pop[testID];

      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      
      TestResult result = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      if (!result.pass) {
        end_program_eval.Trigger(prog_org);
        return false;
      }
    }
    end_program_eval.Trigger(prog_org);
    return true;
  }; 

  // Tell the experiment how to get test phenotypes.
  GetTestPhenotype = [this](size_t testID) -> test_org_phen_t & {
    emp_assert(prob_CountOdds_world->IsOccupied(testID));
    return prob_CountOdds_world->GetOrg(testID).GetPhenotype();
  };

  // Setup how test world updates.
  SetupTestCaseWorldUpdate(prob_CountOdds_world);

  // Setup how test cases mutate.
  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::RANDOM) {
    std::cout << "RANDOM training mode detected, configuring mutation function to RANDOMIZE organisms." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Randomize organism genome on mutate.
      prob_CountOdds_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        test_org.GetGenome() = GenRandomTestInput_CountOdds(*random, 
                                                                  {PROB_COUNT_ODDS__MIN_VEC_LEN, PROB_COUNT_ODDS__MAX_VEC_LEN},
                                                                  {PROB_COUNT_ODDS__MIN_NUM, PROB_COUNT_ODDS__MAX_NUM}); 
        return 1;
      });
    };
  } else {
    std::cout << "Non-RANDOM training mode detected, configuring mutation function normally." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_CountOdds.MIN_VEC_LEN = PROB_COUNT_ODDS__MIN_VEC_LEN;
      prob_utils_CountOdds.MAX_VEC_LEN = PROB_COUNT_ODDS__MAX_VEC_LEN;
      prob_utils_CountOdds.MIN_NUM = PROB_COUNT_ODDS__MIN_NUM;
      prob_utils_CountOdds.MAX_NUM = PROB_COUNT_ODDS__MAX_NUM;
      prob_utils_CountOdds.PER_NUM_SWAP_RATE = PROB_COUNT_ODDS__MUTATION__PER_NUM_SWAP_RATE;
      prob_utils_CountOdds.PER_NUM_DEL_RATE = PROB_COUNT_ODDS__MUTATION__PER_NUM_DEL_RATE;
      prob_utils_CountOdds.PER_NUM_INS_RATE = PROB_COUNT_ODDS__MUTATION__PER_NUM_INS_RATE;
      prob_utils_CountOdds.PER_NUM_SUB_RATE = PROB_COUNT_ODDS__MUTATION__PER_NUM_SUB_RATE;
      // (2) Hook mutator up to world.
      prob_CountOdds_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_CountOdds.Mutate(rnd, test_org.GetGenome());
      });
    };
  }

  // Setup test case fitness function.
  SetupTestFitFun = [this]() {
    prob_CountOdds_world->SetFitFun([](test_org_t & test_org) {
      return (double)test_org.GetPhenotype().num_fails;
    });
  };

  // Tell experiment how to configure hardware inputs when running program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org, emp::Ptr<TestOrg_Base> test_org_base_ptr) {
    // Reset eval stuff
    // Set current test org.
    prob_utils_CountOdds.cur_eval_test_org = test_org_base_ptr.Cast<test_org_t>(); // currently only place need testID for this?
    prob_utils_CountOdds.ResetTestEval();
    prob_utils_CountOdds.MAX_ERROR = prob_utils_CountOdds.cur_eval_test_org->GetGenome().size();
    emp_assert(eval_hardware->GetMemSize() >= 3);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      // Grab some useful references.
      Problem_CountOdds_input_t & input = prob_utils_CountOdds.cur_eval_test_org->GetGenome(); // std::pair<int, double>
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      // Set hardware input.
      wmem.Set(0, input);
    }
  });

  // Tell experiment how to snapshot test population.
  SnapshotTests = [this]() {
    std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate());
    mkdir(snapshot_dir.c_str(), ACCESSPERMS);
    
    emp::DataFile file(snapshot_dir + "/test_pop_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");
    // Test file contents:
    // - test id
    std::function<size_t(void)> get_test_id = [this]() { return stats_util.cur_testID; };
    file.AddFun(get_test_id, "test_id");

    // - test fitness
    std::function<double(void)> get_test_fitness = [this]() { return prob_CountOdds_world->CalcFitnessID(stats_util.cur_testID); };
    file.AddFun(get_test_fitness, "fitness");

    // - num passes
    std::function<size_t(void)> get_test_num_passes = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_passes; };
    file.AddFun(get_test_num_passes, "num_passes");

    // - num fails
    std::function<size_t(void)> get_test_num_fails = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_fails; };
    file.AddFun(get_test_num_fails, "num_fails");

    std::function<size_t(void)> get_num_tested = [this]() { return GetTestPhenotype(stats_util.cur_testID).test_passes.size(); };
    file.AddFun(get_num_tested, "num_programs_tested_against");

    // - test scores by program
    std::function<std::string(void)> get_passes_by_program = [this]() {
      std::string scores = "\"[";
      test_org_phen_t & phen = GetTestPhenotype(stats_util.cur_testID);
      for (size_t i = 0; i < phen.test_passes.size(); ++i) {
        if (i) scores += ",";
        scores += emp::to_string(phen.test_passes[i]);
      }
      scores += "]\"";
      return scores;
    };
    file.AddFun(get_passes_by_program, "passes_by_program");

    // - test
    std::function<std::string(void)> get_test = [this]() {
      std::ostringstream stream;
      stream << "\"";
      prob_CountOdds_world->GetOrg(stats_util.cur_testID).Print(stream);
      stream << "\"";
      return stream.str();
    };
    file.AddFun(get_test, "test");

    file.PrintHeaderKeys();

    // Loop over tests, snapshotting each.
    for (stats_util.cur_testID = 0; stats_util.cur_testID < prob_CountOdds_world->GetSize(); ++stats_util.cur_testID) {
      if (!prob_CountOdds_world->IsOccupied(stats_util.cur_testID)) continue;
      file.Update();
    }
  };

  AddDefaultInstructions({"Add",
                          "Sub",
                          "Mult",
                          "Div",
                          "Mod",
                          "TestNumEqu",
                          "TestNumNEqu",
                          "TestNumLess",
                          "TestNumLessTEqu",
                          "TestNumGreater",
                          "TestNumGreaterTEqu",
                          "Floor",
                          "Not",
                          "Inc",
                          "Dec",
                          "CopyMem",
                          "SwapMem",
                          "Input",
                          "Output",
                          "CommitGlobal",
                          "PullGlobal",
                          "TestMemEqu",
                          "TestMemNEqu",
                          "If",
                          "IfNot",
                          "While",
                          "Countdown",
                          "Foreach",
                          "Close",
                          "Break",
                          "Call",
                          "Routine",
                          "Return",
                          "ModuleDef",
                          "MakeVector",
                          "VecGet",
                          "VecSet",
                          "VecLen",
                          "VecAppend",
                          "VecPop",
                          "VecRemove",
                          "VecReplaceAll",
                          "VecIndexOf",
                          "VecOccurrencesOf",
                          "VecReverse",
                          "VecSwapIfLess",
                          "VecGetFront",
                          "VecGetBack",
                          "IsNum",
                          "IsVec"
  });

  inst_lib->AddInst("LoadVec", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadVec_CountOdds(hw, inst);
  }, 1);

  inst_lib->AddInst("SubmitNum", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitNum_CountOdds(hw, inst);
  }, 1);

  // Add Terminals
  for (size_t i = 0; i <= 10; ++i) {
    inst_lib->AddInst("Set-" + emp::to_string(i),
      [i](hardware_t & hw, const inst_t & inst) {
        hardware_t::CallState & state = hw.GetCurCallState();
        hardware_t::Memory & wmem = state.GetWorkingMem();
        size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
        if (!hw.IsValidMemPos(posA)) return; // Do nothing
        wmem.Set(posA, (double)i);
      });
  }
}

void ProgramSynthesisExperiment::SetupProblem_MirrorImage() { 
  std::cout << "Setting up problem - MirrorImage" << std::endl;

  // A few useful aliases
  using test_org_t = TestOrg_MirrorImage;

  // Load benchmark data for problem.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  std::cout << "Loading training examples." << std::endl;
  prob_utils_MirrorImage.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  std::cout << "Loading testing examples." << std::endl;
  prob_utils_MirrorImage.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  std::cout << "Generating testing set population." << std::endl;
  prob_utils_MirrorImage.GenerateTestingSetPop();
  std::cout << "Loaded training example set size = " << prob_utils_MirrorImage.GetTrainingSet().GetSize() << std::endl;
  std::cout << "Loaded testing example set size = " << prob_utils_MirrorImage.GetTestingSet().GetSize() << std::endl;
  std::cout << "Testing set (non-training examples used to evaluate program accuracy) size = " << prob_utils_MirrorImage.testingset_pop.size() << std::endl;

  // Setup the world
  NewTestCaseWorld(prob_MirrorImage_world, *random, "MirrorImage world");

  // Configure how the population should be initialized
  SetupTestCasePop_Init(prob_MirrorImage_world,
                        prob_utils_MirrorImage.training_set,
                        [this]() { return GenRandomTestInput_MirrorImage(*random, 
                                                                         {PROB_MIRROR_IMAGE__MIN_VEC_LEN, PROB_MIRROR_IMAGE__MAX_VEC_LEN},
                                                                         {PROB_MIRROR_IMAGE__MIN_NUM, PROB_MIRROR_IMAGE__MAX_NUM}); 
                                  }
                        );
  end_setup_sig.AddAction([this]() { std::cout << "TestCase world size= " << prob_MirrorImage_world->GetSize() << std::endl; });

  // Tell the world to calculate the correct test output (given input) on placement.
  prob_MirrorImage_world->OnPlacement([this](size_t pos) { prob_MirrorImage_world->GetOrg(pos).CalcOut(); });

  // How are program results calculated on a test?
  CalcProgramResultOnTest = [this](prog_org_t & prog_org, TestOrg_Base & test_org_base) {
    test_org_t & test_org = static_cast<test_org_t&>(test_org_base);
    TestResult result;
    if (!prob_utils_MirrorImage.submitted) {
      result.score = 0;
      result.pass = false;
      result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_MirrorImage.CalcScorePassFail(test_org.GetCorrectOut(), prob_utils_MirrorImage.submitted_val));
      result.score = r.first;
      result.pass = r.second;
      result.sub = true;
    }
    return result;
  };

  // Setup how evaluation on world test should work.
  EvaluateWorldTest = [this](prog_org_t & prog_org, size_t testID) {
    emp::Ptr<test_org_t> test_org_ptr = prob_MirrorImage_world->GetOrgPtr(testID);
    begin_program_test.Trigger(prog_org, test_org_ptr);
    do_program_test.Trigger(prog_org, test_org_ptr);
    end_program_test.Trigger(prog_org, test_org_ptr);
    return CalcProgramResultOnTest(prog_org, *test_org_ptr);
  };

  // How should we validate programs on testing set?
  prob_utils_MirrorImage.population_validation_outputs.resize(PROG_POP_SIZE);
  DoTestingSetValidation = [this](prog_org_t & prog_org) { 
    // evaluate program on full testing set; update stats utils with results
    begin_program_eval.Trigger(prog_org);
    stats_util.current_program__validation__test_results.resize(prob_utils_MirrorImage.testingset_pop.size());
    stats_util.current_program__validation__total_score = 0;
    stats_util.current_program__validation__total_passes = 0;
    stats_util.current_program__validation__is_solution = false;
    prob_utils_MirrorImage.population_validation_outputs[stats_util.cur_progID].resize(prob_utils_MirrorImage.testingset_pop.size());
    // For each test in validation set, evaluate program.
    for (size_t testID = 0; testID < prob_utils_MirrorImage.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_MirrorImage.testingset_pop[testID];
      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      stats_util.current_program__validation__test_results[testID] = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      stats_util.current_program__validation__total_score += stats_util.current_program__validation__test_results[testID].score;
      stats_util.current_program__validation__total_passes += (size_t)stats_util.current_program__validation__test_results[testID].pass;
      prob_utils_MirrorImage.population_validation_outputs[stats_util.cur_progID][testID] = prob_utils_MirrorImage.submitted_val;
    }
    stats_util.current_program__validation__is_solution = stats_util.current_program__validation__total_passes == prob_utils_MirrorImage.testingset_pop.size();
    end_program_eval.Trigger(prog_org);
  };
  program_stats.get_prog_behavioral_diversity = [this]() { return emp::ShannonEntropy(prob_utils_MirrorImage.population_validation_outputs); };
  program_stats.get_prog_unique_behavioral_phenotypes = [this]() { return emp::UniqueCount(prob_utils_MirrorImage.population_validation_outputs); };

  // How should we screen for a solution?
  ScreenForSolution = [this](prog_org_t & prog_org) {
    begin_program_eval.Trigger(prog_org);
    for (size_t testID = 0; testID < prob_utils_MirrorImage.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_MirrorImage.testingset_pop[testID];

      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      
      TestResult result = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      if (!result.pass) {
        end_program_eval.Trigger(prog_org);
        return false;
      }
    }
    end_program_eval.Trigger(prog_org);
    return true;
  }; 

  // Tell the experiment how to get test phenotypes.
  GetTestPhenotype = [this](size_t testID) -> test_org_phen_t & {
    emp_assert(prob_MirrorImage_world->IsOccupied(testID));
    return prob_MirrorImage_world->GetOrg(testID).GetPhenotype();
  };

  // Setup how test world updates.
  SetupTestCaseWorldUpdate(prob_MirrorImage_world);

  // Setup how test cases mutate.
  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::RANDOM) {
    std::cout << "RANDOM training mode detected, configuring mutation function to RANDOMIZE organisms." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Randomize organism genome on mutate.
      prob_MirrorImage_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        test_org.GetGenome() = GenRandomTestInput_MirrorImage(*random, 
                                                                  {PROB_MIRROR_IMAGE__MIN_VEC_LEN, PROB_MIRROR_IMAGE__MAX_VEC_LEN},
                                                                  {PROB_MIRROR_IMAGE__MIN_NUM, PROB_MIRROR_IMAGE__MAX_NUM}); 
        return 1;
      });
    };
  } else {
    std::cout << "Non-RANDOM training mode detected, configuring mutation function normally." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_MirrorImage.MIN_VEC_LEN = PROB_MIRROR_IMAGE__MIN_VEC_LEN;
      prob_utils_MirrorImage.MAX_VEC_LEN = PROB_MIRROR_IMAGE__MAX_VEC_LEN;
      prob_utils_MirrorImage.MIN_NUM = PROB_MIRROR_IMAGE__MIN_NUM;
      prob_utils_MirrorImage.MAX_NUM = PROB_MIRROR_IMAGE__MAX_NUM;
      prob_utils_MirrorImage.PER_VEC_RANDOMIZE_VAL_RATE = PROB_MIRROR_IMAGE__MUTATION__PER_VEC_RANDOMIZE_VAL_RATE;
      prob_utils_MirrorImage.PER_VEC_MIRROR_RATE = PROB_MIRROR_IMAGE__MUTATION__PER_VEC_MIRROR_RATE;
      prob_utils_MirrorImage.COPY_RATE = PROB_MIRROR_IMAGE__MUTATION__COPY_RATE;
      prob_utils_MirrorImage.INS_RATE = PROB_MIRROR_IMAGE__MUTATION__INS_RATE;
      prob_utils_MirrorImage.DEL_RATE = PROB_MIRROR_IMAGE__MUTATION__DEL_RATE;
      prob_utils_MirrorImage.PER_VEC_SHUFFLE_RATE  = PROB_MIRROR_IMAGE__MUTATION__PER_VEC_SHUFFLE_RATE;
      // (2) Hook mutator up to world.
      prob_MirrorImage_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_MirrorImage.Mutate(rnd, test_org.GetGenome());
      });
    };
  }

  // Setup test case fitness function.
  SetupTestFitFun = [this]() {
    prob_MirrorImage_world->SetFitFun([](test_org_t & test_org) {
      return (double)test_org.GetPhenotype().num_fails;
    });
  };

  // Tell experiment how to configure hardware inputs when running program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org, emp::Ptr<TestOrg_Base> test_org_base_ptr) {
    // Reset eval stuff
    // Set current test org.
    prob_utils_MirrorImage.cur_eval_test_org = test_org_base_ptr.Cast<test_org_t>(); // currently only place need testID for this?
    prob_utils_MirrorImage.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 3);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      // Grab some useful references.
      Problem_MirrorImage_input_t & input = prob_utils_MirrorImage.cur_eval_test_org->GetGenome(); // std::pair<int, double>
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      // Set hardware input.
      wmem.Set(0, input[0]);
      wmem.Set(1, input[1]);
    }
  });

  // Tell experiment how to snapshot test population.
  SnapshotTests = [this]() {
    std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate());
    mkdir(snapshot_dir.c_str(), ACCESSPERMS);
    
    emp::DataFile file(snapshot_dir + "/test_pop_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");
    // Test file contents:
    // - test id
    std::function<size_t(void)> get_test_id = [this]() { return stats_util.cur_testID; };
    file.AddFun(get_test_id, "test_id");

    // - test fitness
    std::function<double(void)> get_test_fitness = [this]() { return prob_MirrorImage_world->CalcFitnessID(stats_util.cur_testID); };
    file.AddFun(get_test_fitness, "fitness");

    // - num passes
    std::function<size_t(void)> get_test_num_passes = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_passes; };
    file.AddFun(get_test_num_passes, "num_passes");

    // - num fails
    std::function<size_t(void)> get_test_num_fails = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_fails; };
    file.AddFun(get_test_num_fails, "num_fails");

    std::function<size_t(void)> get_num_tested = [this]() { return GetTestPhenotype(stats_util.cur_testID).test_passes.size(); };
    file.AddFun(get_num_tested, "num_programs_tested_against");

    // - test scores by program
    std::function<std::string(void)> get_passes_by_program = [this]() {
      std::string scores = "\"[";
      test_org_phen_t & phen = GetTestPhenotype(stats_util.cur_testID);
      for (size_t i = 0; i < phen.test_passes.size(); ++i) {
        if (i) scores += ",";
        scores += emp::to_string(phen.test_passes[i]);
      }
      scores += "]\"";
      return scores;
    };
    file.AddFun(get_passes_by_program, "passes_by_program");

    // - test
    std::function<std::string(void)> get_test = [this]() {
      std::ostringstream stream;
      stream << "\"";
      prob_MirrorImage_world->GetOrg(stats_util.cur_testID).Print(stream);
      stream << "\"";
      return stream.str();
    };
    file.AddFun(get_test, "test");

    file.PrintHeaderKeys();

    // Loop over tests, snapshotting each.
    for (stats_util.cur_testID = 0; stats_util.cur_testID < prob_MirrorImage_world->GetSize(); ++stats_util.cur_testID) {
      if (!prob_MirrorImage_world->IsOccupied(stats_util.cur_testID)) continue;
      file.Update();
    }
  };

  AddDefaultInstructions({"Add",
                          "Sub",
                          "Mult",
                          "Div",
                          "Mod",
                          "TestNumEqu",
                          "TestNumNEqu",
                          "TestNumLess",
                          "TestNumLessTEqu",
                          "TestNumGreater",
                          "TestNumGreaterTEqu",
                          "Floor",
                          "Not",
                          "Inc",
                          "Dec",
                          "CopyMem",
                          "SwapMem",
                          "Input",
                          "Output",
                          "CommitGlobal",
                          "PullGlobal",
                          "TestMemEqu",
                          "TestMemNEqu",
                          "If",
                          "IfNot",
                          "While",
                          "Countdown",
                          "Foreach",
                          "Close",
                          "Break",
                          "Call",
                          "Routine",
                          "Return",
                          "ModuleDef",
                          "MakeVector",
                          "VecGet",
                          "VecSet",
                          "VecLen",
                          "VecAppend",
                          "VecPop",
                          "VecRemove",
                          "VecReplaceAll",
                          "VecIndexOf",
                          "VecOccurrencesOf",
                          "VecReverse",
                          "VecSwapIfLess",
                          "VecGetFront",
                          "VecGetBack",
                          "IsNum",
                          "IsVec"
  });

  inst_lib->AddInst("LoadVec1", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadVec1_MirrorImage(hw, inst);
  }, 1);

  inst_lib->AddInst("LoadVec2", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadVec2_MirrorImage(hw, inst);
  }, 1);

  inst_lib->AddInst("SubmitVal", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitVal_MirrorImage(hw, inst);
  }, 1);

  inst_lib->AddInst("SubmitTrue", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitTrue_MirrorImage(hw, inst);
  }, 1);

  inst_lib->AddInst("SubmitFalse", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitFalse_MirrorImage(hw, inst);
  }, 1);

  // Add Terminals
  for (size_t i = 0; i <= 10; ++i) {
    inst_lib->AddInst("Set-" + emp::to_string(i),
      [i](hardware_t & hw, const inst_t & inst) {
        hardware_t::CallState & state = hw.GetCurCallState();
        hardware_t::Memory & wmem = state.GetWorkingMem();
        size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
        if (!hw.IsValidMemPos(posA)) return; // Do nothing
        wmem.Set(posA, (double)i);
      });
  }
}

void ProgramSynthesisExperiment::SetupProblem_SuperAnagrams() { 
  std::cout << "Problem setup not yet implemented... Exiting." << std::endl;
  exit(-1); 
}

void ProgramSynthesisExperiment::SetupProblem_SumOfSquares() { 
  std::cout << "Setting up problem - SumOfSquares" << std::endl;

  // A few useful aliases
  using test_org_t = TestOrg_SumOfSquares;

  // Load benchmark data for problem.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  std::cout << "Loading training examples." << std::endl;
  prob_utils_SumOfSquares.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  std::cout << "Loading testing examples." << std::endl;
  prob_utils_SumOfSquares.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  std::cout << "Generating testing set population." << std::endl;
  prob_utils_SumOfSquares.GenerateTestingSetPop();
  std::cout << "Loaded training example set size = " << prob_utils_SumOfSquares.GetTrainingSet().GetSize() << std::endl;
  std::cout << "Loaded testing example set size = " << prob_utils_SumOfSquares.GetTestingSet().GetSize() << std::endl;
  std::cout << "Testing set (non-training examples used to evaluate program accuracy) size = " << prob_utils_SumOfSquares.testingset_pop.size() << std::endl;

  // Setup the world
  NewTestCaseWorld(prob_SumOfSquares_world, *random, "SumOfSquares world");

  // Configure how the population should be initialized
  SetupTestCasePop_Init(prob_SumOfSquares_world,
                        prob_utils_SumOfSquares.training_set,
                        [this]() { return GenRandomTestInput_SumOfSquares(*random, 
                                                                          {PROB_SUM_OF_SQUARES__MIN_NUM, PROB_SUM_OF_SQUARES__MAX_NUM}); 
                                  }
                        );
  end_setup_sig.AddAction([this]() { std::cout << "TestCase world size= " << prob_SumOfSquares_world->GetSize() << std::endl; });

  // Tell the world to calculate the correct test output (given input) on placement.
  prob_SumOfSquares_world->OnPlacement([this](size_t pos) { prob_SumOfSquares_world->GetOrg(pos).CalcOut(); });

  // How are program results calculated on a test?
  CalcProgramResultOnTest = [this](prog_org_t & prog_org, TestOrg_Base & test_org_base) {
    test_org_t & test_org = static_cast<test_org_t&>(test_org_base);
    TestResult result;
    if (!prob_utils_SumOfSquares.submitted) {
      result.score = 0;
      result.pass = false;
      result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_SumOfSquares.CalcScoreGradient(test_org.GetCorrectOut(), prob_utils_SumOfSquares.submitted_val));
      result.score = r.first;
      result.pass = r.second;
      result.sub = true;
    }
    return result;
  };

  // Setup how evaluation on world test should work.
  EvaluateWorldTest = [this](prog_org_t & prog_org, size_t testID) {
    emp::Ptr<test_org_t> test_org_ptr = prob_SumOfSquares_world->GetOrgPtr(testID);
    begin_program_test.Trigger(prog_org, test_org_ptr);
    do_program_test.Trigger(prog_org, test_org_ptr);
    end_program_test.Trigger(prog_org, test_org_ptr);
    return CalcProgramResultOnTest(prog_org, *test_org_ptr);
  };

  // How should we validate programs on testing set?
  prob_utils_SumOfSquares.population_validation_outputs.resize(PROG_POP_SIZE);
  DoTestingSetValidation = [this](prog_org_t & prog_org) { 
    // evaluate program on full testing set; update stats utils with results
    begin_program_eval.Trigger(prog_org);
    stats_util.current_program__validation__test_results.resize(prob_utils_SumOfSquares.testingset_pop.size());
    stats_util.current_program__validation__total_score = 0;
    stats_util.current_program__validation__total_passes = 0;
    stats_util.current_program__validation__is_solution = false;
    prob_utils_SumOfSquares.population_validation_outputs[stats_util.cur_progID].resize(prob_utils_SumOfSquares.testingset_pop.size());
    // For each test in validation set, evaluate program.
    for (size_t testID = 0; testID < prob_utils_SumOfSquares.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_SumOfSquares.testingset_pop[testID];
      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      stats_util.current_program__validation__test_results[testID] = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      stats_util.current_program__validation__total_score += stats_util.current_program__validation__test_results[testID].score;
      stats_util.current_program__validation__total_passes += (size_t)stats_util.current_program__validation__test_results[testID].pass;
      prob_utils_SumOfSquares.population_validation_outputs[stats_util.cur_progID][testID] = prob_utils_SumOfSquares.submitted_val;
    }
    stats_util.current_program__validation__is_solution = stats_util.current_program__validation__total_passes == prob_utils_SumOfSquares.testingset_pop.size();
    end_program_eval.Trigger(prog_org);
  };
  program_stats.get_prog_behavioral_diversity = [this]() { return emp::ShannonEntropy(prob_utils_SumOfSquares.population_validation_outputs); };
  program_stats.get_prog_unique_behavioral_phenotypes = [this]() { return emp::UniqueCount(prob_utils_SumOfSquares.population_validation_outputs); };

  // How should we screen for a solution?
  ScreenForSolution = [this](prog_org_t & prog_org) {
    begin_program_eval.Trigger(prog_org);
    for (size_t testID = 0; testID < prob_utils_SumOfSquares.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_SumOfSquares.testingset_pop[testID];

      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      
      TestResult result = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      if (!result.pass) {
        end_program_eval.Trigger(prog_org);
        return false;
      }
    }
    end_program_eval.Trigger(prog_org);
    return true;
  };

  // Tell the experiment how to get test phenotypes.
  GetTestPhenotype = [this](size_t testID) -> test_org_phen_t & {
    emp_assert(prob_SumOfSquares_world->IsOccupied(testID));
    return prob_SumOfSquares_world->GetOrg(testID).GetPhenotype();
  };

  // Setup how test world updates.
  SetupTestCaseWorldUpdate(prob_SumOfSquares_world);

  // todo - test that RANDOM is actually making random things every update..

  // Setup how test cases mutate.
  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::RANDOM) {
    std::cout << "RANDOM training mode detected, configuring mutation function to RANDOMIZE organisms." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Randomize organism genome on mutate.
      prob_SumOfSquares_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        test_org.GetGenome() = GenRandomTestInput_SumOfSquares(*random, {PROB_SUM_OF_SQUARES__MIN_NUM, PROB_SUM_OF_SQUARES__MAX_NUM}); 
        return 1;
      });
    };
  } else {
    std::cout << "Non-RANDOM training mode detected, configuring mutation function normally." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_SumOfSquares.MIN_NUM = PROB_SUM_OF_SQUARES__MIN_NUM;
      prob_utils_SumOfSquares.MAX_NUM = PROB_SUM_OF_SQUARES__MAX_NUM;
      prob_utils_SumOfSquares.NUM_MUT_RATE = PROB_SUM_OF_SQUARES__MUTATION__NUM_MUT_RATE;
      // (2) Hook mutator up to world.
      prob_SumOfSquares_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_SumOfSquares.Mutate(rnd, test_org.GetGenome());
      });
    };
  }

  // Setup test case fitness function.
  SetupTestFitFun = [this]() {
    prob_SumOfSquares_world->SetFitFun([](test_org_t & test_org) {
      return (double)test_org.GetPhenotype().num_fails;
    });
  };

  // Tell experiment how to configure hardware inputs when running program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org, emp::Ptr<TestOrg_Base> test_org_base_ptr) {
    // Reset eval stuff
    // Set current test org.
    prob_utils_SumOfSquares.cur_eval_test_org = test_org_base_ptr.Cast<test_org_t>(); // currently only place need testID for this?
    prob_utils_SumOfSquares.ResetTestEval();
    prob_utils_SumOfSquares.MAX_ERROR = (int)((double)GenCorrectOut_SumOfSquares(prob_utils_SumOfSquares.cur_eval_test_org->GetGenome()) * 0.5);
    emp_assert(eval_hardware->GetMemSize() >= 3);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      // Grab some useful references.
      Problem_SumOfSquares_input_t & input = prob_utils_SumOfSquares.cur_eval_test_org->GetGenome(); // std::pair<int, double>
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      // Set hardware input.
      wmem.Set(0, input);
    }
  });

  // Tell experiment how to snapshot test population.
  SnapshotTests = [this]() {
    std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate());
    mkdir(snapshot_dir.c_str(), ACCESSPERMS);
    
    emp::DataFile file(snapshot_dir + "/test_pop_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");
    // Test file contents:
    // - test id
    std::function<size_t(void)> get_test_id = [this]() { return stats_util.cur_testID; };
    file.AddFun(get_test_id, "test_id");

    // - test fitness
    std::function<double(void)> get_test_fitness = [this]() { return prob_SumOfSquares_world->CalcFitnessID(stats_util.cur_testID); };
    file.AddFun(get_test_fitness, "fitness");

    // - num passes
    std::function<size_t(void)> get_test_num_passes = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_passes; };
    file.AddFun(get_test_num_passes, "num_passes");

    // - num fails
    std::function<size_t(void)> get_test_num_fails = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_fails; };
    file.AddFun(get_test_num_fails, "num_fails");

    std::function<size_t(void)> get_num_tested = [this]() { return GetTestPhenotype(stats_util.cur_testID).test_passes.size(); };
    file.AddFun(get_num_tested, "num_programs_tested_against");

    // - test scores by program
    std::function<std::string(void)> get_passes_by_program = [this]() {
      std::string scores = "\"[";
      test_org_phen_t & phen = GetTestPhenotype(stats_util.cur_testID);
      for (size_t i = 0; i < phen.test_passes.size(); ++i) {
        if (i) scores += ",";
        scores += emp::to_string(phen.test_passes[i]);
      }
      scores += "]\"";
      return scores;
    };
    file.AddFun(get_passes_by_program, "passes_by_program");

    // - test
    std::function<std::string(void)> get_test = [this]() {
      std::ostringstream stream;
      stream << "\"";
      prob_SumOfSquares_world->GetOrg(stats_util.cur_testID).Print(stream);
      stream << "\"";
      return stream.str();
    };
    file.AddFun(get_test, "test");

    file.PrintHeaderKeys();

    // Loop over tests, snapshotting each.
    for (stats_util.cur_testID = 0; stats_util.cur_testID < prob_SumOfSquares_world->GetSize(); ++stats_util.cur_testID) {
      if (!prob_SumOfSquares_world->IsOccupied(stats_util.cur_testID)) continue;
      file.Update();
    }
  };

  // Add default instructions to instruction set.
  AddDefaultInstructions({"Add",
                          "Sub",
                          "Mult",
                          "Div",
                          "Mod",
                          "TestNumEqu",
                          "TestNumNEqu",
                          "TestNumLess",
                          "TestNumLessTEqu",
                          "TestNumGreater",
                          "TestNumGreaterTEqu",
                          "Floor",
                          "Not",
                          "Inc",
                          "Dec",
                          "CopyMem",
                          "SwapMem",
                          "Input",
                          "Output",
                          "CommitGlobal",
                          "PullGlobal",
                          "TestMemEqu",
                          "TestMemNEqu",
                          "If",
                          "IfNot",
                          "While",
                          "Countdown",
                          "Foreach",
                          "Close",
                          "Break",
                          "Call",
                          "Routine",
                          "Return",
                          "ModuleDef",
                          "MakeVector",
                          "VecGet",
                          "VecSet",
                          "VecLen",
                          "VecAppend",
                          "VecPop",
                          "VecRemove",
                          "VecReplaceAll",
                          "VecIndexOf",
                          "VecOccurrencesOf",
                          "VecReverse",
                          "VecSwapIfLess",
                          "VecGetFront",
                          "VecGetBack",
                          "IsNum",
                          "IsVec"
  });

  // -- Custom Instructions --
  inst_lib->AddInst("LoadNum", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadNum_SumOfSquares(hw, inst);
  }, 1);

  inst_lib->AddInst("SubmitNum", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitNum_SumOfSquares(hw, inst);
  }, 1);
  
  // Add Terminals
  for (size_t i = 0; i <= 10; ++i) {
    inst_lib->AddInst("Set-" + emp::to_string(i),
      [i](hardware_t & hw, const inst_t & inst) {
        hardware_t::CallState & state = hw.GetCurCallState();
        hardware_t::Memory & wmem = state.GetWorkingMem();
        size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
        if (!hw.IsValidMemPos(posA)) return; // Do nothing
        wmem.Set(posA, (double)i);
      });
  }

}

void ProgramSynthesisExperiment::SetupProblem_VectorsSummed() { 
  std::cout << "Setting up problem - VectorsSummed" << std::endl;

  // A few useful aliases
  using test_org_t = TestOrg_VectorsSummed;

  // Load benchmark data for problem.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  std::cout << "Loading training examples." << std::endl;
  prob_utils_VectorsSummed.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  std::cout << "Loading testing examples." << std::endl;
  prob_utils_VectorsSummed.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  std::cout << "Generating testing set population." << std::endl;
  prob_utils_VectorsSummed.GenerateTestingSetPop();
  std::cout << "Loaded training example set size = " << prob_utils_VectorsSummed.GetTrainingSet().GetSize() << std::endl;
  std::cout << "Loaded testing example set size = " << prob_utils_VectorsSummed.GetTestingSet().GetSize() << std::endl;
  std::cout << "Testing set (non-training examples used to evaluate program accuracy) size = " << prob_utils_VectorsSummed.testingset_pop.size() << std::endl;

  // Setup the world
  NewTestCaseWorld(prob_VectorsSummed_world, *random, "VectorsSummed world");

  // Configure how the population should be initialized
  SetupTestCasePop_Init(prob_VectorsSummed_world,
                        prob_utils_VectorsSummed.training_set,
                        [this]() { return GenRandomTestInput_VectorsSummed(*random, 
                                                                           {PROB_VECTORS_SUMMED__MIN_VEC_LEN, PROB_VECTORS_SUMMED__MAX_VEC_LEN},
                                                                           {PROB_VECTORS_SUMMED__MIN_NUM, PROB_VECTORS_SUMMED__MAX_NUM}); 
                                  }
                        );
  end_setup_sig.AddAction([this]() { std::cout << "TestCase world size= " << prob_VectorsSummed_world->GetSize() << std::endl; });

  // Tell the world to calculate the correct test output (given input) on placement.
  prob_VectorsSummed_world->OnPlacement([this](size_t pos) { prob_VectorsSummed_world->GetOrg(pos).CalcOut(); });

  // How are program results calculated on a test?
  CalcProgramResultOnTest = [this](prog_org_t & prog_org, TestOrg_Base & test_org_base) {
    test_org_t & test_org = static_cast<test_org_t&>(test_org_base);
    TestResult result;
    if (!prob_utils_VectorsSummed.submitted) {
      result.score = 0;
      result.pass = false;
      result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_VectorsSummed.CalcScoreGradient(test_org.GetCorrectOut(), prob_utils_VectorsSummed.submitted_vec));
      result.score = r.first;
      result.pass = r.second;
      result.sub = true;
    }
    return result;
  };

  // Setup how evaluation on world test should work.
  EvaluateWorldTest = [this](prog_org_t & prog_org, size_t testID) {
    emp::Ptr<test_org_t> test_org_ptr = prob_VectorsSummed_world->GetOrgPtr(testID);
    begin_program_test.Trigger(prog_org, test_org_ptr);
    do_program_test.Trigger(prog_org, test_org_ptr);
    end_program_test.Trigger(prog_org, test_org_ptr);
    return CalcProgramResultOnTest(prog_org, *test_org_ptr);
  };

  // How should we validate programs on testing set?
  prob_utils_VectorsSummed.population_validation_outputs.resize(PROG_POP_SIZE);
  DoTestingSetValidation = [this](prog_org_t & prog_org) { 
    // evaluate program on full testing set; update stats utils with results
    begin_program_eval.Trigger(prog_org);
    stats_util.current_program__validation__test_results.resize(prob_utils_VectorsSummed.testingset_pop.size());
    stats_util.current_program__validation__total_score = 0;
    stats_util.current_program__validation__total_passes = 0;
    stats_util.current_program__validation__is_solution = false;
    prob_utils_VectorsSummed.population_validation_outputs[stats_util.cur_progID].resize(prob_utils_VectorsSummed.testingset_pop.size());
    // For each test in validation set, evaluate program.
    for (size_t testID = 0; testID < prob_utils_VectorsSummed.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_VectorsSummed.testingset_pop[testID];
      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      stats_util.current_program__validation__test_results[testID] = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      stats_util.current_program__validation__total_score += stats_util.current_program__validation__test_results[testID].score;
      stats_util.current_program__validation__total_passes += (size_t)stats_util.current_program__validation__test_results[testID].pass;
      prob_utils_VectorsSummed.population_validation_outputs[stats_util.cur_progID][testID] = prob_utils_VectorsSummed.submitted_vec;
    }
    stats_util.current_program__validation__is_solution = stats_util.current_program__validation__total_passes == prob_utils_VectorsSummed.testingset_pop.size();
    end_program_eval.Trigger(prog_org);
  };
  program_stats.get_prog_behavioral_diversity = [this]() { return emp::ShannonEntropy(prob_utils_VectorsSummed.population_validation_outputs); };
  program_stats.get_prog_unique_behavioral_phenotypes = [this]() { return emp::UniqueCount(prob_utils_VectorsSummed.population_validation_outputs); };

  // How should we screen for a solution?
  ScreenForSolution = [this](prog_org_t & prog_org) {
    begin_program_eval.Trigger(prog_org);
    for (size_t testID = 0; testID < prob_utils_VectorsSummed.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_VectorsSummed.testingset_pop[testID];

      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      
      TestResult result = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      if (!result.pass) {
        end_program_eval.Trigger(prog_org);
        return false;
      }
    }
    end_program_eval.Trigger(prog_org);
    return true;
  }; 

  // Tell the experiment how to get test phenotypes.
  GetTestPhenotype = [this](size_t testID) -> test_org_phen_t & {
    emp_assert(prob_VectorsSummed_world->IsOccupied(testID));
    return prob_VectorsSummed_world->GetOrg(testID).GetPhenotype();
  };

  // Setup how test world updates.
  SetupTestCaseWorldUpdate(prob_VectorsSummed_world);

  // Setup how test cases mutate.
  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::RANDOM) {
    std::cout << "RANDOM training mode detected, configuring mutation function to RANDOMIZE organisms." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Randomize organism genome on mutate.
      prob_VectorsSummed_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        test_org.GetGenome() = GenRandomTestInput_VectorsSummed(*random, 
                                                                  {PROB_VECTORS_SUMMED__MIN_VEC_LEN, PROB_VECTORS_SUMMED__MAX_VEC_LEN},
                                                                  {PROB_VECTORS_SUMMED__MIN_NUM, PROB_VECTORS_SUMMED__MAX_NUM}); 
        return 1;
      });
    };
  } else {
    std::cout << "Non-RANDOM training mode detected, configuring mutation function normally." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_VectorsSummed.MIN_VEC_LEN = PROB_VECTORS_SUMMED__MIN_VEC_LEN;
      prob_utils_VectorsSummed.MAX_VEC_LEN = PROB_VECTORS_SUMMED__MAX_VEC_LEN;
      prob_utils_VectorsSummed.MIN_NUM = PROB_VECTORS_SUMMED__MIN_NUM;
      prob_utils_VectorsSummed.MAX_NUM = PROB_VECTORS_SUMMED__MAX_NUM;
      prob_utils_VectorsSummed.PER_NUM_SUB_RATE = PROB_VECTORS_SUMMED__MUTATION__PER_NUM_SUB_RATE;
      prob_utils_VectorsSummed.COPY_RATE = PROB_VECTORS_SUMMED__MUTATION__COPY_RATE;
      prob_utils_VectorsSummed.INS_RATE = PROB_VECTORS_SUMMED__MUTATION__INS_RATE;
      prob_utils_VectorsSummed.DEL_RATE = PROB_VECTORS_SUMMED__MUTATION__DEL_RATE;
      // (2) Hook mutator up to world.
      prob_VectorsSummed_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_VectorsSummed.Mutate(rnd, test_org.GetGenome());
      });
    };
  }

  // Setup test case fitness function.
  SetupTestFitFun = [this]() {
    prob_VectorsSummed_world->SetFitFun([](test_org_t & test_org) {
      return (double)test_org.GetPhenotype().num_fails;
    });
  };

  // Tell experiment how to configure hardware inputs when running program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org, emp::Ptr<TestOrg_Base> test_org_base_ptr) {
    // Reset eval stuff
    // Set current test org.
    prob_utils_VectorsSummed.cur_eval_test_org = test_org_base_ptr.Cast<test_org_t>(); // currently only place need testID for this?
    prob_utils_VectorsSummed.ResetTestEval();
    prob_utils_VectorsSummed.MAX_ERROR = (2*PROB_VECTORS_SUMMED__MAX_NUM) * prob_utils_VectorsSummed.cur_eval_test_org->GetGenome().size();
    emp_assert(eval_hardware->GetMemSize() >= 3);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      // Grab some useful references.
      Problem_VectorsSummed_input_t & input = prob_utils_VectorsSummed.cur_eval_test_org->GetGenome(); // std::pair<int, double>
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      // Set hardware input.
      wmem.Set(0, input[0]);
      wmem.Set(1, input[1]);
    }
  });

  // Tell experiment how to snapshot test population.
  SnapshotTests = [this]() {
    std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate());
    mkdir(snapshot_dir.c_str(), ACCESSPERMS);
    
    emp::DataFile file(snapshot_dir + "/test_pop_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");
    // Test file contents:
    // - test id
    std::function<size_t(void)> get_test_id = [this]() { return stats_util.cur_testID; };
    file.AddFun(get_test_id, "test_id");

    // - test fitness
    std::function<double(void)> get_test_fitness = [this]() { return prob_VectorsSummed_world->CalcFitnessID(stats_util.cur_testID); };
    file.AddFun(get_test_fitness, "fitness");

    // - num passes
    std::function<size_t(void)> get_test_num_passes = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_passes; };
    file.AddFun(get_test_num_passes, "num_passes");

    // - num fails
    std::function<size_t(void)> get_test_num_fails = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_fails; };
    file.AddFun(get_test_num_fails, "num_fails");

    std::function<size_t(void)> get_num_tested = [this]() { return GetTestPhenotype(stats_util.cur_testID).test_passes.size(); };
    file.AddFun(get_num_tested, "num_programs_tested_against");

    // - test scores by program
    std::function<std::string(void)> get_passes_by_program = [this]() {
      std::string scores = "\"[";
      test_org_phen_t & phen = GetTestPhenotype(stats_util.cur_testID);
      for (size_t i = 0; i < phen.test_passes.size(); ++i) {
        if (i) scores += ",";
        scores += emp::to_string(phen.test_passes[i]);
      }
      scores += "]\"";
      return scores;
    };
    file.AddFun(get_passes_by_program, "passes_by_program");

    // - test
    std::function<std::string(void)> get_test = [this]() {
      std::ostringstream stream;
      stream << "\"";
      prob_VectorsSummed_world->GetOrg(stats_util.cur_testID).Print(stream);
      stream << "\"";
      return stream.str();
    };
    file.AddFun(get_test, "test");

    file.PrintHeaderKeys();

    // Loop over tests, snapshotting each.
    for (stats_util.cur_testID = 0; stats_util.cur_testID < prob_VectorsSummed_world->GetSize(); ++stats_util.cur_testID) {
      if (!prob_VectorsSummed_world->IsOccupied(stats_util.cur_testID)) continue;
      file.Update();
    }
  };

  AddDefaultInstructions({"Add",
                          "Sub",
                          "Mult",
                          "Div",
                          "Mod",
                          "TestNumEqu",
                          "TestNumNEqu",
                          "TestNumLess",
                          "TestNumLessTEqu",
                          "TestNumGreater",
                          "TestNumGreaterTEqu",
                          "Floor",
                          "Not",
                          "Inc",
                          "Dec",
                          "CopyMem",
                          "SwapMem",
                          "Input",
                          "Output",
                          "CommitGlobal",
                          "PullGlobal",
                          "TestMemEqu",
                          "TestMemNEqu",
                          "If",
                          "IfNot",
                          "While",
                          "Countdown",
                          "Foreach",
                          "Close",
                          "Break",
                          "Call",
                          "Routine",
                          "Return",
                          "ModuleDef",
                          "MakeVector",
                          "VecGet",
                          "VecSet",
                          "VecLen",
                          "VecAppend",
                          "VecPop",
                          "VecRemove",
                          "VecReplaceAll",
                          "VecIndexOf",
                          "VecOccurrencesOf",
                          "VecReverse",
                          "VecSwapIfLess",
                          "VecGetFront",
                          "VecGetBack",
                          "IsNum",
                          "IsVec"
  });

  inst_lib->AddInst("LoadVec1", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadVec1_VectorsSummed(hw, inst);
  }, 1);

  inst_lib->AddInst("LoadVec2", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadVec2_VectorsSummed(hw, inst);
  }, 1);

  inst_lib->AddInst("SubmitVec", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitVec_VectorsSummed(hw, inst);
  }, 1);

  // Add Terminals
  for (size_t i = 0; i <= 10; ++i) {
    inst_lib->AddInst("Set-" + emp::to_string(i),
      [i](hardware_t & hw, const inst_t & inst) {
        hardware_t::CallState & state = hw.GetCurCallState();
        hardware_t::Memory & wmem = state.GetWorkingMem();
        size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
        if (!hw.IsValidMemPos(posA)) return; // Do nothing
        wmem.Set(posA, (double)i);
      });
  }
}

void ProgramSynthesisExperiment::SetupProblem_XWordLines() { 
  std::cout << "Problem setup not yet implemented... Exiting." << std::endl;
  exit(-1); 
}

void ProgramSynthesisExperiment::SetupProblem_PigLatin() { 
  std::cout << "Problem setup not yet implemented... Exiting." << std::endl;
  exit(-1); 
}

void ProgramSynthesisExperiment::SetupProblem_NegativeToZero() { 
  std::cout << "Problem setup not yet implemented... Exiting." << std::endl;
  exit(-1); 
}

void ProgramSynthesisExperiment::SetupProblem_ScrabbleScore() { 
  std::cout << "Problem setup not yet implemented... Exiting." << std::endl;
  exit(-1); 
}

void ProgramSynthesisExperiment::SetupProblem_Checksum() { 
  std::cout << "Problem setup not yet implemented... Exiting." << std::endl;
  exit(-1); 
}

void ProgramSynthesisExperiment::SetupProblem_Digits() { 
  std::cout << "Problem setup not yet implemented... Exiting." << std::endl;
  exit(-1); 
}

void ProgramSynthesisExperiment::SetupProblem_Grade() { 
  std::cout << "Setting up problem - Grade" << std::endl;

  // A few useful aliases.
  using test_org_t = TestOrg_Grade;

  // Load benchmark data for problem.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  std::cout << "Loading training examples." << std::endl;
  prob_utils_Grade.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  std::cout << "Loading testing examples." << std::endl;
  prob_utils_Grade.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  std::cout << "Generating testing set population." << std::endl;
  prob_utils_Grade.GenerateTestingSetPop();
  std::cout << "Loaded training example set size = " << prob_utils_Grade.GetTrainingSet().GetSize() << std::endl;
  std::cout << "Loaded testing example set size = " << prob_utils_Grade.GetTestingSet().GetSize() << std::endl;
  std::cout << "Testing set (non-training examples used to evaluate program accuracy) size = " << prob_utils_Grade.testingset_pop.size() << std::endl;

  // Setup the world
  NewTestCaseWorld(prob_Grade_world, *random, "Grade world");

  // Configure how the population should be initialized
  SetupTestCasePop_Init(prob_Grade_world,
                        prob_utils_Grade.training_set,
                        [this]() { return GenRandomTestInput_Grade(*random, 
                                                                   {PROB_GRADE__MIN_NUM, PROB_GRADE__MAX_NUM}); 
                                  }
                        );
  end_setup_sig.AddAction([this]() { std::cout << "TestCase world size= " << prob_Grade_world->GetSize() << std::endl; });

  // Tell the world to calculate the correct test output (given input) on placement.
  prob_Grade_world->OnPlacement([this](size_t pos) { prob_Grade_world->GetOrg(pos).CalcOut(); });

  // How are program results calculated on a test?
  CalcProgramResultOnTest = [this](prog_org_t & prog_org, TestOrg_Base & test_org_base) {
    test_org_t & test_org = static_cast<test_org_t&>(test_org_base);
    TestResult result;
    if (!prob_utils_Grade.submitted) {
      result.score = 0;
      result.pass = false;
      result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_Grade.CalcScorePassFail(test_org.GetCorrectOut(), prob_utils_Grade.submitted_str));
      result.score = r.first;
      result.pass = r.second;
      result.sub = true;
    }
    return result;
  };

  // Setup how evaluation on world test should work.
  EvaluateWorldTest = [this](prog_org_t & prog_org, size_t testID) {
    emp::Ptr<test_org_t> test_org_ptr = prob_Grade_world->GetOrgPtr(testID);
    begin_program_test.Trigger(prog_org, test_org_ptr);
    do_program_test.Trigger(prog_org, test_org_ptr);
    end_program_test.Trigger(prog_org, test_org_ptr);
    return CalcProgramResultOnTest(prog_org, *test_org_ptr);
  };

  // How should we validate programs on testing set?
  prob_utils_Grade.population_validation_outputs.resize(PROG_POP_SIZE);
  DoTestingSetValidation = [this](prog_org_t & prog_org) { 
    // evaluate program on full testing set; update stats utils with results
    begin_program_eval.Trigger(prog_org);
    stats_util.current_program__validation__test_results.resize(prob_utils_Grade.testingset_pop.size());
    stats_util.current_program__validation__total_score = 0;
    stats_util.current_program__validation__total_passes = 0;
    stats_util.current_program__validation__is_solution = false;
    prob_utils_Grade.population_validation_outputs[stats_util.cur_progID].resize(prob_utils_Grade.testingset_pop.size());
    // For each test in validation set, evaluate program.
    for (size_t testID = 0; testID < prob_utils_Grade.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_Grade.testingset_pop[testID];
      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      stats_util.current_program__validation__test_results[testID] = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      stats_util.current_program__validation__total_score += stats_util.current_program__validation__test_results[testID].score;
      stats_util.current_program__validation__total_passes += (size_t)stats_util.current_program__validation__test_results[testID].pass;
      prob_utils_Grade.population_validation_outputs[stats_util.cur_progID][testID] = prob_utils_Grade.submitted_str;
    }
    stats_util.current_program__validation__is_solution = stats_util.current_program__validation__total_passes == prob_utils_Grade.testingset_pop.size();
    end_program_eval.Trigger(prog_org);
  };
  program_stats.get_prog_behavioral_diversity = [this]() { return emp::ShannonEntropy(prob_utils_Grade.population_validation_outputs); };
  program_stats.get_prog_unique_behavioral_phenotypes = [this]() { return emp::UniqueCount(prob_utils_Grade.population_validation_outputs); };

  // How should we screen for a solution?
  ScreenForSolution = [this](prog_org_t & prog_org) {
    begin_program_eval.Trigger(prog_org);
    for (size_t testID = 0; testID < prob_utils_Grade.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_Grade.testingset_pop[testID];

      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      
      TestResult result = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      if (!result.pass) {
        end_program_eval.Trigger(prog_org);
        return false;
      }
    }
    end_program_eval.Trigger(prog_org);
    return true;
  };

  // Tell the experiment how to get test phenotypes.
  GetTestPhenotype = [this](size_t testID) -> test_org_phen_t & {
    emp_assert(prob_Grade_world->IsOccupied(testID));
    return prob_Grade_world->GetOrg(testID).GetPhenotype();
  };

  // Setup how test world updates.
  SetupTestCaseWorldUpdate(prob_Grade_world);

  // todo - test that RANDOM is actually making random things every update..

  // Setup how test cases mutate.
  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::RANDOM) {
    std::cout << "RANDOM training mode detected, configuring mutation function to RANDOMIZE organisms." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Randomize organism genome on mutate.
      prob_Grade_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        test_org.GetGenome() = GenRandomTestInput_Grade(*random, {PROB_GRADE__MIN_NUM, PROB_GRADE__MAX_NUM}); 
        return 1;
      });
    };
  } else {
    std::cout << "Non-RANDOM training mode detected, configuring mutation function normally." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_Grade.MIN_NUM = PROB_GRADE__MIN_NUM;
      prob_utils_Grade.MAX_NUM = PROB_GRADE__MAX_NUM;
      prob_utils_Grade.PER_NUM_ADJUST_RATE = PROB_GRADE__MUTATION__PER_NUM_ADJUST_RATE;
      prob_utils_Grade.PER_NUM_RANDOMIZE_RATE = PROB_GRADE__MUTATION__PER_NUM_RANDOMIZE_RATE;
      // (2) Hook mutator up to world.
      prob_Grade_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_Grade.Mutate(rnd, test_org.GetGenome());
      });
    };
  }

  // Setup test case fitness function.
  SetupTestFitFun = [this]() {
    prob_Grade_world->SetFitFun([](test_org_t & test_org) {
      return (double)test_org.GetPhenotype().num_fails;
    });
  }; 

  // Tell experiment how to configure hardware inputs when running program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org, emp::Ptr<TestOrg_Base> test_org_base_ptr) {
    // Reset eval stuff
    // Set current test org.
    prob_utils_Grade.cur_eval_test_org = test_org_base_ptr.Cast<test_org_t>(); // currently only place need testID for this?
    prob_utils_Grade.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 4);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      // Grab some useful references.
      Problem_Grade_input_t & input = prob_utils_Grade.cur_eval_test_org->GetGenome(); 
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      // std::cout << "Begin program test!" << std::endl;
      // std::cout << "  A thresh: " << input[0] << std::endl;
      // std::cout << "  B thresh: " << input[1] << std::endl;
      // std::cout << "  C thresh: " << input[2] << std::endl;
      // std::cout << "  D thresh: " << input[3] << std::endl;
      // std::cout << "  Grade: " << input[4] << std::endl;
      // Set hardware input.
      wmem.Set(0, input[0]);
      wmem.Set(1, input[1]);
      wmem.Set(2, input[2]);
      wmem.Set(3, input[3]);
      wmem.Set(4, input[4]);
    }
  });

  // Tell experiment how to snapshot test population.
  SnapshotTests = [this]() {
    std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate());
    mkdir(snapshot_dir.c_str(), ACCESSPERMS);
    
    emp::DataFile file(snapshot_dir + "/test_pop_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");
    // Test file contents:
    // - test id
    std::function<size_t(void)> get_test_id = [this]() { return stats_util.cur_testID; };
    file.AddFun(get_test_id, "test_id");

    // - test fitness
    std::function<double(void)> get_test_fitness = [this]() { return prob_Grade_world->CalcFitnessID(stats_util.cur_testID); };
    file.AddFun(get_test_fitness, "fitness");

    // - num passes
    std::function<size_t(void)> get_test_num_passes = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_passes; };
    file.AddFun(get_test_num_passes, "num_passes");

    // - num fails
    std::function<size_t(void)> get_test_num_fails = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_fails; };
    file.AddFun(get_test_num_fails, "num_fails");

    std::function<size_t(void)> get_num_tested = [this]() { return GetTestPhenotype(stats_util.cur_testID).test_passes.size(); };
    file.AddFun(get_num_tested, "num_programs_tested_against");

    // - test scores by program
    std::function<std::string(void)> get_passes_by_program = [this]() {
      std::string scores = "\"[";
      test_org_phen_t & phen = GetTestPhenotype(stats_util.cur_testID);
      for (size_t i = 0; i < phen.test_passes.size(); ++i) {
        if (i) scores += ",";
        scores += emp::to_string(phen.test_passes[i]);
      }
      scores += "]\"";
      return scores;
    };
    file.AddFun(get_passes_by_program, "passes_by_program");

    // - test
    std::function<std::string(void)> get_test = [this]() {
      std::ostringstream stream;
      stream << "\"";
      prob_Grade_world->GetOrg(stats_util.cur_testID).Print(stream);
      stream << "\"";
      return stream.str();
    };
    file.AddFun(get_test, "test", "");

    file.PrintHeaderKeys();

    // Loop over tests, snapshotting each.
    for (stats_util.cur_testID = 0; stats_util.cur_testID < prob_Grade_world->GetSize(); ++stats_util.cur_testID) {
      if (!prob_Grade_world->IsOccupied(stats_util.cur_testID)) continue;
      file.Update();
    }
  };

  // Add default instructions to instruction set.
  AddDefaultInstructions({"Add",
                          "Sub",
                          "Mult",
                          "Div",
                          "Mod",
                          "TestNumEqu",
                          "TestNumNEqu",
                          "TestNumLess",
                          "TestNumLessTEqu",
                          "TestNumGreater",
                          "TestNumGreaterTEqu",
                          "Floor",
                          "Not",
                          "Inc",
                          "Dec",
                          "CopyMem",
                          "SwapMem",
                          "Input",
                          "Output",
                          "CommitGlobal",
                          "PullGlobal",
                          "TestMemEqu",
                          "TestMemNEqu",
                          "If",
                          "IfNot",
                          "While",
                          "Countdown",
                          "Foreach",
                          "Close",
                          "Break",
                          "Call",
                          "Routine",
                          "Return",
                          "ModuleDef",
                          "MakeVector",
                          "VecGet",
                          "VecSet",
                          "VecLen",
                          "VecAppend",
                          "VecPop",
                          "VecRemove",
                          "VecReplaceAll",
                          "VecIndexOf",
                          "VecOccurrencesOf",
                          "VecReverse",
                          "VecSwapIfLess",
                          "VecGetFront",
                          "VecGetBack",
                          "IsNum",
                          "IsVec"
  });

  // -- Custom Instructions --
  inst_lib->AddInst("LoadThreshA", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadThreshA_Grade(hw, inst);
    // this->Inst_LoadNum1_Grade(hw, inst);
  }, 1);

  inst_lib->AddInst("LoadThreshB", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadThreshB_Grade(hw, inst);
  }, 1);

  inst_lib->AddInst("LoadThreshC", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadThreshC_Grade(hw, inst);
  }, 1);

  inst_lib->AddInst("LoadThreshD", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadThreshD_Grade(hw, inst);
  }, 1);

  inst_lib->AddInst("LoadGrade", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadGrade_Grade(hw, inst);
  }, 1);

  inst_lib->AddInst("SubmitA", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitA_Grade(hw, inst);
  }, 1);
  inst_lib->AddInst("SubmitB", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitB_Grade(hw, inst);
  }, 1);
  inst_lib->AddInst("SubmitC", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitC_Grade(hw, inst);
  }, 1);
  inst_lib->AddInst("SubmitD", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitD_Grade(hw, inst);
  }, 1);
  inst_lib->AddInst("SubmitF", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitF_Grade(hw, inst);
  }, 1);
}

void ProgramSynthesisExperiment::SetupProblem_Median() { 
  std::cout << "Setting up problem - Median" << std::endl;

  // A few useful aliases
  using test_org_t = TestOrg_Median;

  // Load benchmark data for problem.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  std::cout << "Loading training examples." << std::endl;
  prob_utils_Median.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  std::cout << "Loading testing examples." << std::endl;
  prob_utils_Median.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  std::cout << "Generating testing set population." << std::endl;
  prob_utils_Median.GenerateTestingSetPop();
  std::cout << "Loaded training example set size = " << prob_utils_Median.GetTrainingSet().GetSize() << std::endl;
  std::cout << "Loaded testing example set size = " << prob_utils_Median.GetTestingSet().GetSize() << std::endl;
  std::cout << "Testing set (non-training examples used to evaluate program accuracy) size = " << prob_utils_Median.testingset_pop.size() << std::endl;

  // Setup the world
  NewTestCaseWorld(prob_Median_world, *random, "Median world");

  // Configure how the population should be initialized
  SetupTestCasePop_Init(prob_Median_world,
                        prob_utils_Median.training_set,
                        [this]() { return GenRandomTestInput_Median(*random, 
                                                                    {PROB_MEDIAN__MIN_NUM, PROB_MEDIAN__MAX_NUM}); 
                                  }
                        );
  end_setup_sig.AddAction([this]() { std::cout << "TestCase world size= " << prob_Median_world->GetSize() << std::endl; });

  // Tell the world to calculate the correct test output (given input) on placement.
  prob_Median_world->OnPlacement([this](size_t pos) { prob_Median_world->GetOrg(pos).CalcOut(); });

  // How are program results calculated on a test?
  CalcProgramResultOnTest = [this](prog_org_t & prog_org, TestOrg_Base & test_org_base) {
    test_org_t & test_org = static_cast<test_org_t&>(test_org_base);
    TestResult result;
    if (!prob_utils_Median.submitted) {
      result.score = 0;
      result.pass = false;
      result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_Median.CalcScorePassFail(test_org.GetCorrectOut(), prob_utils_Median.submitted_val));
      result.score = r.first;
      result.pass = r.second;
      result.sub = true;
    }
    return result;
  };

  // Setup how evaluation on world test should work.
  EvaluateWorldTest = [this](prog_org_t & prog_org, size_t testID) {
    emp::Ptr<test_org_t> test_org_ptr = prob_Median_world->GetOrgPtr(testID);
    begin_program_test.Trigger(prog_org, test_org_ptr);
    do_program_test.Trigger(prog_org, test_org_ptr);
    end_program_test.Trigger(prog_org, test_org_ptr);
    return CalcProgramResultOnTest(prog_org, *test_org_ptr);
  };

  // How should we validate programs on testing set?
  prob_utils_Median.population_validation_outputs.resize(PROG_POP_SIZE);
  DoTestingSetValidation = [this](prog_org_t & prog_org) { 
    // evaluate program on full testing set; update stats utils with results
    begin_program_eval.Trigger(prog_org);
    stats_util.current_program__validation__test_results.resize(prob_utils_Median.testingset_pop.size());
    stats_util.current_program__validation__total_score = 0;
    stats_util.current_program__validation__total_passes = 0;
    stats_util.current_program__validation__is_solution = false;
    prob_utils_Median.population_validation_outputs[stats_util.cur_progID].resize(prob_utils_Median.testingset_pop.size());
    // For each test in validation set, evaluate program.
    for (size_t testID = 0; testID < prob_utils_Median.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_Median.testingset_pop[testID];
      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      stats_util.current_program__validation__test_results[testID] = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      stats_util.current_program__validation__total_score += stats_util.current_program__validation__test_results[testID].score;
      stats_util.current_program__validation__total_passes += (size_t)stats_util.current_program__validation__test_results[testID].pass;
      prob_utils_Median.population_validation_outputs[stats_util.cur_progID][testID] = prob_utils_Median.submitted_val;
    }
    stats_util.current_program__validation__is_solution = stats_util.current_program__validation__total_passes == prob_utils_Median.testingset_pop.size();
    end_program_eval.Trigger(prog_org);
  };
  program_stats.get_prog_behavioral_diversity = [this]() { return emp::ShannonEntropy(prob_utils_Median.population_validation_outputs); };
  program_stats.get_prog_unique_behavioral_phenotypes = [this]() { return emp::UniqueCount(prob_utils_Median.population_validation_outputs); };

  // How should we screen for a solution?
  ScreenForSolution = [this](prog_org_t & prog_org) {
    begin_program_eval.Trigger(prog_org);
    for (size_t testID = 0; testID < prob_utils_Median.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_Median.testingset_pop[testID];

      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      
      TestResult result = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      if (!result.pass) {
        end_program_eval.Trigger(prog_org);
        return false;
      }
    }
    end_program_eval.Trigger(prog_org);
    return true;
  };

  // Tell the experiment how to get test phenotypes.
  GetTestPhenotype = [this](size_t testID) -> test_org_phen_t & {
    emp_assert(prob_Median_world->IsOccupied(testID));
    return prob_Median_world->GetOrg(testID).GetPhenotype();
  };

  // Setup how test world updates.
  SetupTestCaseWorldUpdate(prob_Median_world);

  // todo - test that RANDOM is actually making random things every update..

  // Setup how test cases mutate.
  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::RANDOM) {
    std::cout << "RANDOM training mode detected, configuring mutation function to RANDOMIZE organisms." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Randomize organism genome on mutate.
      prob_Median_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        test_org.GetGenome() = GenRandomTestInput_Median(*random, {PROB_MEDIAN__MIN_NUM, PROB_MEDIAN__MAX_NUM}); 
        return 1;
      });
    };
  } else {
    std::cout << "Non-RANDOM training mode detected, configuring mutation function normally." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_Median.MIN_NUM = PROB_MEDIAN__MIN_NUM;
      prob_utils_Median.MAX_NUM = PROB_MEDIAN__MAX_NUM;
      prob_utils_Median.PER_NUM_COPY_RATE = PROB_MEDIAN__MUTATION__PER_NUM_COPY_RATE;
      prob_utils_Median.PER_NUM_SUB_RATE = PROB_MEDIAN__MUTATION__PER_NUM_SUB_RATE;
      prob_utils_Median.PER_NUM_SWAP_RATE = PROB_MEDIAN__MUTATION__PER_NUM_SWAP_RATE;
      // (2) Hook mutator up to world.
      prob_Median_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_Median.Mutate(rnd, test_org.GetGenome());
      });
    };
  }

  // Setup test case fitness function.
  SetupTestFitFun = [this]() {
    prob_Median_world->SetFitFun([](test_org_t & test_org) {
      return (double)test_org.GetPhenotype().num_fails;
    });
  };

  // Tell experiment how to configure hardware inputs when running program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org, emp::Ptr<TestOrg_Base> test_org_base_ptr) {
    // Reset eval stuff
    // Set current test org.
    prob_utils_Median.cur_eval_test_org = test_org_base_ptr.Cast<test_org_t>(); // currently only place need testID for this?
    prob_utils_Median.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 3);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      // Grab some useful references.
      Problem_Median_input_t & input = prob_utils_Median.cur_eval_test_org->GetGenome(); // std::pair<int, double>
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      // Set hardware input.
      wmem.Set(0, input[0]);
      wmem.Set(1, input[1]);
      wmem.Set(2, input[2]);
    }
  });

  // Tell experiment how to snapshot test population.
  SnapshotTests = [this]() {
    std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate());
    mkdir(snapshot_dir.c_str(), ACCESSPERMS);
    
    emp::DataFile file(snapshot_dir + "/test_pop_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");
    // Test file contents:
    // - test id
    std::function<size_t(void)> get_test_id = [this]() { return stats_util.cur_testID; };
    file.AddFun(get_test_id, "test_id");

    // - test fitness
    std::function<double(void)> get_test_fitness = [this]() { return prob_Median_world->CalcFitnessID(stats_util.cur_testID); };
    file.AddFun(get_test_fitness, "fitness");

    // - num passes
    std::function<size_t(void)> get_test_num_passes = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_passes; };
    file.AddFun(get_test_num_passes, "num_passes");

    // - num fails
    std::function<size_t(void)> get_test_num_fails = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_fails; };
    file.AddFun(get_test_num_fails, "num_fails");

    std::function<size_t(void)> get_num_tested = [this]() { return GetTestPhenotype(stats_util.cur_testID).test_passes.size(); };
    file.AddFun(get_num_tested, "num_programs_tested_against");

    // - test scores by program
    std::function<std::string(void)> get_passes_by_program = [this]() {
      std::string scores = "\"[";
      test_org_phen_t & phen = GetTestPhenotype(stats_util.cur_testID);
      for (size_t i = 0; i < phen.test_passes.size(); ++i) {
        if (i) scores += ",";
        scores += emp::to_string(phen.test_passes[i]);
      }
      scores += "]\"";
      return scores;
    };
    file.AddFun(get_passes_by_program, "passes_by_program");

    // - test
    std::function<std::string(void)> get_test = [this]() {
      std::ostringstream stream;
      stream << "\"";
      prob_Median_world->GetOrg(stats_util.cur_testID).Print(stream);
      stream << "\"";
      return stream.str();
    };
    file.AddFun(get_test, "test");

    file.PrintHeaderKeys();

    // Loop over tests, snapshotting each.
    for (stats_util.cur_testID = 0; stats_util.cur_testID < prob_Median_world->GetSize(); ++stats_util.cur_testID) {
      if (!prob_Median_world->IsOccupied(stats_util.cur_testID)) continue;
      file.Update();
    }
  };

  // Add default instructions to instruction set.
  AddDefaultInstructions({"Add",
                          "Sub",
                          "Mult",
                          "Div",
                          "Mod",
                          "TestNumEqu",
                          "TestNumNEqu",
                          "TestNumLess",
                          "TestNumLessTEqu",
                          "TestNumGreater",
                          "TestNumGreaterTEqu",
                          "Floor",
                          "Not",
                          "Inc",
                          "Dec",
                          "CopyMem",
                          "SwapMem",
                          "Input",
                          "Output",
                          "CommitGlobal",
                          "PullGlobal",
                          "TestMemEqu",
                          "TestMemNEqu",
                          "If",
                          "IfNot",
                          "While",
                          "Countdown",
                          "Foreach",
                          "Close",
                          "Break",
                          "Call",
                          "Routine",
                          "Return",
                          "ModuleDef",
                          "MakeVector",
                          "VecGet",
                          "VecSet",
                          "VecLen",
                          "VecAppend",
                          "VecPop",
                          "VecRemove",
                          "VecReplaceAll",
                          "VecIndexOf",
                          "VecOccurrencesOf",
                          "VecReverse",
                          "VecSwapIfLess",
                          "VecGetFront",
                          "VecGetBack",
                          "IsNum",
                          "IsVec"
  });

  // -- Custom Instructions --
  inst_lib->AddInst("LoadNum1", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadNum1_Median(hw, inst);
  }, 1);

  inst_lib->AddInst("LoadNum2", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadNum2_Median(hw, inst);
  }, 1);

  inst_lib->AddInst("LoadNum3", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadNum3_Median(hw, inst);
  }, 1);


  inst_lib->AddInst("SubmitNum", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitNum_Median(hw, inst);
  }, 1);

}

void ProgramSynthesisExperiment::SetupProblem_Smallest() { 
  std::cout << "Setting up problem - Smallest" << std::endl;

  // A few useful aliases.
  using test_org_t = TestOrg_Smallest;

  // Load benchmark data for problem.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  std::cout << "Loading training examples." << std::endl;
  prob_utils_Smallest.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  std::cout << "Loading testing examples." << std::endl;
  prob_utils_Smallest.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  std::cout << "Generating testing set population." << std::endl;
  prob_utils_Smallest.GenerateTestingSetPop();
  std::cout << "Loaded training example set size = " << prob_utils_Smallest.GetTrainingSet().GetSize() << std::endl;
  std::cout << "Loaded testing example set size = " << prob_utils_Smallest.GetTestingSet().GetSize() << std::endl;
  std::cout << "Testing set (non-training examples used to evaluate program accuracy) size = " << prob_utils_Smallest.testingset_pop.size() << std::endl;

  // Setup the world
  NewTestCaseWorld(prob_Smallest_world, *random, "Smallest world");

  // Configure how the population should be initialized
  SetupTestCasePop_Init(prob_Smallest_world,
                        prob_utils_Smallest.training_set,
                        [this]() { return GenRandomTestInput_Smallest(*random, 
                                                                      {PROB_SMALLEST__MIN_NUM, PROB_SMALLEST__MAX_NUM}); 
                                  }
                        );
  end_setup_sig.AddAction([this]() { std::cout << "TestCase world size= " << prob_Smallest_world->GetSize() << std::endl; });

  // Tell the world to calculate the correct test output (given input) on placement.
  prob_Smallest_world->OnPlacement([this](size_t pos) { prob_Smallest_world->GetOrg(pos).CalcOut(); });

  // How are program results calculated on a test?
  CalcProgramResultOnTest = [this](prog_org_t & prog_org, TestOrg_Base & test_org_base) {
    test_org_t & test_org = static_cast<test_org_t&>(test_org_base);
    TestResult result;
    if (!prob_utils_Smallest.submitted) {
      result.score = 0;
      result.pass = false;
      result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_Smallest.CalcScorePassFail(test_org.GetCorrectOut(), prob_utils_Smallest.submitted_val));
      result.score = r.first;
      result.pass = r.second;
      result.sub = true;
    }
    return result;
  };

  // Setup how evaluation on world test should work.
  EvaluateWorldTest = [this](prog_org_t & prog_org, size_t testID) {
    emp::Ptr<test_org_t> test_org_ptr = prob_Smallest_world->GetOrgPtr(testID);
    begin_program_test.Trigger(prog_org, test_org_ptr);
    do_program_test.Trigger(prog_org, test_org_ptr);
    end_program_test.Trigger(prog_org, test_org_ptr);
    return CalcProgramResultOnTest(prog_org, *test_org_ptr);
  };

  // How should we validate programs on testing set?
  prob_utils_Smallest.population_validation_outputs.resize(PROG_POP_SIZE);
  DoTestingSetValidation = [this](prog_org_t & prog_org) { 
    // evaluate program on full testing set; update stats utils with results
    begin_program_eval.Trigger(prog_org);
    stats_util.current_program__validation__test_results.resize(prob_utils_Smallest.testingset_pop.size());
    stats_util.current_program__validation__total_score = 0;
    stats_util.current_program__validation__total_passes = 0;
    stats_util.current_program__validation__is_solution = false;
    prob_utils_Smallest.population_validation_outputs[stats_util.cur_progID].resize(prob_utils_Smallest.testingset_pop.size());
    // For each test in validation set, evaluate program.
    for (size_t testID = 0; testID < prob_utils_Smallest.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_Smallest.testingset_pop[testID];
      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      stats_util.current_program__validation__test_results[testID] = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      stats_util.current_program__validation__total_score += stats_util.current_program__validation__test_results[testID].score;
      stats_util.current_program__validation__total_passes += (size_t)stats_util.current_program__validation__test_results[testID].pass;
      prob_utils_Smallest.population_validation_outputs[stats_util.cur_progID][testID] = prob_utils_Smallest.submitted_val;
    }
    stats_util.current_program__validation__is_solution = stats_util.current_program__validation__total_passes == prob_utils_Smallest.testingset_pop.size();
    end_program_eval.Trigger(prog_org);
  };
  program_stats.get_prog_behavioral_diversity = [this]() { return emp::ShannonEntropy(prob_utils_Smallest.population_validation_outputs); };
  program_stats.get_prog_unique_behavioral_phenotypes = [this]() { return emp::UniqueCount(prob_utils_Smallest.population_validation_outputs); };

  // How should we screen for a solution?
  ScreenForSolution = [this](prog_org_t & prog_org) {
    begin_program_eval.Trigger(prog_org);
    for (size_t testID = 0; testID < prob_utils_Smallest.testingset_pop.size(); ++testID) {
      stats_util.cur_testID = testID;
      emp::Ptr<test_org_t> test_org_ptr = prob_utils_Smallest.testingset_pop[testID];

      begin_program_test.Trigger(prog_org, test_org_ptr);
      do_program_test.Trigger(prog_org, test_org_ptr);
      end_program_test.Trigger(prog_org, test_org_ptr);
      
      TestResult result = CalcProgramResultOnTest(prog_org, *test_org_ptr);
      if (!result.pass) {
        end_program_eval.Trigger(prog_org);
        return false;
      }
    }
    end_program_eval.Trigger(prog_org);
    return true;
  };

  // Tell the experiment how to get test phenotypes.
  GetTestPhenotype = [this](size_t testID) -> test_org_phen_t & {
    emp_assert(prob_Smallest_world->IsOccupied(testID));
    return prob_Smallest_world->GetOrg(testID).GetPhenotype();
  };

  // Setup how test world updates.
  SetupTestCaseWorldUpdate(prob_Smallest_world);

  // todo - test that RANDOM is actually making random things every update..

  // Setup how test cases mutate.
  if (TRAINING_EXAMPLE_MODE == (size_t)TRAINING_EXAMPLE_MODE_TYPE::RANDOM) {
    std::cout << "RANDOM training mode detected, configuring mutation function to RANDOMIZE organisms." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Randomize organism genome on mutate.
      prob_Smallest_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        test_org.GetGenome() = GenRandomTestInput_Smallest(*random, {PROB_SMALLEST__MIN_NUM, PROB_SMALLEST__MAX_NUM}); 
        return 1;
      });
    };
  } else {
    std::cout << "Non-RANDOM training mode detected, configuring mutation function normally." << std::endl;
    SetupTestMutation = [this]() {
      // (1) Configure mutator.
      prob_utils_Smallest.MIN_NUM = PROB_SMALLEST__MIN_NUM;
      prob_utils_Smallest.MAX_NUM = PROB_SMALLEST__MAX_NUM;
      prob_utils_Smallest.PER_NUM_SUB_RATE = PROB_SMALLEST__MUTATION__PER_NUM_SUB_RATE;
      prob_utils_Smallest.PER_NUM_SWAP_RATE = PROB_SMALLEST__MUTATION__PER_NUM_SWAP_RATE;
      // (2) Hook mutator up to world.
      prob_Smallest_world->SetMutFun([this](test_org_t & test_org, emp::Random & rnd) {
        return prob_utils_Smallest.Mutate(rnd, test_org.GetGenome());
      });
    };
  }

  // Setup test case fitness function.
  SetupTestFitFun = [this]() {
    prob_Smallest_world->SetFitFun([](test_org_t & test_org) {
      return (double)test_org.GetPhenotype().num_fails;
    });
  }; 

  // Tell experiment how to configure hardware inputs when running program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org, emp::Ptr<TestOrg_Base> test_org_base_ptr) {
    // Reset eval stuff
    // Set current test org.
    prob_utils_Smallest.cur_eval_test_org = test_org_base_ptr.Cast<test_org_t>(); // currently only place need testID for this?
    prob_utils_Smallest.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 4);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      // Grab some useful references.
      Problem_Smallest_input_t & input = prob_utils_Smallest.cur_eval_test_org->GetGenome(); // std::pair<int, double>
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      // Set hardware input.
      wmem.Set(0, input[0]);
      wmem.Set(1, input[1]);
      wmem.Set(2, input[2]);
      wmem.Set(3, input[3]);
    }
  });

  // Tell experiment how to snapshot test population.
  SnapshotTests = [this]() {
    std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate());
    mkdir(snapshot_dir.c_str(), ACCESSPERMS);
    
    emp::DataFile file(snapshot_dir + "/test_pop_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");
    // Test file contents:
    // - test id
    std::function<size_t(void)> get_test_id = [this]() { return stats_util.cur_testID; };
    file.AddFun(get_test_id, "test_id");

    // - test fitness
    std::function<double(void)> get_test_fitness = [this]() { return prob_Smallest_world->CalcFitnessID(stats_util.cur_testID); };
    file.AddFun(get_test_fitness, "fitness");

    // - num passes
    std::function<size_t(void)> get_test_num_passes = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_passes; };
    file.AddFun(get_test_num_passes, "num_passes");

    // - num fails
    std::function<size_t(void)> get_test_num_fails = [this]() { return GetTestPhenotype(stats_util.cur_testID).num_fails; };
    file.AddFun(get_test_num_fails, "num_fails");

    std::function<size_t(void)> get_num_tested = [this]() { return GetTestPhenotype(stats_util.cur_testID).test_passes.size(); };
    file.AddFun(get_num_tested, "num_programs_tested_against");

    // - test scores by program
    std::function<std::string(void)> get_passes_by_program = [this]() {
      std::string scores = "\"[";
      test_org_phen_t & phen = GetTestPhenotype(stats_util.cur_testID);
      for (size_t i = 0; i < phen.test_passes.size(); ++i) {
        if (i) scores += ",";
        scores += emp::to_string(phen.test_passes[i]);
      }
      scores += "]\"";
      return scores;
    };
    file.AddFun(get_passes_by_program, "passes_by_program");

    // - test
    std::function<std::string(void)> get_test = [this]() {
      std::ostringstream stream;
      stream << "\"";
      prob_Smallest_world->GetOrg(stats_util.cur_testID).Print(stream);
      stream << "\"";
      return stream.str();
    };
    file.AddFun(get_test, "test", "");

    file.PrintHeaderKeys();

    // Loop over tests, snapshotting each.
    for (stats_util.cur_testID = 0; stats_util.cur_testID < prob_Smallest_world->GetSize(); ++stats_util.cur_testID) {
      if (!prob_Smallest_world->IsOccupied(stats_util.cur_testID)) continue;
      file.Update();
    }
  };

  // Add default instructions to instruction set.
  AddDefaultInstructions({"Add",
                          "Sub",
                          "Mult",
                          "Div",
                          "Mod",
                          "TestNumEqu",
                          "TestNumNEqu",
                          "TestNumLess",
                          "TestNumLessTEqu",
                          "TestNumGreater",
                          "TestNumGreaterTEqu",
                          "Floor",
                          "Not",
                          "Inc",
                          "Dec",
                          "CopyMem",
                          "SwapMem",
                          "Input",
                          "Output",
                          "CommitGlobal",
                          "PullGlobal",
                          "TestMemEqu",
                          "TestMemNEqu",
                          "If",
                          "IfNot",
                          "While",
                          "Countdown",
                          "Foreach",
                          "Close",
                          "Break",
                          "Call",
                          "Routine",
                          "Return",
                          "ModuleDef",
                          "MakeVector",
                          "VecGet",
                          "VecSet",
                          "VecLen",
                          "VecAppend",
                          "VecPop",
                          "VecRemove",
                          "VecReplaceAll",
                          "VecIndexOf",
                          "VecOccurrencesOf",
                          "VecReverse",
                          "VecSwapIfLess",
                          "VecGetFront",
                          "VecGetBack",
                          "IsNum",
                          "IsVec"
  });

  // -- Custom Instructions --
  inst_lib->AddInst("LoadNum1", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadNum1_Smallest(hw, inst);
  }, 1);

  inst_lib->AddInst("LoadNum2", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadNum2_Smallest(hw, inst);
  }, 1);

  inst_lib->AddInst("LoadNum3", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadNum3_Smallest(hw, inst);
  }, 1);

  inst_lib->AddInst("LoadNum4", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_LoadNum4_Smallest(hw, inst);
  }, 1);

  inst_lib->AddInst("SubmitNum", [this](hardware_t & hw, const inst_t & inst) {
    this->Inst_SubmitNum_Smallest(hw, inst);
  }, 1);

}

void ProgramSynthesisExperiment::SetupProblem_Syllables() { 
  std::cout << "Problem setup not yet implemented... Exiting." << std::endl;
  exit(-1); 
}

// ================= PROBLEM-SPECIFIC INSTRUCTION IMPLEMENTATIONS ======================

void ProgramSynthesisExperiment::Inst_LoadInt_NumberIO(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_NumberIO.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_NumberIO.cur_eval_test_org->GetTestIntegerInput());
}

void ProgramSynthesisExperiment::Inst_LoadDouble_NumberIO(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_NumberIO.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_NumberIO.cur_eval_test_org->GetTestDoubleInput());
}

void ProgramSynthesisExperiment::Inst_SubmitNum_NumberIO(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_NumberIO.cur_eval_test_org != nullptr);
  prob_utils_NumberIO.Submit(wmem.AccessVal(posA).GetNum());
}

// ----- SmallOrLarge -----
void ProgramSynthesisExperiment::Inst_LoadInt_SmallOrLarge(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_SmallOrLarge.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_SmallOrLarge.cur_eval_test_org->GetGenome());
}

void ProgramSynthesisExperiment::Inst_SubmitSmall_SmallOrLarge(hardware_t & hw, const inst_t & inst) {
  prob_utils_SmallOrLarge.Submit("small");
}

void ProgramSynthesisExperiment::Inst_SubmitLarge_SmallOrLarge(hardware_t & hw, const inst_t & inst) {
  prob_utils_SmallOrLarge.Submit("large");
}

void ProgramSynthesisExperiment::Inst_SubmitNone_SmallOrLarge(hardware_t & hw, const inst_t & inst) {
  prob_utils_SmallOrLarge.Submit("");
}

// ----- ForLoopIndex -----
void ProgramSynthesisExperiment::Inst_LoadStart_ForLoopIndex(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_ForLoopIndex.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_ForLoopIndex.cur_eval_test_org->GetGenome()[0]);
}

void ProgramSynthesisExperiment::Inst_LoadEnd_ForLoopIndex(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_ForLoopIndex.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_ForLoopIndex.cur_eval_test_org->GetGenome()[1]);
}

void ProgramSynthesisExperiment::Inst_LoadStep_ForLoopIndex(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_ForLoopIndex.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_ForLoopIndex.cur_eval_test_org->GetGenome()[2]);
}

void ProgramSynthesisExperiment::Inst_SubmitNum_ForLoopIndex(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_ForLoopIndex.cur_eval_test_org != nullptr);
  prob_utils_ForLoopIndex.Submit((int)wmem.AccessVal(posA).GetNum());
}

// ----- CompareStringLengths -----
void ProgramSynthesisExperiment::Inst_LoadStr1_CompareStringLengths(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_CompareStringLengths.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_CompareStringLengths.cur_eval_test_org->GetGenome()[0]);
}

void ProgramSynthesisExperiment::Inst_LoadStr2_CompareStringLengths(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_CompareStringLengths.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_CompareStringLengths.cur_eval_test_org->GetGenome()[1]);
}

void ProgramSynthesisExperiment::Inst_LoadStr3_CompareStringLengths(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_CompareStringLengths.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_CompareStringLengths.cur_eval_test_org->GetGenome()[2]);
}

void ProgramSynthesisExperiment::Inst_SubmitTrue_CompareStringLengths(hardware_t & hw, const inst_t & inst) {
  emp_assert(prob_utils_CompareStringLengths.cur_eval_test_org != nullptr);
  prob_utils_CompareStringLengths.Submit(true);
}

void ProgramSynthesisExperiment::Inst_SubmitFalse_CompareStringLengths(hardware_t & hw, const inst_t & inst) {
  emp_assert(prob_utils_CompareStringLengths.cur_eval_test_org != nullptr);
  prob_utils_CompareStringLengths.Submit(false);
}

void ProgramSynthesisExperiment::Inst_SubmitVal_CompareStringLengths(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_CompareStringLengths.cur_eval_test_org != nullptr);
  prob_utils_CompareStringLengths.Submit((bool)wmem.AccessVal(posA).GetNum());
}

// ----- CollatzNumbers -----
void ProgramSynthesisExperiment::Inst_LoadNum_CollatzNumbers(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_CollatzNumbers.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_CollatzNumbers.cur_eval_test_org->GetGenome());
}

void ProgramSynthesisExperiment::Inst_SubmitNum_CollatzNumbers(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_CollatzNumbers.cur_eval_test_org != nullptr);
  prob_utils_CollatzNumbers.Submit((int)wmem.AccessVal(posA).GetNum());
}

// ----- StringLengthsBackwards -----
void ProgramSynthesisExperiment::Inst_LoadStrVec_StringLengthsBackwards(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_StringLengthsBackwards.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_StringLengthsBackwards.cur_eval_test_org->GetGenome());
}

void ProgramSynthesisExperiment::Inst_SubmitVal_StringLengthsBackwards(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_StringLengthsBackwards.cur_eval_test_org != nullptr);
  prob_utils_StringLengthsBackwards.Submit((size_t)wmem.AccessVal(posA).GetNum());
}

// ----- LastIndexOfZero -----
void ProgramSynthesisExperiment::Inst_LoadVec_LastIndexOfZero(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_LastIndexOfZero.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_LastIndexOfZero.cur_eval_test_org->GetGenome());
}

void ProgramSynthesisExperiment::Inst_SubmitNum_LastIndexOfZero(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_LastIndexOfZero.cur_eval_test_org != nullptr);
  prob_utils_LastIndexOfZero.Submit((int)wmem.AccessVal(posA).GetNum());
}

// ---- CountOdds -----
void ProgramSynthesisExperiment::Inst_LoadVec_CountOdds(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_CountOdds.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_CountOdds.cur_eval_test_org->GetGenome());
}

void ProgramSynthesisExperiment::Inst_SubmitNum_CountOdds(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_CountOdds.cur_eval_test_org != nullptr);
  prob_utils_CountOdds.Submit((int)wmem.AccessVal(posA).GetNum());
}

// ----- Mirror Image -----
void ProgramSynthesisExperiment::Inst_LoadVec1_MirrorImage(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_MirrorImage.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_MirrorImage.cur_eval_test_org->GetGenome()[0]);
}

void ProgramSynthesisExperiment::Inst_LoadVec2_MirrorImage(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_MirrorImage.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_MirrorImage.cur_eval_test_org->GetGenome()[1]);
}

void ProgramSynthesisExperiment::Inst_SubmitVal_MirrorImage(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_MirrorImage.cur_eval_test_org != nullptr);
  prob_utils_MirrorImage.Submit((bool)wmem.AccessVal(posA).GetNum());
}

void ProgramSynthesisExperiment::Inst_SubmitTrue_MirrorImage(hardware_t & hw, const inst_t & inst) {
  prob_utils_MirrorImage.Submit(true);
}

void ProgramSynthesisExperiment::Inst_SubmitFalse_MirrorImage(hardware_t & hw, const inst_t & inst) {
  prob_utils_MirrorImage.Submit(false);
}

// ----- VectorsSummed -----
void ProgramSynthesisExperiment::Inst_LoadVec1_VectorsSummed(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_VectorsSummed.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_VectorsSummed.cur_eval_test_org->GetGenome()[0]);
}

void ProgramSynthesisExperiment::Inst_LoadVec2_VectorsSummed(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_VectorsSummed.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_VectorsSummed.cur_eval_test_org->GetGenome()[1]);
}

void ProgramSynthesisExperiment::Inst_SubmitVec_VectorsSummed(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::VEC);
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_VectorsSummed.cur_eval_test_org != nullptr);
  const emp::vector<hardware_t::MemoryValue> & vec = wmem.AccessVec(posA);
  emp::vector<int> output;
  for (size_t i = 0; i < vec.size(); ++i) {
    if (vec[i].GetType() == hardware_t::MemoryValue::MemoryType::NUM) {
      output.emplace_back((int)vec[i].GetNum());
    }
  }
  prob_utils_VectorsSummed.Submit(output); 
}

// ----- Sum Of Squares -----
void ProgramSynthesisExperiment::Inst_LoadNum_SumOfSquares(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_SumOfSquares.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_SumOfSquares.cur_eval_test_org->GetGenome());
}

void ProgramSynthesisExperiment::Inst_SubmitNum_SumOfSquares(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_SumOfSquares.cur_eval_test_org != nullptr);
  prob_utils_SumOfSquares.Submit((int)wmem.AccessVal(posA).GetNum());
}

// ----- VectorAverage ------
void ProgramSynthesisExperiment::Inst_LoadVec_VectorAverage(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_VectorAverage.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_VectorAverage.cur_eval_test_org->GetGenome());
}

void ProgramSynthesisExperiment::Inst_SubmitNum_VectorAverage(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_VectorAverage.cur_eval_test_org != nullptr);
  prob_utils_VectorAverage.Submit((double)wmem.AccessVal(posA).GetNum());
}


// ----- Median -----
void ProgramSynthesisExperiment::Inst_LoadNum1_Median(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_Median.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_Median.cur_eval_test_org->GetGenome()[0]);
}

void ProgramSynthesisExperiment::Inst_LoadNum2_Median(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_Median.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_Median.cur_eval_test_org->GetGenome()[1]);
}

void ProgramSynthesisExperiment::Inst_LoadNum3_Median(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_Median.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_Median.cur_eval_test_org->GetGenome()[2]);
}

void ProgramSynthesisExperiment::Inst_SubmitNum_Median(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_Median.cur_eval_test_org != nullptr);
  prob_utils_Median.Submit((int)wmem.AccessVal(posA).GetNum());
}

// ----- Smallest -----
void ProgramSynthesisExperiment::Inst_LoadNum1_Smallest(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_Smallest.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_Smallest.cur_eval_test_org->GetGenome()[0]);
}

void ProgramSynthesisExperiment::Inst_LoadNum2_Smallest(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_Smallest.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_Smallest.cur_eval_test_org->GetGenome()[1]);
}

void ProgramSynthesisExperiment::Inst_LoadNum3_Smallest(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_Smallest.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_Smallest.cur_eval_test_org->GetGenome()[2]);
}

void ProgramSynthesisExperiment::Inst_LoadNum4_Smallest(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_Smallest.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_Smallest.cur_eval_test_org->GetGenome()[3]);
}

void ProgramSynthesisExperiment::Inst_SubmitNum_Smallest(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_Smallest.cur_eval_test_org != nullptr);
  prob_utils_Smallest.Submit((int)wmem.AccessVal(posA).GetNum());
}

// --- Grade ---
void ProgramSynthesisExperiment::Inst_LoadThreshA_Grade(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_Grade.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_Grade.cur_eval_test_org->GetGenome()[0]);
}

void ProgramSynthesisExperiment::Inst_LoadThreshB_Grade(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_Grade.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_Grade.cur_eval_test_org->GetGenome()[1]);
}

void ProgramSynthesisExperiment::Inst_LoadThreshC_Grade(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_Grade.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_Grade.cur_eval_test_org->GetGenome()[2]);
}

void ProgramSynthesisExperiment::Inst_LoadThreshD_Grade(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_Grade.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_Grade.cur_eval_test_org->GetGenome()[3]);
}

void ProgramSynthesisExperiment::Inst_LoadGrade_Grade(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  emp_assert(prob_utils_Grade.cur_eval_test_org != nullptr);
  wmem.Set(posA, prob_utils_Grade.cur_eval_test_org->GetGenome()[4]);
}

void ProgramSynthesisExperiment::Inst_SubmitA_Grade(hardware_t & hw, const inst_t & inst) {
  prob_utils_Grade.Submit("A");
}

void ProgramSynthesisExperiment::Inst_SubmitB_Grade(hardware_t & hw, const inst_t & inst) {
  prob_utils_Grade.Submit("B");
}

void ProgramSynthesisExperiment::Inst_SubmitC_Grade(hardware_t & hw, const inst_t & inst) {
  prob_utils_Grade.Submit("C");
}

void ProgramSynthesisExperiment::Inst_SubmitD_Grade(hardware_t & hw, const inst_t & inst) {
  prob_utils_Grade.Submit("D");
}

void ProgramSynthesisExperiment::Inst_SubmitF_Grade(hardware_t & hw, const inst_t & inst) {
  prob_utils_Grade.Submit("F");
}

// ------------------------


#endif