#ifndef PROGRAM_SYNTHESIS_CONFIG_H
#define PROGRAM_SYNTHESIS_CONFIG_H

#include "config/config.h"

EMP_BUILD_CONFIG(ProgramSynthesisConfig, 
  GROUP(DEFAULT_GROUP, "General settings"),
  VALUE(SEED, int, 0, "Random number seed (-1 for based on time)"),
  VALUE(GENERATIONS, size_t, 10000, "How many generations should we run for?"),
  VALUE(PROG_POP_SIZE, size_t, 1024, "Population size for programs"),
  VALUE(TEST_POP_SIZE, size_t, 1024, "Population size for tests - how many tests exist at once?"),
  VALUE(EVALUATION_MODE, size_t, 0, "How are programs evaluated? \n0: cohorts \n1: full (on all tests) \n2: program only cohorts \n3: TEST_DOWNSAMPLING treatment"),
  VALUE(PROG_COHORT_SIZE, size_t, 32, "How big should random cohorts be (only relevant when evaluating with random cohort method)?"),
  VALUE(TEST_COHORT_SIZE, size_t, 32, "How big should random cohorts be (only relevant when evaluating with random cohort method)?"),
  VALUE(TRAINING_EXAMPLE_MODE, size_t, 0, "How do training examples change over time? \n0: co-evolution \n1: static \n2: random \n3: Static-gen \n4: STATIC_COEVO "),
  VALUE(PROBLEM, std::string, "number-io", "Which problem to use?"),
  VALUE(BENCHMARK_DATA_DIR, std::string, "../data/prog-synth-examples", "Location to look for problem test case data."),

  GROUP(SELECTION_GROUP, "Settings specific to selection (both tests and programs)."),
  VALUE(PROG_SELECTION_MODE, size_t, 1, "How are selected? \n0: LEXICASE \n1: COHORT_LEXICASE \n2: TOURNAMENT \n3: DRIFT \n4: PROG_ONLY_COHORT_LEXICASE \n5: TEST_DOWNSAMPLING treatment"),
  VALUE(TEST_SELECTION_MODE, size_t, 1, "How are selected? \n0: LEXICASE \n1: COHORT_LEXICASE \n2: TOURNAMENT \n3: DRIFT"),
  VALUE(PROG_LEXICASE_MAX_FUNS, size_t, 0, "Max number of functions to use for program lexicase selection. (0 to use all)"),
  VALUE(PROG_COHORTLEXICASE_MAX_FUNS, size_t, 0, "Max number of functions to use for program cohort lexicase selection. (0 to use all)"),
  VALUE(TEST_LEXICASE_MAX_FUNS, size_t, 0, "Max number of functions to use for test lexicase selection. (0 to use all)"),
  VALUE(TEST_COHORTLEXICASE_MAX_FUNS, size_t, 0, "Max number of functions to use for test cohort lexicase selection. (0 to use all)"),
  VALUE(PROG_TOURNAMENT_SIZE, size_t, 4, "How big should tournaments be during program tournament selection?"),
  VALUE(TEST_TOURNAMENT_SIZE, size_t, 4, "How big should tournaments be during test tournament selection?"),
  VALUE(DISCRIMINATORY_LEXICASE_TESTS, bool, false, "Should we use discriminatory test cases for lexicase selection?"),

  GROUP(PROGRAM_GROUP, "General settings specific to programs."),
  VALUE(MIN_PROG_SIZE, size_t, 1, "Minimum program size"),
  VALUE(MAX_PROG_SIZE, size_t, 128, "Maximum program size"),
  VALUE(PROG_EVAL_TIME, size_t, 256, "How many clock cycles should we give a program during a test?"),
  VALUE(PROG_MUT__PER_BIT_FLIP, double, 0.001, "Program per-bit flip rate."),
  VALUE(PROG_MUT__PER_INST_SUB, double, 0.005, "Program per-instruction substitution mutation rate."),
  VALUE(PROG_MUT__PER_INST_INS, double, 0.005, "Program per-instruction insertion mutation rate."),
  VALUE(PROG_MUT__PER_INST_DEL, double, 0.005, "Program per-instruction deletion mutation rate."),
  VALUE(PROG_MUT__PER_PROG_SLIP, double, 0.05, "Program per-program slip mutation rate."),
  VALUE(PROG_MUT__PER_MOD_DUP, double, 0.05, "Program per-module whole-module duplication rate."),
  VALUE(PROG_MUT__PER_MOD_DEL, double, 0.05, "Program per-module whole-module deletion rate."),

  GROUP(HARDWARE_GROUP, "Settings specific to TagLGP virtual hardware"),
  VALUE(MIN_TAG_SPECIFICITY, double, 0.0, "What is the minimum tag similarity required for a tag to successfully reference another tag?"),
  VALUE(MAX_CALL_DEPTH, size_t, 128, "Maximum depth of hardware's call stack."),

  GROUP(PROB_NUMBER_IO_GROUP, "Settings specific to NumberIO problem."),
  VALUE(PROB_NUMBER_IO__DOUBLE_MIN, double, -100.0, "Min value for input double."),
  VALUE(PROB_NUMBER_IO__DOUBLE_MAX, double, 100.0, "Max value for input double."),
  VALUE(PROB_NUMBER_IO__INT_MIN, int, -100, "Min value for input int."),
  VALUE(PROB_NUMBER_IO__INT_MAX, int, 100, "Max value for input int."),
  VALUE(PROB_NUMBER_IO__MUTATION__PER_INT_RATE, double, 0.1, "Per-integer mutation rate."),
  VALUE(PROB_NUMBER_IO__MUTATION__PER_DOUBLE_RATE, double, 0.1, "Per-double mutation rate."),

  GROUP(PROB_SMALL_OR_LARGE_GROUP, "Settings specific to the small or large problem."),
  VALUE(PROB_SMALL_OR_LARGE__INT_MIN, int, -10000, "Min value for input int."),
  VALUE(PROB_SMALL_OR_LARGE__INT_MAX, int, 10000, "Max value for input int."),
  VALUE(PROB_SMALL_OR_LARGE__MUTATION__PER_INT_RATE, double, 0.1, "Per-integer mutation rate."),

  GROUP(PROB_FOR_LOOP_INDEX_GROUP, "Settings specific to the for loop index problem."),
  VALUE(PROB_FOR_LOOP_INDEX__START_END_MIN, int, -500, ". value for ."),
  VALUE(PROB_FOR_LOOP_INDEX__START_END_MAX, int, 500, ". value for ."),
  VALUE(PROB_FOR_LOOP_INDEX__STEP_MIN, int, 1, ". value for ."),
  VALUE(PROB_FOR_LOOP_INDEX__STEP_MAX, int, 10, ". value for ."),
  VALUE(PROB_FOR_LOOP_INDEX__MUTATION__MUT_RATE, double, 0.1, "Per-integer mutation rate."),
  VALUE(PROB_FOR_LOOP_INDEX__PROMISE_MULTISTEP_TESTCASES, bool, false, "When randomly generating test cases, do we promise to generate test cases whose output is a sequence with length > 1?"),

  GROUP(PROB_COMPARE_STRING_LENGTHS_GROUP, "Settings specific to the compare string lengths problem."),
  VALUE(PROB_COMPARE_STRING_LENGTHS__MIN_STR_LEN, size_t, 0, "."),
  VALUE(PROB_COMPARE_STRING_LENGTHS__MAX_STR_LEN, size_t, 49, "."),
  VALUE(PROB_COMPARE_STRING_LENGTHS__PER_SITE_INS_RATE, double, 0.1, "."),
  VALUE(PROB_COMPARE_STRING_LENGTHS__PER_SITE_DEL_RATE, double, 0.1, "."),
  VALUE(PROB_COMPARE_STRING_LENGTHS__PER_SITE_SUB_RATE, double, 0.1, "."),
  VALUE(PROB_COMPARE_STRING_LENGTHS__PER_STR_SWAP_RATE, double, 0.1, "."),

  GROUP(PROB_COLLATZ_NUMBERS_GROUP, "Settings specific to the collatz numbers problem."),
  VALUE(PROB_COLLATZ_NUMBERS__MIN_NUM, int, 1, "."),
  VALUE(PROB_COLLATZ_NUMBERS__MAX_NUM, int, 10000, "."),
  VALUE(PROB_COLLATZ_NUMBERS__MUTATION__PER_NUM_SUB_RATE, double, 0.1, "."),

  GROUP(PROB_STRING_LENGTHS_BACKWARDS_GROUP, "Settings specific to the string lengths backwards problem"),
  VALUE(PROB_STRING_LENGTHS_BACKWARDS__MIN_STR_LEN, size_t, 0, "."),
  VALUE(PROB_STRING_LENGTHS_BACKWARDS__MAX_STR_LEN, size_t, 50, "."),
  VALUE(PROB_STRING_LENGTHS_BACKWARDS__MIN_STR_CNT, size_t, 0, "."),
  VALUE(PROB_STRING_LENGTHS_BACKWARDS__MAX_STR_CNT, size_t, 50, "."),
  VALUE(PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_CHAR_INS_RATE, double, 0.05, "."),
  VALUE(PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_CHAR_DEL_RATE, double, 0.05, "."),
  VALUE(PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_CHAR_SUB_RATE, double, 0.05, "."),
  VALUE(PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_STR_SWAP_RATE, double, 0.1, "."),
  VALUE(PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_STR_DUP_RATE, double, 0.1, "."),
  VALUE(PROB_STRING_LENGTHS_BACKWARDS__MUTATION__PER_STR_DEL_RATE, double, 0.1, "."),

  GROUP(PROG_LAST_INDEX_OF_ZERO_GROUP, "Settings specific to the last index of zero problem."),
  VALUE(PROB_LAST_INDEX_OF_ZERO__MIN_VEC_LEN, size_t, 1, "."),
  VALUE(PROB_LAST_INDEX_OF_ZERO__MAX_VEC_LEN, size_t, 50, "."),
  VALUE(PROB_LAST_INDEX_OF_ZERO__MIN_NUM, int, -50, "."),
  VALUE(PROB_LAST_INDEX_OF_ZERO__MAX_NUM, int, 50, "."),
  VALUE(PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_SWAP_RATE, double, 0.05, "."),
  VALUE(PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_DEL_RATE, double, 0.05, "."),
  VALUE(PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_INS_RATE, double, 0.05, "."),
  VALUE(PROB_LAST_INDEX_OF_ZERO__MUTATION__PER_NUM_SUB_RATE, double, 0.05, "."),

  GROUP(PROB_COUNT_ODDS_GROUP, "Settings specific to the count odds problem."),
  VALUE(PROB_COUNT_ODDS__MIN_VEC_LEN, size_t, 0, "."),
  VALUE(PROB_COUNT_ODDS__MAX_VEC_LEN, size_t, 50, "."),
  VALUE(PROB_COUNT_ODDS__MIN_NUM, int, -1000, "."),
  VALUE(PROB_COUNT_ODDS__MAX_NUM, int, 1000, "."),
  VALUE(PROB_COUNT_ODDS__MUTATION__PER_NUM_SWAP_RATE, double, 0.05, "."),
  VALUE(PROB_COUNT_ODDS__MUTATION__PER_NUM_DEL_RATE, double, 0.05, "."),
  VALUE(PROB_COUNT_ODDS__MUTATION__PER_NUM_INS_RATE, double, 0.05, "."),
  VALUE(PROB_COUNT_ODDS__MUTATION__PER_NUM_SUB_RATE, double, 0.05, "."),

  GROUP(PROB_MIRROR_IMAGE_GROUP, "Settings specific to the mirror image problem."),
  VALUE(PROB_MIRROR_IMAGE__MIN_VEC_LEN, size_t, 0, "."),
  VALUE(PROB_MIRROR_IMAGE__MAX_VEC_LEN, size_t, 50, "."),
  VALUE(PROB_MIRROR_IMAGE__MIN_NUM, int, -1000, "."),
  VALUE(PROB_MIRROR_IMAGE__MAX_NUM, int, 1000, "."),
  VALUE(PROB_MIRROR_IMAGE__MUTATION__PER_VEC_RANDOMIZE_VAL_RATE, double, 0.1, "."),
  VALUE(PROB_MIRROR_IMAGE__MUTATION__PER_VEC_MIRROR_RATE, double, 0.1, "."),
  VALUE(PROB_MIRROR_IMAGE__MUTATION__COPY_RATE, double, 0.1, "."),
  VALUE(PROB_MIRROR_IMAGE__MUTATION__INS_RATE, double, 0.1, "."),
  VALUE(PROB_MIRROR_IMAGE__MUTATION__DEL_RATE, double, 0.1, "."),
  VALUE(PROB_MIRROR_IMAGE__MUTATION__PER_VEC_SHUFFLE_RATE, double, 0.1, "."),

  GROUP(PROB_VECTORS_SUMMED_GROUP, "Settings specific to the vectors summed problem."),
  VALUE(PROB_VECTORS_SUMMED__MIN_VEC_LEN, size_t, 0, "."),
  VALUE(PROB_VECTORS_SUMMED__MAX_VEC_LEN, size_t, 50, "."),
  VALUE(PROB_VECTORS_SUMMED__MIN_NUM, int, -1000, "."),
  VALUE(PROB_VECTORS_SUMMED__MAX_NUM, int, 1000, "."),
  VALUE(PROB_VECTORS_SUMMED__MUTATION__PER_NUM_SUB_RATE, double, 0.1, "."),
  VALUE(PROB_VECTORS_SUMMED__MUTATION__COPY_RATE, double, 0.1, "."),
  VALUE(PROB_VECTORS_SUMMED__MUTATION__INS_RATE, double, 0.1, "."),
  VALUE(PROB_VECTORS_SUMMED__MUTATION__DEL_RATE, double, 0.1, "."),

  GROUP(PROB_SUM_OF_SQUARES, "Settings specific to the sum of squares problem."),
  VALUE(PROB_SUM_OF_SQUARES__MIN_NUM, int, 1, "."),
  VALUE(PROB_SUM_OF_SQUARES__MAX_NUM, int, 100, "."),
  VALUE(PROB_SUM_OF_SQUARES__MUTATION__NUM_MUT_RATE, double, 0.2, "."),

  GROUP(PROB_VECTOR_AVERAGE_GROUP, "Settings specific to the vector average problem"),
  VALUE(PROB_VECTOR_AVERAGE__EPSILON, double, 0.00005, "."),
  VALUE(PROB_VECTOR_AVERAGE__MIN_VEC_LEN, size_t, 1, "."),
  VALUE(PROB_VECTOR_AVERAGE__MAX_VEC_LEN, size_t, 50, "."),
  VALUE(PROB_VECTOR_AVERAGE__MIN_NUM, double, -1000.0, "."),
  VALUE(PROB_VECTOR_AVERAGE__MAX_NUM, double, 1000.0, "."),
  VALUE(PROB_VECTOR_AVERAGE__MUTATION__INS_RATE, double, 0.1, "."),
  VALUE(PROB_VECTOR_AVERAGE__MUTATION__DEL_RATE, double, 0.1, "."),
  VALUE(PROB_VECTOR_AVERAGE__MUTATION__SUB_RATE, double, 0.1, "."),

  GROUP(PROB_SMALLEST_GROUP, "Settings specific to the smallest problem"),
  VALUE(PROB_SMALLEST__MIN_NUM, int, -100, "."),
  VALUE(PROB_SMALLEST__MAX_NUM, int, 100, "."),
  VALUE(PROB_SMALLEST__MUTATION__PER_NUM_SUB_RATE, double, 0.1, "."),
  VALUE(PROB_SMALLEST__MUTATION__PER_NUM_SWAP_RATE, double, 0.1, "."),

  GROUP(PROB_GRADE_GROUP, "Settings specific to the grade problem"),
  VALUE(PROB_GRADE__MIN_NUM, int, 0, "."),
  VALUE(PROB_GRADE__MAX_NUM, int, 100, "."),
  VALUE(PROB_GRADE__MUTATION__PER_NUM_RANDOMIZE_RATE, double, 0.1, "."),
  VALUE(PROB_GRADE__MUTATION__PER_NUM_ADJUST_RATE, double, 0.1, "."),

  GROUP(PROB_MEDIAN_GROUP, "Settings specific to the median problem"),
  VALUE(PROB_MEDIAN__MIN_NUM, int, -100, "."),
  VALUE(PROB_MEDIAN__MAX_NUM, int, 100, "."),
  VALUE(PROB_MEDIAN__MUTATION__PER_NUM_COPY_RATE, double, 0.1, "."),
  VALUE(PROB_MEDIAN__MUTATION__PER_NUM_SUB_RATE, double, 0.1, "."),
  VALUE(PROB_MEDIAN__MUTATION__PER_NUM_SWAP_RATE, double, 0.1, "."),

  GROUP(DATA_COLLECTION_GROUP, "Settings specific to data collection."),
  VALUE(DATA_DIRECTORY, std::string, "./output", "Where should we dump output files?"),
  VALUE(SNAPSHOT_INTERVAL, size_t, 1000, "How often should we take population snapshots?"),
  VALUE(SUMMARY_STATS_INTERVAL, size_t, 1000, "How often should we output summary stats?"),
  VALUE(SOLUTION_SCREEN_INTERVAL, size_t, 1000, "How often should we screen entire population for solutions?")
  
)

#endif