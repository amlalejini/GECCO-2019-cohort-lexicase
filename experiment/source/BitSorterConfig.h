#ifndef BIT_SORTER_CONFIG_H
#define BIT_SORTER_CONFIG_H

#include "config/config.h"

EMP_BUILD_CONFIG(BitSorterConfig, 
  GROUP(DEFAULT_GROUP, "General experiment settings"),
  VALUE(SEED, int, 2, "Random number seed (-1 for based on time)"),
  VALUE(GENERATIONS, size_t, 10000, "How many generations should we run the experiment for?"),
  VALUE(SORTER_POP_SIZE, size_t, 1024, "Population size for bit sorters"),
  VALUE(TEST_POP_SIZE, size_t, 1024, "Population size for test sorters"),
  VALUE(EVALUATION_MODE, size_t, 0, "How are programs evaluated? \n0: full (on all tests) \n1: cohorts "),
  VALUE(TEST_MODE, size_t, 0, "How do tests change over time? \n0: co-evolution \n1: static (unchanging) \n2: random"),

  GROUP(SELECTION_GROUP, "Settings specific to selection (both tests and sorters)"),
  VALUE(SORTER_SELECTION_MODE, size_t, 0, "How are sorters selected? \n0: LEXICASE \n1: COHORT_LEXICASE, \n2: TOURNAMENT, \n3: DRIFT"),
  VALUE(TEST_SELECTION_MODE, size_t, 0, "How are tests selected? \n0: LEXICASE \n1: COHORT_LEXICASE, \n2: TOURNAMENT, \n3: DRIFT"),
  VALUE(SORTER_COHORT_SIZE, size_t, 32, "."),
  VALUE(TEST_COHORT_SIZE, size_t, 32, "."),
  VALUE(TOURNAMENT_SIZE, size_t, 4, "."),
  VALUE(LEX_MAX_FUNS, size_t, 0, "Max number of lexicase functions to use (0 to use all)."),

  GROUP(SORTER_GROUP, "Settings specific to bit sorters."),
  VALUE(MAX_NETWORK_SIZE, size_t, 128, "Maximum size of a sorting network."),
  VALUE(MIN_NETWORK_SIZE, size_t, 1, "Minimum size of a sorting network."),

  GROUP(SORTER_MUTATION_GROUP, "Settings specific to mutating sorters"),
  VALUE(PER_SORTER_MUTATION, double, 1.0, "??"),
  VALUE(PER_INDEX_SUB, double, 0.001, "."),
  VALUE(PER_PAIR_DUP, double, 0.0005, "."),
  VALUE(PER_PAIR_INS, double, 0.0005, "."),
  VALUE(PER_PAIR_DEL, double, 0.001, "."),
  VALUE(PER_PAIR_SWAP, double, 0.001, "."),

  GROUP(TESTS_GROUP, "Settings specific to tests."),
  VALUE(TEST_SIZE, size_t, 16, "How many bits long are tests?"),

  GROUP(TESTS_MUTATION_GROUP, "Settings specific to tests."),
  VALUE(PER_BIT_FLIP, double, 0.001, "."),
  VALUE(PER_SEQ_INVERSION, double, 0.01, "."),
  VALUE(PER_SEQ_RANDOMIZE, double, 0.01, "."),

  GROUP(DATA_COLLECTION_GROUP, "Settings specific to data collection"),
  VALUE(SNAPSHOT_INTERVAL, size_t, 1000, "How often should we snapshot populations?")
  
)

#endif