#ifndef BIT_SORTER_EXP_H
#define BIT_SORTER_EXP_H

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
#include "control/Signal.h"
#include "Evolve/World.h"
#include "tools/Binomial.h"
#include "tools/Random.h"
#include "tools/random_utils.h"
#include "tools/math.h"
#include "tools/string_utils.h"
#include "tools/stats.h"

#include "BitSorterOrg.h"
#include "BitTestOrg.h"
#include "Selection.h"
#include "Mutators.h"
#include "BitSorterMutators.h"

#include "BitSorterConfig.h"

// TODOs
// - [ ] Setup crossover
// - [ ] Setup data collection

class BitSorterExperiment {
public:

  // Some useful aliases.
  using sorter_org_t = BitSorterOrg;
  using sorter_genome_t = BitSorterOrg::genome_t;

  using test_org_t = BitTestOrg;
  using test_genome_t = BitTestOrg::genome_t;

  using sorter_world_t = emp::World<sorter_org_t>;
  using test_world_t = emp::World<test_org_t>;

  // Experiment toggles
  enum SELECTION_METHODS { LEXICASE=0, COHORT_LEXICASE=1, TOURNAMENT=2, DRIFT=3 };
  enum TEST_MODES { COEVOLVE=0, STATIC=1, RANDOM=2 };
  enum EVALUATION_MODES { FULL=0, COHORT=1 };
  enum SORTER_CROSSOVER_MODES { NONE=0, SINGLE_PT=1, TWO_PT=2 };

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
  
  //////////////////////////////////////////////////////////////////////////////

  // Localized parameters
  int SEED;
  size_t GENERATIONS;
  size_t SORTER_POP_SIZE;
  size_t TEST_POP_SIZE;
  size_t EVALUATION_MODE;
  size_t TEST_MODE;

  size_t SORTER_SELECTION_MODE;
  size_t TEST_SELECTION_MODE;
  size_t SORTER_COHORT_SIZE;
  size_t TEST_COHORT_SIZE;
  size_t TOURNAMENT_SIZE;
  size_t LEX_MAX_FUNS;

  size_t MAX_NETWORK_SIZE;
  size_t MIN_NETWORK_SIZE;

  double PER_SORTER_MUTATION;
  double PER_INDEX_SUB;
  double PER_PAIR_DUP;
  double PER_PAIR_INS;
  double PER_PAIR_DEL;
  double PER_PAIR_SWAP;

  size_t TEST_SIZE;

  double PER_BIT_FLIP;
  double PER_SEQ_INVERSION;
  double PER_SEQ_RANDOMIZE;

  size_t SNAPSHOT_INTERVAL;

  //////////////////////////////////////////////////////////////////////////////

  // Local experiment variables
  bool setup;                      ///< Has setup been run?
  size_t update;
  
  size_t dominant_sorter_id;
  
  size_t smallest_known_sol_size;
  bool solution_found;
  
  size_t MAX_SORTER_PASSES;        ///< Maximum number of test passes for a sorter to achieve
  size_t SORTER_EVAL_TESTCASE_CNT;
  size_t NUM_COHORTS;

  emp::Ptr<emp::Random> random;

  emp::Ptr<sorter_world_t> sorter_world;
  emp::Ptr<test_world_t> test_world;

  emp::Ptr<emp::DataFile> solution_file;

  Cohorts sorter_cohorts;
  Cohorts test_cohorts;

  emp::vector<std::function<double(sorter_org_t &)>> lexicase_sorter_fit_set;
  emp::vector<std::function<double(test_org_t &)>> lexicase_test_fit_set;

  // Mutators
  BitSorterMutator sorter_mutator;
  BitTestMutator test_mutator;

  struct StatsUtil {
    size_t sorterID;
    size_t testID;
  } stats_util;

  // Experiment signals
  emp::Signal<void(void)> do_evaluation_sig;  ///< Trigger network/test evaluations.
  emp::Signal<void(void)> do_selection_sig;   ///< Trigger selection
  emp::Signal<void(void)> do_update_sig;      ///< Trigger world updates

  emp::Signal<void(void)> do_snapshot_sig;
  emp::Signal<void(void)> do_sorter_snapshot_sig; ///< Trigger population snapshot
  emp::Signal<void(void)> do_test_snapshot_sig; ///< Trigger population snapshot
  emp::Signal<void(void)> end_setup_sig;       ///< Triggered at beginning of a run.

  void InitConfigs(const BitSorterConfig & config);  ///< Initialize (localize) configuration settings.

  void InitSorterPop_Random();  ///< Randomly initialize the sorter population.
  void InitTestPop_Random();    ///< Randomly initialize the test population.

  void SetupEvaluation();
  void SetupSelection();
  void SetupMutation();
  void SetupFitFuns();

  void SetupSorterSelection();
  void SetupSorterMutation();
  void SetupSorterFitFun();
  
  void SetupTestSelection();
  void SetupTestMutation();
  void SetupTestFitFun();

  void SetupDataCollection();

  void SnapshotSorters();
  void SnapshotTests();

public:

  BitSorterExperiment()
    : setup(false), update(0),
      sorter_cohorts(), test_cohorts(),
      sorter_mutator(), test_mutator()
  { ; }

  ~BitSorterExperiment() {
    if (setup) {
      #ifndef EMSCRIPTEN
      solution_file.Delete();
      #endif
      sorter_world.Delete();
      test_world.Delete();
      random.Delete();
    }
  }

  /// Configure the experiment using the given config object.
  void Setup(const BitSorterConfig & config);

  /// Run the experiment start->finish
  /// (1) Initialize population(s)
  /// (2) For each generation, RunStep()
  void Run();

  /// Progress experiment by a single time step (generation):
  /// - (1) Evaluate populations
  /// - (2) Select who gets to reproduce
  /// - (3) Update the worlds
  void RunStep();

};


void BitSorterExperiment::Setup(const BitSorterConfig & config) {
  std::cout << "Setting up a BitSorterExperiment." << std::endl;

  // Initialize configs.
  InitConfigs(config);

  if (!setup) { // First call to Setup.
    // Create a new random number generator.
    random = emp::NewPtr<emp::Random>(SEED);
    sorter_world = emp::NewPtr<sorter_world_t>(*random, "Sorter World");
    test_world = emp::NewPtr<test_world_t>(*random, "Sorting test world");
  } else {
    std::cout << "Not allowed to call setup more than once!" << std::endl;
    exit(-1);
  }

  // Configure world pop structures.
  sorter_world->SetPopStruct_Mixed(true);
  test_world->SetPopStruct_Mixed(true);

  smallest_known_sol_size = MAX_NETWORK_SIZE + 1;
  solution_found = false;

  // How does the sorter world update?
  do_update_sig.AddAction([this]() {
    std::cout << "Update: " << update << ", ";
    std::cout << "best sorter (size=" << sorter_world->GetOrg(dominant_sorter_id).GetSize() << "): " << sorter_world->CalcFitnessID(dominant_sorter_id) << ", ";
    std::cout << "solution? " << (size_t)solution_found << ", smallest solution size: " << smallest_known_sol_size << std::endl;

    if (update % SNAPSHOT_INTERVAL) do_snapshot_sig.Trigger();

    sorter_world->Update();
    sorter_world->ClearCache();
  });

  // Setup test world update
  std::cout << "==== EXPERIMENT SETUP => test world update ====" << std::endl;
  switch (TEST_MODE) {
    case (size_t)TEST_MODES::COEVOLVE: {
      std::cout << "TEST MODE = COEVOLVE" << std::endl;
      do_update_sig.AddAction([this]() {
        test_world->Update();
        test_world->ClearCache();
      });
      break;
    }
    case (size_t)TEST_MODES::STATIC: {
      std::cout << "TEST MODE = STATIC" << std::endl;
      do_update_sig.AddAction([this]() {
        test_world->ClearCache();
        // Reset test phenotypes
        for (size_t i = 0; i < test_world->GetSize(); ++i) {
          emp_assert(test_world->IsOccupied(i));
          test_org_t & test_org = test_world->GetOrg(i);
          test_org.GetPhenotype().Reset(test_org.GetPhenotype().test_passes.size());
        }
      });
      break;
    }
    case (size_t)TEST_MODES::RANDOM: {
      std::cout << "TEST MODE = RANDOM" << std::endl;
      do_update_sig.AddAction([this]() {
        test_world->ClearCache();
        // Randomize population
        test_world->DoMutations();
        // Reset test phenotypes
        for (size_t i = 0; i < test_world->GetSize(); ++i) {
          emp_assert(test_world->IsOccupied(i));
          test_org_t & test_org = test_world->GetOrg(i);
          test_org.GetPhenotype().Reset(test_org.GetPhenotype().test_passes.size());
        }
      });
      break;
    }
    default: {
      std::cout << "Unknown TEST_MODE (" << TEST_MODE << "). Exiting." << std::endl;
      exit(-1);
    }
  }

  end_setup_sig.AddAction([this]() {
    // Initialize populations randomly
    InitSorterPop_Random();
    InitTestPop_Random();
  });

  std::cout << "==== EXPERIMENT SETUP => evaluation ====" << std::endl;
  SetupEvaluation();

  // Setup sorter (& test) selection.
  std::cout << "==== EXPERIMENT SETUP => selection ====" << std::endl;
  SetupSelection();

  std::cout << "==== EXPERIMENT SETUP => mutation ====" << std::endl;
  SetupMutation();

  // Setup sorter (& test) fitness calculations.
  std::cout << "==== EXPERIMENT SETUP => fitness functions ====" << std::endl;
  SetupFitFuns();

  // #ifndef EMSCRIPTEN
  // SetupDataCollection();
  // #endif

  // todo - setup sorter snapshots
  // todo - setup test snapshots

  end_setup_sig.Trigger();
  setup = true;
  std::cout << "Finished experiment setup." << std::endl;
}

/// Run the experiment start->finish
/// (1) Initialize population(s)
/// (2) For each generation, RunStep()
void BitSorterExperiment::Run() {
  for (update = 0; update <= GENERATIONS; ++update) {
    RunStep();
  }
}

/// Progress experiment by a single time step (generation):
/// - (1) Evaluate populations
/// - (2) Select who gets to reproduce
/// - (3) Update the worlds
void BitSorterExperiment::RunStep() {
  // std::cout << "-- Doing Evaluation --" << std::endl;
  do_evaluation_sig.Trigger();
  // std::cout << "-- Doing Selection --" << std::endl;
  do_selection_sig.Trigger();
  // std::cout << "-- Doing Update --" << std::endl;
  do_update_sig.Trigger();
}

void BitSorterExperiment::InitConfigs(const BitSorterConfig & config) {
  std::cout << "Initializing configs." << std::endl;
  
  SEED = config.SEED();
  GENERATIONS = config.GENERATIONS();
  SORTER_POP_SIZE = config.SORTER_POP_SIZE();
  TEST_POP_SIZE = config.TEST_POP_SIZE();
  EVALUATION_MODE = config.EVALUATION_MODE();
  TEST_MODE = config.TEST_MODE();

  SORTER_SELECTION_MODE = config.SORTER_SELECTION_MODE();
  TEST_SELECTION_MODE = config.TEST_SELECTION_MODE();
  SORTER_COHORT_SIZE = config.SORTER_COHORT_SIZE();
  TEST_COHORT_SIZE = config.TEST_COHORT_SIZE();
  TOURNAMENT_SIZE = config.TOURNAMENT_SIZE();
  LEX_MAX_FUNS = config.LEX_MAX_FUNS();

  MAX_NETWORK_SIZE = config.MAX_NETWORK_SIZE();
  MIN_NETWORK_SIZE = config.MIN_NETWORK_SIZE();

  PER_SORTER_MUTATION = config.PER_SORTER_MUTATION();
  PER_INDEX_SUB = config.PER_INDEX_SUB();
  PER_PAIR_DUP = config.PER_PAIR_DUP();
  PER_PAIR_INS = config.PER_PAIR_INS();
  PER_PAIR_DEL = config.PER_PAIR_DEL();
  PER_PAIR_SWAP = config.PER_PAIR_SWAP();

  TEST_SIZE = config.TEST_SIZE();

  PER_BIT_FLIP = config.PER_BIT_FLIP();
  PER_SEQ_INVERSION = config.PER_SEQ_INVERSION();
  PER_SEQ_RANDOMIZE = config.PER_SEQ_RANDOMIZE();

  SNAPSHOT_INTERVAL = config.SNAPSHOT_INTERVAL();

  emp_assert(TEST_SIZE <= 32);
}

/// 
void BitSorterExperiment::SetupEvaluation() {
  std::cout << "How do we evaluate sorting networks against tests?" << std::endl;
  switch (EVALUATION_MODE) {
    case (size_t)EVALUATION_MODES::COHORT: {
      std::cout << "  => Evaluate in cohorts." << std::endl;
      // Initialize number of cohorts.
      sorter_cohorts.Setup(SORTER_POP_SIZE, SORTER_COHORT_SIZE);
      std::cout << "Number of sorter cohorts = " << sorter_cohorts.GetCohortCnt() << "; Cohort size = " << sorter_cohorts.GetCohortSize() << std::endl;
      test_cohorts.Setup(TEST_POP_SIZE, TEST_COHORT_SIZE);
      std::cout << "Number of test cohorts = " << test_cohorts.GetCohortCnt() << "; Cohort size = " << test_cohorts.GetCohortSize() << std::endl;
      // Double check that num cohorts will be equal across tests and networks.
      if (test_cohorts.GetCohortCnt() != sorter_cohorts.GetCohortCnt()) {
        // Oh no!
        std::cout << "The number of sorter and test cohorts must be equal. Exiting." << std::endl;
        exit(-1);
      }
      // What's the max number of passes a sorting network can achieve?
      MAX_SORTER_PASSES = TEST_COHORT_SIZE;
      std::cout << "Max number of test passes achievable by a sorting network = " << MAX_SORTER_PASSES << std::endl;
      // Setup world to reset phenotypes on organism placement.
      sorter_world->OnPlacement([this](size_t pos) { sorter_world->GetOrg(pos).GetPhenotype().Reset(TEST_COHORT_SIZE); });
      test_world->OnPlacement([this](size_t pos) { test_world->GetOrg(pos).GetPhenotype().Reset(SORTER_COHORT_SIZE); });
      // What should happen on evaluation?
      do_evaluation_sig.AddAction([this]() {
        // Randomize cohorts.
        sorter_cohorts.Randomize(*random);
        test_cohorts.Randomize(*random);
        // For each cohort, evaluate all sorters against all tests in cohort.
        for (size_t cohortID = 0; cohortID < sorter_cohorts.GetCohortCnt(); ++cohortID) {
          for (size_t sorterID = 0; sorterID < SORTER_COHORT_SIZE; ++sorterID) {
            // Test this sorter against all tests in associated cohort.
            sorter_org_t & sorter_org = sorter_world->GetOrg(sorter_cohorts.GetWorldID(cohortID, sorterID)); 
            for (size_t testID = 0; testID < TEST_COHORT_SIZE; ++testID) {
              // Evaluate this test against current sorter org.
              test_org_t & test_org = test_world->GetOrg(test_cohorts.GetWorldID(cohortID, testID));
              // Evaluate!
              bool can_sort = sorter_org.GetGenome().TestSortable(test_org.GetGenome());
              // Update sorter phenotype
              sorter_org.GetPhenotype().RecordPassFail(testID, can_sort);
              sorter_org.GetPhenotype().RecordScore(testID, (double)can_sort);
              // Update test phenotype
              test_org.GetPhenotype().RecordPassFail(sorterID, can_sort);
              test_org.GetPhenotype().RecordScore(sorterID, !can_sort); // Score here is 1 if not sorted, 0 if sorted.
            }
          }
        } 
      });
      break;
    }
    case (size_t)EVALUATION_MODES::FULL: {
      std::cout << "  => Evaluate against full test population." << std::endl;
      // Setup world to reset phenotypes on organism placement
      sorter_world->OnPlacement([this](size_t pos) { sorter_world->GetOrg(pos).GetPhenotype().Reset(TEST_POP_SIZE); });
      test_world->OnPlacement([this](size_t pos) { test_world->GetOrg(pos).GetPhenotype().Reset(SORTER_POP_SIZE); });
      MAX_SORTER_PASSES = TEST_POP_SIZE;
      // What should happen on evaluation?
      do_evaluation_sig.AddAction([this]() {
        for (size_t sorterID = 0; sorterID < sorter_world->GetSize(); ++sorterID) {
          sorter_org_t & sorter_org = sorter_world->GetOrg(sorterID);
          for (size_t testID = 0; testID < test_world->GetSize(); ++testID) {
            test_org_t & test_org = test_world->GetOrg(testID);
            // Evaluate sorter and test.
            bool can_sort = sorter_org.GetGenome().TestSortable(test_org.GetGenome());
            // Update sorter phenotype
            sorter_org.GetPhenotype().RecordPassFail(testID, can_sort);
            sorter_org.GetPhenotype().RecordScore(testID, (double)can_sort);
            // Update test phenotype
            test_org.GetPhenotype().RecordPassFail(sorterID, can_sort);
            test_org.GetPhenotype().RecordScore(sorterID, !can_sort); // Score here is 1 if not sorted, 0 if sorted.
          }
        }
      });
      break;
    }
    default: {
      std::cout << "Unknown EVALUATION_MODE (" << EVALUATION_MODE << "). Exiting." << std::endl;
    }
  }

  // Add evaluation action to screen for solutions
  do_evaluation_sig.AddAction([this]() {
    // Find a 'dominant' sorting network.
    double cur_best_score = 0;
    for (size_t sID = 0; sID < sorter_world->GetSize(); ++sID) {
      emp_assert(sorter_world->IsOccupied(sID));
      sorter_org_t & sorter_org = sorter_world->GetOrg(sID);
      const size_t pass_total = sorter_org.GetPhenotype().num_passes;
      const double score_total = sorter_org.GetPhenotype().total_score;

      // Is this the highest score for this generation?
      if (score_total > cur_best_score || sID == 0) {
        dominant_sorter_id = sID;
        cur_best_score = score_total;
      }
      // At this point, sorting networks have been evaluated against all tests.
      if (pass_total == MAX_SORTER_PASSES && sorter_org.GetGenome().GetSize() < smallest_known_sol_size) {
        stats_util.sorterID = sID;
        if (sorter_org.GetGenome().IsCorrect(TEST_SIZE)) {
          solution_found = true;
          smallest_known_sol_size = sorter_org.GetGenome().GetSize();
          // Add to solutions file.
          // solution_file->Update(); // TODO - setup solution file
        }
      }
    }
  });
}

void BitSorterExperiment::SetupSelection() {
  SetupSorterSelection();
  if (TEST_MODE == TEST_MODES::COEVOLVE) { SetupTestSelection(); }
}

void BitSorterExperiment::SetupSorterSelection() {
  std::cout << "Setting up sorter selection." << std::endl;
  switch (SORTER_SELECTION_MODE) {
    case (size_t)SELECTION_METHODS::LEXICASE: {
      std::cout << "  Sorter selection method: LEXICASE" << std::endl;
      emp_assert(EVALUATION_MODE == EVALUATION_MODES::FULL);
      // For lexicase selection, one function for every test.
      for (size_t i = 0; i < TEST_POP_SIZE; ++i) {
        lexicase_sorter_fit_set.push_back([this, i](sorter_org_t & sorter) {
          return sorter.GetPhenotype().test_scores[i];
        });
      }
      // Add one function for size.
      lexicase_sorter_fit_set.push_back([this](sorter_org_t & sorter) {
        if (sorter.GetPhenotype().num_passes == MAX_SORTER_PASSES) return (double)MAX_NETWORK_SIZE - (double)sorter.GetGenome().GetSize();
        else return 0.0;
      });
      // Setup selection signal action.
      do_selection_sig.AddAction([this]() {
        std::cout << "=>Enter: Sorter Lexicase selection" << std::endl;
        emp::LexicaseSelect_NAIVE(*sorter_world, 
                                  lexicase_sorter_fit_set, 
                                  SORTER_POP_SIZE,
                                  LEX_MAX_FUNS);
        std::cout << "=>Exit: Sorter Lexicase selection" << std::endl;
      });
      break;
    }
    case (size_t)SELECTION_METHODS::COHORT_LEXICASE: {
      std::cout << "  Sorter selection method: COHORT_LEXICASE" << std::endl;
      emp_assert(EVALUATION_MODE == EVALUATION_MODES::COHORT);
      // For lexicase selection, one function for every test.
      for (size_t i = 0; i < TEST_COHORT_SIZE; ++i) {
        lexicase_sorter_fit_set.push_back([this, i](sorter_org_t & sorter) {
          return sorter.GetPhenotype().test_scores[i];
        });
      }
      // Add one function for size.
      lexicase_sorter_fit_set.push_back([this](sorter_org_t & sorter) {
        if (sorter.GetPhenotype().num_passes == MAX_SORTER_PASSES) return (double)MAX_NETWORK_SIZE - (double)sorter.GetGenome().GetSize();
        else return 0.0;
      });
      // Setup selection signal action.
      do_selection_sig.AddAction([this]() {
        // For each cohort, run selection
        for (size_t cID = 0; cID < sorter_cohorts.GetCohortCnt(); ++cID) {
          emp::CohortLexicaseSelect_NAIVE(*sorter_world, 
                                          lexicase_sorter_fit_set,
                                          sorter_cohorts.GetCohort(cID),
                                          SORTER_COHORT_SIZE,
                                          LEX_MAX_FUNS);
        }
      });
      break;
    }
    case (size_t)SELECTION_METHODS::TOURNAMENT: {
      std::cout << "  Sorter selection method: TOURNAMENT" << std::endl;
      do_selection_sig.AddAction([this]() {
        emp::TournamentSelect(*sorter_world, TOURNAMENT_SIZE, SORTER_POP_SIZE);
      });
      break;
    }
    case (size_t)SELECTION_METHODS::DRIFT: {
      std::cout << "  Sorter selection method: DRIFT" << std::endl;
      do_selection_sig.AddAction([this]() {
        emp::RandomSelect(*sorter_world, SORTER_POP_SIZE);
      });
      break;
    }
    default: {
      std::cout << "Unknown SORTER_SELECTION_MODE (" << SORTER_SELECTION_MODE << "). Exiting." << std::endl;
      exit(-1);
    }
  }
}

void BitSorterExperiment::SetupTestSelection() {
  std::cout << "Setting up test selection." << std::endl;
  switch (TEST_SELECTION_MODE) {
    case (size_t)SELECTION_METHODS::LEXICASE: {
      std::cout << "  Test selection method: LEXICASE" <<  std::endl;
      emp_assert(EVALUATION_MODE == EVALUATION_MODES::FULL);
      for (size_t i = 0; i < SORTER_POP_SIZE; ++i) {
        lexicase_test_fit_set.push_back([i](test_org_t & test) { 
          return 1 - test.GetPhenotype().test_passes[i]; // Max if test_pass[i] = 0
        });
      }
      do_selection_sig.AddAction([this]() {
        std::cout << "=>Enter: Test Lexicase selection" << std::endl;
        emp::LexicaseSelect_NAIVE(*test_world, lexicase_test_fit_set, TEST_POP_SIZE, LEX_MAX_FUNS);
        std::cout << "=>Exit: Test Lexicase selection" << std::endl;
      });
      break;
    }
    case (size_t)SELECTION_METHODS::COHORT_LEXICASE: {
      std::cout << "  Test selection method: COHORT_LEXICASE" <<  std::endl;
      emp_assert(EVALUATION_MODE == EVALUATION_MODES::COHORT);
      for (size_t i = 0; i < SORTER_COHORT_SIZE; ++i) {
        lexicase_test_fit_set.push_back([i](test_org_t & test) { 
          return 1 - test.GetPhenotype().test_passes[i]; // Max if test_pass[i] = 0
        });
      }
      do_selection_sig.AddAction([this]() {
        for (size_t cID = 0; cID < test_cohorts.GetCohortCnt(); ++cID) {
          emp::CohortLexicaseSelect_NAIVE(*test_world,
                                          lexicase_test_fit_set,
                                          test_cohorts.GetCohort(cID),
                                          TEST_COHORT_SIZE,
                                          LEX_MAX_FUNS);
        }
      });
      break;
    }
    case (size_t)SELECTION_METHODS::TOURNAMENT: {
      std::cout << "  Test selection method: TOURNAMENT" <<  std::endl;
      do_selection_sig.AddAction([this]() {
        emp::TournamentSelect(*test_world, TOURNAMENT_SIZE, TEST_POP_SIZE);
      });
      break;
    }
    case (size_t)SELECTION_METHODS::DRIFT: {
      std::cout << "  Test selection method: DRIFT" <<  std::endl;
      do_selection_sig.AddAction([this]() {
        emp::RandomSelect(*test_world, TEST_POP_SIZE);
      });
      break;
    }
    default: {
      std::cout << "Unknown TEST_SELECTION_MODE (" << TEST_SELECTION_MODE << "). Exiting." << std::endl;
      exit(-1);
    }
  }
}

void BitSorterExperiment::SetupMutation() {
  // (1) Setup sorter mutations
  SetupSorterMutation();
  // (2) Setup test mutations
  SetupTestMutation();
  // (3) Setup world(s) to auto mutate
  end_setup_sig.AddAction([this]() {
    sorter_world->SetAutoMutate();
    if (TEST_MODE == (size_t)TEST_MODES::COEVOLVE) { test_world->SetAutoMutate(); }
  });
}

void BitSorterExperiment::SetupSorterMutation() {
  std::cout << "Setting up sorter mutation." << std::endl;
  // Configure sorter mutator.
  sorter_mutator.MAX_NETWORK_SIZE = MAX_NETWORK_SIZE;
  sorter_mutator.MIN_NETWORK_SIZE = MIN_NETWORK_SIZE;
  sorter_mutator.SORT_SEQ_SIZE = TEST_SIZE;
  sorter_mutator.PER_INDEX_SUB = PER_INDEX_SUB;
  sorter_mutator.PER_PAIR_DUP = PER_PAIR_DUP;
  sorter_mutator.PER_PAIR_INS = PER_PAIR_INS;
  sorter_mutator.PER_PAIR_DEL = PER_PAIR_DEL;
  sorter_mutator.PER_PAIR_SWAP = PER_PAIR_SWAP;

  if (PER_SORTER_MUTATION) {
    sorter_world->SetMutFun([this](sorter_org_t & sorter_org, emp::Random & rnd) {
      // return 1.0;
      return sorter_mutator.Mutate(rnd, sorter_org.GetGenome());
    });
  } else {
    sorter_world->SetMutFun([this](sorter_org_t & sorter_org, emp::Random & rnd) {
      if (rnd.P(PER_SORTER_MUTATION)) {
        return (double)sorter_mutator.Mutate(rnd, sorter_org.GetGenome());
      } else {
        return 0.0;
      }
    });
  }

  // TODO - configure crossover (with do_crossover_sig(?))
}

void BitSorterExperiment::SetupTestMutation() { // todo - if test mode == random, ...
  std::cout << "Setting up test mutation." << std::endl;
  test_mutator.NUM_BITS = TEST_SIZE;
  test_mutator.PER_BIT_FLIP = PER_BIT_FLIP;
  test_mutator.PER_SEQ_INVERSION = PER_SEQ_INVERSION;
  test_mutator.PER_SEQ_RANDOMIZE = PER_SEQ_RANDOMIZE;
  if (TEST_MODE == TEST_MODES::RANDOM) {
    test_world->SetMutFun([this](test_org_t & test, emp::Random & rnd) {
      test.SetGenome(rnd.GetUInt(0, 1 << TEST_SIZE));  // Randomize genome
      return 1;
    });
  } else {
    test_world->SetMutFun([this](test_org_t & test, emp::Random & rnd) {
      return test_mutator.Mutate(rnd, test.GetGenome());
    });
  }

}

void BitSorterExperiment::SetupFitFuns() {
  // (1) Setup sorter fitness function
  SetupSorterFitFun();
  // (2) Setup test fitness function
  SetupTestFitFun();
}

void BitSorterExperiment::SetupSorterFitFun() {
  sorter_world->SetFitFun([this](sorter_org_t & sorter_org) {
    if (sorter_org.GetPhenotype().num_passes == MAX_SORTER_PASSES) {
      const double test_score = (double)sorter_org.GetPhenotype().total_score;
      const double size_bonus = ((double)MAX_NETWORK_SIZE - (double)sorter_org.GetGenome().GetSize())/(double)MAX_NETWORK_SIZE;
      return test_score + size_bonus;
    } else {
      return (double)sorter_org.GetPhenotype().total_score;
    }
  });
}

void BitSorterExperiment::SetupTestFitFun() {
  test_world->SetFitFun([this](test_org_t & test_org) {
    return (double)test_org.GetPhenotype().num_fails;
  });
}

/// Initialize bit sorter population randomly
void BitSorterExperiment::InitSorterPop_Random() {
  std::cout << "SORTER POP INIT - randomly initializing sorter population." << std::endl;
  // Inject random networks into the sorter world up to population size.
  for (size_t i = 0; i < SORTER_POP_SIZE; ++i) {
    sorter_world->Inject(sorter_mutator.GenRandomBitSorter(*random));
  }
  std::cout << "Done with sorter pop init (" << sorter_world->GetSize() << ")." << std::endl;
}

/// Initialize test population randomly. If possible, guarantee test uniqueness.
void BitSorterExperiment::InitTestPop_Random() {
  std::cout << "TEST POP INIT - randomly initializing test population." << std::endl;
  // Generate all possible tests.
  size_t total_tests = emp::Pow2(TEST_SIZE);
  emp::vector<uint32_t> all_tests(total_tests, 0);
  for (uint32_t val = 0; val < all_tests.size(); ++val) {
    all_tests[val] = val;
  }
  std::cout << "  Total possible tests = " << all_tests.size() << "([" << all_tests.front() << "," << all_tests.back() << "])" << std::endl;
  // Inject tests into the population.
  while (test_world->GetSize() < TEST_POP_SIZE) {
    emp::Shuffle(*random, all_tests);
    for (size_t i = 0; i < all_tests.size() && test_world->GetSize() < TEST_POP_SIZE; ++i) {
      test_world->Inject(all_tests[i]);
    }
  }
  std::cout << "Done with test population init (" << test_world->GetSize() << "). " << std::endl;
}



#endif