#ifndef SORTING_NETWORK_EXP_H
#define SORTING_NETWORK_EXP_H

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

#include "SortingNetworkConfig.h"
#include "SortingNetworkOrg.h"
#include "SortingTestOrg.h"
#include "Selection.h"
#include "Mutators.h"

/*

SortingNetworkOrgs - Phenotype:

- [test org 0 perf, test org 1 perf, ...]
  - where each performance score is a pass count
- Num passes
- Network size

SortingTestOrgs - Phenotype

- [network org 0 perf, ]

Major components:

- SortingNetworkWorld
- SortingTestWorld

Experiment configuration/setup:

- InitConfigs

...
- SetupMut, SetupFitness

Experiment execution flow:

- (1) Evaluation
- (2) Selection
- (3) World update

*/

class SortingNetworkExperiment {
public:

  // Setup some convenient aliases.
  using network_org_t = SortingNetworkOrg;
  using network_genome_t = SortingNetworkOrg::genome_t;
  
  using test_org_t = SortingTestOrg;
  using test_genome_t = SortingTestOrg::genome_t;
  
  using network_world_t = emp::World<network_org_t>;
  using test_world_t = emp::World<test_org_t>;

  // Experiment toggles
  enum SELECTION_METHODS { LEXICASE=0, COHORT_LEXICASE=1, TOURNAMENT=2 };
  enum TEST_MODES { COEVOLVE=0, STATIC=1, RANDOM=2, DRIFT=3 };
  enum NETWORK_CROSSOVER_MODES { NONE=0, SINGLE_PT=1, TWO_PT=2 };
  
protected:

  // Localized experiment parameters
  int SEED;
  size_t GENERATIONS;
  size_t NETWORK_POP_SIZE;
  size_t TEST_POP_SIZE;
  size_t TEST_MODE;
  size_t COHORT_SIZE;

  size_t SELECTION_MODE;
  size_t LEX_MAX_FUNS;
  size_t COHORTLEX_MAX_FUNS;
  size_t TOURNAMENT_SIZE;
  bool DISCRIMINATORY_LEXICASE_TESTS;

  size_t MAX_NETWORK_SIZE;
  size_t MIN_NETWORK_SIZE;

  double PER_INDEX_SUB;
  double PER_PAIR_DUP;
  double PER_PAIR_INS;
  double PER_PAIR_DEL;
  double PER_PAIR_SWAP;
  size_t NETWORK_CROSSOVER_MODE;
  double PER_ORG_CROSSOVER;
  double PER_ORG_MUTATION;

  size_t SORT_SIZE;
  size_t SORTS_PER_TEST;

  double PER_SITE_SUB;
  double PER_SEQ_INVERSION;
  double PER_SEQ_RANDOMIZE;

  std::string DATA_DIRECTORY;
  size_t SNAPSHOT_INTERVAL;
  size_t DOMINANT_STATS_INTERVAL;
  size_t AGGREGATE_STATS_INTERVAL;
  size_t CORRECTNESS_SAMPLE_SIZE;
  size_t SOLUTION_SCREEN_INTERVAL;
  bool COLLECT_TEST_PHYLOGENIES;

  // Experiment variables
  bool setup;                 ///< Has setup been run?
  size_t update;
  size_t dominant_network_id;
  size_t dominant_test_id;
  size_t smallest_known_sol_size;
  size_t MAX_PASSES;

  emp::Ptr<emp::Random> random;

  emp::Ptr<network_world_t> network_world;
  emp::Ptr<test_world_t> test_world;

  // emp::Ptr<emp::Systematics<network_org_t, network_genome_t>> network_genotypic_systematics;
  emp::Ptr<emp::Systematics<test_org_t, test_genome_t>> test_genotypic_systematics;


  struct Cohorts {
    emp::vector<size_t> population_ids;
    emp::vector<emp::vector<size_t>> cohorts;

    size_t cohort_size;

    /// Used internally to access population_ids
    size_t GetOrgID(size_t cohortID, size_t memberID) const {
      return (cohortID * cohort_size) + memberID;
    }

    /// Used to get world ID of organism #memberID belonging to given cohortID
    size_t GetWorldID(size_t cohortID, size_t memberID) const {
      emp_assert(cohortID < cohorts.size());
      emp_assert(memberID < cohorts[cohortID].size());
      return cohorts[cohortID][memberID];
    }

    size_t GetNumCohorts() const { return cohorts.size(); }
    size_t GetCohortSize() const { return cohort_size; }
    const emp::vector<size_t> & GetCohort(size_t id) const { 
      emp_assert(id < cohorts.size());
      return cohorts[id]; 
    }

    void Init(size_t _pop_size, size_t _cohort_size) {
      emp_assert(_cohort_size > 0);
      population_ids.clear();
      cohorts.clear();
      cohort_size = _cohort_size;
      const size_t num_cohorts = _pop_size / _cohort_size;
      // Initialize population ids
      for (size_t i = 0; i < _pop_size; ++i) population_ids.emplace_back(i);
      // Initialize cohort vectors
      for (size_t cID = 0; cID < num_cohorts; ++cID) {
        cohorts.emplace_back(cohort_size);
        for (size_t i = 0; i < cohort_size; ++i) cohorts[cID][i] = population_ids[GetOrgID(cID, i)];
      }
    }

    void Randomize(emp::Random & rnd) {
      // Shuffle population ids
      emp::Shuffle(rnd, population_ids);
      // Reassign cohorts
      for (size_t cID = 0; cID < cohorts.size(); ++cID) {
        for (size_t i = 0; i < cohorts[cID].size(); ++i) {
          cohorts[cID][i] = population_ids[GetOrgID(cID, i)];
        }
      }
    }

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

  Cohorts network_cohorts;
  Cohorts test_cohorts;

  emp::vector<std::function<double(network_org_t &)>> lexicase_network_fit_set;
  emp::vector<std::function<double(test_org_t &)>> lexicase_test_fit_set;

  // Mutators
  SortingNetworkMutator network_mutator;
  SortingTestMutator test_mutator;

  struct StatID {
    size_t networkID;
    size_t testID;
    StatID(size_t nID=0, size_t tID=0) : networkID(nID), testID(tID) { ; }
  } curIDs;

  struct PopIDs {
    emp::vector<size_t> popIDs;
    void Init(size_t pop_size) {
      popIDs.clear();
      popIDs.resize(pop_size);
      for (size_t i = 0; i < pop_size; ++i) popIDs[i] = i;
    }
    void Shuffle(emp::Random & rnd) {
      emp::Shuffle(rnd, popIDs);
    }
  } network_pop_ids;

  emp::Ptr<emp::Binomial> network_cross_binomial;

  struct CompleteTestSet {
    emp::vector<size_t> testIDs;
    emp::vector<SortingTest> tests;
    size_t swap_id;

    size_t GetSize() const { return tests.size(); }

    void Generate(size_t test_size) {
      testIDs.clear();
      tests.clear();
      size_t total_tests = emp::Pow2(test_size);
      tests.resize(total_tests, test_size);
      size_t interval = 1;
      bool val = false;
      for (size_t pos = 0; pos < test_size; ++pos) {
        for (size_t tID = 0; tID < tests.size(); ++tID) {
          // Set test at this position to val.
          tests[tID][pos] = (int)val;
          // Change value on interval
          if ((tID % interval) == 0) {
            val = (val) ? false : true;
          }
        }
        interval *= 2;
      }
      // Print!
      // std::unordered_set<std::string> set;
      // for (size_t tID = 0; tID < tests.size(); ++tID) {
      //   tests[tID].Print(); std::cout << std::endl;
      //   set.insert(emp::to_string(tests[tID].GetTest()));
      // }
      for (size_t i = 0; i < tests.size(); ++i) testIDs.emplace_back(i);
    }

    void SuffleTestIDs(emp::Random & rnd) {
      emp::Shuffle(rnd, testIDs);
    }

    size_t EvaluateAll(const SortingNetwork & network) {
      size_t sort_cnt = 0;
      for (size_t i = 0; i < tests.size(); ++i) {
        sort_cnt += (size_t)tests[i].Evaluate(network);
      }
      return sort_cnt;
    }

    size_t EvaluateN(const SortingNetwork & network, size_t N) {
      size_t sort_cnt = 0;
      for (size_t i = 0; i < N; ++i) {
        sort_cnt += (size_t)tests[testIDs[i]].Evaluate(network);
      }
      return sort_cnt;
    }

    // TODO - push up tests likely to cause failure
    bool Correct(const SortingNetwork & network) {
      for (size_t i = 0; i < tests.size(); ++i) {
        if (!tests[i].Evaluate(network)) {
          std::swap(tests[i], tests[swap_id]);
          ++swap_id; 
          if (swap_id >= tests.size()) swap_id = 0;
          return false;
        }
      }
      return true;
    }

  } complete_test_set;

  // Network stats
  std::function<size_t(void)> get_networkID;
  std::function<double(void)> get_network_fitness;
  std::function<double(void)> get_network_true_correct;
  std::function<double(void)> get_network_possible_correct;
  std::function<double(void)> get_network_sample_correct;
  std::function<double(void)> get_network_possible_sample_correct;
  std::function<size_t(void)> get_network_pass_total;
  std::function<size_t(void)> get_network_size;
  std::function<size_t(void)> get_network_antagonist_cnt;
  std::function<size_t(void)> get_network_sorts_per_antagonist;
  std::function<std::string(void)> get_network_passes_by_antagonist;
  std::function<std::string(void)> get_network;
  // Test stats
  std::function<size_t(void)> get_testID;
  std::function<double(void)> get_test_fitness;
  std::function<size_t(void)> get_test_pass_total;
  std::function<size_t(void)> get_test_fail_total;
  std::function<size_t(void)> get_test_size;
  std::function<std::string(void)> get_test_passes_by_antagonist;
  std::function<size_t(void)> get_test_sorts_per_antagonist;
  std::function<std::string(void)> get_test;

  emp::Ptr<emp::DataFile> sol_file;
  emp::Ptr<emp::DataFile> small_sol_file;

  // Experiment signals
  emp::Signal<void(void)> do_evaluation_sig;  ///< Trigger network/test evaluations.
  emp::Signal<void(void)> do_selection_sig;   ///< Trigger selection
  emp::Signal<void(void)> do_update_sig;      ///< Trigger world updates

  emp::Signal<void(void)> do_pop_snapshot_sig; ///< Trigger population snapshot
  emp::Signal<void(void)> do_sol_screen_sig;  ///< Trigger a screen for solutions
  emp::Signal<void(size_t)> do_small_sol_screen_sig;

  emp::Signal<void(void)> end_setup_sig;    ///< Triggered at beginning of a run.

  std::function<void(size_t)> on_lex_test_sel;
  std::function<void(size_t)> on_lex_repro;

  struct SelectionInfo {
    bool use_network_size;
  } sel_info;

  void InitConfigs(const SortingNetworkConfig & config);  ///< Initialize (localize) configuration settings.

  void InitNetworkPop_Random(); ///< Randomly initialize network population.
  void InitTestPop_Random();    ///< Randomly initialize testing population.

  void SetupEvaluation();         ///< Setup evaluation for networks and tests.
  void SetupNetworkSelection();   ///< Configure network selection.
  void SetupNetworkMutation();    ///< Configure network mutations.
  void SetupNetworkTesting();     ///< Configure sorting network testing (coevolving vs. static vs. random)
  void SetupTestSelection();      ///< Configure selection for sorting network tests. Only called if co-evolving tests.
  void SetupTestMutation();       ///< Configure mutation for sorting network tests. Only called if coevolving tests.
  void SetupDataCollection();

  void SetupNetworkStats();
  void SetupTestStats();

  void SnapshotNetworks();  ///< Output a snapshot of network population.
  void SnapshotTests();     ///< Output a snapshot of test population.

  void SetupSolutionsFile();

  /// Evaluate SortingNetworkOrg network against SortingTestOrg test,
  /// return number of passes.
  size_t EvaluateNetworkOrg(const SortingNetworkOrg & network, const SortingTestOrg & test) const;

public:

  SortingNetworkExperiment() 
    : setup(false), update(0), 
      network_cohorts(), test_cohorts(),
      network_mutator(), test_mutator(),
      curIDs(0, 0)
    { ; }

  ~SortingNetworkExperiment() {
    if (setup) {
      network_cross_binomial.Delete();

      #ifndef EMSCRIPTEN
      sol_file.Delete();
      small_sol_file.Delete();
      #endif 

      network_world.Delete();
      test_world.Delete();
      random.Delete();
    }
  }

  /// Configure the experiment
  void Setup(const SortingNetworkConfig & config);

  /// Run the experiment start->finish
  /// - Initialize population
  /// - For each generation: RunStep()
  void Run();     

  /// Progress experiment by single time step (generation):
  /// - (1) evaluation population(s)
  /// - (2) select who gets to reproduce
  /// - (3) update world(s)
  void RunStep(); 

  // -- Accessors --
  emp::Random & GetRandom() { return *random; }
  network_world_t & GetNetworkWorld() { return *network_world; }
  test_world_t & GetTestWorld() { return *test_world; }

};

/// Setup is, in theory (but untested), able to be called multiple times.
void SortingNetworkExperiment::Setup(const SortingNetworkConfig & config) {
  std::cout << "Running SortingNetworkExperiment setup." << std::endl;
  // Initialize localized configs
  InitConfigs(config);

  if (!setup) {
    // Create a random number generator
    random = emp::NewPtr<emp::Random>(SEED);
    // Create the world(s)
    network_world = emp::NewPtr<network_world_t>(*random, "Network World");
    test_world = emp::NewPtr<test_world_t>(*random, "Sorting Test World");
  } else {
    network_world->Reset();
    test_world->Reset();
    lexicase_network_fit_set.clear();
    lexicase_test_fit_set.clear();
    do_evaluation_sig.Clear();
    do_selection_sig.Clear();
    do_update_sig.Clear();
    do_pop_snapshot_sig.Clear();
    do_sol_screen_sig.Clear();
    do_small_sol_screen_sig.Clear();
    end_setup_sig.Clear();

    network_cross_binomial.Delete();
    #ifndef EMSCRIPTEN
    sol_file.Delete();
    small_sol_file.Delete();
    #endif 
  }

  network_world->SetPopStruct_Mixed(true);
  test_world->SetPopStruct_Mixed(true);

  network_pop_ids.Init(NETWORK_POP_SIZE);
  network_cross_binomial = emp::NewPtr<emp::Binomial>(PER_ORG_CROSSOVER, NETWORK_POP_SIZE);

  complete_test_set.Generate(SORT_SIZE);
  smallest_known_sol_size = MAX_NETWORK_SIZE + 1;

  // Add network world updates
  do_update_sig.AddAction([this]() {
    std::cout << "Update: " << update << ", ";
    std::cout << "best-network (size=" << network_world->GetOrg(dominant_network_id).GetSize() << "): " << network_world->CalcFitnessID(dominant_network_id) << ", ";
    // std::cout << "best-test: " << test_world->CalcFitnessID(dominant_test_id) << std::endl;
    std::cout << "solution found? " << (smallest_known_sol_size < (MAX_NETWORK_SIZE + 1)) << "; Solution size: " << smallest_known_sol_size << std::endl;
    
    if (update % SNAPSHOT_INTERVAL == 0) do_pop_snapshot_sig.Trigger();
    // if (update % SOLUTION_SCREEN_INTERVAL == 0) do_sol_screen_sig.Trigger();

    network_world->Update();
    network_world->ClearCache();

    test_world->Update();
    test_world->ClearCache();

  });

  end_setup_sig.AddAction([this]() {
    // Initialize populations
    InitNetworkPop_Random();
    InitTestPop_Random();
  });
  
  SetupEvaluation();       // Setup population evaluation
  SetupNetworkSelection(); // Setup network selection
  SetupNetworkMutation();  // Setup network mutations
  SetupNetworkTesting();   // Setup testing mode 
  #ifndef EMSCRIPTEN
  SetupDataCollection();
  #endif 

  // Setup fitness function
  network_world->SetFitFun([this](network_org_t & network) {
    if (network.GetPhenotype().num_passes == MAX_PASSES) {
      return (double)network.GetPhenotype().num_passes + ((double)(MAX_NETWORK_SIZE - network.GetSize())/(double)MAX_NETWORK_SIZE);
    } else {
      return (double)network.GetPhenotype().num_passes;
    }
  });

  test_world->SetFitFun([](test_org_t & test) {
    return (double)test.GetPhenotype().num_fails;
  });
  
  end_setup_sig.Trigger();
  setup = true;
  std::cout << "Done with experiment setup." << std::endl;
}

void SortingNetworkExperiment::SetupDataCollection() {
  std::cout << "Setting up data collection!" << std::endl;
  // Make a data directory
  mkdir(DATA_DIRECTORY.c_str(), ACCESSPERMS);
  if (DATA_DIRECTORY.back() != '/') DATA_DIRECTORY += '/';  
  
  // Setup stats getters
  SetupNetworkStats();
  SetupTestStats();

  // Setup solutions file
  SetupSolutionsFile();

  // Setup fitness files
  // SetupFitnessFile
  auto & network_fit_file = network_world->SetupFitnessFile(DATA_DIRECTORY + "network_stats.csv", false);
  network_fit_file.SetTimingRepeat(AGGREGATE_STATS_INTERVAL);
  // Add shannon entropy to fitness 
  network_fit_file.template AddFun<double>([this]() -> double {
    std::function<network_genome_t(emp::Ptr<network_org_t>)> get_network_org_genome = [](emp::Ptr<network_org_t> org) -> network_genome_t { return org->GetGenome(); };
    auto vec = emp::ApplyFunction(get_network_org_genome, network_world->GetFullPop());
    return emp::ShannonEntropy(vec);
  }, "diversity", "Shannon diversity of genotypes in population.");
  network_fit_file.PrintHeaderKeys();

  auto & test_fit_file = test_world->SetupFitnessFile(DATA_DIRECTORY + "test_stats.csv", false);
  test_fit_file.SetTimingRepeat(AGGREGATE_STATS_INTERVAL);
  test_fit_file.template AddFun<double>([this]() -> double {
    std::function<test_genome_t(emp::Ptr<test_org_t>)> get_test_org_genome = [](emp::Ptr<test_org_t> org) -> test_genome_t { return org->GetGenome(); };
    auto vec = emp::ApplyFunction(get_test_org_genome, test_world->GetFullPop());
    return emp::ShannonEntropy(vec);
  }, "diversity", "Shannon diversity of genotypes in population.");
  test_fit_file.PrintHeaderKeys();

  // Setup systematics managers
  // - Network systematics manager
  // network_genotypic_systematics = emp::NewPtr<emp::Systematics<network_org_t, network_genome_t>>([](const network_org_t & org) { return org.GetGenome(); });
  // network_genotypic_systematics->AddEvolutionaryDistinctivenessDataNode();
  // network_genotypic_systematics->AddPairwiseDistanceDataNode();
  // network_genotypic_systematics->AddPhylogeneticDiversityDataNode();
  // network_world->AddSystematics(network_genotypic_systematics, "network_genotype");
  // auto & network_gen_sys_file = network_world->SetupSystematicsFile("network_genotype", DATA_DIRECTORY + "/network_gen_sys.csv", false);
  // network_gen_sys_file.SetTimingRepeat(AGGREGATE_STATS_INTERVAL);
  // Functions to add:
  // network_gen_sys_file.AddStats(*network_genotypic_systematics->GetDataNode("evolutionary_distinctiveness") , "evolutionary_distinctiveness", "evolutionary distinctiveness for a single update", true, true);
  // network_gen_sys_file.AddStats(*network_genotypic_systematics->GetDataNode("pairwise_distances"), "pairwise_distance", "pairwise distance for a single update", true, true);
  // - GetPhylogeneticDiversity
  // network_gen_sys_file.AddCurrent(*network_genotypic_systematics->GetDataNode("phylogenetic_diversity"), "current_phylogenetic_diversity", "current phylogenetic_diversity", true, true);
  // - GetTreeSize
  // network_gen_sys_file.template AddFun<size_t>([this]() { return network_genotypic_systematics->GetTreeSize(); }, "tree_size", "Phylogenetic tree size");
  // network_gen_sys_file.PrintHeaderKeys();
  // Add function(s) to program systematics snapshot
  // using network_taxon_t = typename emp::Systematics<network_org_t, network_genome_t>::taxon_t;
  // network_genotypic_systematics->AddSnapshotFun([](const network_taxon_t & t) {
  //   std::ostringstream stream;
  //   stream << "\"";
  //   t.GetInfo().Print(stream, ",");
  //   stream << "\"";
  //   return stream.str();
  // }, "network", "sorting network");

  if (COLLECT_TEST_PHYLOGENIES) {
    std::cout << "Collecting test phylogenies!" << std::endl;
    test_genotypic_systematics = emp::NewPtr<emp::Systematics<test_org_t, test_genome_t>>([](const test_org_t & org) { return org.GetGenome(); });
    test_genotypic_systematics->AddEvolutionaryDistinctivenessDataNode();
    test_genotypic_systematics->AddPairwiseDistanceDataNode();
    test_genotypic_systematics->AddPhylogeneticDiversityDataNode();
    test_world->AddSystematics(test_genotypic_systematics, "test_genotype");
    auto & test_gen_sys_file = test_world->SetupSystematicsFile("test_genotype", DATA_DIRECTORY + "/test_gen_sys.csv", false);
    test_gen_sys_file.SetTimingRepeat(AGGREGATE_STATS_INTERVAL);
    // Functions to add:
    test_gen_sys_file.AddStats(*test_genotypic_systematics->GetDataNode("evolutionary_distinctiveness") , "evolutionary_distinctiveness", "evolutionary distinctiveness for a single update", true, true);
    test_gen_sys_file.AddStats(*test_genotypic_systematics->GetDataNode("pairwise_distances"), "pairwise_distance", "pairwise distance for a single update", true, true);
    // - GetPhylogeneticDiversity
    test_gen_sys_file.AddCurrent(*test_genotypic_systematics->GetDataNode("phylogenetic_diversity"), "current_phylogenetic_diversity", "current phylogenetic_diversity", true, true);
    // - GetTreeSize
    test_gen_sys_file.template AddFun<size_t>([this]() { return test_genotypic_systematics->GetTreeSize(); }, "tree_size", "Phylogenetic tree size");
    test_gen_sys_file.PrintHeaderKeys();
    // Add function(s) to program systematics snapshot
    using test_taxon_t = typename emp::Systematics<test_org_t, test_genome_t>::taxon_t;
    test_genotypic_systematics->AddSnapshotFun([](const test_taxon_t & t) {
      std::ostringstream stream;
      // t.GetInfo().PrintCSVEntry(stream);
      const test_genome_t & test_genome = t.GetInfo();
      stream << "\"[";
      for (size_t i = 0; i < test_genome.test_set.size(); ++i) {
        if (i) stream << ",";
        test_genome.test_set[i].Print(stream);
      }
      stream << "]\"";

      return stream.str();
    }, "test", "Test org");
  }


  
  // Setup network/test snapshots
  do_pop_snapshot_sig.AddAction([this]() { 
    SnapshotNetworks(); 
    SnapshotTests();
    // network_genotypic_systematics->Snapshot(DATA_DIRECTORY + "pop_" + emp::to_string(network_world->GetUpdate()) + "/network_phylogeny_" + emp::to_string((int)network_world->GetUpdate()) + ".csv");
    if (TEST_MODE == (size_t)TEST_MODES::COEVOLVE && COLLECT_TEST_PHYLOGENIES) {
      test_genotypic_systematics->Snapshot(DATA_DIRECTORY + "pop_" + emp::to_string(network_world->GetUpdate()) + "/test_phylogeny_" + emp::to_string((int)network_world->GetUpdate()) + ".csv");
    }
  });
}

void SortingNetworkExperiment::SetupNetworkStats() {
  
  get_networkID = [this]() { return curIDs.networkID; };
  
  get_network_fitness = [this]() { return network_world->CalcFitnessID(curIDs.networkID); };
  
  get_network_true_correct = [this]() {
    return complete_test_set.EvaluateAll(network_world->GetOrg(curIDs.networkID).GetGenome());
  };

  get_network_possible_correct = [this]() {
    return emp::Pow2(SORT_SIZE);
  };

  get_network_sample_correct = [this]() {
    const size_t n = complete_test_set.EvaluateN(network_world->GetOrg(curIDs.networkID).GetGenome(), CORRECTNESS_SAMPLE_SIZE);
    // if (n == CORRECTNESS_SAMPLE_SIZE) networkIDs_to_test_for_correctness.emplace_back(curIDs.networkID);
    return n;
  };
  
  get_network_possible_sample_correct = [this]() {
    return CORRECTNESS_SAMPLE_SIZE;
  };

  get_network_pass_total = [this]() { 
    return network_world->GetOrg(curIDs.networkID).GetPhenotype().num_passes; 
  };
  
  get_network_size = [this]() {
    return network_world->GetOrg(curIDs.networkID).GetSize();
  };
  
  get_network_antagonist_cnt = [this]() {
    return network_world->GetOrg(curIDs.networkID).GetPhenotype().test_results.size();
  };
  
  get_network_sorts_per_antagonist = [this]() {
    return SORTS_PER_TEST;
  };
  
  get_network_passes_by_antagonist = [this]() {
    network_org_t & network = network_world->GetOrg(curIDs.networkID);
    std::string scores = "\"[";
    for (size_t i = 0; i < network.GetPhenotype().test_results.size(); ++i) {
      if (i) scores += ",";
      scores += emp::to_string(network.GetPhenotype().test_results[i]);
    }
    scores += "]\"";
    return scores;
  };
  
  get_network = [this]() {
    std::ostringstream stream;
    stream << "\"";
    network_world->GetOrg(curIDs.networkID).GetGenome().Print(stream, ",");
    stream << "\"";
    return stream.str();
  };
}

void SortingNetworkExperiment::SetupTestStats() {
  get_testID = [this]() {
    return curIDs.testID;
  };
  get_test_fitness = [this]() {
    return test_world->CalcFitnessID(curIDs.testID);
  };
  get_test_pass_total = [this]() {
    return test_world->GetOrg(curIDs.testID).GetPhenotype().num_passes;
  };
  get_test_fail_total = [this]() {
    return test_world->GetOrg(curIDs.testID).GetPhenotype().num_fails;
  };
  get_test_size = [this]() {
    return test_world->GetOrg(curIDs.testID).GetTestSize();
  };
  get_test_passes_by_antagonist = [this]() {
    test_org_t & test = test_world->GetOrg(curIDs.testID);
    std::string scores = "\"[";
    for (size_t i = 0; i < test.GetPhenotype().test_results.size(); ++i) {
      if (i) scores += ",";
      scores += emp::to_string(test.GetPhenotype().test_results[i]);
    }
    scores += "]\"";
    return scores;
  };
  get_test_sorts_per_antagonist = [this]() {
    return test_world->GetOrg(curIDs.testID).GetNumTests();
  };
  get_test = [this]() {
    std::ostringstream stream;
    stream << "\"";
    test_world->GetOrg(curIDs.testID).PrintMin(stream);
    stream << "\"";
    return stream.str();
  };
}

void SortingNetworkExperiment::SetupEvaluation() {
  // Setup population evaluation.
  const bool cohort_eval = SELECTION_MODE == SELECTION_METHODS::COHORT_LEXICASE;
  if (cohort_eval) {  // We're evaluating networks with tests in cohorts.
    // Make sure settings abide by expectations.
    emp_assert(NETWORK_POP_SIZE == TEST_POP_SIZE, "Network and test population sizes must match in random cohort evaluation mode.");
    emp_assert(NETWORK_POP_SIZE % COHORT_SIZE == 0, "Population sizes must be evenly divisible by cohort size in random cohort evaluation mode.");
    // Setup cohorts
    network_cohorts.Init(NETWORK_POP_SIZE, COHORT_SIZE);
    test_cohorts.Init(TEST_POP_SIZE, COHORT_SIZE);
    MAX_PASSES = SORTS_PER_TEST * COHORT_SIZE;
    // Setup world to reset phenotypes on organism placement      
    network_world->OnPlacement([this](size_t pos){ network_world->GetOrg(pos).GetPhenotype().Reset(COHORT_SIZE); });
    test_world->OnPlacement([this](size_t pos){ test_world->GetOrg(pos).GetPhenotype().Reset(COHORT_SIZE); });
    // What should happen on evaluation?
    do_evaluation_sig.AddAction([this]() {
      // Randomize the cohorts.
      network_cohorts.Randomize(*random);
      test_cohorts.Randomize(*random);
      // For each cohort, evaluate all networks in cohort against all tests in cohort.
      for (size_t cID = 0; cID < network_cohorts.GetNumCohorts(); ++cID) {
        for (size_t nID = 0; nID < COHORT_SIZE; ++nID) {
          network_org_t & network = network_world->GetOrg(network_cohorts.GetWorldID(cID, nID));
          for (size_t tID = 0; tID < COHORT_SIZE; ++tID) {
            test_org_t & test = test_world->GetOrg(test_cohorts.GetWorldID(cID, tID));
            // Evaluate network, nID, on test, tID.
            const size_t passes = EvaluateNetworkOrg(network, test);
            network.GetPhenotype().test_results[tID] = passes;
            test.GetPhenotype().test_results[nID] = passes;
          }
        }
      }
    });
  } else { // Evaluate all networks on all tests.
    // Setup world to reset phenotypes on organism placement (each phenotype will need spots to store scores against antagonist's pop size)
    network_world->OnPlacement([this](size_t pos){ network_world->GetOrg(pos).GetPhenotype().Reset(TEST_POP_SIZE); });
    test_world->OnPlacement([this](size_t pos){ test_world->GetOrg(pos).GetPhenotype().Reset(NETWORK_POP_SIZE); });
    MAX_PASSES = SORTS_PER_TEST * TEST_POP_SIZE;
    // What should happen on evaluation?
    do_evaluation_sig.AddAction([this]() {
      for (size_t nID = 0; nID < network_world->GetSize(); ++nID) {
        network_org_t & network = network_world->GetOrg(nID);
        for (size_t tID = 0; tID < test_world->GetSize(); ++tID) {
          test_org_t & test = test_world->GetOrg(tID);
          // Evaluate network, nID, on test, tID.
          const size_t passes = EvaluateNetworkOrg(network, test);
          network.GetPhenotype().test_results[tID] = passes;

          test.GetPhenotype().test_results[nID] = passes;
        }
      }
    });
  }

  // Add evaluation action to calculate num_passes for tests and networks.
  do_evaluation_sig.AddAction([this]() {
    // Sum pass totals for networks.
    dominant_network_id = 0;
    double cur_best = 0;
    complete_test_set.swap_id = 0;
    for (size_t nID = 0; nID < network_world->GetSize(); ++nID) {
      if (!network_world->IsOccupied(nID)) continue;
      network_org_t & network = network_world->GetOrg(nID);
      const size_t pass_total = emp::Sum(network.GetPhenotype().test_results);
      network.GetPhenotype().num_passes = pass_total;
      // Is this highest fitness network this generation?
      if (pass_total > cur_best || nID == 0) {
        dominant_network_id = nID;
        cur_best = pass_total;
      }
      // At this point, network has been evaluated against all tests. Screen
      // for possible smallest known solution.
      if (pass_total == MAX_PASSES && network.GetSize() < smallest_known_sol_size) {
        do_small_sol_screen_sig.Trigger(nID);
      }
    }
    // Sum pass totals for tests.
    dominant_test_id = 0;
    cur_best = 0;
    for (size_t tID = 0; tID < test_world->GetSize(); ++tID) {
      if (!test_world->IsOccupied(tID)) continue;
      test_org_t & test = test_world->GetOrg(tID);
      const size_t pass_total = emp::Sum(test.GetPhenotype().test_results);
      const size_t fail_total = (test.GetPhenotype().test_results.size() * SORTS_PER_TEST) - pass_total;
      test.GetPhenotype().num_passes = pass_total;
      test.GetPhenotype().num_fails = fail_total;
      if (fail_total > cur_best || tID == 0) { 
        dominant_test_id = tID;
        cur_best = fail_total;
      }
    }
  });

}

void SortingNetworkExperiment::SetupNetworkSelection() {
  // Configure network selection.

  sel_info.use_network_size = false;

  switch (SELECTION_MODE) {
    
    case (size_t)SELECTION_METHODS::COHORT_LEXICASE: {
      on_lex_test_sel = [this](size_t fit_id) {
        // sel_info.use_network_size = random->P(0.05);
        sel_info.use_network_size = false;
      };
      // Setup network fit funs
      // - 1 function for every cohort member
      for (size_t i = 0; i < COHORT_SIZE; ++i) {
        lexicase_network_fit_set.push_back([this, i](network_org_t & network) {
          double score = network.GetPhenotype().test_results[i];
          // if (sel_info.use_network_size && score == SORTS_PER_TEST) {
          //   score += ((double)(MAX_NETWORK_SIZE - network.GetSize())/(double)MAX_NETWORK_SIZE);
          // }
          return score;
        });
      }
      // - 1 test case for being small
      lexicase_network_fit_set.push_back([this](network_org_t & network) {
        if (network.GetPhenotype().num_passes == (COHORT_SIZE * SORTS_PER_TEST)) return (double)MAX_NETWORK_SIZE - (double)network.GetSize();
        return 0.0;
      });
      // Add selection action
      emp_assert(COHORT_SIZE * network_cohorts.GetNumCohorts() == NETWORK_POP_SIZE);
      do_selection_sig.AddAction([this]() {
        // For each cohort, run selection
        for (size_t cID = 0; cID < network_cohorts.GetNumCohorts(); ++cID) {
          emp::CohortLexicaseSelect_NAIVE(*network_world, 
                                          lexicase_network_fit_set,
                                          network_cohorts.GetCohort(cID),
                                          COHORT_SIZE,
                                          COHORTLEX_MAX_FUNS,
                                          on_lex_test_sel);
        }
      });
      break;
    }
    
    case (size_t)SELECTION_METHODS::LEXICASE: {
      on_lex_test_sel = [this](size_t fit_id) {
        // std::cout << "FitID: " << fit_id << std::endl;
        // sel_info.use_network_size = random->P(0.05);
        sel_info.use_network_size = false;
      };
      // For lexicase selection, one function for every test organism.
      for (size_t i = 0; i < TEST_POP_SIZE; ++i) {
        lexicase_network_fit_set.push_back([this, i](network_org_t & network) {
          
          double score = network.GetPhenotype().test_results[i];
          // if (sel_info.use_network_size && score == SORTS_PER_TEST) {
          //   score += ((double)(MAX_NETWORK_SIZE - network.GetSize())/(double)MAX_NETWORK_SIZE);
          // }
          // std::cout << "Use size? " << sel_info.use_network_size << "  score=" << score << std::endl;
          return score;
        });
      }
      lexicase_network_fit_set.push_back([this](network_org_t & network) {
        if (network.GetPhenotype().num_passes == (TEST_POP_SIZE * SORTS_PER_TEST)) return (double)MAX_NETWORK_SIZE - (double)network.GetSize();
        return 0.0;
      });
      do_selection_sig.AddAction([this]() {
        emp::LexicaseSelect_NAIVE(*network_world,
                            lexicase_network_fit_set,
                            NETWORK_POP_SIZE,
                            LEX_MAX_FUNS,
                            on_lex_test_sel);
      });
      break;
    }
    
    case (size_t)SELECTION_METHODS::TOURNAMENT: {
      do_selection_sig.AddAction([this]() {
        emp::TournamentSelect(*network_world, TOURNAMENT_SIZE, NETWORK_POP_SIZE);
      });
      break;
    }

    default: {
      std::cout << "Unrecognized SELECTION_MODE(" << SELECTION_MODE << "). Exiting..." << std::endl;
      exit(-1);
    }
  }

  // Do crossover
  // GetRandBinomial(const double n, const double p)
  switch (NETWORK_CROSSOVER_MODE) {
    case NETWORK_CROSSOVER_MODES::NONE: {
      break;
    }
    case NETWORK_CROSSOVER_MODES::SINGLE_PT: {
      do_selection_sig.AddAction([this]() {
        network_pop_ids.Shuffle(*random);
        // const uint32_t num_crosses = random->GetRandBinomial((double)NETWORK_POP_SIZE, (double)PER_ORG_CROSSOVER);
        const size_t num_crosses = network_cross_binomial->PickRandom(*random);
        for (size_t i = 0; (int)i < (int)num_crosses-1; ++i) {
          emp_assert(i+1 < NETWORK_POP_SIZE, NETWORK_POP_SIZE, num_crosses);
          network_mutator.Crossover1Pt(*random,
                                      network_world->GetNextOrg(network_pop_ids.popIDs[i]).GetGenome(), 
                                      network_world->GetNextOrg(network_pop_ids.popIDs[i+1]).GetGenome());
        }
      });
      break;
    }
    case NETWORK_CROSSOVER_MODES::TWO_PT: {
      do_selection_sig.AddAction([this]() {
        network_pop_ids.Shuffle(*random);
        const size_t num_crosses = network_cross_binomial->PickRandom(*random);
        for (size_t i = 0; (int)i < (int)num_crosses-1; ++i) {
          emp_assert(i+1 < NETWORK_POP_SIZE, NETWORK_POP_SIZE, num_crosses);
          network_mutator.Crossover2Pt(*random,
                                      network_world->GetNextOrg(network_pop_ids.popIDs[i]).GetGenome(), 
                                      network_world->GetNextOrg(network_pop_ids.popIDs[i+1]).GetGenome());
        }
      });
      break;
    }
    default: {
      std::cout << "Unrecognized crossover mode (" << NETWORK_CROSSOVER_MODE << "). Exiting..." << std::endl;
      exit(-1);
      break;
    }
  }
}

void SortingNetworkExperiment::SetupNetworkTesting() {
  
  switch (TEST_MODE) {
    case (size_t)TEST_MODES::COEVOLVE: {
      std::cout << "Setting up network testing - CO-EVOLUTION MODE" << std::endl;
      // Sorting tests coevolve with sorting networks.
      // - Evaluation already done. 
      // - Setup selection
      SetupTestSelection();
      SetupTestMutation();
      break;
    }
    case (size_t)TEST_MODES::DRIFT: {
      std::cout << "Setting up network testing - DRIFT MODE" << std::endl;
      do_selection_sig.AddAction([this]() {
        emp::RandomSelect(*test_world, TEST_POP_SIZE);
      });
      SetupTestMutation();
      break;
    }
    case (size_t)TEST_MODES::STATIC: {
      std::cout << "Setting up network testing - STATIC MODE" << std::endl;
      // Sorting tests are static over time.
      test_world->SetPopStruct_Mixed(false);
      // ... do nothing ...
      break;
    }
    case (size_t)TEST_MODES::RANDOM: {
      std::cout << "Setting up network testing - RANDOM MODE" << std::endl;
      test_world->SetPopStruct_Mixed(false);
      // On world update, randomize test cases.
      do_update_sig.AddAction([this]() {
        for (size_t tID = 0; tID < test_world->GetSize(); ++tID) {
          if (!test_world->IsOccupied(tID)) continue;
          test_world->GetOrg(tID).GetGenome().Randomize(*random);
        }
      });
      break;
    }
    default: {
      std::cout << "Unrecognized TEST_MODE (" << TEST_MODE << "). Exiting..." << std::endl;
      exit(-1);
    }
  }
}

void SortingNetworkExperiment::SetupNetworkMutation() {
  // Setup network mutation functions
  network_mutator.MAX_NETWORK_SIZE = MAX_NETWORK_SIZE;
  network_mutator.MIN_NETWORK_SIZE = MIN_NETWORK_SIZE;
  network_mutator.SORT_SEQ_SIZE = SORT_SIZE;
  network_mutator.PER_INDEX_SUB = PER_INDEX_SUB;
  network_mutator.PER_PAIR_DUP = PER_PAIR_DUP;
  network_mutator.PER_PAIR_INS = PER_PAIR_INS;
  network_mutator.PER_PAIR_DEL = PER_PAIR_DEL;
  network_mutator.PER_PAIR_SWAP = PER_PAIR_SWAP;

  if (PER_ORG_MUTATION == 1.0) {
    network_world->SetMutFun([this](network_org_t & network, emp::Random & rnd) {
      return network_mutator.Mutate(rnd, network.GetGenome());
    });
  } else {
    network_world->SetMutFun([this](network_org_t & network, emp::Random & rnd) {
      if (rnd.P(PER_ORG_MUTATION)) {
        return (double)network_mutator.Mutate(rnd, network.GetGenome());
      } else {
        return 0.0;
      }
    });
  }

  end_setup_sig.AddAction([this]() {
    network_world->SetAutoMutate(); // After we've initialized populations, turn auto-mutate on.
  });
}

void SortingNetworkExperiment::SetupTestSelection() {
  std::cout << "Setting up test selection!" << std::endl;
  // Configure test selection scheme (will match network selection mode).
  switch (SELECTION_MODE) {
    case (size_t)SELECTION_METHODS::COHORT_LEXICASE: {
      // Setup test fit funs
      // - 1 function for every cohort member
      if (DISCRIMINATORY_LEXICASE_TESTS) {
        std::cout << "  Configuring tests to be DISCRIMINATORY" << std::endl;
        for (size_t i = 0; i < COHORT_SIZE; ++i) {
          lexicase_test_fit_set.push_back([i, this](test_org_t & test) {
            if (test.GetPhenotype().test_results[i] == SORTS_PER_TEST) {
              return (double)test.GetPhenotype().num_fails;
            } else if (test.GetPhenotype().num_passes == 0) {
              return 0.5;
            } else {
              return 0.0;
            } 
          });
        }
      } else {
        for (size_t i = 0; i < COHORT_SIZE; ++i) {
          lexicase_test_fit_set.push_back([i](test_org_t & test) {
            return test.GetNumTests() - test.GetPhenotype().test_results[i];
          });
        }
      }
      // Add selection action.
      emp_assert(COHORT_SIZE * test_cohorts.GetNumCohorts() == TEST_POP_SIZE);
      do_selection_sig.AddAction([this]() {
        // For each cohort, run selection.
        for (size_t cID = 0; cID < test_cohorts.GetNumCohorts(); ++cID) {
          emp::CohortLexicaseSelect(*test_world,
                                    lexicase_test_fit_set,
                                    test_cohorts.GetCohort(cID),
                                    COHORT_SIZE,
                                    COHORTLEX_MAX_FUNS);
        }
      });
      break;
    }
    case (size_t)SELECTION_METHODS::LEXICASE: {
      // For lexicase selection, one function for every network organism.
      if (DISCRIMINATORY_LEXICASE_TESTS) {
        std::cout << "  Configuring tests to be DISCRIMINATORY" << std::endl;
        for (size_t i = 0; i < NETWORK_POP_SIZE; ++i) {
          lexicase_test_fit_set.push_back([i, this](test_org_t & test) {
            if (test.GetPhenotype().test_results[i] == SORTS_PER_TEST) {
              return (double)test.GetPhenotype().num_fails;
            } else if (test.GetPhenotype().num_passes == 0) {
              return 0.5;
            } else {
              return 0.0;
            } 
          });
        }
      } else {
        for (size_t i = 0; i < NETWORK_POP_SIZE; ++i) {
          lexicase_test_fit_set.push_back([i](test_org_t & test) {
            return test.GetNumTests() - test.GetPhenotype().test_results[i]; // Incorrect!
          });
        }
      }

      do_selection_sig.AddAction([this]() {
        emp::LexicaseSelect_NAIVE(*test_world, lexicase_test_fit_set, TEST_POP_SIZE, LEX_MAX_FUNS);
      });
      break;
    }
    case (size_t)SELECTION_METHODS::TOURNAMENT: {
      do_selection_sig.AddAction([this]() {
        emp::TournamentSelect(*test_world, TOURNAMENT_SIZE, TEST_POP_SIZE);
      });
      break;
    }
    default: {
      std::cout << "Unrecognized SELECTION_MODE (" << SELECTION_MODE << "). Exiting..." << std::endl;
      exit(-1);
    }
  }
}

void SortingNetworkExperiment::SetupTestMutation() {
  std::cout << "Setting up test mutation!" << std::endl;
  test_mutator.bit_mode = true;
  test_mutator.MAX_VALUE = 1;
  test_mutator.MIN_VALUE = 0;
  test_mutator.PER_SITE_SUB = PER_SITE_SUB;
  test_mutator.PER_SEQ_INVERSION = PER_SEQ_INVERSION;
  test_mutator.PER_SEQ_RANDOMIZE = PER_SEQ_RANDOMIZE;
  test_world->SetMutFun([this](test_org_t & test, emp::Random & rnd) {
    return test_mutator.Mutate(rnd, test.GetGenome());
  });
  end_setup_sig.AddAction([this]() {
    test_world->SetAutoMutate();    // After populations have been initialized, set world to auto mutate.
  });
}

void SortingNetworkExperiment::Run() {
  for (update = 0; update <= GENERATIONS; ++update) {
    RunStep();
  }
}

void SortingNetworkExperiment::RunStep() {
  do_evaluation_sig.Trigger();
  do_selection_sig.Trigger();
  do_update_sig.Trigger();
}

void SortingNetworkExperiment::InitConfigs(const SortingNetworkConfig & config) {
  // Default group
  SEED = config.SEED();
  GENERATIONS = config.GENERATIONS();
  NETWORK_POP_SIZE = config.NETWORK_POP_SIZE();
  TEST_POP_SIZE = config.TEST_POP_SIZE();
  TEST_MODE = config.TEST_MODE();
  COHORT_SIZE = config.COHORT_SIZE();

  SELECTION_MODE = config.SELECTION_MODE();
  LEX_MAX_FUNS = config.LEX_MAX_FUNS();
  COHORTLEX_MAX_FUNS = config.COHORTLEX_MAX_FUNS();
  TOURNAMENT_SIZE = config.TOURNAMENT_SIZE();
  DISCRIMINATORY_LEXICASE_TESTS = config.DISCRIMINATORY_LEXICASE_TESTS();
  
  MAX_NETWORK_SIZE = config.MAX_NETWORK_SIZE();
  MIN_NETWORK_SIZE = config.MIN_NETWORK_SIZE();

  PER_INDEX_SUB = config.PER_INDEX_SUB();
  PER_PAIR_DUP = config.PER_PAIR_DUP();
  PER_PAIR_INS = config.PER_PAIR_INS();
  PER_PAIR_DEL = config.PER_PAIR_DEL();
  PER_PAIR_SWAP = config.PER_PAIR_SWAP();
  NETWORK_CROSSOVER_MODE = config.NETWORK_CROSSOVER_MODE();
  PER_ORG_CROSSOVER = config.PER_ORG_CROSSOVER();
  PER_ORG_MUTATION = config.PER_ORG_MUTATION();

  SORT_SIZE = config.SORT_SIZE();
  SORTS_PER_TEST = config.SORTS_PER_TEST();

  PER_SITE_SUB = config.PER_SITE_SUB();
  PER_SEQ_INVERSION = config.PER_SEQ_INVERSION();
  PER_SEQ_RANDOMIZE = config.PER_SEQ_RANDOMIZE();

  DATA_DIRECTORY = config.DATA_DIRECTORY();
  SNAPSHOT_INTERVAL = config.SNAPSHOT_INTERVAL();
  DOMINANT_STATS_INTERVAL = config.DOMINANT_STATS_INTERVAL();
  AGGREGATE_STATS_INTERVAL = config.AGGREGATE_STATS_INTERVAL();
  CORRECTNESS_SAMPLE_SIZE = config.CORRECTNESS_SAMPLE_SIZE();
  SOLUTION_SCREEN_INTERVAL = config.SOLUTION_SCREEN_INTERVAL();
  COLLECT_TEST_PHYLOGENIES = config.COLLECT_TEST_PHYLOGENIES();

}

void SortingNetworkExperiment::InitNetworkPop_Random() {
  std::cout << "Randomly initializing sorting network population...";
  // Inject random networks into network world up to population size.
  for (size_t i = 0; i < NETWORK_POP_SIZE; ++i) {
    network_world->Inject(network_genome_t(*random, SORT_SIZE, MIN_NETWORK_SIZE, MAX_NETWORK_SIZE), 1);
  }
  std::cout << " Done." << std::endl;
}

void SortingNetworkExperiment::InitTestPop_Random() {
  std::cout << "Random initializing sorting test population...";
  for (size_t i = 0; i < TEST_POP_SIZE; ++i) {
    test_world->Inject(test_genome_t(*random, SORT_SIZE, SORTS_PER_TEST), 1);
  }
  std::cout << " Done." << std::endl;
}

void SortingNetworkExperiment::SnapshotNetworks() {
  std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(network_world->GetUpdate());
  mkdir(snapshot_dir.c_str(), ACCESSPERMS);
  
  emp::DataFile file(snapshot_dir + "/network_pop_" + emp::to_string((int)network_world->GetUpdate()) + ".csv");

  // Add functions to data file
  // - networkID
  file.AddFun(get_networkID, "network_id", "Network ID");
  // - Fitness
  file.AddFun(get_network_fitness, "fitness", "");
  // - pass total
  file.AddFun(get_network_pass_total, "pass_total", "");
  
  file.AddFun(get_network_sample_correct, "sample_passes");
  file.AddFun(get_network_possible_sample_correct, "sample_size");

  // - network size
  file.AddFun(get_network_size, "network_size");
  // - num_antagonists
  file.AddFun(get_network_antagonist_cnt, "num_antagonists");
  // - tests per antagonist
  file.AddFun(get_network_sorts_per_antagonist, "sorts_per_antagonist");
  // - antagonist scores]
  file.AddFun(get_network_passes_by_antagonist, "scores_by_antagonist");
  // - network
  file.AddFun(get_network, "network");

  // Output file headers
  file.PrintHeaderKeys();

  // For each network in the population, dump the network and anything we want to know about it.
  complete_test_set.SuffleTestIDs(*random);
  for (curIDs.networkID = 0; curIDs.networkID < network_world->GetSize(); ++curIDs.networkID) {
    if (!network_world->IsOccupied(curIDs.networkID)) continue;
    file.Update();
  }
  
}

void SortingNetworkExperiment::SnapshotTests() {
  std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(test_world->GetUpdate());
  mkdir(snapshot_dir.c_str(), ACCESSPERMS);

  emp::DataFile file(snapshot_dir + "/test_pop_" + emp::to_string((int)test_world->GetUpdate()) + ".csv");  

  file.AddFun(get_testID, "test_id");
  file.AddFun(get_test_fitness, "fitness");
  file.AddFun(get_test_pass_total, "pass_total");
  file.AddFun(get_test_fail_total, "fail_total");
  file.AddFun(get_test_sorts_per_antagonist, "sorts_per_antagonist");
  file.AddFun(get_test_passes_by_antagonist, "passes_by_antagonist");
  file.AddFun(get_test_size, "test_size");
  file.AddFun(get_test, "test");

  // Output file headers
  file.PrintHeaderKeys();

  // For each test in the population, dump the test and anything we want to know about it.
  for (curIDs.testID = 0; curIDs.testID < test_world->GetSize(); ++curIDs.testID) {
    if (!test_world->IsOccupied(curIDs.testID)) continue;
    file.Update();
  }
}

size_t SortingNetworkExperiment::EvaluateNetworkOrg(const SortingNetworkOrg & network,
                                                    const SortingTestOrg & test) const {
  size_t passes = 0;
  for (size_t i = 0; i < test.GetNumTests(); ++i) {
    passes += (size_t)(test.GetGenome().test_set[i].Evaluate(network.GetGenome()));
  }
  return passes;                                                    
}

void SortingNetworkExperiment::SetupSolutionsFile() {
  sol_file = emp::NewPtr<emp::DataFile>(DATA_DIRECTORY + "/solutions.csv");

  do_sol_screen_sig.AddAction([this]() {
    // - For each potential solution -> is_correct? -> if so, sol_file.update
    complete_test_set.swap_id = 0;
    for (curIDs.networkID = 0; curIDs.networkID < network_world->GetSize(); ++curIDs.networkID) {
      // Is network a candidate for solution-checking?
      network_org_t & network = network_world->GetOrg(curIDs.networkID);
      if (network.GetPhenotype().num_passes == MAX_PASSES) {
        bool correct = complete_test_set.Correct(network_world->GetOrg(curIDs.networkID).GetGenome());
        if (correct) sol_file->Update();
      }
    }

  });
  
  std::function<size_t(void)> get_update = [this]() { return network_world->GetUpdate(); };
  
  std::function<double(void)> get_evaluations = [this]() { // --bookmark--
    if (SELECTION_MODE == SELECTION_METHODS::COHORT_LEXICASE) {
      // evals = update * (test cohort size * program cohort size * num cohorts)
      return network_world->GetUpdate() * COHORT_SIZE * COHORT_SIZE * network_cohorts.GetNumCohorts();
    } else {
      return network_world->GetUpdate() * network_world->GetSize() * test_world->GetSize();
    }
  };


  sol_file->AddFun(get_update, "update");
  sol_file->AddFun(get_evaluations, "evaluations");

  sol_file->AddFun(get_networkID, "network_id", "Network ID");
  
  sol_file->AddFun(get_network_fitness, "fitness");
  sol_file->AddFun(get_network_pass_total, "pass_total");
  sol_file->AddFun(get_network_size, "network_size");
  sol_file->AddFun(get_network_antagonist_cnt, "num_antagonists");
  sol_file->AddFun(get_network_sorts_per_antagonist, "sorts_per_antagonist");
  sol_file->AddFun(get_network, "network");
  
  sol_file->PrintHeaderKeys();

  // Setup small sol file (will only have 1 solution per size found)
  small_sol_file = emp::NewPtr<emp::DataFile>(DATA_DIRECTORY + "/small_solutions.csv");

  do_small_sol_screen_sig.AddAction([this](size_t id) {
    curIDs.networkID = id;
    network_org_t & network = network_world->GetOrg(id);
    if (network.GetPhenotype().num_passes == MAX_PASSES && network.GetSize() < smallest_known_sol_size) {
      bool correct = complete_test_set.Correct(network.GetGenome());
      if (correct) {
        smallest_known_sol_size = network.GetSize();
        small_sol_file->Update();
      }
    }
  });

  small_sol_file->AddFun(get_update, "update");
  small_sol_file->AddFun(get_evaluations, "evaluations");
  small_sol_file->AddFun(get_networkID, "network_id", "Network ID");
  small_sol_file->AddFun(get_network_fitness, "fitness");
  small_sol_file->AddFun(get_network_pass_total, "pass_total");
  small_sol_file->AddFun(get_network_size, "network_size");
  small_sol_file->AddFun(get_network_antagonist_cnt, "num_antagonists");
  small_sol_file->AddFun(get_network_sorts_per_antagonist, "sorts_per_antagonist");
  small_sol_file->AddFun(get_network, "network");
  small_sol_file->PrintHeaderKeys();

}

#endif