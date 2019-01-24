#ifndef ALEX_SELECTION_H
#define ALEX_SELECTION_H

#include <functional>

#include "base/assert.h"
#include "base/vector.h"
#include "tools/Random.h"
#include "tools/random_utils.h"
#include "tools/vector_utils.h"

namespace emp {
  template<typename ORG> class World;

  /// NOTE: This method was a copy-pasta version of the Empirical lexicase select.
  /// - Needed to modify internals to support experiment
  template<typename ORG>
  void LexicaseSelect_NAIVE(World<ORG> & world,
                      const emp::vector< std::function<double(ORG &)> > & fit_funs,
                      size_t repro_count=1,
                      size_t max_funs=0,
                      const std::function<void(size_t)> & on_lex_test_sel = [](size_t){ ; },
                      const std::function<void(size_t)> & on_lex_repro = [](size_t) { ; })
  {
    emp_assert(world.GetSize() > 0);
    emp_assert(fit_funs.size() > 0);

    // std::cout << "=======Lexicase selection! (update: "<< world.GetUpdate() <<")=======" << std::endl;

    emp::vector<size_t> org_list;

    // list all orgs in pop
    for (size_t org_id = 0; org_id < world.GetSize(); ++org_id) {
      if (world.IsOccupied(org_id)) { org_list.emplace_back(org_id); }
    }

    emp::vector<size_t> all_orgs(org_list.size()), cur_orgs, next_orgs;

    for (size_t i = 0; i < org_list.size(); ++i) {
      all_orgs[i] = org_list[i];  
    }

    if (!max_funs) max_funs = fit_funs.size();

    // Go through a new ordering of fitness functions for each selections.
    for (size_t repro = 0; repro < repro_count; ++repro) {
      // std::cout << "--- REPRO: " << repro << " ---" << std::endl;
      // Determine the current ordering of the functions.
      emp::vector<size_t> order;

      if (max_funs == fit_funs.size()) {
        order = GetPermutation(world.GetRandom(), fit_funs.size());
      } else {
        order.resize(max_funs); // We want to limit the total numebr of tests done.
        for (auto & x : order) x = world.GetRandom().GetUInt(fit_funs.size());
      }
      
      // Step through the functions in the proper order.
      cur_orgs = all_orgs;  // Start with all of the organisms.
      int depth = -1;
      for (size_t fit_id : order) {
        // -- OnLexTest.Trigger --
        on_lex_test_sel(fit_id);
        depth++;

        double max_fit = fit_funs[fit_id](world.GetOrg(cur_orgs[0]));
        next_orgs.push_back(cur_orgs[0]);
        
        for (size_t i = 1; i < cur_orgs.size(); ++i) {

          const double cur_fit = fit_funs[fit_id](world.GetOrg(cur_orgs[i]));
          if (cur_fit > max_fit) {
            max_fit = cur_fit;             // This is a the NEW maximum fitness for this function
            next_orgs.resize(0);           // Clear out orgs with former maximum fitness
            next_orgs.push_back(cur_orgs[i]);   // Add this org as only one with new max fitness
          }
          else if (cur_fit == max_fit) {
            next_orgs.push_back(cur_orgs[i]);   // Same as cur max fitness -- save this org too.
          }
        }
        // Make next_orgs into new cur_orgs; make cur_orgs allocated space for next_orgs.
        std::swap(cur_orgs, next_orgs);
        next_orgs.resize(0);

        if (cur_orgs.size() == 1) break;  // Stop if we're down to just one organism.
      }
      // Place a random survivor (all equal) into the next generation!
      emp_assert(cur_orgs.size() > 0, cur_orgs.size(), fit_funs.size(), all_orgs.size());
      const size_t winner = world.GetRandom().GetUInt(cur_orgs.size());
      const size_t reproID = cur_orgs[winner];
      // std::cout << "Winner id = " << reproID << " Depth = " << depth << std::endl;
      // - OnSel.Trigger -
      on_lex_repro(reproID);
      world.DoBirth(world.GetGenomeAt(reproID), reproID);
    }
  }

  /// ==COHORT-LEXICASE== Selection runs through multiple fitness functions in a random order for
  /// EACH offspring produced. Only run select from population IDs specified by cohort vector.
  /// @param world The emp::World object with the organisms to be selected.
  /// @param fit_funs The set of fitness functions to shuffle for each organism reproduced.
  /// @param cohort The set of organism IDs (corresponding to world pop) to select from.
  /// @param repro_count How many rounds of repliction should we do. (default 1)
  /// @param max_funs The maximum number of fitness functions to use. (use 0 for all; default)
  template<typename ORG>
  void CohortLexicaseSelect_NAIVE(World<ORG> & world,
                            const emp::vector< std::function<double(ORG &)> > & fit_funs,
                            const emp::vector<size_t> & cohort,
                            size_t repro_count=1,
                            size_t max_funs=0,
                            const std::function<void(size_t)> & on_lex_test_sel = [](size_t){ ; },
                            const std::function<void(size_t)> & on_lex_repro = [](size_t) { ; }) 
  {
    emp_assert(world.GetSize() > 0);
    emp_assert(fit_funs.size() > 0);
    emp_assert(cohort.size() > 0);

    if (!max_funs) max_funs = fit_funs.size();

    // Reminder: these index into cohort, cohort has proper worldID.
    emp::vector<size_t> all_orgs, cur_orgs, next_orgs; 
    for (size_t i = 0; i < cohort.size(); ++i) {
      if (!world.IsOccupied(cohort[i])) continue;
      all_orgs.emplace_back(cohort[i]);
    }
    
    // Do selection!
    for (size_t repro = 0; repro < repro_count; ++repro) {
      // Get a random ordering of fitness functions.
      emp::vector<size_t> order(GetPermutation(world.GetRandom(), fit_funs.size()));
      if (max_funs < fit_funs.size()) { // If we don't use them all, toss some.
        order.resize(max_funs);
      }

      // Step through the functions in the proper order.
      cur_orgs = all_orgs; // Start with all of the organisms.
      int depth = -1;
      for (size_t fit_id : order) {
        on_lex_test_sel(fit_id);
        ++depth;
        
        double max_fit = fit_funs[fit_id](world.GetOrg(cur_orgs[0]));
        next_orgs.emplace_back(cur_orgs[0]);

        for (size_t i = 1; i < cur_orgs.size(); ++i) {
          const double cur_fit = fit_funs[fit_id](world.GetOrg(cur_orgs[i]));
          if (cur_fit > max_fit) {
            max_fit = cur_fit;      // This is a NEW maximum fitness for this function.
            next_orgs.resize(1);    // Clear out orgs with former max fitness.
            next_orgs[0] = cur_orgs[i];  // Add this org as the only one with the new max fitness.
          } else if (cur_fit == max_fit) {
            next_orgs.emplace_back(cur_orgs[i]); // Same as current max fitness, save this one too.
          }
        }
        // Make next_orgs into new cur_orgs; make cur_orgs allocated space for next_orgs.
        std::swap(cur_orgs, next_orgs);
        next_orgs.resize(0);
        if (cur_orgs.size() == 1) break; // Stop if we're down to just one organism.
      }

      // Place a random survivor (all equal) into the next generation.
      emp_assert(cur_orgs.size() > 0, cur_orgs.size(), fit_funs.size(), all_orgs.size());
      const size_t winner = world.GetRandom().GetUInt(cur_orgs.size());
      const size_t reproID = cur_orgs[winner];
      on_lex_repro(reproID);
      world.DoBirth(world.GetGenomeAt(reproID), reproID);
    }
  }

  /// ==COHORT-LEXICASE== Selection runs through multiple fitness functions in a random order for
  /// EACH offspring produced. Only run select from population IDs specified by cohort vector.
  /// @param world The emp::World object with the organisms to be selected.
  /// @param fit_funs The set of fitness functions to shuffle for each organism reproduced.
  /// @param cohort The set of organism IDs (corresponding to world pop) to select from.
  /// @param repro_count How many rounds of repliction should we do. (default 1)
  /// @param max_funs The maximum number of fitness functions to use. (use 0 for all; default)
  template<typename ORG>
  void CohortLexicaseSelect(World<ORG> & world,
                            const emp::vector< std::function<double(ORG &)> > & fit_funs,
                            const emp::vector<size_t> & cohort,
                            size_t repro_count=1,
                            size_t max_funs=0) 
  {
    emp_assert(world.GetSize() > 0);
    emp_assert(fit_funs.size() > 0);
    emp_assert(cohort.size() > 0);

    if (!max_funs) max_funs = fit_funs.size();

    // Reminder: these index into cohort, cohort has proper worldID.
    emp::vector<size_t> all_orgs(cohort.size()), cur_orgs, next_orgs; 
    for (size_t i = 0; i < cohort.size(); ++i) all_orgs[i] = i;

    // Get fitnesses for each cohort member on all fitness functions.
    emp::vector< emp::vector<double> > fitnesses(fit_funs.size());
    for (size_t fit_id = 0; fit_id < fit_funs.size(); ++fit_id) {
      fitnesses[fit_id].resize(cohort.size());
      for (size_t i = 0; i < cohort.size(); ++i) {
        emp_assert(world.IsOccupied(cohort[i]));
        fitnesses[fit_id][i] = fit_funs[fit_id](world.GetOrg(cohort[i]));
      }
    }

    // Do selection!
    for (size_t repro = 0; repro < repro_count; ++repro) {
      // Get a random ordering of fitness functions.
      emp::vector<size_t> order(GetPermutation(world.GetRandom(), fit_funs.size()));
      if (max_funs < fit_funs.size()) { // If we don't use them all, toss some.
        order.resize(max_funs);
      }

      // Step through the functions in the proper order.
      cur_orgs = all_orgs; // Start with all of the organisms.
      int depth = -1;
      for (size_t fit_id : order) {
        ++depth;
        double max_fit = fitnesses[fit_id][cur_orgs[0]];
        next_orgs.emplace_back(cur_orgs[0]);
        for (size_t org_id = 1; org_id < cur_orgs.size(); ++org_id) {
          const double cur_fit = fitnesses[fit_id][org_id];
          if (cur_fit > max_fit) {
            max_fit = cur_fit;      // This is a NEW maximum fitness for this function.
            next_orgs.resize(1);    // Clear out orgs with former max fitness.
            next_orgs[0] = org_id;  // Add this org as the only one with the new max fitness.
          } else if (cur_fit == max_fit) {
            next_orgs.emplace_back(org_id); // Same as current max fitness, save this one too.
          }
        }
        // Make next_orgs into new cur_orgs; make cur_orgs allocated space for next_orgs.
        std::swap(cur_orgs, next_orgs);
        next_orgs.resize(0);
        if (cur_orgs.size() == 1) break; // Stop if we're down to just one organism.
      }

      // Place a random survivor (all equal) into the next generation.
      emp_assert(cur_orgs.size() > 0, cur_orgs.size(), fit_funs.size(), all_orgs.size());
      const size_t winner = world.GetRandom().GetUInt(cur_orgs.size());
      const size_t reproID = cohort[cur_orgs[winner]];
      world.DoBirth(world.GetGenomeAt(reproID), reproID);
    }
  }

  /// Assumes fitness [0:BIG NUMBER] (i.e., non-negative to qualify for resource)
  template<typename ORG>
  void CohortEcoSelect_NAIVE(World<ORG> & world,
                             const emp::vector<std::function<double(ORG &)> > & extra_funs,
                             const emp::vector<double> & pool_sizes,
                             const emp::vector<size_t> & cohort,
                             size_t t_size,
                             size_t tourny_count=1) 
  {
    emp_assert(world.GetFitFun(), "Must define a base fitness function");
    // TODO - more asserts for safety!
    // emp_assert(world.GetSize() > 0);
    // emp_assert(t_size > 0 && t_size <= world.GetSize(), t_size, world.GetSize());

    // Setup info to track fitnesses
    emp::vector<double> base_fitness(cohort.size());
    emp::vector< emp::vector<double> > extra_fitnesses(extra_funs.size());
    emp::vector<double> max_extra_fit(extra_funs.size(), 0.0);
    emp::vector<size_t> max_count(extra_funs.size(), 0);
    for (size_t i = 0; i < extra_funs.size(); ++i) {
      extra_fitnesses[i].resize(cohort.size());
    }

    // Collect all fitness info (from each individual in cohort)
    for (size_t cID = 0; cID < cohort.size(); ++cID) {
      const size_t worldID = cohort[cID];
      base_fitness[cID] = world.CalcFitnessID(worldID);
      for (size_t exID = 0; exID < extra_funs.size(); ++exID) {
        double cur_fit = extra_funs[exID](world.GetOrg(worldID));
        extra_fitnesses[exID][cID] = cur_fit;
        if (cur_fit > max_extra_fit[exID]) {
          max_extra_fit[exID] = cur_fit;
          max_count[exID] = 1;
        } else if (cur_fit == max_extra_fit[exID]) {
          max_count[exID]++;
        }
      }
    }

    // Readjust base fitnesses to reflect extra resources.
    for (size_t exID = 0; exID < extra_funs.size(); ++exID) {
      if (max_count[exID] == 0) continue; // No one gets this reward.
      // The current bonus is divided up among the organisms that earned it.
      const double cur_bonus = pool_sizes[exID] / max_count[exID];
      for (size_t cID = 0; cID < cohort.size(); ++cID) {
        // If this organism is the best at the current resource, give it the bonus.
        if (extra_fitnesses[exID][cID] == max_extra_fit[exID]) {
          base_fitness[cID] += cur_bonus;
        }
      }
    }

    // Do tournament selection w/adjusted fitnesses
    emp::vector<size_t> entries;
    for (size_t T = 0; T < tourny_count; ++T) {
      entries.clear();
      for (size_t i = 0; i < t_size; ++i) {
        entries.push_back(world.GetRandom().GetUInt(0, cohort.size()));
      }
      double best_fit = base_fitness[entries[0]];
      size_t best_id = entries[0];
      // Search for a higher fit org in tournament.
      for (size_t i = 1; i < t_size; ++i) {
        const double cur_fit = base_fitness[entries[i]];
        if (cur_fit > best_fit) {
          best_fit = cur_fit;
          best_id = entries[i];
        }
      }

      // Place highest fitness into the next generation.
      const size_t worldID = cohort[best_id];
      world.DoBirth(world.GetGenomeAt(worldID), worldID, 1);
    }
  }

}

#endif