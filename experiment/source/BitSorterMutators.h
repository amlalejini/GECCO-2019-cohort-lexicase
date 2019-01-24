#ifndef ALEX_BITSORTER_MUTATORS_H
#define ALEX_BITSORTER_MUTATORS_H

#include "base/vector.h"
#include "tools/Random.h"
#include "tools/random_utils.h"

#include "SortingNetworkOrg.h"
#include "SortingTestOrg.h"
#include "SortingTest.h"

#include "TagLinearGP.h"
#include "TagLinearGP_Utilities.h"

#include "BitSorterOrg.h"
#include "BitTestOrg.h"

/// Bit sorter organism mutator
struct BitSorterMutator {
  using genome_t = BitSorterOrg::genome_t;

  size_t MAX_NETWORK_SIZE;  ///< Maximum size network can grow
  size_t MIN_NETWORK_SIZE;  ///< Minimum size network can shrink
  size_t SORT_SEQ_SIZE;     ///< Sort input size (defines range for i,j values)

  double PER_INDEX_SUB;
  double PER_PAIR_DUP;
  double PER_PAIR_INS;
  double PER_PAIR_DEL;
  double PER_PAIR_SWAP;

  BitSorterMutator() 
    : MAX_NETWORK_SIZE(64),
      MIN_NETWORK_SIZE(1),
      SORT_SEQ_SIZE(4),
      PER_INDEX_SUB(0.001),
      PER_PAIR_DUP(0.001),
      PER_PAIR_INS(0.001),
      PER_PAIR_DEL(0.001),
      PER_PAIR_SWAP(0.001)
  { ; }

  size_t Mutate(emp::Random & rnd, genome_t & genome) {
    size_t muts = 0;
    emp::vector<std::pair<size_t, size_t>> sorting_network; // Extract bit sorter into more convenient form for mutating?
    size_t expected_size = genome.GetSize();
    for (size_t i = 0; i < genome.GetSize(); ++i) {
      // Deletions!
      if (rnd.P(PER_PAIR_DEL) && (expected_size > MIN_NETWORK_SIZE)) {
        ++muts;
        --expected_size;
        continue;
      }
      // Copy over.
      const size_t whead = sorting_network.size();
      sorting_network.emplace_back(genome.GetComparator(i));
      // Do we insert?
      if (rnd.P(PER_PAIR_INS) && (expected_size < MAX_NETWORK_SIZE)) {
        ++muts;
        ++expected_size;
        // Insert randomly
        sorting_network.emplace_back(std::pair<size_t,size_t>{rnd.GetUInt(0, SORT_SEQ_SIZE), rnd.GetUInt(0, SORT_SEQ_SIZE)});
      }
      // Do we duplicate?
      if (rnd.P(PER_PAIR_DUP) && (expected_size < MAX_NETWORK_SIZE)) {
        ++muts;
        ++expected_size;
        // Duplicate!
        sorting_network.emplace_back(genome.GetComparator(i));
      }
      // Per-index substitutions?
      if (rnd.P(PER_INDEX_SUB)) {
        sorting_network[whead].first = rnd.GetUInt(0, SORT_SEQ_SIZE);
        ++muts;
      }
      if (rnd.P(PER_INDEX_SUB)) {
        sorting_network[whead].second = rnd.GetUInt(0, SORT_SEQ_SIZE);
        ++muts;
      }
    }
    // How about swaps?
    if (PER_PAIR_SWAP > 0) {
      for (size_t i = 0; i < sorting_network.size(); ++i) {
        if (rnd.P(PER_PAIR_SWAP)) {
          // Select two random positions
          const size_t pos = rnd.GetUInt(sorting_network.size());
          if (pos == i) continue;
          std::swap(sorting_network[i], sorting_network[pos]);
          ++muts;
        }
      }
    }
    // Update genome with mutated network!
    genome.Clear();
    for (std::pair<size_t,size_t> & comp : sorting_network) genome.AddCompare(comp.first, comp.second);
    return muts;
  }

  genome_t GenRandomBitSorter(emp::Random & rnd) {
    emp::BitSorter rando_sorter;
    // How big?
    const size_t sorter_size = rnd.GetUInt(MAX_NETWORK_SIZE+1);
    // Build!
    for (size_t i = 0; i < sorter_size; ++i) {
      rando_sorter.AddCompare(rnd.GetUInt(0, SORT_SEQ_SIZE),rnd.GetUInt(0, SORT_SEQ_SIZE));
    }
    return rando_sorter;
  }
};

// Note - look at the BitSorter ToString method to build bitvectors to mutate
// - then go back to uint32
// Generate random - limit = 1 << num_bits
struct BitTestMutator {
  using genome_t = uint32_t;

  size_t NUM_BITS;
  double PER_BIT_FLIP;
  double PER_SEQ_INVERSION;
  double PER_SEQ_RANDOMIZE;

  BitTestMutator()
    : NUM_BITS(16),
      PER_BIT_FLIP(0.001),
      PER_SEQ_INVERSION(0.01),
      PER_SEQ_RANDOMIZE(0.01)
  { ; }

  size_t Mutate(emp::Random & rnd, genome_t & genome) {
    size_t muts = 0;
    
    // Randomize?
    if (rnd.P(PER_SEQ_RANDOMIZE)) {
      genome = rnd.GetUInt(0, 1 << NUM_BITS);
      return 1;
    }
    // Past this point, it's easier to work with a bit vector.
    emp::BitVector bit_vec(emp::BitSorter::ToBitVector(genome, NUM_BITS));
    
    // Per-site bit flips
    for (size_t i = 0; i < bit_vec.GetSize(); ++i) {
      if (rnd.P(PER_BIT_FLIP)) {
        ++muts;
        bit_vec.Set(i, !bit_vec.Get(i));
      }
    }

    // TODO - test inversions!
    // Inversions?
    if (rnd.P(PER_SEQ_INVERSION)) {
      ++muts;
      int p0 = (int)rnd.GetUInt(0, bit_vec.GetSize());
      int p1 = (int)rnd.GetUInt(0, bit_vec.GetSize());
      if (p1 < p0) std::swap(p0, p1);
      emp_assert(p0<=p1);
      while (p0 < p1) {
        // Swap p0 and p1 values
        bool p0_val = bit_vec.Get(p0);
        bool p1_val = bit_vec.Get(p1);
        bit_vec.Set(p0, p1_val);
        bit_vec.Set(p1, p0_val);
        // std::swap(bit_vec[p0], bit_vec[p1]);
        ++p0; --p1;
      }
    }

    // Update genome being mutated
    genome = bit_vec.GetUInt(0);
    return muts;

  }

};


#endif