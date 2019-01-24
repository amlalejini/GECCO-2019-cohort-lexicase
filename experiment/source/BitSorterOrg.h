#ifndef BIT_SORTER_ORG_H
#define BIT_SORTER_ORG_H

#include "base/array.h"
#include "base/vector.h"
#include "tools/Random.h"

#include "hardware/BitSorter.h"

class BitSorterOrg {
public:

  struct Phenotype {
    emp::vector<double> test_scores;
    double total_score;
    
    emp::vector<size_t> test_passes;
    
    size_t num_passes;
    size_t num_fails;

    void Reset(size_t s=0) {
      test_scores.clear();
      test_scores.resize(s, 0);
      total_score = 0;

      test_passes.clear();
      test_passes.resize(s, 0);

      num_passes = 0;
      num_fails = 0;
    }

    void RecordPassFail(size_t testID, bool pass) {
      test_passes[testID] = (size_t)pass;
      if (pass) ++num_passes;
      else ++num_fails;
    }

    void RecordScore(size_t testID, double score) {
      test_scores[testID] = score;
      total_score += score;
    }
  };

  using phenotype_t = Phenotype;
  using genome_t = emp::BitSorter;

protected:

  genome_t genome;
  phenotype_t phenotype;

public:

  BitSorterOrg() : genome(), phenotype() { ; }
  BitSorterOrg(const genome_t & _g) : genome(_g), phenotype() { ; }

  size_t GetSize() const { return genome.GetSize(); }

  genome_t & GetGenome() { return genome; }
  const genome_t & GetGenome() const { return genome; }

  phenotype_t & GetPhenotype() { return phenotype; }
  const phenotype_t & GetPhenotype() const { return phenotype; }

  void SetGenome(const genome_t & in) { genome = in; }

  void RandomizeGenome(emp::Random & rnd, size_t input_size, size_t network_size);

  void Print(std::ostream & out=std::cout) const;
};

void BitSorterOrg::RandomizeGenome(emp::Random & rnd, size_t input_size, size_t network_size) {
  emp::BitSorter rand_sorter;
  for (size_t i = 0; i < network_size; ++i) {
    rand_sorter.AddCompare(rnd.GetUInt(0, input_size), rnd.GetUInt(0, input_size));
  }
  genome = rand_sorter;
}

void BitSorterOrg::Print(std::ostream & out) const {
  out << genome.AsString();
}

#endif