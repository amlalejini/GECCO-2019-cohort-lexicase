#ifndef BIT_SORTER_TEST_ORG_H
#define BIT_SORTER_TEST_ORG_H

#include "base/array.h"
#include "base/vector.h"
#include "tools/Random.h"

class BitTestOrg {

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
  
  using test_t = uint32_t;
  using genome_t = test_t;
  using phenotype_t = Phenotype;

protected:

  genome_t genome;
  phenotype_t phenotype;

public:

  BitTestOrg() : genome(0), phenotype()  { ; }
  BitTestOrg(const genome_t & _g) : genome(_g), phenotype() { ; }

  genome_t & GetGenome() { return genome; }
  const genome_t & GetGenome() const { return genome; }  

  void SetGenome(uint32_t in) { genome = in; } 

  phenotype_t & GetPhenotype() { return phenotype; }
  const phenotype_t & GetPhenotype() const { return phenotype; }

  void Print(size_t num_bits, std::ostream & out=std::cout) const;
};

void BitTestOrg::Print(size_t num_bits, std::ostream & out) const {
  for (size_t id = num_bits; id > 0; id--) {
    if (genome & 1 << (id-1)) out << "1";
    else out << "0";
  }
}

#endif