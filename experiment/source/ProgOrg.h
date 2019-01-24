#ifndef PROG_ORG_H
#define PROG_ORG_H

#include "base/array.h"
#include "base/vector.h"
#include "tools/Random.h"

#include "TagLinearGP.h"

template<size_t TAG_WIDTH>
class ProgOrg {
public:
  struct Phenotype;

  using phenotype_t = Phenotype;
  using genome_t = typename TagLGP::TagLinearGP_TW<TAG_WIDTH>::Program;

  struct Phenotype {
    // emp::vector<double> test_results;
    emp::vector<double> test_scores;
    double total_score;
    
    emp::vector<bool> test_passes;
    size_t num_passes;
    size_t num_fails;

    size_t total_submissions;

    void Reset(size_t s=0) {
      // test_results.clear();
      // test_results.resize(s, 0);

      test_scores.clear();
      test_scores.resize(s, 0);
      total_score = 0;
      
      test_passes.clear();
      test_passes.resize(s, false);
      num_passes = 0;
      num_fails = 0;

      total_submissions = 0;
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

    void RecordSubmission(bool sub) {
      total_submissions += (size_t)sub;
    }

  };

protected:

  phenotype_t phenotype;
  genome_t genome;

public:

  ProgOrg(const genome_t & _g) 
    : phenotype(), genome(_g) { ; }

  ProgOrg(const ProgOrg &) = default;
  ProgOrg(ProgOrg &&) = default;

  genome_t & GetGenome() { return genome; }
  const genome_t & GetGenome() const { return genome; }

  phenotype_t & GetPhenotype() { return phenotype; }
  const phenotype_t & GetPhenotype() const { return phenotype; }

  void SetGenome(const genome_t & in) { genome = in; }

  // Todo: randomize genome

  // Todo: print linear (print such that program can be put on single line)
};

#endif