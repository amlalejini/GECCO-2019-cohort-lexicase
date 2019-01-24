#ifndef SORTING_TEST_ORG_H
#define SORTING_TEST_ORG_H

#include "SortingTest.h"

class SortingTestOrg {
public:

  struct Genome {
    emp::vector<SortingTest> test_set;
    size_t test_size;

    Genome(size_t test_size, size_t num_tests=1) 
      : test_set(num_tests, test_size), test_size(test_size) { ; }

    Genome(emp::Random & rnd, size_t test_size, size_t num_tests=1)
      : test_set(num_tests, test_size), test_size(test_size) 
    {
        for (size_t i = 0; i < test_set.size(); ++i) test_set[i].RandomizeTest(rnd);
    }

    Genome(Genome &&) = default;
    Genome(const Genome &) = default;
    
    Genome & operator=(const Genome &) = default;
    Genome & operator=(Genome &&) = default;

    bool operator==(const Genome & in) const { return in.test_set == test_set; }
    bool operator!=(const Genome & in) const { return !(in == *this); }
    bool operator<(const Genome & in) const { return test_set < in.test_set; }

    void Randomize(emp::Random & rnd) {
      for (size_t i = 0; i < test_set.size(); ++i) {
        test_set[i].RandomizeTest(rnd);
      }
    }

    bool Validate(int min_val, int max_val) {
      for (size_t tID = 0; tID < test_set.size(); ++tID) {
        if (!test_set[tID].Validate(test_size, min_val, max_val)) return false;
      }
      return true;
    }

  };

  struct Phenotype {
    emp::vector<size_t> test_results; ///< Correspond to per-organism passes, not per-test passes!
    size_t num_passes;
    size_t num_fails;

    void Reset(size_t s=0) { 
      test_results.clear(); 
      test_results.resize(s, 0); 
      num_passes = 0;
      num_fails = 0;
    }
  };

  using phenotype_t = Phenotype;
  using genome_t = Genome;

protected:

  genome_t genome;
  phenotype_t phenotype;

public:

  SortingTestOrg(size_t test_size, size_t num_tests=1) 
    : genome(test_size, num_tests), phenotype()
  { ;  } 

  SortingTestOrg(emp::Random & rnd, size_t test_size, size_t num_tests=1)
    : genome(rnd, test_size, num_tests), phenotype()
  { ; }
  
  SortingTestOrg(const genome_t & _g) : genome(_g), phenotype() { ; }

  size_t GetNumTests() const { return genome.test_set.size(); }
  size_t GetTestSize() const { return genome.test_size; }
  emp::vector<SortingTest> & GetTestSet() { return genome.test_set; }

  genome_t & GetGenome() { return genome; }
  const genome_t & GetGenome() const { return genome; }
  
  phenotype_t & GetPhenotype() { return phenotype; }

  void Print(std::ostream & out=std::cout) const;
  void PrintMin(std::ostream & out=std::cout) const;

};

void SortingTestOrg::Print(std::ostream & out) const {
  out << "TestOrg(seqsize=" << GetTestSize() << "," << "numtests=" << GetNumTests() << "):\n";
  for (size_t i = 0; i < GetNumTests(); ++i) {
    out << "  Test[" << i << "]: "; 
    genome.test_set[i].Print();
    out << "\n";
  }
}

void SortingTestOrg::PrintMin(std::ostream & out) const {
  out << "[";
  for (size_t i = 0; i < GetNumTests(); ++i) {
    if (i) out << ",";
    genome.test_set[i].Print(out);
  }
  out << "]";
}

// todo - flat print

#endif