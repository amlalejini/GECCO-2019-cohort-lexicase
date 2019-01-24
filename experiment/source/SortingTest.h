#ifndef SORTING_TEST_H
#define SORTING_TEST_H

#include <algorithm>
#include <iterator>

#include "base/array.h"
#include "base/vector.h"
#include "tools/Random.h"

#include "SortingNetwork.h"

class SortingTest {
public:
  using test_t = emp::vector<int>;

protected:
  test_t test;

public:
  SortingTest()
    : test() { ; }
  
  SortingTest(emp::Random & rnd, size_t test_size, int min_val=0, int max_val=1)
    : test(test_size) { RandomizeTest(rnd, test_size, min_val, max_val); }

  SortingTest(size_t test_size) : test(test_size, 0) { ; }

  SortingTest(SortingTest &&) = default;
  SortingTest(const SortingTest &) = default;

  SortingTest & operator=(const SortingTest &) = default;
  SortingTest & operator=(SortingTest &&) = default;

  bool operator==(const SortingTest & in) const { return in.test == test; }
  bool operator!=(const SortingTest & in) const { return !(in == *this); }
  bool operator<(const SortingTest & in) const { return test < in.test; }

  int & operator[](size_t id) {
    emp_assert(id < test.size());
    return test[id];
  }

  const int & operator[](size_t id) const {
    emp_assert(id < test.size());
    return test[id];
  }
  
  size_t GetSize() const { return test.size(); }
  test_t & GetTest() { return test; }

  void RandomizeTest(emp::Random & rnd, size_t test_size, int min_val=0, int max_val=1);
  void RandomizeTest(emp::Random & rnd, int min_val=0, int max_val=1);

  bool Evaluate(const SortingNetwork & network) const;

  bool Validate(size_t test_size, int min_val, int max_val);

  void Print(std::ostream & out=std::cout) const;

};

void SortingTest::RandomizeTest(emp::Random & rnd, size_t test_size, int min_val, int max_val) {
  test.resize(test_size);
  RandomizeTest(rnd, min_val, max_val);
}

void SortingTest::RandomizeTest(emp::Random & rnd, int min_val, int max_val) {
  for (size_t i = 0; i < test.size(); ++i) test[i] = rnd.GetInt(min_val, max_val+1);
}

bool SortingTest::Evaluate(const SortingNetwork & network) const {
  emp_assert(network.Validate(test.size()));
  // Make copy of test to be sorted.
  test_t eval_test(test);
  // Evaluate the sorting network.
  for (size_t ni = 0; ni < network.GetSize(); ++ni) {
    // if [j] < [i] => swap([i],[j])
    const size_t i = network[ni][0];
    const size_t j = network[ni][1];
    
    size_t head;
    size_t tail;
    if (i < j) { head = i; tail = j; } 
    else { head = j; tail = i; }

    if (eval_test[tail] < eval_test[head]) {
      std::swap(eval_test[tail], eval_test[head]);
    }
  }
  // Is eval test sorted?
  return std::is_sorted(std::begin(eval_test), std::end(eval_test));
}

bool SortingTest::Validate(size_t test_size, int min_val, int max_val) {
  if (test.size() != test_size) return false; 
  for (size_t i = 0; i < test.size(); ++i) {
    if (test[i] < min_val || test[i] > max_val) return false;
  }
  return true;
}

void SortingTest::Print(std::ostream & out) const {
  out << "[";
  for (size_t i = 0; i < test.size(); ++i) {
    if (i) out << ",";
    out << test[i];
  }
  out << "]";
}




#endif