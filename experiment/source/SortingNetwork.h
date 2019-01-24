#ifndef SORTING_NETWORK_H
#define SORTING_NETWORK_H

#include "base/array.h"
#include "base/vector.h"
#include "tools/Random.h"

class SortingNetwork {
public:
  using op_t = emp::array<size_t, 2>;
  using network_t = emp::vector<op_t>;

protected:
  network_t network;

public:
  SortingNetwork() 
    : network() { ; }
  
  SortingNetwork(emp::Random & rnd, size_t input_size, size_t min_network_size, size_t max_network_size)
    : network() 
  { RandomizeNetwork(rnd, input_size, min_network_size, max_network_size); }
  
  SortingNetwork(emp::Random & rnd, size_t input_size, size_t network_size) 
    : network()
  { RandomizeNetwork(rnd, input_size, network_size); }

  SortingNetwork(size_t network_size) 
    : network(network_size) 
  { 
    for (size_t i = 0; i < network.size(); ++i) {
      network[i][0] = 0; network[i][1] = 0;
    }
  }

  SortingNetwork(SortingNetwork &&) = default;
  SortingNetwork(const SortingNetwork &) = default;

  SortingNetwork & operator=(const SortingNetwork & in) { network = in.network; return *this; }

  bool operator==(const SortingNetwork & in) const { return in.network == network; }
  bool operator!=(const SortingNetwork & in) const { return !(in == *this); }
  bool operator<(const SortingNetwork & in) const { return network < in.network; }

  op_t & operator[](size_t id) {
    emp_assert(id < network.size());
    return network[id];
  }

  const op_t & operator[](size_t id) const {
    emp_assert(id < network.size());
    return network[id];
  }

  size_t GetSize() const { return network.size(); }
  network_t & GetNetwork() { return network; }

  void RandomizeNetwork(emp::Random & rnd, size_t input_size, size_t network_size);
  void RandomizeNetwork(emp::Random & rnd, size_t input_size, size_t min_network_size, size_t max_network_size);

  bool Validate(size_t input_size, size_t min_network_size=0, size_t max_network_size=(size_t)-1) const;

  void Print(std::ostream & out=std::cout, std::string op_sep="=>") const;
  void PrintVert(std::ostream & out=std::cout) const;

};

void SortingNetwork::RandomizeNetwork(emp::Random & rnd, size_t input_size, 
                                        size_t network_size) {
  network.resize(network_size);
  for (size_t i = 0; i < network.size(); ++i) {
    network[i][0] = rnd.GetUInt(0, input_size);
    network[i][1] = rnd.GetUInt(0, input_size);
  }
}

void SortingNetwork::RandomizeNetwork(emp::Random & rnd, size_t input_size,
                                     size_t min_network_size, size_t max_network_size) {
  const size_t s = rnd.GetUInt(min_network_size, max_network_size+1);
  RandomizeNetwork(rnd, input_size, s);
}

bool SortingNetwork::Validate(size_t input_size, size_t min_network_size, 
                                     size_t max_network_size) const {
  if (GetSize() > max_network_size || GetSize() < min_network_size) return false;
  for (size_t i = 0; i < network.size(); ++i) {
    if (network[i][0] >= input_size || network[i][1] >= input_size) return false;
  }
  return true;
}

void SortingNetwork::Print(std::ostream & out, std::string op_sep) const {
  out << "[";
  for (size_t i=0; i < network.size(); ++i) {
    if (i) out << ",";
    out << "(" << network[i][0] << op_sep << network[i][1] << ")";
  }
  out << "]";
}

void SortingNetwork::PrintVert(std::ostream & out) const {
  for (size_t i=0; i < network.size(); ++i) {
    if (i) out << "\n";
    out << network[i][0] << " <=> " << network[i][1];
  }
} 

#endif