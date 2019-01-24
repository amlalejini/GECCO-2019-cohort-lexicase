#ifndef ALEX_MUTATORS_H
#define ALEX_MUTATORS_H

#include "base/vector.h"
#include "tools/Random.h"
#include "tools/random_utils.h"

#include "SortingNetworkOrg.h"
#include "SortingTestOrg.h"
#include "SortingTest.h"

#include "TagLinearGP.h"
#include "TagLinearGP_Utilities.h"

// For now, all mutation types are integrated into Mutate function.
// To make cleaner, could break each out into their own function.
struct SortingNetworkMutator {
  using genome_t = SortingNetworkOrg::genome_t;

  size_t MAX_NETWORK_SIZE;  ///< Maximum size network can grow
  size_t MIN_NETWORK_SIZE;  ///< Minimum size network can shrink
  size_t SORT_SEQ_SIZE;     ///< Sort input size (defines range for i,j values)

  double PER_INDEX_SUB;
  double PER_PAIR_DUP;
  double PER_PAIR_INS;
  double PER_PAIR_DEL;
  double PER_PAIR_SWAP;

  SortingNetworkMutator() 
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
    using gene_t = emp::array<size_t, 2>;

    genome_t new_genome;
    size_t mut_cnt = 0;

    size_t expected_size = genome.GetSize();

    // For gene (compare-exchange operation) in genome:
    // - copy gene to new genome, applying mutations
    for (size_t geneID = 0; geneID < genome.GetSize(); ++geneID) {

      // Do we delete?
      if (rnd.P(PER_PAIR_DEL) && (expected_size > MIN_NETWORK_SIZE)) { 
        --expected_size;
        ++mut_cnt; 
        continue; 
      }
      const size_t rhead = new_genome.GetSize(); // Where in the new genome are we copying to?
      new_genome.GetNetwork().emplace_back(gene_t{genome[geneID][0], genome[geneID][1]});
      // gene_t & gene_copy = new_genome.GetNetwork()[rhead];

      // Do we insert?
      if (rnd.P(PER_PAIR_INS) && (expected_size < MAX_NETWORK_SIZE)) {
        new_genome.GetNetwork().emplace_back(gene_t{rnd.GetUInt(0, SORT_SEQ_SIZE), rnd.GetUInt(0, SORT_SEQ_SIZE)});
        ++expected_size;
        ++mut_cnt;
      }

      // Do we duplicate?
      if (rnd.P(PER_PAIR_DUP) && (expected_size < MAX_NETWORK_SIZE)) {
        new_genome.GetNetwork().emplace_back(gene_t{genome[geneID][0], genome[geneID][1]});
        ++expected_size;
        ++mut_cnt;
      }
      
      // Do index changes?
      if (rnd.P(PER_INDEX_SUB)) {
        new_genome.GetNetwork()[rhead][0] = rnd.GetUInt(0, SORT_SEQ_SIZE);
        ++mut_cnt;
      }

      if (rnd.P(PER_INDEX_SUB)) {
        new_genome.GetNetwork()[rhead][1] = rnd.GetUInt(0, SORT_SEQ_SIZE);
        ++mut_cnt;
      }
    }

    // Do swaps?
    if (PER_PAIR_SWAP > 0.0) {
      for (size_t geneID = 0; geneID < new_genome.GetSize(); ++geneID) {
        if (rnd.P(PER_PAIR_SWAP)) {
          // Select two random positions
          const size_t pos = rnd.GetUInt(new_genome.GetSize());
          if (pos == geneID) continue;
          std::swap(new_genome.GetNetwork()[geneID], new_genome.GetNetwork()[pos]);
          ++mut_cnt;
        }
      }
    }

    genome = new_genome;
    return mut_cnt;
  }

  /// Given genome A and B, crossover.
  /// If (Valid(AB)) A = AB; If(Valid(BA)) B = BA;
  void Crossover1Pt(emp::Random & rnd, genome_t & genomeA, genome_t & genomeB) {
    const size_t min_size = emp::Min(genomeA.GetSize(), genomeB.GetSize());
    const size_t max_size = emp::Max(genomeA.GetSize(), genomeB.GetSize());
    const size_t xpoint = rnd.GetUInt(min_size);
    genome_t genomeAB;
    genome_t genomeBA;

    for (size_t r = 0; r < xpoint; ++r) {
      genomeAB.GetNetwork().emplace_back(genomeA[r]);
      genomeBA.GetNetwork().emplace_back(genomeB[r]);
    }

    for (size_t r = xpoint; r < max_size; ++r) {
      if (r < genomeA.GetSize()) { genomeBA.GetNetwork().emplace_back(genomeA[r]); } 
      if (r < genomeB.GetSize()) { genomeAB.GetNetwork().emplace_back(genomeB[r]); }
    }

    genomeA = genomeAB;
    genomeB = genomeBA;
  }

  /// Given genome A and B, crossover.
  /// If (Valid(ABA)) A = ABA; If(Valid(BAB)) B = BAB;
  void Crossover2Pt(emp::Random & rnd, genome_t & genomeA, genome_t & genomeB) {
    // std::cout << "Genome A: "; genomeA.Print(); std::cout << std::endl;
    // std::cout << "Genome B: "; genomeB.Print(); std::cout << std::endl;

    double pct1 = rnd.GetDouble();
    double pct2 = rnd.GetDouble();
    if (pct2 < pct1) std::swap(pct1, pct2); // pct1 < pct2

    size_t pos1A = (size_t)genomeA.GetSize() * pct1;
    size_t pos1B = (size_t)genomeB.GetSize() * pct1;

    size_t pos2A = (size_t)genomeA.GetSize() * pct2;
    size_t pos2B = (size_t)genomeB.GetSize() * pct2;

    size_t ABA_size = pos1A + (pos2B - pos1B) + (genomeA.GetSize() - pos2A);
    size_t BAB_size = pos1B + (pos2A - pos1A) + (genomeB.GetSize() - pos2B);

    // std::cout << "pct1=" << pct1 << std::endl;
    // std::cout << "pct2=" << pct2 << std::endl;

    // std::cout << "pA = [" << pos1A << "," << pos2A << "]" << std::endl;
    // std::cout << "pB = [" << pos1B << "," << pos2B << "]" << std::endl;

    // std::cout << "ABA size = " << ABA_size << std::endl;
    // std::cout << "BAB size = " << BAB_size << std::endl;

    // ABA: |---A---|---B---|---A---|
    // BAB: |---B---|---A---|---B---|
    const bool build_aba = ABA_size <= MAX_NETWORK_SIZE && ABA_size >= MIN_NETWORK_SIZE;
    const bool build_bab = BAB_size <= MAX_NETWORK_SIZE && BAB_size >= MIN_NETWORK_SIZE;
    genome_t genomeABA;
    genome_t genomeBAB;

    if (build_aba) {
      // std::cout << "Buildinb ABA" << std::endl;    
      // Build ABA
      for (size_t r = 0; r < pos1A; ++r) { 
        // std::cout << "  Copying A[" << r << "]: " << genomeA[r][0] << "=>" << genomeA[r][1] << std::endl; 
        genomeABA.GetNetwork().emplace_back(genomeA[r]); 
      }
      for (size_t r = pos1B; r < pos2B; ++r) { 
        // std::cout << "  Copying B[" << r << "]: " << genomeB[r][0] << "=>" << genomeB[r][1] << std::endl;
        genomeABA.GetNetwork().emplace_back(genomeB[r]); 
      }
      for (size_t r = pos2A; r < genomeA.GetSize(); ++r) {
        // std::cout << "  Copying A[" << r << "]: " << genomeA[r][0] << "=>" << genomeA[r][1] << std::endl; 
        genomeABA.GetNetwork().emplace_back(genomeA[r]); 
      }
      emp_assert(genomeABA.GetSize() <= MAX_NETWORK_SIZE && genomeABA.GetSize() >= MIN_NETWORK_SIZE, ABA_size, genomeABA.GetSize(), pos1A, pos1B, pos2A, pos2B);
    }

    if (build_bab) {
      // std::cout << "Building BAB" << std::endl;
      // Build BAB
      for (size_t r = 0; r < pos1B; ++r) { 
        // std::cout << "  Copying B[" << r << "]: " << genomeB[r][0] << "=>" << genomeB[r][1] << std::endl;
        genomeBAB.GetNetwork().emplace_back(genomeB[r]); 
      }
      for (size_t r = pos1A; r < pos2A; ++r) { 
        // std::cout << "  Copying A[" << r << "]: " << genomeA[r][0] << "=>" << genomeA[r][1] << std::endl; 
        genomeBAB.GetNetwork().emplace_back(genomeA[r]); 
      }
      for (size_t r = pos2B; r < genomeB.GetSize(); ++r) {
        // std::cout << "  Copying B[" << r << "]: " << genomeB[r][0] << "=>" << genomeB[r][1] << std::endl;
        genomeBAB.GetNetwork().emplace_back(genomeB[r]); 
      }
      emp_assert(genomeBAB.GetSize() <= MAX_NETWORK_SIZE && genomeBAB.GetSize() >= MIN_NETWORK_SIZE, BAB_size, genomeBAB.GetSize(), pos1A, pos1B, pos2A, pos2B);
    }

    if (build_aba) genomeA = genomeABA;
    if (build_bab) genomeB = genomeBAB;

  }

};

/// Sorting network test mutator
struct SortingTestMutator {
  using genome_t = SortingTestOrg::genome_t;

  bool bit_mode;

  int MAX_VALUE;
  int MIN_VALUE;
  
  double PER_SITE_SUB;
  double PER_SEQ_INVERSION;
  double PER_SEQ_RANDOMIZE;

  SortingTestMutator() 
    : bit_mode(true),
      MAX_VALUE(1),
      MIN_VALUE(0),
      PER_SITE_SUB(0.001),
      PER_SEQ_INVERSION(0.01),
      PER_SEQ_RANDOMIZE(0.0)
  { ; }

  size_t Mutate(emp::Random & rnd, genome_t & genome) {
    size_t mut_cnt = 0;
    
    // For each sorting test in test set.
    for (size_t testID = 0; testID < genome.test_set.size(); ++testID) {
      SortingTest & test = genome.test_set[testID];
      
      if (rnd.P(PER_SEQ_RANDOMIZE)) {
        test.RandomizeTest(rnd);
        ++mut_cnt;
        continue;
      }
      
      // For each site in test
      for (size_t i = 0; i < test.GetSize(); ++i) {
        if (rnd.P(PER_SITE_SUB)) {
          if (bit_mode) test[i] = (int)!((bool)test[i]);
          else test[i] = rnd.GetInt(MIN_VALUE, MAX_VALUE+1);
          ++mut_cnt;
        }
      }
      // Inversions?
      if (rnd.P(PER_SEQ_INVERSION)) {
        int p0 = (int)rnd.GetUInt(0, test.GetSize());
        int p1 = (int)rnd.GetUInt(0, test.GetSize());
        if (p1 < p0) std::swap(p0, p1);
        emp_assert(p0 <= p1);
        while (p0 < p1) {
          std::swap(test[p0], test[p1]);
          ++p0; --p1;
        }
      }

    }
    return mut_cnt;
  }

};

/// Mutator for tag LGP programs
template<size_t TAG_WIDTH>
struct TagLGPMutator {
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>; 
  using program_t = typename hardware_t::Program;
  using tag_t = typename hardware_t::tag_t;
  using inst_t = typename hardware_t::inst_t;
  using inst_lib_t = typename hardware_t::inst_lib_t;
  using inst_prop_t = typename hardware_t::inst_prop_t;
  
  size_t MAX_PROGRAM_LEN;
  size_t MIN_PROGRAM_LEN;

  double PER_BIT_FLIP;
  double PER_INST_SUB;
  double PER_INST_INS;
  double PER_INST_DEL;
  double PER_PROG_SLIP;
  double PER_MOD_DUP;
  double PER_MOD_DEL;

  enum ModuleMutType { DUP=0, DEL, NONE };

  struct ModuleInfo {
    int begin;  // Includes ModDef instruction
    int end;
    ModuleMutType mut;
    emp::vector<size_t> positions;
    ModuleInfo(int _b, int _e) : begin(_b), end(_e), mut(ModuleMutType::NONE), positions() { ; }
    size_t GetSize() const { return positions.size(); }
  };

  size_t Mutate(emp::Random & rnd, program_t & program) {
    const inst_lib_t & ilib = program.GetInstLib();
    
    int expected_size = program.GetSize();
    size_t mut_cnt = 0;
    int rhead = 0;

    // assert that expected size must be in correct range
    /////////////////////////////////////////// 
    // Whole module duplications/deletions
    if (PER_MOD_DEL != 0.0 || PER_MOD_DUP != 0.0) {
      
      // Scan for modules
      size_t mod_mut_cnt = 0;
      emp::vector<ModuleInfo> modules;  // Modules we're going to dup or delete.
      emp::vector<size_t> module_membership(program.GetSize(), 0);
      emp::vector<bool> dels(program.GetSize(), false);
      emp::vector<size_t> danglers;
      
      for (size_t i = 0; i < program.GetSize(); ++i) {
        // inst_t & inst = program[i];
        if (ilib.HasProperty(program[i].id, inst_prop_t::MODULE)) {
          if (modules.size()) { modules.back().end = i; }
          const size_t mod_id = modules.size();
          modules.emplace_back(i, -1);
          modules.back().positions.emplace_back(i); // Add start position to module
          module_membership[i] = mod_id;
        } else {
          // Not a new module definition
          if (modules.size()) { modules.back().positions.emplace_back(i); module_membership[i] = modules.size()-1; }
          else { danglers.emplace_back(i); }
        }
      }

      // Take care of danglers.
      if (modules.size()) {
        if (modules[0].begin == 0) { modules.back().end = program.GetSize(); }  // Case where last module does not wrap.
        else { modules.back().end = modules[0].begin-1; }                       // Last module wraps around.
        for (size_t i = 0; i < danglers.size(); ++i) {
          modules.back().positions.emplace_back(danglers[i]); 
          module_membership[i] = modules.size()-1;
        }
      }

      // Did we do the above thing correctly?
      // std::cout << "Modules found: " << modules.size() << std::endl;
      // for (size_t i = 0; i < modules.size(); ++i) {
      //   std::cout << "  Module " << i << "["<<modules[i].begin << "," << modules[i].end << ")" << std::endl;
      //   for (size_t p = 0; p < modules[i].positions.size(); ++p) {
      //     if (p == 0) std::cout << "    => [";
      //     else std::cout << ",";
      //     std::cout << modules[i].positions[p];
      //   }
      //   std::cout << "]" << std::endl;
      // }

      // Are we mutating any of the modules?
      for (size_t i = 0; i < modules.size(); ++i) {
        bool mod_dup = rnd.P(PER_MOD_DUP);
        bool mod_del = rnd.P(PER_MOD_DEL);
        if (mod_dup && mod_del) { mod_dup = false; mod_dup = false; } // If we would both dup and delete module, do nothing instead.
        if (mod_dup && ((expected_size + modules[i].GetSize()) <= MAX_PROGRAM_LEN) ) { 
          mod_mut_cnt++; 
          mut_cnt++; 
          modules[i].mut = ModuleMutType::DUP; 
          expected_size += modules[i].GetSize(); 
        }
        if (mod_del && (((int)expected_size - (int)modules[i].GetSize()) >= (int)MIN_PROGRAM_LEN) ) { 
          mod_mut_cnt++; 
          mut_cnt++; 
          modules[i].mut = ModuleMutType::DEL; 
          expected_size -= modules[i].GetSize(); 
          for (size_t p = 0; p < modules[i].GetSize(); ++p) dels[modules[i].positions[p]] = true;
        }
      }
      
      // Do all of the dups/deletions
      if (mod_mut_cnt > 0) {
        program_t new_program(program.GetInstLibPtr()); // Program we're writing to. (will be copied over.) 
        // for (rhead = 0; rhead < program.GetSize(); ++rhead) {
        rhead = 0;
        while ((int)new_program.GetSize() < expected_size) {
          const size_t cur_module = module_membership[rhead];
          // Did we find a new module?
          if (ilib.HasProperty(program[rhead].id, inst_prop_t::MODULE)) {
            // Should we duplicate cur_module?
            if (modules[cur_module].mut == ModuleMutType::DUP) {
              // We're supposed to duplicate current module.
              for (size_t i = 0; i < modules[cur_module].GetSize(); ++i) {
                new_program.PushInst(program[modules[cur_module].positions[i]]);
              }
            }
          }
          if (!dels[rhead]) {
            new_program.PushInst(program[rhead]);
          }
          ++rhead;
        }
        program = new_program;
      }

    }



    // Slip? -> If so, where?
    bool slip = false;
    bool slip_dup = false;
    bool slip_del = false;
    size_t slip_begin=0;
    size_t slip_end=0;
    // int slip_dup_size = 0;
    // int slip_del_size = 0;
    if (rnd.P(PER_PROG_SLIP)) {
      slip_begin = rnd.GetUInt(program.GetSize());
      slip_end = rnd.GetUInt(program.GetSize());
      slip_dup = slip_begin < slip_end;
      slip_del = slip_begin > slip_end;
      slip = slip_dup || slip_del; // If we may dup or del, we're slipping! (well, maybe - need to check that slip would't break constraints)
      if (slip) {
        const int slip_dup_size = 1 + (int)slip_end - (int)slip_begin;
        const int slip_del_size = 1 + (int)slip_begin - (int)slip_end;

        if (slip_dup && ((expected_size+slip_dup_size) > (int)MAX_PROGRAM_LEN)) { 
          // If slip-dup would break constraints, don't slip.
          slip = false; slip_dup=false; 
        } 
        if (slip_dup) { expected_size += slip_dup_size; }

        if (slip_del && ((expected_size-slip_del_size) < (int)MIN_PROGRAM_LEN)) { 
          // If slip-del would break constraints, don't slip.
          slip = false; slip_del = false; 
        } 
        if (slip_del) { expected_size -= slip_del_size; }

      }
    }

    program_t new_program(program.GetInstLibPtr());
    for (rhead = 0; rhead < (int)program.GetSize(); ++rhead) {
      // Check for slip.
      if (slip_dup && rhead == (int)slip_end) {
        // Copy over inst.
        new_program.PushInst(program[rhead]);
        // Duplicate slip segment.
        for (size_t si = slip_begin; si <= slip_end; ++si) {
          new_program.PushInst(program[si]);
        }
        mut_cnt++;
        continue;
      } else if (slip_del && rhead == (int)slip_end) {
        mut_cnt++;
        rhead = slip_begin;
        continue;
      }
      
      // Instruction deletion
      if (rnd.P(PER_INST_DEL) && ((expected_size-1)>=(int)MIN_PROGRAM_LEN)) {
        --expected_size;
        ++mut_cnt;
        continue;
      }

      // Copy instruction at read head
      const size_t whead = new_program.GetSize();
      new_program.PushInst(program[rhead]);

      // Instruction insertion
      if (rnd.P(PER_INST_INS) && ((expected_size+1)<=(int)MAX_PROGRAM_LEN)) {
        ++expected_size;
        ++mut_cnt;
        new_program.PushInst(TagLGP::GenRandTagGPInst(rnd, ilib));
      }

      // Instruction substitution
      if (rnd.P(PER_INST_SUB)) {
        ++mut_cnt;
        new_program[whead].id = rnd.GetUInt(ilib.GetSize());
      }

      // Instruction argument bit flips
      for (size_t arg = 0; arg < new_program[whead].arg_tags.size(); ++arg) {
        tag_t & tag = new_program[whead].arg_tags[arg];
        for (size_t k = 0; k < tag.GetSize(); ++k) {
          if (rnd.P(PER_BIT_FLIP)) {
            ++mut_cnt;
            tag.Toggle(k);
          }
        }
      }
    }

    program = new_program;
    return mut_cnt;
  }



};

#endif