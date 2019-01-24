#ifndef TAG_LINEAR_GP_UTILITIES_H
#define TAG_LINEAR_GP_UTILITIES_H

#include "tools/Random.h"
#include "TagLinearGP.h"
#include "Utilities.h"

namespace TagLGP {

  /// ====== Generators ======

  template<size_t TAG_WIDTH, size_t MAX_ARG_CNT>
  typename TagLinearGP_TW<TAG_WIDTH>::Program GenRandTagGPProgram(emp::Random & rnd,
                                                                  emp::Ptr<InstLib<TagLinearGP_TW<TAG_WIDTH>, MAX_ARG_CNT>> inst_lib,
                                                                  size_t MIN_LEN=1, size_t MAX_LEN=100) 
  {
    emp_assert(MIN_LEN < MAX_LEN);
    emp_assert(inst_lib->GetSize());
    using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
    using program_t = typename hardware_t::program_t;
    
    program_t prog(inst_lib);
    size_t N = rnd.GetUInt(MIN_LEN, MAX_LEN+1);
    for (size_t i = 0; i < N; ++i) {
      prog.PushInst(GenRandTagGPInst(rnd, *inst_lib));
    }
    return prog;
  }


  template<size_t TAG_WIDTH, size_t MAX_ARG_CNT> 
  typename TagLinearGP_TW<TAG_WIDTH>::Instruction GenRandTagGPInst(emp::Random & rnd, 
                                                                   const InstLib<TagLinearGP_TW<TAG_WIDTH>, MAX_ARG_CNT> & inst_lib) 
  {
    emp_assert(inst_lib.GetSize() > 0, "Instruction library must have at least one instruction definition before being used to generate a random instruction.");
    
    // using inst_lib_t = InstLib<TagLinearGP_TW<TAG_WIDTH>, MAX_ARG_CNT>;
    using hardware_t = TagLinearGP_TW<TAG_WIDTH>;
    using inst_t = typename hardware_t::inst_t;
    using tag_t = typename hardware_t::tag_t;

    // Generate tag arguments.
    emp::vector<tag_t> tags;
    for (size_t i = 0; i < MAX_ARG_CNT; ++i) {
      tags.emplace_back(rnd, 0.5);
    }
    
    return inst_t(rnd.GetUInt(inst_lib.GetSize()), tags);
  }

  // /// ====== Print ======
  // void PrintInst(const inst_t & inst, std::ostream & os=std::cout) {
  //   os << inst_lib->GetName(inst.id);
  //   os << "(";
  //   for (size_t t = 0; t < inst.arg_tags.size(); ++t) {
  //     if (t) os << ",";
  //     inst.arg_tags[t].Print(os);
  //   }
  // }

}

#endif