#ifndef TAG_LINEAR_GP_H
#define TAG_LINEAR_GP_H

#include <string>
#include <functional>
#include <unordered_set>
#include <cmath>

#include "base/vector.h"
#include "tools/BitSet.h"
#include "tools/map_utils.h"
#include "tools/Random.h"
#include "tools/string_utils.h"

#include "TagLinearGP_InstLib.h"
#include "Utilities.h"

namespace TagLGP {
  /////////////
  // TODO:
  // - Instruction asserts for instruction argument counts.
  // NOTES:
  // - WARNING - Can currently have infinitely recursing routine flows...
  // - FUTURE:
  //   - To solve problems with looping: make Flows condition-aware. Never need
  //     to close/reopen a flow for a loop. Instruction just passes lambda that
  //     defines the condition for the loop.
  //     - Should flows have an 'in-flow' ip map?
  /////////////



  /// TagLinearGP class.
  ///  - Description: 'Simple' lgp designed to handle program synthesis benchmarking
  ///    problems (i.e., designed to facilitate *more* arbitrary data types in memory).
  ///  - Simplifying non-features:
  ///   - No threads
  ///   - No events
  template<size_t TAG_WIDTH>
  class TagLinearGP_TW {
  public:
    struct Module;
    struct Instruction;
    struct Flow;
    struct CallState;
    class Memory;
    class MemoryValue;

    struct Program;

    using hardware_t = TagLinearGP_TW<TAG_WIDTH>;
    using tag_t = emp::BitSet<TAG_WIDTH>;
    using memory_t = Memory;
    using mem_type_t = typename MemoryValue::MemoryType;
    // using mem_pos_t = MemoryPosition;

    using module_t = Module;
    using inst_t = Instruction;
    using inst_lib_t = InstLib<hardware_t>;
    using inst_prop_t = typename inst_lib_t::InstProperty;
    using inst_seq_t = emp::vector<Instruction>;
    using args_t = emp::vector<tag_t>;

    using program_t = Program;

    // using fun_get_modules_t = std::function<emp::vector<module_t>(const program_t &)>;

    static constexpr size_t DEFAULT_MEM_SIZE = 16;
    static constexpr size_t DEFAULT_MAX_CALL_DEPTH = 128;
    static constexpr double DEFAULT_MIN_TAG_SPECIFICITY = 0.0;
    static constexpr size_t DEFAULT_MAX_STR_LEN = 4096;

    enum MemPosType { NUM=0, STR, VEC, ANY};
    enum FlowType { BASIC=0, LOOP, ROUTINE, CALL };

    struct Flow {
      FlowType type;
      size_t begin;
      size_t end;
      size_t iptr;
      size_t mptr;
      size_t iter;

      Flow(FlowType _type, size_t _begin, size_t _end,
           size_t _mptr, size_t _iptr) 
        : type(_type), 
          begin(_begin), 
          end(_end), 
          iptr(_iptr),
          mptr(_mptr),
          iter(0)
      { ; }
    };

    /// Struct - State
    /// - Description: maintains information about local call state.
    struct CallState {
      size_t mem_size;

      memory_t working_mem;
      memory_t input_mem;
      memory_t output_mem;

      emp::vector<Flow> flow_stack; ///< Stack of Flow (top is current).

      bool returnable;  ///< Can we return from this state? (or, are we trapped here forever!)
      bool circular;    ///< Does call have an implicit return at EOM? Or, is it circular?
      // {mem_size, returnable, circular, default_mem_val}
      CallState(size_t _mem_size=DEFAULT_MEM_SIZE, bool _returnable=true, bool _circular=false,
                const MemoryValue & _def_mem=MemoryValue())
        : mem_size(_mem_size),
          working_mem(_mem_size, _def_mem),
          input_mem(_mem_size, _def_mem),
          output_mem(_mem_size, _def_mem),
          flow_stack(),
          returnable(_returnable),
          circular(_circular)
      { ; }

      CallState(const CallState &) = default;
      CallState(CallState &&) = default;

      void Reset(const MemoryValue & def_mem=MemoryValue()) {
        working_mem.Reset(); // TODO - use default memory value?
        input_mem.Reset();
        output_mem.Reset();
      }

      bool IsReturnable() const { return returnable; }
      bool IsCircular() const { return circular; }
      bool IsFlow() const { return !flow_stack.empty(); }

      size_t GetMP() const { 
        if (flow_stack.size()) return flow_stack.back().mptr; 
        else return (size_t)-1;
      }
      
      size_t GetIP() const { 
        if (flow_stack.size()) return flow_stack.back().iptr;
        else return (size_t)-1;
      }

      emp::vector<Flow> & GetFlowStack() { return flow_stack; }
      Flow & GetTopFlow() { emp_assert(flow_stack.size()); return flow_stack.back(); }

      memory_t & GetWorkingMem() { return working_mem; }
      memory_t & GetInputMem() { return input_mem; }
      memory_t & GetOutputMem() { return output_mem; }

      void SetMP(size_t mp) { 
        if (flow_stack.size()) flow_stack.back().mptr = mp;
      }

      void SetIP(size_t ip) { 
        if (flow_stack.size()) flow_stack.back().iptr = ip; 
      }
      
      void AdvanceIP(size_t inc=1) { 
        if (flow_stack.size()) flow_stack.back().iptr += inc; 
      }
    
    };
    
    /// Structure of this class hails from C++ Primer 5th Ed. By Lippman, Lajoie, and Moo.
    class MemoryValue {
      public:
        enum MemoryType { NUM=0, STR=1};
        
      protected:
        MemoryType type;
        
        std::string default_str;
        double default_num;

        union {
          double num;
          std::string str;
        };

        void CopyUnion(const MemoryValue & in) {
          switch (in.type) {
            case MemoryType::NUM: {
              num = in.num;
              break;
            }
            case MemoryType::STR: {
              new (&str) std::string(in.str);
              break;
            }
            default:
              emp_assert(false, "Uknown memory value type.", in.type);
          }
        }

      public:
        MemoryValue() 
          : type(MemoryType::NUM), 
            default_str(""), 
            default_num(0),
            num(default_num) 
        { ; }

        MemoryValue(const MemoryValue & p) 
          : type(p.type),
            default_str(p.default_str),
            default_num(p.default_num)
        { 
          CopyUnion(p); // Copy union doesn't delete anything
        }

        MemoryValue & operator=(const MemoryValue & in) {
          using std::string;
          if (type == MemoryType::STR && in.type != MemoryType::STR) str.~string();
          if (type == MemoryType::STR && in.type == MemoryType::STR) {
            str = in.str;
          } else {
            CopyUnion(in);
          }
          type = in.type;
          return *this;
        }

        MemoryValue & operator=(const std::string & in) {
          if (type == MemoryType::STR) {
            str = in;
          } else {
            new (&str) std::string(in);
          }
          type = MemoryType::STR;
          return *this;
        }

        MemoryValue & operator=(double in) {
          using std::string;
          if (type == MemoryType::STR) str.~string();
          num = in;
          type = MemoryType::NUM;
          return *this;
        }

        bool operator==(const MemoryValue & in) const {
          if (in.type == type) {
            switch (type) {
              case MemoryType::STR: return str == in.str;
              case MemoryType::NUM: return num == in.num;
            }
          }
          return false;
        }

        bool operator!=(const MemoryValue & in) const {
          return !(*this == in);
        }

        // If the union holds a non built-in type, we need to destroy it.
        ~MemoryValue() { 
          using std::string;
          if (type == MemoryType::STR) str.~string(); 
        }

        const std::string & GetDefaultStr() const { return default_str; }
        double GetDefaultNum() const { return default_num; }
        
        std::string & GetStr() { 
          if (type == MemoryType::STR) {
            return str;
          } else {
            emp_assert(false, "Requesting string from non-string memory value.");
            return default_str;
          }
        }

        double GetNum() const {
          if (type == MemoryType::NUM) {
            return num;
          } else {
            // TODO - Can str be converted to number?
            emp_assert(false, "Requesting num from non-num memory value.");
            return default_num;
          }
        }

        MemoryType GetType() const { return type; }

        void SetDefaultStr(const std::string & ds) { default_str = ds; }
        void SetDefaultNum(double dn) { default_num = dn; }

        void Print(std::ostream & os=std::cout) const {
          switch (type) {
            case MemoryType::NUM: os << num; break;
            case MemoryType::STR: os << "\"" << str << "\""; break;
            default: os << "UNKOWN TYPE"; break;
          }
        }
    };

    /// Memory class to manage TagLGP memory.
    class Memory {
      public:
        
        struct MemoryPosition {
          bool is_vector;
          bool set;
          emp::vector<MemoryValue> pos;

          MemoryPosition(const MemoryValue & val=MemoryValue(), bool is_vec=false)
            : is_vector(is_vec), 
              set(false),
              pos(1, val) { ; }
          MemoryPosition(const MemoryPosition &) = default;

          MemoryPosition & operator=(const MemoryPosition &) = default;
          MemoryPosition & operator=(MemoryPosition &&) = default;

          bool operator==(const MemoryPosition & in) const {
            if (in.is_vector == is_vector) { return pos == in.pos; }
            return false;
          }

          bool operator!=(const MemoryPosition & in) const {
            return !(*this == in);
          }

          void Print(std::ostream & os=std::cout) const {
            if (is_vector) {
              os << "[";
              for (size_t i = 0; i < pos.size(); ++i) {
                if (i) os << ",";
                pos[i].Print(os);
              }
              os << "]";
            } else {
              pos[0].Print(os);
            }
          }
        };

      protected:
        emp::vector<MemoryPosition> memory;
        
      public:
        Memory(size_t size, const MemoryValue & default_value=MemoryValue()) 
          : memory(size, {default_value})
        { ; }

        Memory(const Memory &) = default;
        Memory(Memory &&) = default;

        void Reset(const MemoryValue & default_value=MemoryValue()) {
          const size_t size = memory.size();
          memory.clear();
          memory.resize(size, {default_value});
        }

        void Resize(size_t size, const MemoryValue & default_value=MemoryValue()) {
          memory.resize(size, {default_value});
        }

        bool IsSet(size_t id) const {
          emp_assert(id < memory.size());
          return memory[id].set;
        }

        bool IsVec(size_t id) const {
          emp_assert(id < memory.size());
          return memory[id].is_vector;
        }

        MemPosType GetPosType(size_t id) const {
          if (memory[id].is_vector) return MemPosType::VEC;
          else if (memory[id].pos[0].GetType() == mem_type_t::NUM) return MemPosType::NUM;
          else if (memory[id].pos[0].GetType() == mem_type_t::STR) return MemPosType::STR;
          return MemPosType::ANY;
        }

        size_t GetSize() const { return memory.size(); }

        const MemoryPosition & GetPos(size_t id) const {
          emp_assert(id < memory.size());
          return memory[id];
        }

        MemoryValue & AccessVal(size_t id) {
          emp_assert(id < memory.size());
          emp_assert(memory[id].pos.size());
          memory[id].set = true;          // Memory becomes set on access by reference.
          return memory[id].pos[0];          
        }

        emp::vector<MemoryValue> & AccessVec(size_t id) {
          emp_assert(id < memory.size());
          emp_assert(memory[id].is_vector);
          memory[id].set = true;          // Memory becomes set on access by reference.
          return memory[id].pos;
        }

        void Swap(size_t a, size_t b) {
          emp_assert(a < memory.size());
          emp_assert(b < memory.size());
          std::swap(memory[a], memory[b]);
        }

        /// Set memory[id] = other
        void Set(size_t id, const MemoryPosition & other) {
          emp_assert(id < memory.size());
          memory[id] = other;
          memory[id].set = true;
        }
        
        /// Set memory[id] = vector of MemoryValues
        void Set(size_t id, const emp::vector<MemoryValue> & val_vec) {
          emp_assert(id < memory.size());
          // emp_assert(val_vec.size());
          memory[id].pos = val_vec;
          memory[id].set = true;
          memory[id].is_vector = true;
        }

        /// Set memory[id] = given MemoryValue
        void Set(size_t id, const MemoryValue & val) {
          emp_assert(id < memory.size());
          memory[id].pos.resize(1); // Memory position is now just a value, shrink appropriately.
          memory[id].pos[0] = val;
          memory[id].set = true;
          memory[id].is_vector = false;
        }

        void Set(size_t id, double val) {
          emp_assert(id < memory.size());
          memory[id].pos.resize(1);
          memory[id].pos[0] = val;
          memory[id].set = true;
          memory[id].is_vector = false;
        }

        void Set(size_t id, const emp::vector<double> & val_vec) {
          emp_assert(id < memory.size());
          // emp_assert(val_vec.size());
          memory[id].pos.resize(val_vec.size());
          for (size_t i = 0; i < memory[id].pos.size(); ++i) memory[id].pos[i] = val_vec[i];
          memory[id].set = true;
          memory[id].is_vector = true;
        }

        void Set(size_t id, const emp::vector<int> & val_vec) {
          emp_assert(id < memory.size());
          // emp_assert(val_vec.size());
          memory[id].pos.resize(val_vec.size());
          for (size_t i = 0; i < memory[id].pos.size(); ++i) memory[id].pos[i] = (double)val_vec[i];
          memory[id].set = true;
          memory[id].is_vector = true;
        }

        void Set(size_t id, const std::string & str) {
          emp_assert(id < memory.size());
          memory[id].pos.resize(1);
          memory[id].pos[0] = str;
          memory[id].set = true;
          memory[id].is_vector = false;
        }

        void Set(size_t id, const emp::vector<std::string> & str_vec) {
          emp_assert(id < memory.size());
          // emp_assert(str_vec.size());
          memory[id].pos.resize(str_vec.size());
          for (size_t i = 0; i < memory[id].pos.size(); ++i) memory[id].pos[i] = str_vec[i];
          memory[id].set = true;
          memory[id].is_vector = true;
        }

        // TODO - Append, Resize - (case - to zero - unset, resize 1 default), Clear()

        void Print(std::ostream & os=std::cout) const {
          os << "{";
          for (size_t i = 0; i < memory.size(); ++i) {
            if (i) os << ",";
            os << i << ":";
            memory[i].Print(os);
          }
          os << "}";
        }

        void PrintSetMem(std::ostream & os=std::cout) const {
          os << "{";
          for (size_t i = 0; i < memory.size(); ++i) {
            if (!memory[i].set) { continue; }
            if (i) os << ",";
            os << i << ":";
            memory[i].Print(os);
          }
          os << "}";
        }

    };

    /// Module definition.
    struct Module {
      size_t id;
      size_t begin;   ///< First instruction in module (will be executed first).
      size_t end;     ///< Instruction pointer value this module returns (or loops back) on (1 past last instruction that is part of this module).
      tag_t tag;      ///< Module tag. Used to call/reference module.
      std::unordered_set<size_t> in_module; ///< instruction positions belonging to this module.

      Module(size_t _id, size_t _begin=0, size_t _end=0, tag_t _tag=tag_t())
        : id(_id), begin(_begin), end(_end), tag(_tag) { ; }

      size_t GetLen() const { return in_module.size(); }

      bool InModule(size_t ip) const {
        return emp::Has(in_module, ip);
      }
    };
    
    struct Instruction {
      size_t id;
      emp::vector<tag_t> arg_tags;

      Instruction(size_t _id=0, const emp::vector<tag_t> & _arg_tags=emp::vector<tag_t>())
        : id(_id), arg_tags(_arg_tags) { ; }
      
      Instruction(const Instruction &) = default;
      Instruction(Instruction &&) = default;

      Instruction & operator=(const Instruction &) = default;
      Instruction & operator=(Instruction &&) = default;

      bool operator==(const Instruction & in) const {
        return id == in.id && arg_tags == in.arg_tags;
      }
      bool operator!=(const Instruction & in) const { return !(*this == in); }

      bool operator<(const Instruction & other) const {
          return std::tie(id, arg_tags) < std::tie(other.id, other.arg_tags);
      }

      size_t GetNumArgs() const { return arg_tags.size(); }

      emp::vector<tag_t> & GetArgTags() {
        return arg_tags;
      }

      tag_t & GetArgTag(size_t i) {
        emp_assert(i < arg_tags.size());
        return arg_tags[i];
      }

      void Set(size_t _id, const emp::vector<tag_t> & _args) {
        id = _id;
        arg_tags = _args;
      }

      void Set(const Instruction & other) {
        id = other.id;
        arg_tags = other.arg_tags;
      } 

    };

    struct Program {
      emp::Ptr<const inst_lib_t> inst_lib;  ///< Pointer to instruction library associated with this program.
      inst_seq_t program;     ///< Programs are linear sequences of instructions.

      Program(emp::Ptr<const inst_lib_t> _ilib, const inst_seq_t & _prog=inst_seq_t())
        : inst_lib(_ilib), program(_prog){ ; }
      
      Program(const Program &) = default;

      bool operator==(const Program & in) const { return program == in.program; }
      bool operator!=(const Program & in) const { return !(*this == in); }
      bool operator<(const Program & other) const { return program < other.program; }

      /// Allow program's instruction sequence to be indexed as if a vector.
      Instruction & operator[](size_t id) { 
        emp_assert(id < program.size());
        return program[id]; 
      }

      /// Allow program's instruction sequence to be indexed as if a vector.
      const Instruction & operator[](size_t id) const {
        emp_assert(id < program.size());
        return program[id];
      }

      // Clear program's instruction sequence.
      void Clear() { program.clear(); }

      /// Get program size.
      size_t GetSize() const { return program.size(); }

      /// Get a pointer to const instruction library.
      emp::Ptr<const inst_lib_t> GetInstLibPtr() const { return inst_lib; }

      /// Get const reference to instruction library.
      const inst_lib_t & GetInstLib() const { return *inst_lib; }

      /// Is given position a valid position in this program?
      bool ValidPosition(size_t pos) const { return pos < GetSize(); }

      /// Set program's instruction sequence to the one given.
      void SetProgram(const inst_seq_t & p) { program = p; }

      /// Push new instruction to program by instruction ID.
      void PushInst(size_t id, const args_t & args) {
        program.emplace_back(id, args);
      }

      /// Push new instruction to program by instruction name.
      void PushInst(const std::string & name, const args_t & args) {
        emp_assert(inst_lib->IsInst(name), "Uknown instruction name", name);
        PushInst(inst_lib->GetID(name), args);
      }

      /// Push given instruction onto program.
      void PushInst(const inst_t & inst) {
        program.emplace_back(inst);
      }

      /// Overwrite instruction in sequence (@ position pos).
      void SetInst(size_t pos, size_t id, const args_t & args) {
        emp_assert(pos < GetSize());
        program[pos].Set(id, args);
      }

      /// Overwrite instruction in sequence (@ position pos).
      void SetInst(size_t pos, std::string & name, const args_t & args) {
        SetInst(pos, inst_lib->GetID(name), args);
      }

      /// Overwrite instruction in sequence (@ position pos).
      void SetInst(size_t pos, const inst_t & inst) {
        program[pos] = inst;
      }

      void PrintInst(const inst_t & inst, std::ostream & os=std::cout) const {
        // Print instruction
        os << inst_lib->GetName(inst.id);
        os << "(";
        for (size_t t = 0; t < inst.arg_tags.size(); ++t) {
          if (t) os << ",";
          inst.arg_tags[t].Print(os);
        }
        os << ")";
      }

      /// Plain-print program.
      void Print(std::ostream & os=std::cout) const {
        for (size_t i = 0; i < program.size(); ++i) {
          const inst_t & inst = program[i];
          PrintInst(inst, os);
          os << "\n";
        }
      }

      /// Print full program as valid CSV entry. 
      /// "[InstName(1111,1111,1111),InstName(1111,1111,1111),...]"
      void PrintCSVEntry(std::ostream & os=std::cout) const {
        os << "\"[";
        for (size_t i = 0; i < program.size(); ++i) {
          if (i) os << ",";
          const inst_t & inst = program[i];
          PrintInst(inst, os);
        }
        os << "]\"";
      }

    };

  protected:

    emp::Ptr<emp::Random> random_ptr;
    bool random_owner;

    program_t program;
    
    emp::vector<module_t> modules;

    MemoryValue default_mem_val;

    size_t mem_size;
    emp::vector<tag_t> mem_tags;
    // emp::vector<tag_t> working_mem_tags;
    // emp::vector<tag_t> input_mem_tags;
    // emp::vector<tag_t> output_mem_tags;

    memory_t global_mem;

    emp::vector<CallState> call_stack;

    size_t max_call_depth;
    double min_tag_specificity;

    bool is_executing;

    void OpenFlow_CALL(CallState & state, const Module & module) {
      emp_assert(state.GetFlowStack().size() == 0); // Should only put call flows at bottom of stack.
      OpenFlow(state, {FlowType::CALL, module.begin, module.end, module.id, module.begin});
    }

    void CloseFlow_BASIC(CallState & state, bool implicit) {
      emp_assert(state.IsFlow());
      // Closing BASIC flow:
      // - Pop basic flow from flow stack, passing IP and MP down.
      const size_t ip = state.GetTopFlow().iptr;
      const size_t mp = state.GetTopFlow().mptr;
      state.GetFlowStack().pop_back();
      if (state.IsFlow()) {
        Flow & top = state.GetTopFlow();
        top.iptr = ip;
        top.mptr = mp;
      }
    }

    void CloseFlow_LOOP(CallState & state, bool implicit) {
      emp_assert(state.IsFlow());
      // Closing a LOOP flow:
      const size_t loop_begin = state.GetTopFlow().begin;
      const size_t mp = state.GetTopFlow().mptr;
      
      state.GetFlowStack().pop_back();

      if (state.IsFlow()) {
        Flow & top = state.GetTopFlow();
        top.iptr = loop_begin; // Loop parent flow back to beginning of while loop.
        top.mptr = mp;
        ++top.iter;
      }
    }

    void CloseFlow_ROUTINE(CallState & state, bool implicit) {
      emp_assert(state.IsFlow());
      // Closing a ROUTINE flow:
      // - Pop ROUTINE flow from flow stack.
      // - We don't pass IP and MP down (those on the lower stack are where we
      //   should return to).
      state.GetFlowStack().pop_back();
    }

    void CloseFlow_CALL(CallState & state, bool implicit) {
      emp_assert(state.IsFlow());
      // Closing a CALL flow:
      // - Pop call flow from flow stack.
      // - No need to pass IP and MP down (presumably, this was the bottom of the
      //   flow stack).
      if (implicit && state.IsCircular()) {
        Flow & top = state.GetTopFlow();
        top.iptr = top.begin;
      } else if (!state.IsReturnable() && state.GetFlowStack().size() == 1) {
        Flow & top = state.GetTopFlow();
        emp_assert(top.type == FlowType::CALL, "Expected flow type here to be CALL.");
        top.iptr = top.begin; // Wrap
      } else {
        state.GetFlowStack().pop_back();
      }
    }


    void BreakFlow_BASIC(CallState & state) {
      emp_assert(state.IsFlow());
      // Break out of basic flow:
      const size_t flow_end = state.GetTopFlow().end;
      state.GetFlowStack().pop_back();
      if (state.IsFlow()) {
        Flow & top = state.GetTopFlow();
        top.iptr = flow_end;
        top.iter = 0;
        if (ValidPosition(top.mptr, top.iptr)) state.AdvanceIP(); // Skip over Close
      }
    }

    void BreakFlow_LOOP(CallState & state) {
      emp_assert(state.IsFlow());
      // Break out of loop flow:
      const size_t flow_end = state.GetTopFlow().end;
      state.GetFlowStack().pop_back();
      if (state.IsFlow()) {
        Flow & top = state.GetTopFlow();
        top.iptr = flow_end;
        top.iter = 0;
        if (ValidPosition(top.mptr, top.iptr)) state.AdvanceIP(); // Skip over Close
      }
    }

  public:

    TagLinearGP_TW(emp::Ptr<const inst_lib_t> _ilib,
                  emp::Ptr<emp::Random> rnd=nullptr)
      : random_ptr(rnd), random_owner(false),
        program(_ilib),
        modules(),
        default_mem_val(),
        mem_size(DEFAULT_MEM_SIZE),
        // global_mem_tags(mem_size),
        // working_mem_tags(mem_size),
        // input_mem_tags(mem_size),
        // output_mem_tags(mem_size),
        mem_tags(mem_size),
        global_mem(mem_size, default_mem_val),
        call_stack(),
        max_call_depth(DEFAULT_MAX_CALL_DEPTH),
        min_tag_specificity(DEFAULT_MIN_TAG_SPECIFICITY),
        is_executing(false)
    { 
      // If no random number generator is provided, create one (taking ownership).
      if (!rnd) NewRandom(); 
    }

    TagLinearGP_TW(inst_lib_t & _ilib, 
                  emp::Ptr<emp::Random> rnd=nullptr)
      : TagLinearGP_TW(&_ilib, rnd) { ; }
    
    TagLinearGP_TW(const hardware_t & in) 
      : random_ptr(nullptr), random_owner(false),
        program(in.program),
        modules(in.modules),
        default_mem_val(in.default_mem_val),
        mem_size(in.mem_size),
        // global_mem_tags(in.global_mem_tags),
        // working_mem_tags(in.working_mem_tags),
        // input_mem_tags(in.input_mem_tags),
        // output_mem_tags(in.output_mem_tags),
        mem_tags(in.mem_tags),
        global_mem(in.global_mem),
        call_stack(in.call_stack),
        max_call_depth(in.max_call_depth),
        min_tag_specificity(in.min_tag_specificity),
        is_executing(in.is_executing)
    {
      if (in.random_owner) NewRandom();
      else random_ptr = in.random_ptr;
    }

    ~TagLinearGP_TW() { if (random_owner) random_ptr.Delete(); }

    // ---------------------------- Hardware configuration ----------------------------

    /// Set program for this hardware object.
    /// After updating hardware's program, run UpdateModules to update module definitions.
    void SetProgram(const program_t & _program) { 
      emp_assert(!is_executing);
      program = _program; 
      UpdateModules();
    }

    /// Set hardware memory size (number of memory positions per memory buffer)
    void SetMemSize(size_t size) {
      mem_size = size;
      global_mem.Resize(mem_size, default_mem_val);
      mem_tags.resize(mem_size, tag_t());
      for (size_t i = 0; i < call_stack.size(); ++i) {
        CallState & state = call_stack[i];
        
        state.input_mem.Resize(mem_size, default_mem_val);
        state.working_mem.Resize(mem_size, default_mem_val);
        state.output_mem.Resize(mem_size, default_mem_val);
      }
    }

    void SetMemTags(const emp::vector<tag_t> & tags) {
      emp_assert(tags.size() == mem_size);
      mem_tags = tags;
    }

    void SetMinTagSpecificity(double val) { min_tag_specificity = val; }

    void SetMaxCallDepth(size_t depth) { max_call_depth = depth; }

    // ---------------------------- Hardware control ----------------------------
    /// Reset everything, including the program.
    void Reset() {
      emp_assert(!is_executing);
      ResetHardware();
      modules.clear();
      program.Clear();
    }

    /// Reset only hardware, not program.
    /// Not allowed to reset hardware during execution.
    void ResetHardware() {
      emp_assert(!is_executing);
      global_mem.Reset(default_mem_val);
      call_stack.clear();
      is_executing = false;
    }

    /// Update modules.
    void UpdateModules() {
      // Grab reference to program's instruction library.
      const inst_lib_t & ilib = program.GetInstLib();
      // Clear out current module definitions.
      modules.clear();
      // Scan program for module definitions.
      std::unordered_set<size_t> dangling_instructions;
      for (size_t pos = 0; pos < program.GetSize(); ++pos) {
        inst_t & inst = program[pos];
        // Is this a module definition?
        if (ilib.HasProperty(inst.id, inst_prop_t::MODULE)) {
          if (modules.size()) { modules.back().end = pos; } 
          const size_t mod_id = modules.size();
          modules.emplace_back(mod_id, (pos+1)%program.GetSize(), -1, inst.arg_tags[0]);
        } else {
          // Not a new module definition.
          if (modules.size()) { modules.back().in_module.emplace(pos); }
          else { dangling_instructions.emplace(pos); }
        }
      }
      // Take care of dangling instructions and add default module if necessary.
      if (modules.size()) {                               // if we've found at least one module, set it's end point.
        if (modules[0].begin == 0) modules.back().end = program.GetSize();
        else modules.back().end = modules[0].begin-1;
      } else {                                            // found no modules, add a default module.
        // Default module starts at beginning of program and ends at the end.
        modules.emplace_back(0, 0, program.GetSize(), tag_t());
      }
      for (size_t val : dangling_instructions) modules.back().in_module.emplace(val);
    }

    // ---------------------------- Hardware execution ----------------------------
    /// Process a single instruction, provided by the caller.
    void ProcessInst(const inst_t & inst) { program.GetInstLib().ProcessInst(*this, inst); }

    /// Advance hardware by a single instruction.
    void SingleProcess() {
      emp_assert(program.GetSize()); // Must have a non-empty program to advance the hardware.

      is_executing = true;

      // If there's a call state on the call stack, execute an instruction in
      // that call.

      // Repeat:
      while (call_stack.size()) {
        CallState & state = call_stack.back();
        // todo - check validity of ip/mp
        if (state.IsFlow()) { // TODO - may need to switch this to flag (e.g., state.done)
          Flow & top_flow = state.GetTopFlow();
          size_t ip = top_flow.iptr;
          size_t mp = top_flow.mptr;
          emp_assert(mp < modules.size());
          if (modules[mp].InModule(ip)) { // Valid IP
            // First, increment flow's IP by 1. This may invalidate the IP, but
            // that's okay.
            ++top_flow.iptr;
            // Next, run instruction @ ip.
            GetInstLib().ProcessInst(*this, program[ip]);
          } else if (ip >= program.GetSize() && modules[mp].InModule(0) && modules[mp].end < modules[mp].begin) {  // If module wraps.
            ip = 0;
            top_flow.iptr = 1;
            GetInstLib().ProcessInst(*this, program[ip]);
          } else { 
            CloseFlow(state, true);
            continue;
          }
        } else {
          // Return from CallState
          ReturnCall(true);
          continue; // Implicit returns are free for now...
        }
        break;
      }
      is_executing = false;
    }

    void OpenFlowRoutine(CallState & state, const Module & module) {
      OpenFlow(state, {FlowType::ROUTINE, module.begin, module.end, module.id, module.begin});
    }

    void OpenFlowRoutine(CallState & state, size_t mID) {
      emp_assert(mID < modules.size());
      const Module & module = modules[mID];
      OpenFlow(state, {FlowType::ROUTINE, module.begin, module.end, module.id, module.begin});
    }

    void OpenFlow(CallState & state, FlowType type, 
                  size_t begin, size_t end, size_t mptr, size_t iptr) {
      OpenFlow(state, {type, begin, end, mptr, iptr});
    }

    void OpenFlow(CallState & state, const Flow & new_flow) {
      emp_assert(!(state.GetFlowStack().size() > 0 && new_flow.type == FlowType::CALL), "Call flows are only allowed at bottom of flow stack.");
      state.GetFlowStack().emplace_back(new_flow);
    }

    void CloseFlow(CallState & state, bool implicit=false) {
      if (!state.IsFlow()) return;
      switch (state.GetTopFlow().type) {
        case BASIC: CloseFlow_BASIC(state, implicit); break;
        case LOOP: CloseFlow_LOOP(state, implicit); break;
        case ROUTINE: CloseFlow_ROUTINE(state, implicit); break;
        case CALL: CloseFlow_CALL(state, implicit); break;
        default:
          std::cout << "Uknown flow type (" << state.GetTopFlow().type << ")!" << std::endl;
      }
    }

    void BreakFlow(CallState & state) {
      if (state.IsFlow()) {
        switch (state.GetTopFlow().type) {
          case BASIC: BreakFlow_BASIC(state); break;
          case LOOP: BreakFlow_LOOP(state); break;
          case ROUTINE: break;
          case CALL: break;
          default:
            std::cout << "Uknown flow type (" << state.GetTopFlow().type << ")!" << std::endl;
        }
      }
    }

    void CallModule(const tag_t & tag, double threshold=0.0, bool returnable=true, bool circular=false) {
      const size_t target_module = FindBestModuleMatch(tag, threshold);
      if (target_module < modules.size()) {
        CallModule(target_module, returnable, circular);
      }
    }

    void CallModule(size_t mID, bool returnable=true, bool circular=false) {
      emp_assert(mID < modules.size());
      // Are we at max depth? If so, call fails.
      if (call_stack.size() >= max_call_depth) return;
      // Push new state onto stack.
      call_stack.emplace_back(mem_size, returnable, circular, default_mem_val);
      // Open call flow on stack w/called module.
      OpenFlow_CALL(call_stack.back(), modules[mID]);
      // If there's at least one call state before this one, configure new
      // state's memory appropriately; otherwise, leave defaults.
      if (call_stack.size() > 1) {
        CallState & new_state = call_stack.back();
        CallState & caller_state = call_stack[call_stack.size() - 2]; 
        Memory & new_input_mem = new_state.GetInputMem();
        Memory & caller_working_mem = caller_state.GetWorkingMem();
        for (size_t i = 0; i < mem_size; ++i) {
          if (caller_working_mem.IsSet(i)) { new_input_mem.Set(i, caller_working_mem.GetPos(i)); }
        }
      }
    }

    void ReturnCall(bool implicit=false) {
      // bool returnable;  ///< Can we return from this state? (or, are we trapped here forever!)
      // bool circular;    ///< Does call have an implicit return at EOM? Or, is it circular?
      if (call_stack.empty()) return; // Nothing to return from.
      CallState & returning_state = GetCurCallState();

      // // No returning from non-returnable call state.
      if (!returning_state.IsReturnable()) {
        while (returning_state.GetFlowStack().size() > 1) {
          CloseFlow(returning_state, implicit);
        }
        emp_assert(returning_state.GetTopFlow().type == FlowType::CALL, "Bottom flow should be the call.", returning_state.GetTopFlow().type);
        CloseFlow(returning_state, implicit);
        return;
      }  

      // Is there anything to return to?
      if (call_stack.size() > 1) {
        // If so, copy returning state's output memory into caller state's local memory.
        CallState & caller_state = call_stack[call_stack.size() - 2];
        Memory & out_mem = returning_state.GetOutputMem();
        Memory & working_mem = caller_state.GetWorkingMem();
        for (size_t i = 0; i < mem_size; ++i) {
          if (out_mem.IsSet(i)) { working_mem.Set(i, out_mem.GetPos(i)); } 
        }
      }

      // Pop returning state from call stack.
      call_stack.pop_back();
    }
    
    // ---------------------------- Accessors ----------------------------
    /// Get instruction library associated with hardware's program.
    emp::Ptr<const inst_lib_t> GetInstLibPtr() const { return program.GetInstLibPtr(); }
    const inst_lib_t & GetInstLib() const { return program.GetInstLib(); }

    /// Get random number generator.
    emp::Random & GetRandom() { return *random_ptr; }
    
    /// Get pointer to random number generator.
    emp::Ptr<emp::Random> GetRandomPtr() { return random_ptr; }

    /// Get program loaded on this hardware.
    const program_t & GetProgram() const { return program; }
    program_t & GetProgram() { return program; }

    /// Get the minimum tag specificity (i.e., required similarity between two tags
    /// for a successful reference to occur).
    double GetMinTagSpecificity() const { return min_tag_specificity; }

    /// Get maximum allowed call depth (maximum number of call states allowed on
    /// the call stack at any one time).
    size_t GetMaxCallDepth() const { return max_call_depth; }

    /// Get memory size. How many memory positions are available in input, output,
    /// working, and global memory.
    size_t GetMemSize() const { return mem_size; }

    size_t GetModuleCnt() const { return modules.size(); }

    /// Get global memory vector.
    memory_t & GetGlobalMem() { return global_mem; }
    const memory_t & GetGlobalMem() const { return global_mem; }

    /// memory tag accessors
    emp::vector<tag_t> & GetMemTags() { return mem_tags; }
    const emp::vector<tag_t> & GetMemTags() const { return mem_tags; }
    // emp::vector<tag_t> & GetGlobalMemTags() { return global_mem_tags; }
    // const emp::vector<tag_t> & GetGlobalMemTags() const { return global_mem_tags; }
    // emp::vector<tag_t> & GetWorkingMemTags() { return working_mem_tags; }
    // const emp::vector<tag_t> & GetWorkingMemTags() const { return working_mem_tags; }
    // emp::vector<tag_t> & GetInputMemTags() { return input_mem_tags; }
    // const emp::vector<tag_t> & GetInputMemTags() const { return input_mem_tags; }
    // emp::vector<tag_t> & GetOutputMemTags() { return output_mem_tags; }
    // const emp::vector<tag_t> & GetOutputMemTags() const { return output_mem_tags; }

    CallState & GetCurCallState() {
      emp_assert(call_stack.size(), "Cannot query for current call state if call stack is empty.");
      return call_stack.back();
    }

    const CallState & GetCurCallState() const { 
      emp_assert(call_stack.size(), "Cannot query for current call state if call stack is empty.");
      return call_stack.back();
    }

    /// Get size of call stack.
    size_t GetCallStackSize() const { return call_stack.size(); }

    // ---------------------------- Hardware utilities ----------------------------
    void NewRandom(int seed=-1) {
      if (random_owner) random_ptr.Delete();
      else random_ptr = nullptr;
      random_ptr = emp::NewPtr<emp::Random>(seed);
      random_owner = true;
    }

    size_t FindBestMemoryMatch_SMC(const Memory & mem, const tag_t & tag, double threshold=0.0, MemPosType mem_type=MemPosType::ANY) { 
      emp::vector<size_t> best_matches;
      for (size_t i = 0; i < mem_tags.size(); ++i) {
        // ANY, VEC, STR, NUM
        if (mem_type == MemPosType::ANY || mem_type == mem.GetPosType(i)) {
          double match = emp::SimpleMatchCoeff(tag, mem_tags[i]);
          if (match == threshold) best_matches.emplace_back(i);
          else if (match > threshold) { // If distance is closer.
            best_matches.resize(1);
            best_matches[0] = i;
            threshold = match;
          }
        }
      }
      if (best_matches.size()) {
        return best_matches[0];
      } else {
        return (size_t)-1;
      }
    }

    size_t FindBestMemoryMatch_SMC(const Memory & mem, const tag_t & tag, std::unordered_set<MemPosType> mem_types, double threshold=0.0) { 
      emp::vector<size_t> best_matches;
      for (size_t i = 0; i < mem_tags.size(); ++i) {
        // ANY, VEC, STR, NUM
        if (emp::Has(mem_types, mem.GetPosType(i))) {
          double match = emp::SimpleMatchCoeff(tag, mem_tags[i]);
          if (match == threshold) best_matches.emplace_back(i);
          else if (match > threshold) { // If distance is closer.
            best_matches.resize(1);
            best_matches[0] = i;
            threshold = match;
          }
        }
      }
      if (best_matches.size()) {
        return best_matches[0];
      } else {
        return (size_t)-1;
      }
    }
  
    /// Return best matching memory
    /// TODO - configurable tie-breaking procedure
    size_t FindBestMemoryMatch(const Memory & mem, const tag_t & tag, double threshold=0.0, MemPosType mem_type=MemPosType::ANY) { 
      double dist_thresh = TAG_WIDTH - (threshold * (double)TAG_WIDTH);
      emp::vector<size_t> best_matches;
      for (size_t i = 0; i < mem_tags.size(); ++i) {
        // ANY, VEC, STR, NUM
        if (mem_type == MemPosType::ANY || mem_type == mem.GetPosType(i)) {
          // double match = emp::SimpleMatchCoeff(tag, mem_tags[i]);
          double dist = (double)HammingDist(tag, mem_tags[i]);
          if (dist == dist_thresh) best_matches.emplace_back(i);
          else if (dist < dist_thresh) {
            best_matches.resize(1);
            best_matches[0] = i;
            dist_thresh = dist;
          }
        }
      }
      if (best_matches.size()) {
        return best_matches[0];
      } else {
        return (size_t)-1;
      }
    }

    size_t FindBestMemoryMatch(const Memory & mem, const tag_t & tag, std::unordered_set<MemPosType> mem_types, double threshold=0.0) { 
      double dist_thresh = TAG_WIDTH - (threshold * (double)TAG_WIDTH);
      // TODO - type checking
      emp::vector<size_t> best_matches;
      for (size_t i = 0; i < mem_tags.size(); ++i) {
        // ANY, VEC, STR, NUM
        if (emp::Has(mem_types, mem.GetPosType(i))) {
          // double match = emp::SimpleMatchCoeff(tag, mem_tags[i]);
          double dist = (double)HammingDist(tag, mem_tags[i]);
          if (dist == dist_thresh) best_matches.emplace_back(i);
          else if (dist < dist_thresh) { // If distance is closer.
            best_matches.resize(1);
            best_matches[0] = i;
            dist_thresh = dist;
          }
        }
      }
      if (best_matches.size()) {
        return best_matches[0];
      } else {
        return (size_t)-1;
      }
    }

    size_t FindBestModuleMatch(const tag_t & tag, double threshold=0.0) { 
      emp::vector<size_t> best_matches;
      for (size_t i = 0; i < modules.size(); ++i) {
        double match = emp::SimpleMatchCoeff(tag, modules[i].tag);
        if (match == threshold) best_matches.emplace_back(i);
        else if (match > threshold) {
          best_matches.resize(1);
          best_matches[0] = i;
          threshold = match;
        }
      }
      if (best_matches.size()) {
        return best_matches[0];
      } else {
        return (size_t)-1;
      }
    }

    bool IsValidMemPos(size_t pos) const { return pos < mem_size; }

    bool ValidPosition(size_t mp, size_t ip) const {
      if (mp < modules.size()) {
        if (modules[mp].InModule(ip)) return true;
      }
      return false;
    }

    size_t FindEndOfFlow(size_t mp, size_t ip) {
      emp_assert(mp < modules.size());
      // std::cout << "EOF(mp=" << mp << "," << "ip=" << ip << ")" << std::endl;
      const inst_lib_t & ilib = program.GetInstLib();
      int depth_counter = 1;
      std::unordered_set<size_t> seen;
      while (true) {
        // std::cout << ">> ip=" << ip << "; d=" << depth_counter << std::endl;
        if (!ValidPosition(mp, ip)) break;
        const inst_t & inst = program[ip];
        if (ilib.HasProperty(inst.id, inst_prop_t::BEGIN_FLOW)) {
          ++depth_counter;
        } else if (ilib.HasProperty(inst.id, inst_prop_t::END_FLOW)) {
          --depth_counter;
          // If depth counter is ever 0 after subtracting, we've found the end
          // to the initial flow.
          if (depth_counter == 0) break;
        }
        seen.emplace(ip); 
        ++ip; 
        if (ip >= program.GetSize() && seen.size() < modules[mp].GetLen()) ip %= program.GetSize(); // Wrap ip around
      }
      // std::cout << "EOF=" << ip << std::endl;
      return ip;
    }

    // ---------------------------- Printing ----------------------------
    std::string FlowTypeToString(FlowType type) {
      switch (type) {
        case FlowType::BASIC: return "BASIC";
        case FlowType::LOOP: return "LOOP";
        case FlowType::ROUTINE: return "ROUTINE";
        case FlowType::CALL: return "CALL";
        default: return "UNKNOWN";
      }
    }
    
    void PrintModules(std::ostream & os=std::cout) {
      os << "Modules: {";
      for (size_t i = 0; i < modules.size(); ++i) {
        if (i) os << ",";
        os << "[" << modules[i].begin << ":" << modules[i].end << "]";
        os << "(";
        modules[i].tag.Print(os);
        os << ")";
      }
      os << "}";
    }

    void PrintModuleSequences(std::ostream & os=std::cout) {
      for (size_t i = 0; i < modules.size(); ++i) {
        Module & module = modules[i];
        os << "Module["<<module.id<<"](";
        module.tag.Print(os);
        os << ")\n";
        size_t msize = (module.begin < module.end) ? module.end - module.begin : (program.GetSize() - module.begin) + module.end;
        std::cout << "#size=" << msize << std::endl;
        for (size_t mip = 0; mip < msize; ++mip) {
          const size_t ip = (module.begin+mip)%program.GetSize();
          os << "  ip[" << ip << "]: ";
          program.PrintInst(program[ip], os);
          os << "\n";
        }
      }
    }

    void PrintMemoryVerbose(const Memory & memory, const emp::vector<tag_t> & mem_tags, std::ostream & os=std::cout, const std::string & indent="") {
      emp_assert(memory.GetSize() == mem_tags.size());
      for (size_t i = 0; i < memory.GetSize(); ++i) {
        os << indent << "mem[" << i << "](";
        mem_tags[i].Print(os);
        os << "): ";
        memory.GetPos(i).Print(os);
        if (!memory.IsSet(i)) { os << " (unset)"; }
        os << "\n";
      } 
    }

    void PrintHardwareMemoryVerbose(std::ostream & os=std::cout) {
      os << "-- Global memory --\n";
      for (size_t gi = 0; gi < mem_size; ++gi) {
        os << "  mem[" << gi << "](";
        mem_tags[gi].Print(os);
        os << "): ";
        global_mem.GetPos(gi).Print(os);
        if (!global_mem.IsSet(gi)) { os << " (unset)"; }
        os << "\n";
      }
      for (int ci = (int)call_stack.size()-1; ci >= 0; --ci) {
        CallState & state = call_stack[ci];
        os << "Local Memory (stack id=" << ci << ")\n";
        
        os << "  -- Input Memory -- \n";
        for (size_t i = 0; i < mem_size; ++i) {
          os << "    mem[" << i << "](";
          mem_tags[i].Print(os);
          os << "): ";
          state.GetInputMem().GetPos(i).Print(os);
          if (!state.GetInputMem().IsSet(i)) { os << " (unset)"; }
          os << "\n";
        } 

        os << "  -- Working Memory -- \n";
        for (size_t i = 0; i < mem_size; ++i) {
          os << "    mem[" << i << "](";
          mem_tags[i].Print(os);
          os << "): ";
          state.GetWorkingMem().GetPos(i).Print(os);
          if (!state.GetWorkingMem().IsSet(i)) { os << " (unset)"; }
          os << "\n";
        } 

        os << "  -- Output Memory -- \n";
        for (size_t i = 0; i < mem_size; ++i) {
          os << "    mem[" << i << "](";
          mem_tags[i].Print(os);
          os << "): ";
          state.GetOutputMem().GetPos(i).Print(os);
          if (!state.GetOutputMem().IsSet(i)) { os << " (unset)"; }
          os << "\n";
        } 

      }
    }

    void PrintHardwareState(std::ostream & os=std::cout) {
      // Print state of global memory
      os << "Global memory (size=" << global_mem.GetSize() << "):\n";
      PrintMemoryVerbose(global_mem, mem_tags, os, "  ");
      // Print call stack
      os << "Call stack (size=" << call_stack.size() << "):\n";
      for (int ci = (int)call_stack.size()-1; ci >= 0; --ci) {
        CallState & state = call_stack[ci];
        os << "----------ID="<<ci<<"----------\n";
        os << "  Attributes: {" << "returnable:" << state.IsReturnable() << "," << "circular:" << state.IsCircular() << "}\n";
        os << "  Flow stack: [";
        for (int fi = (int)state.GetFlowStack().size()-1; fi >= 0; --fi) {
          Flow & flow = state.GetFlowStack()[fi];
          os << "{";
          os << "type:" << FlowTypeToString(flow.type) << ",";
          os << "mp:" << flow.mptr << ",";
          os << "ip:" << flow.iptr << ",";
          os << "begin:" << flow.begin << ",";
          os << "end:" << flow.end << ",";
          os << "iter:" << flow.iter; 
          os << "}";
          if (fi) os << ",";
        }
        os << "]\n";

        if (state.IsFlow()) {
          os << "  Inst: ";
          if (state.GetTopFlow().iptr < program.GetSize()) {
            program.PrintInst(program[state.GetTopFlow().iptr], os);
          } else {
            os << " INVALID";
          }
          os << "\n";
        }

        os << "  Input memory (size=" << state.GetInputMem().GetSize() << "):\n";
        PrintMemoryVerbose(state.GetInputMem(), mem_tags, os, "    ");
        os << "  Working memory (size=" << state.GetWorkingMem().GetSize() << "):\n";
        PrintMemoryVerbose(state.GetWorkingMem(), mem_tags, os, "    ");
        os << "  Output memory (size=" << state.GetOutputMem().GetSize() << "):\n";
        PrintMemoryVerbose(state.GetOutputMem(), mem_tags, os, "    ");
      }
      os << "----------BOTTOM----------\n"; 
    }

    void PrintHardwareStateFlat(std::ostream & os=std::cout) {
      // os << "Global memory tags: ";
      // os << "Input memory tags: ";
      // os << "Working memory tags: ";
      // os << "Output memory tags: ";   
      // Print state of global memory
      os << "Global memory (size=" << global_mem.GetSize() << "):";
      global_mem.Print(os);
      os << "\n";
      
      // Print call stack
      os << "Call stack (size=" << call_stack.size() << "):\n";
      for (int ci = (int)call_stack.size()-1; ci >= 0; --ci) {
        CallState & state = call_stack[ci];
        os << "----------ID="<<ci<<"----------\n";
        os << "  Attributes: {" << "returnable:" << state.IsReturnable() << "," << "circular:" << state.IsCircular() << "}\n";
        os << "  Flow stack: [";
        for (int fi = (int)state.GetFlowStack().size()-1; fi >= 0; --fi) {
          Flow & flow = state.GetFlowStack()[fi];
          os << "{";
          os << "type:" << FlowTypeToString(flow.type) << ",";
          os << "mp:" << flow.mptr << ",";
          os << "ip:" << flow.iptr << ",";
          os << "begin:" << flow.begin << ",";
          os << "end:" << flow.end;
          os << "}";
          if (fi) os << ",";
        }
        os << "]\n";

        os << "  Input memory (size=" << state.GetInputMem().GetSize() << "):";
        state.GetInputMem().Print(os);
        os << "\n";
        os << "  Working memory (size=" << state.GetWorkingMem().GetSize() << "):";
        state.GetWorkingMem().Print(os);
        os << "\n";
        os << "  Output memory (size=" << state.GetOutputMem().GetSize() << "):";
        state.GetWorkingMem().Print(os);
        os << "\n";
      }
      os << "----------BOTTOM----------\n"; 
    }

    // ------- Instructions --------

    // - Numeric -

    static void Inst_Add(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();
      // Find arguments.
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;
      //  mem[C] = mem[A] + mem[B]
      const double A = wmem.AccessVal(posA).GetNum();
      const double B = wmem.AccessVal(posB).GetNum();
      wmem.Set(posC, A + B);
    }

    static void Inst_Sub(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments.
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;

      //  mem[C] = mem[A] - mem[B]
      const double A = wmem.AccessVal(posA).GetNum();
      const double B = wmem.AccessVal(posB).GetNum();
      wmem.Set(posC, A - B);
    }

    static void Inst_Mult(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments.
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;

      //  mem[C] = mem[A] * mem[B]
      const double A = wmem.AccessVal(posA).GetNum();
      const double B = wmem.AccessVal(posB).GetNum();
      wmem.Set(posC, A * B);
    }

    static void Inst_Div(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments.
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;

      //  mem[C] = mem[A] / mem[B]
      const double A = wmem.AccessVal(posA).GetNum();
      const double B = wmem.AccessVal(posB).GetNum();
      if (B != 0.0) { wmem.Set(posC, A / B); }
    }

    static void Inst_Mod(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments.
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;

      //  mem[C] = mem[A] % mem[B]
      const int A = (int)wmem.AccessVal(posA).GetNum();
      const int B = (int)wmem.AccessVal(posB).GetNum();
      if (B != 0) { wmem.Set(posC, static_cast<int64_t>(A) % static_cast<int64_t>(B)); }
    }

    static void Inst_Inc(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posA)) return;

      // mem[A] += 1;
      const double A = wmem.AccessVal(posA).GetNum();
      wmem.Set(posA, A+1.0);
    }

    static void Inst_Dec(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posA)) return;

      // mem[A] -= 1;
      const double A = wmem.AccessVal(posA).GetNum();
      wmem.Set(posA, A-1.0);
    }

    static void Inst_Not(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posA)) return;

      // mem[A] = !mem[A];
      const double A = wmem.AccessVal(posA).GetNum();
      wmem.Set(posA, (double)(!(bool)A));
    }

    static void Inst_TestNumEqu(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments.
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;

      const double A = wmem.AccessVal(posA).GetNum();
      const double B = wmem.AccessVal(posB).GetNum();

      wmem.Set(posC, A == B);
    }

    static void Inst_TestNumNEqu(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments.
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;

      const double A = wmem.AccessVal(posA).GetNum();
      const double B = wmem.AccessVal(posB).GetNum();

      wmem.Set(posC, A != B);
    }

    static void Inst_TestNumLess(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments.
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;

      const double A = wmem.AccessVal(posA).GetNum();
      const double B = wmem.AccessVal(posB).GetNum();
      wmem.Set(posC, A < B);
    }

    static void Inst_TestNumLessTEqu(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments.
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;

      const double A = wmem.AccessVal(posA).GetNum();
      const double B = wmem.AccessVal(posB).GetNum();
      wmem.Set(posC, A <= B);
    }

    static void Inst_TestNumGreater(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments.
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;

      const double A = wmem.AccessVal(posA).GetNum();
      const double B = wmem.AccessVal(posB).GetNum();
      wmem.Set(posC, A > B);
    }

    static void Inst_TestNumGreaterTEqu(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments.
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;

      const double A = wmem.AccessVal(posA).GetNum();
      const double B = wmem.AccessVal(posB).GetNum();
      wmem.Set(posC, A >= B);
    }

    static void Inst_Floor(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posA)) return;

      // floor(mem[A]);
      const double A = wmem.AccessVal(posA).GetNum();
      wmem.Set(posA, std::floor(A));
    }

    // - Memory -
    static void Inst_CopyMem(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments.
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing

      wmem.Set(posB, wmem.GetPos(posA));
    }

    static void Inst_SwapMem(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments.
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing

      wmem.Swap(posA, posB);
    }

    static void Inst_Input(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();
      memory_t & imem = state.GetInputMem();

      // Find arguments.
      size_t iposA = hw.FindBestMemoryMatch(imem, inst.arg_tags[0], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(iposA)) return; // Do nothing
      size_t wposB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(wposB)) return; // Do nothing

      wmem.Set(wposB, imem.GetPos(iposA));
    }

    static void Inst_Output(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();
      memory_t & omem = state.GetOutputMem();

      // Find arguments.
      size_t wposA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(wposA)) return; // Do nothing

      size_t oposB = hw.FindBestMemoryMatch(omem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(oposB)) return; // Do nothing

      omem.Set(oposB, wmem.GetPos(wposA));
    }

    static void Inst_CommitGlobal(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();
      memory_t & gmem = hw.GetGlobalMem();

      // Find arguments.
      size_t wposA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(wposA)) return; // Do nothing
      size_t gposB = hw.FindBestMemoryMatch(gmem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(gposB)) return; // Do nothing

      gmem.Set(gposB, wmem.GetPos(wposA));
    }

    static void Inst_PullGlobal(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();
      memory_t & gmem = hw.GetGlobalMem();

      // Find arguments.
      size_t gposA = hw.FindBestMemoryMatch(gmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(gposA)) return; // Do nothing
      size_t wposB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(wposB)) return; // Do nothing

      wmem.Set(wposB, gmem.GetPos(gposA));
    }

    static void Inst_TestMemEqu(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments.
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;

      wmem.Set(posC, wmem.GetPos(posA) == wmem.GetPos(posB));
    }

    static void Inst_TestMemNEqu(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find arguments.
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;

      wmem.Set(posC, wmem.GetPos(posA) != wmem.GetPos(posB));
    }

    // - Vectors -
    /// Instruction: MakeVector (MemPos A, MemPos B, MemPos C)
    /// - mem[C] = {mem[A:B]} || {mem[B:A]}
    static void Inst_MakeVector(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;

      posA = (posA <= posB) ? posA : posB;
      emp_assert(posA <= posB);

      emp::vector<MemoryValue> vec;
      for (size_t i = posA; i <= posB; ++i) {
        if (wmem.IsVec(i)) continue;
        vec.emplace_back(wmem.AccessVal(i));
      }
      // We'll only set mem[posc] if there were values between [posA:posB]
      // if (vec.size()) { // If vec has members, assign mem[posc] = vec
      //   wmem.Set(posC, vec);
      // }
      wmem.Set(posC, vec);
    }

    /// Instruction: VecGet (VEC, NUM, POS)
    /// mem[C] = mem[A][mem[B]]
    static void Inst_VecGet(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::VEC);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;
      
      emp::vector<MemoryValue> & vec = wmem.AccessVec(posA);
      size_t i = (size_t)wmem.AccessVal(posB).GetNum();
      if (i < vec.size()) {
        wmem.Set(posC, vec[i]);
      }
    }

    /// Instruction: VecSet (VEC, NUM, POS)
    /// mem[A][mem[B]] = mem[C]
    static void Inst_VecSet(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::VEC);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], {MemPosType::NUM, MemPosType::STR}, hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;
      
      emp::vector<MemoryValue> & vec = wmem.AccessVec(posA);
      size_t i = (size_t)wmem.AccessVal(posB).GetNum();
      if (i < vec.size()) {
        vec[i] = wmem.AccessVal(posC);
      }
    }

    /// Instruction: VecLen (VEC, POS)
    static void Inst_VecLen(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();
      
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::VEC);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing

      wmem.Set(posB, wmem.AccessVec(posA).size());
    }

    /// Instruction: VecAppend (VEC, NUM/STR)
    static void Inst_VecAppend(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();
      
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::VEC);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], {MemPosType::NUM, MemPosType::STR}, hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing

      emp::vector<MemoryValue> & vec = wmem.AccessVec(posA);
      vec.emplace_back(wmem.AccessVal(posB));
    }

    /// Instruction VecPop (VEC, POS)
    static void Inst_VecPop(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();
      
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::VEC);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing

      emp::vector<MemoryValue> & vec = wmem.AccessVec(posA);
      if (vec.size()) {
        MemoryValue val(vec.back()); // Need to do things this way because vector could be overwritten with pop.
        vec.pop_back();
        wmem.Set(posB, val);
      }
    }

    /// Instruction: VecRemove (VEC, NUM)
    /// 
    static void Inst_VecRemove(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();
      
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::VEC);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posB)) return; // Do nothing

      emp::vector<MemoryValue> & vec = wmem.AccessVec(posA);
      size_t i = (size_t)wmem.AccessVal(posB).GetNum();
      if (i < vec.size()) {
        vec.erase(vec.begin()+i);
      }
    }

    static void Inst_VecReplaceAll(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::VEC);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], {MemPosType::NUM, MemPosType::STR}, hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], {MemPosType::NUM, MemPosType::STR}, hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;
      
      emp::vector<MemoryValue> & vec = wmem.AccessVec(posA);
      MemoryValue remove_val(wmem.AccessVal(posB));
      MemoryValue replace_val(wmem.AccessVal(posC));
      for (size_t i = 0; i < vec.size(); ++i) {
        if (vec[i] == remove_val) vec[i] = replace_val; 
      }
    }

    /// Instruction: VecIndexOf (VEC, VAL, POS)
    static void Inst_VecIndexOf(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::VEC);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], {MemPosType::NUM, MemPosType::STR}, hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;
      
      emp::vector<MemoryValue> & vec = wmem.AccessVec(posA);
      MemoryValue & val = wmem.AccessVal(posB);
      bool found = false;
      for (size_t i = 0; i < vec.size(); ++i) {
        if (vec[i] == val) {
          found = true;
          wmem.Set(posC, i);
          break;
        }
      }
      if (!found) wmem.Set(posC, vec.size());
    }

    static void Inst_VecOccurrencesOf(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::VEC);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], {MemPosType::NUM, MemPosType::STR}, hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;
      
      emp::vector<MemoryValue> & vec = wmem.AccessVec(posA);
      MemoryValue & val = wmem.AccessVal(posB);
      size_t cnt = 0;
      for (size_t i = 0; i < vec.size(); ++i) {
        if (vec[i] == val) {
          cnt += 1;
        }
      }
      wmem.Set(posC, cnt);
    }

    static void Inst_VecReverse(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::VEC);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      
      emp::vector<MemoryValue> & vec = wmem.AccessVec(posA);
      std::reverse(std::begin(vec), std::end(vec));
    }

    static void Inst_VecSwapIfLess(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::VEC);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity(), MemPosType::NUM);
      if (!hw.IsValidMemPos(posC)) return;
      
      emp::vector<MemoryValue> & vec = wmem.AccessVec(posA);
      
      const size_t ai = (size_t)wmem.AccessVal(posB).GetNum();
      const size_t bi = (size_t)wmem.AccessVal(posC).GetNum();

      if (ai < vec.size() && bi < vec.size()) {
        // Want numbers
        double a = (vec[ai].GetType() == MemoryValue::MemoryType::NUM) ? vec[ai].GetNum() : 0;
        double b = (vec[bi].GetType() == MemoryValue::MemoryType::NUM) ? vec[bi].GetNum() : 0;
        if (a < b) {
          std::swap(vec[ai], vec[bi]);
        } 
        // else {
        // }
      }
    }

    static void Inst_VecGetFront(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::VEC);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing

      emp::vector<MemoryValue> & vec = wmem.AccessVec(posA);
      if (vec.size()) {
        wmem.Set(posB, vec.front());
      }
    }

    static void Inst_VecGetBack(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::VEC);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing

      emp::vector<MemoryValue> & vec = wmem.AccessVec(posA);
      if (vec.size()) {
        wmem.Set(posB, vec.back());
      }
    }

    // - str instructions -
    static void Inst_StrLength(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();
      
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::STR);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing

      wmem.Set(posB, wmem.AccessVal(posA).GetStr().size());
    } // todo - test

    static void Inst_StrConcat(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();
      // Find arguments.
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::STR);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity(), MemPosType::STR);
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      size_t posC = hw.FindBestMemoryMatch(wmem, inst.arg_tags[2], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posC)) return;
      //  mem[C] = mem[A] + mem[B]
      const std::string & A = wmem.AccessVal(posA).GetStr();
      const std::string & B = wmem.AccessVal(posB).GetStr();
      if ((A.size() + B.size()) <= DEFAULT_MAX_STR_LEN) { // Don't want exponential string growth...
        wmem.Set(posC, A + B);
      }
    } // todo - test

    static void Inst_StrToVec(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();
      // Find arguments.
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), MemPosType::STR);
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing
      //  mem[C] = mem[A] + mem[B]
      const std::string & A = wmem.AccessVal(posA).GetStr();
      emp::vector<std::string> vec;
      for (size_t i = 0; i < A.size(); ++i) {
        vec.emplace_back(A[i]);
      }
      wmem.Set(posB, vec);
    } // todo - test

    // - types -
    static void Inst_IsStr(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing

      wmem.Set(posB, (double)(wmem.GetPosType(posA) == MemPosType::STR));
    }

    static void Inst_IsNum(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing

      wmem.Set(posB, (double)(wmem.GetPosType(posA) == MemPosType::NUM));
    }

    static void Inst_IsVec(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posA)) return; // Do nothing
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posB)) return; // Do nothing

      wmem.Set(posB, (double)(wmem.GetPosType(posA) == MemPosType::VEC));
    }

    // static void Inst_ToStr(hardware_t & hw, const inst_t & inst) {

    // }

    // static void Inst_ToNum(hardware_t & hw, const inst_t & inst) {

    // }

    // - Flow control -
    /// Instruction: If ()
    /// If (mem[A]) { }
    static void Inst_If(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find the end of flow
      const size_t cur_ip = state.GetIP();
      const size_t cur_mp = state.GetMP();
      size_t eof = hw.FindEndOfFlow(cur_mp, cur_ip);
      size_t bof = (cur_ip == 0) ? hw.GetProgram().GetSize() - 1 : cur_ip - 1;
      
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());

      bool skip = false;
      if (!hw.IsValidMemPos(posA)) skip = true;                       // SKip if failed to find valid mem pos
      else if (wmem.GetPosType(posA) != MemPosType::NUM) skip = true;    // Skip if best match is not a number
      else if (wmem.AccessVal(posA).GetNum() == 0) skip = true;       // Skip if best match = 0 

      if (skip) {
        state.SetIP(eof);
        // Advance past the flow close if not at end of module
        if (hw.ValidPosition(state.GetMP(), eof)) state.AdvanceIP();
      } else {
        // Open flow
        hw.OpenFlow(state, FlowType::BASIC, bof, eof, cur_mp, cur_ip);
      }
    }

    static void Inst_IfNot(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find the end of flow
      const size_t cur_ip = state.GetIP();
      const size_t cur_mp = state.GetMP();
      size_t eof = hw.FindEndOfFlow(cur_mp, cur_ip);
      size_t bof = (cur_ip == 0) ? hw.GetProgram().GetSize() - 1 : cur_ip - 1;
      
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());

      bool skip = false;
      if (!hw.IsValidMemPos(posA)) skip = true;                       // SKip if failed to find valid mem pos
      else if (wmem.GetPosType(posA) != MemPosType::NUM) skip = true; // Skip if best match is not a number
      else if (wmem.AccessVal(posA).GetNum() != 0) skip = true;       // Skip if best match != 0 

      if (skip) {
        state.SetIP(eof);
        // Advance past the flow close if not at end of module
        if (hw.ValidPosition(state.GetMP(), eof)) state.AdvanceIP();
      } else {
        // Open flow
        hw.OpenFlow(state, FlowType::BASIC, bof, eof, cur_mp, cur_ip);
      }
    }

    static void Inst_While(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find the end of flow
      const size_t cur_ip = state.GetIP();
      const size_t cur_mp = state.GetMP();
      size_t eof = hw.FindEndOfFlow(cur_mp, cur_ip);
      size_t bof = (cur_ip == 0) ? hw.GetProgram().GetSize() - 1 : cur_ip - 1;
      
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());

      bool skip = false;
      if (!hw.IsValidMemPos(posA)) skip = true;                       // SKip if failed to find valid mem pos
      else if (wmem.GetPosType(posA) != MemPosType::NUM) skip = true;    // Skip if best match is not a number
      else if (wmem.AccessVal(posA).GetNum() == 0) skip = true;       // Skip if best match = 0 

      if (skip) {
        state.SetIP(eof);
        state.GetTopFlow().iter = 0;
        // Advance past the flow close if not at end of module
        if (hw.ValidPosition(state.GetMP(), eof)) state.AdvanceIP();
      } else {
        // Open flow
        hw.OpenFlow(state, FlowType::LOOP, bof, eof, cur_mp, cur_ip);
      }
    }

    static void Inst_Countdown(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find the end of flow
      const size_t cur_ip = state.GetIP();
      const size_t cur_mp = state.GetMP();
      size_t eof = hw.FindEndOfFlow(cur_mp, cur_ip);
      size_t bof = (cur_ip == 0) ? hw.GetProgram().GetSize() - 1 : cur_ip - 1;
      
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());

      bool skip = false;
      if (!hw.IsValidMemPos(posA)) skip = true;                       // SKip if failed to find valid mem pos
      else if (wmem.GetPosType(posA) != MemPosType::NUM) skip = true;    // Skip if best match is not a number
      else if (wmem.AccessVal(posA).GetNum() == 0) skip = true;       // Skip if best match = 0 

      if (skip) {
        state.SetIP(eof);
        state.GetTopFlow().iter = 0;
        // Advance past the flow close if not at end of module
        if (hw.ValidPosition(state.GetMP(), eof)) state.AdvanceIP();
      } else {
        emp_assert(wmem.GetPosType(posA) == MemPosType::NUM);
        // Decrement
        const double val = wmem.AccessVal(posA).GetNum();
        wmem.Set(posA, val-1);
        // Open flow
        hw.OpenFlow(state, FlowType::LOOP, bof, eof, cur_mp, cur_ip);
      }
    }


    /// For each in vec
    /// Foreach (POS, VEC)
    static void Inst_Foreach(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      memory_t & wmem = state.GetWorkingMem();

      // Find the end of flow
      const size_t cur_ip = state.GetIP();
      const size_t cur_mp = state.GetMP();
      size_t eof = hw.FindEndOfFlow(cur_mp, cur_ip);
      size_t bof = (cur_ip == 0) ? hw.GetProgram().GetSize() - 1 : cur_ip - 1;
      
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
      size_t posB = hw.FindBestMemoryMatch(wmem, inst.arg_tags[1], hw.GetMinTagSpecificity(), MemPosType::VEC);

      bool skip = false;
      if (!hw.IsValidMemPos(posA) || !hw.IsValidMemPos(posB)) skip = true;                       // SKip if failed to find valid mem pos
      else if (state.GetTopFlow().iter >= wmem.AccessVec(posB).size()) skip = true;
  
      if (skip) {
        state.SetIP(eof);
        state.GetTopFlow().iter = 0; // Loop's parent flow
        // Advance past the flow close if not at end of module
        if (hw.ValidPosition(state.GetMP(), eof)) state.AdvanceIP();
      } else {
        // emp_assert(wmem.GetPosType(posA) == MemPosType::NUM);
        // mem[posA] = vec[iter]
        wmem.Set(posA, wmem.AccessVec(posB)[state.GetTopFlow().iter]);
        // Open flow
        hw.OpenFlow(state, FlowType::LOOP, bof, eof, cur_mp, cur_ip);
      }
    }

    static void Inst_Close(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      FlowType cur_flow_type = state.GetTopFlow().type;
      if (cur_flow_type == FlowType::BASIC || cur_flow_type == FlowType::LOOP) {
        hw.CloseFlow(state);
      }
    }

    static void Inst_Break(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      FlowType cur_flow_type = state.GetTopFlow().type;
      if (cur_flow_type == FlowType::BASIC || cur_flow_type == FlowType::LOOP) {
        hw.BreakFlow(state);
      }
    }

    static void Inst_Call(hardware_t & hw, const inst_t & inst) {
      hw.CallModule(inst.arg_tags[0], hw.GetMinTagSpecificity());
    }

    static void Inst_Routine(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      const size_t target_module = hw.FindBestModuleMatch(inst.arg_tags[0], hw.GetMinTagSpecificity());
      if (target_module < hw.GetModuleCnt()) {
        hw.OpenFlowRoutine(state, target_module);
      }
    }

    static void Inst_Return(hardware_t & hw, const inst_t & inst) {
      CallState & state = hw.GetCurCallState();
      while (state.IsFlow()) {
        Flow & top = state.GetTopFlow();
        if (top.type == FlowType::CALL || top.type == FlowType::ROUTINE) {
          hw.CloseFlow(state);
          break;
        } else {
          hw.CloseFlow(state);
        }
      }
    }

    
    // - MISC -
    static void Inst_Nop(hardware_t & hw, const inst_t & inst) { return; }
    

  };

  // TODO - make writing loop instructions way easier
  // todo - make writing instructions in general way easier

}



#endif