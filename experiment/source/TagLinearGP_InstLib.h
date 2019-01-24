#ifndef TAG_LINEAR_GP_INSTLIB_H
#define TAG_LINEAR_GP_INSTLIB_H

#include <map>
#include <string>
#include <unordered_set>

#include "base/array.h"
#include "base/Ptr.h"
#include "base/vector.h"
#include "tools/map_utils.h"
#include "tools/string_utils.h"

namespace TagLGP {

  /// Instruction library class used by TagLinearGP class.
  /// Original version pulled from InstLib.h in Empirical library.
  /// - Modified to potentially facilitate some experimental functionality.
  template<typename HARDWARE_T, size_t MAX_ARG_CNT=3>
  class InstLib {
  public:
    using hardware_t = HARDWARE_T;
    using inst_t = typename hardware_t::inst_t;
    using fun_t = std::function<void(hardware_t &, const inst_t &)>; // Provide arguments, too?
    
    enum InstProperty { BEGIN_FLOW, END_FLOW, MODULE };

    using inst_properties_t = std::unordered_set<InstProperty>;

    struct InstDef {
      std::string name;                   ///< Name of this instruction.
      fun_t fun_call;                     ///< Function to call when executing.
      size_t num_args;                    ///< Number of args needed by function.
      std::string desc;                   ///< Description of function.
      inst_properties_t properties;

      InstDef(const std::string & _n, fun_t _fun, 
              size_t _nargs, const std::string & _d,
              const inst_properties_t & _properties = inst_properties_t())
        : name(_n), fun_call(_fun), num_args(_nargs), desc(_d), properties(_properties) { ; }
      InstDef(const InstDef &) = default;
    };
    
  protected:
    emp::vector<InstDef> inst_lib;           ///< Full definitions for instructions.
    emp::vector<fun_t> inst_funs;            ///< Map of instruction IDs to their functions.
    std::map<std::string, size_t> name_map;  ///< How do names link to instructions?

    
  public:
    InstLib() : inst_lib(), inst_funs(), name_map() { ; }  ///< Default Constructor
    InstLib(const InstLib &) = default;                    ///< Copy Constructor
    InstLib(InstLib &&) = default;                         ///< Move Constructor
    ~InstLib() { ; }

    InstLib & operator=(const InstLib &) = default;        ///< Copy Operator
    InstLib & operator=(InstLib &&) = default;       

    /// Return the name associated with the specified instruction ID.
    const std::string & GetName(size_t id) const { return inst_lib[id].name; }

    /// Return the function associated with the specified instruction ID.
    const fun_t & GetFunction(size_t id) const { return inst_lib[id].fun_call; }

    /// Return the number of arguments expected for the specified instruction ID.
    size_t GetNumArgs(size_t id) const { return inst_lib[id].num_args; }

    /// Return the provided description for the provided instruction ID.
    const std::string & GetDesc(size_t id) const { return inst_lib[id].desc; }

    /// Return the set of properties for the provided instruction ID.
    const inst_properties_t & GetProperties(size_t id) const { return inst_lib[id].properties; }

    /// Does the given instruction ID have the given property value?
    bool HasProperty(size_t id, InstProperty property) const { return inst_lib[id].properties.count(property); }

    /// 
    std::string GetPropertyStr(InstProperty property) {
      switch (property) {
        case InstProperty::MODULE: return "MODULE";
        default: return "UNKNOWN";
      }
    }

    /// Get the number of instructions in this set.
    size_t GetSize() const { return inst_lib.size(); }

    bool IsInst(const std::string & name) const {
        return emp::Has(name_map, name);
    }

    /// Return the ID of the instruction that has the specified name.
    size_t GetID(const std::string & name) const {
      emp_assert(emp::Has(name_map, name), name);
      return emp::Find(name_map, name, (size_t) -1);
    }

    /// @brief Add a new instruction to the set.
    /// @param name A unique string name for this instruction.
    /// @param fun_call The function that should be called when this instruction is executed.
    /// @param num_args How many arguments does this function require? (default=0)
    /// @param desc A description of how this function operates. (default="")
    void AddInst(const std::string & name,
                 const fun_t & fun_call,
                 size_t num_args=0,
                 const std::string & desc="",
                 const inst_properties_t & inst_properties=inst_properties_t())
    {
      const size_t id = inst_lib.size();
      inst_lib.emplace_back(name, fun_call, num_args, desc, inst_properties);
      inst_funs.emplace_back(fun_call);
      name_map[name] = id;
    }

    /// Process a specified instruction in the provided hardware.
    void ProcessInst(hardware_t & hw, const inst_t & inst) const {
      inst_funs[inst.id](hw, inst);
    }

    /// Process a specified instruction on hardware that can be converted to the correct type.
    template <typename IN_HW>
    void ProcessInst(emp::Ptr<IN_HW> hw, const inst_t & inst) const {
      emp_assert( dynamic_cast<hardware_t*>(hw.Raw()) );
      inst_funs[inst.id](*(hw.template Cast<hardware_t>()), inst);
    }


    void Print(std::ostream & os=std::cout) {
      os << "{\n";
      for (size_t i = 0; i < inst_lib.size(); ++i) {
        os << "inst: {\n";
        os << "  name: " <<  inst_lib[i].name  << "\n"; 
        os << "  num args: " << inst_lib[i].num_args << "\n";
        os << "  desc: \"" << inst_lib[i].desc << "\"\n";
        os << "  properties: [";
        size_t cnt = 0;
        for (const auto & property : inst_lib[i].properties) {
          if (cnt) os << ",";
          os << property;
          ++cnt;
        } 
        os << "]\n";
      }
      os << "}";
    }
  };

}

#endif