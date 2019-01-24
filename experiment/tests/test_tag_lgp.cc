#define CATCH_CONFIG_MAIN
#include "third-party/Catch/single_include/catch.hpp"

#include <iostream>
#include <cmath>

#include "base/Ptr.h"
#include "base/vector.h"
#include "tools/Random.h"

#include "TagLinearGP.h"
#include "TagLinearGP_InstLib.h"
#include "TagLinearGP_Utilities.h"
#include "Utilities.h"

// String => number; number=>string conversion
//  - If str.numeric: return double(str)
//    else: return 0.0
// Number to string
//  return to_string(number)

TEST_CASE("Inst_Add", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 2;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("Add", hardware_t::Inst_Add, 3, "mem[C] = mem[A] + mem[B]");
  // RND + RND = (RND + RND)
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());

    prog.PushInst("Add", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();
    const double A = (double)random->GetInt(-1000, 1000);
    const double B = (double)random->GetInt(-1000, 1000);
    const double C = (posA != posB) ? A + B : B + B;
    wmem.Set(posA, A);
    wmem.Set(posB, B);
    cpu.SingleProcess();
    double memC = wmem.AccessVal(posC).GetNum();
    if (memC != C) {
      cpu.PrintHardwareState();
      std::cout << "posA = " << posA << "; posB = " << posB << "; posC = " << posC << std::endl;
      std::cout << "A = " << A << "; B = " << B << "; C = " << C << std::endl;
    }
    REQUIRE(memC == C);
  }
  /////////////////////////////////////

  // Clean up
  inst_lib.Delete();
  random.Delete();
}


TEST_CASE("Inst_Sub", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 3;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("Sub", hardware_t::Inst_Sub, 3, "mem[C] = mem[A] - mem[B]");
  // RND + RND = (RND + RND)
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());

    prog.PushInst("Sub", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();
    const double A = (double)random->GetInt(-1000, 1000);
    const double B = (double)random->GetInt(-1000, 1000);
    const double C = (posA != posB) ? A - B : B - B;
    wmem.Set(posA, A);
    wmem.Set(posB, B);
    cpu.SingleProcess();
    double memC = wmem.AccessVal(posC).GetNum();
    if (memC != C) {
      cpu.GetProgram().Print();
      cpu.PrintHardwareState();
      std::cout << "posA = " << posA << "; posB = " << posB << "; posC = " << posC << std::endl;
      std::cout << "A = " << A << "; B = " << B << "; C = " << C << std::endl;
    }
    REQUIRE(memC == C);
  }
}


TEST_CASE("Inst_Mult", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 4;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("Mult", hardware_t::Inst_Mult, 3, "mem[C] = mem[A] * mem[B]");
  // RND + RND = (RND + RND)
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());

    prog.PushInst("Mult", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();
    const double A = (double)random->GetDouble(-1000, 1000);
    const double B = (double)random->GetDouble(-1000, 1000);
    const double C = (posA != posB) ? A * B : B * B;
    wmem.Set(posA, A);
    wmem.Set(posB, B);
    cpu.SingleProcess();
    double memC = wmem.AccessVal(posC).GetNum();
    if (memC != C) {
      cpu.GetProgram().Print();
      cpu.PrintHardwareState();
      std::cout << "posA = " << posA << "; posB = " << posB << "; posC = " << posC << std::endl;
      std::cout << "A = " << A << "; B = " << B << "; C = " << C << std::endl;
    }
    REQUIRE(memC == C);
  }
}


TEST_CASE("Inst_Inc", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 4;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("Inc", hardware_t::Inst_Inc, 1, "++mem-NUM[A]");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());

    prog.PushInst("Inc", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();
    const double A = (double)random->GetDouble(-1000, 1000);
    wmem.Set(posA, A);
    cpu.SingleProcess();
    double memA = wmem.AccessVal(posA).GetNum();
    REQUIRE(memA == A+1);
  }
}


TEST_CASE("Inst_Dec", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 5;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("Dec", hardware_t::Inst_Dec, 1, "--mem-NUM[A]");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());

    prog.PushInst("Dec", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();
    const double A = (double)random->GetDouble(-1000, 1000);
    wmem.Set(posA, A);
    cpu.SingleProcess();
    double memA = wmem.AccessVal(posA).GetNum();
    REQUIRE(memA == A-1);
  }
}


TEST_CASE("Inst_Div", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 6;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("Div", hardware_t::Inst_Div, 3, "mem-ANY[C] = mem-NUM[A] / mem-NUM[B]");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());

    prog.PushInst("Div", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();
    const double A = (double)random->GetDouble(-1000, 1000);
    const double B = (double)random->GetDouble(-1000, 1000);
    double C = 0;

    wmem.Set(posA, A);
    wmem.Set(posB, B);
    
    if (B != 0.0) C = (posA != posB) ? A / B : B / B;
    else C = wmem.AccessVal(posC).GetNum();

    cpu.SingleProcess();
    double memC = wmem.AccessVal(posC).GetNum();
    if (memC != C) {
      cpu.GetProgram().Print();
      cpu.PrintHardwareState();
      std::cout << "posA = " << posA << "; posB = " << posB << "; posC = " << posC << std::endl;
      std::cout << "A = " << A << "; B = " << B << "; C = " << C << std::endl;
    }
    REQUIRE(memC == C);
  }
}


TEST_CASE("Inst_Mod", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 7;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("Mod", hardware_t::Inst_Mod, 3, "mem-ANY[C] = mem-NUM[A] % mem-NUM[B]");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());

    prog.PushInst("Mod", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();
    const double A = (double)random->GetDouble(-1000, 1000);
    const double B = (double)random->GetDouble(-1000, 1000);
    double C = 0;

    wmem.Set(posA, A);
    wmem.Set(posB, B);
    
    if ((int)B != 0.0) C = (posA != posB) ? (int)A % (int)B : (int)B % (int)B;
    else C = wmem.AccessVal(posC).GetNum();

    cpu.SingleProcess();
    double memC = wmem.AccessVal(posC).GetNum();
    if (memC != C) {
      cpu.GetProgram().Print();
      cpu.PrintHardwareState();
      std::cout << "posA = " << posA << "; posB = " << posB << "; posC = " << posC << std::endl;
      std::cout << "A = " << A << "; B = " << B << "; C = " << C << std::endl;
    }
    REQUIRE(memC == C);
  }
}


TEST_CASE("Inst_TestNumEqu", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 8;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("TestNumEqu", hardware_t::Inst_TestNumEqu, 3, "mem[C] = mem[A] == mem[B]");
  // RND + RND = (RND + RND)
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());

    prog.PushInst("TestNumEqu", {matrix[posA], matrix[posB], matrix[posC]});
    prog.PushInst("TestNumEqu", {matrix[posC], matrix[posC], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();
    const double A = (double)random->GetDouble(-1000, 1000);
    const double B = (double)random->GetDouble(-1000, 1000);
    const double C = (posA != posB) ? A == B : B == B;

    wmem.Set(posA, A);
    wmem.Set(posB, B);

    cpu.SingleProcess();
    double memC = wmem.AccessVal(posC).GetNum();
    if (memC != C) {
      cpu.GetProgram().Print();
      cpu.PrintHardwareState();
      std::cout << "posA = " << posA << "; posB = " << posB << "; posC = " << posC << std::endl;
      std::cout << "A = " << A << "; B = " << B << "; C = " << C << std::endl;
    }
    REQUIRE(memC == C);

    cpu.SingleProcess();
    REQUIRE(wmem.AccessVal(posC).GetNum() == 1.0);
  }
}


TEST_CASE("Inst_TestNumNEqu", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 8;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("TestNumNEqu", hardware_t::Inst_TestNumNEqu, 3, "mem[C] = mem[A] != mem[B]");
  // RND + RND = (RND + RND)
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());

    prog.PushInst("TestNumNEqu", {matrix[posA], matrix[posB], matrix[posC]});
    prog.PushInst("TestNumNEqu", {matrix[posC], matrix[posC], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();
    const double A = (double)random->GetDouble(-1000, 1000);
    const double B = (double)random->GetDouble(-1000, 1000);
    const double C = (posA != posB) ? A != B : B != B;

    wmem.Set(posA, A);
    wmem.Set(posB, B);

    cpu.SingleProcess();
    double memC = wmem.AccessVal(posC).GetNum();
    if (memC != C) {
      cpu.GetProgram().Print();
      cpu.PrintHardwareState();
      std::cout << "posA = " << posA << "; posB = " << posB << "; posC = " << posC << std::endl;
      std::cout << "A = " << A << "; B = " << B << "; C = " << C << std::endl;
    }
    REQUIRE(memC == C);

    cpu.SingleProcess();
    REQUIRE(wmem.AccessVal(posC).GetNum() == 0.0);
  }
}


TEST_CASE("Inst_TestNumLess", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 8;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("TestNumLess", hardware_t::Inst_TestNumLess, 3, "mem[C] = mem[A] < mem[B]");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());

    prog.PushInst("TestNumLess", {matrix[posA], matrix[posB], matrix[posC]});
    prog.PushInst("TestNumLess", {matrix[posC], matrix[posC], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();
    const double A = (double)random->GetDouble(-1000, 1000);
    const double B = (double)random->GetDouble(-1000, 1000);
    const double C = (posA != posB) ? A < B : B < B;

    wmem.Set(posA, A);
    wmem.Set(posB, B);

    cpu.SingleProcess();
    double memC = wmem.AccessVal(posC).GetNum();
    if (memC != C) {
      cpu.GetProgram().Print();
      cpu.PrintHardwareState();
      std::cout << "posA = " << posA << "; posB = " << posB << "; posC = " << posC << std::endl;
      std::cout << "A = " << A << "; B = " << B << "; C = " << C << std::endl;
    }
    REQUIRE(memC == C);

    cpu.SingleProcess();
    REQUIRE(wmem.AccessVal(posC).GetNum() == 0.0);
  }
}


TEST_CASE("Inst_Floor", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 10;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("Floor", hardware_t::Inst_Floor, 1, "mem-NUM[A] = FLOOR(mem-NUM[A])");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());

    prog.PushInst("Floor", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();
    const double A = (double)random->GetDouble(-1000, 1000);
    wmem.Set(posA, A);
    cpu.SingleProcess();
    double memA = wmem.AccessVal(posA).GetNum();
    if (memA != std::floor(A)) {
      cpu.PrintHardwareState();
      std::cout << "posA = " << posA << std::endl;
      std::cout << "A = " << A << std::endl;
    }
    REQUIRE(memA == std::floor(A));
  }
}


TEST_CASE("Inst_CopyMem", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 11;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("CopyMem", hardware_t::Inst_CopyMem, 3, "mem-ANY[B] = mem-ANY[A]");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());
    prog.PushInst("CopyMem", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    size_t sA = random->GetUInt(0, 3);
    size_t sB = random->GetUInt(0, 3);
    size_t sC = random->GetUInt(0, 3);
    size_t sD = random->GetUInt(0, 3);

    double baseA = random->GetDouble(-1000, 1000);
    double baseB = random->GetDouble(-1000, 1000);
    double baseC = random->GetDouble(-1000, 1000);
    double baseD = random->GetDouble(-1000, 1000);

    // Set wmem to base values
    (sA == 1) ? wmem.Set(0, emp::to_string(baseA)) : wmem.Set(0, baseA);
    (sB == 1) ? wmem.Set(1, emp::to_string(baseB)) : wmem.Set(1, baseB);
    (sC == 1) ? wmem.Set(2, emp::to_string(baseC)) : wmem.Set(2, baseC);
    (sD == 1) ? wmem.Set(3, emp::to_string(baseD)) : wmem.Set(3, baseD);
    // drop some vectors in
    if (sA == 2) wmem.Set(0, emp::vector<std::string>{"hello", "world"});
    if (sB == 2) wmem.Set(1, emp::vector<mem_val_t>{wmem.AccessVal(0)});
    if (sC == 2) wmem.Set(2, emp::vector<mem_val_t>{wmem.AccessVal(1), wmem.AccessVal(2), wmem.AccessVal(3)});
    if (sD == 2) wmem.Set(3, emp::vector<std::string>{"blah", "blah", "blah"});

    cpu.SingleProcess();
    REQUIRE(wmem.GetPos(posA) == wmem.GetPos(posB));
  }
}


TEST_CASE("Inst_SwapMem", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 12;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("SwapMem", hardware_t::Inst_SwapMem, 3, "SWAP(mem-ANY[A], mem-ANY[B])");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());
    prog.PushInst("SwapMem", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    size_t sA = random->GetUInt(0, 3);
    size_t sB = random->GetUInt(0, 3);
    size_t sC = random->GetUInt(0, 3);
    size_t sD = random->GetUInt(0, 3);

    double baseA = random->GetDouble(-1000, 1000);
    double baseB = random->GetDouble(-1000, 1000);
    double baseC = random->GetDouble(-1000, 1000);
    double baseD = random->GetDouble(-1000, 1000);

    // Set wmem to base values
    (sA == 1) ? wmem.Set(0, emp::to_string(baseA)) : wmem.Set(0, baseA);
    (sB == 1) ? wmem.Set(1, emp::to_string(baseB)) : wmem.Set(1, baseB);
    (sC == 1) ? wmem.Set(2, emp::to_string(baseC)) : wmem.Set(2, baseC);
    (sD == 1) ? wmem.Set(3, emp::to_string(baseD)) : wmem.Set(3, baseD);
    // drop some vectors in
    if (sA == 2) wmem.Set(0, emp::vector<std::string>{"hello", "world"});
    if (sB == 2) wmem.Set(1, emp::vector<mem_val_t>{wmem.AccessVal(0)});
    if (sC == 2) wmem.Set(2, emp::vector<mem_val_t>{wmem.AccessVal(1), wmem.AccessVal(2), wmem.AccessVal(3)});
    if (sD == 2) wmem.Set(3, emp::vector<std::string>{"blah", "blah", "blah"});
    
    auto A(wmem.GetPos(posA));
    auto B(wmem.GetPos(posB));

    cpu.SingleProcess();

    REQUIRE(wmem.GetPos(posA) == B);
    REQUIRE(wmem.GetPos(posB) == A);
  }
}



TEST_CASE("Inst_Input", "[taglgp]") {
 // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 14;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("Input", hardware_t::Inst_Input, 2, "wmem-ANY[B] = imem-ANY[A]");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());
    prog.PushInst("Input", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();
    memory_t & imem = state.GetInputMem();

    size_t sA = random->GetUInt(0, 3);
    size_t sB = random->GetUInt(0, 3);
    
    double baseA = random->GetDouble(-1000, 1000);
    double baseB = random->GetDouble(-1000, 1000);
    
    // Set wmem to base values
    (sA == 1) ? imem.Set(posA, emp::to_string(baseA)) : imem.Set(posA, baseA);
    (sB == 1) ? wmem.Set(posB, emp::to_string(baseB)) : wmem.Set(posB, baseB);

    // drop some vectors in
    if (sA == 2) imem.Set(posA, emp::vector<std::string>{"hello", "world"});
    if (sB == 2) wmem.Set(posB, emp::vector<double>{random->GetDouble(-100, 100), random->GetDouble(-100, 100), random->GetDouble(-100, 100), random->GetDouble(-100, 100), random->GetDouble(-100, 100)});

    cpu.SingleProcess();

    REQUIRE(imem.GetPos(posA) == wmem.GetPos(posB));
  }
}



TEST_CASE("Inst_Output", "[taglgp]") {
 // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 14;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("Output", hardware_t::Inst_Output, 2, "omem-ANY[B] = wmem-ANY[A]");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());
    prog.PushInst("Output", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & omem = state.GetOutputMem();
    memory_t & wmem = state.GetWorkingMem();

    size_t sA = random->GetUInt(0, 3);
    size_t sB = random->GetUInt(0, 3);
    
    double baseA = random->GetDouble(-1000, 1000);
    double baseB = random->GetDouble(-1000, 1000);
    
    // Set wmem to base values
    (sA == 1) ? wmem.Set(posA, emp::to_string(baseA)) : wmem.Set(posA, baseA);
    (sB == 1) ? omem.Set(posB, emp::to_string(baseB)) : omem.Set(posB, baseB);

    // drop some vectors in
    if (sA == 2) wmem.Set(posA, emp::vector<std::string>{"hello", "world"});
    if (sB == 2) omem.Set(posB, emp::vector<double>{random->GetDouble(-100, 100), random->GetDouble(-100, 100), random->GetDouble(-100, 100), random->GetDouble(-100, 100), random->GetDouble(-100, 100)});

    cpu.SingleProcess();

    if (wmem.GetPos(posA) != omem.GetPos(posB)) {
      cpu.GetProgram().Print();
      cpu.PrintHardwareState();
      break;
    }
    
    REQUIRE(wmem.GetPos(posA) == omem.GetPos(posB));
  }
}



TEST_CASE("Inst_CommitGlobal", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 15;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("CommitGlobal", hardware_t::Inst_CommitGlobal, 2, "gmem-ANY[B] = wmem-ANY[A]");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());
    prog.PushInst("CommitGlobal", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & gmem = cpu.GetGlobalMem();
    memory_t & wmem = state.GetWorkingMem();

    size_t sA = random->GetUInt(0, 3);
    size_t sB = random->GetUInt(0, 3);
    
    double baseA = random->GetDouble(-1000, 1000);
    double baseB = random->GetDouble(-1000, 1000);
    
    // Set wmem to base values
    (sA == 1) ? wmem.Set(posA, emp::to_string(baseA)) : wmem.Set(posA, baseA);
    (sB == 1) ? gmem.Set(posB, emp::to_string(baseB)) : gmem.Set(posB, baseB);

    // drop some vectors in
    if (sA == 2) wmem.Set(posA, emp::vector<std::string>{"hello", "world"});
    if (sB == 2) gmem.Set(posB, emp::vector<double>{random->GetDouble(-100, 100), random->GetDouble(-100, 100), random->GetDouble(-100, 100), random->GetDouble(-100, 100), random->GetDouble(-100, 100)});

    cpu.SingleProcess();

    if (wmem.GetPos(posA) != gmem.GetPos(posB)) {
      cpu.GetProgram().Print();
      cpu.PrintHardwareState();
      break;
    }
    
    REQUIRE(wmem.GetPos(posA) == gmem.GetPos(posB));
  }
}



TEST_CASE("Inst_PullGlobal", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 21;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("PullGlobal", hardware_t::Inst_PullGlobal, 2, "wmem-ANY[B] = gmem-ANY[A]");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());
    prog.PushInst("PullGlobal", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();
    memory_t & gmem = cpu.GetGlobalMem();

    size_t sA = random->GetUInt(0, 3);
    size_t sB = random->GetUInt(0, 3);
    
    double baseA = random->GetDouble(-1000, 1000);
    double baseB = random->GetDouble(-1000, 1000);
    
    // Set wmem to base values
    (sA == 1) ? gmem.Set(posA, emp::to_string(baseA)) : gmem.Set(posA, baseA);
    (sB == 1) ? wmem.Set(posB, emp::to_string(baseB)) : wmem.Set(posB, baseB);

    // drop some vectors in
    if (sA == 2) gmem.Set(posA, emp::vector<std::string>{"hello", "world"});
    if (sB == 2) wmem.Set(posB, emp::vector<double>{random->GetDouble(-100, 100), random->GetDouble(-100, 100), random->GetDouble(-100, 100), random->GetDouble(-100, 100), random->GetDouble(-100, 100)});

    cpu.SingleProcess();

    REQUIRE(gmem.GetPos(posA) == wmem.GetPos(posB));
  }
}



TEST_CASE("Inst_TestMemEqu", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 12;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("TestMemEqu", hardware_t::Inst_TestMemEqu, 3, "mem-ANY[C] = mem-ANY[A] == mem-ANY[B])");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());
    prog.PushInst("TestMemEqu", {matrix[posA], matrix[posB], matrix[posC]});
    prog.PushInst("TestMemEqu", {matrix[posC], matrix[posC], matrix[posA]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    size_t sA = random->GetUInt(0, 3);
    size_t sB = random->GetUInt(0, 3);
    size_t sC = random->GetUInt(0, 3);
    size_t sD = random->GetUInt(0, 3);

    double baseA = random->GetDouble(-1000, 1000);
    double baseB = random->GetDouble(-1000, 1000);
    double baseC = random->GetDouble(-1000, 1000);
    double baseD = random->GetDouble(-1000, 1000);

    // Set wmem to base values
    (sA == 1) ? wmem.Set(0, emp::to_string(baseA)) : wmem.Set(0, baseA);
    (sB == 1) ? wmem.Set(1, emp::to_string(baseB)) : wmem.Set(1, baseB);
    (sC == 1) ? wmem.Set(2, emp::to_string(baseC)) : wmem.Set(2, baseC);
    (sD == 1) ? wmem.Set(3, emp::to_string(baseD)) : wmem.Set(3, baseD);
    // drop some vectors in
    if (sA == 2) wmem.Set(0, emp::vector<std::string>{"hello", "world"});
    if (sB == 2) wmem.Set(1, emp::vector<mem_val_t>{wmem.AccessVal(0)});
    if (sC == 2) wmem.Set(2, emp::vector<mem_val_t>{wmem.AccessVal(1), wmem.AccessVal(2), wmem.AccessVal(3)});
    if (sD == 2) wmem.Set(3, emp::vector<std::string>{"blah", "blah", "blah"});
    
    double correct = (double)(wmem.GetPos(posA) == wmem.GetPos(posB));

    cpu.SingleProcess();
    REQUIRE(correct == wmem.AccessVal(posC).GetNum());
    
    cpu.SingleProcess();
    REQUIRE(wmem.AccessVal(posA).GetNum() == 1.0);
  }
}



TEST_CASE("Inst_TestMemNEqu", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 12;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("TestMemNEqu", hardware_t::Inst_TestMemNEqu, 3, "mem-ANY[C] = mem-ANY[A] != mem-ANY[B])");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());
    prog.PushInst("TestMemNEqu", {matrix[posA], matrix[posB], matrix[posC]});
    prog.PushInst("TestMemNEqu", {matrix[posC], matrix[posC], matrix[posA]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    size_t sA = random->GetUInt(0, 3);
    size_t sB = random->GetUInt(0, 3);
    size_t sC = random->GetUInt(0, 3);
    size_t sD = random->GetUInt(0, 3);

    double baseA = random->GetDouble(-1000, 1000);
    double baseB = random->GetDouble(-1000, 1000);
    double baseC = random->GetDouble(-1000, 1000);
    double baseD = random->GetDouble(-1000, 1000);

    // Set wmem to base values
    (sA == 1) ? wmem.Set(0, emp::to_string(baseA)) : wmem.Set(0, baseA);
    (sB == 1) ? wmem.Set(1, emp::to_string(baseB)) : wmem.Set(1, baseB);
    (sC == 1) ? wmem.Set(2, emp::to_string(baseC)) : wmem.Set(2, baseC);
    (sD == 1) ? wmem.Set(3, emp::to_string(baseD)) : wmem.Set(3, baseD);
    // drop some vectors in
    if (sA == 2) wmem.Set(0, emp::vector<std::string>{"hello", "world"});
    if (sB == 2) wmem.Set(1, emp::vector<mem_val_t>{wmem.AccessVal(0)});
    if (sC == 2) wmem.Set(2, emp::vector<mem_val_t>{wmem.AccessVal(1), wmem.AccessVal(2), wmem.AccessVal(3)});
    if (sD == 2) wmem.Set(3, emp::vector<std::string>{"blah", "blah", "blah"});
    
    double correct = (double)(wmem.GetPos(posA) != wmem.GetPos(posB));

    cpu.SingleProcess();
    REQUIRE(correct == wmem.AccessVal(posC).GetNum());
    
    cpu.SingleProcess();
    REQUIRE(wmem.AccessVal(posA).GetNum() == 0.0);
  }
}



TEST_CASE("Inst_MakeVector", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 12;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("MakeVector", hardware_t::Inst_MakeVector, 3, "mem-ANY[C] = mem-ANY[A] == mem-ANY[B])");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    const size_t posA = random->GetUInt(0, matrix.size());
    const size_t posB = random->GetUInt(0, matrix.size());
    const size_t posC = random->GetUInt(0, matrix.size());
    prog.PushInst("MakeVector", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    size_t sA = random->GetUInt(0, 3);
    size_t sB = random->GetUInt(0, 3);
    size_t sC = random->GetUInt(0, 3);
    size_t sD = random->GetUInt(0, 3);

    double baseA = random->GetDouble(-1000, 1000);
    double baseB = random->GetDouble(-1000, 1000);
    double baseC = random->GetDouble(-1000, 1000);
    double baseD = random->GetDouble(-1000, 1000);

    // Set wmem to base values
    (sA == 1) ? wmem.Set(0, emp::to_string(baseA)) : wmem.Set(0, baseA);
    (sB == 1) ? wmem.Set(1, emp::to_string(baseB)) : wmem.Set(1, baseB);
    (sC == 1) ? wmem.Set(2, emp::to_string(baseC)) : wmem.Set(2, baseC);
    (sD == 1) ? wmem.Set(3, emp::to_string(baseD)) : wmem.Set(3, baseD);
    // drop some vectors in
    if (sA == 2) wmem.Set(0, emp::vector<std::string>{"hello", "world"});
    if (sB == 2) wmem.Set(1, emp::vector<double>{random->GetDouble(-100, 100), random->GetDouble(-100, 100)});
    if (sC == 2) wmem.Set(2, emp::vector<double>{0.0, 42.0, 256.0});
    if (sD == 2) wmem.Set(3, emp::vector<std::string>{"blah", "blah", "blah"});
    
    // What should the vector @ posC look like?
    emp::vector<mem_val_t> vec;
    size_t iA = (posA <= posB) ? posA : posB;
    for (size_t k = iA; k <= posB; ++k) {
      if (wmem.IsVec(k)) continue;
      vec.emplace_back(wmem.AccessVal(k));
    }
    cpu.SingleProcess();

    REQUIRE(wmem.IsVec(posC));
    REQUIRE(wmem.AccessVec(posC) == vec);  
  }
}



TEST_CASE("Inst_VecGet", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 12;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("VecGet", hardware_t::Inst_VecGet, 3, "mem-ANY[C] = mem-VEC[A][mem-NUM[B]])");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    size_t posA = random->GetUInt(0, matrix.size());
    size_t posB = random->GetUInt(0, matrix.size());
    while (posB == posA) posB = random->GetUInt(0, matrix.size());
    size_t posC = random->GetUInt(0, matrix.size());
    prog.PushInst("VecGet", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    emp::vector<mem_val_t> vec(random->GetUInt(10));
    for (size_t vi = 0; vi < vec.size(); ++vi) {
      if (random->P(0.5)) {
        vec[vi] = random->GetDouble(-100, 100);
      } else {
        vec[vi] = emp::to_string(random->GetDouble(-100, 100));
      }
    }

    size_t B = random->GetUInt(vec.size());
    
    wmem.Set(posA, vec);
    wmem.Set(posB, B);

    if (vec.size()) {
      mem_val_t C(vec[B]);
      cpu.SingleProcess();
      REQUIRE(wmem.IsVec(posC) == false);
      
      if (C != wmem.AccessVal(posC)) {
        cpu.GetProgram().Print();
        cpu.PrintHardwareState();
      }
      
      REQUIRE(C == wmem.AccessVal(posC));

    } else { // Vec is empty, instruction should do nothing.
      auto pos = wmem.GetPos(posC);
      cpu.SingleProcess();
      REQUIRE(wmem.GetPos(posC) == pos);
    }    
  }
}



TEST_CASE("Inst_VecSet", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 12;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("VecSet", hardware_t::Inst_VecSet, 3, "mem-VEC[A][mem-NUM[B]] = mem-NUM,STR[C]");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    size_t posA = random->GetUInt(0, matrix.size());
    size_t posB = random->GetUInt(0, matrix.size());
    while (posB == posA) posB = random->GetUInt(0, matrix.size());
    size_t posC = random->GetUInt(0, matrix.size());
    while (posC == posA || posC == posB) posC = random->GetUInt(0, matrix.size());
    prog.PushInst("VecSet", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    emp::vector<mem_val_t> vec(random->GetUInt(10));
    for (size_t vi = 0; vi < vec.size(); ++vi) {
      if (random->P(0.5)) {
        vec[vi] = random->GetDouble(-100, 100);
      } else {
        vec[vi] = emp::to_string(random->GetDouble(-100, 100));
      }
    }

    size_t B = random->GetUInt(vec.size());
    mem_val_t C;
    if (random->P(0.5)) {
      C = emp::to_string(random->GetDouble(-100, 100));
    } else {
      C = random->GetDouble(-100, 100);
    }
    
    wmem.Set(posA, vec);
    wmem.Set(posB, B);
    wmem.Set(posC, C);

    if (vec.size()) {
      cpu.SingleProcess();
      REQUIRE(wmem.IsVec(posA) == true);
      REQUIRE(wmem.GetPosType(posB) == hardware_t::MemPosType::NUM);
      REQUIRE(wmem.AccessVec(posA)[B] == C);
    } else { // Vec is empty, instruction should do nothing.
      auto pos = wmem.GetPos(posA);
      cpu.SingleProcess();
      REQUIRE(wmem.GetPos(posA) == pos);
    }    
  }
}



TEST_CASE("Inst_VecLen", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 12;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("VecLen", hardware_t::Inst_VecLen, 2, "mem-ANY[B] = LEN(mem-VEC[A])");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    size_t posA = random->GetUInt(0, matrix.size());
    size_t posB = random->GetUInt(0, matrix.size());
    size_t posC = random->GetUInt(0, matrix.size());
    
    prog.PushInst("VecLen", {matrix[posA], matrix[posB], matrix[posC]});
    prog.PushInst("VecLen", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    emp::vector<mem_val_t> vec(random->GetUInt(128));
    for (size_t vi = 0; vi < vec.size(); ++vi) {
      if (random->P(0.5)) {
        vec[vi] = random->GetDouble(-100, 100);
      } else {
        vec[vi] = emp::to_string(random->GetDouble(-100, 100));
      }
    }

    wmem.Set(posA, vec);
    cpu.SingleProcess();
    REQUIRE(wmem.GetPosType(posB) == hardware_t::MemPosType::NUM);
    REQUIRE(wmem.AccessVal(posB).GetNum() == vec.size());

    // Get rid of vectors and make sure instruction did nothing
    wmem.Set(0, 0);
    wmem.Set(1, 0);
    wmem.Set(2, 0);
    wmem.Set(3, 0);
    cpu.SingleProcess();
    REQUIRE(wmem.AccessVal(0).GetNum() == 0);
    REQUIRE(wmem.AccessVal(1).GetNum() == 0);
    REQUIRE(wmem.AccessVal(2).GetNum() == 0);
    REQUIRE(wmem.AccessVal(3).GetNum() == 0);
  }
}



TEST_CASE("Inst_VecAppend", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 12;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("VecAppend", hardware_t::Inst_VecAppend, 2, "mem-VEC[A].append(mem-NUM,STR[B])");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    size_t posA = random->GetUInt(0, matrix.size());
    size_t posB = random->GetUInt(0, matrix.size());
    size_t posC = random->GetUInt(0, matrix.size());
    
    prog.PushInst("VecAppend", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    emp::vector<mem_val_t> vec(random->GetUInt(128));
    for (size_t vi = 0; vi < vec.size(); ++vi) {
      if (random->P(0.5)) {
        vec[vi] = random->GetDouble(-100, 100);
      } else {
        vec[vi] = emp::to_string(random->GetDouble(-100, 100));
      }
    }

    mem_val_t B;
    if (random->P(0.5)) {
      B = emp::to_string(random->GetDouble(-100, 100));
    } else {
      B = random->GetDouble(-100, 100);
    }

    wmem.Set(posA, vec);
    wmem.Set(posB, B);

    vec.emplace_back(B);
    cpu.SingleProcess();

    if (posA == posB) {
      // No vectors, nothing should happen.
      REQUIRE(wmem.IsVec(posA) == false);
      REQUIRE(wmem.AccessVal(posA) == wmem.AccessVal(posB));

    } else {
      // Did an append happen?

      if (wmem.AccessVec(posA) != vec) {
        cpu.GetProgram().Print();
        cpu.PrintHardwareState();
      }

      REQUIRE(wmem.IsVec(posA));
      REQUIRE(wmem.AccessVec(posA) == vec);  
    }
  }
}



TEST_CASE("Inst_VecPop", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 12;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("VecPop", hardware_t::Inst_VecPop, 2, "mem-ANY[B] = mem-VEC[A].pop()");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    size_t posA = random->GetUInt(0, matrix.size());
    size_t posB = random->GetUInt(0, matrix.size());
    while (posB == posA) posB = random->GetUInt(0, matrix.size());
    size_t posC = random->GetUInt(0, matrix.size());
    
    prog.PushInst("VecPop", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    emp::vector<mem_val_t> vec(random->GetUInt(128));
    for (size_t vi = 0; vi < vec.size(); ++vi) {
      if (random->P(0.5)) {
        vec[vi] = random->GetDouble(-100, 100);
      } else {
        vec[vi] = emp::to_string(random->GetDouble(-100, 100));
      }
    }

    wmem.Set(posA, vec);
    wmem.Set(posB, 0);
    mem_val_t B;
    if (vec.size()) {
      B = vec.back();  
    } else {
      B = 0;
    }

    cpu.SingleProcess();

    if (vec.size()) {
      vec.pop_back();
      REQUIRE(wmem.IsVec(posA));
      REQUIRE(wmem.AccessVec(posA) == vec);
    } 
    REQUIRE(!wmem.IsVec(posB));
    REQUIRE(wmem.AccessVal(posB) == B);
  }
}



TEST_CASE("Inst_VecRemove", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 12;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("VecRemove", hardware_t::Inst_VecRemove, 2, "mem-VEC[A].Remove(mem-NUM[A])");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    size_t posA = random->GetUInt(0, matrix.size());
    size_t posB = random->GetUInt(0, matrix.size());
    while (posB == posA) posB = random->GetUInt(0, matrix.size());
    size_t posC = random->GetUInt(0, matrix.size());
    
    prog.PushInst("VecRemove", {matrix[posA], matrix[posB], matrix[posC]});
    
    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    emp::vector<mem_val_t> vec(random->GetUInt(128));
    for (size_t vi = 0; vi < vec.size(); ++vi) {
      if (random->P(0.5)) {
        vec[vi] = random->GetDouble(-100, 100);
      } else {
        vec[vi] = emp::to_string(random->GetDouble(-100, 100));
      }
    }

    double B = (double)random->GetUInt(0, 128);
    wmem.Set(posA, vec);
    wmem.Set(posB, B);

    cpu.SingleProcess();

    if (B < vec.size()) {
      // VecRemove should have actually removed element B.
      vec.erase(vec.begin()+(size_t)B);
      REQUIRE(wmem.IsVec(posA));
      REQUIRE(wmem.AccessVec(posA) == vec);
    } 
    REQUIRE(wmem.IsVec(posA));
    REQUIRE(wmem.AccessVec(posA) == vec);
  }
}

TEST_CASE("Inst_VecReplaceAll", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 12;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("VecReplaceAll", hardware_t::Inst_VecReplaceAll, 3, "mem-VEC[A].Replace(mem-NUM,STR[B], mem-NUM[C])");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    size_t posA = random->GetUInt(0, matrix.size());
    size_t posB = random->GetUInt(0, matrix.size());
    while (posB == posA) posB = random->GetUInt(0, matrix.size());
    size_t posC = random->GetUInt(0, matrix.size());
    while (posC == posA || posC == posB) posC = random->GetUInt(0, matrix.size());
    
    prog.PushInst("VecReplaceAll", {matrix[posA], matrix[posB], matrix[posC]});

    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    emp::vector<mem_val_t> vec(random->GetUInt(128));
    for (size_t vi = 0; vi < vec.size(); ++vi) {
      if (random->P(0.5)) {
        vec[vi] = random->GetUInt(0, 5);
      } else {
        vec[vi] = emp::to_string(random->GetUInt(0, 5));
      }
    }

    mem_val_t B; 
    B = random->GetUInt(0, 5);
    mem_val_t C;
    C = "REPLACE VAL";
    
    wmem.Set(posA, vec);
    wmem.Set(posB, B);
    wmem.Set(posC, C);

    cpu.SingleProcess();

    for (size_t vi = 0; vi < vec.size(); ++vi) {
      if (vec[vi] == B) vec[vi] = C;
    }

    REQUIRE(wmem.IsVec(posA));
    REQUIRE(vec == wmem.AccessVec(posA));
  }
}

TEST_CASE("Inst_VecIndexOf", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 12;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("VecIndexOf", hardware_t::Inst_VecIndexOf, 3, "");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    size_t posA = random->GetUInt(0, matrix.size());
    size_t posB = random->GetUInt(0, matrix.size());
    while (posB == posA) posB = random->GetUInt(0, matrix.size());
    size_t posC = random->GetUInt(0, matrix.size());
    while (posC == posA || posC == posB) posC = random->GetUInt(0, matrix.size());
    
    prog.PushInst("VecIndexOf", {matrix[posA], matrix[posB], matrix[posC]});

    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    emp::vector<mem_val_t> vec(random->GetUInt(1,128));
    for (size_t vi = 0; vi < vec.size(); ++vi) {
      if (random->P(0.5)) {
        vec[vi] = random->GetUInt(0, 5);
      } else {
        vec[vi] = emp::to_string(random->GetUInt(0, 5));
      }
    }

    size_t pos = random->GetUInt(vec.size());
    mem_val_t B;
    B = "FIND ME";
    vec[pos] = B;

    wmem.Set(posA, vec);
    wmem.Set(posB, B);

    cpu.SingleProcess();

    REQUIRE(wmem.AccessVal(posC).GetNum() == pos);
  }
}

TEST_CASE("Inst_VecOccurrencesOf", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 12;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("VecOccurrencesOf", hardware_t::Inst_VecOccurrencesOf, 3, "");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    size_t posA = random->GetUInt(0, matrix.size());
    size_t posB = random->GetUInt(0, matrix.size());
    while (posB == posA) posB = random->GetUInt(0, matrix.size());
    size_t posC = random->GetUInt(0, matrix.size());
    while (posC == posA || posC == posB) posC = random->GetUInt(0, matrix.size());
    
    prog.PushInst("VecOccurrencesOf", {matrix[posA], matrix[posB], matrix[posC]});

    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    emp::vector<mem_val_t> vec(random->GetUInt(128));
    for (size_t vi = 0; vi < vec.size(); ++vi) {
      if (random->P(0.5)) {
        vec[vi] = random->GetUInt(0, 5);
      } else {
        vec[vi] = emp::to_string(random->GetUInt(0, 5));
      }
    }

    mem_val_t B; 
    B = random->GetUInt(0, 5);
    
    
    wmem.Set(posA, vec);
    wmem.Set(posB, B);
    
    cpu.SingleProcess();

    size_t cnt = 0;
    for (size_t vi = 0; vi < vec.size(); ++vi) {
      if (vec[vi] == B) cnt++;
    }
    REQUIRE(cnt == wmem.AccessVal(posC).GetNum());
  }
}



TEST_CASE("Inst_VecReverse", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 12;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("VecReverse", hardware_t::Inst_VecReverse, 1, "");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    size_t posA = random->GetUInt(0, matrix.size());
    size_t posB = random->GetUInt(0, matrix.size());
    size_t posC = random->GetUInt(0, matrix.size());
    
    prog.PushInst("VecReverse", {matrix[posA], matrix[posB], matrix[posC]});

    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    emp::vector<mem_val_t> vec(random->GetUInt(128));
    for (size_t vi = 0; vi < vec.size(); ++vi) {
      if (random->P(0.5)) {
        vec[vi] = random->GetUInt(0, 5);
      } else {
        vec[vi] = emp::to_string(random->GetUInt(0, 5));
      }
    }
    
    wmem.Set(posA, vec);
    
    cpu.SingleProcess();

    std::reverse(std::begin(vec), std::end(vec));

    REQUIRE(wmem.AccessVec(posA) == vec);
    
  }
}

TEST_CASE("Inst_VecSwapIfLess", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 12;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("VecSwapIfLess", hardware_t::Inst_VecSwapIfLess, 1, "");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    size_t posA = random->GetUInt(0, matrix.size());
    size_t posB = random->GetUInt(0, matrix.size());
    while (posB == posA) posB = random->GetUInt(0, matrix.size());
    size_t posC = random->GetUInt(0, matrix.size());
    while (posC == posA || posC == posB) posC = random->GetUInt(0, matrix.size());
    
    prog.PushInst("VecSwapIfLess", {matrix[posA], matrix[posB], matrix[posC]});

    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    emp::vector<mem_val_t> vec(random->GetUInt(1,128));
    for (size_t vi = 0; vi < vec.size(); ++vi) {
      if (random->P(0.5)) {
        vec[vi] = random->GetUInt(0, 5);
      } else {
        vec[vi] = emp::to_string(random->GetUInt(0, 5));
      }
    }

    size_t B = random->GetUInt(vec.size());
    size_t C = random->GetUInt(vec.size());

    wmem.Set(posA, vec);
    wmem.Set(posB, B);
    wmem.Set(posC, C);
  
    // cpu.GetProgram().Print();
    // cpu.PrintHardwareState();

    cpu.SingleProcess();

    if (((vec[B].GetType() == mem_val_t::MemoryType::NUM) ? vec[B].GetNum() : 0) <
        ((vec[C].GetType() == mem_val_t::MemoryType::NUM) ? vec[C].GetNum() : 0)) {
      std::swap(vec[B], vec[C]);
    }

    // if (vec != wmem.AccessVec(posA)) {
    // cpu.PrintHardwareState();
    // }
    REQUIRE(vec == wmem.AccessVec(posA));
  }
}



TEST_CASE("Inst_VecGetFront", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 12;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("VecGetFront", hardware_t::Inst_VecGetFront, 2, "");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    size_t posA = random->GetUInt(0, matrix.size());
    size_t posB = random->GetUInt(0, matrix.size());
    size_t posC = random->GetUInt(0, matrix.size());
    
    prog.PushInst("VecGetFront", {matrix[posA], matrix[posB], matrix[posC]});

    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    emp::vector<mem_val_t> vec(random->GetUInt(1,128));
    for (size_t vi = 0; vi < vec.size(); ++vi) {
      if (random->P(0.5)) {
        vec[vi] = random->GetUInt(0, 5);
      } else {
        vec[vi] = emp::to_string(random->GetUInt(0, 5));
      }
    }

    mem_val_t B(vec.front());
    
    wmem.Set(posA, vec);
    
    cpu.SingleProcess();

    REQUIRE(wmem.AccessVal(posB) == B);
  }
}

TEST_CASE("Inst_VecGetBack", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 13;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;
  using mem_val_t = typename hardware_t::MemoryValue;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("VecGetBack", hardware_t::Inst_VecGetBack, 2, "");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    size_t posA = random->GetUInt(0, matrix.size());
    size_t posB = random->GetUInt(0, matrix.size());
    size_t posC = random->GetUInt(0, matrix.size());
    
    prog.PushInst("VecGetBack", {matrix[posA], matrix[posB], matrix[posC]});

    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    emp::vector<mem_val_t> vec(random->GetUInt(1,128));
    for (size_t vi = 0; vi < vec.size(); ++vi) {
      if (random->P(0.5)) {
        vec[vi] = random->GetUInt(0, 5);
      } else {
        vec[vi] = emp::to_string(random->GetUInt(0, 5));
      }
    }

    mem_val_t B(vec.back());
    
    wmem.Set(posA, vec);
    
    cpu.SingleProcess();

    REQUIRE(wmem.AccessVal(posB) == B);
  }
}

TEST_CASE("Inst_IsStr", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 13;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("IsStr", hardware_t::Inst_IsStr, 2, "");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    size_t posA = random->GetUInt(0, matrix.size());
    size_t posB = random->GetUInt(0, matrix.size());
    size_t posC = random->GetUInt(0, matrix.size());
    
    prog.PushInst("IsStr", {matrix[posA], matrix[posB], matrix[posC]});

    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    size_t sA = random->GetUInt(3);
    if (sA == 0) {
      wmem.Set(posA, random->GetUInt(0, 5));
    } else if (sA == 1) {
      wmem.Set(posA, emp::to_string(random->GetUInt(0, 5)));
    } else {
      wmem.Set(posA, emp::vector<std::string>{emp::to_string(random->GetUInt(0, 5)), emp::to_string(random->GetUInt(0, 5)), emp::to_string(random->GetUInt(0, 5)), emp::to_string(random->GetUInt(0, 5)), emp::to_string(random->GetUInt(0, 5))});
    }

    cpu.SingleProcess();

    if (sA == 1) {
      REQUIRE(wmem.AccessVal(posB).GetNum() == 1.0);
    } else {
      REQUIRE(wmem.AccessVal(posB).GetNum() == 0.0);
    } 
  }
}

TEST_CASE("Inst_IsNum", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 13;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("IsNum", hardware_t::Inst_IsNum, 2, "");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    size_t posA = random->GetUInt(0, matrix.size());
    size_t posB = random->GetUInt(0, matrix.size());
    size_t posC = random->GetUInt(0, matrix.size());
    
    prog.PushInst("IsNum", {matrix[posA], matrix[posB], matrix[posC]});

    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    size_t sA = random->GetUInt(3);
    if (sA == 0) {
      wmem.Set(posA, random->GetUInt(0, 5));
    } else if (sA == 1) {
      wmem.Set(posA, emp::to_string(random->GetUInt(0, 5)));
    } else {
      wmem.Set(posA, emp::vector<std::string>{emp::to_string(random->GetUInt(0, 5)), emp::to_string(random->GetUInt(0, 5)), emp::to_string(random->GetUInt(0, 5)), emp::to_string(random->GetUInt(0, 5)), emp::to_string(random->GetUInt(0, 5))});
    }

    cpu.SingleProcess();

    if (sA == 0) {
      REQUIRE(wmem.AccessVal(posB).GetNum() == 1.0);
    } else {
      REQUIRE(wmem.AccessVal(posB).GetNum() == 0.0);
    } 
  }
}

TEST_CASE("Inst_IsVec", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 13;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  inst_lib->AddInst("IsVec", hardware_t::Inst_IsVec, 2, "");
  for (size_t i = 0; i < 1000; ++i) {
    cpu.Reset(); // Hard reset on virtual CPU
    prog.Clear();
    
    size_t posA = random->GetUInt(0, matrix.size());
    size_t posB = random->GetUInt(0, matrix.size());
    size_t posC = random->GetUInt(0, matrix.size());
    
    prog.PushInst("IsVec", {matrix[posA], matrix[posB], matrix[posC]});

    cpu.SetProgram(prog);
    cpu.CallModule(0);

    callstate_t & state = cpu.GetCurCallState();
    memory_t & wmem = state.GetWorkingMem();

    size_t sA = random->GetUInt(3);
    if (sA == 0) {
      wmem.Set(posA, random->GetUInt(0, 5));
    } else if (sA == 1) {
      wmem.Set(posA, emp::to_string(random->GetUInt(0, 5)));
    } else {
      wmem.Set(posA, emp::vector<std::string>{emp::to_string(random->GetUInt(0, 5)), emp::to_string(random->GetUInt(0, 5)), emp::to_string(random->GetUInt(0, 5)), emp::to_string(random->GetUInt(0, 5)), emp::to_string(random->GetUInt(0, 5))});
    }

    cpu.SingleProcess();

    if (sA == 2) {
      REQUIRE(wmem.AccessVal(posB).GetNum() == 1.0);
    } else {
      REQUIRE(wmem.AccessVal(posB).GetNum() == 0.0);
    } 
  }
}


TEST_CASE("Inst_If", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 13;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;
  using callstate_t = typename hardware_t::CallState;
  using memory_t = typename hardware_t::Memory;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  // InstProperty BEGIN_FLOW END_FLOW MODULE
  inst_lib->AddInst("Inc", hardware_t::Inst_Inc, 1, "");
  inst_lib->AddInst("Dec", hardware_t::Inst_Dec, 1, "");
  inst_lib->AddInst("Nop", hardware_t::Inst_Nop, 0, "");
  inst_lib->AddInst("If", hardware_t::Inst_If, 1, "", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("Close", hardware_t::Inst_Close, 0, "", {inst_lib_t::InstProperty::END_FLOW});
  inst_lib->AddInst("Break", hardware_t::Inst_Break, 0, "");
  
  // ---------------------------------------
  // - TEST: If(true) ... Close ...
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("If", {matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Close", {matrix[0]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 9; ++i) {
  //   std::cout << "\n\n--- AFTER CYCLE " << i << std::endl;
    cpu.SingleProcess();
  //   cpu.PrintHardwareState();
  }
  callstate_t & state = cpu.GetCurCallState();
  memory_t & wmem = state.GetWorkingMem();
  REQUIRE(wmem.AccessVal(0).GetNum() == 1);
  REQUIRE(wmem.AccessVal(1).GetNum() == 2);
  REQUIRE(wmem.AccessVal(2).GetNum() == 4);
  // ---------------------------------------


  // ---------------------------------------
  // - TEST: If(false) ... Close ...
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("Nop", {matrix[0]});
  prog.PushInst("If", {matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Close", {matrix[0]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 6; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  callstate_t & state0 = cpu.GetCurCallState();
  memory_t & wmem0 = state0.GetWorkingMem();
  REQUIRE(wmem0.AccessVal(0).GetNum() == 0);
  REQUIRE(wmem0.AccessVal(1).GetNum() == 0);
  REQUIRE(wmem0.AccessVal(2).GetNum() == 4);
  // ---------------------------------------


  // ---------------------------------------
  // - TEST: If(true) ... - Implicit close
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("If", {matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 8; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << "---" << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  callstate_t & state1 = cpu.GetCurCallState();
  memory_t & wmem1 = state1.GetWorkingMem();
  REQUIRE(wmem1.AccessVal(0).GetNum() == 1);
  REQUIRE(wmem1.AccessVal(1).GetNum() == 2);
  REQUIRE(wmem1.AccessVal(2).GetNum() == 4);
  // ---------------------------------------


  // ---------------------------------------
  // - TEST: If(false) ... - Implicit close
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("Nop", {matrix[0]});
  prog.PushInst("Nop", {matrix[0]});
  prog.PushInst("If", {matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0, true, true);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 16; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << "---" << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(0).GetNum() == 0);
  // ---------------------------------------

  // ---------------------------------------
  // - TEST: If(true) .. break .. - Implicit close
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("If", {matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Break", {});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 4; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << "---" << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(0).GetNum() == 1);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(1).GetNum() == 1);
  // ---------------------------------------

  // ---------------------------------------
  // - TEST: If(true) .. break .. Close
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("If", {matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Break", {});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Close", {});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 8; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << "---" << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(0).GetNum() == 1);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(1).GetNum() == 1);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(2).GetNum() == 4);
  // ---------------------------------------

  // ---------------------------------------
  // - TEST: nested ifs
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("If", {matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("If", {matrix[1]});
  prog.PushInst("Nop", {});
  prog.PushInst("Break", {});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Close", {});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Close", {});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 12; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << "---" << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(0).GetNum() == 1);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(1).GetNum() == 2);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(2).GetNum() == 4);
  // ---------------------------------------

  /////////////////////////////////////
}

TEST_CASE("Inst_While", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 13;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  // InstProperty BEGIN_FLOW END_FLOW MODULE
  inst_lib->AddInst("Inc", hardware_t::Inst_Inc, 1, "");
  inst_lib->AddInst("Dec", hardware_t::Inst_Dec, 1, "");
  inst_lib->AddInst("Nop", hardware_t::Inst_Nop, 0, "");
  inst_lib->AddInst("If", hardware_t::Inst_If, 1, "", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("While", hardware_t::Inst_While, 1, "", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("Close", hardware_t::Inst_Close, 0, "", {inst_lib_t::InstProperty::END_FLOW});
  inst_lib->AddInst("Break", hardware_t::Inst_Break, 0, "");

  // ---------------------------------------
  // - TEST: If(true) ... Close ...
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("While", {matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Dec", {matrix[0]});
  prog.PushInst("Close", {matrix[0]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 20; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(0).GetNum() == 0);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(1).GetNum() == 3);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(2).GetNum() == 4);
  // ---------------------------------------

  // ---------------------------------------
  // - TEST: If(false) ... Close ...
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("Nop", {matrix[0]});
  prog.PushInst("Nop", {matrix[0]});
  prog.PushInst("Nop", {matrix[0]});
  prog.PushInst("While", {matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Dec", {matrix[0]});
  prog.PushInst("Close", {matrix[0]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 8; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(0).GetNum() == 0);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(1).GetNum() == 0);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(2).GetNum() == 4);
  // ---------------------------------------

  // ---------------------------------------
  // - TEST: If(true) ... Implicit close
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("While", {matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Dec", {matrix[0]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 25; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(0).GetNum() == 0);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(1).GetNum() == 3);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(2).GetNum() == 12);
  // ---------------------------------------

  // ---------------------------------------
  // - TEST: If(false) ... Implicit close
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("Nop", {matrix[0]});
  prog.PushInst("Nop", {matrix[0]});
  prog.PushInst("Nop", {matrix[0]});
  prog.PushInst("While", {matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Dec", {matrix[0]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 4; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(0).GetNum() == 0);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(1).GetNum() == 0);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(2).GetNum() == 0);
  // ---------------------------------------

  // ---------------------------------------
  // - TEST: while(true) ... break ... Close ...
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("While", {matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Dec", {matrix[0]});
  prog.PushInst("Break", {});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Close", {matrix[0]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 11; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(0).GetNum() == 2);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(1).GetNum() == 1);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(2).GetNum() == 4);
  // ---------------------------------------

  // ---------------------------------------
  // - TEST: while(true) ... break ... Close ...
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("While", {matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Dec", {matrix[0]});
  prog.PushInst("While", {matrix[1]});
  prog.PushInst("Nop", {});
  prog.PushInst("Dec", {matrix[1]});
  prog.PushInst("Break", {});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Close", {matrix[0]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 37; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(0).GetNum() == 0);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(1).GetNum() == 0);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(2).GetNum() == 12);
  // ---------------------------------------

  ///////////////////////////////////////////
}

TEST_CASE("Inst_Countdown", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 13;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  // InstProperty BEGIN_FLOW END_FLOW MODULE
  inst_lib->AddInst("Inc", hardware_t::Inst_Inc, 1, "");
  inst_lib->AddInst("Dec", hardware_t::Inst_Dec, 1, "");
  inst_lib->AddInst("Nop", hardware_t::Inst_Nop, 0, "");
  inst_lib->AddInst("If", hardware_t::Inst_If, 1, "", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("Countdown", hardware_t::Inst_Countdown, 1, "", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("Close", hardware_t::Inst_Close, 0, "", {inst_lib_t::InstProperty::END_FLOW});
  inst_lib->AddInst("Break", hardware_t::Inst_Break, 0, "");

  // ---------------------------------------
  // - TEST: If(true) ... Close ...
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Countdown", {matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Close", {matrix[0]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 17; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(0).GetNum() == 0);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(1).GetNum() == 3);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(2).GetNum() == 4);
  // ---------------------------------------

}

TEST_CASE("Inst_Foreach", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 13;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  // InstProperty BEGIN_FLOW END_FLOW MODULE
  inst_lib->AddInst("Inc", hardware_t::Inst_Inc, 1, "");
  inst_lib->AddInst("Dec", hardware_t::Inst_Dec, 1, "");
  inst_lib->AddInst("Nop", hardware_t::Inst_Nop, 0, "");
  inst_lib->AddInst("If", hardware_t::Inst_If, 1, "", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("MakeVector", hardware_t::Inst_MakeVector, 3, "");
  inst_lib->AddInst("Foreach", hardware_t::Inst_Foreach, 1, "", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("Close", hardware_t::Inst_Close, 0, "", {inst_lib_t::InstProperty::END_FLOW});
  inst_lib->AddInst("Break", hardware_t::Inst_Break, 0, "");

  // ---------------------------------------
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();
  prog.PushInst("Dec", {matrix[0]});
  prog.PushInst("Dec", {matrix[1]});
  prog.PushInst("Dec", {matrix[1]});
  prog.PushInst("Dec", {matrix[2]});
  prog.PushInst("Dec", {matrix[2]});
  prog.PushInst("Dec", {matrix[2]});
  prog.PushInst("Dec", {matrix[3]});
  prog.PushInst("Dec", {matrix[3]});
  prog.PushInst("Dec", {matrix[3]});
  prog.PushInst("Dec", {matrix[3]});
  prog.PushInst("MakeVector", {matrix[0], matrix[3], matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Foreach", {matrix[3], matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Close", {});
  prog.PushInst("Nop", {});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 27; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(1).GetNum() ==  4);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(2).GetNum() == -3);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(3).GetNum() == -4);
  // ---------------------------------------
  // exit(-1);
  // ---------------------------------------
  // - TEST: If(true) ... Close ...
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("Dec", {matrix[0]});
  prog.PushInst("Dec", {matrix[1]});
  prog.PushInst("Dec", {matrix[1]});
  prog.PushInst("Dec", {matrix[2]});
  prog.PushInst("Dec", {matrix[2]});
  prog.PushInst("Dec", {matrix[2]});
  prog.PushInst("Dec", {matrix[3]});
  prog.PushInst("Dec", {matrix[3]});
  prog.PushInst("Dec", {matrix[3]});
  prog.PushInst("Dec", {matrix[3]});
  prog.PushInst("MakeVector", {matrix[0], matrix[3], matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Foreach", {matrix[0], matrix[0]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Close", {});
  prog.PushInst("Nop", {});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 18; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(0).GetNum() ==  -1);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(1).GetNum() ==  1);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(2).GetNum() == -3);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(3).GetNum() == -4);
  // ---------------------------------------
  

}

TEST_CASE("Inst_Call", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 13;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  // InstProperty BEGIN_FLOW END_FLOW MODULE
  inst_lib->AddInst("Inc", hardware_t::Inst_Inc, 1, "");
  inst_lib->AddInst("Dec", hardware_t::Inst_Dec, 1, "");
  inst_lib->AddInst("Nop", hardware_t::Inst_Nop, 0, "");
  inst_lib->AddInst("If", hardware_t::Inst_If, 1, "", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("Call", hardware_t::Inst_Call, 1, "");
  inst_lib->AddInst("Return", hardware_t::Inst_Return, 0, "");
  inst_lib->AddInst("Output", hardware_t::Inst_Output, 2, "");
  inst_lib->AddInst("ModuleDef", hardware_t::Inst_Nop, 1, "", {inst_lib_t::InstProperty::MODULE});
  inst_lib->AddInst("Close", hardware_t::Inst_Close, 0, "", {inst_lib_t::InstProperty::END_FLOW});
  inst_lib->AddInst("Break", hardware_t::Inst_Break, 0, "");

  // ---------------------------------------
  // - TEST: If(true) ... Close ...
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("Nop", {});
  prog.PushInst("Return", {matrix[2]});
  prog.PushInst("Output", {matrix[2], matrix[3]});
  prog.PushInst("ModuleDef", {matrix[0]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Call", {matrix[1]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Call", {matrix[2]});
  prog.PushInst("Nop", {matrix[2]});
  prog.PushInst("ModuleDef", {matrix[1]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("Output", {matrix[1], matrix[1]});
  prog.PushInst("ModuleDef", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Output", {matrix[2], matrix[2]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 12; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(0).GetNum() == 2);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(1).GetNum() == 1);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(2).GetNum() == 2);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(3).GetNum() == 0);
  // ---------------------------------------

}

TEST_CASE("Inst_Routine", "[taglgp]") {
  // Constants
  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 13;

  // Convenient aliases
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;

  // Create new random number generator
  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  // Create new instruction library
  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  
  // Create virtual hardware w/inst_lib
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  // Create new program
  program_t prog(inst_lib);

  /////////////////////////////////////
  // Instruction testing
  // InstProperty BEGIN_FLOW END_FLOW MODULE
  inst_lib->AddInst("Inc", hardware_t::Inst_Inc, 1, "");
  inst_lib->AddInst("Dec", hardware_t::Inst_Dec, 1, "");
  inst_lib->AddInst("Nop", hardware_t::Inst_Nop, 0, "");
  inst_lib->AddInst("If", hardware_t::Inst_If, 1, "", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("IfNot", hardware_t::Inst_IfNot, 1, "", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("Call", hardware_t::Inst_Call, 1, "");
  inst_lib->AddInst("Routine", hardware_t::Inst_Routine, 1, "");
  inst_lib->AddInst("While", hardware_t::Inst_While, 1, "");
  inst_lib->AddInst("Return", hardware_t::Inst_Return, 0, "");
  inst_lib->AddInst("Output", hardware_t::Inst_Output, 2, "");
  inst_lib->AddInst("ModuleDef", hardware_t::Inst_Nop, 1, "", {inst_lib_t::InstProperty::MODULE});
  inst_lib->AddInst("Close", hardware_t::Inst_Close, 0, "", {inst_lib_t::InstProperty::END_FLOW});
  inst_lib->AddInst("Break", hardware_t::Inst_Break, 0, "");

  // ---------------------------------------
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("Nop", {});
  prog.PushInst("Return", {matrix[2]});
  prog.PushInst("Dec", {matrix[2]});
  prog.PushInst("ModuleDef", {matrix[0]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Routine", {matrix[1]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Routine", {matrix[2]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Nop", {matrix[2]});
  prog.PushInst("ModuleDef", {matrix[1]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("ModuleDef", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 11; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(0).GetNum() == 3);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(1).GetNum() == 1);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(2).GetNum() == 2);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(3).GetNum() == 0);
  // ---------------------------------------

  // ---------------------------------------
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();

  prog.PushInst("If", {matrix[2]});
  prog.PushInst("Nop", {matrix[2]});
  prog.PushInst("While", {matrix[2]});
  prog.PushInst("Nop", {matrix[2]});
  prog.PushInst("Return", {matrix[2]});
  prog.PushInst("Close", {matrix[2]});
  prog.PushInst("Close", {matrix[2]});
  prog.PushInst("Routine", {matrix[2]});
  prog.PushInst("ModuleDef", {matrix[0]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Routine", {matrix[1]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Routine", {matrix[2]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Nop", {matrix[2]});
  prog.PushInst("ModuleDef", {matrix[1]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("ModuleDef", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  prog.PushInst("Inc", {matrix[2]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 14; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(0).GetNum() == 3);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(1).GetNum() == 1);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(2).GetNum() == 2);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(3).GetNum() == 0);
  // ---------------------------------------

  // ---------------------------------------
  cpu.Reset(); // Hard reset on virtual CPU
  prog.Clear();


  prog.PushInst("IfNot", {matrix[3]});
  prog.PushInst("Return", {});
  prog.PushInst("Close", {matrix[0]});
  prog.PushInst("Dec", {matrix[3]});
  prog.PushInst("Close", {matrix[0]});
  prog.PushInst("Nop", {matrix[0]});
  prog.PushInst("Routine", {matrix[1]});

  prog.PushInst("ModuleDef", {matrix[0]});
  prog.PushInst("Inc", {matrix[3]});
  prog.PushInst("Inc", {matrix[3]});
  prog.PushInst("Inc", {matrix[3]});
  prog.PushInst("Inc", {matrix[3]});
  prog.PushInst("Inc", {matrix[3]});
  prog.PushInst("Routine", {matrix[1]});
  prog.PushInst("Inc", {matrix[0]});
  prog.PushInst("Nop", {matrix[0]});

  prog.PushInst("ModuleDef", {matrix[1]});
  prog.PushInst("Inc", {matrix[1]});
  prog.PushInst("If", {matrix[3]});
  
  cpu.SetProgram(prog);
  cpu.CallModule(0);

  // std::cout << "\n\n---PROGRAM---\n\n";
  // cpu.GetProgram().Print();
  // std::cout << "\n\n--- INITIAL STATE ---" << std::endl;
  // cpu.PrintHardwareState();
  for (size_t i = 0; i < 47; ++i) {
    // std::cout << "\n\n--- AFTER CYCLE " << i << std::endl;
    cpu.SingleProcess();
    // cpu.PrintHardwareState();
  }
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(0).GetNum() == 1);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(1).GetNum() == 6);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(2).GetNum() == 0);
  REQUIRE(cpu.GetCurCallState().GetWorkingMem().AccessVal(3).GetNum() == 0);
  // ---------------------------------------

}


TEST_CASE("RandomPrograms", "[taglgp]") {

  constexpr size_t TAG_WIDTH = 4;
  constexpr int seed = 2;

  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using program_t = typename hardware_t::program_t;
  using inst_t = typename hardware_t::inst_t;
  using inst_lib_t = TagLGP::InstLib<hardware_t>;

  emp::Ptr<emp::Random> random = emp::NewPtr<emp::Random>(seed);

  emp::Ptr<inst_lib_t> inst_lib = emp::NewPtr<inst_lib_t>();
  hardware_t cpu(inst_lib, random);

  // Configure CPU
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  cpu.SetMemSize(TAG_WIDTH);
  cpu.SetMemTags(matrix);

  std::cout << "Empty instruction library:" << std::endl;
  inst_lib->Print(); std::cout << std::endl;

  inst_lib->AddInst("Add", hardware_t::Inst_Add, 3, "");
  inst_lib->AddInst("Sub", hardware_t::Inst_Sub, 3, "");
  inst_lib->AddInst("Mult", hardware_t::Inst_Mult, 3, "");
  inst_lib->AddInst("Inc", hardware_t::Inst_Inc, 1, "");
  inst_lib->AddInst("Dec", hardware_t::Inst_Dec, 1, "");
  inst_lib->AddInst("Div", hardware_t::Inst_Div, 3, "");
  inst_lib->AddInst("Mod", hardware_t::Inst_Mod, 3, "");
  inst_lib->AddInst("Not", hardware_t::Inst_Not, 3, ""); 
  inst_lib->AddInst("TestNumEqu", hardware_t::Inst_TestNumEqu, 3, "");
  inst_lib->AddInst("TestNumNEqu", hardware_t::Inst_TestNumNEqu, 3, "");
  inst_lib->AddInst("TestNumLess", hardware_t::Inst_TestNumLess, 3, "");
  inst_lib->AddInst("Floor", hardware_t::Inst_Floor, 3, "");
  
  inst_lib->AddInst("CopyMem", hardware_t::Inst_CopyMem, 3, "");
  inst_lib->AddInst("SwapMem", hardware_t::Inst_SwapMem, 3, "");
  inst_lib->AddInst("Input", hardware_t::Inst_Input, 3, "");
  inst_lib->AddInst("Output", hardware_t::Inst_Output, 3, "");
  inst_lib->AddInst("CommitGlobal", hardware_t::Inst_CommitGlobal, 3, "");
  inst_lib->AddInst("PullGlobal", hardware_t::Inst_PullGlobal, 3, "");
  inst_lib->AddInst("TestMemEqu", hardware_t::Inst_TestMemEqu, 3, "");
  inst_lib->AddInst("TestMemNEqu", hardware_t::Inst_TestMemNEqu, 3, "");

  inst_lib->AddInst("MakeVector", hardware_t::Inst_MakeVector, 3, "");
  inst_lib->AddInst("VecGet", hardware_t::Inst_VecGet, 3, "");
  inst_lib->AddInst("VecSet", hardware_t::Inst_VecSet, 3, "");
  inst_lib->AddInst("VecLen", hardware_t::Inst_VecLen, 3, "");
  inst_lib->AddInst("VecAppend", hardware_t::Inst_VecAppend, 3, "");
  inst_lib->AddInst("VecPop", hardware_t::Inst_VecPop, 3, "");
  inst_lib->AddInst("VecRemove", hardware_t::Inst_VecRemove, 3, "");
  inst_lib->AddInst("VecReplaceAll", hardware_t::Inst_VecReplaceAll, 3, "");
  inst_lib->AddInst("VecIndexOf", hardware_t::Inst_VecIndexOf, 3, "");
  inst_lib->AddInst("VecOccurrencesOf", hardware_t::Inst_VecOccurrencesOf, 3, "");
  inst_lib->AddInst("VecReverse", hardware_t::Inst_VecReverse, 3, "");
  inst_lib->AddInst("VecSwapIfLess", hardware_t::Inst_VecSwapIfLess, 3, "");
  inst_lib->AddInst("VecGetFront", hardware_t::Inst_VecGetFront, 3, "");
  inst_lib->AddInst("VecGetBack", hardware_t::Inst_VecGetBack, 3, "");

  inst_lib->AddInst("IsStr", hardware_t::Inst_IsStr, 3, "");
  inst_lib->AddInst("IsNum", hardware_t::Inst_IsNum, 3, "");
  inst_lib->AddInst("IsVec", hardware_t::Inst_IsVec, 3, "");

  inst_lib->AddInst("If", hardware_t::Inst_If, 3, "", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("IfNot", hardware_t::Inst_IfNot, 3, "", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("While", hardware_t::Inst_While, 3, "", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("Countdown", hardware_t::Inst_Countdown, 3, "", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("Foreach", hardware_t::Inst_Foreach, 3, "", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("Close", hardware_t::Inst_Close, 3, "", {inst_lib_t::InstProperty::END_FLOW});
  inst_lib->AddInst("Break", hardware_t::Inst_Break, 3, "");
  inst_lib->AddInst("Call", hardware_t::Inst_Call, 3, "");
  inst_lib->AddInst("Routine", hardware_t::Inst_Routine, 3, "");
  inst_lib->AddInst("Return", hardware_t::Inst_Return, 3, "");

  inst_lib->AddInst("ModuleDef", hardware_t::Inst_Nop, 3, "", {inst_lib_t::InstProperty::MODULE});

  inst_lib->AddInst("Nop", hardware_t::Inst_Nop, 3, "");

  // Add Terminals
  for (size_t i = 0; i <= 16; ++i) {
    inst_lib->AddInst("Set-" + emp::to_string(i),
      [i](hardware_t & hw, const inst_t & inst) {
        hardware_t::CallState & state = hw.GetCurCallState();
        hardware_t::Memory & wmem = state.GetWorkingMem();
        
        size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
        if (!hw.IsValidMemPos(posA)) return; // Do nothing

        wmem.Set(posA, (double)i);

      });
  }

  inst_lib->AddInst("Set-A",
    [](hardware_t & hw, const inst_t & inst) {
      hardware_t::CallState & state = hw.GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posA)) return; // Do nothing

      wmem.Set(posA, "A");

    });

  inst_lib->AddInst("Set-B",
    [](hardware_t & hw, const inst_t & inst) {
      hardware_t::CallState & state = hw.GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posA)) return; // Do nothing

      wmem.Set(posA, "A");

    });

  inst_lib->AddInst("Set-C",
    [](hardware_t & hw, const inst_t & inst) {
      hardware_t::CallState & state = hw.GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();
      
      size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
      if (!hw.IsValidMemPos(posA)) return; // Do nothing

      wmem.Set(posA, "A");
    });
  // Print populated instruction library
  // std::cout << "Populated instruction library:" << std::endl;
  // inst_lib->Print(); std::cout << std::endl;

  // Run a bunch of randomly generated programs for 1024 updates
  for (size_t p = 0; p < 10000; ++p) {
    cpu.Reset(); // Hard reset on virtual CPU

    std::cout << "Testing random program #" << p << "..." << std::endl;
    program_t prg(inst_lib);
    size_t N = random->GetUInt(512, 2046);
    for (size_t i = 0; i < N; ++i) {
      prg.PushInst(TagLGP::GenRandTagGPInst(*random, *inst_lib));
    } 
    cpu.SetProgram(prg);
    cpu.CallModule(0); 

    // if (p==252) {
      // std::cout << "\n\n----- PROGRAM ----" << std::endl;
      // cpu.GetProgram().Print();
      // std::cout << "\n\n----- PROGRAM (by modules) ----" << std::endl;
      // cpu.PrintModuleSequences();
    // }

    for (size_t i = 0; i < 4096; ++i) {
      // std::cout << ">> Cycle: " << i << std::endl;
      cpu.SingleProcess();
    }
    std::cout << " ...Done." << std::endl;
  }
}

