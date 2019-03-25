# Genetic Programming System

Here, we provide a more detailed description of the linear GP system we used in this work. This document also reports configuration details used in this work. Our exact configuration files along with experiment source code can be found in our [online GitHub repository](https://github.com/amlalejini/GECCO-2019-cohort-lexicase).

**Contents**

<!-- TOC -->

- [GP Representation](#gp-representation)
  - [Tag-based Referencing](#tag-based-referencing)
  - [Programs](#programs)
  - [Virtual CPU](#virtual-cpu)
  - [Tag-accessed Memory](#tag-accessed-memory)
  - [Module Execution](#module-execution)
- [Instruction Set](#instruction-set)
  - [Default Instructions](#default-instructions)
  - [String-specific Instructions](#string-specific-instructions)
  - [Problem-specific Instructions](#problem-specific-instructions)
    - [Problem - Number IO](#problem---number-io)
    - [Problem - For Loop Index](#problem---for-loop-index)
    - [Problem - Compare String Lengths](#problem---compare-string-lengths)
    - [Problem - Grade](#problem---grade)
    - [Problem - Median](#problem---median)
    - [Problem - Smallest](#problem---smallest)
- [Evolution](#evolution)
- [References](#references)

<!-- /TOC -->

## GP Representation

In this work, we used a tag-based linear GP representation.
Our representation's extensive use of tag-based referencing makes it unique relative to traditional linear GP representations.

### Tag-based Referencing

Tags are evolvable labels that can be mutated, and the similarity (or dissimilarity) between any two tags is quantifiable.
Spector et al. initially demonstrated tags in GP as a mechanism for labeling and referring to program modules (Spector et al., 2011a; Spector et al., 2011b; Spector et al., 2012b).
The use of tags in GP was further expanded in SignalGP (Lalejini and Ofria, 2018) where tag-based referencing was used to facilitate signal-driven program execution.
Tags allow for _inexact_ referencing; a referring tag can always link to the program module labeled with the most similar tag, ensuring that all possible tags are valid references.
Here, tags are represented as length-16 bit strings, and the Hamming distance between two tags specifies their dissimilarity.
As in Spector et al.'s previous work (Spector et al., 2011a; Spector et al., 2011b; Spector et al., 2012b) and in SignalGP (Lalejini and Ofria, 2018), we use tags to label and reference program modules; however, we take the use of tags one step further, using them to specify locations in memory.

### Programs

Programs in our representation are linear sequences of instructions, and each instruction has up to three tag-based arguments that may modify its behavior.
Our instruction set supports the declaration of modules, which are named (tagged) segments of code that can be triggered (i.e., called by other instructions) during execution.
Modules are declared and labeled by `Module-Def` instructions whose tag-based arguments are used to tag the declared module.

### Virtual CPU

Programs are executed in the context of a virtual CPU, which maintains a call stack that stores state information about the program's active module calls (e.g., local memory contents and a pointer to the current position in the program).
The virtual CPU also gives programs access to four types of memory buffers, each consisting of 16 statically tagged positions: input memory, output memory, working memory, and global memory.
Input, output, and working memory are local to a particular module call, while changes to global memory are shared across module calls.
Input memory is used to pass information to a triggered module from the calling state, and output memory is used to return information from triggered modules back to the calling state.
Working memory is used for performing local operations (e.g., addition, multiplication, etc.).

### Tag-accessed Memory

Many traditional GP systems that give genetic programs access to memory (_e.g._, indexable memory registers) use rigid naming schemes where memory is numerically indexed, and mutation operators must guarantee the validity of memory-referencing instructions. Tag-accessed memory allows programs to use tag-based referencing to index into memory registers. Tags are evolvable labels that give genetic programs a flexible mechanism for specification. Tags allow for _inexact_ referencing; a referring tag references the _closest matching_ referent. To facilitate inexact referencing, the similarity (or dissimilarity) between any two tags must be quantifiable; thus, a referring tag can _always_ reference the closest matching referent tag. This ensures that all possible tags are valid references.
When using a tag-accessed memory model, each of the 16 memory registers in any of the four types of memory buffers (input, output, working, and global) are
statically tagged with length-16 bit strings.
Tags used for memory registers were generated using the [Hadamard matrix](https://en.wikipedia.org/wiki/Hadamard_matrix) and were as follows:

- Register 0:  `1111111111111111`
- Register 1:  `0101010101010101`
- Register 2:  `0011001100110011`
- Register 3:  `1001100110011001`
- Register 4:  `0000111100001111`
- Register 5:  `1010010110100101`
- Register 6:  `1100001111000011`
- Register 7:  `0110100101101001`
- Register 8:  `0000000011111111`
- Register 9:  `1010101001010101`
- Register 10: `1100110000110011`
- Register 11: `0110011010011001`
- Register 12: `1111000000001111`
- Register 13: `0101101010100101`
- Register 14: `0011110011000011`
- Register 15: `1001011001101001`

Note that register 0 has the same tag in each of the four types of memory buffers (input, output, working, and global).
The same is true for all other registers.

Each program instruction has three _tag_ arguments (_i.e._, each instruction argument is a length-16 bit string). Tag-based instruction arguments reference the memory position with the closest matching tag; as such, argument tags need not _exactly_ match any of the tags with memory
positions.

Because many problems in the general program synthesis benchmark suite (Helmuth and Spector, 2015) require the use of a variety of data types,
our virtual CPU allows all positions in memory to contain any one of three data types: numbers, strings, or lists (which can contain a mix of numbers and strings).
Instruction arguments specify positions in memory using tag-based referencing, searching for and retrieving the contents of the closest matching memory position _of the required type_.
For example, the `Add` instruction's first two arguments specify which two numbers in memory to add together, and the third argument specifies where to store the result.
If there are no numbers in memory, the instruction does nothing.
The same is true for instruction arguments that specify strings or lists.

In this work, we mutated all tag arguments at a per-bit rate (of 0.001). The tags on memory registers never changed.

### Module Execution

When a module is triggered by a `Call` instruction (using the instruction's first argument tag to reference the module with the closest matching tag-based label), a new call state is created and pushed onto the virtual CPU's call stack.
The working memory of the calling state is copied into the input memory of the new call state (i.e., the arguments to the called module are the full contents of the calling state's working memory); the working and output memory of the new call state are initially empty.
Module calls may return by either executing a `Return` instruction or by reaching the end of the module's instruction sequence.
When a module call returns, its call state is popped from the call stack, and anything stored in the output memory of the returning call state is copied to the working memory of the caller state (otherwise leaving the state of the caller's working memory unchanged).
Alternatively, modules can be triggered by executing a `Routine` instruction, which behaves in the same way as the `Call` instruction, except the called module shares its input, working, and output memory with the caller.

## Instruction Set

Below, we describe our default instruction set (used across all problems), and all problem-specific instructions.

### Default Instructions

Note:

- EOP: End of program
- wReg: Working memory register
- iReg: Input memory register
- oReg: Output memory register
- gReg: Global memory register
- wReg[0]/iReg[0]/oReg[0]/gReg[0] indicates the value at the register specified by an instruction's first _argument_ (either tag-based or numeric), Reg[1] indicates the value at the register specified by an instruction's second argument, and Reg[2] indicates the value at the register specified by the instruction's third argument.
    - wReg[0], wReg[1], _etc_: working register 0, working register 1, _etc._

Instructions that would produce undefined behavior (e.g., division by zero) are treated as no operations.

| Instruction | Arguments Used | Description |
| :--- | :---: | :--- |
| `Add` | 3 | wReg[2] = wReg[0] + wReg[1] |
| `Sub` | 3 | wReg[2] = wReg[0] - wReg[1] |
| `Mult`  | 3 | wReg[2] = wReg[0] * wReg[1] |
| `Div` | 3 | wReg[2] = wReg[0] / wReg[1] |
| `Mod` | 3 | wReg[2] = wReg[0] % wReg[1] |
| `TestNumEqu`  | 3 | wReg[2] = wReg[0] == wReg[1] |
| `TestNumNEqu` | 3 | wReg[2] = wReg[0] != wReg[1] |
| `TestNumLess` | 3 | wReg[2] = wReg[0] < wReg[1] |
| `TestNumLessTEqu` | 3 | wReg[2] = wReg[0] <= wReg[1] |
| `TestNumGreater`  | 3 | wReg[2] = wReg[0] > wReg[1] |
| `TestNumGreaterTEqu`  | 3 | wReg[2] = wReg[0] >= wReg[1] |
| `Floor` | 1 | Floor(wReg[0]) |
| `Not` | 1 | wReg[0] = !wReg[0] |
| `Inc` | 1 | wReg[0] = wReg[0] + 1 |
| `Dec` | 1 | wReg[0] = wReg[0] - 1 |
| `CopyMem` | 2 | wReg[0] = wReg[1] |
| `SwapMem` | 2 | Swap(wReg[0], wReg[1]) |
| `Input` |  |  |
| `Output` |    |  |
| `CommitGlobal` |    |  |
| `PullGlobal` |    |  |
| `TestMemEqu` |    |  |
| `TestmemNEqu` |    |  |
| `If`      | 1   | If wReg[0] != 0, proceed. Otherwise skip to the next `Close` or EOP. |
| `IfNot` | 1     | If wReg[0] == 0, proceed. Otherwise skip to the next `Close` or EOP. |
| `While` | 1     | While wReg[0] != 0, loop. Otherwise skip to next `Close` or EOP. |
| `Countdown` | 1 | Same as `While`, but decrement wReg[0] if wReg[0] != 0. |
| `Foreach` |    |  |
| `Close` | 0     | Indicate the end of a control block of code (e.g., loop, if). |
| `Break` | 0     | Break out of current control flow (e.g., loop). |
| `Call` |    |  |
| `Routine` |    |  |
| `Return`  | 0    | Return from program execution (exit program execution). |
| `ModuleDef` |    |  |
| `MakeVector` |    |  |
| `VecGet` |    |  |
| `VecSet` |    |  |
| `VecLen` |    |  |
| `VecAppend` |    |  |
| `VecPop` |    |  |
| `VecRemove` |    |  |
| `VecReplaceAll` |    |  |
| `VecIndexOf` |    |  |
| `VecOccurencesOf` |    |  |
| `VecReverse` |    |  |
| `VecSwapIfLess` |    |  |
| `VecGetFront` |    |  |
| `VecGetBack` |    |  |
| `IsNum` |    |  |
| `IsVec` |    |  |
| `Set-0` | 1 | wReg[0] = 0|
| `Set-1` | 1 | wReg[0] = 1|
| `Set-2` | 1 | wReg[0] = 2|
| `Set-3` | 1 | wReg[0] = 3|
| `Set-4` | 1 | wReg[0] = 4|
| `Set-5` | 1 | wReg[0] = 5|
| `Set-6` | 1 | wReg[0] = 6|
| `Set-7` | 1 | wReg[0] = 7|
| `Set-8` | 1 | wReg[0] = 8|
| `Set-9` | 1 | wReg[0] = 9|
| `Set-10` | 1 | wReg[0] = 10|
| `Set-11` | 1 | wReg[0] = 11|
| `Set-12` | 1 | wReg[0] = 12 |
| `Set-13` | 1 | wReg[0] = 13 |
| `Set-14` | 1 | wReg[0] = 14 |
| `Set-15` | 1 | wReg[0] = 15 |
| `Set-16` | 1 | wReg[0] = 16 |

### String-specific Instructions

These instructions were included only in problems where programs needed to manipulate
strings.

| Instruction | Arguments Used | Description |
| :--- | :---: | :--- |
| `StrLength` |  |  |
| `StrConcat` |  |  |
| `IsStr` |  |  |

### Problem-specific Instructions

#### Problem - Number IO

| Instruction | # Arguments Used | Description |
| :--- | :---: | :--- |
| `LoadInt` | 1 | wReg[0] = integer input |
| `LoadDouble` | 1 | wReg[0] = double input |
| `SubmitNum` | 1 | Output wReg[0] |


#### Problem - For Loop Index

| Instruction | # Arguments Used | Description |
| :--- | :---: | :--- |
| `LoadStart` | 1 | wReg[0] = start input|
| `LoadEnd` | 1 | wReg[0] = end input |
| `LoadStep` | 1 | wReg[0] = step input |
| `SubmitNum` | 1 | Output wReg[0] |

#### Problem - Compare String Lengths

| Instruction | # Arguments Used | Description |
| :--- | :---: | :--- |
| `LoadStr1` | 1 | wReg[0] = string-1 |
| `LoadStr2` | 1 | wReg[0] = string-2 |
| `LoadStr3` | 1 | wReg[0] = string-3 |
| `SubmitTrue` | 0 | Output true |
| `SubmitFalse` | 0 | Output false |
| `SubmitVal` | 1 | Output logical value at wReg[0] |

#### Problem - Grade

| Instruction | # Arguments Used | Description |
| :--- | :---: | :--- |
| `SubmitA` | 0 | Classify grade as 'A' |
| `SubmitB` | 0 | Classify grade as 'B' |
| `SubmitC` | 0 | Classify grade as 'C' |
| `SubmitD` | 0 | Classify grade as 'D' |
| `SubmitF` | 0 | Classify grade as 'F' |
| `LoadThreshA` | 1 | wReg[0] = input threshold for 'A' |
| `LoadThreshB` | 1 | wReg[0] = input threshold for 'B' |
| `LoadThreshC` | 1 | wReg[0] = input threshold for 'C' |
| `LoadThreshD` | 1 | wReg[0] = input threshold for 'D' |
| `LoadGrade` | 1 | wReg[0] = input grade to classify|

#### Problem - Median

| Instruction | # Arguments Used | Description |
| :--- | :---: | :--- |
| `LoadNum1` | 1 | wReg[0] = input 1|
| `LoadNum2` | 1 | wReg[0] = input 2|
| `LoadNum3` | 1 | wReg[0] = input 3|
| `SubmitNum` | 1 | Output wReg[0]|

#### Problem - Smallest

| Instruction | # Arguments Used | Description |
| :--- | :---: | :--- |
| `LoadNum1` | 1 | wReg[0] = input 1|
| `LoadNum2` | 1 | wReg[0] = input 2|
| `LoadNum3` | 1 | wReg[0] = input 3|
| `LoadNum4` | 1 | wReg[0] = input 4|
| `SubmitNum` | 1 | Output wReg[0] |

## Evolution

Across all of our experiments, we evolved populations of 512 programs for varied numbers of generations depending on the experimental treatment.
We propagated programs asexually and applied mutations to offspring.
We applied single-instruction substitutions, insertions, and deletions at a per-instruction rate of 0.005 each and multi-instruction sequence duplications and deletions (slip mutations~\citep{Lalejini2017}) at a per-program rate of 0.05.
We applied single-module duplications and deletions at a per-module rate of 0.05, and we mutated argument tags at a per-bit rate of 0.001.
During mutation, we limited the maximum length of programs (which varied by problem).

## References

Helmuth, T., & Spector, L. (2015). General Program Synthesis Benchmark Suite. In Proceedings of the 2015 on Genetic and Evolutionary Computation Conference - GECCO ’15 (pp. 1039–1046). New York, New York, USA: ACM Press. https://doi.org/10.1145/2739480.2754769

Lalejini, A., & Ofria, C. (2018). Evolving event-driven programs with SignalGP. In Proceedings of the Genetic and Evolutionary Computation Conference on - GECCO ’18 (pp. 1135–1142). New York, New York, USA: ACM Press. https://doi.org/10.1145/3205455.3205523

Spector, L., Martin, B., Harrington, K., & Helmuth, T. (2011a). Tag-based modules in genetic programming. In Proceedings of the 13th annual conference on Genetic and evolutionary computation - GECCO ’11 (p. 1419). New York, New York, USA: ACM Press. https://doi.org/10.1145/2001576.2001767

Spector, L., Harrington, K., Martin, B., & Helmuth, T. (2011b). What’s in an Evolved Name? The Evolution of Modularity via Tag-Based Reference. In R. Riolo, E. Vladislavleva, & J. H. Moore (Eds.), Genetic Programming Theory and Practice IX (pp. 1–16). New York, NY: Springer New York. https://doi.org/10.1007/978-1-4614-1770-5_1

Spector, L., Harrington, K., & Helmuth, T. (2012). Tag-based modularity in tree-based genetic programming. In Proceedings of the fourteenth international conference on Genetic and evolutionary computation conference - GECCO ’12 (p. 815). New York, New York, USA: ACM Press. https://doi.org/10.1145/2330163.2330276