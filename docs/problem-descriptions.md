# Program Synthesis Problems

We used five test-based problems from Helmuth and Spector's general program synthesis
benchmark suite (Helmuth and Spector, 2015) to evaluate the utility of our lexicase
selection variants:

- [Small or Large](#problem---small-or-large)
- [For Loop Index](#problem---for-loop-index)
- [Compare String Lengths](#problem---compare-string-lengths)
- [Median](#problem---median)
- [Smallest](#problem---smallest)

All of the problems in the general program synthesis benchmark suite were selected
from sources for introductory computer science programming problems; while not particularly
challenging for experienced human programmers, they can be challenging for current
GP systems (Helmuth and Spector, 2015; Forstenlechner et al., 2018).
The suite of benchmarks was previously used to compare lexicase selection against
other, more traditional selection schemes (Helmuth and Spector, 2015) using PushGP.

Each problem is defined by a set of test cases in which programs are given specified
input data and are scored on how close their output is to the correct output; depending
on the problem, we measured scores on a gradient or on a binary pass-fail basis.
During an evaluation, we limited the number of steps (instructions) a program could
execute; this limit varied by problem. Programs were not required to stop on their
own as long as they output their results before reaching their execution limit.
For each problem, we added [problem-specific instructions](./gp-system.md#problem-specific-instructions) to the instruction set to
allow programs to lead test cases and to submit output.

Below, we provide a brief description of each of the [five problems](#problem-descriptions)
as well as description of [how we configured the training and testing sets for each](#trainingtesting-sets).

For a more fully detailed description of these problems, see (Helmuth and Spector, 2015).

## Problem Descriptions

### Problem - Small or Large

Programs receive an input _n_ (-10,000 &le; _n_ &le; 10,000) and must classify _n_
as 'small' (_n_ < 1,000), 'large' (_n_ &ge; 2,000), or 'neither' (1,000 &le; _n_ < 2,000).

During evolution, we limited program length to a maximum of 64 instructions and also
limited programs' maximum instruction-execution steps to 64.

For each test case, we evaluated programs on a pass-fail basis.

### Problem - For Loop Index

Programs receive three integer inputs _start_ (), _end_ (), and _step_ ().
Programs must output the following sequence:

_n_<sub>0</sub> = _start_ <br/>
_n_<sub>_i_</sub> = _n_<sub>_i_-1</sub> + _step_

for each _n_<sub>_i_</sub> < _end_.

In python:

```python
print([n for n in range(start, end, step)])
```

During evolution, we limited program length to a maximum of 128 instructions and also limited the
maximum number of instruction-execution steps to 256.

Program performance against a test case was measured on a gradient, using the Levenshtein
distance between the program's output and the correct output sequence.

### Problem - Compare String Lengths

Given three string inputs _s_<sub>1</sub>, _s_<sub>2</sub>, and _s_<sub>3</sub>,
programs must output _true_ if length(_s_<sub>1</sub>) < length(_s_<sub>2</sub>)
< length(_s_<sub>3</sub>) and _false_ otherwise.

During evolution, we limited program length to a maximum of 64 instructions and
also limited the maximum number of instruction-execution steps to 64.

We measured program performance against test cases on a pass-fail basis.

### Problem - Median

Programs are given three integer inputs (-100 &le; _input_<sub>_i_</sub> &le; 100)
and must output the median value.

During evolution, we limited program length to a maximum of 64 instructions and
also limited the maximum number of instruction-execution steps to 64.

We measured program performance against test cases on a pass-fail basis.

### Problem - Smallest

Programs are given four integer inputs (-100 &le; _input_<sub>_i_</sub> &le; 100)
and must output the smallest value.

During evolution, we limited program length to a maximum of 64 instructions and
also limited the maximum number of instruction-execution steps to 64.

We measured program performance against test cases on a pass-fail basis.

## Training/Testing Sets

During evolution, programs were assessed using a training set of test cases, which
defined the selection criteria used for lexicase selection. To qualify as a solution,
a program needed to perfectly pass all test cases in a separate testing set
(withheld generalization examples) in addition to passing all tests in the
training set used during evaluation. For all problems, we used an identical testing
set (1,000 tests) and the same input constraints as in (Helmuth and Spector, 2015).
In runs that used standard lexicase, we used identical training sets (100 tests)
as in (Helmuth and Spector, 2015), and in runs using down-sampled and cohort lexicase,
we always randomly subsampled from Helmuth and Spector's training sets.

The exact training and testing sets used can also be found in this repository:
[../data/prog-synth-examples/](../data/prog-synth-examples/).

## References

Helmuth, T., & Spector, L. (2015). General Program Synthesis Benchmark
Suite. In Proceedings of the 2015 on Genetic and Evolutionary Computation
Conference - GECCO ’15 (pp. 1039–1046). New York, New York, USA: ACM Press.
[https://doi.org/10.1145/2739480.2754769](https://doi.org/10.1145/2739480.2754769)

Forstenlechner, S., Fagan, D., Nicolau, M., & O’Neill, M. (2018). Towards Understanding and Refining the General Program Synthesis Benchmark Suite with Genetic Programming. In 2018 IEEE Congress on Evolutionary Computation (CEC) (pp. 1–6). IEEE. https://doi.org/10.1109/CEC.2018.8477953