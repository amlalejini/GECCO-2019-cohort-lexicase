# Cohorts enhance speed and solution quality of lexicase selection

This repository contains supplementary material (e.g., source code, data, documentation) 
for our 2019 GECCO submission, Cohorts Enhance Speed and Solution Quality of Lexicase
Selection.

**Navigation**

<!-- TOC -->

- [Project Overview](#project-overview)
  - [Contribution Authors](#contribution-authors)
- [Repository Guide](#repository-guide)

<!-- /TOC -->

## Project Overview

Here, we introduce cohort lexicase selection, and new variant on the standard lexicase
parent selection algorithm.

Lexicase selection has proven to be highly successful for finding solutions in genetic
programming, especially for test-based problems where there are many distinct test cases that must all be passed.
However, lexicase selection (like many other selection schemes) requires prospective
solutions to be evaluated against most test cases each generation, which can be
computationally demanding.
In this work, we reduce the number of required evaluations per generation by dividing
the population and set of test cases into cohorts.
Cohorts are randomly reassigned each generation, and candidate solutions are only ever
evaluated on test cases from within their own cohort, reducing the total number of
evaluations needed while ensuring that each lineage eventually encounters all test
cases.
We test a range of cohort numbers/sizes on five program synthesis problems.
We demonstrate that for a moderate number of cohorts, cohort lexicase does at least as well, and often better than standard lexicase, given the same number of evaluations.
Overall, this technique allows us to evolve populations for more generations, use larger population or test set sizes, or simply to find solutions using fewer total evaluations.

### Contribution Authors

- Jose Hernandez
- [Alexander Lalejini](lalejini.com)
- [Emily Dolson](emilyldolson.com)
- [Charles Ofria](ofria.com)

## Repository Guide

- [analysis/](https://github.com/amlalejini/GECCO-2019-cohort-lexicase/tree/master/analysis/)
  - Contains R scripts used for data analyses and generating graphs.
- [data/](https://github.com/amlalejini/GECCO-2019-cohort-lexicase/tree/master/data/)
  - Contains raw data for preliminary and published experiments as well as the
    training and testing examples used for the programming synthesis benchmark 
    problems (taken from [Tom Helmuth's example repository](https://github.com/thelmuth/Program-Synthesis-Benchmark-Data)).
- [docs/](https://github.com/amlalejini/GECCO-2019-cohort-lexicase/tree/master/docs/)
  - Contains miscellaneous documentation associated with this work.
- [experiment/](https://github.com/amlalejini/GECCO-2019-cohort-lexicase/tree/master/experiment/)
  - Contains the source code (C++) for our simple linear GP representation and for
    running the experiments discussed in our contribution.
- [hpcc/](https://github.com/amlalejini/GECCO-2019-cohort-lexicase/tree/master/hpcc/)
  - Contains scripts used to submit experiment jobs to MSU's HPCC.
- [media/](https://github.com/amlalejini/GECCO-2019-cohort-lexicase/tree/master/media/)
  - Contains media (images) associated with this work.
- [scripts/](https://github.com/amlalejini/GECCO-2019-cohort-lexicase/tree/master/scripts/)
  - Contains utility scripts used for managing experiments on the HPCC and for aggregating
    and manipulating experiment data.