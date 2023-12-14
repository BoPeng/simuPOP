simuPOP is a general-purpose individual-based forward-time population genetics
simulation environment. Please refer to the simuPOP homepage http://simupop.sourceforge.net for details.

## Installation

simuPOP is distributed under a GPL3 license. Starting from simuPOP 1.1.8, it supports only Python 3 (3.5 and up) on windows, mac and Linux systems.

simuPOP is part of the [conda-forge](https://conda-forge.github.io/), if you use [Anaconda Python 3](https://www.continuum.io/downloads), you can install simuPOP with command

```
conda install -c conda-forge simuPOP
```

If you would like to use simuPOP with Python 2.5 - 2.7, please compile simuPOP from source, using either [simuPOP 1.1.7](https://pypi.python.org/pypi/simuPOP/1.1.7),
or the [Python 2.x branch of simuPOP](https://github.com/BoPeng/simuPOP/tree/python2). Note that features that has been marked deprecated
(e.g. `simuOpt.Param`, `simuPOP.plotting`) in simuPOP 1.1.7 and earlier are removed in simuPOP 1.1.8+, so simuPOP 1.1.7 would be your best
bet for the execution of legendary simuPOP scripts.

## Documenation

The user guide and reference manual of simuPOP is available at http://bopeng.github.io/simuPOP/. simuPOP is also introduced in the following two books:

1. **Forward-Time Population Genetics Simulations: Methods, Implementation, and Applications** by Bo Peng, Marek Kimmel and Christopher I Amos, published by Wiley & Sons Inc, and available at [Amazon](http://www.amazon.com/gp/product/0470503483/?tag=wwwwileycom-20) and other bookstores.
2. **Bioinformatics with Python Cookbook** by Tiago Antao, available at [Amazon](https://www.amazon.com/Bioinformatics-Python-Cookbook-Tiago-Antao/dp/1782175113)

## Change Log since 1.1.7

### simuPOP 1.1.12
* [#102](https://github.com/BoPeng/simuPOP/issues/102): Fix compatibility issues with Python 3.11
* [#98](https://github.com/BoPeng/simuPOP/issues/98): Write banner messaged to stderr insted of stdout

### simuPOP 1.1.11
* [#94](https://github.com/BoPeng/simuPOP/issues/94): Fix `HeteroMating` when being used in a `ConditionalMating` mating scheme.
* [#93](https://github.com/BoPeng/simuPOP/issues/93): Add function form of selection operators such as `maSelect` (for `MapSelector`) and `mlSelect` (for `MlSelector`).


### simuPOP 1.1.10
* [#70](https://github.com/BoPeng/simuPOP/issues/70): Fix compatibility with Mac OS using `libc++` instead of `libstdc++`.

### simuPOP 1.1.9
* [#28](https://github.com/BoPeng/simuPOP/issues/28): Fix a segmentation fault when providing a non-existent VSP index
* [#31](https://github.com/BoPeng/simuPOP/issues/31): Allow operator DiscardIf to accept a probability in addition to `True`/`False`
* [#35](https://github.com/BoPeng/simuPOP/issues/35): Add a `weightBy` parameter to allow `HeteroMating` to produce offspring subpopulation with weights determined by not only the size of the parental subpopulation, but also by for example number of mating pairs.
* [#49](https://github.com/BoPeng/simuPOP/issues/49): Fix output of loci positions in MS and other formats.
