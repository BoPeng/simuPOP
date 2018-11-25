simuPOP is a general-purpose individual-based forward-time population genetics
simulation environment. Please refer to the simuPOP homepage

    http://simupop.sourceforge.net

for details.

## Installation

simuPOP is distributed under a GPL3 license. Starting from simuPOP 1.1.8, it supports only Python 3 (3.5 and up) on windows, mac and Linux systems.

simuPOP is part of the [conda-forge](https://conda-forge.github.io/), if you use [Anaconda Python 3](https://www.continuum.io/downloads),
you can install simuPOP with command

```
conda config --add channels conda-forge
conda install simuPOP
```

The first command only needs to be run only so you can upgrade simuPOP later using command

```
conda install simuPOP --upgrade
```

If you would like to use simuPOP with Python 2.5 - 2.7, please compile simuPOP from source, using either [simuPOP 1.1.7](https://pypi.python.org/pypi/simuPOP/1.1.7), 
or the [Python 2.x branch of simuPOP](https://github.com/BoPeng/simuPOP/tree/python2). Note that features that has been marked deprecated 
(e.g. `simuOpt.Param`, `simuPOP.plotting`) in simuPOP 1.1.7 and earlier are removed in simuPOP 1.1.8+, so simuPOP 1.1.7 would be your best
bet for the execution of legendary simuPOP scripts.

## Change Log since 1.1.7

### simuPOP 1.1.9
* [#28](https://github.com/BoPeng/simuPOP/issues/28): Fix a segmentation fault when providing a non-existent VSP index
* [#31](https://github.com/BoPeng/simuPOP/issues/31): Allow operator DiscardIf to accept a probability in addition to `True`/`False`
* [#35](https://github.com/BoPeng/simuPOP/issues/35): Add a `weightBy` parameter to allow `HeteroMating` to produce offspring subpopulation with weights determined by not only the size of the parental subpopulation, but also by for example number of mating pairs.
* [#49](https://github.com/BoPeng/simuPOP/issues/49): Fix output of loci positions in MS and other formats.
