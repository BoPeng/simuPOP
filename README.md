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

