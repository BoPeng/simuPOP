simuPOP is a general-purpose individual-based forward-time population genetics
simulation environment. Please refer to the simuPOP homepage [https://bopeng.github.io/simuPOP/](https://bopeng.github.io/simuPOP/) for details.

## Installation

simuPOP is distributed under a GPL3 license. Starting from simuPOP 1.1.8, it supports only Python 3 (3.9 and up) on Windows, macOS and Linux systems.

### Install from conda-forge

simuPOP is part of the [conda-forge](https://conda-forge.github.io/). If you use [Anaconda Python 3](https://www.continuum.io/downloads), you can install simuPOP with command

```
conda install -c conda-forge simuPOP
```

### Install from PyPI

```
pip install simuPOP
```

### Install from source

If you are working with a development version of simuPOP or an unsupported platform, you can install from source:

```
git clone https://github.com/BoPeng/simuPOP.git
cd simuPOP
pip install .
```

Please refer to [INSTALL](https://github.com/BoPeng/simuPOP/blob/master/INSTALL) for detailed instructions and platform-specific requirements.

## Documentation

The user guide and reference manual of simuPOP is available at http://bopeng.github.io/simuPOP/. simuPOP is also introduced in the following two books:

1. **Forward-Time Population Genetics Simulations: Methods, Implementation, and Applications** by Bo Peng, Marek Kimmel and Christopher I Amos, published by Wiley & Sons Inc, and available at [Amazon](http://www.amazon.com/gp/product/0470503483/?tag=wwwwileycom-20) and other bookstores.
2. **Bioinformatics with Python Cookbook** by Tiago Antao, available at [Amazon](https://www.amazon.com/Bioinformatics-Python-Cookbook-Tiago-Antao/dp/1782175113)

## Change Log since 1.1.7

### simuPOP 1.1.18
* Modernize build system: migrate from setup.py to CMake + scikit-build-core (PEP 517 compliant)
* Add multi-platform CI testing (Linux, macOS, Windows)
* Move `simuOpt` module into the simuPOP package. Use `from simuPOP.simuOpt import setOptions` instead of `from simuOpt import setOptions`. The old import path still works but is deprecated.

### simuPOP 1.1.17
* Fix compatibility with Python 3.13 and latest version of Xcode on mac.

### simuPOP 1.1.16
* [#103](https://github.com/BoPeng/simuPOP/issues/103) Fix compatibility with recent compilers with c++17 support, and fix conda-forge release.

### simuPOP 1.1.15
* [#122](https://github.com/BoPeng/simuPOP/issues/122) Fix compatibility with Visual Studio 2022 under windows.

### simuPOP 1.1.14
* [#114](https://github.com/BoPeng/simuPOP/issues/117) Fix compatibility with Python 3.11.

### simuPOP 1.1.13
* [#114](https://github.com/BoPeng/simuPOP/issues/114) Allow negative weight of a heterogeneous mating scheme to generate less individuals than calculated.

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
