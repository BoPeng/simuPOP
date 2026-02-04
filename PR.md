# Modernize Build System: Migrate from setup.py to CMake + scikit-build-core

## Summary

This PR modernizes the simuPOP build system by migrating from the legacy `setup.py` distutils-based build to a modern CMake + scikit-build-core system. This provides better cross-platform support, automatic compiler feature detection, and compatibility with modern Python packaging standards (PEP 517/518).

## Changes

### New Files

- **CMakeLists.txt**: Complete CMake build configuration
  - Automatic detection of compiler features (isfinite, isnan, isinf, etc.)
  - Automatic Boost download if not present
  - Builds all 10 module variants (std, op, la, laop, ba, baop, mu, muop, lin, linop)
  - OpenMP support detection
  - Uses pre-generated SWIG wrapper files to avoid SWIG version compatibility issues

- **cmake/SimuPOPConfig.h.in**: Template for auto-generated config header
  - Replaces manual platform detection in config.h
  - Uses CMake's `check_function_exists()` and `check_symbol_exists()` for portable feature detection

- **pyproject.toml**: Modern Python packaging configuration
  - Uses scikit-build-core as build backend
  - Dynamic version from `src/_version.py`
  - Proper metadata and classifiers

### Modified Files

- **src/simuPOP_cfg.h**: Updated to use auto-generated config header from CMake build
- **src/_version.py**: Bumped version to 1.1.18

## Build Instructions

```bash
# Using pip (recommended)
pip install .

# Or using python -m build
pip install build
python -m build

# Direct CMake build (for development)
cmake -B build -S .
cmake --build build -j4
```

## Key Features

1. **Platform-independent feature detection**: No more manual `#ifdef` blocks for different platforms
2. **Modern Python packaging**: Compatible with pip, build, and other PEP 517 tools
3. **Pre-generated SWIG wrappers**: Avoids SWIG version compatibility issues
4. **Automatic Boost handling**: Downloads Boost if not present
5. **Single source of truth for version**: `src/_version.py` feeds both CMake and pyproject.toml

## Compatibility

- Python 3.9+
- CMake 3.20+
- C++17 compatible compiler
- macOS, Linux, Windows (with Visual Studio 2022)
