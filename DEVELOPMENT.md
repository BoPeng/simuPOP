# simuPOP Development Guide

This document contains notes for developers working on simuPOP.

## Build System Overview

simuPOP uses a modern Python build system:
- **Build backend**: scikit-build-core (PEP 517/518 compliant)
- **Build tool**: CMake 3.20+
- **SWIG**: 4.0+ for Python/C++ bindings

## Dependencies

### Vendored Dependencies

simuPOP vendors the following libraries to simplify building:

#### Boost (src/boost)

A minimal subset of Boost 1.86.0 is vendored in `src/boost/`. This includes only the components actually used by simuPOP:

- **boost::serialization** - Object serialization for saving/loading populations
- **boost::iostreams** - Compressed file I/O (gzip support)
- **boost::regex** - Regular expression support
- **boost::format** - String formatting
- **boost::numeric/ublas** - Linear algebra operations
- **boost::tuple**, **boost::lambda**, **boost::lexical_cast** - Utility libraries

**Why vendor Boost?**

1. **Self-contained builds**: `python -m build` works without network access or `--no-isolation`
2. **Reproducible builds**: Consistent Boost version across all platforms
3. **Smaller footprint**: 25MB subset vs 180MB full Boost distribution
4. **No system dependency**: Users don't need to install Boost separately

**Updating the vendored Boost:**

If you need to update Boost or add new Boost components:

1. Download the full Boost distribution:
   ```bash
   curl -L -O https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz
   tar -xzf boost_1_86_0.tar.gz
   ```

2. Install the `bcp` tool (Boost copy):
   ```bash
   # macOS
   brew install boost-bcp

   # Or build from Boost source
   cd boost_1_86_0/tools/bcp
   b2
   ```

3. Extract only the needed components by scanning source files:
   ```bash
   rm -rf src/boost
   mkdir -p src/boost
   bcp --boost=./boost_1_86_0 --scan \
       src/simuPOP/*.cpp src/simuPOP/*.h src/simuPOP/*.hpp \
       src/boost
   ```

4. Clean up unnecessary files:
   ```bash
   cd src/boost
   rm -rf doc boost.css boost.png rst.css Jamroot tools
   find libs -type d \( -name "test" -o -name "doc" -o -name "example" \) -exec rm -rf {} +
   ```

5. Verify the build works:
   ```bash
   python -m build
   ```

#### GSL (gsl/)

A subset of the GNU Scientific Library is vendored in `gsl/` for random number generation and statistical functions.

## Building simuPOP

### Standard Build

```bash
# Build and install
pip install .

# Build sdist and wheel
python -m build

# Editable install for development
pip install -e .
```

### Build Options

Pass CMake options via pip:

```bash
# Disable OpenMP
pip install . --config-settings=cmake.define.USE_OPENMP=OFF

# Verbose build
pip install -v .
```

### Module Variants

simuPOP builds multiple module variants for different use cases:

| Module | Description |
|--------|-------------|
| `std` | Standard module with debug checks |
| `op` | Optimized module (NDEBUG, no debug checks) |
| `la` | Long allele support |
| `laop` | Long allele, optimized |
| `ba` | Binary allele support |
| `baop` | Binary allele, optimized |
| `mu` | Mutant allele support |
| `muop` | Mutant allele, optimized |
| `lin` | Lineage tracking |
| `linop` | Lineage tracking, optimized |

## SWIG Bindings

The Python/C++ bindings are generated using SWIG. Pre-generated wrapper files are included in the repository (`src/simuPOP/*_wrap.cpp`) to avoid SWIG version compatibility issues during builds.

To regenerate SWIG wrappers (if needed):

```bash
cd src/simuPOP
swig -python -c++ -o simuPOP_std_wrap.cpp simuPOP_std.i
# Repeat for other variants
```

## Testing

```bash
# Install test dependencies
pip install pytest

# Run tests
pytest test/
```

## Release Process

1. Update version in `src/simuPOP/_version.py`
2. Build and test:
   ```bash
   python -m build
   pip install dist/*.whl
   pytest test/
   ```
3. Upload to PyPI:
   ```bash
   twine upload dist/*
   ```

## Platform-Specific Notes

### macOS

- OpenMP requires `libomp` from Homebrew:
  ```bash
  brew install libomp
  export OpenMP_ROOT=$(brew --prefix libomp)
  pip install .
  ```

### Windows

- Requires Visual Studio 2022 with C++ build tools
- SWIG can be installed via Chocolatey: `choco install swig`

### Linux

- Install build dependencies:
  ```bash
  sudo apt-get install zlib1g-dev cmake swig ninja-build
  ```

## Code Style

- C++ code follows the existing style in the codebase
- Python code should be compatible with Python 3.9+

## Troubleshooting

### Build fails with missing Boost headers

Ensure `src/boost/` exists and contains the vendored Boost subset. If headers are missing, re-extract using the bcp procedure above.

### OpenMP not found on macOS

Install libomp and set the OpenMP_ROOT environment variable:
```bash
brew install libomp
export OpenMP_ROOT=$(brew --prefix libomp)
```

### SWIG version mismatch

The repository includes pre-generated SWIG wrappers. If you encounter SWIG-related issues, ensure you're using SWIG 4.0 or later, or use the pre-generated wrappers.
