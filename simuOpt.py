#!/usr/bin/env python
"""
Backward compatibility stub for simuOpt module.

DEPRECATED: Please use 'from simuPOP.simuOpt import setOptions' instead.
            This import path will be removed in a future version.
"""

import warnings

warnings.warn(
    "Importing 'simuOpt' directly is deprecated. "
    "Please use 'from simuPOP.simuOpt import setOptions' instead. "
    "This import path will be removed in a future version.",
    DeprecationWarning,
    stacklevel=2
)

from simuPOP.simuOpt import *
