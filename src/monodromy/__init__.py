"""monodromy/__init__.py.

Top-level imports for `monodromy`.
"""

import monodromy.backend
import monodromy.coverage
import monodromy.polytopes
from monodromy.backend.lrs import LRSBackend

monodromy.backend.backend = LRSBackend()


import warnings

# 抑制 NumPy 中行列式计算的警告
warnings.filterwarnings("ignore", message="divide by zero encountered in det")
warnings.filterwarnings("ignore", message="invalid value encountered in det")
