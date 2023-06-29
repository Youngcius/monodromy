"""monodromy/__init__.py.

Top-level imports for `monodromy`.
"""

import monodromy.backend
import monodromy.coverage
import monodromy.polytopes
from monodromy.backend.lrs import LRSBackend

monodromy.backend.backend = LRSBackend()
