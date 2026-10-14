"""Feature modules. Importing this package registers every feature in REGISTRY.

Build order (Code_Architecture_Plan.md §7): rate -> ventricular -> atrial -> quality.
Each module is added only after its tests exist.
"""

from . import rate, ventricular  # noqa: F401
from . import atrial  # noqa: F401  (after ventricular: build order, column order)
