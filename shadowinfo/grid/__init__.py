"""Grid-world simulator producing shadow event streams (DESIGN §5).

Mirrored bit for bit by ``docs/js/grid.js``.
"""

from .gridmap import MAPS_DIR, GridMap, expand_path, load_office, load_path
from .shadows import ShadowTracker, Transition, label_components
from .simulate import MODES, MOVES, GridSimulator, initial_bounds, label_hash, simulate

__all__ = [
    "MAPS_DIR", "GridMap", "expand_path", "load_office", "load_path",
    "ShadowTracker", "Transition", "label_components",
    "MODES", "MOVES", "GridSimulator", "initial_bounds", "label_hash", "simulate",
]
