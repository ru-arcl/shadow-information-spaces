"""Filters over shadow information spaces (Yu & LaValle, IEEE T-RO 2012)."""

from .events import *  # noqa: F401,F403
from .bipartite import BipartiteIState, Pseudo, vertex_name
from .fov import FovBatch, batch_fov, fov_to_component
from .lp import ShadowLP, lp_all_bounds, lp_bounds, lp_feasible
from .maxflow import FlowNetwork, InfeasibleError
from .nondeterministic import (CombinatorialFilter, MultiTeamFilter, RedOrBlueFilter, TwoPassBounds, Witness,
                               counting_filter, evader_status, pursuit_evasion_filter)
from .probabilistic import *  # noqa: F401,F403
from . import grid, polygon, rng  # noqa: F401

try:
    from importlib.metadata import PackageNotFoundError, version

    __version__ = version("shadowinfo")
except PackageNotFoundError:  # running from a source tree without installation
    __version__ = "0.1.0"

__all__ = [
    *events.__all__,  # noqa: F405
    "BipartiteIState", "Pseudo", "vertex_name", "FovBatch", "batch_fov", "fov_to_component",
    "ShadowLP", "lp_all_bounds", "lp_bounds", "lp_feasible", "FlowNetwork", "InfeasibleError",
    "CombinatorialFilter", "MultiTeamFilter", "RedOrBlueFilter", "TwoPassBounds", "Witness",
    "counting_filter", "evader_status", "pursuit_evasion_filter",
    *probabilistic.__all__,  # noqa: F405
    "grid", "polygon", "rng",
]
