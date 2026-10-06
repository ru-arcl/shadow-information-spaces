"""Polygon front end: a port of the original Java implementation (2008-2012, no applets).

The original program (https://github.com/arc-l/shadow-information-space) drives one omnidirectional
robot along a path inside a simple polygon, finds where the path crosses the polygon's critical lines
(inflection, bitangent and single-tangent rays), samples the visible gaps just after every crossing
and turns them into shadow component events for the combinatorial filter (T-RO 2012 Secs. II, V).
See ``docs/notes/original_java.md`` (inventory, quirks B1-B18, porting plan §5) and
``docs/notes/original_maps.md`` (which map and path produced which paper figure).

Modules (``docs/notes/original_java.md`` §5.1):

* :mod:`.io` -- the 14 original maps and every robot path (``data/``), Java constants and formatting;
* :mod:`.geometry` -- ``Algorithm`` primitives, ``Polygon``, ``Path`` and the ``java.awt.geom``
  semantics they rely on, bit for bit;
* :mod:`.cuts` -- inflection, bitangent and single-tangent cuts and the path's critical points;
* :mod:`.visibility` -- physical gaps and visibility polygons;
* :mod:`.gaps` -- gap tracking along a path (``getGaps``): shadow IDs, links and times;
* :mod:`.events` -- gap histories as :mod:`shadowinfo.events` component events, never-see-evader;
* :mod:`.oracle` -- ``java.util.Random`` and the single-type-agent oracle (simulated targets);
* :mod:`.stagent` -- the Java filter: equations, bipartite I-state, Java's max-flow bounds, exact bounds;
* :mod:`.simulate` -- new (no Java counterpart): a robot on a path plus random-walk targets, shadow pockets,
  component and FOV events and ground truth, the polygon counterpart of :mod:`shadowinfo.grid.simulate`;
* ``python -m shadowinfo.polygon`` (:mod:`.__main__`) -- what the ``ProjectPanel*`` classes print, plus
  ``list``, ``cuts N`` and ``vis N x y``.

A full run, as ``ProjectPanel5`` (T-RO Fig. 15(b)) but with the bugs fixed::

    poly, path = load_path("P14b")
    history = get_gaps(poly, path)                                   # 386 gap sets
    oracle = SingleTypeAgentOracle(1_000_000, rng=1).initialize(history)
    exact_bounds(oracle)                                             # {shadow: (lo, hi)} at the end
    oracle.sequence()                                                # observations for any filter

All coordinates are the y-up, counter-clockwise ``.dat`` frame; nothing is ever flipped.  Where the
Java code has a bug, functions take ``compat="java"`` (reproduce it) or ``compat="fixed"`` (default).
"""

from .cuts import (Bitangent, CriticalPoint, Cut, CutType, GeneralInflection, InflectionType, all_critical_points,
                   all_cuts, bitangent_cuts, bitangent_lines, cut_intersect_points, general_inflection_cuts,
                   inflections, single_tangent_cuts)
from .geometry import (COMPAT_MODES, GeneralPath, GeometryError, Path, Polygon, purturb_point_along_seg,
                       seg_in_polygon)
from .io import (EPSILON, PANELS, PERTURB, demo_path_record, demo_paths, java_double_str, load_demo_path, load_path,
                 load_paths, load_polygon, path_labels, path_record)
from .visibility import PhysicalGap, physical_gaps, vertex_visibility, visibility_polygon
from .gaps import (Gap, GapHistory, GapTrackingError, PathSamples, chain_graph, format_gap_set, format_gap_sets,
                   get_gaps, java_hashset_order, sample_path)
from .events import (GapEvents, NSEState, gap_history_to_events, never_see_evader, never_see_evader_via_filter,
                     transition_events)
from .simulate import (AT_END, PolygonSimulator, Shadow, pocket_polygon, segment_in_polygon, shadow_pockets, simulate,
                       validate_path)
from .oracle import EventType, JavaRandom, OracleError, SingleTypeAgentEvent, SingleTypeAgentOracle
from .stagent import (Equation, JavaBipartiteGraph, JavaBounds, bounds_lines, derive_gap_evolving_equations,
                      derive_shadow_info_state, exact_bounds, java_bounds, project_panel_lines)

__all__ = [
    "Bitangent", "CriticalPoint", "Cut", "CutType", "GeneralInflection", "InflectionType", "all_critical_points",
    "all_cuts", "cut_intersect_points", "bitangent_cuts", "bitangent_lines", "general_inflection_cuts", "inflections",
    "single_tangent_cuts", "COMPAT_MODES", "GeneralPath", "GeometryError", "Path", "Polygon",
    "purturb_point_along_seg", "seg_in_polygon", "EPSILON", "PANELS", "PERTURB", "java_double_str", "load_path",
    "load_paths", "load_polygon", "path_labels", "path_record", "PhysicalGap", "physical_gaps",
    "vertex_visibility", "visibility_polygon", "Gap", "GapHistory", "GapTrackingError", "PathSamples",
    "chain_graph", "format_gap_set", "format_gap_sets", "get_gaps", "java_hashset_order", "sample_path",
    "GapEvents", "NSEState", "gap_history_to_events", "never_see_evader", "never_see_evader_via_filter",
    "transition_events", "EventType", "JavaRandom", "OracleError", "SingleTypeAgentEvent", "SingleTypeAgentOracle",
    "Equation", "JavaBipartiteGraph", "JavaBounds", "bounds_lines", "derive_gap_evolving_equations",
    "derive_shadow_info_state", "exact_bounds", "java_bounds", "project_panel_lines", "demo_path_record",
    "demo_paths", "load_demo_path", "AT_END", "PolygonSimulator", "Shadow", "pocket_polygon", "segment_in_polygon",
    "shadow_pockets", "simulate", "validate_path",
]
