"""Loading the original polygons and robot paths; Java constants and number formatting.

Port of the loading side of the original Java implementation (``drawable.Polygon.read`` /
``Polygon.getPolygon``, the ``ProjectPanel*`` path literals and the ``Algorithm`` constants).
See ``docs/notes/original_java.md`` §2.2 and §3 and ``docs/notes/original_maps.md``.

The 14 maps ``data/{1..14}.dat`` are byte-identical copies of ``source/polygons/N.dat`` of
https://github.com/arc-l/shadow-information-space (BSD-3-Clause, Rutgers ARC group).  A file holds
one vertex ``x y`` per line in the y-up, counter-clockwise ``.dat`` frame; no transform is ever
applied (the Java code only flips y when drawing).  ``data/paths.json`` holds every robot path of
the Java panels plus the paths digitised from the paper figures, each with its provenance.
"""

from __future__ import annotations

import functools
import json
import math
from decimal import Decimal
from pathlib import Path as _FsPath
from typing import Dict, Iterable, List, Sequence, Tuple, Union

DATA_DIR = _FsPath(__file__).resolve().parent / "data"

EPSILON = 5e-5
"""``Algorithm.epsilon`` (= ``0.00005``): absolute tolerance, in map units, of every geometric test."""

PERTURB = 0.005
"""``Algorithm.purturb``: step used to sample just after a critical point (``purturbPointAlongSeg``)."""

Y_MAX = 1000.0
"""``scaling-factor`` of ``config/geometry.properties``: the maps live in ``[0, 1000]^2``."""

N_POLYGONS = 14

__all__ = [
    "DATA_DIR", "EPSILON", "PERTURB", "Y_MAX", "N_POLYGONS", "PANELS", "parse_dat", "read_dat", "format_dat",
    "write_dat", "NUMBER_OF_POLYGON_DATA_FILES", "DATA_FILE_POSTFIX", "random_polygon_from_data", "random_polygon",
    "load_polygon", "load_paths", "path_record", "load_path", "path_labels", "demo_paths", "demo_path_record",
    "load_demo_path", "java_double_str",
]


def parse_dat(text: str) -> List[Tuple[float, float]]:
    """Parse the text of a ``.dat`` polygon file into ``[(x, y), ...]``.

    ``Polygon.read``: each line is ``line.trim().split(" ")`` and the first two fields are parsed
    with ``Double.parseDouble`` (Python ``float`` parses decimal text to the same double).
    Fix of quirk B16: blank lines are skipped (Java throws ``NumberFormatException`` on them), and
    runs of whitespace separate fields (Java's ``split(" ")`` fails on two spaces).
    """
    pts: List[Tuple[float, float]] = []
    for lineno, line in enumerate(text.splitlines(), 1):
        f = line.split()
        if not f:
            continue
        if len(f) < 2:
            raise ValueError(f"line {lineno}: expected 'x y', got {line!r}")
        pts.append((float(f[0]), float(f[1])))
    return pts


def read_dat(path: Union[str, _FsPath]) -> List[Tuple[float, float]]:
    """Read a ``.dat`` polygon file (see :func:`parse_dat`)."""
    return parse_dat(_FsPath(path).read_text())


def format_dat(vertices: Iterable[Sequence[float]]) -> str:
    """``Polygon.write``: ``"x y\\n"`` per vertex (Java ``Double.toString``), then one blank line.

    Java's own ``read`` throws on that trailing blank line (B16); :func:`parse_dat` skips it, so
    ``parse_dat(format_dat(v)) == v`` for every finite vertex list.  Accepts a
    :class:`~shadowinfo.polygon.geometry.Polygon` (its ``vertices``) or a list of points.
    """
    vs = getattr(vertices, "vertices", vertices)
    return "".join(f"{java_double_str(float(v[0]))} {java_double_str(float(v[1]))}\n" for v in vs) + "\n"


def write_dat(vertices: Iterable[Sequence[float]], path: Union[str, _FsPath]) -> None:
    """Write a ``.dat`` polygon file in ``Polygon.write`` format (see :func:`format_dat`)."""
    _FsPath(path).write_text(format_dat(vertices))


NUMBER_OF_POLYGON_DATA_FILES = 7
"""``number`` of ``polygons/polygons.properties`` (``Polygon.getNumberOfPolygonDataFiles``): only maps
1..7 can be drawn by :func:`random_polygon_from_data`; 8..14 are reachable by number only."""

DATA_FILE_POSTFIX = ".dat"
"""``postfix`` of ``polygons/polygons.properties`` (``Polygon.getDataFilePostfix``)."""


def _java_random(rng):
    from .oracle import JavaRandom

    return rng if isinstance(rng, JavaRandom) else JavaRandom(rng)


def random_polygon_from_data(rng=None, number: int = NUMBER_OF_POLYGON_DATA_FILES):
    """``Polygon.randomPolygonFromData``: map ``(int) (Math.random() * number + 1)``, i.e. one of 1..``number``.

    ``rng`` is a :class:`~shadowinfo.polygon.oracle.JavaRandom` or a seed for one (None: unseeded, like
    ``Math.random()``); one ``next_double`` is drawn, as in Java.
    """
    r = _java_random(rng)
    return load_polygon(int(r.next_double() * number + 1))


def random_polygon(num_vertices: int, width: float = Y_MAX, height: float = Y_MAX, rng=None):
    """``Polygon.randomPolygon``: ``num_vertices`` points uniform in ``[0, width) x [0, height)``, x then y
    drawn per vertex from ``Math.random()`` (``rng`` as in :func:`random_polygon_from_data`).

    Like the Java (which never calls it) the result is generally **not** a simple polygon, so it is built
    with ``validate=False``; the cut and visibility functions assume a simple counter-clockwise polygon.
    """
    from .geometry import Polygon

    r = _java_random(rng)
    pts = []
    for _ in range(num_vertices):
        x = width * r.next_double()
        pts.append((x, height * r.next_double()))
    return Polygon(pts, name="random", validate=False)


@functools.lru_cache(maxsize=None)
def _load_numbered(n: int):
    from .geometry import Polygon

    if not 1 <= n <= N_POLYGONS:
        raise ValueError(f"polygon number must be 1..{N_POLYGONS}, got {n}")
    return Polygon(read_dat(DATA_DIR / f"{n}.dat"), name=str(n))


def load_polygon(n_or_path: Union[int, str, _FsPath]):
    """``Polygon.getPolygon(n, dc)``: load ``data/<n>.dat`` (``n`` in 1..14) or a ``.dat`` file path.

    Numbered polygons are cached, so repeated loads return the same (immutable) object and share its
    cached cuts.  The vertices are used exactly as read (y-up, counter-clockwise).
    """
    from .geometry import Polygon

    if isinstance(n_or_path, int) and not isinstance(n_or_path, bool):
        return _load_numbered(n_or_path)
    p = _FsPath(n_or_path)
    return Polygon(read_dat(p), name=p.stem)


@functools.lru_cache(maxsize=1)
def _paths_json() -> dict:
    return json.loads((DATA_DIR / "paths.json").read_text())


def load_paths() -> Dict[str, dict]:
    """All robot paths, keyed by label (``P12``, ``P1``, ``P5``, ``P14a``, ``P14b``, ``P1s``, ``P13``,
    ``P13+``, ``P13+tail``, ``P5s``, ``P1L`` and the digitised ``fig_*`` paths).

    Each record has ``polygon``, ``waypoints`` (``.dat`` frame), ``kind`` (``code``/``digitised``),
    ``status``, ``panels``, ``figure``, ``provenance`` and ``harness`` (the golden-fixture run).
    Returns a fresh copy.
    """
    return json.loads(json.dumps(_paths_json()["paths"]))


PANELS: Dict[str, str] = {"ProjectPanel": "P12", "ProjectPanel2": "P1", "ProjectPanel3": "P5",
                          "ProjectPanel4": "P14a", "ProjectPanel5": "P14b"}
"""The path each ``ProjectPanel*`` class runs (original_java.md §3)."""


def path_labels() -> List[str]:
    """Labels of all stored paths, in file order."""
    return list(_paths_json()["paths"])


def path_record(label: str) -> dict:
    """The ``paths.json`` record of one path (accepts a label, a panel class name or a harness name
    such as ``PP5@14``)."""
    paths = _paths_json()["paths"]
    if label in PANELS:
        label = PANELS[label]
    if label not in paths and "@" in label:
        name, poly = label.split("@")
        for k, rec in paths.items():
            if rec["harness"]["name"] == name and rec["polygon"] == int(poly):
                label = k
                break
    if label not in paths:
        raise KeyError(f"unknown path {label!r}; known: {', '.join(paths)}")
    return json.loads(json.dumps(paths[label]))


def load_path(label: str):
    """``(polygon, path)`` for a stored path label (see :func:`path_record`)."""
    from .geometry import Path

    rec = path_record(label)
    return load_polygon(rec["polygon"]), Path([tuple(p) for p in rec["waypoints"]])


def demo_paths() -> Dict[str, dict]:
    """The default demo path of every map, keyed by map number (``"1"`` .. ``"14"``).

    Section ``demo`` of ``data/paths.json``.  Each record has ``polygon``, ``waypoints`` (``.dat``
    frame), ``at_end`` (what :class:`~shadowinfo.polygon.simulate.PolygonSimulator` does at the end of
    the path: ``"loop"`` for a closed patrol, ``"reverse"`` to run an open path back and forth),
    ``source`` (the stored path it copies, or ``None`` for a path designed for the demo),
    ``provenance`` (``"demo"``) and ``note``.  Returns a fresh copy.
    """
    return {k: v for k, v in json.loads(json.dumps(_paths_json()["demo"])).items() if not k.startswith("_")}


def demo_path_record(n: int) -> dict:
    """The demo-path record of map ``n`` (see :func:`demo_paths`)."""
    recs = demo_paths()
    if str(n) not in recs:
        raise KeyError(f"no demo path for polygon {n!r}")
    return recs[str(n)]


def load_demo_path(n: int):
    """``(path, at_end)``: the default demo path of map ``n`` as a :class:`~shadowinfo.polygon.geometry.Path`."""
    from .geometry import Path

    rec = demo_path_record(n)
    return Path([tuple(p) for p in rec["waypoints"]]), rec["at_end"]


def java_double_str(x: float) -> str:
    """Java ``Double.toString(x)`` layout of the shortest round-trip digits of ``x``.

    ``1e-3 <= |x| < 1e7`` prints as plain decimal with at least one fractional digit (``0.0``,
    ``100.0``, ``0.05623``); other values use ``d.dddE<exp>`` (``1.0E-4``, ``1.2345E7``).  The digits
    are Python's shortest repr, which is what JDK >= 19 prints.  JDK <= 18 (the fixtures were made with
    17) occasionally printed one or two extra digits (JDK-4511638), and subnormals may differ in the
    second digit (Java prints ``Double.MIN_VALUE`` as ``4.9E-324``); none of the golden outputs is
    affected.
    """
    if math.isnan(x):
        return "NaN"
    if math.isinf(x):
        return "Infinity" if x > 0 else "-Infinity"
    sign = "-" if math.copysign(1.0, x) < 0 else ""
    ax = abs(x)
    if ax == 0.0:
        return sign + "0.0"
    _, digits, exp = Decimal(repr(ax)).as_tuple()
    ds = list(digits)
    while len(ds) > 1 and ds[-1] == 0:
        ds.pop()
        exp += 1
    s = "".join(map(str, ds))
    e10 = len(s) - 1 + exp  # value = s[0].s[1:] * 10**e10
    if 1e-3 <= ax < 1e7:
        if e10 >= 0:
            ip = s[:e10 + 1].ljust(e10 + 1, "0")
            fp = s[e10 + 1:] or "0"
        else:
            ip, fp = "0", "0" * (-e10 - 1) + s
        return f"{sign}{ip}.{fp}"
    return f"{sign}{s[0]}.{s[1:] or '0'}E{e10}"
