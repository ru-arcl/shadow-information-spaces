"""ASCII grid maps and integer cell-to-cell visibility (DESIGN §5).

The grid world is a discrete stand-in for the polygonal environments of the
paper (T-RO 2012, Sec. II): ``#`` cells are obstacles, ``.`` cells are free.
A robot standing at a free cell sees a free cell ``(r, c)`` iff the segment
joining the two cell centres does not pass through the interior of an obstacle
cell (a segment slipping exactly through the shared corner of two diagonal
obstacle cells is blocked).  The test is an integer supercover traversal, so
the JS port (``docs/js/grid.js``) agrees bit for bit.

Cells are addressed either as ``(r, c)`` or by the row-major index
``r * cols + c``; masks and label grids are flat lists in row-major order.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

Cell = Tuple[int, int]

MAPS_DIR = Path(__file__).resolve().parent / "maps"  # shipped as package data


class GridMap:
    """Static occupancy grid with cached visibility masks."""

    def __init__(self, rows: Sequence[str]) -> None:
        rows = [r.rstrip() for r in rows if r.strip()]  # same trailing-whitespace rule as JS
        if not rows:
            raise ValueError("empty map")
        self.rows = len(rows)
        self.cols = max(len(r) for r in rows)
        self.blocked: List[int] = []
        for line in rows:
            line = line.ljust(self.cols, "#")
            for ch in line:
                if ch not in "#.":
                    raise ValueError(f"bad map character {ch!r}")
                self.blocked.append(1 if ch == "#" else 0)
        self.free_cells: List[int] = [i for i, b in enumerate(self.blocked) if not b]
        self._vis_cache: Dict[Tuple[int, Optional[int]], List[int]] = {}

    @classmethod
    def from_text(cls, text: str) -> "GridMap":
        return cls(text.splitlines())

    @classmethod
    def load(cls, path: Union[str, Path]) -> "GridMap":
        return cls.from_text(Path(path).read_text())

    def to_text(self) -> str:
        return "\n".join(
            "".join("#" if self.blocked[r * self.cols + c] else "." for c in range(self.cols))
            for r in range(self.rows)
        ) + "\n"

    def index(self, r: int, c: int) -> int:
        return r * self.cols + c

    def cell(self, i: int) -> Cell:
        return divmod(i, self.cols)

    def is_free(self, r: int, c: int) -> bool:
        return 0 <= r < self.rows and 0 <= c < self.cols and not self.blocked[r * self.cols + c]

    # -- visibility -------------------------------------------------------
    def line_of_sight(self, r0: int, c0: int, r1: int, c1: int) -> bool:
        """Integer supercover walk from the centre of ``(r0, c0)`` to ``(r1, c1)``.

        At each step the next boundary crossed is decided by comparing
        ``(1 + 2 ix) * ny`` with ``(1 + 2 iy) * nx``; equality means the
        segment passes exactly through a cell corner, which is blocked only
        when both side cells are obstacles.
        """
        nx, ny = abs(c1 - c0), abs(r1 - r0)
        sx = 1 if c1 > c0 else -1
        sy = 1 if r1 > r0 else -1
        cols, blocked = self.cols, self.blocked
        r, c, ix, iy = r0, c0, 0, 0
        while ix < nx or iy < ny:
            decision = (1 + 2 * ix) * ny - (1 + 2 * iy) * nx
            if decision == 0:
                if blocked[r * cols + c + sx] and blocked[(r + sy) * cols + c]:
                    return False
                r += sy
                c += sx
                ix += 1
                iy += 1
            elif decision < 0:
                c += sx
                ix += 1
            else:
                r += sy
                iy += 1
            if blocked[r * cols + c]:
                return False
        return True

    def visibility(self, robot: Cell, radius: Optional[int] = None) -> List[int]:
        """Flat 0/1 mask of free cells visible from ``robot`` (cached).

        ``radius`` (in cells) limits sensing to cells whose centre lies within
        ``radius + 1/2`` of the robot's, i.e. ``dr**2 + dc**2 <= radius**2 + radius``
        (a round integer disc; plain ``<= radius**2`` leaves one-cell nubs at
        the axis extremes that spawn spurious one-cell shadows along walls).
        """
        key = (self.index(*robot), radius)
        mask = self._vis_cache.get(key)
        if mask is None:
            r0, c0 = robot
            if not self.is_free(r0, c0):
                raise ValueError(f"robot cell {robot} is not free")
            r2 = None if radius is None else radius * radius + radius
            mask = [0] * (self.rows * self.cols)
            for i in self.free_cells:
                r, c = divmod(i, self.cols)
                if r2 is not None and (r - r0) ** 2 + (c - c0) ** 2 > r2:
                    continue
                if self.line_of_sight(r0, c0, r, c):
                    mask[i] = 1
            self._vis_cache[key] = mask
        return mask

    def shadow_mask(self, robot: Cell, radius: Optional[int] = None) -> List[int]:
        """Flat 0/1 mask of the shadow region: free and not visible."""
        vis = self.visibility(robot, radius)
        return [0 if (b or v) else 1 for b, v in zip(self.blocked, vis)]


def expand_path(waypoints: Sequence[Sequence[int]], loop: bool = True) -> List[Cell]:
    """Expand axis-aligned waypoints into a per-cell robot path.

    With ``loop`` the last waypoint is joined back to the first and the
    returned list does not repeat the starting cell, so one lap takes exactly
    ``len(path)`` unit moves.
    """
    pts = [(int(p[0]), int(p[1])) for p in waypoints]
    if loop and len(pts) > 1 and pts[-1] != pts[0]:
        pts.append(pts[0])
    path: List[Cell] = [pts[0]]
    for (r0, c0), (r1, c1) in zip(pts, pts[1:]):
        if r0 != r1 and c0 != c1:
            raise ValueError(f"segment {(r0, c0)}->{(r1, c1)} is not axis aligned")
        dr = (r1 > r0) - (r1 < r0)
        dc = (c1 > c0) - (c1 < c0)
        r, c = r0, c0
        while (r, c) != (r1, c1):
            r += dr
            c += dc
            path.append((r, c))
    if loop and len(path) > 1 and path[-1] == path[0]:
        path.pop()
    return path


def load_path(path: Union[str, Path]) -> List[Cell]:
    """Load ``{"waypoints": [[r, c], ...], "loop": true}`` and expand it."""
    d = json.loads(Path(path).read_text())
    return expand_path(d["waypoints"], d.get("loop", True))


def load_office() -> Tuple[GridMap, List[Cell]]:
    """The fixed demo map ``shadowinfo/grid/maps/office.txt`` and its robot path."""
    return GridMap.load(MAPS_DIR / "office.txt"), load_path(MAPS_DIR / "office_path.json")
