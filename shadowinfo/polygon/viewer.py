"""Live visibility viewer: the robot follows the mouse (a stand-in for the original applet windows).

``python -m shadowinfo.polygon show [MAP]`` opens a matplotlib window on one of the 14 original maps.  As
the mouse moves inside the polygon, the visibility polygon :math:`V(q)` (yellow) and the shadow pockets
(T-RO 2012, Sec. II, Fig. 1) are recomputed with :func:`~shadowinfo.polygon.visibility.physical_gaps`,
:func:`~shadowinfo.polygon.visibility.visibility_polygon` and
:func:`~shadowinfo.polygon.simulate.pocket_polygon`.  Each pocket is coloured by the reflex vertex that casts
its gap, so a shadow keeps its colour while the robot moves.  A click pins/unpins the robot; the left/right
arrow keys change the map.

Requires matplotlib (``pip install -e '.[plot]'``).
"""

from __future__ import annotations

from typing import Optional, Tuple

from .io import N_POLYGONS, load_demo_path, load_polygon
from .simulate import pocket_polygon
from .visibility import PhysicalGap, physical_gaps, visibility_polygon

Point = Tuple[float, float]


def gap_anchor(g: PhysicalGap) -> int:
    """Index of the reflex vertex casting gap ``g`` (see :func:`pocket_polygon`); stable while the robot moves."""
    return g.start_edge if g.start_point is None else g.end_edge + 1


class Viewer:
    """Matplotlib view of map ``map_no`` with the robot at ``q`` (default: start of the map's demo path)."""

    def __init__(self, ax, map_no: int = 12, q: Optional[Point] = None):
        import matplotlib.pyplot as plt

        self.ax = ax
        self.colors = plt.get_cmap("tab20").colors
        self.pinned = False
        self.load(map_no, q)

    def load(self, map_no: int, q: Optional[Point] = None) -> None:
        from matplotlib.patches import Polygon as MplPolygon

        self.map_no = map_no
        self.poly = load_polygon(map_no)
        if q is None:
            q = tuple(load_demo_path(map_no)[0].points[0])
        elif not self.poly.contains(q):
            raise ValueError(f"point {q} is outside map {map_no}")
        ax = self.ax
        ax.clear()
        ax.add_patch(MplPolygon(self.poly.vertices, closed=True, fc="#d9dde2", ec="#222", lw=1.2, zorder=0))
        self.vis = MplPolygon([(0.0, 0.0)], closed=True, fc="#fff3b0", ec="#c9a800", lw=1, zorder=2)
        ax.add_patch(self.vis)
        self.pockets = []
        (self.robot,) = ax.plot([], [], "o", ms=9, mec="white", zorder=3)
        xs, ys = zip(*self.poly.vertices)
        ax.set_xlim(min(xs) - 20, max(xs) + 20)
        ax.set_ylim(min(ys) - 20, max(ys) + 20)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        self.update(q)

    def update(self, q: Point) -> None:
        from matplotlib.patches import Polygon as MplPolygon

        self.q = q
        for p in self.pockets:
            p.remove()
        gaps = physical_gaps(self.poly, q)
        self.pockets = []
        for g in gaps:
            patch = MplPolygon(pocket_polygon(self.poly, g), closed=True, ec="none", alpha=0.55, zorder=1,
                               fc=self.colors[gap_anchor(g) % len(self.colors)])
            self.ax.add_patch(patch)
            self.pockets.append(patch)
        self.vis.set_xy(visibility_polygon(self.poly, q))
        self.robot.set_data([q[0]], [q[1]])
        self.robot.set_color("#d33" if self.pinned else "#1f5fff")
        self.ax.set_title(f"Map {self.map_no}  ·  robot ({q[0]:.0f}, {q[1]:.0f})  ·  {len(gaps)} shadows"
                          f"{'  ·  pinned' if self.pinned else ''}\n"
                          "move the mouse · click to pin/unpin · ←/→ change map", fontsize=10)
        self.ax.figure.canvas.draw_idle()

    def on_move(self, ev) -> None:
        if self.pinned or ev.inaxes is not self.ax or ev.xdata is None:
            return
        p = (ev.xdata, ev.ydata)
        if self.poly.contains(p):
            self.update(p)

    def on_click(self, ev) -> None:
        if ev.inaxes is not self.ax or ev.xdata is None:
            return
        p = (ev.xdata, ev.ydata)
        self.pinned = not self.pinned
        self.update(p if self.poly.contains(p) else self.q)

    def on_key(self, ev) -> None:
        step = {"right": 1, "left": -1}.get(ev.key)
        if step:
            self.pinned = False
            self.load((self.map_no - 1 + step) % N_POLYGONS + 1)


def show(map_no: int = 12, q: Optional[Point] = None, save: Optional[str] = None) -> Viewer:
    """Open the live viewer (or, with ``save``, write one frame to that image file and return)."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.5, 7.8))
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.92)
    v = Viewer(ax, map_no, q)
    if save:
        fig.savefig(save, dpi=110)
        plt.close(fig)
        return v
    fig.canvas.mpl_connect("motion_notify_event", v.on_move)
    fig.canvas.mpl_connect("button_press_event", v.on_click)
    fig.canvas.mpl_connect("key_press_event", v.on_key)
    plt.show()
    return v


__all__ = ["Viewer", "gap_anchor", "show"]
