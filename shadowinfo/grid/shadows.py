"""Shadow components and frame-to-frame component events (DESIGN §5).

The shadow region (free cells not seen by the robot) is split into
4-connected components, numbered in row-major scan order of their first cell.
Between two consecutive frames, old shadows (labelled) and new components are
joined whenever they share a cell; every connected group of this overlap
bipartite graph with ``m`` old and ``n`` new components yields the component
events of T-RO 2012, Sec. III-A:

* ``(0, 1)`` appear, ``(1, 0)`` disappear, ``(1, 1)`` the label lives on,
* ``(1, 2)`` split, ``(2, 1)`` merge,
* general ``(m, n)``: merge the old shadows one by one (ascending label), then
  split off the new components one by one (scan order).

Labels are integers allocated by a global counter and never reused.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Tuple

from ..events import Appear, Disappear, Event, Merge, Split


def label_components(mask: Sequence[int], rows: int, cols: int) -> Tuple[List[int], List[List[int]]]:
    """4-connected components of the cells where ``mask`` is nonzero.

    Returns ``(comp, cells)``: ``comp[i]`` is ``-1`` or the 0-based component
    index of cell ``i`` (components ordered by their first cell in row-major
    order) and ``cells[k]`` lists the cells of component ``k`` in ascending
    order.
    """
    n = rows * cols
    comp = [-1] * n
    cells: List[List[int]] = []
    for start in range(n):
        if not mask[start] or comp[start] >= 0:
            continue
        k = len(cells)
        comp[start] = k
        members = [start]
        stack = [start]
        while stack:
            i = stack.pop()
            r, c = divmod(i, cols)
            for j, ok in ((i - cols, r > 0), (i + cols, r < rows - 1), (i - 1, c > 0), (i + 1, c < cols - 1)):
                if ok and mask[j] and comp[j] < 0:
                    comp[j] = k
                    members.append(j)
                    stack.append(j)
        members.sort()
        cells.append(members)
    return comp, cells


@dataclass
class Transition:
    """Result of :meth:`ShadowTracker.update` for one frame change.

    ``events`` are the component events in order; ``labels`` is the new flat
    label grid (0 = not shadow); ``appeared`` / ``disappeared`` list the labels
    of ``(0, 1)`` / ``(1, 0)`` groups; ``covers`` maps each split-created label
    to the new components (indices into ``cells``) it covers, so that the
    ground-truth count of such a shadow is the sum over them, and
    ``merged_from`` maps each merge-created label to the old labels it absorbs.
    """

    events: List[Event]
    labels: List[int]
    cells: List[List[int]]
    comp_labels: List[int]
    appeared: List[int] = field(default_factory=list)
    disappeared: List[int] = field(default_factory=list)
    covers: Dict[int, List[int]] = field(default_factory=dict)
    merged_from: Dict[int, List[int]] = field(default_factory=dict)


class ShadowTracker:
    """Maintains shadow labels across frames and emits component events."""

    def __init__(self, rows: int, cols: int) -> None:
        self.rows = rows
        self.cols = cols
        self.labels: List[int] = [0] * (rows * cols)
        self.next_label = 1

    def _new_label(self) -> int:
        s = self.next_label
        self.next_label += 1
        return s

    def reset(self, mask: Sequence[int]) -> List[int]:
        """Label the initial shadows ``1..k`` in scan order; return the labels."""
        comp, cells = label_components(mask, self.rows, self.cols)
        self.next_label = 1
        comp_labels = [self._new_label() for _ in cells]
        self.labels = [comp_labels[k] if k >= 0 else 0 for k in comp]
        return comp_labels

    def alive(self) -> List[int]:
        return sorted(set(self.labels) - {0})

    def update(self, mask: Sequence[int]) -> Transition:
        """Advance to a new shadow mask and return the component events.

        ``Appear`` / ``Disappear`` events carry ``lo = hi = 0`` here; the
        simulator fills in the revealed target counts.
        """
        old = self.labels
        comp, cells = label_components(mask, self.rows, self.cols)
        n_new = len(cells)
        old_labels = sorted(set(old) - {0})
        node_of = {lab: n_new + j for j, lab in enumerate(old_labels)}
        parent = list(range(n_new + len(old_labels)))

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for i, lab in enumerate(old):
            if lab and comp[i] >= 0:
                rx, ry = find(comp[i]), find(node_of[lab])
                if rx != ry:
                    parent[max(rx, ry)] = min(rx, ry)

        news_of: Dict[int, List[int]] = {}
        olds_of: Dict[int, List[int]] = {}
        for k in range(n_new):
            news_of.setdefault(find(k), []).append(k)
        for lab in old_labels:
            olds_of.setdefault(find(node_of[lab]), []).append(lab)

        events: List[Event] = []
        comp_labels = [0] * n_new
        tr = Transition(events, [], cells, comp_labels)
        for root in sorted(olds_of):
            if root not in news_of:
                lab = olds_of[root][0]
                events.append(Disappear(lab))
                tr.disappeared.append(lab)
        for root in sorted(news_of):
            news = news_of[root]
            olds = olds_of.get(root, [])
            if not olds:
                s = self._new_label()
                events.append(Appear(s))
                tr.appeared.append(s)
                comp_labels[news[0]] = s
                continue
            cur = olds[0]
            for j in range(1, len(olds)):
                merged = self._new_label()
                events.append(Merge(cur, olds[j], merged))
                tr.merged_from[merged] = olds[:j + 1]
                cur = merged
            for j in range(len(news) - 1):
                a = self._new_label()
                b = self._new_label()
                events.append(Split(cur, a, b))
                comp_labels[news[j]] = a
                tr.covers[a] = [news[j]]
                tr.covers[b] = news[j + 1:]
                cur = b
            comp_labels[news[-1]] = cur
        self.labels = [comp_labels[k] if k >= 0 else 0 for k in comp]
        tr.labels = self.labels
        return tr
