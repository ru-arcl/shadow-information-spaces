"""Field-of-view events in the nondeterministic setting (T-RO 2012, Sec. V-D).

There is no null FOV event here.  By Observation 2 an enter or exit event can
be rewritten as component events:

* ``Enter(s, k)``  ==  ``Appear(t, k, k)`` + ``Merge(s, t, s')``;
* ``Exit(s, k)``   ==  ``Split(s, s', t)`` + ``Disappear(t, k, k)``.

Doing this per event creates two component events per FOV event.  Sec. V-D
instead keeps, per shadow, the running net flow ``d_j`` (``+k`` per enter,
``-k`` per exit), its minimum ``d_min`` (the "surplus" the shadow must have
held) and the total ``d_tot``.  The whole run of FOV events on a shadow between
two component events is then equivalent to at most one batch exit of
``|d_min|`` targets followed by one batch enter of ``d_tot - d_min`` targets.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

from .events import Appear, Disappear, Enter, Event, Exit, Merge, ShadowSequence, Split


@dataclass
class FovBatch:
    """Running ``(d_tot, d_min)`` counter for the FOV events of one shadow (Sec. V-D)."""

    d_tot: int = 0
    d_min: int = 0

    def add(self, e: Event) -> None:
        if isinstance(e, Enter):
            self.d_tot += e.k
        elif isinstance(e, Exit):
            self.d_tot -= e.k
            self.d_min = min(self.d_min, self.d_tot)
        else:
            raise TypeError(f"not a FOV event: {e}")

    def batch_events(self, s: int) -> List[Event]:
        """At most one batch exit then one batch enter, equivalent to the run (four cases)."""
        out: List[Event] = []
        if self.d_min < 0:
            out.append(Exit(s, -self.d_min))
        if self.d_tot - self.d_min > 0:
            out.append(Enter(s, self.d_tot - self.d_min))
        return out


def batch_fov(events: List[Event]) -> List[Event]:
    """Replace every run of FOV events per shadow by its batch form.

    Pending batches are emitted right before the shadow takes part in a
    component event, and at the end of the sequence for shadows still alive.
    """
    pending: dict = {}
    out: List[Event] = []

    def flush(s: int) -> None:
        b = pending.pop(s, None)
        if b is not None:
            out.extend(b.batch_events(s))

    for e in events:
        if isinstance(e, (Enter, Exit)):
            pending.setdefault(e.s, FovBatch()).add(e)
            continue
        if isinstance(e, (Split, Disappear)):
            flush(e.s)
        elif isinstance(e, Merge):
            flush(e.a)
            flush(e.b)
        out.append(e)
    for s in list(pending):
        flush(s)
    return out


def fov_to_component(seq: ShadowSequence, first_label: Optional[int] = None) -> Tuple[ShadowSequence, dict]:
    """Rewrite all FOV events of ``seq`` as component events (naive conversion).

    Fresh labels start at ``first_label`` (default: one above the largest label
    in use).  Returns the new sequence and a map from each original label to the
    label that carries it at the end (``s`` is renamed after each FOV event).
    """
    nxt = (max(seq.all_labels(), default=0) + 1) if first_label is None else first_label
    current = {s: s for s in seq.all_labels()}
    out: List[Event] = []

    def fresh() -> int:
        nonlocal nxt
        nxt += 1
        return nxt - 1

    for e in seq.events:
        if isinstance(e, Enter):
            t, s2 = fresh(), fresh()
            out += [Appear(t, e.k, e.k), Merge(current[e.s], t, s2)]
            current[e.s] = s2
        elif isinstance(e, Exit):
            s2, t = fresh(), fresh()
            out += [Split(current[e.s], s2, t), Disappear(t, e.k, e.k)]
            current[e.s] = s2
        elif isinstance(e, Appear):
            out.append(e)
        elif isinstance(e, Disappear):
            out.append(Disappear(current[e.s], e.lo, e.hi))
        elif isinstance(e, Split):
            out.append(Split(current[e.s], e.a, e.b))
        else:
            out.append(Merge(current[e.a], current[e.b], e.s))
    return ShadowSequence(dict(seq.initial), out), current


__all__ = ["FovBatch", "batch_fov", "fov_to_component"]
