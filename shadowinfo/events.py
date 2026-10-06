"""Critical events over shadows (T-RO 2012, Secs. II-III).

Two kinds of critical events change what can be known about hidden targets:

* component events -- a shadow appears, disappears, splits or merges;
* field-of-view (FOV) events -- targets enter or exit a shadow through the
  boundary between the shadow and the sensors' field of view.

Shadows are identified by integer labels that are never reused.  A
:class:`ShadowSequence` bundles the initial condition (bounds on the number of
targets in each shadow at ``t = t0``) with a chronologically ordered list of
events, and validates that the labels are used consistently.

Upper bounds may be ``math.inf`` (serialised as JSON ``null``).
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple, Union

Bound = Tuple[int, float]  # (lo, hi); hi may be math.inf


@dataclass(frozen=True)
class Appear:
    """New shadow ``s`` appears holding between ``lo`` and ``hi`` targets."""

    s: int
    lo: int = 0
    hi: float = 0


@dataclass(frozen=True)
class Disappear:
    """Shadow ``s`` disappears, revealing between ``lo`` and ``hi`` targets."""

    s: int
    lo: int = 0
    hi: float = 0


@dataclass(frozen=True)
class Split:
    """Shadow ``s`` splits into shadows ``a`` and ``b``."""

    s: int
    a: int
    b: int


@dataclass(frozen=True)
class Merge:
    """Shadows ``a`` and ``b`` merge into shadow ``s``."""

    a: int
    b: int
    s: int


@dataclass(frozen=True)
class Enter:
    """``k`` targets enter shadow ``s`` from the field of view."""

    s: int
    k: int = 1


@dataclass(frozen=True)
class Exit:
    """``k`` targets exit shadow ``s`` into the field of view."""

    s: int
    k: int = 1


Event = Union[Appear, Disappear, Split, Merge, Enter, Exit]
ComponentEvent = (Appear, Disappear, Split, Merge)
FovEvent = (Enter, Exit)

_TYPES = {
    "appear": Appear,
    "disappear": Disappear,
    "split": Split,
    "merge": Merge,
    "enter": Enter,
    "exit": Exit,
}
_NAMES = {cls: name for name, cls in _TYPES.items()}


class InvalidSequenceError(ValueError):
    """Raised when an event refers to labels inconsistently."""


def _hi_to_json(hi: float) -> Optional[float]:
    return None if hi == math.inf else int(hi)


def _hi_from_json(hi: Optional[float]) -> float:
    return math.inf if hi is None else int(hi)


def event_to_dict(e: Event) -> dict:
    d = {"type": _NAMES[type(e)]}
    for k, v in e.__dict__.items():
        d[k] = _hi_to_json(v) if k == "hi" else v
    return d


def event_from_dict(d: dict) -> Event:
    cls = _TYPES[d["type"]]
    kw = {k: v for k, v in d.items() if k != "type"}
    if "hi" in kw:
        kw["hi"] = _hi_from_json(kw["hi"])
    return cls(**kw)


def event_labels(e: Event) -> Tuple[Tuple[int, ...], Tuple[int, ...]]:
    """Return ``(consumed, created)`` labels of an event.

    FOV events neither consume nor create labels; they require ``s`` alive.
    """
    if isinstance(e, Appear):
        return (), (e.s,)
    if isinstance(e, Disappear):
        return (e.s,), ()
    if isinstance(e, Split):
        return (e.s,), (e.a, e.b)
    if isinstance(e, Merge):
        return (e.a, e.b), (e.s,)
    return (), ()


@dataclass
class ShadowSequence:
    """Initial condition plus an ordered list of critical events."""

    initial: Dict[int, Bound]
    events: List[Event] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.initial = {int(k): (int(v[0]), float(v[1]) if v[1] == math.inf else int(v[1]))
                        for k, v in self.initial.items()}
        self.validate()

    # -- validation -----------------------------------------------------
    def validate(self) -> None:
        alive = set(self.initial)
        seen = set(self.initial)
        for lo, hi in self.initial.values():
            if lo < 0 or hi < lo:
                raise InvalidSequenceError(f"bad initial bound {(lo, hi)}")
        for i, e in enumerate(self.events):
            consumed, created = event_labels(e)
            if isinstance(e, FovEvent):
                if e.s not in alive:
                    raise InvalidSequenceError(f"event {i} {e}: shadow {e.s} not alive")
                if e.k < 1:
                    raise InvalidSequenceError(f"event {i} {e}: k must be >= 1")
                continue
            if isinstance(e, (Appear, Disappear)) and (e.lo < 0 or e.hi < e.lo):
                raise InvalidSequenceError(f"event {i} {e}: bad bound")
            if len(set(consumed)) != len(consumed) or len(set(created)) != len(created):
                raise InvalidSequenceError(f"event {i} {e}: repeated label")
            for s in consumed:
                if s not in alive:
                    raise InvalidSequenceError(f"event {i} {e}: shadow {s} not alive")
            for s in created:
                if s in seen:
                    raise InvalidSequenceError(f"event {i} {e}: label {s} reused")
            alive.difference_update(consumed)
            alive.update(created)
            seen.update(created)

    # -- queries --------------------------------------------------------
    def alive_at_end(self) -> List[int]:
        alive = list(self.initial)
        for e in self.events:
            consumed, created = event_labels(e)
            alive = [s for s in alive if s not in consumed] + list(created)
        return sorted(alive)

    def all_labels(self) -> List[int]:
        labels = set(self.initial)
        for e in self.events:
            labels.update(event_labels(e)[1])
        return sorted(labels)

    def append(self, e: Event) -> None:
        self.events.append(e)
        try:
            self.validate()
        except InvalidSequenceError:
            self.events.pop()
            raise

    def extend(self, events: Iterable[Event]) -> None:
        for e in events:
            self.append(e)

    # -- (de)serialisation ------------------------------------------------
    def to_dict(self) -> dict:
        return {
            "initial": {str(k): [lo, _hi_to_json(hi)] for k, (lo, hi) in self.initial.items()},
            "events": [event_to_dict(e) for e in self.events],
        }

    @classmethod
    def from_dict(cls, d: dict) -> "ShadowSequence":
        initial = {int(k): (int(v[0]), _hi_from_json(v[1])) for k, v in d["initial"].items()}
        return cls(initial, [event_from_dict(e) for e in d["events"]])

    def to_json(self, **kw) -> str:
        return json.dumps(self.to_dict(), **kw)

    @classmethod
    def from_json(cls, s: str) -> "ShadowSequence":
        return cls.from_dict(json.loads(s))


__all__ = [
    "Bound", "Appear", "Disappear", "Split", "Merge", "Enter", "Exit", "Event", "ComponentEvent", "FovEvent",
    "InvalidSequenceError", "ShadowSequence", "event_to_dict", "event_from_dict", "event_labels",
]
