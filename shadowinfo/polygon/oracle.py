"""The simulated targets: ``java.util.Random`` and the single-type-agent oracle.

Port of ``pe.is.Oracle``, ``pe.is.Event`` and ``pe.is.singleTypeAgent.SingleTypeAgentOracle`` /
``SingleTypeAgentEvent`` of the original Java implementation (``docs/notes/original_java.md`` §2.9).
The oracle walks a gap history (:func:`~shadowinfo.polygon.gaps.get_gaps`) and distributes ``N``
indistinguishable targets at random as the shadows evolve.  This produces what the robot observes
(T-RO 2012 Secs. IV-B and V): the number of targets hidden at the start, the number that enter each
newly appeared shadow (taken from the visible ones) and the number revealed by each disappearing
shadow.  Splits share a shadow's targets between its two parts; merges add them up.

All randomness is ``Math.random()``, i.e. one ``java.util.Random``.  :class:`JavaRandom` is an
exact clone of it (48-bit LCG), so a seeded oracle reproduces the Java runs of the golden fixtures
(``tests/fixtures/java/*.json``, ``oracle_runs``; the harness seeds ``Math.random()``) number for
number.

Two modes (``compat``):

* ``"java"`` -- ``SingleTypeAgentOracle.distributeAgents`` statement by statement: it compares the
  sizes of consecutive gap sets and assumes one event per transition (B9), stores events in a
  ``TreeMap`` keyed by relative time (equal times overwrite), wraps ``int`` arithmetic, leaves the
  type ``None`` where no event is recognised and raises :class:`OracleError` where Java throws.
  ``merge_bug`` (default ``True`` here) reproduces B7: the second merging parent sets the merged
  count to *its own count twice*, so the first parent's targets are lost and later observations are
  inconsistent.  ``merge_bug=False`` is the harness's ``MergeFixedOracle``
  (``tools/java_reference/src/MergeFixedOracle.java``, the one-line fix).
* ``"fixed"`` (default) -- the same random draws in the same order, driven by the component events of
  the history (:func:`~shadowinfo.polygon.events.transition_events`).  It therefore handles histories
  with several events per transition and splits/merges of more than two shadows (``compat="fixed"``
  gap histories have those), conserves targets (``merge_bug`` defaults to ``False``), and never loops
  forever (B8).  On every Java-shaped history (one event per transition, distinct times) it gives
  exactly the events of ``compat="java"`` with the same ``merge_bug``.

:meth:`SingleTypeAgentOracle.sequence` turns the observations into a
:class:`~shadowinfo.events.ShadowSequence` for the filters of the package; the hidden total at the
start is :attr:`SingleTypeAgentOracle.hidden_total` (the Java filter's pooled left vertex 0, see
:mod:`.stagent`).
"""

from __future__ import annotations

import enum
import math
import os
import struct
import time as _time
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple, Union

from ..events import Appear, Disappear, Event, Merge, ShadowSequence, Split
from .events import transition_events
from .gaps import Gap, _check_compat

__all__ = [
    "JavaRandom", "strict_log", "OracleError", "EventType", "SingleTypeAgentEvent", "SingleTypeAgentOracle", "java_int",
]

_NPE = "java.lang.NullPointerException"
_AIOOBE = "java.lang.ArrayIndexOutOfBoundsException"


def java_int(x: int) -> int:
    """Wrap an integer to Java ``int`` (32-bit two's complement), as Java ``int`` arithmetic does."""
    x &= 0xFFFFFFFF
    return x - 0x100000000 if x & 0x80000000 else x


def _java_long(x: int) -> int:
    x &= 0xFFFFFFFFFFFFFFFF
    return x - 0x10000000000000000 if x & 0x8000000000000000 else x


def _d2i(x: float) -> int:
    """Java ``(int) x`` for a double: truncation toward zero, NaN -> 0, saturating at the ``int`` range."""
    if math.isnan(x):
        return 0
    if x >= 2147483647.0:
        return 2147483647
    if x <= -2147483648.0:
        return -2147483648
    return int(x)


# --------------------------------------------------------------------------- java.util.Random

_LN2_HI = 6.93147180369123816490e-01
_LN2_LO = 1.90821492927058770002e-10
_TWO54 = 1.80143985094819840000e+16
_LG = (6.666666666666735130e-01, 3.999999999940941908e-01, 2.857142874366239149e-01, 2.222219843214978396e-01,
       1.818357216161805012e-01, 1.531383769920937332e-01, 1.479819860511658591e-01)


def _words(x: float) -> Tuple[int, int]:
    b = struct.unpack("<q", struct.pack("<d", x))[0]
    return java_int(b >> 32), b & 0xFFFFFFFF


def _from_words(hi: int, lo: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", ((hi & 0xFFFFFFFF) << 32) | lo))[0]


def strict_log(x: float) -> float:
    """``StrictMath.log`` (fdlibm ``__ieee754_log``), bit for bit.  ``math.log`` may differ in the last bit."""
    hx, lx = _words(x)
    k = 0
    if hx < 0x00100000:
        if ((hx & 0x7FFFFFFF) | lx) == 0:
            return -math.inf
        if hx < 0:
            return math.nan
        k -= 54
        x *= _TWO54
        hx = _words(x)[0]
    if hx >= 0x7FF00000:
        return x + x
    k += (hx >> 20) - 1023
    hx &= 0x000FFFFF
    i = (hx + 0x95F64) & 0x100000
    x = _from_words(hx | (i ^ 0x3FF00000), _words(x)[1])
    k += i >> 20
    f = x - 1.0
    if (0x000FFFFF & (2 + hx)) < 3:
        if f == 0.0:
            if k == 0:
                return 0.0
            dk = float(k)
            return dk * _LN2_HI + dk * _LN2_LO
        r = f * f * (0.5 - 0.33333333333333333 * f)
        if k == 0:
            return f - r
        dk = float(k)
        return dk * _LN2_HI - ((r - dk * _LN2_LO) - f)
    lg1, lg2, lg3, lg4, lg5, lg6, lg7 = _LG
    s = f / (2.0 + f)
    dk = float(k)
    z = s * s
    i = hx - 0x6147A
    w = z * z
    j = 0x6B851 - hx
    t1 = w * (lg2 + w * (lg4 + w * lg6))
    t2 = z * (lg1 + w * (lg3 + w * (lg5 + w * lg7)))
    i |= j
    r = t2 + t1
    if i > 0:
        hfsq = 0.5 * f * f
        if k == 0:
            return f - (hfsq - s * (hfsq + r))
        return dk * _LN2_HI - ((hfsq - (s * (hfsq + r) + dk * _LN2_LO)) - f)
    if k == 0:
        return f - s * (f - r)
    return dk * _LN2_HI - ((s * (f - r) - dk * _LN2_LO) - f)



class JavaRandom:
    """Exact clone of ``java.util.Random`` (the generator behind ``Math.random()``).

    48-bit linear congruential generator ``seed = seed * 0x5DEECE66D + 0xB mod 2^48``; ``next(bits)``
    returns the top ``bits`` bits as a Java ``int``.  ``Random(seed)`` scrambles the seed with
    ``0x5DEECE66D``.  :meth:`next_double` is ``((long) next(26) << 27) + next(27)) * 2^-53``, which
    is what ``Math.random()`` returns; the other methods follow OpenJDK 17 ``java.util.Random``.

    ``JavaRandom()`` without a seed is seeded from the OS (like Java's unseeded constructor it is not
    reproducible).  :meth:`next_gaussian` uses :func:`strict_log` (fdlibm, as ``StrictMath.log``) and
    ``math.sqrt`` (correctly rounded, as ``StrictMath.sqrt``).
    """

    MULTIPLIER = 0x5DEECE66D
    ADDEND = 0xB
    MASK = (1 << 48) - 1

    __slots__ = ("_seed", "_next_next_gaussian", "_have_next_next_gaussian")

    def __init__(self, seed: Optional[int] = None):
        if seed is None:
            seed = int.from_bytes(os.urandom(8), "little") ^ _time.perf_counter_ns()
        self._seed = 0
        self._next_next_gaussian = 0.0
        self._have_next_next_gaussian = False
        self.set_seed(seed)

    def set_seed(self, seed: int) -> None:
        """``setSeed(seed)``: ``(seed ^ 0x5DEECE66D) & (2^48 - 1)``; clears the cached Gaussian."""
        self._seed = (int(seed) ^ self.MULTIPLIER) & self.MASK
        self._have_next_next_gaussian = False

    @property
    def state(self) -> int:
        """The current 48-bit internal state."""
        return self._seed

    def next(self, bits: int) -> int:
        """``next(bits)``: advance the LCG and return its top ``bits`` bits as a Java ``int``."""
        self._seed = (self._seed * self.MULTIPLIER + self.ADDEND) & self.MASK
        return java_int(self._seed >> (48 - bits))

    def next_int(self, bound: Optional[int] = None) -> int:
        """``nextInt()`` (any ``int``) or ``nextInt(bound)`` (uniform in ``[0, bound)``, rejection sampling)."""
        if bound is None:
            return self.next(32)
        if bound <= 0:
            raise ValueError("bound must be positive")
        r = self.next(31)
        m = bound - 1
        if bound & m == 0:  # power of two
            return java_int((bound * r) >> 31)
        u = r
        while True:
            r = u % bound
            if java_int(u - r + m) >= 0:
                return r
            u = self.next(31)

    def next_long(self) -> int:
        """``nextLong()``: ``((long) next(32) << 32) + next(32)``."""
        hi = self.next(32)
        return _java_long((hi << 32) + self.next(32))

    def next_boolean(self) -> bool:
        """``nextBoolean()``: ``next(1) != 0``."""
        return self.next(1) != 0

    def next_float(self) -> float:
        """``nextFloat()``: ``next(24) / 2^24``."""
        return self.next(24) / float(1 << 24)

    def next_double(self) -> float:
        """``nextDouble()`` (= ``Math.random()``): ``(((long) next(26) << 27) + next(27)) * 2^-53``."""
        return ((self.next(26) << 27) + self.next(27)) * (1.0 / (1 << 53))

    def next_gaussian(self) -> float:
        """``nextGaussian()``: Marsaglia's polar method with ``StrictMath``, caching the second value."""
        if self._have_next_next_gaussian:
            self._have_next_next_gaussian = False
            return self._next_next_gaussian
        while True:
            v1 = 2 * self.next_double() - 1
            v2 = 2 * self.next_double() - 1
            s = v1 * v1 + v2 * v2
            if s < 1 and s != 0:
                break
        multiplier = math.sqrt(-2 * strict_log(s) / s)
        self._next_next_gaussian = v2 * multiplier
        self._have_next_next_gaussian = True
        return v1 * multiplier

    random = next_double  # Math.random()


# --------------------------------------------------------------------------- events


class OracleError(RuntimeError):
    """The Java oracle or filter throws (``compat="java"``), or an input cannot be processed.

    ``java_exception`` is the Java exception class (``java.lang.NullPointerException`` ...) when the
    error reproduces a Java crash, else ``""``.
    """

    def __init__(self, message: str, *, java_exception: str = ""):
        super().__init__(message)
        self.java_exception = java_exception


class EventType(enum.Enum):
    """``Oracle.EVENT_TYPE``."""

    INIT = "INIT"
    SPLIT = "SPLIT"
    MERGE = "MERGE"
    APPEAR = "APPEAR"
    DISAPPEAR = "DISAPPEAR"
    AGENT_APPEAR = "AGENT_APPEAR"
    AGENT_DISAPPEAR = "AGENT_DISAPPEAR"


@dataclass
class SingleTypeAgentEvent:
    """``SingleTypeAgentEvent``: one observation of the oracle.

    * ``type`` -- :class:`EventType` (``None``: Java recognised no event in that transition, B9);
    * ``time`` -- relative time (arc-length fraction) of the gap set the event leads to;
    * ``from_gaps`` / ``to_gaps`` -- shadow IDs before / after (``None`` = unset, Java ``null``).  INIT:
      ``to`` = the initial shadows.  SPLIT ``s -> a b``: ``from=[s]``, ``to=[a, b]``.  MERGE
      ``a b -> s``: ``from=[a, b]``, ``to=[s]``.  APPEAR: ``to=[s]``.  DISAPPEAR: ``from=[s]``;
    * ``moved`` -- targets that entered the appeared shadow, or were revealed by the disappeared one;
    * ``visible`` -- targets visible to the robot after the event.

    ``str()`` is Java's ``toString()`` (``INIT EVENT 1 2 3 4``, ``APPEAR EVENT 5 <- 2``,
    ``DISAPPEAR EVENT 14 -> 36``, ``SPLIT EVENT 1 -> 6 7``, ``MERGE EVENT 7 2 -> 15``).
    """

    type: Optional[EventType] = None
    time: float = 0.0
    from_gaps: Optional[List[int]] = None
    to_gaps: Optional[List[int]] = None
    visible: int = 0
    moved: int = 0

    def __str__(self) -> str:
        t, f, to = self.type, self.from_gaps, self.to_gaps
        if t is None:
            raise OracleError("SingleTypeAgentEvent.toString: event type is null", java_exception=_NPE)
        try:
            if t is EventType.AGENT_APPEAR:
                return f"AGENT_APPEAR: {self.moved} appeared from gap {f[0]}"  # type: ignore[index]
            if t is EventType.AGENT_DISAPPEAR:
                return f"AGENT_DISAPPEAR: {self.moved} disappeared from gap {to[0]}"  # type: ignore[index]
            if t is EventType.INIT:
                return "INIT EVENT" + "".join(f" {g}" for g in to)  # type: ignore[union-attr]
            if t is EventType.APPEAR:
                return f"APPEAR EVENT {to[0]} <- {self.moved}"  # type: ignore[index]
            if t is EventType.DISAPPEAR:
                return f"DISAPPEAR EVENT {f[0]} -> {self.moved}"  # type: ignore[index]
            if t is EventType.SPLIT:
                return f"SPLIT EVENT {f[0]} -> {to[0]} {to[1]}"  # type: ignore[index]
            return f"MERGE EVENT {f[0]} {f[1]} -> {to[0]}"  # type: ignore[index]
        except TypeError as ex:  # a null array
            raise OracleError(f"SingleTypeAgentEvent.toString: {ex}", java_exception=_NPE) from None

    def to_row(self) -> list:
        """``[time, type, fromGap, toGap, movedAgents, visibleAgents]`` (golden-fixture ``events`` format)."""
        return [self.time, None if self.type is None else self.type.value,
                None if self.from_gaps is None else list(self.from_gaps),
                None if self.to_gaps is None else list(self.to_gaps), self.moved, self.visible]

    def shadow_event(self) -> Optional[Event]:
        """The observation as a :mod:`shadowinfo.events` event (``None`` for INIT and untyped events).

        Counts are exact: ``Appear(s, moved, moved)``, ``Disappear(s, moved, moved)``.
        """
        t = self.type
        if t is EventType.APPEAR:
            return Appear(self.to_gaps[0], self.moved, self.moved)  # type: ignore[index]
        if t is EventType.DISAPPEAR:
            return Disappear(self.from_gaps[0], self.moved, self.moved)  # type: ignore[index]
        if t is EventType.SPLIT:
            return Split(self.from_gaps[0], self.to_gaps[0], self.to_gaps[1])  # type: ignore[index]
        if t is EventType.MERGE:
            return Merge(self.from_gaps[0], self.from_gaps[1], self.to_gaps[0])  # type: ignore[index]
        if t in (EventType.AGENT_APPEAR, EventType.AGENT_DISAPPEAR):
            raise OracleError(f"{t.value} events are not produced by this oracle (B13)")
        return None


def _time_key(t: float) -> Tuple[int, float, int]:
    """Sort/equality key of ``Double.compareTo`` (``TreeMap<Double, ...>``): ``-0.0 < 0.0``, NaN last, NaN == NaN."""
    if math.isnan(t):
        return (1, 0.0, 0)
    return (0, t, -1 if t == 0 and math.copysign(1.0, t) < 0 else 0)


# --------------------------------------------------------------------------- the oracle


def _history_times(history: Sequence[Sequence[Gap]]) -> List[float]:
    times = getattr(history, "times", None)
    if times is not None and len(times) == len(history):
        return list(times)
    return [gs[0].relative_time if gs else 0.0 for gs in history]


class SingleTypeAgentOracle:
    """``SingleTypeAgentOracle(N)``: distributes ``N`` targets at random as the gaps evolve.

    ``rng`` is a :class:`JavaRandom`, a seed for one, or ``None`` (an unseeded one, like
    ``Math.random()``).  ``merge_bug`` reproduces quirk B7 (default: ``True`` for ``compat="java"``,
    ``False`` for ``compat="fixed"``).  See the module docstring for the two modes.

    :meth:`initialize` runs the distribution on a gap history (``Oracle.initialize``).  Afterwards:

    * :attr:`events` -- the observations in time order (``getEvent(i)``; Java ``TreeMap`` order),
      :meth:`event_by_time` (``getEventByTime``), :attr:`number_of_events`;
    * :attr:`visible0` / :attr:`hidden_total` -- targets visible / hidden at the start (``N - visible0``);
    * :attr:`counts` -- the ground truth: for every gap set, ``{id: targets}`` (``None`` where Java left a
      gap without a state); :attr:`final_counts` and :meth:`truth` for the last set; :attr:`label_counts`
      -- the targets of every shadow label when it was created (including the intermediate labels of
      ``compat="fixed"`` multi-event transitions), against which every equation holds unless the merge
      bug is on;
    * :meth:`sequence` -- the observations as a :class:`~shadowinfo.events.ShadowSequence`.
    """

    def __init__(self, n_targets: int, rng: Union[JavaRandom, int, None] = None, merge_bug: Optional[bool] = None,
                 compat: str = "fixed"):
        _check_compat(compat)
        if n_targets < 0:
            raise ValueError("the number of targets must be non-negative")
        self.n_targets = int(n_targets)
        self.rng = rng if isinstance(rng, JavaRandom) else JavaRandom(rng)
        self.compat = compat
        self.merge_bug = (compat == "java") if merge_bug is None else bool(merge_bug)
        self.events: List[SingleTypeAgentEvent] = []
        self.counts: List[Dict[int, Optional[int]]] = []
        self.label_counts: Dict[int, Optional[int]] = {}
        self.set_ids: List[List[int]] = []
        self._by_time: Dict[Tuple[int, float, int], SingleTypeAgentEvent] = {}
        self.visible0 = 0

    # -- Oracle interface ------------------------------------------------------------
    def get_total_agent_number(self) -> int:
        """``getTotalAgentNumber()``: ``N``."""
        return self.n_targets

    @property
    def number_of_events(self) -> int:
        """``getNumberOfEvents()``."""
        return len(self.events)

    def get_event(self, i: int) -> SingleTypeAgentEvent:
        """``getEvent(i)``."""
        return self.events[i]

    def event_by_time(self, t: float) -> Optional[SingleTypeAgentEvent]:
        """``getEventByTime(t)``: the event stored under relative time ``t`` (``None`` if there is none).

        With ``compat="fixed"`` several events can share a time; this returns the last of them.
        """
        return self._by_time.get(_time_key(t))

    @property
    def hidden_total(self) -> int:
        """Targets hidden in the initial shadows: ``N - visible0`` (Java int arithmetic in java mode)."""
        h = self.n_targets - self.visible0
        return java_int(h) if self.compat == "java" else h

    @property
    def final_counts(self) -> Dict[int, Optional[int]]:
        """``{id: targets}`` of the last gap set (the shadows alive at the end of the path)."""
        return dict(self.counts[-1]) if self.counts else {}

    def truth(self, ids: Optional[Sequence[int]] = None) -> List[Optional[int]]:
        """Targets in the shadows ``ids`` (default: the last gap set, in its order) at the end.

        This is the golden fixtures' ``truth_final`` (aligned with ``get_gaps.final_ids``).
        """
        fc = self.counts[-1] if self.counts else {}
        return [fc.get(g) for g in (self.set_ids[-1] if ids is None and self.set_ids else ids or [])]

    def initialize(self, history: Sequence[Sequence[Gap]]) -> "SingleTypeAgentOracle":
        """``Oracle.initialize(gapss)``: distribute the targets along a gap history.  Returns ``self``."""
        self.events = []
        self.counts = []
        self.label_counts = {}
        self._by_time = {}
        self.set_ids = [[g.id for g in gs] for gs in history]
        if self.compat == "java":
            self._distribute_java(history)
        else:
            self._distribute_fixed(history)
        return self

    # -- distributeOneBatch --------------------------------------------------------------
    def distribute_one_batch(self, t: int, g: int) -> List[int]:
        """``distributeOneBatch(t, g)``: split ``t`` targets into ``g`` random parts.

        Draws ``g - 1`` *distinct* cut points ``(int)(rand * t)`` in ``[0, t)`` (a ``TreeSet``), adds ``t``
        and returns the differences, so the last part is ``>= 1`` whenever ``t > 0`` (B8: an appearing
        shadow never takes all visible targets; with one visible target it always takes 0).  ``g == 1``
        returns ``[t]`` and ``t == 0`` returns zeros without drawing.

        If ``t < g - 1`` Java loops forever (B8): ``compat="java"`` raises :class:`OracleError`,
        ``compat="fixed"`` draws the cut points in ``[0, t]`` with repetition instead.
        """
        if g == 1:
            return [t]
        if t == 0:
            return [0] * max(g, 0)
        rand = self.rng.next_double
        if self.compat == "fixed":
            if g <= 0:
                return []
            if t < g - 1:
                cuts = sorted(_d2i(rand() * (t + 1)) for _ in range(g - 1)) + [t]
                return [cuts[0]] + [cuts[i] - cuts[i - 1] for i in range(1, g)]
        elif 0 < t < g - 1:
            raise OracleError(f"distributeOneBatch({t}, {g}) never terminates in Java (needs {g - 1} distinct "
                              f"cut points in [0, {t}))")
        sp = set()
        while len(sp) < g - 1:
            sp.add(_d2i(rand() * t))
        sp.add(t)
        ret_i = sorted(sp)
        if g <= 0:
            raise OracleError(f"distributeOneBatch({t}, {g}): index 0 out of bounds for length {g}",
                              java_exception=_AIOOBE)
        ret = [0] * g
        ret[0] = ret_i[0]
        for i in range(1, g):
            ret[i] = java_int(ret_i[i] - ret_i[i - 1])
        return ret

    # -- compat="java" -----------------------------------------------------------------
    def _distribute_java(self, gapss: Sequence[Sequence[Gap]]) -> None:
        """``SingleTypeAgentOracle.distributeAgents`` (``MergeFixedOracle`` with ``merge_bug=False``)."""
        n = self.n_targets
        prev_vis = 0

        def state(cnt: Dict[int, Optional[int]], gid: int, k: int) -> int:
            v = cnt.get(gid)
            if v is None:
                raise OracleError(f"distributeAgents, set {k}: gap {gid} has no state "
                                  "(SingleTypeAgentState.getNumberOfAgents on null)", java_exception=_NPE)
            return v

        def gap(cnt: Dict[int, Optional[int]], gid: int, k: int) -> None:
            if gid not in cnt:
                raise OracleError(f"distributeAgents, set {k}: gap {gid} is not in the set "
                                  "(Map.get returns null)", java_exception=_NPE)

        for i, gaps in enumerate(gapss):
            e = SingleTypeAgentEvent()
            cnt: Dict[int, Optional[int]] = {g.id: None for g in gaps}  # Gap.state (getGapMap: last wins)
            if i == 0:
                prev_vis = _d2i(self.rng.next_double() * 0.5 * n)
                self.visible0 = prev_vis
                e.visible = prev_vis
                e.type = EventType.INIT
                b = self.distribute_one_batch(java_int(n - prev_vis), len(gaps))
                e.to_gaps = [g.id for g in gaps]
                for j, g in enumerate(gaps):
                    cnt[g.id] = b[j]
            else:
                pgaps, pcnt = gapss[i - 1], self.counts[i - 1]
                if len(gaps) > len(pgaps):  # split or appear
                    gm = dict.fromkeys(cnt)
                    for pg in pgaps:
                        to = list(pg.to_gaps)
                        if len(to) == 2:
                            e.type = EventType.SPLIT
                            t = state(pcnt, pg.id, i)
                            b = self.distribute_one_batch(t, 2)
                            gap(cnt, to[0], i)
                            cnt[to[0]] = b[0]
                            gap(cnt, to[1], i)
                            cnt[to[1]] = b[1]
                            e.from_gaps = [pg.id]
                            e.to_gaps = [to[0], to[1]]
                        else:
                            gap(cnt, pg.id, i)
                            cnt[pg.id] = pcnt.get(pg.id)
                        gm.pop(pg.id, None)
                    if len(gm) == 1:
                        e.type = EventType.APPEAR
                        gn = next(iter(gm))
                        b = self.distribute_one_batch(prev_vis, 2)
                        cnt[gn] = b[0]
                        prev_vis = b[1]
                        e.visible = prev_vis
                        e.to_gaps = [gn]
                        e.moved = b[0]
                    else:
                        e.visible = prev_vis
                else:  # merge or disappear
                    for pg in pgaps:
                        to = list(pg.to_gaps)
                        if len(to) == 1:
                            t = state(pcnt, pg.id, i)
                            gap(cnt, to[0], i)
                            if cnt[to[0]] is None:
                                e.type = EventType.MERGE
                                e.visible = prev_vis
                                cnt[to[0]] = t
                                e.to_gaps = [to[0]]
                                e.from_gaps = [pg.id, 0]
                            else:
                                # B7: the original adds this parent's count to itself
                                cnt[to[0]] = java_int(t + t if self.merge_bug else cnt[to[0]] + t)  # type: ignore[operator]
                                if e.from_gaps is None:
                                    raise OracleError(f"distributeAgents, set {i}: fromGap is null",
                                                      java_exception=_NPE)
                                if len(e.from_gaps) < 2:  # a disappearance in between set fromGap = {id}
                                    raise OracleError(f"distributeAgents, set {i}: fromGap has length "
                                                      f"{len(e.from_gaps)} (index 1 out of bounds)",
                                                      java_exception=_AIOOBE)
                                e.from_gaps[1] = pg.id
                        elif pg.id not in cnt:  # disappear
                            e.type = EventType.DISAPPEAR
                            ng = state(pcnt, pg.id, i)
                            prev_vis = java_int(prev_vis + ng)
                            e.visible = prev_vis
                            e.from_gaps = [pg.id]
                            e.moved = ng
                        else:
                            cnt[pg.id] = pcnt.get(pg.id)
            self.counts.append(cnt)
            if not gaps:
                raise OracleError(f"distributeAgents: gap set {i} is empty (gapss[{i}][0])", java_exception=_AIOOBE)
            e.time = gaps[0].relative_time
            self._by_time[_time_key(e.time)] = e  # TreeMap.put: equal times overwrite
        self.events = [self._by_time[k] for k in sorted(self._by_time)]
        for c in self.counts:
            self.label_counts.update(c)

    # -- compat="fixed" -----------------------------------------------------------------
    def _distribute_fixed(self, history: Sequence[Sequence[Gap]]) -> None:
        """The same draws, driven by the component events of every transition (any number per transition)."""
        if not history:
            return
        times = _history_times(history)
        n = self.n_targets
        ids0 = [g.id for g in history[0]]
        vis = _d2i(self.rng.next_double() * 0.5 * n)
        if not ids0:  # nothing hidden: every target is visible
            vis = n
        self.visible0 = vis
        b = self.distribute_one_batch(n - vis, len(ids0))
        cnt: Dict[int, int] = dict(zip(ids0, b))
        lc = self.label_counts
        lc.update(cnt)
        self._add(SingleTypeAgentEvent(EventType.INIT, times[0], None, ids0, vis, 0))
        self.counts.append(dict(cnt))
        for k, evs in enumerate(transition_events(history)):
            t = times[k + 1]
            for ev in evs:
                if isinstance(ev, Split):
                    b = self.distribute_one_batch(cnt.pop(ev.s), 2)
                    cnt[ev.a] = lc[ev.a] = b[0]
                    cnt[ev.b] = lc[ev.b] = b[1]
                    self._add(SingleTypeAgentEvent(EventType.SPLIT, t, [ev.s], [ev.a, ev.b], vis, 0))
                elif isinstance(ev, Merge):
                    ca, cb = cnt.pop(ev.a), cnt.pop(ev.b)
                    cnt[ev.s] = lc[ev.s] = cb + cb if self.merge_bug else ca + cb
                    self._add(SingleTypeAgentEvent(EventType.MERGE, t, [ev.a, ev.b], [ev.s], vis, 0))
                elif isinstance(ev, Appear):
                    b = self.distribute_one_batch(vis, 2)
                    cnt[ev.s] = lc[ev.s] = b[0]
                    vis = b[1]
                    self._add(SingleTypeAgentEvent(EventType.APPEAR, t, None, [ev.s], vis, b[0]))
                elif isinstance(ev, Disappear):
                    m = cnt.pop(ev.s)
                    vis += m
                    self._add(SingleTypeAgentEvent(EventType.DISAPPEAR, t, [ev.s], None, vis, m))
                else:  # pragma: no cover - transition_events only yields component events
                    raise OracleError(f"unexpected event {ev!r}")
            self.counts.append({g.id: cnt[g.id] for g in history[k + 1]})

    def _add(self, e: SingleTypeAgentEvent) -> None:
        self.events.append(e)
        self._by_time[_time_key(e.time)] = e

    @classmethod
    def from_sequence(cls, seq: ShadowSequence, n_targets: int, visible0: int,
                      times: Optional[Sequence[float]] = None) -> "SingleTypeAgentOracle":
        """An oracle that reports given observations instead of drawing them (``compat="fixed"``).

        ``seq`` lists the initial shadows and the events; appear / disappear counts must be exact
        (``lo == hi``).  ``visible0`` targets are visible at the start, so ``N - visible0`` are hidden.
        ``times`` (one per event) default to ``(i + 1) / (len(events) + 1)``.  There is no ground truth
        (:attr:`counts` is empty).  Useful to replay a recorded run or a crafted case through
        :func:`~shadowinfo.polygon.stagent.derive_shadow_info_state`, :func:`~shadowinfo.polygon.stagent.java_bounds`
        and :func:`~shadowinfo.polygon.stagent.exact_bounds`.
        """
        o = cls(n_targets, 0, merge_bug=False, compat="fixed")
        o.visible0 = vis = visible0
        evs = list(seq.events)
        ts = list(times) if times is not None else [(i + 1) / (len(evs) + 1) for i in range(len(evs))]
        if len(ts) != len(evs):
            raise ValueError("one time per event is needed")
        o._add(SingleTypeAgentEvent(EventType.INIT, 0.0, None, list(seq.initial), vis, 0))
        for ev, t in zip(evs, ts):
            if isinstance(ev, (Appear, Disappear)):
                if ev.lo != ev.hi:
                    raise ValueError(f"{ev!r}: the oracle's counts are exact (lo == hi)")
                k = int(ev.lo)
                if isinstance(ev, Appear):
                    vis -= k
                    o._add(SingleTypeAgentEvent(EventType.APPEAR, t, None, [ev.s], vis, k))
                else:
                    vis += k
                    o._add(SingleTypeAgentEvent(EventType.DISAPPEAR, t, [ev.s], None, vis, k))
            elif isinstance(ev, Split):
                o._add(SingleTypeAgentEvent(EventType.SPLIT, t, [ev.s], [ev.a, ev.b], vis, 0))
            elif isinstance(ev, Merge):
                o._add(SingleTypeAgentEvent(EventType.MERGE, t, [ev.a, ev.b], [ev.s], vis, 0))
            else:
                raise ValueError(f"{ev!r}: only component events (the Java code has no FOV events, B13)")
        return o

    # -- conversions -------------------------------------------------------------------
    def shadow_events(self) -> List[Event]:
        """The observations after INIT as :mod:`shadowinfo.events` events, in order (see
        :meth:`SingleTypeAgentEvent.shadow_event`).  Raises :class:`OracleError` on an untyped event."""
        out: List[Event] = []
        for e in self.events[1:]:
            se = e.shadow_event()
            if se is None:
                raise OracleError(f"event at t={e.time} has no type (B9)", java_exception=_NPE)
            out.append(se)
        return out

    def initial_shadows(self) -> List[int]:
        """The shadows of the INIT event."""
        if not self.events or self.events[0].type is not EventType.INIT:
            raise OracleError("the oracle has no INIT event (call initialize first)")
        return list(self.events[0].to_gaps or [])

    def sequence(self) -> ShadowSequence:
        """The observations as a :class:`~shadowinfo.events.ShadowSequence`.

        Initial shadows get ``[0, inf)`` each; the known hidden total goes to the filter separately as
        ``total=(hidden_total, hidden_total)`` (:func:`~shadowinfo.polygon.stagent.exact_bounds`), which is
        what the Java filter's pooled left vertex 0 encodes.  Appear/disappear counts are exact.
        """
        return ShadowSequence({g: (0, math.inf) for g in self.initial_shadows()}, self.shadow_events())

    def event_strings(self) -> List[str]:
        """``event.toString()`` of every event (the last block of the ``ProjectPanel`` output)."""
        return [str(e) for e in self.events]

    def __repr__(self) -> str:
        return (f"SingleTypeAgentOracle(n_targets={self.n_targets}, compat={self.compat!r}, "
                f"merge_bug={self.merge_bug}, events={len(self.events)})")
