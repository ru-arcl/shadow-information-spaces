"""Probabilistic shadow information spaces (T-RO 2012, Sec. VI; ICRA 2010).

The state is a joint probability mass function ``P(s_1, ..., s_n)`` over the
number of targets in each alive shadow, stored as an ordered label tuple plus a
table ``{(x_1, ..., x_n): p}`` (Sec. VI-A).  :class:`ExactFilter` implements
Algorithm 1 (PROCESSPROBABILITYMASS) and Algorithm 2 (PROCESSFOVEVENT) without
approximation; :class:`TruncatedFilter` adds the truncation heuristics TR,
RT and RT-LA of Sec. VI-E / VII-B; :func:`monte_carlo` is the rejection
sampling baseline of Sec. VI-E.

The observation queue ``Q`` of Algorithm 1 mixes component events from
:mod:`shadowinfo.events`, noisy FOV observations (:class:`FovObservation`) and,
optionally, noise-free FOV events (:class:`~shadowinfo.events.Enter`,
:class:`~shadowinfo.events.Exit`).  Count distributions for appear/disappear
events and per-event split probabilities are carried by the subclasses
:class:`ProbAppear`, :class:`ProbDisappear` and :class:`ProbSplit`, which remain
valid members of a :class:`~shadowinfo.events.ShadowSequence`.

All arithmetic is generic: with :class:`fractions.Fraction` inputs the filter
computes exact rationals (e.g. Table III ends at exactly 1/13, 2/13, 10/13).
"""

from __future__ import annotations

import heapq
import itertools
import math
from dataclasses import dataclass, field
from fractions import Fraction
from typing import (
    Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union,
)

from .events import (
    Appear, Disappear, Enter, Event, Exit, Merge, Split, event_from_dict, event_to_dict,
)
from .rng import Rng

Key = Tuple[int, ...]
Table = Dict[Key, float]
ObsModel = Mapping[str, Mapping[str, float]]
SplitRule = Callable[[int], Mapping[Tuple[int, int], float]]

FOV_EVENTS = ("enter", "exit", "null")
_DELTA = {"enter": 1, "null": 0, "exit": -1}


class InconsistentObservationError(ValueError):
    """No probability mass is consistent with the observations (Sec. VI-B)."""


class TruncationFailure(InconsistentObservationError):
    """A truncated filter was left with no mass after it had truncated.

    This is the "failure" outcome of Tables IV and V (Sec. VII-B).  It is
    raised whenever the mass runs out after at least one truncation, so it
    does not prove that truncation was the cause: observations that are
    inconsistent even for :class:`ExactFilter` raise it too once a truncation
    has happened.  Run the exact filter to tell the two apart.
    """


# -- observations ------------------------------------------------------------

@dataclass(frozen=True)
class FovObservation:
    """Noisy FOV observation ``y`` in ``{"enter", "exit", "null"}`` on shadow ``s``.

    The true event is drawn from ``P(e | y)`` (Sec. VI-A, assumption 4).
    """

    s: int
    y: str

    def __post_init__(self) -> None:
        if self.y not in FOV_EVENTS:
            raise ValueError(f"FOV observation must be one of {FOV_EVENTS}, got {self.y!r}")


def _dist_tuple(dist: Mapping[int, float]) -> Tuple[Tuple[int, float], ...]:
    items = tuple(sorted((int(n), p) for n, p in dist.items() if p))
    if not items or any(n < 0 or p < 0 for n, p in items):
        raise ValueError(f"bad count distribution {dict(dist)}")
    return items


@dataclass(frozen=True)
class ProbAppear(Appear):
    """Appear event whose count follows ``dist`` (``P(s_e = n_e)``, Table II).

    Build with :meth:`of`; ``lo``/``hi`` are the support bounds so the event is
    also usable by the nondeterministic filter.  Weights that do not sum to
    one are normalised by the filter.
    """

    dist: Tuple[Tuple[int, float], ...] = ()

    @classmethod
    def of(cls, s: int, dist: Mapping[int, float]) -> "ProbAppear":
        d = _dist_tuple(dist)
        return cls(s, d[0][0], d[-1][0], d)


@dataclass(frozen=True)
class ProbDisappear(Disappear):
    """Disappear event revealing a count distributed as ``dist`` (``P(s_v = n_v)``)."""

    dist: Tuple[Tuple[int, float], ...] = ()

    @classmethod
    def of(cls, s: int, dist: Mapping[int, float]) -> "ProbDisappear":
        d = _dist_tuple(dist)
        return cls(s, d[0][0], d[-1][0], d)


@dataclass(frozen=True)
class ProbSplit(Split):
    """Split event with its own binomial parameter ``p`` (probability that each
    target goes to child ``a``), e.g. proportional to the children's areas."""

    p: float = 0.5


Observation = Union[Event, FovObservation]


# -- models ------------------------------------------------------------------

class BinomialSplit:
    """Split rule: each target independently enters the first child with
    probability ``p`` (the rule of the Fig. 12 example uses ``p = 0.5``)."""

    def __init__(self, p: float = 0.5) -> None:
        if not 0 <= p <= 1:
            raise ValueError("p must lie in [0, 1]")
        self.p = p
        self._cache: Dict[int, Dict[Tuple[int, int], float]] = {}

    def __call__(self, n: int) -> Mapping[Tuple[int, int], float]:
        out = self._cache.get(n)
        if out is None:
            p, q = self.p, 1 - self.p
            try:
                out = {(a, n - a): math.comb(n, a) * p ** a * q ** (n - a) for a in range(n + 1)}
            except OverflowError:  # float p and comb(n, a) > 1.8e308, i.e. n >= 1030
                out = self._log_pmf(n)
            out = {k: v for k, v in out.items() if v}
            self._cache[n] = out
        return out

    def _log_pmf(self, n: int) -> Dict[Tuple[int, int], float]:
        """The pmf computed in log space (float), renormalised to sum to one."""
        p = float(self.p)
        if p in (0.0, 1.0):
            return {(n, 0) if p else (0, n): 1.0}
        lp, lq, c = math.log(p), math.log1p(-p), math.lgamma(n + 1)
        out = {(a, n - a): math.exp(c - math.lgamma(a + 1) - math.lgamma(n - a + 1)
                                    + a * lp + (n - a) * lq) for a in range(n + 1)}
        z = math.fsum(out.values())
        return {k: v / z for k, v in out.items()}

    def __repr__(self) -> str:
        return f"BinomialSplit({self.p!r})"


def perfect_obs_model() -> Dict[str, Dict[str, float]]:
    """``P(e | y) = 1`` iff ``e == y``: noise-free FOV sensing."""
    return {y: {e: int(e == y) for e in FOV_EVENTS} for y in FOV_EVENTS}


def symmetric_obs_model(p: float = 0.9, null: bool = False) -> Dict[str, Dict[str, float]]:
    """Observation model with true-positive rate ``p`` (Sec. VI-D uses ``p = 0.9``).

    Without null events (the paper's example) an enter observation is an exit
    with probability ``1 - p`` and vice versa.  With ``null=True`` the error
    mass ``1 - p`` is shared equally by the two other events.
    """
    if not null:
        return {
            "enter": {"enter": p, "exit": 1 - p, "null": 0},
            "exit": {"enter": 1 - p, "exit": p, "null": 0},
            "null": {"enter": 0, "exit": 0, "null": 1},
        }
    q = (1 - p) / 2
    return {y: {e: (p if e == y else q) for e in FOV_EVENTS} for y in FOV_EVENTS}


_SUM_TOL = 1e-9


def _normalized(weights: Dict, what: str) -> Dict:
    """``weights`` rescaled to sum to one.  Negative weights raise
    :class:`ValueError`; all-zero weights and weights that already sum to one
    up to rounding (``_SUM_TOL``) are returned unchanged."""
    if any(w < 0 for w in weights.values()):
        raise ValueError(f"{what}: negative weight in {weights}")
    z = sum(weights.values(), 0)
    if z and abs(z - 1) > _SUM_TOL:
        return {k: w / z for k, w in weights.items()}
    return weights


def obs_model_from_dict(d: Mapping[str, Mapping[str, float]]) -> Dict[str, Dict[str, float]]:
    """Parse ``obs_model[y][e] = P(e | y)`` (fixture/JS format), filling zeros.

    Unknown ``y``/``e`` keys and negative weights raise :class:`ValueError`.
    A row that does not sum to one is normalised.  An omitted or all-zero row
    means that observation ``y`` is impossible: observing it raises
    :class:`InconsistentObservationError`.
    """
    unknown = [repr(y) for y in d if y not in FOV_EVENTS]
    unknown += [f"{y}/{e!r}" for y, row in d.items() if y in FOV_EVENTS for e in (row or {})
                if e not in FOV_EVENTS]
    if unknown:
        raise ValueError(f"observation model: unknown FOV event(s) {', '.join(unknown)}; "
                         f"expected {FOV_EVENTS}")
    return {y: _normalized({e: (d.get(y) or {}).get(e, 0) for e in FOV_EVENTS},
                           f"observation model row {y!r}") for y in FOV_EVENTS}


def _split_dist(rule: SplitRule, n: int) -> Mapping[Tuple[int, int], float]:
    """``rule(n)``; a custom rule's output is checked and normalised."""
    dist = rule(n)
    if isinstance(rule, BinomialSplit):
        return dist
    if any(na < 0 or nb < 0 or na + nb != n for na, nb in dist):
        raise ValueError(f"split rule: {dict(dist)} is not a split of {n} targets")
    if not any(dist.values()):
        raise ValueError(f"split rule: no probability mass for {n} targets")
    return _normalized(dict(dist), f"split rule output for {n} targets")


def split_rule_from_dict(d: Mapping) -> BinomialSplit:
    if d.get("type", "binomial") != "binomial":
        raise ValueError(f"unknown split rule {d!r}")
    return BinomialSplit(d.get("p", 0.5))


# -- (de)serialisation -------------------------------------------------------

_FIELDS = {"appear": ("s", "lo", "hi"), "disappear": ("s", "lo", "hi"), "split": ("s", "a", "b"),
           "merge": ("a", "b", "s"), "enter": ("s", "k"), "exit": ("s", "k")}


def observation_from_dict(d: Mapping) -> Observation:
    """Parse one queue item; extends :func:`events.event_from_dict` with
    ``{"type": "fov", "s", "y"}``, ``"dist": {count: p}`` on appear/disappear
    and ``"p"`` on split.  Unknown keys (annotations) are ignored."""
    t = d["type"]
    if t == "fov":
        return FovObservation(int(d["s"]), d["y"])
    if t in ("appear", "disappear") and "dist" in d:
        dist = {int(n): p for n, p in d["dist"].items()}
        return (ProbAppear if t == "appear" else ProbDisappear).of(int(d["s"]), dist)
    if t == "split" and "p" in d:
        return ProbSplit(int(d["s"]), int(d["a"]), int(d["b"]), d["p"])
    return event_from_dict({"type": t, **{k: d[k] for k in _FIELDS[t] if k in d}})


def _json_num(x):
    """JSON has no rationals: :class:`~fractions.Fraction` is written as float."""
    return float(x) if isinstance(x, Fraction) else x


def observation_to_dict(o: Observation) -> dict:
    """Inverse of :func:`observation_from_dict`; Fraction weights become floats."""
    if isinstance(o, FovObservation):
        return {"type": "fov", "s": o.s, "y": o.y}
    if isinstance(o, (ProbAppear, ProbDisappear)):
        base = event_to_dict(Appear(o.s, o.lo, o.hi) if isinstance(o, Appear)
                             else Disappear(o.s, o.lo, o.hi))
        return {**base, "dist": {str(n): _json_num(p) for n, p in o.dist}}
    if isinstance(o, ProbSplit):
        return {"type": "split", "s": o.s, "a": o.a, "b": o.b, "p": _json_num(o.p)}
    return event_to_dict(o)


@dataclass
class ProbSequence:
    """Initial joint pmf plus the observation queue ``Q`` of Algorithm 1.

    JSON format (shared with the fixtures and the JS demo)::

        {"labels": [1, 2], "joint": [[[2, 2], 1.0]],
         "events": [{"type": "fov", "s": 1, "y": "exit"}, ...]}

    Repeated joint keys are summed on input.  JSON is float valued:
    :class:`~fractions.Fraction` masses and weights are written as floats.
    """

    labels: Tuple[int, ...]
    joint: Table
    observations: List[Observation] = field(default_factory=list)

    @classmethod
    def from_dict(cls, d: Mapping) -> "ProbSequence":
        labels = tuple(int(s) for s in d["labels"])
        joint: Table = {}
        for k, p in d["joint"]:
            k = tuple(int(x) for x in k)
            joint[k] = joint.get(k, 0) + p
        return cls(labels, joint, [observation_from_dict(e) for e in d["events"]])

    def to_dict(self) -> dict:
        return {"labels": list(self.labels),
                "joint": [[list(k), _json_num(p)] for k, p in self.joint.items()],
                "events": [observation_to_dict(o) for o in self.observations]}

    @classmethod
    def from_counts(cls, counts: Mapping[int, int],
                    observations: Iterable[Observation] = ()) -> "ProbSequence":
        labels = tuple(counts)
        return cls(labels, {tuple(int(counts[s]) for s in labels): 1}, list(observations))


# -- the joint pmf -----------------------------------------------------------

class JointPmf:
    """A joint pmf over the target counts of the shadows in ``labels``."""

    def __init__(self, labels: Sequence[int], table: Mapping[Key, float]) -> None:
        self.labels: Tuple[int, ...] = tuple(labels)
        self.table: Table = {tuple(k): p for k, p in table.items() if p}
        if any(len(k) != len(self.labels) for k in self.table):
            raise ValueError("joint table keys must match the labels")
        if any(p < 0 for p in self.table.values()):
            raise ValueError("negative probability mass")

    @property
    def entries(self) -> int:
        """Number of stored probability mass entries."""
        return len(self.table)

    def total(self) -> float:
        return sum(self.table.values())

    def index(self, s: int) -> int:
        try:
            return self.labels.index(s)
        except ValueError:
            raise KeyError(f"shadow {s} is not alive") from None

    def marginal(self, s: int) -> Dict[int, float]:
        """``P(s = x)`` for each count ``x`` with positive mass."""
        i = self.index(s)
        m: Dict[int, float] = {}
        for k, p in self.table.items():
            m[k[i]] = m.get(k[i], 0) + p
        return dict(sorted(m.items()))

    def expected_counts(self) -> Dict[int, float]:
        """``E[s]`` for every alive shadow."""
        out = {s: 0 for s in self.labels}
        for k, p in self.table.items():
            for s, x in zip(self.labels, k):
                out[s] += x * p
        return out

    def probability(self, assignment: Mapping[int, int]) -> float:
        """Probability that every shadow in ``assignment`` holds the given count."""
        idx = [(self.index(s), x) for s, x in assignment.items()]
        return sum((p for k, p in self.table.items() if all(k[i] == x for i, x in idx)), 0)

    def as_dict(self) -> Dict[Key, float]:
        return dict(sorted(self.table.items()))

    def reordered(self, labels: Sequence[int]) -> "JointPmf":
        """The same pmf with its key columns in the order of ``labels``, a
        permutation of :attr:`labels` (e.g. to compare a :func:`monte_carlo`
        result, whose labels are sorted, with an :class:`ExactFilter`)."""
        labels = tuple(labels)
        if sorted(labels) != sorted(self.labels):
            raise KeyError(f"labels {labels} are not a permutation of {self.labels}")
        idx = [self.index(s) for s in labels]
        return JointPmf(labels, {tuple(k[i] for i in idx): p for k, p in self.table.items()})

    def _normalize(self) -> None:
        z = self.total()
        if not z:
            raise InconsistentObservationError("no probability mass left")
        if z != 1:
            self.table = {k: v for k, v in ((k, p / z) for k, p in self.table.items()) if v}

    def __repr__(self) -> str:
        return f"{type(self).__name__}(labels={self.labels}, entries={self.entries})"


class ExactFilter(JointPmf):
    """Exact Bayesian filter over joint target counts (Algorithms 1 and 2).

    ``split_rule(n)`` gives ``{(n_a, n_b): prob}`` (default binomial, p = 0.5);
    a :class:`ProbSplit` overrides it with its own binomial ``p``.
    ``obs_model[y][e] = P(e | y)`` (default: noise-free).
    """

    def __init__(self, labels: Sequence[int], joint: Mapping[Key, float],
                 split_rule: Optional[SplitRule] = None,
                 obs_model: Optional[ObsModel] = None) -> None:
        super().__init__(labels, joint)
        self._normalize()
        self.split_rule: SplitRule = split_rule or BinomialSplit(0.5)
        self.obs_model = obs_model_from_dict(obs_model) if obs_model else perfect_obs_model()
        self.peak_entries = self.entries
        self._split_rules: Dict[float, BinomialSplit] = {}

    @classmethod
    def from_counts(cls, counts: Mapping[int, int], **kw) -> "ExactFilter":
        """Known initial counts: ``P(s_1 = c_1, ..., s_n = c_n) = 1``."""
        labels = tuple(counts)
        return cls(labels, {tuple(int(counts[s]) for s in labels): 1}, **kw)

    @classmethod
    def from_independent(cls, dists: Mapping[int, Mapping[int, float]], **kw) -> "ExactFilter":
        """Independent initial marginals ``{label: {count: prob}}``."""
        labels = tuple(dists)
        table: Table = {}
        for combo in itertools.product(*(dists[s].items() for s in labels)):
            p = 1
            for _, q in combo:
                p *= q
            if p:
                table[tuple(n for n, _ in combo)] = p
        return cls(labels, table, **kw)

    @classmethod
    def from_sequence(cls, seq: ProbSequence, **kw) -> "ExactFilter":
        return cls(seq.labels, seq.joint, **kw)

    def copy(self) -> "ExactFilter":
        other = object.__new__(type(self))
        other.__dict__.update(self.__dict__)
        other.table = dict(self.table)
        return other

    # -- Algorithm 1 --------------------------------------------------------
    def run(self, observations: Sequence[Observation]) -> "ExactFilter":
        """Process the whole queue ``Q`` (Algorithm 1) and return ``self``."""
        observations = list(observations)
        for t, o in enumerate(observations):
            self.apply(o, observations[t + 1:])
        return self

    def apply(self, o: Observation, upcoming: Optional[Sequence[Observation]] = None) -> None:
        """Process one queue item.  ``upcoming`` (the rest of the queue) is
        only used by the RT-LA truncation heuristic."""
        if isinstance(o, FovObservation):
            self._fov(o.s, o.y)
        elif isinstance(o, Split):
            self._split(o)
        elif isinstance(o, Merge):
            self._merge(o.a, o.b, o.s)
        elif isinstance(o, Appear):
            self._appear(o)
        elif isinstance(o, Disappear):
            self._disappear(o)
        elif isinstance(o, Enter):
            self._shift(o.s, o.k)
        elif isinstance(o, Exit):
            self._shift(o.s, -o.k)
        else:
            raise TypeError(f"unknown observation {o!r}")
        self.peak_entries = max(self.peak_entries, self.entries)
        self._after_update(upcoming)

    def _after_update(self, upcoming: Optional[Sequence[Observation]]) -> None:
        pass

    # -- component events (Sec. VI-B) ----------------------------------------
    def _appear(self, e: Appear) -> None:
        if e.s in self.labels:
            raise KeyError(f"shadow {e.s} already alive")
        if isinstance(e, ProbAppear):
            dist = e.dist
        else:
            if e.hi == math.inf:
                raise ValueError("an appear event needs a finite upper bound or a dist")
            w = Fraction(1, int(e.hi) - e.lo + 1) if e.hi > e.lo else 1
            dist = tuple((n, w) for n in range(e.lo, int(e.hi) + 1))
        self.table = {k + (n,): v for k, p in self.table.items() for n, w in dist
                      for v in (p * w,) if v}
        self.labels += (e.s,)
        self._normalize()

    def _disappear(self, e: Disappear) -> None:
        i = self.index(e.s)
        if isinstance(e, ProbDisappear):
            dist = dict(e.dist)
            like = dist.get
        else:
            lo, hi = e.lo, e.hi
            like = lambda x: 1 if lo <= x <= hi else 0  # noqa: E731
        out: Table = {}
        for k, p in self.table.items():
            w = like(k[i])
            v = p * w if w else 0
            if v:
                nk = k[:i] + k[i + 1:]
                out[nk] = out.get(nk, 0) + v
        self.table = out
        self.labels = self.labels[:i] + self.labels[i + 1:]
        self._normalize_or_fail(f"disappear of shadow {e.s}")

    def _split(self, e: Split) -> None:
        i = self.index(e.s)
        rest_labels = self.labels[:i] + self.labels[i + 1:]
        if e.a == e.b or e.a in rest_labels or e.b in rest_labels:
            raise KeyError(f"split of shadow {e.s}: children {e.a}, {e.b} must be distinct "
                           f"and not alive")
        if isinstance(e, ProbSplit):
            rule = self._split_rules.get(e.p)
            if rule is None:
                rule = self._split_rules[e.p] = BinomialSplit(e.p)
        else:
            rule = self.split_rule
        dists: Dict[int, Mapping[Tuple[int, int], float]] = {}
        out: Table = {}
        for k, p in self.table.items():
            rest = k[:i] + k[i + 1:]
            dist = dists.get(k[i])
            if dist is None:
                dist = dists[k[i]] = _split_dist(rule, k[i])
            for (na, nb), q in dist.items():
                v = p * q
                if v:
                    nk = rest + (na, nb)
                    out[nk] = out.get(nk, 0) + v
        self.table = out
        self.labels = rest_labels + (e.a, e.b)

    def _merge(self, a: int, b: int, s: int) -> None:
        i, j = self.index(a), self.index(b)
        if i == j:
            raise ValueError("cannot merge a shadow with itself")
        if s in self.labels and s not in (a, b):
            raise KeyError(f"merge into shadow {s}: already alive")
        keep = [q for q in range(len(self.labels)) if q not in (i, j)]
        out: Table = {}
        for k, p in self.table.items():
            nk = tuple(k[q] for q in keep) + (k[i] + k[j],)
            out[nk] = out.get(nk, 0) + p
        self.table = out
        self.labels = tuple(self.labels[q] for q in keep) + (s,)

    # -- FOV events (Sec. VI-C) -----------------------------------------------
    def _fov(self, s: int, y: str) -> None:
        """Algorithm 2.  When ``x_s = 0`` the exit branch is impossible and the
        enter/null branches of *that entry* are rescaled to its mass ``p_j``;
        an entry with no feasible branch is dropped (global renormalisation).
        An observation whose model row is all zero is impossible."""
        i = self.index(s)
        model = self.obs_model[y]
        branches = [(_DELTA[e], model[e]) for e in FOV_EVENTS if model[e]]
        z0 = sum((w for d, w in branches if d >= 0), 0)
        dropped = False
        out: Table = {}
        for k, p in self.table.items():
            x = k[i]
            if x == 0 and not z0:
                dropped = True
                continue
            for d, w in branches:
                if x == 0:
                    if d < 0:
                        continue
                    w = w / z0
                v = p * w
                if v:
                    nk = k[:i] + (x + d,) + k[i + 1:]
                    out[nk] = out.get(nk, 0) + v
        self.table = out
        if dropped or not out:
            self._normalize_or_fail(f"FOV observation {y!r} on shadow {s}")

    def _shift(self, s: int, k: int) -> None:
        """Noise-free FOV events: ``Enter(s, k)`` adds ``k`` targets to every
        entry; ``Exit(s, k)`` removes ``k``, ruling out entries with ``x_s < k``."""
        i = self.index(s)
        out: Table = {}
        for key, p in self.table.items():
            x = key[i] + k
            if x >= 0:
                out[key[:i] + (x,) + key[i + 1:]] = p
        dropped = len(out) < len(self.table)
        self.table = out
        if dropped:
            self._normalize_or_fail(f"exit of {-k} target(s) from shadow {s}")

    def _normalize_or_fail(self, what: str) -> None:
        if not self.table:
            raise InconsistentObservationError(f"{what}: no probability mass left")
        self._normalize()


# -- truncation heuristics (Sec. VI-E, VII-B) --------------------------------

TRUNCATION_MODES = ("TR", "RT", "RT-LA")


class TruncatedFilter(ExactFilter):
    """:class:`ExactFilter` that keeps at most ``max_entries`` entries.

    After every processed observation (component or FOV) with more than
    ``max_entries`` entries:

    * ``"TR"`` (basic truncation, called RT in ICRA 2010) keeps the entries
      of largest mass ``p_j``;
    * ``"RT"`` (random truncation) keeps the largest ``p_j * u_j`` with
      ``u_j ~ U(0, 1)`` drawn from a seeded :class:`~shadowinfo.rng.Rng`;
    * ``"RT-LA"`` is RT, but skips truncation when a disappear is within the
      next ``lookahead[0]`` queue items or a merge within the next
      ``lookahead[1]`` (paper: 4 and 2).  With ``lookahead_counts_fov=False``
      only component events are counted.  It needs the future queue, so use
      :meth:`run` or pass ``upcoming`` to :meth:`apply`; without it RT-LA
      behaves as RT.

    The kept entries are renormalised.  An observation that leaves no mass
    after a truncation raises :class:`TruncationFailure` (which may also mean
    that the observations are inconsistent outright).
    """

    def __init__(self, labels: Sequence[int], joint: Mapping[Key, float], max_entries: int,
                 mode: str = "TR", seed: int = 1, lookahead: Tuple[int, int] = (4, 2),
                 lookahead_counts_fov: bool = True, **kw) -> None:
        if mode not in TRUNCATION_MODES:
            raise ValueError(f"mode must be one of {TRUNCATION_MODES}")
        if max_entries < 1:
            raise ValueError("max_entries must be positive")
        super().__init__(labels, joint, **kw)
        self.max_entries = max_entries
        self.mode = mode
        self.rng = Rng(seed)
        self.lookahead = lookahead
        self.lookahead_counts_fov = lookahead_counts_fov
        self.truncations = 0
        self.truncated_mass = 0.0

    def copy(self) -> "TruncatedFilter":
        other = super().copy()
        other.rng = Rng()
        other.rng.state = self.rng.state
        return other

    def _skip_for_lookahead(self, upcoming: Optional[Sequence[Observation]]) -> bool:
        if self.mode != "RT-LA" or not upcoming:
            return False
        n_dis, n_mer = self.lookahead
        items: Iterable[Observation] = upcoming
        if not self.lookahead_counts_fov:
            items = (o for o in upcoming if isinstance(o, (Appear, Disappear, Split, Merge)))
        for pos, o in enumerate(itertools.islice(items, max(n_dis, n_mer)), 1):
            if (isinstance(o, Disappear) and pos <= n_dis) or (isinstance(o, Merge) and pos <= n_mer):
                return True
        return False

    def _after_update(self, upcoming: Optional[Sequence[Observation]]) -> None:
        if self.entries <= self.max_entries or self._skip_for_lookahead(upcoming):
            return
        items = sorted(self.table.items())
        if self.mode == "TR":
            kept = heapq.nlargest(self.max_entries, items, key=lambda kp: kp[1])
        else:
            scores = [p * self.rng.random() for _, p in items]
            top = heapq.nlargest(self.max_entries, range(len(items)), key=scores.__getitem__)
            kept = [items[j] for j in top]
        before = self.total()
        self.table = dict(kept)
        self.truncated_mass += float(before - self.total())
        self.truncations += 1
        self._normalize()

    def _normalize_or_fail(self, what: str) -> None:
        if not self.table and self.truncations:
            raise TruncationFailure(f"{what}: no probability mass left after "
                                    f"{self.truncations} truncation(s)")
        super()._normalize_or_fail(what)


# -- Monte Carlo baseline (Sec. VI-E) ----------------------------------------

class MonteCarloResult(JointPmf):
    """Empirical joint pmf of the final counts over ``successes`` accepted
    trials out of ``attempts``."""

    def __init__(self, labels: Sequence[int], table: Mapping[Key, float],
                 successes: int, attempts: int) -> None:
        super().__init__(labels, table)
        self.successes = successes
        self.attempts = attempts


def _sample(rng: Rng, items: Sequence[Tuple[object, float]]) -> object:
    """Inverse-CDF draw from ``[(value, weight), ...]`` (weights need not sum to 1)."""
    total = sum(w for _, w in items)
    u = rng.random() * total
    acc = 0
    for v, w in items:
        acc += w
        if u < acc:
            return v
    return items[-1][0]


def monte_carlo(seq: ProbSequence, trials: int = 1000, seed: int = 1,
                split_rule: Optional[SplitRule] = None, obs_model: Optional[ObsModel] = None,
                max_attempts: Optional[int] = None) -> MonteCarloResult:
    """Rejection-sampling baseline of Sec. VI-E, run until ``trials`` successes.

    Each trial samples initial counts from the joint pmf and propagates them
    through the queue: binomial splits route every target independently
    (other rules sample ``(n_a, n_b)`` from ``split_rule(n)``); a FOV
    observation samples its event from ``P(e | y)`` restricted to the events
    feasible at the current count (matching Algorithm 2); an appear samples
    its count; a disappear samples the revealed count and *discards* the
    trial if it differs from the simulated one.  Trials whose noise-free exit
    or FOV observation is infeasible are discarded too.

    The result's labels are *sorted*, whereas :class:`ExactFilter` keeps them
    in event order; use :meth:`JointPmf.reordered` before comparing tables.
    """
    split_rule = split_rule or BinomialSplit(0.5)
    model = obs_model_from_dict(obs_model) if obs_model else perfect_obs_model()
    max_attempts = max_attempts if max_attempts is not None else 1000 * trials
    rng = Rng(seed)
    init = sorted((k, p) for k, p in seq.joint.items() if p)
    fov_items = {y: [(e, model[y][e]) for e in FOV_EVENTS if model[y][e]] for y in FOV_EVENTS}
    final_labels: Optional[Tuple[int, ...]] = None
    counts: Dict[Key, int] = {}
    successes = attempts = 0
    while successes < trials:
        if attempts >= max_attempts:
            raise InconsistentObservationError(
                f"only {successes} of {attempts} Monte Carlo trials were consistent")
        attempts += 1
        x = dict(zip(seq.labels, _sample(rng, init)))
        ok = True
        for o in seq.observations:
            if isinstance(o, FovObservation):
                n = x[o.s]
                items = [(e, w) for e, w in fov_items[o.y] if n > 0 or e != "exit"]
                if not items:
                    ok = False
                    break
                x[o.s] = n + _DELTA[_sample(rng, items)]
            elif isinstance(o, Split):
                n = x.pop(o.s)
                if o.a == o.b or o.a in x or o.b in x:
                    raise KeyError(f"split of shadow {o.s}: children {o.a}, {o.b} must be "
                                   f"distinct and not alive")
                rule = BinomialSplit(o.p) if isinstance(o, ProbSplit) else split_rule
                if isinstance(rule, BinomialSplit):
                    na = sum(rng.random() < rule.p for _ in range(n))
                else:
                    na = _sample(rng, [(k[0], w) for k, w in _split_dist(rule, n).items()])
                x[o.a], x[o.b] = na, n - na
            elif isinstance(o, Merge):
                if o.a == o.b:
                    raise ValueError("cannot merge a shadow with itself")
                m = x.pop(o.a) + x.pop(o.b)
                if o.s in x:
                    raise KeyError(f"merge into shadow {o.s}: already alive")
                x[o.s] = m
            elif isinstance(o, Appear):
                if o.s in x:
                    raise KeyError(f"shadow {o.s} already alive")
                if isinstance(o, ProbAppear):
                    x[o.s] = _sample(rng, o.dist)
                elif o.hi == math.inf:
                    raise ValueError("an appear event needs a finite upper bound or a dist")
                else:
                    x[o.s] = o.lo + rng.randint(int(o.hi) - o.lo + 1)
            elif isinstance(o, Disappear):
                n = x.pop(o.s)
                if isinstance(o, ProbDisappear):
                    ok = _sample(rng, o.dist) == n
                else:
                    ok = o.lo <= n <= o.hi
                if not ok:
                    break
            elif isinstance(o, (Enter, Exit)):
                x[o.s] += o.k if isinstance(o, Enter) else -o.k
                if x[o.s] < 0:
                    ok = False
                    break
            else:
                raise TypeError(f"unknown observation {o!r}")
        if not ok:
            continue
        if final_labels is None:
            final_labels = tuple(sorted(x))
        key = tuple(x[s] for s in final_labels)
        counts[key] = counts.get(key, 0) + 1
        successes += 1
    table = {k: c / successes for k, c in counts.items()}
    return MonteCarloResult(final_labels or (), table, successes, attempts)


__all__ = [
    "FOV_EVENTS", "TRUNCATION_MODES", "BinomialSplit", "ExactFilter", "FovObservation",
    "InconsistentObservationError", "JointPmf", "MonteCarloResult", "ProbAppear",
    "ProbDisappear", "ProbSequence", "ProbSplit", "TruncatedFilter", "TruncationFailure",
    "monte_carlo", "obs_model_from_dict", "observation_from_dict", "observation_to_dict",
    "perfect_obs_model", "split_rule_from_dict", "symmetric_obs_model",
]
