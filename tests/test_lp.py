"""LP baseline on the shadow sequence (Sec. IV, Eqs. (2)–(4))."""

from __future__ import annotations

import math

import pytest

from shadowinfo.events import Appear, Disappear, Enter, Exit, Merge, ShadowSequence, Split
from shadowinfo.lp import ShadowLP, lp_all_bounds, lp_bounds, lp_feasible
from shadowinfo.maxflow import InfeasibleError


def test_split_merge_disappear():
    seq = ShadowSequence({1: (2, 4), 2: (1, 1)}, [Split(1, 3, 4), Disappear(3, 1, 1), Merge(4, 2, 5)])
    assert lp_all_bounds(seq) == {5: (2, 4)}
    assert lp_bounds(seq, [1], at="start") == (2, 4)


def test_fov_events_as_components():
    seq = ShadowSequence({1: (0, 5)}, [Exit(1, 3), Enter(1, 1)])
    assert lp_all_bounds(seq) == {1: (1, 3)}
    assert lp_bounds(seq, [1], at="start") == (3, 5)


def test_unbounded_and_total():
    seq = ShadowSequence({1: (0, math.inf), 2: (0, math.inf)}, [Appear(3, 1, 2)])
    assert lp_all_bounds(seq) == {1: (0, math.inf), 2: (0, math.inf), 3: (1, 2)}
    assert lp_all_bounds(seq, total=(3, 3)) == {1: (0, 3), 2: (0, 3), 3: (1, 2)}
    assert lp_bounds(seq, [1, 2], total=(3, 3)) == (3, 3)


def test_infeasible():
    seq = ShadowSequence({1: (2, 2)}, [Disappear(1, 0, 1)])
    assert not lp_feasible(seq)
    with pytest.raises(InfeasibleError):
        lp_bounds(seq, [1], at="start")
    assert not lp_feasible(ShadowSequence({1: (0, 2)}, [Exit(1, 3)]))
    assert not lp_feasible(ShadowSequence({1: (0, 2)}, []), total=(3, 4))


def test_variable_count():
    seq = ShadowSequence({1: (0, 1)}, [Enter(1), Split(1, 2, 3)])
    lp = ShadowLP(seq)
    assert len(lp.lo) == 5 and len(lp.eq) == 2
