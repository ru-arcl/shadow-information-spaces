"""Incremental bipartite I-state (Fig. 9, Fig. 11(c))."""

from __future__ import annotations

import math

import pytest

from shadowinfo.bipartite import BipartiteIState, Pseudo, vertex_name
from shadowinfo.events import Appear, Disappear, Enter, Exit, InvalidSequenceError, Merge, Split


def test_fig9_operations():
    st = BipartiteIState({1: (1, 2), 2: (0, 1)})
    st.apply(Appear(3, 2, 2))
    assert st.left == {1: (1, 2), 2: (0, 1), 3: (2, 2)}
    assert st.reach[3] == {3}
    st.apply(Split(1, 4, 5))
    assert st.reach[4] == st.reach[5] == {1}
    st.apply(Merge(5, 2, 6))
    assert st.reach[6] == {1, 2}
    st.apply(Disappear(4, 1, 1))
    assert st.disappeared == {4: (1, 1)} and st.disappeared_reach[4] == {1}
    assert sorted(st.alive()) == [3, 6]
    assert st.edges() == [(1, ("gone", 4)), (3, ("alive", 3)), (1, ("alive", 6)), (2, ("alive", 6))]
    assert st.right()[("alive", 6)] == (0, math.inf)


def test_fig11_bipartite_graph(fig11):
    st = BipartiteIState.from_sequence(fig11["seq"])
    b = fig11["bipartite"]
    assert list(st.left) == b["left"]
    assert sorted(v for _, v in st.right()) == b["right"]
    assert sorted((u, v) for u, (_, v) in st.edges()) == sorted(map(tuple, b["edges"]))
    assert 18 in st.left and ("alive", 18) in st.right()


def test_naive_fov_creates_pseudo_vertices():
    st = BipartiteIState({1: (0, 2)}, fov="naive")
    st.apply(Enter(1, 2))
    st.apply(Exit(1, 1))
    p_in, p_out = Pseudo("enter", 1, 1), Pseudo("exit", 1, 2)
    assert st.left[p_in] == (2, 2)
    assert st.reach[1] == {1, p_in}
    assert st.disappeared[p_out] == (1, 1) and st.disappeared_reach[p_out] == {1, p_in}
    assert st.real_left() == [1]
    assert vertex_name(p_in) == "+1#1" and vertex_name(p_out) == "-1#2"


def test_pseudo_labels_never_collide_with_real_labels():
    st = BipartiteIState({1: (0, 2)}, fov="naive")
    st.apply(Enter(1))
    st.apply(Appear(2, 1, 1))
    st.apply(Appear(3, 0, 0))
    assert len(st.left) == 4
    assert Pseudo("enter", 1, 1) != 1 and Pseudo("enter", 1, 1) not in (1, 2, 3)


def test_batch_pending_until_component_event_or_flush():
    st = BipartiteIState({1: (0, 5), 2: (0, 5)})
    for e in (Exit(1), Exit(1), Enter(1), Enter(1), Enter(2)):
        st.apply(e)
    assert len(st.left) == 2 and not st.disappeared
    st.apply(Split(1, 3, 4))
    assert [v for v in st.disappeared] == [Pseudo("exit", 1, 1)]
    assert st.disappeared[Pseudo("exit", 1, 1)] == (2, 2)
    assert st.left[Pseudo("enter", 1, 2)] == (2, 2)
    assert st.reach[3] == {1, Pseudo("enter", 1, 2)}
    assert 2 in st.pending
    d = st.to_dict()
    assert not st.pending
    assert ["+2#3", 1, 1] in d["left"]


def test_invalid_events():
    st = BipartiteIState({1: (0, 1)})
    with pytest.raises(InvalidSequenceError):
        st.apply(Split(7, 8, 9))
    with pytest.raises(InvalidSequenceError):
        st.apply(Appear(1))
    with pytest.raises(InvalidSequenceError):
        st.apply(Exit(5))
    with pytest.raises(InvalidSequenceError):
        st.apply(Disappear(1, 2, 1))
    with pytest.raises(InvalidSequenceError, match="repeated label"):
        st.apply(Split(1, 2, 2))
    st.apply(Disappear(1, 0, 0))
    with pytest.raises(InvalidSequenceError):
        st.apply(Appear(1))
    with pytest.raises(ValueError):
        BipartiteIState({}, fov="lazy")
