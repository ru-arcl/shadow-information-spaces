"""FOV batching (Sec. V-D) and FOV -> component conversion."""

from __future__ import annotations

from shadowinfo.events import Appear, Disappear, Enter, Exit, Merge, ShadowSequence, Split
from shadowinfo.fov import FovBatch, batch_fov, fov_to_component
from shadowinfo.lp import lp_all_bounds


def _batch(*events):
    b = FovBatch()
    for e in events:
        b.add(e)
    return b


def test_four_cases():
    assert _batch().batch_events(1) == []
    assert _batch(Exit(1), Enter(1)).batch_events(1) == [Exit(1, 1), Enter(1, 1)]
    assert _batch(Enter(1), Enter(1, 2)).batch_events(1) == [Enter(1, 3)]
    assert _batch(Exit(1, 2), Enter(1)).batch_events(1) == [Exit(1, 2), Enter(1, 1)]
    assert _batch(Enter(1), Exit(1, 3)).batch_events(1) == [Exit(1, 2)]
    b = _batch(Exit(1), Exit(1), Enter(1), Enter(1))
    assert (b.d_min, b.d_tot) == (-2, 0)
    assert b.batch_events(1) == [Exit(1, 2), Enter(1, 2)]


def test_batch_fov_flushes_before_component_events():
    ev = [Enter(1), Exit(1), Exit(2), Merge(1, 2, 3), Enter(3), Exit(3, 2), Split(3, 4, 5), Enter(4)]
    assert batch_fov(ev) == [Exit(2, 1), Merge(1, 2, 3), Exit(3, 1), Split(3, 4, 5), Enter(4, 1)]


def test_batched_and_converted_sequences_preserve_lp_bounds(instances):
    for inst in instances[:120]:
        seq = inst.seq
        expect = lp_all_bounds(seq, inst.total)
        batched = ShadowSequence(seq.initial, batch_fov(seq.events))
        assert lp_all_bounds(batched, inst.total) == expect
        conv, current = fov_to_component(seq)
        assert not any(isinstance(e, (Enter, Exit)) for e in conv.events)
        got = lp_all_bounds(conv, inst.total)
        assert {s: got[current[s]] for s in expect} == expect


def test_fov_to_component_shape():
    seq = ShadowSequence({1: (0, 3)}, [Enter(1, 2), Exit(1, 1), Appear(9), Disappear(1, 0, 4)])
    conv, current = fov_to_component(seq)
    assert conv.events == [Appear(10, 2, 2), Merge(1, 10, 11), Split(11, 12, 13), Disappear(13, 1, 1),
                           Appear(9, 0, 0), Disappear(12, 0, 4)]
    assert current[1] == 12
