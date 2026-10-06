"""Live viewer (``python -m shadowinfo.polygon show``): headless checks with the Agg backend."""

from types import SimpleNamespace

import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from shadowinfo.polygon import physical_gaps  # noqa: E402
from shadowinfo.polygon.__main__ import main  # noqa: E402
from shadowinfo.polygon.viewer import Viewer, gap_anchor  # noqa: E402


@pytest.fixture
def viewer():
    fig, ax = plt.subplots()
    yield Viewer(ax, 12)
    plt.close(fig)


def _event(v, x, y, key=None):
    return SimpleNamespace(inaxes=v.ax, xdata=x, ydata=y, key=key)


@pytest.mark.parametrize("map_no", range(1, 15))
def test_show_saves_every_map(tmp_path, map_no):
    out = tmp_path / f"m{map_no}.png"
    assert main(["show", str(map_no), "--save", str(out)]) == 0
    assert out.stat().st_size > 1000


def test_show_rejects_bad_input(tmp_path, capsys):
    assert main(["show", "15"]) == 2
    assert main(["show", "12", "--at", "-50", "-50", "--save", str(tmp_path / "x.png")]) == 2
    assert "outside" in capsys.readouterr().err


def test_mouse_moves_robot_and_updates_pockets(viewer):
    viewer.on_move(_event(viewer, 300, 500))
    assert viewer.q == (300, 500)
    assert len(viewer.pockets) == len(physical_gaps(viewer.poly, (300, 500)))
    viewer.on_move(_event(viewer, -50, -50))  # outside the polygon: ignored
    assert viewer.q == (300, 500)


def test_click_pins_and_unpins(viewer):
    viewer.on_click(_event(viewer, 120, 450))
    assert viewer.pinned and viewer.q == (120, 450)
    viewer.on_move(_event(viewer, 300, 500))
    assert viewer.q == (120, 450)
    viewer.on_click(_event(viewer, 300, 500))
    assert not viewer.pinned and viewer.q == (300, 500)


def test_arrow_keys_cycle_maps(viewer):
    viewer.on_key(_event(viewer, None, None, key="left"))
    assert viewer.map_no == 11
    for _ in range(4):
        viewer.on_key(_event(viewer, None, None, key="right"))
    assert viewer.map_no == 1


def test_pocket_colours_follow_casting_vertex(viewer):
    """Small moves keep each gap's anchor vertex, hence its colour."""
    a = {gap_anchor(g) for g in physical_gaps(viewer.poly, (300, 500))}
    b = {gap_anchor(g) for g in physical_gaps(viewer.poly, (302, 501))}
    assert a == b
