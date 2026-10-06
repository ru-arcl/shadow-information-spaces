"""Smoke tests: every script in ``examples/`` runs with small parameters and prints its result."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


def run(script: str, *args: str, timeout: float = 240) -> str:
    p = subprocess.run([sys.executable, str(EXAMPLES / script), *args], capture_output=True, text=True,
                       timeout=timeout, cwd=EXAMPLES.parent)
    assert p.returncode == 0, p.stdout + p.stderr
    return p.stdout


def test_fig11_maxflow():
    out = run("fig11_maxflow.py")
    assert "feasible flow with lower bounds : [10, 24]" in out
    assert "paper's recipe, Eqs. (7)/(8)    : [10, 24]" in out
    assert "LP relaxation, Eqs. (2)-(4)     : [10, 24]" in out
    assert "agree with the paper: yes" in out
    # the geometry behind the figure: map 12 + the ProjectPanel path (polygon front end)
    assert "30 critical points, 14 component events, labels 1-19; initial shadows [1, 2, 3, 4]" in out
    assert "Split(s=1, a=7, b=6)   (Java labels: Split(s=1, a=6, b=7))" in out
    assert "the same 13 events as Fig. 11(b) after swapping labels 6 and 7: yes" in out
    assert "bipartite I-state = Fig. 11(c) (11 edges): yes" in out
    assert "bounds on s19: [10, 24]; final shadows [13, 15, 18, 19]" in out
    assert "reproduces Fig. 11(b)-(c) and the bounds: yes" in out


def test_table3_probabilistic():
    out = run("table3_probabilistic.py", "--trials", "100", "--seeds", "1")
    assert "P(s4=0) =      1/13" in out and "P(s4=2) =     10/13" in out
    assert "y_x, s1  (2 entries)" in out and "merge s1, s3 -> s5  (8 entries)" in out
    assert "P(s4 = 0, 1, 2) = 0, 1/6, 5/6" in out
    assert "100/" in out


def test_truncation_heuristics():
    out = run("truncation_heuristics.py", "--fov", "4", "--max-entries", "200", "2000", "--runs", "2",
              "--trials", "30")
    for name in ("exact", "TR-2000", "RT-2000", "RT-LA-2000", "Monte Carlo"):
        assert f"  {name} " in out
    assert "30/" in out


def test_office_grid():
    out = run("office_grid.py", "--ticks", "40")
    assert "40 ticks" in out and "truth outside the bounds at any tick: 0" in out
    assert "E[count]" in out


@pytest.mark.parametrize("args", [["--mode", "unknown", "--targets", "20"], ["--mode", "evader", "--targets", "1"]])
def test_office_grid_modes(args):
    out = run("office_grid.py", "--ticks", "60", *args)
    assert "truth outside the bounds at any tick: 0" in out


def test_office_grid_plot(tmp_path):
    pytest.importorskip("matplotlib")
    out = run("office_grid.py", "--ticks", "40", "--plot", "--out", str(tmp_path))
    for name in ("office_snapshot.png", "shadow_sequence.png", "bipartite.png"):
        assert (tmp_path / name).stat().st_size > 10_000
        assert name in out


def test_counting_and_pursuit():
    out = run("counting_and_pursuit.py", "--laps", "1", "--every", "98")
    assert "n (counting)" in out and "contaminated" in out
    assert "\n   98 " in out and "\n  196 " in out


def test_readme_figures(tmp_path):
    pytest.importorskip("matplotlib")
    run("figures.py", "--out", str(tmp_path))
    assert sorted(p.name for p in tmp_path.iterdir()) == ["bipartite.png", "events.png", "office_snapshot.png",
                                                          "shadow_sequence.png"]


def test_event_illustrations(tmp_path):
    pytest.importorskip("matplotlib")
    out = run("event_illustrations.py", "--out", str(tmp_path / "events.png"))
    assert "appear at critical point 2, disappear at critical point 15, split at critical point 4, " \
           "merge at critical point 21" in out
    assert (tmp_path / "events.png").stat().st_size > 50_000
    listing = run("event_illustrations.py", "--list")
    assert "Split(s=1, a=6, b=7)" in listing and "Merge(a=11, b=5, s=16)" in listing


def row(out: str, name: str) -> list:
    """The cells of a table row of ``paper_fig15.py`` (``name``, paper, compat=java, default, note)."""
    m = re.search(rf"^  {re.escape(name)} +(\S+) +(\S+) +(\S+)", out, re.M)
    assert m, name
    return list(m.groups())


def test_paper_fig15():
    out = run("paper_fig15.py", "--seed", "2")
    a, b = out.split("T-RO Fig. 15(b)")[0], out.split("T-RO Fig. 15(b) = ICRA'08 Fig. 8(b): ")[1]
    assert "map 13 (140 vertices, 907 critical lines), path P13 (9 waypoints, 192 critical points)" in a
    assert row(a, "component events") == ["85", "85", "85"]
    assert row(a, "  of which Java prints bounds for") == ["18", "18", "18"]
    assert row(a, "BG vertices, Java (pooled t0)") == ["41", "12+29=41", "12+29=41"]
    assert row(a, "BG edges, Java (pooled t0)") == ["60", "80", "80"]
    assert "map 14 (267 vertices, 1877 critical lines), path P14b (21 waypoints, 671 critical points)" in b
    assert row(b, "component events") == ["385", "385", "385"]
    assert row(b, "shadow labels (IDs handed out)") == ["491", "491", "491"]
    assert row(b, "  of which Java prints bounds for") == ["12", "12", "12"]
    assert row(b, "BG vertices, one per t0 shadow") == ["124", "61+63=124", "61+63=124"]
    assert row(b, "BG edges, Java (pooled t0)")[:2] == ["339", "339"]
    assert out.count("truth inside the exact bounds: yes") == 4
    assert b.count("original oracle (merge bug): infeasible (B7)") == 2
    assert "Fig. 15(a)'s 41 vertices" in out


def test_original_maps(tmp_path):
    pytest.importorskip("matplotlib")
    out = run("original_maps.py", "--out", str(tmp_path))
    assert re.search(r"^ 14   267     134   1877   727   78  190  882  P14a, P14b +2  T-RO Fig. 15\(b\)", out, re.M)
    assert re.search(r"^ 12    34      15    117 .* P12 +3  T-RO Fig. 11", out, re.M)
    assert re.search(r"^P13\+ +13  code +12  commented: ProjectPanel +throws", out, re.M)
    assert "demo: closed tours on maps 2, 3, 4, 6, 7, 8, 9, 10, 11" in out
    assert (tmp_path / "original_maps.png").stat().st_size > 50_000
