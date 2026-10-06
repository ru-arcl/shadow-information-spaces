"""The JS port of the polygon world (``docs/js/polygon.js``, ``polygon_sim.js``, ``maps_polygon.js``): generated
files are up to date.

``npm test`` checks the JS itself: ``js-tests/polygon.test.js`` against the Java golden fixtures
(``tests/fixtures/java``) and ``js-tests/polygon_parity.test.js`` against ``tests/fixtures/polygon_parity.json``,
which ``tools/gen_polygon_parity_fixture.py`` writes from :mod:`shadowinfo.polygon`.  Here we check that the
embedded maps and the parity fixture still equal what the generators produce from the current Python code and
data, so that a change on the Python side cannot silently leave the JS tests comparing against stale numbers.
"""

from __future__ import annotations

import importlib.util
import json
import warnings
from pathlib import Path

from shadowinfo.polygon import demo_paths, load_paths, load_polygon

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "polygon_parity.json"
JS = ROOT / "docs" / "js"


def _load_tool_from(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _load_tool(name):
    return _load_tool_from(ROOT / "tools" / f"{name}.py")


def test_js_polygon_maps_up_to_date():
    expected = _load_tool("gen_js_polygon_maps").render()
    assert (JS / "maps_polygon.js").read_text() == expected, "run: python tools/gen_js_polygon_maps.py"


def test_js_polygon_maps_content():
    """The embedded data is the package data: vertices as loaded, every path and demo path."""
    text = (JS / "maps_polygon.js").read_text()
    for n in range(1, 15):
        verts = json.dumps([list(v) for v in load_polygon(n).vertices], separators=(", ", ": "))
        assert verts in text, n
    for label in load_paths():
        assert json.dumps(label) + ": {" in text
    assert len(demo_paths()) == 14


def test_parity_fixture_covers_every_map():
    fx = json.loads(FIXTURE.read_text())
    assert [h["map"] for h in fx["histories"]] == list(range(1, 15))
    seeds = {(r["map"], r["seed"]) for r in fx["runs"] if not r["script"]}
    assert seeds == {(n, s) for n in range(1, 15) for s in (1, 2)}
    assert sum(1 for r in fx["runs"] if r["script"]) >= 3
    kinds = {e["type"] for r in fx["runs"] for f in r["frames"] for e in f.get("events", ())}
    assert kinds == {"appear", "disappear", "split", "merge", "enter", "exit"}
    assert any(r["stats"]["legs"] >= 2 for r in fx["runs"] if not r["script"])  # loop / reverse turnarounds
    assert sum(len(r.get("geometry", {})) for r in fx["runs"]) >= 3 * 14


def test_parity_fixture_up_to_date():
    """The fixture equals what the generator writes now.

    If its ``inputs_sha256`` equals the fingerprint of the current inputs (``shadowinfo`` sources, polygon data,
    generator), it was generated from exactly this code and data.  Otherwise the fixture is regenerated (about
    15 s; tracks are cached per polygon and path, shared with ``test_polygon_simulate``) and compared; if only
    the fingerprint differs, the test passes with a warning to refresh it.
    """
    gen = _load_tool("gen_polygon_parity_fixture")
    stored = json.loads(FIXTURE.read_text())
    if stored.get("inputs_sha256") == gen.input_fingerprint():
        return
    fresh = json.loads(json.dumps(gen.build_fixture()))
    hint = "run: python tools/gen_polygon_parity_fixture.py"
    assert fresh.keys() == stored.keys(), hint
    assert fresh["_doc"] == stored["_doc"], hint
    assert len(fresh["histories"]) == len(stored["histories"]), hint
    for a, b in zip(fresh["histories"], stored["histories"]):
        assert a == b, (hint, a["map"])
    assert len(fresh["runs"]) == len(stored["runs"]), hint
    for a, b in zip(fresh["runs"], stored["runs"]):
        key = (a["map"], a["seed"], bool(a["script"]))
        assert {k: v for k, v in a.items() if k != "frames"} == {k: v for k, v in b.items() if k != "frames"}, (hint, key)
        assert len(a["frames"]) == len(b["frames"]), (hint, key)
        for fa, fb in zip(a["frames"], b["frames"]):
            assert fa == fb, (hint, key, fa["t"])
    warnings.warn("tests/fixtures/polygon_parity.json is up to date but its inputs_sha256 is stale; " + hint)


def test_parity_fingerprint_covers_inputs():
    """The fingerprint covers the shadowinfo sources and the polygon data, so a change triggers the full check."""
    gen = _load_tool("gen_polygon_parity_fixture")
    files = {f.relative_to(ROOT).as_posix() for f in gen.input_files()}
    for rel in ("shadowinfo/polygon/gaps.py", "shadowinfo/polygon/simulate.py", "shadowinfo/rng.py",
                "shadowinfo/events.py", "shadowinfo/grid/simulate.py", "shadowinfo/polygon/data/14.dat",
                "shadowinfo/polygon/data/paths.json", "tools/gen_polygon_parity_fixture.py"):
        assert rel in files, rel
    assert not any("__pycache__" in f for f in files)
    fp = gen.input_fingerprint()
    assert len(fp) == 64 and fp == gen.input_fingerprint()
