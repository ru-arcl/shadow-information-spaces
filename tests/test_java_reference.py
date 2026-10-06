"""Cross-checks of the golden dumps of the original Java implementation (tests/fixtures/java).

The dumps are produced by tools/java_reference (see its README). These tests pin down
what the dumps say about the Java behaviour, and check that the dumps hold everything
needed to reproduce it in Python. The polygon front end of the port should add its own
parity tests against the same files.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools" / "java_reference"))

import convert  # noqa: E402

FIX = ROOT / "tests" / "fixtures" / "java"
RUNS = sorted(p.stem for p in FIX.glob("*_*.json"))
OK_RUNS = [r for r in RUNS if json.loads((FIX / f"{r}.json").read_text())["get_gaps"]["ok"]]


def load(name):
    return json.loads((FIX / f"{name}.json").read_text())


def test_index_and_random_clone():
    idx = load("index")
    rc = idx["java_random_check"]
    assert rc["equals_new_Random_seed_nextDouble"]
    assert convert.java_random_doubles(rc["seed"], len(rc["math_random_after_seeding"])) == rc["math_random_after_seeding"]
    assert [int(x * 100) for x in convert.java_random_doubles(42, 6)] == rc["new_Random_42_int_nextDouble_x100"]
    assert idx["non_round_trip_doubles"] == 0
    assert {r["file"][:-5] for r in idx["runs"]} == set(RUNS)


# cut counts (docs/notes/original_java.md §2.5): n, reflex, single-tangent, non-general, general, bitangent rays
CUT_COUNTS = {
    1: (31, 13, 41, 14, 12, 10), 2: (68, 32, 104, 34, 30, 80), 3: (18, 9, 17, 10, 8, 20),
    4: (124, 58, 181, 70, 46, 116), 5: (253, 127, 403, 124, 130, 394), 6: (18, 10, 16, 12, 8, 12),
    7: (22, 8, 31, 6, 10, 18), 8: (30, 12, 35, 12, 12, 24), 9: (33, 13, 44, 12, 14, 20),
    10: (54, 25, 87, 30, 20, 56), 11: (45, 21, 58, 24, 18, 56), 12: (34, 15, 47, 18, 12, 40),
    13: (140, 68, 407, 74, 62, 364), 14: (267, 134, 727, 78, 190, 882),
}


@pytest.mark.parametrize("n", range(1, 15))
def test_polygon_dump(n):
    d = load(f"poly{n}")
    c, cuts = d["counts"], d["cuts"]
    got = (d["n"], c["reflex"], c["single_tangent"], c["inflection_nongeneral"], c["general_inflection_cut"],
           c["bitangent_cut"])
    assert got == CUT_COUNTS[n]
    assert len(d["vertices"]) == d["n"] and d["reflex"] == [i for i, t in enumerate(d["turn"]) if t == 1]
    assert cuts["inflection_general_equals_general_inflection_cut_lines"]
    assert cuts["bitangent_lines_equal_bitangent_cut_lines"]
    assert [k for k, _ in cuts["get_cuts_layout"]] == ["SINGLETANGENT", "NONGENERAL_INFLECTION",
                                                      "GENERAL_INFLECTION", "BITANGENT"]
    # y-up, counter-clockwise data (shoelace area > 0); the reflex test depends on it
    v = d["vertices"]
    area = sum(v[i][0] * v[(i + 1) % len(v)][1] - v[(i + 1) % len(v)][0] * v[i][1] for i in range(len(v)))
    assert area > 0
    assert "errors" not in d


@pytest.mark.parametrize("name", OK_RUNS)
def test_gap_history_events_bg_nse(name):
    run = load(name)
    rep = convert.check_run(run)
    assert rep["sets_equal_samples"] and rep["getGaps_repeatable"]
    assert rep["gap_history_anomalies"] == []
    assert rep["gap_history_sequence_valid"] is True
    assert rep["nse"]["equal"], rep["nse"]
    assert rep["bg_structure_seed_independent"]
    for o in rep["oracle_runs"]:
        assert o["problems"] == [] and o["events_equal_gap_history"] and o["times_equal_set_times"]
        bg = o["bg"]
        assert bg["edges_equal"] and bg["left_weights_equal"] and bg["right_weights_equal"], bg
        b = o["bounds"]
        if o["variant"] == "merge_fixed":
            # consistent observations: the oracle's hidden counts satisfy the exact bounds
            assert b["feasible"] and b["truth_in_exact"]
        if b["feasible"]:
            # observed on every run: some identity-hash order makes Java print the exact bounds
            assert b["any_outcome_equals_exact"]


@pytest.mark.parametrize("name", OK_RUNS)
def test_oracle_replica_and_text_formats(name):
    run = load(name)
    for k, o in enumerate(run["oracle_runs"]):
        rows = convert.replay_oracle(run, o["seed"], o["n_targets"], merge_bug=o["variant"] == "original")
        assert rows == [e[:6] for e in o["events"]], (o["variant"], o["seed"])
        assert convert.format_equations(o) == o["equations"]
        if k == 0:
            assert [convert.format_event(e) for e in o["events"]] == [e[6] for e in o["events"]]
    for s in run["get_gaps"]["sets"]:
        body = "".join(f"[{i}, {ie}" + (" | " + "".join(f"{x} " for x in to) if to else "") + "] "
                       for i, ie, to in s["gaps"])
        assert s["line"] == f"{s['t']!r} {body}"


def test_java_bounds_are_order_dependent_and_b7_breaks_consistency():
    rep = convert.check_run(load("12_PP"))
    assert all(o["bounds"]["n_java_outcomes"] > 1 for o in rep["oracle_runs"])
    # B11: s18 appeared and survived, Java prints no bound for it
    assert all(o["bounds"]["shadows_without_java_bound"] == [18] for o in rep["oracle_runs"])
    pp5 = convert.check_run(load("14_PP5"))
    assert [o["bounds"]["feasible"] for o in pp5["oracle_runs"] if o["variant"] == "original"] == [False] * 5


EXPECTED_FAILURES = {
    "14_PP4": ("java.lang.NullPointerException", 654, "BITANGENT", 286),
    "13_c_second": ("java.lang.ArrayIndexOutOfBoundsException", 239, "GENERAL_INFLECTION", 138),
    "13_c_tail": ("java.lang.ArrayIndexOutOfBoundsException", 239, "GENERAL_INFLECTION", 138),
    "14_fig_ICRA08-Fig8b": ("java.lang.ArrayIndexOutOfBoundsException", 596, "GENERAL_INFLECTION", 138),
}


def test_get_gaps_failures():
    failing = sorted(set(RUNS) - set(OK_RUNS))
    assert failing == sorted(EXPECTED_FAILURES)
    for name, (cls, idx, typ, line) in EXPECTED_FAILURES.items():
        g = load(name)["get_gaps"]
        assert g["exception"]["class"] == cls
        assert (g["failure"]["critical_index"], g["failure"]["critical_type"], g["failure"]["getGaps_line"]) == (idx, typ, line)


def test_icra08_fig4_is_the_polygon12_panel_run():
    chk = convert.icra08_fig4_check(load("12_PP"))
    assert chk["bipartite_edges_equal"] and chk["same_event_multiset"]


def test_paper_fig15_counts():
    a = convert.paper_counts(load("13_c_first"), load("13_c_first")["oracle_runs"][0])
    assert a["component_events"] == convert.PAPER["fig15a"]["events"] == 85
    assert a["final_shadows_bounded_by_java"] == convert.PAPER["fig15a"]["final_shadows"] == 18
    assert a["java_bg_vertices"] == convert.PAPER["fig15a"]["bg_vertices"] == 41
    assert a["java_bg_edges"] == 80 != convert.PAPER["fig15a"]["bg_edges"]  # paper: 60
    b = convert.paper_counts(load("14_PP5"), load("14_PP5")["oracle_runs"][0])
    assert b["component_events"] == 385 and b["shadow_ids"] == 491
    assert b["final_shadows_bounded_by_java"] == 12
    assert b["java_bg_edges"] == 339
    assert b["java_bg_vertices"] == 114 and b["shadowinfo_bg_vertices"] == 124  # paper: 124
