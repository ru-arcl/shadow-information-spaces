"""Package namespace: only the public API is exported; the README snippet prints what it claims."""

from __future__ import annotations

import contextlib
import io
import re
from pathlib import Path

import shadowinfo
from shadowinfo import events

README = Path(__file__).resolve().parents[1] / "README.md"


def test_no_stdlib_or_typing_names_leak():
    for name in ("json", "math", "dataclass", "field", "annotations", "Dict", "List", "Optional", "Tuple",
                 "Union", "Iterable"):
        assert not hasattr(shadowinfo, name), name
    assert set(events.__all__) <= set(shadowinfo.__all__)
    assert all(hasattr(shadowinfo, n) for n in shadowinfo.__all__)
    ns: dict = {}
    exec("from shadowinfo import *", ns)
    assert "json" not in ns and "CombinatorialFilter" in ns and "RedOrBlueFilter" in ns


def test_readme_snippet_output():
    text = README.read_text()
    code = re.search(r"```python\n(.*?)```", text, re.S).group(1)
    expected = re.search(r"```text\n(.*?)```", text, re.S).group(1)
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        exec(code, {})
    assert out.getvalue() == expected
