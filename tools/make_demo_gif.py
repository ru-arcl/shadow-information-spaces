#!/usr/bin/env python3
"""Regenerate docs/assets/demo.gif, the animated hero image of the README.

The GIF is a sequence of snapshots of the interactive demo (docs/index.html)
on its default setting: map 12 (T-RO Fig. 11(a)) with its default path,
default seed and ground truth on. Each frame is the page fast-forwarded to a
later tick with the ``?ticks=N`` URL parameter, so the animation shows shadows
appearing, splitting, merging and disappearing while the count bounds in the
Shadows table tighten.

Requirements: ``google-chrome`` (or ``chromium``; override with $CHROME) on
PATH, and Pillow. Run from anywhere::

    python tools/make_demo_gif.py            # writes docs/assets/demo.gif
    python tools/make_demo_gif.py --keep-png /tmp/frames   # also keep the PNGs

Each frame is rendered headless at 1400x900 CSS px with device scale factor 2,
cropped to the map card and the Shadows table, then downscaled for crisp text.
The Shadows table's 300 px scroll limit is lifted in the capture so no row is cut.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
PAGE = ROOT / "docs" / "index.html"
OUT = ROOT / "docs" / "assets" / "demo.gif"

# Ticks to snapshot along the first pass of map 12's default path (it ends near t = 235).
TICKS = [0, 14, 21, 36, 49, 59, 68, 85, 96, 110, 131, 150, 170, 200, 232]
VIEW_W, VIEW_H, SCALE = 1400, 900, 2
# Crop box in CSS px: map card (left) + Shadows card (right), below the intro text.
CROP = (24, 214, 1376, 864)
OUT_W = 1100
FRAME_MS, FIRST_MS, LAST_MS = 800, 1400, 2600
COLOURS = 256
EXACT_SLOTS = 48  # most frequent exact colours kept verbatim


def find_chrome() -> str:
    for name in (os.environ.get("CHROME"), "google-chrome", "google-chrome-stable", "chromium", "chromium-browser"):
        if name and shutil.which(name):
            return shutil.which(name)
    sys.exit("make_demo_gif: google-chrome (or chromium) not found; set $CHROME")


# Let the Shadows table show all rows (it scrolls past 300 px on the page) so no row is cut.
UNCLIP = ("var d=this.contentDocument,s=d.createElement('style');"
          "s.textContent='.table-wrap{max-height:none!important}';d.head.appendChild(s)")


def snapshot(chrome: str, ticks: int, out: Path, tmp: Path) -> None:
    # Headless Chrome's --window-size does not give an exact viewport (part of the window
    # height is lost), so the page is loaded into a VIEW_W x VIEW_H iframe of a wrapper page
    # rendered in a taller window.
    url = f"{PAGE.as_uri()}?map=12&ticks={ticks}&theme=light"
    wrap = tmp / "wrap.html"
    wrap.write_text(
        "<!doctype html><html><head><style>html,body{margin:0;background:#fff}"
        f"iframe{{border:0;display:block;width:{VIEW_W}px;height:{VIEW_H}px}}</style></head>"
        f'<body><iframe src="{url}" onload="{UNCLIP}"></iframe></body></html>')
    cmd = [chrome, "--headless=new", "--no-sandbox", "--disable-gpu", "--hide-scrollbars",
           "--allow-file-access-from-files", f"--user-data-dir={tmp / 'profile'}",
           f"--force-device-scale-factor={SCALE}", f"--window-size={VIEW_W},{VIEW_H + 200}",
           "--virtual-time-budget=3000", f"--screenshot={out}", wrap.as_uri()]
    subprocess.run(cmd, check=True, timeout=120, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if not out.exists():
        sys.exit(f"make_demo_gif: chrome produced no screenshot for ticks={ticks}")


def build_frame(png: Path) -> Image.Image:
    im = Image.open(png).convert("RGB")
    s = im.width / VIEW_W  # actual device scale
    box = tuple(round(v * s) for v in CROP)
    im = im.crop(box)
    h = round(im.height * OUT_W / im.width)
    return im.resize((OUT_W, h), Image.LANCZOS)


def shared_palette(frames: list[Image.Image]) -> Image.Image:
    """One adaptive palette for all frames. The most frequent exact colours (flat UI fills,
    shadow regions) are kept verbatim; the other slots come from a median cut in which
    saturated pixels (legend swatches, target dots, coloured labels) are over-weighted, so
    small coloured marks do not collapse onto the greys of anti-aliased text."""
    peak: dict[tuple, int] = {}  # colour -> largest pixel count in any single frame
    for f in frames:
        for n, c in f.getcolors(maxcolors=1 << 24):
            peak[c] = max(peak.get(c, 0), n)
    exact = sorted(peak, key=lambda c: -peak[c])[:EXACT_SLOTS]
    sat = [c for c in peak if max(c) - min(c) >= 80 and c not in exact]
    w, h = frames[0].size
    train = [px for f in frames for px in f.getdata()]
    train += sat * max(1, len(train) // (3 * max(1, len(sat))))  # saturated colours ~ 1/4 of the sample
    side = int(len(train) ** 0.5) + 1
    img = Image.new("RGB", (side, side))
    img.putdata(train + [train[0]] * (side * side - len(train)))
    rest = img.quantize(colors=COLOURS - len(exact), method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
    flat = [v for c in exact for v in c] + rest.getpalette()[: 3 * (COLOURS - len(exact))]
    pal = Image.new("P", (1, 1))
    pal.putpalette(flat + [0] * (768 - len(flat)))
    return pal


def assemble(frames: list[Image.Image], out: Path) -> None:
    # Shared palette, no dithering: flat UI colours stay flat and frames diff well.
    pal = shared_palette(frames)
    q = [f.quantize(palette=pal, dither=Image.Dither.NONE) for f in frames]
    durations = [FIRST_MS] + [FRAME_MS] * (len(q) - 2) + [LAST_MS]
    out.parent.mkdir(parents=True, exist_ok=True)
    q[0].save(out, save_all=True, append_images=q[1:], duration=durations, loop=0,
              optimize=True, disposal=1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT, help="output GIF (default: docs/assets/demo.gif)")
    ap.add_argument("--keep-png", type=Path, help="directory to keep the raw screenshots in")
    args = ap.parse_args()

    chrome = find_chrome()
    with tempfile.TemporaryDirectory(prefix="demo-gif-") as tmp:
        tmp = Path(tmp)
        shots_dir = args.keep_png or tmp
        shots_dir.mkdir(parents=True, exist_ok=True)
        frames = []
        for t in TICKS:
            png = shots_dir / f"t{t:04d}.png"
            png.unlink(missing_ok=True)
            snapshot(chrome, t, png, tmp)
            frames.append(build_frame(png))
            print(f"  t = {t:4d}  captured", file=sys.stderr)
        assemble(frames, args.out)
    size = args.out.stat().st_size
    print(f"wrote {args.out} ({len(frames)} frames, {frames[0].size[0]}x{frames[0].size[1]}, {size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
