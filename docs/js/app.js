/* Demo page (DESIGN §6): canvas rendering, side panel and controls around SI.Controller.
 * Two worlds: the grid office map and the 14 polygon maps of the original Java code (y-up .dat frame,
 * drawn y-up like the Java window, never flipped).
 * URL parameters: map (office, 1..14), path (default, a stored label, or x1,y1;x2,y2;... drawn waypoints),
 * autoplay=1, ticks=N (fast-forward), n, seed, mode, radius, gt=0, speed, theme. */
(function (SI) {
  "use strict";

  const $ = (id) => document.getElementById(id);
  const INF = Infinity;
  const LOG_TICKS = 120;
  const DEFAULT_MAP = "12";
  const params = new URLSearchParams(location.search);

  const els = {
    wrap: $("canvas-wrap"), canvas: $("map"), play: $("btn-play"), step: $("btn-step"), reset: $("btn-reset"),
    speed: $("speed"), speedOut: $("speed-out"), tick: $("tick"), n: $("opt-n"), seed: $("opt-seed"),
    mode: $("opt-mode"), radius: $("opt-radius"), gt: $("opt-gt"), rows: $("shadow-rows"),
    aggregate: $("aggregate"), notes: $("notes"), graph: $("graph"), log: $("log"), warning: $("warning"),
    map: $("opt-map"), path: $("opt-path"), draw: $("btn-draw"), drawBar: $("draw-bar"), drawMsg: $("draw-msg"),
    drawUndo: $("draw-undo"), drawCancel: $("draw-cancel"), drawStart: $("draw-start"), pathInfo: $("path-info"),
    scene: $("scene"),
  };
  const ctx = els.canvas.getContext("2d");

  const state = {
    ctrl: null,
    playing: false,
    gt: true,
    speed: 6,
    stepAt: 0,
    animDur: 1,
    prev: null,
    dirtyCells: true,
    needsDraw: true,
    cell: 10,
    dpr: 1,
    theme: {},
    anchors: [],
    graphW: 0,
    view: null,
    custom: new Map(),
    draw: null,
  };

  const cells = document.createElement("canvas");
  const cellsCtx = cells.getContext("2d");

  /* ---- options ------------------------------------------------------------ */

  function clampInt(v, lo, hi, dflt) {
    const x = parseInt(v, 10);
    return Number.isFinite(x) ? Math.min(hi, Math.max(lo, x)) : dflt;
  }

  const isPolygon = (map) => SI.demo.isPolygonMap(map);

  function buildMapSelect() {
    const groups = [["Grid world", "grid"], ["Original maps (Java code, polygons)", "polygon"]];
    let html = "";
    groups.forEach(([label, world]) => {
      html += '<optgroup label="' + esc(label) + '">';
      SI.demo.mapList().filter((m) => m.world === world).forEach((m) => {
        html += '<option value="' + m.id + '" title="' + esc(m.title) + '">' + esc(m.name) + "</option>";
      });
      html += "</optgroup>";
    });
    els.map.innerHTML = html;
    els.map.value = DEFAULT_MAP;
  }

  /* Path choices of the selected map: the demo path, every stored path (invalid ones disabled),
   * and the path drawn on this map, if any. */
  function buildPathSelect(selected) {
    const map = els.map.value;
    let html = "";
    if (!isPolygon(map)) {
      html = '<option value="default">Patrol loop</option>';
    } else {
      SI.demo.pathOptions(map).forEach((o) => {
        html += '<option value="' + esc(o.id) + '"' + (o.valid ? "" : " disabled") + ' title="' +
          esc(o.valid ? o.info : o.error) + '">' + esc(o.name + (o.valid ? "" : " — leaves the map")) + "</option>";
      });
      const c = state.custom.get(map);
      if (c) html += '<option value="custom">Drawn path (' + c.length + " waypoints)</option>";
    }
    els.path.innerHTML = html;
    const ok = Array.from(els.path.options).some((o) => o.value === selected && !o.disabled);
    els.path.value = ok ? selected : "default";
    els.path.disabled = !isPolygon(map);
    els.draw.hidden = !isPolygon(map);
    if (isPolygon(map)) {
      if (els.radius.value) els.radius.dataset.prev = els.radius.value;
      els.radius.value = "";
      els.radius.disabled = true;
      els.radius.title = "Polygon maps: unlimited range, as in the paper";
    } else if (els.radius.disabled) {
      els.radius.disabled = false;
      els.radius.title = "";
      if (els.radius.dataset.prev) els.radius.value = els.radius.dataset.prev;
      delete els.radius.dataset.prev;
    }
  }

  function parseWaypoints(text) {
    const pts = String(text).split(";").map((p) => p.split(",").map(Number));
    return pts.length >= 2 && pts.every((p) => p.length === 2 && p.every(Number.isFinite)) ? pts : null;
  }

  function applyParams() {
    if (params.has("n")) els.n.value = clampInt(params.get("n"), 1, 40, 10);
    if (params.has("seed")) els.seed.value = clampInt(params.get("seed"), 1, 99999, 1);
    if (params.has("mode") && SI.grid.MODES.includes(params.get("mode"))) els.mode.value = params.get("mode");
    if (params.has("radius")) els.radius.value = params.get("radius");
    if (params.get("gt") === "0") els.gt.checked = false;
    if (params.has("speed")) els.speed.value = clampInt(params.get("speed"), 1, 60, 6);
    const theme = params.get("theme");
    if (theme === "dark" || theme === "light") document.documentElement.dataset.theme = theme;
    const map = params.get("map");
    if (map && SI.demo.mapList().some((m) => m.id === map)) els.map.value = map;
    let path = params.get("path") || "default";
    if (path.includes(";") && isPolygon(els.map.value)) {
      const pts = parseWaypoints(path);
      let ok = !!pts;
      if (ok) {
        try {
          SI.polygonSim.validatePath(SI.polygon.loadPolygon(+els.map.value), pts);
        } catch (err) {
          ok = false;
        }
      }
      if (ok) state.custom.set(els.map.value, pts);
      path = ok ? "custom" : "default";
    }
    buildPathSelect(path);
  }

  function options() {
    const mode = els.mode.value, evader = mode === "evader";
    if (evader && !els.n.disabled) {
      els.n.dataset.prev = els.n.value;
    } else if (!evader && els.n.disabled && els.n.dataset.prev) {
      els.n.value = els.n.dataset.prev;
      delete els.n.dataset.prev;
    }
    els.n.disabled = evader;
    if (evader) els.n.value = 1;
    const n = clampInt(els.n.value, 1, 40, 10);
    const seed = clampInt(els.seed.value, 1, 99999, 1);
    els.n.value = n;
    els.seed.value = seed;
    const map = els.map.value;
    const opts = { map: map, nTargets: n, seed: seed, mode: mode, radius: els.radius.value ? +els.radius.value : null };
    if (isPolygon(map)) {
      const p = els.path.value;
      opts.path = p === "custom" ? state.custom.get(map) : p;
    }
    return opts;
  }

  /* ---- theme and colours -------------------------------------------------- */

  function readTheme() {
    const cs = getComputedStyle(document.documentElement);
    const v = (k) => cs.getPropertyValue(k).trim();
    state.theme = {
      wall: v("--wall"), floor: v("--floor"), visible: v("--visible"), path: v("--path"), robot: v("--robot"),
      ring: v("--robot-ring"), target: v("--target"), hidden: v("--hidden-target"), range: v("--range"),
      pillBg: v("--pill-bg"), pillText: v("--pill-text"), l: v("--shadow-l") || "58%",
      a: parseFloat(v("--shadow-a")) || 0.55, outside: v("--outside"), edge: v("--wall-edge"),
      draw: v("--draw"), bad: v("--danger"),
    };
  }

  /* Stable per-label colours: a split's larger child and a merge's larger input keep
   * the parent's colour, new labels take the least used palette slot; then resolveColours
   * separates live shadows that ended up with the same colour. */
  const PALETTE = [205, 140, 340, 265, 180, 48, 305, 95, 2, 230, 160, 285];
  const colours = { of: new Map(), next: 0 };

  function freshColour(label, live) {
    const used = new Array(PALETTE.length).fill(0);
    live.forEach((s) => { if (colours.of.has(s)) used[colours.of.get(s)] += 1; });
    let best = -1;
    for (let k = 0; k < PALETTE.length; k++) {
      const idx = (colours.next + k) % PALETTE.length;
      if (best < 0 || used[idx] < used[best]) best = idx;
    }
    colours.next = (best + 1) % PALETTE.length;
    colours.of.set(label, best);
  }

  function assignColours(events, prevAreas) {
    const ctrl = state.ctrl;
    const live = new Set(ctrl.shadows.map((s) => s.label));
    const ps = SI.demo.splitProbabilities(events, ctrl.areas());
    const area = new Map(prevAreas);
    const inherit = (child, parent) => {
      if (colours.of.has(parent)) colours.of.set(child, colours.of.get(parent));
      else freshColour(child, live);
    };
    events.forEach((e, i) => {
      if (e.type === "appear") {
        freshColour(e.s, live);
      } else if (e.type === "merge") {
        const a = area.get(e.a) || 0, b = area.get(e.b) || 0;
        inherit(e.s, a >= b ? e.a : e.b);
        area.set(e.s, a + b);
      } else if (e.type === "split") {
        const big = ps.get(i) >= 0.5 ? e.a : e.b, small = big === e.a ? e.b : e.a;
        inherit(big, e.s);
        freshColour(small, new Set(Array.from(live).concat([big])));
      }
    });
    if (events.length) resolveColours(live, ctrl.areas());
    if (colours.of.size > 4000) {
      Array.from(colours.of.keys()).forEach((s) => { if (!live.has(s)) colours.of.delete(s); });
    }
  }

  /* Inheritance can hand a live shadow's colour to another one, and with more shadows than palette
   * slots some must share.  Going from the largest shadow down, each of the PALETTE.length largest
   * shadows that has the colour of a larger one moves to the least used slot that no larger shadow
   * uses; smaller shadows keep theirs.  So the largest shadows always have distinct colours (all of
   * them while a slot is free), and a colour changes only on a tick with events. */
  function resolveColours(live, areaOf) {
    const L = PALETTE.length;
    const order = Array.from(live).sort((a, b) => (areaOf.get(b) || 0) - (areaOf.get(a) || 0) || a - b);
    const cnt = new Array(L).fill(0), seen = new Array(L).fill(0);
    order.forEach((s) => { if (!colours.of.has(s)) freshColour(s, live); cnt[colours.of.get(s)] += 1; });
    order.forEach((s, rank) => {
      const c = colours.of.get(s);
      if (rank < L && seen[c] > 0) {
        let best = -1;
        for (let k = 0; k < L; k++) {
          const d = (colours.next + k) % L;
          if (d === c || seen[d] > 0) continue;
          if (best < 0 || cnt[d] < cnt[best]) best = d;
        }
        if (best >= 0) {
          cnt[c] -= 1;
          cnt[best] += 1;
          colours.of.set(s, best);
          colours.next = (best + 1) % L;
        }
      }
      seen[colours.of.get(s)] += 1;
    });
  }

  function hue(label) {
    return colours.of.has(label) ? PALETTE[colours.of.get(label)] : (label * 137.508) % 360;
  }

  function shadowColor(label, alpha) {
    return "hsla(" + hue(label) + ", 68%, " + state.theme.l + ", " +
      (alpha === undefined ? state.theme.a : alpha) + ")";
  }

  function solidColor(label) {
    return "hsl(" + hue(label) + ", 62%, 48%)";
  }

  /* ---- formatting ----------------------------------------------------------- */

  function fmtBounds(lo, hi) {
    return "[" + lo + ", " + (hi === INF ? "∞)" : hi + "]");
  }

  function fmtShort(lo, hi) {
    if (lo === hi) return String(lo);
    if (hi === INF) return "≥" + lo;
    return lo + "–" + hi;
  }

  function esc(s) {
    return String(s).replace(/[&<>"]/g, (ch) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[ch]);
  }

  /* ---- simulation control ---------------------------------------------------- */

  function build() {
    stopDrawing();
    let ctrl;
    try {
      ctrl = new SI.Controller(options());
    } catch (err) {
      showWarning("Could not start: " + err.message);
      if (els.path.value === "default") throw err;
      els.path.value = "default";
      ctrl = new SI.Controller(options());
    }
    state.ctrl = ctrl;
    colours.of.clear();
    colours.next = 0;
    ctrl.shadows.forEach((sh) => freshColour(sh.label, new Set()));
    resolveColours(new Set(ctrl.shadows.map((sh) => sh.label)), ctrl.areas());
    state.prev = null;
    state.dirtyCells = true;
    els.log.textContent = "";
    logTicks.length = 0;
    updatePathInfo();
    sizeCanvas();
    updatePanel();
    schedulePrefetch(ctrl);
  }

  /* The gap track of the reversed path is computed while the page is idle, not at the turnaround. */
  function schedulePrefetch(ctrl) {
    const run = () => { if (state.ctrl === ctrl) ctrl.prefetch(); };
    if (window.requestIdleCallback) window.requestIdleCallback(run, { timeout: 1500 });
    else setTimeout(run, 200);
  }

  function updatePathInfo() {
    const ctrl = state.ctrl;
    let text;
    if (ctrl.world === "grid") {
      text = "Grid world of this repository: the robot patrols a fixed loop through the office, one cell per tick.";
    } else if (els.path.value === "custom") {
      const closed = ctrl.sim.atEnd === "loop";
      text = "Your drawn path (" + (closed ? ctrl.path.length - 1 + " waypoints, closed" : ctrl.path.length + " waypoints") +
        "), driven " + (closed ? "in a loop." : "back and forth.");
    } else {
      const o = SI.demo.pathOptions(ctrl.opts.map).find((p) => p.id === els.path.value);
      text = o ? o.info : "";
    }
    const m = SI.demo.mapList().find((x) => x.id === ctrl.opts.map);
    if (ctrl.world === "polygon") {
      const rec = SI.polygonMaps.polygons[ctrl.opts.map];
      text = rec.file + " (" + rec.n_vertices + " vertices" + (m.figure ? ", " + m.figure : "") + "). " + text;
    }
    els.pathInfo.textContent = text;
  }

  function snapshotPositions() {
    const sim = state.ctrl.sim;
    return { robot: sim.robot.slice(), targets: sim.targets.map((p) => (Array.isArray(p) ? p.slice() : p)) };
  }

  function doStep(now, record) {
    const before = snapshotPositions();
    const areas = state.ctrl.shadows.map((sh) => [sh.label, sh.area]);
    try {
      state.ctrl.step();
    } catch (err) {
      setPlaying(false);
      showWarning("Filter error: " + err.message);
      return false;
    }
    assignColours(state.ctrl.lastEvents, areas);
    state.prev = before;
    state.stepAt = now;
    state.animDur = state.playing ? 1000 / state.speed : Math.min(1000 / state.speed, 260);
    state.dirtyCells = true;
    state.needsDraw = true;
    if (record !== false) logEvents(state.ctrl.t, state.ctrl.lastEvents);
    return true;
  }

  function setPlaying(on) {
    if (on && state.draw) return;
    state.playing = on;
    els.play.classList.toggle("playing", on);
    els.play.querySelector(".lbl").textContent = on ? "Pause" : "Play";
    els.play.setAttribute("aria-pressed", on ? "true" : "false");
    if (on) state.stepAt = performance.now();
  }

  function stepOnce() {
    if (state.draw) return;
    setPlaying(false);
    if (doStep(performance.now())) updatePanel();
  }

  /* ---- canvas geometry -------------------------------------------------------- */

  /* World <-> canvas (CSS px).  Grid: cell (r, c) centre at ((c + .5) cell, (r + .5) cell).
   * Polygon: the .dat frame, y up, framed by the polygon's bounding box plus a margin. */
  function makeView(ctrl, w) {
    if (ctrl.world === "grid") {
      const map = ctrl.map;
      return { world: "grid", aspect: map.cols / map.rows, cell: w / map.cols };
    }
    const vx = ctrl.map.vertices.map((p) => p[0]), vy = ctrl.map.vertices.map((p) => p[1]);
    const x0 = Math.min(...vx), x1 = Math.max(...vx), y0 = Math.min(...vy), y1 = Math.max(...vy);
    const pad = 0.025 * Math.max(x1 - x0, y1 - y0);
    const bw = x1 - x0 + 2 * pad, bh = y1 - y0 + 2 * pad;
    const scale = w / bw;
    return { world: "polygon", aspect: bw / bh, x0: x0 - pad, y1: y1 + pad, scale: scale, cell: w / 38 };
  }

  const sx = (p) => (p[0] - state.view.x0) * state.view.scale;
  const sy = (p) => (state.view.y1 - p[1]) * state.view.scale;

  function toWorld(X, Y) {
    const v = state.view, r = (x) => Math.round(x * 10) / 10;
    return [r(X / v.scale + v.x0), r(v.y1 - Y / v.scale)];
  }

  function sizeCanvas() {
    const ctrl = state.ctrl;
    const probe = makeView(ctrl, 1);
    els.wrap.style.setProperty("--ar", probe.aspect.toFixed(4));
    const w = els.wrap.clientWidth;
    const h = Math.round(w / probe.aspect);
    const dpr = Math.max(1, window.devicePixelRatio || 1);
    state.dpr = dpr;
    state.view = makeView(ctrl, w);
    state.cell = state.view.cell;
    els.canvas.width = Math.round(w * dpr);
    els.canvas.height = Math.round(h * dpr);
    cells.width = els.canvas.width;
    cells.height = els.canvas.height;
    state.dirtyCells = true;
    state.needsDraw = true;
  }

  /* ---- grid world -------------------------------------------------------------- */

  function drawCells() {
    if (state.ctrl.world === "polygon") {
      drawPolygonLayer();
      return;
    }
    const ctrl = state.ctrl, map = ctrl.map, th = state.theme;
    const c = state.cell * state.dpr;
    const g = cellsCtx;
    g.setTransform(1, 0, 0, 1, 0, 0);
    g.fillStyle = th.wall;
    g.fillRect(0, 0, cells.width, cells.height);
    const labels = ctrl.sim.labels;
    const vis = map.visibility(ctrl.sim.robot, ctrl.opts.radius);
    const x = (col) => Math.round(col * c), y = (row) => Math.round(row * c);
    for (let i = 0; i < labels.length; i++) {
      if (map.blocked[i]) continue;
      const r = Math.floor(i / map.cols), col = i % map.cols;
      const lab = labels[i];
      if (lab) {
        g.fillStyle = th.floor;
        g.fillRect(x(col), y(r), x(col + 1) - x(col), y(r + 1) - y(r));
        g.fillStyle = shadowColor(lab);
      } else {
        g.fillStyle = vis[i] ? th.visible : th.floor;
      }
      g.fillRect(x(col), y(r), x(col + 1) - x(col), y(r + 1) - y(r));
    }
    g.lineWidth = Math.max(1, 1.2 * state.dpr);
    g.lineCap = "round";
    for (let i = 0; i < labels.length; i++) {
      const lab = labels[i];
      if (!lab) continue;
      const r = Math.floor(i / map.cols), col = i % map.cols;
      g.strokeStyle = shadowColor(lab, 0.95);
      g.beginPath();
      if (r === 0 || labels[i - map.cols] !== lab) { g.moveTo(x(col), y(r)); g.lineTo(x(col + 1), y(r)); }
      if (r === map.rows - 1 || labels[i + map.cols] !== lab) { g.moveTo(x(col), y(r + 1)); g.lineTo(x(col + 1), y(r + 1)); }
      if (col === 0 || labels[i - 1] !== lab) { g.moveTo(x(col), y(r)); g.lineTo(x(col), y(r + 1)); }
      if (col === map.cols - 1 || labels[i + 1] !== lab) { g.moveTo(x(col + 1), y(r)); g.lineTo(x(col + 1), y(r + 1)); }
      g.stroke();
    }
    state.anchors = computeAnchors(labels, map);
    state.dirtyCells = false;
  }

  /* Label position: the shadow cell nearest to the shadow's centroid; fitW/fitH bound the pill size. */
  function computeAnchors(labels, map) {
    const acc = new Map(), cell = state.cell;
    for (let i = 0; i < labels.length; i++) {
      const lab = labels[i];
      if (!lab) continue;
      const r = Math.floor(i / map.cols), c = i % map.cols;
      const a = acc.get(lab) || { n: 0, r: 0, c: 0, cells: [], r0: r, r1: r, c0: c, c1: c };
      a.n += 1;
      a.r0 = Math.min(a.r0, r);
      a.r1 = Math.max(a.r1, r);
      a.c0 = Math.min(a.c0, c);
      a.c1 = Math.max(a.c1, c);
      a.r += r;
      a.c += c;
      a.cells.push(i);
      acc.set(lab, a);
    }
    const out = [];
    acc.forEach((a, lab) => {
      const cr = a.r / a.n, cc = a.c / a.n;
      let best = a.cells[0], bd = INF;
      a.cells.forEach((i) => {
        const r = Math.floor(i / map.cols), c = i % map.cols;
        const d = (r - cr) * (r - cr) + (c - cc) * (c - cc);
        if (d < bd) { bd = d; best = i; }
      });
      out.push({ label: lab, x: ((best % map.cols) + 0.5) * cell, y: (Math.floor(best / map.cols) + 0.5) * cell,
                 size: a.n, fitW: (a.c1 - a.c0 + 2) * cell, fitH: (a.r1 - a.r0 + 2) * cell });
    });
    out.sort((p, q) => q.size - p.size);
    return out;
  }

  /* ---- polygon world ------------------------------------------------------------ */

  function tracePath(g, pts, close) {
    g.beginPath();
    pts.forEach((p, k) => (k ? g.lineTo(sx(p), sy(p)) : g.moveTo(sx(p), sy(p))));
    if (close) g.closePath();
  }

  /* Map, visible region and the shadow pockets of the current tick (cached between ticks). */
  function drawPolygonLayer() {
    const ctrl = state.ctrl, sim = ctrl.sim, th = state.theme, g = cellsCtx, dpr = state.dpr;
    const poly = ctrl.map.vertices;
    g.setTransform(1, 0, 0, 1, 0, 0);
    g.fillStyle = th.outside;
    g.fillRect(0, 0, cells.width, cells.height);
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.lineJoin = "round";
    tracePath(g, poly, true);
    g.fillStyle = th.floor;
    g.fill();
    tracePath(g, sim.visibilityPolygon(), true);
    g.fillStyle = th.visible;
    g.fill();
    const shadows = sim.shadows();
    shadows.forEach((sh) => {
      tracePath(g, sh.pocket, true);
      g.fillStyle = shadowColor(sh.label);
      g.fill();
      g.lineWidth = 1.2;
      g.strokeStyle = shadowColor(sh.label, 0.95);
      g.stroke();
    });
    tracePath(g, poly, true);
    g.lineWidth = 2;
    g.strokeStyle = th.edge;
    g.stroke();
    shadows.forEach((sh) => {
      g.beginPath();
      g.moveTo(sx(sh.window[0]), sy(sh.window[0]));
      g.lineTo(sx(sh.window[1]), sy(sh.window[1]));
      g.lineWidth = 1.6;
      g.setLineDash([5, 3]);
      g.strokeStyle = solidColor(sh.label);
      g.stroke();
      g.setLineDash([]);
    });
    state.anchors = shadows.map((sh) => pocketAnchor(sh)).filter(Boolean);
    state.anchors.sort((p, q) => q.size - p.size);
    state.dirtyCells = false;
  }

  /* Crossings of the horizontal line y (or the vertical line x when vertical) with a closed polyline. */
  function crossings(pts, c, vertical) {
    const out = [], n = pts.length, i0 = vertical ? 1 : 0, i1 = vertical ? 0 : 1;
    for (let i = 0; i < n; i++) {
      const a = pts[i], b = pts[(i + 1) % n];
      if ((a[i1] > c) !== (b[i1] > c)) out.push(a[i0] + ((c - a[i1]) * (b[i0] - a[i0])) / (b[i1] - a[i1]));
    }
    return out.sort((u, v) => u - v);
  }

  /* Pill position inside a pocket: over a few horizontal scan lines, the midpoint of the widest
   * inside run, scored by the free width and height around it (pockets need not be convex). */
  function pocketAnchor(sh) {
    const pts = sh.pocket.map((p) => [sx(p), sy(p)]);
    const ys = pts.map((p) => p[1]);
    const top = Math.min(...ys), bottom = Math.max(...ys);
    const area = Math.abs(SI.polygonSim.polygonArea(pts));
    let best = null;
    const K = 9;
    for (let k = 0; k < K; k++) {
      const y = top + ((k + 0.5) / K) * (bottom - top);
      const xs = crossings(pts, y, false);
      for (let i = 0; i + 1 < xs.length; i += 2) {
        const x = (xs[i] + xs[i + 1]) / 2, w = xs[i + 1] - xs[i];
        const yv = crossings(pts, x, true);
        let hgt = 0, yc = y;
        for (let j = 0; j + 1 < yv.length; j += 2) {
          if (yv[j] <= y && y <= yv[j + 1]) {
            hgt = yv[j + 1] - yv[j];
            yc = (yv[j] + yv[j + 1]) / 2;
          }
        }
        const score = Math.min(w, 3 * hgt);
        if (!best || score > best.score) best = { score: score, x: x, y: yc, fitW: w, fitH: hgt };
      }
    }
    if (!best) {
      const xs = pts.map((p) => p[0]);
      best = { x: (Math.min(...xs) + Math.max(...xs)) / 2, y: (top + bottom) / 2, fitW: 0, fitH: 0 };
    }
    return { label: sh.label, x: best.x, y: best.y, fitW: best.fitW, fitH: best.fitH, size: area };
  }

  /* ---- drawing -------------------------------------------------------------------- */

  function lerp(a, b, t) {
    return a + (b - a) * t;
  }

  function ease(t) {
    return t < 0.5 ? 2 * t * t : 1 - Math.pow(-2 * t + 2, 2) / 2;
  }

  function draw(now) {
    const ctrl = state.ctrl, th = state.theme;
    if (state.dirtyCells) drawCells();
    const dpr = state.dpr, cell = state.cell, grid = ctrl.world === "grid";
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.drawImage(cells, 0, 0);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const t = state.prev ? ease(Math.min(1, (now - state.stepAt) / state.animDur)) : 1;
    /* screen position of a robot/target position (grid: [row, col]; polygon: [x, y]) */
    const pos = grid ? (p) => [(p[1] + 0.5) * cell, (p[0] + 0.5) * cell] : (p) => [sx(p), sy(p)];

    const sim = ctrl.sim;
    ctx.save();
    ctx.strokeStyle = th.path;
    ctx.lineWidth = Math.max(1, grid ? cell * 0.12 : 1.6);
    ctx.setLineDash(grid ? [cell * 0.35, cell * 0.35] : [5, 5]);
    ctx.lineJoin = "round";
    ctx.beginPath();
    const path = grid ? ctrl.path : sim.path.points;
    path.forEach((p, k) => {
      const q = pos(p);
      if (k) ctx.lineTo(q[0], q[1]);
      else ctx.moveTo(q[0], q[1]);
    });
    if (grid) ctx.closePath();
    ctx.stroke();
    ctx.restore();

    const robotW = state.prev
      ? [lerp(state.prev.robot[0], sim.robot[0], t), lerp(state.prev.robot[1], sim.robot[1], t)]
      : sim.robot;
    const robot = pos(robotW);

    if (grid && ctrl.opts.radius) {
      ctx.save();
      ctx.strokeStyle = th.range;
      ctx.lineWidth = 1.5;
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      ctx.arc(robot[0], robot[1], (ctrl.opts.radius + 0.5) * cell, 0, 2 * Math.PI);
      ctx.stroke();
      ctx.restore();
    }

    drawPills(cell);

    const rad = Math.max(3, cell * 0.3);
    const cols = grid ? ctrl.map.cols : 0;
    const hiddenOf = grid ? null : sim.targetLabels;
    sim.targets.forEach((tg, j) => {
      const hidden = grid ? sim.labels[tg] !== 0 : hiddenOf[j] !== 0;
      if (hidden && !state.gt) return;
      let p = grid ? [Math.floor(tg / cols), tg % cols] : tg;
      if (state.prev) {
        const k = state.prev.targets[j];
        const q = grid ? [Math.floor(k / cols), k % cols] : k;
        p = [lerp(q[0], p[0], t), lerp(q[1], p[1], t)];
      }
      const s = pos(p);
      ctx.beginPath();
      ctx.arc(s[0], s[1], rad, 0, 2 * Math.PI);
      if (hidden) {
        ctx.lineWidth = Math.max(1.4, cell * 0.11);
        ctx.strokeStyle = th.hidden;
        ctx.stroke();
      } else {
        ctx.fillStyle = th.target;
        ctx.fill();
        ctx.lineWidth = 1;
        ctx.strokeStyle = th.ring;
        ctx.stroke();
      }
    });

    const rr = Math.max(4, cell * 0.55);
    ctx.beginPath();
    ctx.arc(robot[0], robot[1], rr + 2, 0, 2 * Math.PI);
    ctx.fillStyle = th.ring;
    ctx.fill();
    ctx.beginPath();
    ctx.arc(robot[0], robot[1], rr, 0, 2 * Math.PI);
    ctx.fillStyle = th.robot;
    ctx.fill();

    if (state.draw) drawDrawing(now);

    state.needsDraw = t < 1 || !!(state.draw && state.draw.flash);
  }

  /* Label pills, drawn under the targets and the robot.  A pill that would not fit inside its
   * shadow (anchor.fitW x anchor.fitH) shows only the label (the bounds are in the table). */
  function drawPills(cell) {
    const th = state.theme, w = els.canvas.width / state.dpr, h = els.canvas.height / state.dpr;
    const M = 4; /* margin: keeps pills clear of the wrapper's rounded corners */
    const byLabel = new Map(state.ctrl.shadows.map((s) => [s.label, s]));
    const fs = Math.max(9, Math.min(13, cell * 0.72));
    ctx.font = "600 " + fs.toFixed(1) + "px system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif";
    ctx.textBaseline = "middle";
    const placed = [];
    state.anchors.forEach((a) => {
      const sh = byLabel.get(a.label);
      if (!sh) return;
      const ph = fs + 6;
      const width = (t) => ctx.measureText(t).width + fs * 0.9 + 6;
      let text = "s" + a.label + " · " + fmtShort(sh.lo, sh.hi);
      let pw = width(text);
      if (pw > a.fitW || ph > a.fitH) {
        text = "s" + a.label;
        pw = width(text);
      }
      let x = Math.min(Math.max(M, a.x - pw / 2), w - pw - M);
      let y = Math.min(Math.max(M, a.y - ph / 2), h - ph - M);
      for (let k = 0; k < 6; k++) {
        const hit = placed.find((p) => x < p.x + p.w && p.x < x + pw && y < p.y + p.h && p.y < y + ph);
        if (!hit) break;
        y = k % 2 ? hit.y - ph - 2 : hit.y + hit.h + 2;
        y = Math.min(Math.max(M, y), h - ph - M);
      }
      placed.push({ x: x, y: y, w: pw, h: ph });
      ctx.beginPath();
      roundRect(ctx, x, y, pw, ph, ph / 2);
      ctx.fillStyle = th.pillBg;
      ctx.globalAlpha = 0.85;
      ctx.fill();
      ctx.globalAlpha = 1;
      ctx.lineWidth = 1.5;
      ctx.strokeStyle = solidColor(a.label);
      ctx.stroke();
      ctx.beginPath();
      ctx.arc(x + ph / 2, y + ph / 2, fs * 0.28, 0, 2 * Math.PI);
      ctx.fillStyle = solidColor(a.label);
      ctx.fill();
      ctx.fillStyle = th.pillText;
      ctx.fillText(text, x + ph / 2 + fs * 0.45, y + ph / 2 + 0.5);
    });
  }

  function roundRect(g, x, y, w, h, r) {
    g.moveTo(x + r, y);
    g.arcTo(x + w, y, x + w, y + h, r);
    g.arcTo(x + w, y + h, x, y + h, r);
    g.arcTo(x, y + h, x, y, r);
    g.arcTo(x, y, x + w, y, r);
    g.closePath();
  }

  /* ---- draw your own path ----------------------------------------------------------- */

  const DRAW_HELP = "Click inside the map to add waypoints; each leg must stay inside. " +
    "Click the start point again to close a loop. Double-click or press Enter to start, Esc to cancel.";

  function startDrawing() {
    if (state.ctrl.world !== "polygon" || state.draw) return;
    setPlaying(false);
    state.draw = { points: [], hover: null, flash: null };
    els.drawBar.style.minHeight = els.scene.offsetHeight + "px"; /* takes the map/path row's place without a layout shift */
    els.scene.hidden = true;
    els.drawBar.hidden = false;
    els.wrap.classList.add("drawing");
    els.draw.setAttribute("aria-pressed", "true");
    drawMessage(DRAW_HELP, false);
    state.needsDraw = true;
  }

  function stopDrawing() {
    if (!state.draw) return;
    state.draw = null;
    els.drawBar.hidden = true;
    els.scene.hidden = false;
    els.wrap.classList.remove("drawing");
    els.draw.setAttribute("aria-pressed", "false");
    state.needsDraw = true;
  }

  function drawMessage(text, bad) {
    els.drawMsg.textContent = text;
    els.drawMsg.classList.toggle("bad", !!bad);
    const n = state.draw ? state.draw.points.length : 0;
    els.drawStart.disabled = n < 2;
    els.drawUndo.disabled = n < 1;
  }

  function canvasPoint(ev) {
    const r = els.canvas.getBoundingClientRect();
    return [ev.clientX - r.left, ev.clientY - r.top];
  }

  const REASONS = {
    outside: "That point is outside the map. Waypoints must lie strictly inside the polygon.",
    leaves: "That leg leaves the polygon. Pick a point the robot can reach in a straight line.",
    along: "That leg runs along a wall. Move it off the wall a little.",
    near: "That is the last waypoint: click farther away to add one, or double-click / Enter / Start to run.",
  };

  const SNAP_PX = 9; /* a click this close to the start point (with 3+ waypoints) closes the loop */

  function addWaypoint(ev) {
    const d = state.draw;
    const s = canvasPoint(ev);
    let p = toWorld(s[0], s[1]);
    const last = d.points[d.points.length - 1] || null, first = d.points[0];
    const closing = d.points.length >= 3 && Math.hypot(sx(first) - s[0], sy(first) - s[1]) < SNAP_PX;
    if (closing) p = first.slice(); /* an exact copy: the path is then driven as a loop */
    else if (last && Math.hypot(sx(last) - s[0], sy(last) - s[1]) < 5) {
      if (ev.detail < 2) drawMessage(REASONS.near, false); /* a hint, as the click may start a double-click */
      return;
    }
    const why = SI.demo.checkWaypoint(state.ctrl.map, last, p);
    if (why === "repeat") return;
    if (why) {
      d.flash = { a: last, b: p, until: performance.now() + 1600 };
      drawMessage(REASONS[why], true);
    } else {
      d.points.push(p);
      d.flash = null;
      if (closing) {
        finishDrawing();
        return;
      }
      drawMessage(d.points.length < 2 ? "Start point set. Click to add the next waypoint." :
        d.points.length + " waypoints. Keep clicking, click the start to close a loop, or double-click / Enter / Start to run.", false);
    }
    state.needsDraw = true;
  }

  function finishDrawing() {
    const d = state.draw;
    if (!d || d.points.length < 2) return;
    const map = state.ctrl.opts.map;
    state.custom.set(map, d.points.slice());
    stopDrawing();
    buildPathSelect("custom");
    showWarning("");
    build();
  }

  function drawDrawing(now) {
    const d = state.draw, th = state.theme;
    ctx.save();
    ctx.fillStyle = "rgba(0, 0, 0, 0.18)";
    ctx.fillRect(0, 0, els.canvas.width / state.dpr, els.canvas.height / state.dpr);
    ctx.lineJoin = "round";
    ctx.lineCap = "round";
    const last = d.points[d.points.length - 1];
    if (last && d.hover && !d.flash) {
      const ok = !SI.demo.checkWaypoint(state.ctrl.map, last, d.hover);
      ctx.beginPath();
      ctx.moveTo(sx(last), sy(last));
      ctx.lineTo(sx(d.hover), sy(d.hover));
      ctx.setLineDash([6, 4]);
      ctx.lineWidth = 2;
      ctx.strokeStyle = ok ? th.draw : th.bad;
      ctx.stroke();
      ctx.setLineDash([]);
    }
    if (d.points.length) {
      tracePath(ctx, d.points, false);
      ctx.lineWidth = 3;
      ctx.strokeStyle = th.draw;
      ctx.stroke();
      d.points.forEach((p, k) => {
        ctx.beginPath();
        ctx.arc(sx(p), sy(p), k ? 4 : 6, 0, 2 * Math.PI);
        ctx.fillStyle = k ? th.draw : th.ring;
        ctx.fill();
        ctx.lineWidth = 2;
        ctx.strokeStyle = th.draw;
        ctx.stroke();
      });
    }
    if (d.flash) {
      if (now > d.flash.until) {
        d.flash = null;
      } else {
        const b = d.flash.b;
        ctx.globalAlpha = Math.min(1, (d.flash.until - now) / 400);
        if (d.flash.a) {
          ctx.beginPath();
          ctx.moveTo(sx(d.flash.a), sy(d.flash.a));
          ctx.lineTo(sx(b), sy(b));
          ctx.lineWidth = 3;
          ctx.setLineDash([6, 4]);
          ctx.strokeStyle = th.bad;
          ctx.stroke();
          ctx.setLineDash([]);
        }
        const X = sx(b), Y = sy(b);
        ctx.beginPath();
        ctx.moveTo(X - 6, Y - 6);
        ctx.lineTo(X + 6, Y + 6);
        ctx.moveTo(X + 6, Y - 6);
        ctx.lineTo(X - 6, Y + 6);
        ctx.lineWidth = 3;
        ctx.strokeStyle = th.bad;
        ctx.stroke();
        ctx.globalAlpha = 1;
      }
    }
    ctx.restore();
  }

  /* ---- side panel ------------------------------------------------------------ */

  function showWarning(msg) {
    els.warning.hidden = !msg;
    els.warning.textContent = msg || "";
  }

  function updatePanel() {
    const ctrl = state.ctrl;
    const lap = ctrl.lapTicks, back = ctrl.world === "polygon" && ctrl.sim.atEnd === "reverse";
    els.tick.textContent = "t = " + ctrl.t + " · " + (back ? "pass " : "lap ") + (Math.floor(ctrl.t / lap) + 1);
    document.body.classList.toggle("no-gt", !state.gt);

    let html = "";
    ctrl.shadows.forEach((s) => {
      const exact = s.lo === s.hi;
      const e = s.expected === null ? "—" : s.expected.toFixed(2);
      html += "<tr" + (s.ok ? "" : ' class="bad"') + ">" +
        '<td><span class="lab"><span class="chip" style="background:' + solidColor(s.label) + '"></span>s' + s.label + "</span></td>" +
        '<td class="col-true">' + s.truth + "</td>" +
        "<td" + (exact ? ' class="exact"' : "") + ">" + fmtBounds(s.lo, s.hi) + "</td>" +
        "<td>" + e + "</td></tr>";
    });
    if (!ctrl.shadows.length) html = '<tr class="empty"><td colspan="4">no shadows</td></tr>';
    els.rows.innerHTML = html;

    const h = ctrl.hidden;
    const rng = h.lo === h.hi ? String(h.lo) : h.hi === INF ? "≥ " + h.lo : h.lo + "–" + h.hi;
    els.aggregate.innerHTML = "Hidden targets: <b>" + rng + "</b>" +
      (state.gt ? ' <span class="muted">(true ' + h.truth + " of " + h.total + ")</span>" : "");

    const notes = [];
    if (ctrl.prob) {
      notes.push("E[n]: probabilistic filter, area-proportional binomial splits, " + ctrl.opts.truncation +
        " truncation to " + ctrl.opts.maxEntries + " entries (" + ctrl.prob.entries + " now).");
    } else {
      notes.push("E[n] unavailable: " + ctrl.probNote + ".");
    }
    let iso = "I-state: " + ctrl.edges + " edges.";
    const cp = ctrl.compactions;
    if (cp.exact || cp.relaxed) {
      iso += " Re-rooted at the current shadows " + (cp.exact + cp.relaxed) + "× (" +
        (cp.relaxed ? cp.relaxed + " relaxed, " : "") + cp.exact + " exact; last at t = " + cp.lastT + ").";
    }
    notes.push(iso);
    els.notes.textContent = notes.join(" ");

    const v = ctrl.violations;
    if (v.length) {
      const x = v[v.length - 1];
      showWarning("Bound violated (" + v.length + "×): t = " + x.t + ", " +
        (x.label === "all" ? "all shadows" : "s" + x.label) + " holds " + x.truth + " ∉ " + fmtBounds(x.lo, x.hi));
    } else if (!ctrl.errors.length) {
      showWarning("");
    }
    drawGraph();
  }

  /* ---- bipartite graph -------------------------------------------------------- */

  /* Text widths for the graph layout; GRAPH_FONT matches "#graph text" in style.css. */
  const GRAPH_FS = 11;
  const GRAPH_FONT = "system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif";
  const measureCtx = document.createElement("canvas").getContext("2d");

  function textWidth(text, bold) {
    measureCtx.font = (bold ? "600 " : "400 ") + GRAPH_FS + "px " + GRAPH_FONT;
    return measureCtx.measureText(text).width;
  }

  /* Drawn in CSS pixels (viewBox width = card width) so the text keeps its size; the
   * viewBox only grows past the card, and the graph scales down, when the nodes need it. */
  function drawGraph() {
    const view = state.ctrl.bipartiteView(14, 5);
    const NH = 18, GAP = 22, TOP = 26, PAD = 7, LINE = 16;
    const rooted = view.rootT !== null;
    const tagOf = (v) => (v.kind === "initial" ? (rooted ? "root" : "t₀") :
      v.kind === "appeared" ? "new" : v.kind === "gone" ? "gone" : "");
    const parts = (v) => (v.kind === "fov"
      ? { name: v.id === "enter" ? "+ enters" : "− exits", b: "Σ " + v.lo }
      : { name: "s" + v.id, b: fmtShort(v.lo, v.hi) });
    let NW = 96;
    view.left.concat(view.right).forEach((v) => {
      const p = parts(v), tag = tagOf(v);
      v.nameW = textWidth(p.name, true);
      NW = Math.max(NW, PAD + v.nameW + (tag ? 5 + textWidth(tag) : 0) + 10 + textWidth(p.b, true) + PAD);
    });
    NW = Math.ceil(NW);
    const headL = (rooted ? "root (t = " + view.rootT + ")" : "t₀") + " / appeared";
    const headR = "disappeared / now";
    const notes = [];
    if (view.hiddenLeft) notes.push(["+ " + view.hiddenLeft + " older left vertices not shown", false]);
    if (view.hiddenGone) notes.push(["+ " + view.hiddenGone + " earlier disappeared not shown", true]);
    const noteW = notes.map((n) => textWidth(n[0]));
    const minW = Math.ceil(Math.max(2 * NW + 48, textWidth(headL) + textWidth(headR) + 24, Math.max(0, ...noteW) + 8));
    const W = Math.max(Math.floor(els.graph.parentNode.clientWidth) || 360, minW);
    state.graphW = els.graph.parentNode.clientWidth;
    const oneLine = notes.length === 2 && noteW[0] + noteW[1] + 24 <= W;
    const noteLines = oneLine ? 1 : notes.length;
    const nL = view.left.length, nR = view.right.length;
    const rows = Math.max(nL, nR, 1);
    const H = TOP + rows * GAP + 4 + noteLines * LINE - (noteLines ? 4 : 0);
    const xL = 4, xR = W - NW - 4;
    const yOf = (i, n) => TOP + (rows - n) * GAP / 2 + i * GAP;
    let s = "";
    s += '<text x="' + xL + '" y="12" class="muted">' + esc(headL) + "</text>";
    s += '<text x="' + (xR + NW) + '" y="12" class="muted" text-anchor="end">' + headR + "</text>";
    view.edges.forEach((e) => {
      const r = view.right[e[1]];
      const y1 = yOf(e[0], nL) + NH / 2, y2 = yOf(e[1], nR) + NH / 2;
      const x1 = xL + NW, x2 = xR, mx = (x1 + x2) / 2;
      const style = r.kind === "alive" ? ' style="stroke:' + solidColor(r.id) + ';stroke-opacity:.55"' : "";
      s += '<path class="edge"' + style + ' d="M' + x1 + " " + y1 + "C" + mx + " " + y1 + " " + mx + " " + y2 + " " + x2 + " " + y2 + '"/>';
    });
    const node = (v, x, y) => {
      let cls = "node", fill = "";
      const p = parts(v), tag = tagOf(v);
      if (v.kind === "fov") cls += " fov";
      if (v.kind === "alive") fill = ' style="fill:' + shadowColor(v.id, 0.5) + ";stroke:" + solidColor(v.id) + '"';
      if (v.kind === "gone") cls += " gone";
      return '<rect class="' + cls + '"' + fill + ' x="' + x + '" y="' + y + '" width="' + NW + '" height="' + NH + '" rx="6"/>' +
        '<text x="' + (x + PAD) + '" y="' + (y + 13) + '">' + esc(p.name) + "</text>" +
        (tag ? '<text x="' + (x + PAD + v.nameW + 5).toFixed(1) + '" y="' + (y + 13) + '" class="muted">' + tag + "</text>" : "") +
        '<text x="' + (x + NW - PAD) + '" y="' + (y + 13) + '" text-anchor="end">' + esc(p.b) + "</text>";
    };
    view.left.forEach((v, i) => { s += node(v, xL, yOf(i, nL)); });
    view.right.forEach((v, i) => { s += node(v, xR, yOf(i, nR)); });
    notes.forEach((n, k) => {
      const y = H - 4 - (oneLine ? 0 : (notes.length - 1 - k) * LINE);
      s += n[1] ? '<text x="' + (W - 4) + '" y="' + y + '" class="muted" text-anchor="end">' + n[0] + "</text>"
        : '<text x="' + xL + '" y="' + y + '" class="muted">' + n[0] + "</text>";
    });
    els.graph.setAttribute("viewBox", "0 0 " + W + " " + H);
    els.graph.innerHTML = s;
  }

  /* ---- event log --------------------------------------------------------------- */

  const ICONS = {
    appear: '<circle cx="8" cy="8" r="5.5"/><path d="M8 5.2v5.6M5.2 8h5.6"/>',
    disappear: '<circle cx="8" cy="8" r="5.5" stroke-dasharray="2.4 2"/><path d="M5.2 8h5.6"/>',
    split: '<path d="M1.5 8h6M7.5 8l6.5-5M7.5 8l6.5 5"/>',
    merge: '<path d="M2 3l6.5 5L2 13M8.5 8h6"/>',
    enter: '<path d="M1.5 8h8.5M7 4.8L10.2 8 7 11.2M13.5 2.5v11"/>',
    exit: '<path d="M5.5 8h9M11.5 4.8L14.7 8l-3.2 3.2M2.5 2.5v11"/>',
  };
  const logTicks = [];

  function logEvents(t, events) {
    if (!events.length) return;
    const frag = document.createDocumentFragment();
    const head = document.createElement("li");
    head.className = "tick-head";
    head.textContent = "t = " + t;
    frag.appendChild(head);
    const items = [head];
    events.forEach((e) => {
      const li = document.createElement("li");
      li.className = "ev ev-" + e.type;
      li.innerHTML = '<span class="ev-ico" title="' + e.type + '"><svg viewBox="0 0 16 16">' + ICONS[e.type] + "</svg></span>" +
        "<span>" + esc(SI.demo.describeEvent(e)) + "</span>";
      frag.appendChild(li);
      items.push(li);
    });
    els.log.insertBefore(frag, els.log.firstChild);
    logTicks.push(items);
    while (logTicks.length > LOG_TICKS) logTicks.shift().forEach((li) => li.remove());
  }

  /* ---- loop and wiring -------------------------------------------------------- */

  function frame(now) {
    if (state.playing) {
      const interval = 1000 / state.speed;
      let steps = 0;
      while (now - state.stepAt >= interval && steps < 3) {
        const at = now - state.stepAt > 2 * interval ? now : state.stepAt + interval;
        if (!doStep(at)) break;
        steps++;
      }
      if (steps) updatePanel();
    }
    if (state.needsDraw || state.dirtyCells || state.playing) draw(now);
    requestAnimationFrame(frame);
  }

  function restart() {
    const was = state.playing;
    build();
    setPlaying(was);
  }

  const more = $("intro-more");
  more.addEventListener("click", () => {
    const open = $("intro").classList.toggle("open");
    more.textContent = open ? "Show less" : "Read more";
    more.setAttribute("aria-expanded", open ? "true" : "false");
  });
  els.play.addEventListener("click", () => setPlaying(!state.playing));
  els.step.addEventListener("click", stepOnce);
  els.reset.addEventListener("click", restart);
  els.speed.addEventListener("input", () => {
    state.speed = +els.speed.value;
    els.speedOut.textContent = state.speed + "/s";
  });
  [els.n, els.seed, els.mode, els.radius, els.path].forEach((el) => el.addEventListener("change", restart));
  els.map.addEventListener("change", () => {
    showWarning("");
    buildPathSelect("default");
    restart();
  });
  els.gt.addEventListener("change", () => {
    state.gt = els.gt.checked;
    state.needsDraw = true;
    updatePanel();
  });

  els.draw.addEventListener("click", () => (state.draw ? stopDrawing() : startDrawing()));
  els.drawUndo.addEventListener("click", () => {
    if (!state.draw) return;
    state.draw.points.pop();
    state.draw.flash = null;
    drawMessage(state.draw.points.length ? state.draw.points.length + " waypoints." : DRAW_HELP, false);
    state.needsDraw = true;
  });
  els.drawCancel.addEventListener("click", stopDrawing);
  els.drawStart.addEventListener("click", finishDrawing);
  els.canvas.addEventListener("click", (ev) => {
    if (state.draw) addWaypoint(ev);
  });
  els.canvas.addEventListener("dblclick", (ev) => {
    if (!state.draw) return;
    ev.preventDefault();
    finishDrawing();
  });
  els.canvas.addEventListener("pointermove", (ev) => {
    if (!state.draw) return;
    const s = canvasPoint(ev);
    state.draw.hover = toWorld(s[0], s[1]);
    state.needsDraw = true;
  });
  els.canvas.addEventListener("pointerleave", () => {
    if (!state.draw) return;
    state.draw.hover = null;
    state.needsDraw = true;
  });

  /* Space and ArrowRight are page shortcuts unless the focused element uses the key itself:
   * Space activates buttons, links and checkboxes; arrows move sliders, selects and carets.
   * While drawing a path, Enter starts it, Escape cancels and Backspace removes the last waypoint. */
  function ownsKey(target, key) {
    if (!target || !target.closest) return false;
    const sel = key === " " || key === "Enter" ? "input, select, textarea, button, a[href], summary, [contenteditable]"
      : "input, select, textarea, [contenteditable]";
    return !!target.closest(sel);
  }

  document.addEventListener("keydown", (ev) => {
    if (ev.altKey || ev.ctrlKey || ev.metaKey) return;
    if (state.draw) {
      if (ev.key === "Escape") {
        ev.preventDefault();
        stopDrawing();
      } else if (ev.key === "Enter" && !ownsKey(ev.target, "Enter")) {
        ev.preventDefault();
        finishDrawing();
      } else if (ev.key === "Backspace" && !ownsKey(ev.target, ev.key)) {
        ev.preventDefault();
        els.drawUndo.click();
      }
      return;
    }
    const space = ev.code === "Space" || ev.key === " ";
    if (!space && ev.key !== "ArrowRight") return;
    if (ownsKey(ev.target, space ? " " : ev.key)) return;
    ev.preventDefault();
    if (space) {
      if (!ev.repeat) setPlaying(!state.playing);
    } else {
      stepOnce();
    }
  });

  const dark = window.matchMedia ? window.matchMedia("(prefers-color-scheme: dark)") : null;
  if (dark && dark.addEventListener) {
    dark.addEventListener("change", () => {
      readTheme();
      state.dirtyCells = true;
      updatePanel();
    });
  }
  if (window.ResizeObserver) {
    new ResizeObserver(() => sizeCanvas()).observe(els.wrap);
    new ResizeObserver(() => {
      if (els.graph.parentNode.clientWidth !== state.graphW) drawGraph();
    }).observe(els.graph.parentNode);
  }
  window.addEventListener("resize", sizeCanvas);

  buildMapSelect();
  applyParams();
  readTheme();
  state.gt = els.gt.checked;
  state.speed = +els.speed.value;
  els.speedOut.textContent = state.speed + "/s";
  build();
  const ff = clampInt(params.get("ticks"), 0, 100000, 0);
  for (let k = 0; k < ff; k++) {
    if (!doStep(0, k >= ff - 40)) break;
  }
  state.prev = null;
  updatePanel();
  setPlaying(params.get("autoplay") === "1");
  draw(performance.now());
  requestAnimationFrame(frame);
})(globalThis.SI);
