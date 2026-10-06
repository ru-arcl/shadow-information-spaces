import geometry.Algorithm;
import geometry.FileHelper;
import geometry.Geometry;
import geometry.cut.Bitangent;
import geometry.cut.Cut;
import geometry.cut.GeneralInflection;
import geometry.cut.PathCutIntersectPoint;
import geometry.drawable.Path;
import geometry.drawable.Polygon;
import geometry.gap.Gap;
import geometry.gap.GapState;
import geometry.gap.PhysicalGap;
import geometry.graph.BipartiteGraph;
import geometry.graph.Edge;
import geometry.graph.Vertex;
import geometry.pe.is.Equation;
import geometry.pe.is.Oracle;
import geometry.pe.is.basic.NSEState;
import geometry.pe.is.singleTypeAgent.SingleTypeAgentAlgorithm;
import geometry.pe.is.singleTypeAgent.SingleTypeAgentEvent;
import geometry.pe.is.singleTypeAgent.SingleTypeAgentOracle;
import geometry.pe.is.singleTypeAgent.SingleTypeAgentState;

import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintStream;
import java.io.Writer;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Headless golden-output harness for the original Java implementation of
 * Yu & LaValle, "Shadow Information Spaces" (T-RO 2012, ICRA 2008).
 *
 * It calls the non-UI classes of the original code (package geometry, unmodified)
 * the way the ProjectPanel* constructors do, but with no AWT window, and
 * writes JSON fixtures for the Python port. See README.md in this folder.
 *
 * Usage:
 *   java -Djava.awt.headless=true --add-opens java.base/java.lang=ALL-UNNAMED \
 *        -cp <classes>:<log4j.jar> GoldenDump <orig/source dir> <paths.tsv> <out dir> [options]
 * Options:
 *   --polys 1,2,...       polygons to dump (default 1..14)
 *   --runs NAME@POLY,...  path runs to dump (default: every row of paths.tsv)
 *   --seeds 1,2,3,4,5     seeds for the original oracle
 *   --fixed-seeds 1,2,3   seeds for the merge-fixed oracle variant
 *   --burnins 0,1,2,...   identity-hash burn-ins for the bounds order check
 *   --no-polys / --no-runs
 *
 * What the harness adds beyond calling the original code:
 *   1. It seeds Math.random(). Math.random() delegates to one static
 *      java.util.Random held in java.lang.Math$RandomNumberGeneratorHolder.
 *      The harness takes that instance by reflection, which is why the java
 *      command needs --add-opens, and calls setSeed(seed). The original
 *      SingleTypeAgentOracle then runs unmodified with a reproducible
 *      stream: java.util.Random(seed).nextDouble().
 *   2. It re-implements the sampling loop at the top of Algorithm.getGaps
 *      (getAllCriticalPoints, then purturbPointAlongSeg, then
 *      getPhysicalGaps). All three are public. Doing this exposes the
 *      per-sample physical gaps, which getGaps keeps to itself. The ID
 *      propagation is not re-implemented: Algorithm.getGaps itself is called
 *      for that.
 *   3. RecordingPath, a subclass of Path, logs the points passed to
 *      ptDistFromStart(Point2D). Path.java is not modified. getGaps calls this
 *      method at i == 0 and at every GI/BT event, so the log locates the
 *      critical point at which a crashing getGaps fails.
 *   4. It runs deriveShadowBounds several times, each after a different
 *      number of System.identityHashCode "burn-in" calls. Its result depends
 *      on HashMap<Vertex,...> iteration order, which follows identity
 *      hashes. The distinct outcomes are recorded.
 */
public class GoldenDump {

    static final String[] TYPE = {"GENERAL_INFLECTION", "NONGENERAL_INFLECTION", "SINGLETANGENT", "BITANGENT", "NONE"};

    // ------------------------------------------------------------------ options
    static int[] polys = range(1, 14);
    static List<String> runFilter = null;
    static long[] seeds = {1, 2, 3, 4, 5};
    static long[] fixedSeeds = {1, 2, 3};
    static int[] burnins = range(0, 23);
    static boolean doPolys = true, doRuns = true;

    static int[] range(int a, int b) { int[] r = new int[b - a + 1]; for (int i = a; i <= b; i++) r[i - a] = i; return r; }

    static int[] ints(String s) { String[] p = s.split(","); int[] r = new int[p.length]; for (int i = 0; i < p.length; i++) r[i] = Integer.parseInt(p[i].trim()); return r; }

    static long[] longs(String s) { if (s.isEmpty()) return new long[0]; String[] p = s.split(","); long[] r = new long[p.length]; for (int i = 0; i < p.length; i++) r[i] = Long.parseLong(p[i].trim()); return r; }

    public static void main(String[] a) throws Exception {
        if (a.length < 3) {
            System.err.println("usage: GoldenDump <orig/source dir> <paths.tsv> <out dir> [options]");
            System.exit(2);
        }
        File src = new File(a[0]).getAbsoluteFile();
        File pathsFile = new File(a[1]);
        File out = new File(a[2]);
        for (int i = 3; i < a.length; i++) {
            if (a[i].equals("--polys")) polys = ints(a[++i]);
            else if (a[i].equals("--runs")) runFilter = Arrays.asList(a[++i].split(","));
            else if (a[i].equals("--seeds")) seeds = longs(a[++i]);
            else if (a[i].equals("--fixed-seeds")) fixedSeeds = longs(a[++i]);
            else if (a[i].equals("--burnins")) burnins = ints(a[++i]);
            else if (a[i].equals("--no-polys")) doPolys = false;
            else if (a[i].equals("--no-runs")) doRuns = false;
            else throw new IllegalArgumentException(a[i]);
        }
        // Must happen before the static initialisers of Geometry and Polygon run.
        FileHelper.setBaseUrl(src.toURI().toString());
        out.mkdirs();

        Map<String, Object> index = Json.obj();
        index.put("generator", "tools/java_reference/src/GoldenDump.java");
        index.put("java_version", System.getProperty("java.version"));
        index.put("java_vm", System.getProperty("java.vm.name") + " " + System.getProperty("java.vm.version"));
        index.put("epsilon", Algorithm.epsilon);
        index.put("purturb", Algorithm.purturb);
        index.put("seeds_original", seeds);
        index.put("seeds_merge_fixed", fixedSeeds);
        index.put("burnins", burnins);
        // RNG check: after seedMathRandom(s), Math.random() must equal new java.util.Random(s).nextDouble()
        {
            Map<String, Object> rc = Json.obj();
            seedMathRandom(1);
            double[] m = new double[6];
            for (int i = 0; i < m.length; i++) m[i] = Math.random();
            Random ref = new Random(1);
            boolean same = true;
            for (int i = 0; i < m.length; i++) same &= m[i] == ref.nextDouble();
            rc.put("seed", 1);
            rc.put("math_random_after_seeding", m);
            rc.put("equals_new_Random_seed_nextDouble", same);
            int[] ints = new int[6];
            Random r2 = new Random(42);
            for (int i = 0; i < ints.length; i++) ints[i] = (int) (r2.nextDouble() * 100);
            rc.put("new_Random_42_int_nextDouble_x100", ints);
            index.put("java_random_check", rc);
            if (!same) throw new IllegalStateException("seeding Math.random() did not work");
        }

        if (doPolys) {
            List<Object> pl = new ArrayList<Object>();
            for (int n : polys) {
                long t0 = System.currentTimeMillis();
                Map<String, Object> d = dumpPolygon(n);
                String fn = "poly" + n + ".json";
                writeFile(new File(out, fn), d);
                pl.add(fn);
                System.err.println("polygon " + n + ": " + (System.currentTimeMillis() - t0) + " ms");
            }
            index.put("polygon_files", pl);
        }
        if (doRuns) {
            List<Object> rl = new ArrayList<Object>();
            BufferedReader br = new BufferedReader(new FileReader(pathsFile));
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty() || line.startsWith("#")) continue;
                String[] f = line.split("\t");
                String name = f[0];
                int poly = Integer.parseInt(f[1]);
                int nTargets = Integer.parseInt(f[2]);
                if (runFilter != null && !runFilter.contains(name + "@" + poly)) continue;
                List<double[]> wps = new ArrayList<double[]>();
                for (String xy : f[4].trim().split(" +")) {
                    String[] c = xy.split(",");
                    wps.add(new double[]{Double.parseDouble(c[0]), Double.parseDouble(c[1])});
                }
                long t0 = System.currentTimeMillis();
                Map<String, Object> d = dumpRun(name, poly, nTargets, f[3], wps);
                String fn = poly + "_" + name + ".json";
                writeFile(new File(out, fn), d);
                Map<String, Object> r = Json.obj();
                r.put("file", fn);
                r.put("name", name);
                r.put("polygon", poly);
                r.put("source", f[3]);
                r.put("get_gaps_ok", ((Map<?, ?>) d.get("get_gaps")).get("ok"));
                rl.add(r);
                System.err.println("run " + name + "@" + poly + ": " + (System.currentTimeMillis() - t0) + " ms");
            }
            br.close();
            index.put("runs", rl);
        }
        index.put("non_round_trip_doubles", Json.nonRoundTrip);
        if (doPolys && doRuns && runFilter == null) writeFile(new File(out, "index.json"), index);
        System.err.println("non-round-trip doubles: " + Json.nonRoundTrip);
    }

    static void writeFile(File f, Object d) throws Exception {
        Writer w = new FileWriter(f);
        w.write(Json.write(d));
        w.close();
    }

    // ------------------------------------------------------------------ helpers
    static double[] pt(Point2D p) { return p == null ? null : new double[]{p.getX(), p.getY()}; }

    static double[] ln(Line2D l) { return l == null ? null : new double[]{l.getX1(), l.getY1(), l.getX2(), l.getY2()}; }

    static boolean sameLine(Line2D a, Line2D b) {
        return a.getX1() == b.getX1() && a.getY1() == b.getY1() && a.getX2() == b.getX2() && a.getY2() == b.getY2();
    }

    static List<Object> lines(Line2D[] ls) { List<Object> r = new ArrayList<Object>(); for (Line2D l : ls) r.add(ln(l)); return r; }

    static List<Object> points(Point2D[] ps) { List<Object> r = new ArrayList<Object>(); for (Point2D p : ps) r.add(pt(p)); return r; }

    static Map<String, Object> exc(Throwable t) {
        Map<String, Object> m = Json.obj();
        m.put("class", t.getClass().getName());
        m.put("message", t.getMessage());
        List<Object> st = new ArrayList<Object>();
        StackTraceElement[] s = t.getStackTrace();
        for (int i = 0; i < Math.min(8, s.length); i++) st.add(s[i].toString());
        m.put("stack", st);
        return m;
    }

    static List<Object> pgap(PhysicalGap g) {
        return Json.arr(g.getStartEdge(), g.getEndEdge(), pt(g.getStartPoint()), pt(g.getEndPoint()));
    }

    static List<Object> pgaps(PhysicalGap[] gs) { List<Object> r = new ArrayList<Object>(); for (PhysicalGap g : gs) r.add(pgap(g)); return r; }

    static String bits(boolean[] b) { StringBuilder sb = new StringBuilder(); for (boolean x : b) sb.append(x ? '1' : '0'); return sb.toString(); }

    static String typeName(Algorithm.CUT_TYPE t) { return t == null ? null : t.name(); }

    /** Path that records the points passed to the public ptDistFromStart(Point2D). */
    static class RecordingPath extends Path {
        final List<Point2D> calls = new ArrayList<Point2D>();
        RecordingPath() { super(Geometry.DEFAULT_DC); }
        @Override public double ptDistFromStart(Point2D p) { calls.add(p); return super.ptDistFromStart(p); }
    }

    /** Seed the java.util.Random behind Math.random(), see the class comment. */
    static void seedMathRandom(long seed) throws Exception {
        Class<?> h = Class.forName("java.lang.Math$RandomNumberGeneratorHolder");
        Field f = h.getDeclaredField("randomNumberGenerator");
        f.setAccessible(true);
        ((Random) f.get(null)).setSeed(seed);
    }

    static final PrintStream REAL_OUT = System.out;

    interface Body { void run() throws Exception; }

    /** Run body with System.out captured; returns the captured text. */
    static String capture(Body b) throws Exception {
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        PrintStream ps = new PrintStream(bos, true, "UTF-8");
        System.setOut(ps);
        try { b.run(); } finally { System.setOut(REAL_OUT); ps.flush(); }
        return bos.toString("UTF-8");
    }

    static void burn(int n) { for (int i = 0; i < n; i++) System.identityHashCode(new Object()); }

    // ------------------------------------------------------------------ polygons
    static Map<String, Object> dumpPolygon(int n) {
        Map<String, Object> d = Json.obj();
        Polygon p = Polygon.getPolygon(n, Geometry.DEFAULT_DC);
        Line2D[] L = p.getLineSegmentArray();
        int nv = L.length;
        d.put("polygon", n);
        d.put("file", "polygons/" + n + ".dat");
        d.put("n", nv);
        d.put("vertices", points(p.getPointArray()));
        // turn[i] = lines[i-1].relativeCCW(v_{i+1}); reflex <=> +1 (test used by every cut function)
        int[] turn = new int[nv];
        List<Object> reflex = new ArrayList<Object>();
        for (int i = 0; i < nv; i++) {
            Line2D prev = L[i == 0 ? nv - 1 : i - 1];
            turn[i] = prev.relativeCCW(L[i].getP2());
            if (turn[i] == 1) reflex.add(i);
        }
        d.put("turn", turn);
        d.put("reflex", reflex);

        Map<String, Object> cuts = Json.obj();
        Map<String, Object> errs = Json.obj();
        Line2D[] infAll = null, infGen = null, infNon = null, bitLines = null;
        Cut[] st = null;
        GeneralInflection[] gic = null;
        Bitangent[] bt = null;
        try { st = Algorithm.getSingletangentCut(p); List<Object> l = new ArrayList<Object>(); for (Cut c : st) l.add(ln(c.getCut())); cuts.put("single_tangent", l); } catch (Throwable t) { errs.put("single_tangent", exc(t)); }
        try { infNon = Algorithm.getInflection(p, Algorithm.INFLECTION_TYPE.NONGENERAL); cuts.put("inflection_nongeneral", lines(infNon)); } catch (Throwable t) { errs.put("inflection_nongeneral", exc(t)); }
        try { infGen = Algorithm.getInflection(p, Algorithm.INFLECTION_TYPE.GENERAL); } catch (Throwable t) { errs.put("inflection_general", exc(t)); }
        try { infAll = Algorithm.getInflection(p, Algorithm.INFLECTION_TYPE.ALL); cuts.put("inflection_all", lines(infAll)); } catch (Throwable t) { errs.put("inflection_all", exc(t)); }
        try {
            gic = Algorithm.getGeneralInflectionCut(p);
            List<Object> l = new ArrayList<Object>();
            for (GeneralInflection g : gic) l.add(Json.arr(ln(g.getCut()), g.getFromLine(), g.isCounterClockWise()));
            cuts.put("general_inflection_cut", l);
        } catch (Throwable t) { errs.put("general_inflection_cut", exc(t)); }
        try {
            bt = Algorithm.getBitangentCut(p);
            List<Object> l = new ArrayList<Object>();
            for (Bitangent b : bt) l.add(Json.arr(ln(b.getCut()), b.getThisPoint(), b.getOppositePoint(), ln(b.getOppositeSegment()), pt(b.getCurveToPoint())));
            cuts.put("bitangent_cut", l);
        } catch (Throwable t) { errs.put("bitangent_cut", exc(t)); }
        try { bitLines = Algorithm.getBitangent(p); } catch (Throwable t) { errs.put("bitangent_lines", exc(t)); }
        cuts.put("_format", "single_tangent/inflection_*: [x1,y1,x2,y2] (ray from the tangent/reflex vertex P1 to the boundary hit P2); "
                + "general_inflection_cut: [[x1,y1,x2,y2], fromLine, counterClockWise]; "
                + "bitangent_cut: [[cut], thisPoint, oppositePoint, [oppositeSegment], [curveToPoint]]");
        // getInflection(GENERAL) and getBitangent() are expected to repeat the lines of
        // getGeneralInflectionCut / getBitangentCut; store them only if they differ.
        if (infGen != null && gic != null) {
            boolean same = infGen.length == gic.length;
            for (int i = 0; same && i < gic.length; i++) same = sameLine(infGen[i], gic[i].getCut());
            cuts.put("inflection_general_equals_general_inflection_cut_lines", same);
            if (!same) cuts.put("inflection_general", lines(infGen));
        }
        if (bitLines != null && bt != null) {
            boolean same = bitLines.length == bt.length;
            for (int i = 0; same && i < bt.length; i++) same = sameLine(bitLines[i], bt[i].getCut());
            cuts.put("bitangent_lines_equal_bitangent_cut_lines", same);
            if (!same) cuts.put("bitangent_lines", lines(bitLines));
        }
        try {
            Cut[] all = Algorithm.getCuts(p);
            List<Object> runs = new ArrayList<Object>();
            String cur = null; int cnt = 0;
            for (Cut c : all) {
                String t = typeName(c.getCutType());
                if (!t.equals(cur)) { if (cur != null) runs.add(Json.arr(cur, cnt)); cur = t; cnt = 0; }
                cnt++;
            }
            if (cur != null) runs.add(Json.arr(cur, cnt));
            cuts.put("get_cuts_layout", runs);
            cuts.put("get_cuts_total", all.length);
        } catch (Throwable t) { errs.put("get_cuts", exc(t)); }
        Map<String, Object> counts = Json.obj();
        counts.put("reflex", reflex.size());
        counts.put("single_tangent", st == null ? null : st.length);
        counts.put("inflection_nongeneral", infNon == null ? null : infNon.length);
        counts.put("inflection_general", infGen == null ? null : infGen.length);
        counts.put("inflection_all", infAll == null ? null : infAll.length);
        counts.put("general_inflection_cut", gic == null ? null : gic.length);
        counts.put("bitangent_cut", bt == null ? null : bt.length);
        counts.put("bitangent_lines", bitLines == null ? null : bitLines.length);
        d.put("counts", counts);
        d.put("cuts", cuts);

        // point-in-polygon (Path2D.Float, WIND_NON_ZERO) on a lattice, at the vertices and at edge midpoints
        Map<String, Object> pip = Json.obj();
        pip.put("lattice", "x,y in 50,100,...,950; rows by y ascending, chars by x ascending");
        List<Object> rows = new ArrayList<Object>();
        for (int yi = 1; yi <= 19; yi++) {
            boolean[] b = new boolean[19];
            for (int xi = 1; xi <= 19; xi++) b[xi - 1] = p.pointInPolygon(new Point2D.Double(50 * xi, 50 * yi));
            rows.add(bits(b));
        }
        pip.put("rows", rows);
        boolean[] atV = new boolean[nv], atMid = new boolean[nv];
        for (int i = 0; i < nv; i++) {
            atV[i] = p.pointInPolygon(L[i].getP1());
            atMid[i] = p.pointInPolygon(new Point2D.Double((L[i].getX1() + L[i].getX2()) / 2, (L[i].getY1() + L[i].getY2()) / 2));
        }
        pip.put("at_vertices", bits(atV));
        pip.put("at_edge_midpoints", bits(atMid));
        d.put("point_in_polygon", pip);

        // visibility from the interior points of a 100-unit lattice
        List<Object> vis = new ArrayList<Object>();
        List<Point2D> inside = new ArrayList<Point2D>();
        for (int yi = 1; yi <= 9; yi++)
            for (int xi = 1; xi <= 9; xi++) {
                Point2D q = new Point2D.Double(100 * xi, 100 * yi);
                if (p.pointInPolygon(q)) inside.add(q);
            }
        int every = Math.max(1, (inside.size() + 11) / 12);
        for (int k = 0; k < inside.size(); k++) {
            Point2D q = inside.get(k);
            Map<String, Object> v = Json.obj();
            v.put("q", pt(q));
            boolean[] sip = new boolean[nv];
            for (int j = 0; j < nv; j++) sip[j] = Algorithm.segInPolygon(p, new Line2D.Double(q, L[j].getP1()));
            v.put("seg_in_polygon_to_vertices", bits(sip));
            try { v.put("physical_gaps", pgaps(Algorithm.getPhysicalGaps(p, q))); } catch (Throwable t) { v.put("physical_gaps_error", exc(t)); }
            if (k % every == 0) {
                try { v.put("visibility_polygon", points(Algorithm.getVisibilityPolygon(p, q, Geometry.DEFAULT_DC).getPointArray())); }
                catch (Throwable t) { v.put("visibility_polygon_error", exc(t)); }
            }
            vis.add(v);
        }
        d.put("visibility_format", "physical_gaps: [startEdge, endEdge, startPoint|null, endPoint|null] in Java output order; "
                + "seg_in_polygon_to_vertices: char j = segInPolygon(q -> v_j); visibility_polygon (every few points only): getVisibilityPolygon vertex list incl. duplicates");
        d.put("visibility", vis);
        if (!errs.isEmpty()) d.put("errors", errs);
        return d;
    }

    // ------------------------------------------------------------------ runs
    static RecordingPath makePath(List<double[]> wps) {
        RecordingPath path = new RecordingPath();
        for (double[] w : wps) path.addPoint(new Point2D.Double(w[0], w[1]));
        return path;
    }

    static List<String> gapLines(Gap[][] gapss, boolean withTime) {
        List<String> r = new ArrayList<String>();
        for (Gap[] gs : gapss) {
            StringBuilder sb = new StringBuilder();
            if (withTime) sb.append(gs[0].getRelativeTime()).append(' ');
            for (Gap g : gs) sb.append(g.toString()).append(' ');
            r.add(sb.toString());
        }
        return r;
    }

    static Map<String, Object> dumpRun(String name, int polyNo, int nTargets, String source, List<double[]> wps) throws Exception {
        // Runs on digitised figure paths are for comparison with the paper only:
        // they skip the per-sample physical gaps and use one seed per oracle variant.
        boolean digitised = name.startsWith("fig_");
        boolean withSampleGaps = !digitised;
        Map<String, Object> d = Json.obj();
        d.put("name", name);
        d.put("polygon", polyNo);
        d.put("source", source);
        d.put("n_targets", nTargets);
        d.put("waypoints", wps);
        Polygon poly = Polygon.getPolygon(polyNo, Geometry.DEFAULT_DC);
        RecordingPath path = makePath(wps);
        boolean[] inside = new boolean[wps.size()];
        for (int i = 0; i < wps.size(); i++) inside[i] = poly.pointInPolygon(new Point2D.Double(wps.get(i)[0], wps.get(i)[1]));
        d.put("waypoints_inside", bits(inside));
        d.put("path_length", path.getLength());

        // ---- critical points + per-sample physical gaps (the sampling loop of getGaps)
        Cut[] cuts = Algorithm.getCuts(poly);
        Line2D[] segs = path.getLineArray();
        PathCutIntersectPoint[] pc = null;
        try { pc = path.getAllCriticalPoints(cuts); } catch (Throwable t) { d.put("critical_points_error", exc(t)); }
        List<Object> cps = new ArrayList<Object>();
        List<Point2D> samplePts = new ArrayList<Point2D>();
        List<PhysicalGap[]> samplePg = new ArrayList<PhysicalGap[]>();
        int nEvents = 0;
        if (pc != null) {
            for (PathCutIntersectPoint c : pc) {
                int ci = -1, si = -1;
                for (int k = 0; k < cuts.length; k++) if (cuts[k] == c.getCut()) { ci = k; break; }
                for (int k = 0; k < segs.length; k++) if (segs[k] == c.getSeg()) { si = k; break; }
                if (c.getCutType() == Algorithm.CUT_TYPE.GENERAL_INFLECTION || c.getCutType() == Algorithm.CUT_TYPE.BITANGENT) nEvents++;
                List<Object> row = Json.arr(c.getDistance(), pt(c.getPoint()), typeName(c.getCutType()), ci, si);
                Point2D pp = null;
                try { pp = Algorithm.purturbPointAlongSeg(c.getPoint(), c.getSeg(), Algorithm.purturb, true); row.add(pt(pp)); }
                catch (Throwable t) { row.add(null); row.add(exc(t)); }
                PhysicalGap[] spg = null;
                if (pp != null) {
                    try { spg = Algorithm.getPhysicalGaps(poly, pp); if (withSampleGaps) row.add(pgaps(spg)); } catch (Throwable t) { row.add(null); row.add(exc(t)); }
                }
                samplePts.add(pp);
                samplePg.add(spg);
                cps.add(row);
            }
        }
        d.put("critical_points_format", "[distance, point, cutType, index into getCuts(polygon) (-1 for NONE), path segment index, "
                + "perturbed sample point = purturbPointAlongSeg(point, seg, purturb, far=true)"
                + (withSampleGaps ? ", getPhysicalGaps(polygon, sample) as [startEdge,endEdge,startPoint,endPoint]]" : "] (per-sample physical gaps omitted for this run to save space)"));
        d.put("n_critical_points", cps.size());
        d.put("n_gi_bt_crossings", nEvents);
        d.put("critical_points", cps);

        // ---- what the panels compute for drawing
        Map<String, Object> pi = Json.obj();
        try { pi.put("bitangent_lines", points(path.getIntersectPoints(Algorithm.getBitangent(poly)))); } catch (Throwable t) { pi.put("bitangent_lines_error", exc(t)); }
        try { pi.put("general_inflection_lines", points(path.getIntersectPoints(Algorithm.getInflection(poly, Algorithm.INFLECTION_TYPE.GENERAL)))); } catch (Throwable t) { pi.put("general_inflection_lines_error", exc(t)); }
        d.put("panel_intersect_points", pi);

        // ---- visibility at waypoints (inside only) and at <= 5 samples
        List<Object> vis = new ArrayList<Object>();
        List<Point2D> vq = new ArrayList<Point2D>();
        List<String> vlab = new ArrayList<String>();
        for (int i = 0; i < wps.size(); i++) if (inside[i]) { vq.add(new Point2D.Double(wps.get(i)[0], wps.get(i)[1])); vlab.add("waypoint " + i); }
        int ns = samplePts.size();
        for (int k = 0; k < 5 && ns > 0; k++) {
            int idx = (int) ((long) k * (ns - 1) / 4);
            if (samplePts.get(idx) != null && poly.pointInPolygon(samplePts.get(idx))) { vq.add(samplePts.get(idx)); vlab.add("sample " + idx); }
        }
        for (int i = 0; i < vq.size(); i++) {
            Map<String, Object> v = Json.obj();
            v.put("where", vlab.get(i));
            v.put("q", pt(vq.get(i)));
            try { v.put("physical_gaps", pgaps(Algorithm.getPhysicalGaps(poly, vq.get(i)))); } catch (Throwable t) { v.put("physical_gaps_error", exc(t)); }
            try { v.put("visibility_polygon", points(Algorithm.getVisibilityPolygon(poly, vq.get(i), Geometry.DEFAULT_DC).getPointArray())); }
            catch (Throwable t) { v.put("visibility_polygon_error", exc(t)); }
            vis.add(v);
        }
        d.put("visibility", vis);

        // ---- Algorithm.getGaps
        Map<String, Object> gg = Json.obj();
        d.put("get_gaps", gg);
        Gap[][] gapss;
        path.calls.clear();
        try {
            gapss = Algorithm.getGaps(poly, path);
        } catch (Throwable t) {
            gg.put("ok", false);
            gg.put("exception", exc(t));
            gg.put("failure", locateFailure(t, pc, path.calls));
            return d;
        }
        gg.put("ok", true);
        List<String> lines0 = gapLines(gapss, true);
        gg.put("n_sets", gapss.length);
        int maxId = 0;
        boolean onlyFirstTimed = true;
        // Set k is the physical-gap array of the last sample of group k, i.e. the
        // sample at critical index closeIdx[k] (the index before the (k+1)-th GI/BT
        // crossing, or the last index). Check that, then store only id/iEdge/links.
        List<Integer> closeIdx = new ArrayList<Integer>();
        for (int i = 0; i < pc.length; i++) {
            if (i == pc.length - 1 || pc[i + 1].getCutType() == Algorithm.CUT_TYPE.GENERAL_INFLECTION
                    || pc[i + 1].getCutType() == Algorithm.CUT_TYPE.BITANGENT) closeIdx.add(i);
        }
        int sampleMismatch = closeIdx.size() == gapss.length ? 0 : -1;
        List<Object> sets = new ArrayList<Object>();
        for (int k = 0; k < gapss.length; k++) {
            Map<String, Object> s = Json.obj();
            s.put("t", gapss[k][0].getRelativeTime());
            if (sampleMismatch >= 0) {
                s.put("sample_index", closeIdx.get(k));
                PhysicalGap[] spg = samplePg.get(closeIdx.get(k));
                boolean same = spg != null && spg.length == gapss[k].length;
                for (int j = 0; same && j < spg.length; j++) {
                    PhysicalGap g = (PhysicalGap) gapss[k][j];
                    same = g.getStartEdge() == spg[j].getStartEdge() && g.getEndEdge() == spg[j].getEndEdge()
                            && java.util.Objects.equals(g.getStartPoint(), spg[j].getStartPoint())
                            && java.util.Objects.equals(g.getEndPoint(), spg[j].getEndPoint());
                }
                if (!same) sampleMismatch++;
            }
            List<Object> gl = new ArrayList<Object>();
            for (int j = 0; j < gapss[k].length; j++) {
                PhysicalGap g = (PhysicalGap) gapss[k][j];
                if (j > 0 && g.getRelativeTime() != 0.0) onlyFirstTimed = false;
                maxId = Math.max(maxId, g.getId());
                List<Object> to = new ArrayList<Object>(g.getToGapSet()); // Java HashSet iteration order
                gl.add(Json.arr(g.getId(), g.getIEdge(), to));
            }
            s.put("gaps", gl);
            s.put("line", lines0.get(k));
            sets.add(s);
        }
        gg.put("max_id", maxId);
        gg.put("only_first_gap_carries_time", onlyFirstTimed);
        gg.put("final_ids", idsOf(gapss[gapss.length - 1]));
        gg.put("sets_equal_samples", sampleMismatch == 0);
        gg.put("sets_format", "t = gaps[0].relativeTime; sample_index = critical point whose sample (perturbed point) this set is: its gaps are, "
                + "in the same order, exactly getPhysicalGaps at that sample (sets_equal_samples); "
                + "gaps: [id, iEdge, toGapSet in Java HashSet iteration order]; "
                + "line = ProjectPanel/ProjectPanel5 stdout line: Double.toString(relativeTime) + ' ' + each Gap.toString() + ' '");
        gg.put("sets", sets);
        int repeats = 0, repeatMismatch = 0;

        // ---- NSE (processGapHistoryInformation, NEVER_SEE_EVADER)
        Map<String, Object> nse = Json.obj();
        try {
            Gap[][] g2 = Algorithm.getGaps(poly, makePath(wps));
            repeats++;
            if (!gapLines(g2, true).equals(lines0)) repeatMismatch++;
            geometry.pe.Algorithm.processGapHistoryInformation(g2, geometry.pe.Algorithm.SCENARIO.NEVER_SEE_EVADER, null);
            List<Object> st = new ArrayList<Object>();
            for (Gap[] gs : g2) {
                StringBuilder sb = new StringBuilder();
                for (Gap g : gs) { GapState s = g.getState(); sb.append(s == null ? '-' : (((NSEState) s).isClear() ? '0' : '1')); }
                st.add(sb.toString());
            }
            nse.put("states_format", "one string per gap set, one char per gap in set order: 0 = CLEAR, 1 = CONTAMINATED, - = null");
            nse.put("states", st);
            nse.put("lines_note", "ProjectPanel2/3 print each set as Gap.toString() + ' ' per gap, where toString inserts '| 0' (clear) or '| 1' (contaminated) after the iEdge");
        } catch (Throwable t) { nse.put("exception", exc(t)); }
        d.put("nse", nse);

        // ---- oracle + filter
        List<Object> oruns = new ArrayList<Object>();
        for (int variant = 0; variant < 2; variant++) {
            long[] ss = variant == 0 ? seeds : fixedSeeds;
            if (digitised && ss.length > 1) ss = new long[]{ss[0]};
            for (long seed : ss) {
                Gap[][] g3 = Algorithm.getGaps(poly, makePath(wps));
                repeats++;
                if (!gapLines(g3, true).equals(lines0)) repeatMismatch++;
                oruns.add(oracleRun(g3, nTargets, seed, variant == 0 ? "original" : "merge_fixed", oruns.isEmpty()));
            }
        }
        gg.put("repeat_calls", repeats);
        gg.put("repeat_calls_identical", repeatMismatch == 0);
        d.put("oracle_runs", oruns);
        return d;
    }

    static int[] idsOf(Gap[] gs) { int[] r = new int[gs.length]; for (int i = 0; i < gs.length; i++) r[i] = gs[i].getId(); return r; }

    /** Locate the critical point at which getGaps failed (see the class comment, item 3). */
    static Map<String, Object> locateFailure(Throwable t, PathCutIntersectPoint[] pc, List<Point2D> calls) {
        Map<String, Object> m = Json.obj();
        int getGapsLine = -1;
        String inner = null;
        for (StackTraceElement e : t.getStackTrace()) {
            if (e.getClassName().equals("geometry.Algorithm") && e.getMethodName().equals("getGaps")) { getGapsLine = e.getLineNumber(); break; }
            if (inner == null) inner = e.getClassName() + "." + e.getMethodName() + ":" + e.getLineNumber();
        }
        m.put("getGaps_line", getGapsLine);
        m.put("innermost_frame", inner);
        if (pc == null) return m;
        // calls[0] is the end point (from getAllCriticalPoints); then i == 0 and each GI/BT index
        int last = -1, from = 0;
        for (int c = 1; c < calls.size(); c++) {
            Point2D p = calls.get(c);
            for (int i = from; i < pc.length; i++) if (pc[i].getPoint().equals(p)) { last = i; from = i + 1; break; }
        }
        int next = -1;
        for (int i = last + 1; i < pc.length; i++) {
            Algorithm.CUT_TYPE ty = pc[i].getCutType();
            if (ty == Algorithm.CUT_TYPE.GENERAL_INFLECTION || ty == Algorithm.CUT_TYPE.BITANGENT) { next = i; break; }
        }
        m.put("last_timed_critical_index", last);
        String phase;
        int at;
        // Line ranges of Algorithm.getGaps in the original Algorithm.java:
        // 102-166 general-inflection branch, 168-290 bitangent branch, 302-310 group close.
        if (getGapsLine >= 102 && getGapsLine <= 166) { phase = "general-inflection branch of the event at critical_index"; at = next; }
        else if (getGapsLine >= 168 && getGapsLine <= 290) { phase = "bitangent branch of the event at critical_index"; at = next; }
        else if (getGapsLine >= 302 && getGapsLine <= 310) { phase = "group close (collapsePhysicalGaps over the samples of the group ending at critical_index)"; at = next == -1 ? pc.length - 1 : next - 1; }
        else { phase = "other"; at = -1; }
        m.put("phase", phase);
        m.put("critical_index", at);
        m.put("critical_type", at >= 0 ? typeName(pc[at].getCutType()) : null);
        return m;
    }

    static final Pattern BOUND = Pattern.compile("(Max|Min) flow in shadow (-?\\d+) is (-?\\d+)");

    static Map<String, Object> oracleRun(final Gap[][] gapss, int nTargets, long seed, String variant, boolean full) throws Exception {
        Map<String, Object> r = Json.obj();
        r.put("variant", variant);
        r.put("seed", seed);
        r.put("n_targets", nTargets);
        final Oracle o = variant.equals("original") ? new SingleTypeAgentOracle(nTargets) : new MergeFixedOracle(nTargets);
        try {
            seedMathRandom(seed);
            o.initialize(gapss);
        } catch (Throwable t) {
            r.put("initialize_exception", exc(t));
            return r;
        }
        // events
        List<Object> evs = new ArrayList<Object>();
        for (int i = 0; i < o.getNumberOfEvents(); i++) {
            SingleTypeAgentEvent e = (SingleTypeAgentEvent) o.getEvent(i);
            String s;
            try { s = e.toString(); } catch (Throwable t) { s = "EXC " + t; }
            List<Object> row = Json.arr(e.getEventRelativeTime(), e.getEventType() == null ? null : e.getEventType().name(),
                    e.getFromGap(), e.getToGap(), e.getMovedAgents(), e.getVisibleAgents());
            if (full) row.add(s);
            evs.add(row);
        }
        r.put("n_events", o.getNumberOfEvents());
        r.put("events_format", "[relativeTime, type, fromGap, toGap, movedAgents, visibleAgents" + (full ? ", toString()" : "")
                + "] in TreeMap (time) order" + (full ? "" : "; toString() is stored only in the first oracle run of the file"));
        r.put("events", evs);
        // ground truth of the final gap set, aligned with get_gaps.final_ids
        List<Object> truth = new ArrayList<Object>();
        for (Gap g : gapss[gapss.length - 1]) truth.add(g.getState() == null ? null : ((SingleTypeAgentState) g.getState()).getNumberOfAgents());
        r.put("truth_final_format", "number of targets the oracle put in each final gap, aligned with get_gaps.final_ids");
        r.put("truth_final", truth);
        // equations
        try {
            Equation[] eqs = SingleTypeAgentAlgorithm.deriveGapEvolvingEquations(gapss, o);
            List<Object> el = new ArrayList<Object>();
            for (Equation eq : eqs) el.add(eq.toString());
            r.put("equations", el);
        } catch (Throwable t) { r.put("equations_exception", exc(t)); }
        // bipartite graph (first build) + bounds under several identity-hash orders
        Map<String, Object> outcomes = new LinkedHashMap<String, Object>();
        Map<String, List<Object>> outcomeBurn = new LinkedHashMap<String, List<Object>>();
        for (int bi = 0; bi < burnins.length; bi++) {
            burn(burnins[bi]);
            final BipartiteGraph[] bg = new BipartiteGraph[1];
            String bgText;
            try { bgText = capture(new Body() { public void run() { bg[0] = SingleTypeAgentAlgorithm.deriveShadowInfoState(gapss, o); } }); }
            catch (Throwable t) { r.put("bg_exception", exc(t)); break; }
            if (bi == 0) {
                if (full) r.put("bg_text_example", bgText);
                r.put("bg", bgStructure(bg[0]));
            }
            String bText;
            try { bText = capture(new Body() { public void run() { SingleTypeAgentAlgorithm.deriveShadowBounds(bg[0]); } }); }
            catch (Throwable t) { r.put("bounds_exception", exc(t)); break; }
            Map<Integer, int[]> b = new TreeMap<Integer, int[]>();
            Matcher m = BOUND.matcher(bText);
            while (m.find()) {
                int id = Integer.parseInt(m.group(2));
                int v = Integer.parseInt(m.group(3));
                int[] mm = b.get(id);
                if (mm == null) { mm = new int[]{Integer.MIN_VALUE, Integer.MIN_VALUE}; b.put(id, mm); }
                if (m.group(1).equals("Min")) mm[0] = v; else mm[1] = v;
            }
            Map<String, Object> bj = Json.obj();
            for (Map.Entry<Integer, int[]> e : b.entrySet()) bj.put(String.valueOf(e.getKey()), e.getValue());
            String key = Json.write(bj);
            if (!outcomes.containsKey(key)) { outcomes.put(key, bj); outcomeBurn.put(key, new ArrayList<Object>()); }
            outcomeBurn.get(key).add(burnins[bi]);
        }
        List<Object> ol = new ArrayList<Object>();
        for (String k : outcomes.keySet()) {
            Map<String, Object> x = Json.obj();
            x.put("burnins", outcomeBurn.get(k));
            x.put("bounds", outcomes.get(k));
            ol.add(x);
        }
        r.put("java_bounds_format", "{shadow id: [min, max]} parsed from the 'Min/Max flow in shadow i is v' lines of deriveShadowBounds; "
                + "one entry per distinct outcome over the identity-hash burn-ins");
        r.put("java_bounds_outcomes", ol);
        r.put("java_bounds_order_dependent", ol.size() > 1);
        return r;
    }

    static Map<String, Object> bgStructure(BipartiteGraph bg) {
        Map<String, Object> m = Json.obj();
        List<Object> left = new ArrayList<Object>();
        int edges = 0;
        for (Vertex v : bg.leftVertexMap.values()) {
            List<Integer> to = new ArrayList<Integer>();
            for (Edge e : v.veMap.values()) to.add(e.toVertex.id);
            java.util.Collections.sort(to);
            edges += to.size();
            left.add(Json.arr(v.id, v.minWeight, v.maxWeight, new ArrayList<Object>(to)));
        }
        List<Object> right = new ArrayList<Object>();
        for (Vertex v : bg.rightVertexMap.values()) right.add(Json.arr(v.id, v.minWeight, v.maxWeight));
        m.put("format", "left: [id, minWeight, maxWeight, sorted ids of right vertices it has edges to]; right: [id, minWeight, maxWeight] (-1 = unknown, i.e. alive)");
        m.put("n_left", bg.leftVertexMap.size());
        m.put("n_right", bg.rightVertexMap.size());
        m.put("n_edges", edges);
        m.put("left", left);
        m.put("right", right);
        return m;
    }
}
