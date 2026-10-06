import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Minimal JSON writer, so the harness needs nothing beyond the JDK and log4j.
 *
 * Values may be Map (written in insertion order), List, int[], double[],
 * Integer, Long, Double, Boolean, String or null.
 *
 * Doubles are written with Java's Double.toString. That text is valid JSON,
 * e.g. "1.0E-4", and Python's float() parses it back to the same double. Every
 * double is checked with Double.parseDouble(Double.toString(d)) == d. A value
 * that fails this check, and any NaN or infinity, is written as a string:
 * "NaN", "Infinity", "-Infinity", or "hex:<Double.toHexString>". The counter
 * nonRoundTrip records how many such strings were written; it is expected
 * to stay 0.
 *
 * Layout, to keep files small but diffable: maps up to depth 3 that hold
 * containers, lists of maps up to depth 2, and lists of containers at
 * depth 0 or 1 are written one element per line with a one-space indent.
 * Everything else goes on one line.
 */
public final class Json {
    public static int nonRoundTrip = 0;

    private Json() {}

    public static Map<String, Object> obj() { return new LinkedHashMap<String, Object>(); }

    public static List<Object> arr(Object... xs) { return new ArrayList<Object>(Arrays.asList(xs)); }

    public static String write(Object v) {
        StringBuilder sb = new StringBuilder();
        write(sb, v, 0);
        sb.append('\n');
        return sb.toString();
    }

    private static boolean isContainer(Object v) {
        return v instanceof Map || v instanceof List || v instanceof int[] || v instanceof long[] || v instanceof double[];
    }

    private static boolean hasContainerChild(Object v) {
        if (v instanceof Map) {
            for (Object c : ((Map<?, ?>) v).values()) if (isContainer(c)) return true;
        } else if (v instanceof List) {
            for (Object c : (List<?>) v) if (isContainer(c)) return true;
        }
        return false;
    }

    private static boolean hasMapChild(List<?> l) {
        for (Object c : l) if (c instanceof Map) return true;
        return false;
    }

    private static void nl(StringBuilder sb, int depth) {
        sb.append('\n');
        for (int i = 0; i < depth; i++) sb.append(' ');
    }

    @SuppressWarnings("unchecked")
    private static void write(StringBuilder sb, Object v, int depth) {
        if (v == null) {
            sb.append("null");
        } else if (v instanceof String) {
            str(sb, (String) v);
        } else if (v instanceof Boolean || v instanceof Integer || v instanceof Long) {
            sb.append(v.toString());
        } else if (v instanceof Double || v instanceof Float) {
            dbl(sb, ((Number) v).doubleValue());
        } else if (v instanceof int[]) {
            int[] a = (int[]) v;
            sb.append('[');
            for (int i = 0; i < a.length; i++) { if (i > 0) sb.append(','); sb.append(a[i]); }
            sb.append(']');
        } else if (v instanceof long[]) {
            long[] a = (long[]) v;
            sb.append('[');
            for (int i = 0; i < a.length; i++) { if (i > 0) sb.append(','); sb.append(a[i]); }
            sb.append(']');
        } else if (v instanceof double[]) {
            double[] a = (double[]) v;
            sb.append('[');
            for (int i = 0; i < a.length; i++) { if (i > 0) sb.append(','); dbl(sb, a[i]); }
            sb.append(']');
        } else if (v instanceof Map) {
            Map<String, Object> m = (Map<String, Object>) v;
            boolean multi = depth <= 3 && hasContainerChild(v);
            sb.append('{');
            boolean first = true;
            for (Map.Entry<String, Object> e : m.entrySet()) {
                if (!first) sb.append(',');
                first = false;
                if (multi) nl(sb, depth + 1);
                str(sb, e.getKey());
                sb.append(':');
                write(sb, e.getValue(), depth + 1);
            }
            if (multi && !m.isEmpty()) nl(sb, depth);
            sb.append('}');
        } else if (v instanceof List) {
            List<Object> l = (List<Object>) v;
            boolean multi = hasContainerChild(v) && (depth <= 1 || (depth <= 2 && hasMapChild(l)));
            sb.append('[');
            for (int i = 0; i < l.size(); i++) {
                if (i > 0) sb.append(',');
                if (multi) nl(sb, depth + 1);
                write(sb, l.get(i), depth + 1);
            }
            if (multi && !l.isEmpty()) nl(sb, depth);
            sb.append(']');
        } else {
            throw new IllegalArgumentException("cannot serialise " + v.getClass());
        }
    }

    private static void dbl(StringBuilder sb, double d) {
        if (Double.isNaN(d)) { sb.append("\"NaN\""); return; }
        if (Double.isInfinite(d)) { sb.append(d > 0 ? "\"Infinity\"" : "\"-Infinity\""); return; }
        String s = Double.toString(d);
        if (Double.parseDouble(s) != d || (d == 0 && (1 / d) != (1 / Double.parseDouble(s)))) {
            nonRoundTrip++;
            sb.append("\"hex:").append(Double.toHexString(d)).append('"');
            return;
        }
        sb.append(s);
    }

    private static void str(StringBuilder sb, String s) {
        sb.append('"');
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            switch (c) {
                case '"': sb.append("\\\""); break;
                case '\\': sb.append("\\\\"); break;
                case '\n': sb.append("\\n"); break;
                case '\r': sb.append("\\r"); break;
                case '\t': sb.append("\\t"); break;
                default:
                    if (c < 0x20) sb.append(String.format("\\u%04x", (int) c));
                    else sb.append(c);
            }
        }
        sb.append('"');
    }
}
