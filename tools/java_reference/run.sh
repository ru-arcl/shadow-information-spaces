#!/usr/bin/env bash
# Build the original Java sources and the harness, then write tests/fixtures/java/*.json.
#
#   ORIG_SRC=/path/to/shadow-information-space/source \
#   LOG4J_JAR=/path/to/log4j-1.2.12.jar \
#   JAVA_HOME=/path/to/jdk \
#   tools/java_reference/run.sh [extra GoldenDump options]
#
# See README.md in this folder for details.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
: "${ORIG_SRC:?set ORIG_SRC to the 'source' folder of the original code (contains geometry/, polygons/, config/)}"
: "${LOG4J_JAR:?set LOG4J_JAR to log4j-1.2.12.jar}"
JAVA="${JAVA_HOME:+$JAVA_HOME/bin/}java"
JAVAC="${JAVA_HOME:+$JAVA_HOME/bin/}javac"
BUILD="${BUILD:-$HERE/build}"
OUT="${OUT:-$REPO/tests/fixtures/java}"

mkdir -p "$BUILD/orig" "$BUILD/harness"
# 1. the original code, unmodified (the pe/ui applet classes compile headless too)
find "$ORIG_SRC/geometry" -name '*.java' > "$BUILD/orig-sources.txt"
if ! "$JAVAC" -nowarn -encoding ISO-8859-1 -cp "$LOG4J_JAR" -d "$BUILD/orig" @"$BUILD/orig-sources.txt" > "$BUILD/javac-orig.log" 2>&1; then
  cat "$BUILD/javac-orig.log" >&2; exit 1
fi
# 2. the harness
"$JAVAC" -encoding UTF-8 -cp "$BUILD/orig:$LOG4J_JAR" -d "$BUILD/harness" "$HERE"/src/*.java
# 3. run
"$JAVA" -Djava.awt.headless=true -XX:-OmitStackTraceInFastThrow --add-opens java.base/java.lang=ALL-UNNAMED -Xss16m \
  -cp "$BUILD/harness:$BUILD/orig:$LOG4J_JAR" \
  GoldenDump "$ORIG_SRC" "$HERE/paths.tsv" "$OUT" "$@"
