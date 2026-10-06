import geometry.gap.Gap;
import geometry.pe.is.Event;
import geometry.pe.is.Oracle;
import geometry.pe.is.singleTypeAgent.SingleTypeAgentEvent;
import geometry.pe.is.singleTypeAgent.SingleTypeAgentState;

import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * A copy of geometry.pe.is.singleTypeAgent.SingleTypeAgentOracle with one
 * intended change: the merge bug (B7 in docs/notes/original_java.md) is fixed.
 *
 * The original code handles the second merging parent like this:
 *
 *   toGap.setState(new SingleTypeAgentState(
 *       ((SingleTypeAgentState) pGaps[j].getState()).getNumberOfAgents() + t));
 *
 * Here t is the count of that same parent, so the target gets
 * count2 + count2, and the first parent's count1 is lost. This copy uses
 * count1 + count2 instead (see the FIX line below), so targets are conserved
 * and the observations stay consistent.
 *
 * Everything else is copied verbatim: the Math.random() draws (the harness
 * seeds them, as for the original), distributeOneBatch, event types, the
 * TreeMap keyed by relative time, and the unused random visibility events.
 * Two harmless edits: Double.valueOf replaces the deprecated new Double(...),
 * and Integer.valueOf replaces new Integer(...).
 *
 * The harness runs this oracle only as a separate variant ("merge_fixed").
 * The golden "original" runs use the unmodified class from the original
 * sources.
 */
public class MergeFixedOracle extends Oracle {

    private int totalAgents;

    private Map<Double, SingleTypeAgentEvent> eventMap = new TreeMap<Double, SingleTypeAgentEvent>();

    private SingleTypeAgentEvent stae[] = null;

    public MergeFixedOracle(int totalAgents) {
        super();
        this.totalAgents = totalAgents;
    }

    public int getTotalAgentNumber() { return totalAgents; }

    @Override
    public int getNumberOfEvents() { return stae.length; }

    @Override
    public Event getEvent(int i) { return stae[i]; }

    @Override
    public void initialize(Gap[][] gapss) {
        eventMap.clear();
        distributeAgents(gapss);
        stae = eventMap.values().toArray(new SingleTypeAgentEvent[0]);
    }

    @Override
    public Event getEventByTime(double t) { return eventMap.get(Double.valueOf(t)); }

    private void distributeAgents(Gap[][] gapss) {
        int previousVisibleAgents = 0;
        for (int i = 0; i < gapss.length; i++) {
            Gap[] gaps = gapss[i];
            SingleTypeAgentEvent e = new SingleTypeAgentEvent();
            if (i == 0) {
                previousVisibleAgents = (int) (Math.random() * 0.5 * totalAgents);
                e.setVisibleAgents(previousVisibleAgents);
                e.setEventType(EVENT_TYPE.INIT);
                int[] b = distributeOneBatch(totalAgents - previousVisibleAgents, gaps.length);
                int toGaps[] = new int[gaps.length];
                e.setToGap(toGaps);
                for (int j = 0; j < gaps.length; j++) {
                    gaps[j].setState(new SingleTypeAgentState(b[j]));
                    toGaps[j] = gaps[j].getId();
                }
            } else {
                Map<Integer, Gap> gm = Gap.getGapMap(gaps);
                Gap[] pGaps = gapss[i - 1];
                Map<Integer, Gap> pgm = Gap.getGapMap(pGaps);
                if (gaps.length > pGaps.length) {
                    for (int j = 0; j < pGaps.length; j++) {
                        Integer[] toGaps = pGaps[j].getToGapSet().toArray(new Integer[0]);
                        if (toGaps.length == 2) {
                            e.setEventType(EVENT_TYPE.SPLIT);
                            int t = ((SingleTypeAgentState) pGaps[j].getState()).getNumberOfAgents();
                            int[] b = distributeOneBatch(t, 2);
                            gm.get(toGaps[0]).setState(new SingleTypeAgentState(b[0]));
                            gm.get(toGaps[1]).setState(new SingleTypeAgentState(b[1]));
                            e.setFromGap(new int[]{pGaps[j].getId()});
                            e.setToGap(new int[]{toGaps[0].intValue(), toGaps[1].intValue()});
                        } else {
                            gm.get(pGaps[j].getId()).setState(pGaps[j].getState());
                        }
                        gm.remove(pGaps[j].getId());
                    }
                    if (gm.size() == 1) {
                        e.setEventType(EVENT_TYPE.APPEAR);
                        int gn = gm.keySet().toArray(new Integer[0])[0].intValue();
                        int[] b = distributeOneBatch(previousVisibleAgents, 2);
                        gm.get(gn).setState(new SingleTypeAgentState(b[0]));
                        previousVisibleAgents = b[1];
                        e.setVisibleAgents(previousVisibleAgents);
                        e.setToGap(new int[]{gm.get(gn).getId()});
                        e.setMovedAgents(b[0]);
                    } else {
                        e.setVisibleAgents(previousVisibleAgents);
                    }
                } else {
                    for (int j = 0; j < pGaps.length; j++) {
                        Integer[] toGaps = pGaps[j].getToGapSet().toArray(new Integer[0]);
                        if (toGaps.length == 1) {
                            int t = ((SingleTypeAgentState) pGaps[j].getState()).getNumberOfAgents();
                            Gap toGap = gm.get(toGaps[0]);
                            if (toGap.getState() == null) {
                                e.setEventType(EVENT_TYPE.MERGE);
                                e.setVisibleAgents(previousVisibleAgents);
                                toGap.setState(new SingleTypeAgentState(t));
                                e.setToGap(new int[]{toGap.getId()});
                                e.setFromGap(new int[]{pGaps[j].getId(), 0});
                            } else {
                                // FIX (the only semantic change): count1 + count2.
                                toGap.setState(new SingleTypeAgentState(
                                        ((SingleTypeAgentState) toGap.getState()).getNumberOfAgents() + t));
                                e.getFromGap()[1] = pGaps[j].getId();
                            }
                        } else {
                            Gap toGap = gm.get(pGaps[j].getId());
                            if (toGap == null) {
                                e.setEventType(EVENT_TYPE.DISAPPEAR);
                                int ng = ((SingleTypeAgentState) pGaps[j].getState()).getNumberOfAgents();
                                previousVisibleAgents += ng;
                                e.setVisibleAgents(previousVisibleAgents);
                                e.setFromGap(new int[]{pGaps[j].getId()});
                                e.setMovedAgents(ng);
                            } else {
                                toGap.setState(pGaps[j].getState());
                            }
                        }
                        pgm.remove(pGaps[j].getId());
                    }
                }
            }
            e.setEventRelativeTime(gapss[i][0].getRelativeTime());
            eventMap.put(gapss[i][0].getRelativeTime(), e);
        }
    }

    private int[] distributeOneBatch(int t, int g) {
        if (g == 1) return new int[]{t};
        if (t == 0) return new int[g];
        Set<Integer> spSet = new TreeSet<Integer>();
        while (spSet.size() < g - 1) {
            spSet.add(Integer.valueOf((int) (Math.random() * t)));
        }
        spSet.add(t);
        Integer[] retI = spSet.toArray(new Integer[0]);
        int[] ret = new int[g];
        ret[0] = retI[0].intValue();
        for (int i = 1; i < ret.length; i++) ret[i] = retI[i] - retI[i - 1];
        return ret;
    }
}
