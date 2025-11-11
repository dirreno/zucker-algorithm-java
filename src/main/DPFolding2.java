package main;

import java.util.Arrays;

public class DPFolding2 {
    private static final int MIN_LOOP = 3;
    private static final double INF = 1e9;
    private static final double EPS = 1e-6;

    private final EnergyModel model;

    public DPFolding2(EnergyModel model) {
        this.model = model;
    }

    public Result fold(String seqIn) {
        String seq = seqIn.toUpperCase().replaceAll("\\s+", "");
        final int n = seq.length();
        if (n == 0) return new Result("", 0.0, new int[0]);

        double[][] W = new double[n][n];
        double[][] V = new double[n][n];
        double[][] WM = new double[n][n];
        

        for (int i = 0; i < n; i++) {
            Arrays.fill(W[i], INF);
            Arrays.fill(V[i], INF);
            Arrays.fill(WM[i], INF);
        }

        // length-1 intervals: empty / single base
        for (int i = 0; i < n; i++) {
            W[i][i] = 0.0;
            // V and WM remain INF for singletons
        }

        // fill by increasing length
        for (int len = 2; len <= n; len++) {
            for (int i = 0; i + len - 1 < n; i++) {
                int j = i + len - 1;
                if (j - i < MIN_LOOP) {
                    // small intervals: allow W to be zero or carry previous
                    W[i][j] = 0.0;
                    continue;
                }

                // ---- V[i][j] ----
                // only compute V if i and j can pair
                if (model.canPair(seq.charAt(i), seq.charAt(j))) {
                    double bestV = INF;

                    // hairpin: loop size = j - i - 1
                    int loopSize = j - i - 1;
                    if (loopSize >= MIN_LOOP) {
                        // use hairpinEnergy(int)
                        double hairpin = model.hairpinEnergy(loopSize);
                        bestV = Math.min(bestV, hairpin);
                    }

                    // stacked pair: (i+1, j-1) exists and paired
                    if (i + 1 < j - 1 && V[i + 1][j - 1] < INF/2
                        && model.canPair(seq.charAt(i + 1), seq.charAt(j - 1))) {
                        double stack = V[i + 1][j - 1]
                                + model.stackEnergy(seq.charAt(i), seq.charAt(i + 1), seq.charAt(j - 1), seq.charAt(j));
                        bestV = Math.min(bestV, stack);
                    }

                    // internal/bulge loops: find an inner closing pair (i2,j2)
                    double bestInternal = INF;
                    for (int i2 = i + 1; i2 <= j - MIN_LOOP - 1; i2++) {
                        for (int j2 = Math.max(i2 + MIN_LOOP + 1, i + 2); j2 < j; j2++) {
                            if (!model.canPair(seq.charAt(i2), seq.charAt(j2))) continue;
                            if (V[i2][j2] >= INF/2) continue;
                            double eL = internalLoopEnergyApprox(i, j, i2, j2);
                            bestInternal = Math.min(bestInternal, V[i2][j2] + eL);
                        }
                    }
                    bestV = Math.min(bestV, bestInternal);

                    // multi-loop closed by (i,j): WM[i+1][j-1] + a
                    if (i + 1 <= j - 1 && WM[i + 1][j - 1] < INF/2) {
                        bestV = Math.min(bestV, WM[i + 1][j - 1] + model.multiInit());
                    }

                    V[i][j] = bestV;
                } // end V

                // ---- WM[i][j] ----
                double bestWM = INF;
                // i unpaired inside multi-loop
                if (i + 1 <= j && WM[i + 1][j] < INF/2) {
                    bestWM = Math.min(bestWM, WM[i + 1][j] + model.multiUnpaired());
                }
                // j unpaired inside multi-loop
                if (i <= j - 1 && WM[i][j - 1] < INF/2) {
                    bestWM = Math.min(bestWM, WM[i][j - 1] + model.multiUnpaired());
                }
                // single branch: closing pair (i,j) forms a branch
                if (V[i][j] < INF/2) {
                    bestWM = Math.min(bestWM, V[i][j] + model.multiBranch());
                }
                // split multi-loop into two WM regions
                for (int k = i + 1; k < j; k++) {
                    if (WM[i][k] < INF/2 && WM[k + 1][j] < INF/2) {
                        bestWM = Math.min(bestWM, WM[i][k] + WM[k + 1][j]);
                    }
                }
                // base case: allow a WM single branch starting at an inner V
                // if nothing found, keep INF
                WM[i][j] = bestWM;

                // ---- W[i][j] ----
                // j unpaired
                double bestW = W[i][j - 1];
                // split: W[i][k-1] + V[k][j]
                for (int k = i; k <= j - MIN_LOOP - 1; k++) {
                    if (V[k][j] < INF/2) {
                        double left = (k - 1 >= i) ? W[i][k - 1] : 0.0;
                        bestW = Math.min(bestW, left + V[k][j]);
                    }
                }
                W[i][j] = bestW;
            } // end i loop
        } // end len loop

        // traceback using W,V,WM to build pairTo
        int[] pairTo = new int[n];
        Arrays.fill(pairTo, -1);
        tracebackW(0, n - 1, W, V, WM, seq, pairTo);

        String dotBracket = toDotBracket(pairTo);
        double mfe = W[0][n - 1];
        return new Result(dotBracket, mfe, pairTo);
    }

    /* ---------------------- helper functions ------------------------- */

    // simple approximate internal loop energy: use bulgePenalty if one side 0,
    // otherwise combine hairpin penalties as a cheap approximation.
    private double internalLoopEnergyApprox(int i, int j, int i2, int j2) {
        // outer pair (i,j) and inner pair (i2,j2) with i<i2<j2<j
        int a = i2 - i - 1;    // unpaired on left
        int b = j - j2 - 1;    // unpaired on right
        int total = a + b;
        if (total < 0) return INF;
        if (a == 0 || b == 0) {
            // bulge
            int bulgeSize = total;
            if (bulgeSize < 1) return INF;
            if (bulgeSize > 30) bulgeSize = 30;
            return model.InternalLoopPenalty(bulgeSize);
        } else {
            // internal loop: approximate by average hairpin penalties + mild penalty
            int ta = Math.min(a, 30);
            int tb = Math.min(b, 30);
            double ha = model.hairpinEnergy(Math.max(3, ta));
            double hb = model.hairpinEnergy(Math.max(3, tb));
            return 0.5 * (ha + hb) + 0.2 * total;
        }
    }

    // Traceback for W
    private void tracebackW(int i, int j, double[][] W, double[][] V, double[][] WM,
                            String seq, int[] pairTo) {
        if (i > j) return;
        if (i == j) return;
        double cur = W[i][j];
        // j unpaired?
        if (Math.abs(cur - W[i][j - 1]) < EPS) {
            tracebackW(i, j - 1, W, V, WM, seq, pairTo);
            return;
        }
        // otherwise find k where cur == left + V[k][j]
        for (int k = i; k <= j - MIN_LOOP - 1; k++) {
            if (V[k][j] >= INF/2) continue;
            double left = (k - 1 >= i) ? W[i][k - 1] : 0.0;
            if (Math.abs(cur - (left + V[k][j])) < 1e-3) {
                // reconstruct the V[k][j] region and left part
                tracebackV(k, j, V, WM, seq, pairTo);
                if (k - 1 >= i) tracebackW(i, k - 1, W, V, WM, seq, pairTo);
                return;
            }
        }
        // fallback: nothing matched, try shrinking the interval
        tracebackW(i, j - 1, W, V, WM, seq, pairTo);
    }

    // Traceback for V (closed by (i,j))
    private void tracebackV(int i, int j, double[][] V, double[][] WM,
                            String seq, int[] pairTo) {
        if (i >= j) return;
        if (V[i][j] >= INF/2) return;
        pairTo[i] = j;
        pairTo[j] = i;

        double cur = V[i][j];

        // hairpin?
        int loopSize = j - i - 1;
        if (loopSize >= MIN_LOOP) {
            double hairpin = model.hairpinEnergy(loopSize);
            if (Math.abs(cur - hairpin) < 1e-3) {
                // hairpin closed by (i,j) — nothing inside
                return;
            }
        }

        // stacked?
        if (i + 1 < j - 1 && V[i + 1][j - 1] < INF/2) {
            double stack = V[i + 1][j - 1]
                    + model.stackEnergy(seq.charAt(i), seq.charAt(i + 1),
                                        seq.charAt(j - 1), seq.charAt(j));
            if (Math.abs(cur - stack) < 1e-3) {
                // stacked pair: continue inside
                tracebackV(i + 1, j - 1, V, WM, seq, pairTo);
                return;
            }
        }

        // internal loop case: find inner pair (i2,j2)
        for (int i2 = i + 1; i2 <= j - MIN_LOOP - 1; i2++) {
            for (int j2 = Math.max(i2 + MIN_LOOP + 1, i + 2); j2 < j; j2++) {
                if (!model.canPair(seq.charAt(i2), seq.charAt(j2))) continue;
                if (V[i2][j2] >= INF/2) continue;
                double eL = internalLoopEnergyApprox(i, j, i2, j2);
                if (Math.abs(cur - (V[i2][j2] + eL)) < 1e-3) {
                    // reconstruct inner pair
                    tracebackV(i2, j2, V, WM, seq, pairTo);
                    return;
                }
            }
        }

        // multi-loop closure: V[i][j] == WM[i+1][j-1] + a
        if (i + 1 <= j - 1 && WM[i + 1][j - 1] < INF/2) {
            double multi = WM[i + 1][j - 1] + model.multiInit();
            if (Math.abs(cur - multi) < 1e-3) {
                tracebackWM(i + 1, j - 1, WM, V, seq, pairTo);
                return;
            }
        }
    }

    // Traceback for WM (multi-loop interior)
    private void tracebackWM(int i, int j, double[][] WM, double[][] V,
                             String seq, int[] pairTo) {
        if (i > j) return;
        if (WM[i][j] >= INF/2) return;
        double cur = WM[i][j];

        // V[i][j] + b ?
        if (V[i][j] < INF/2) {
            double cand = V[i][j] + model.multiBranch();
            if (Math.abs(cur - cand) < 1e-3) {
                tracebackV(i, j, V, WM, seq, pairTo);
                return;
            }
        }

        // i unpaired inside multi-loop
        if (i + 1 <= j && WM[i + 1][j] < INF/2) {
            double cand = WM[i + 1][j] + model.multiUnpaired();
            if (Math.abs(cur - cand) < 1e-3) {
                tracebackWM(i + 1, j, WM, V, seq, pairTo);
                return;
            }
        }

        // j unpaired inside multi-loop
        if (i <= j - 1 && WM[i][j - 1] < INF/2) {
            double cand = WM[i][j - 1] + model.multiUnpaired();
            if (Math.abs(cur - cand) < 1e-3) {
                tracebackWM(i, j - 1, WM, V, seq, pairTo);
                return;
            }
        }

        // split
        for (int k = i + 1; k < j; k++) {
            if (WM[i][k] < INF/2 && WM[k + 1][j] < INF/2) {
                if (Math.abs(cur - (WM[i][k] + WM[k + 1][j])) < 1e-3) {
                    tracebackWM(i, k, WM, V, seq, pairTo);
                    tracebackWM(k + 1, j, WM, V, seq, pairTo);
                    return;
                }
            }
        }
    }

    private String toDotBracket(int[] pairTo) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < pairTo.length; i++) {
            int p = pairTo[i];
            sb.append(p == -1 ? '.' : (p > i ? '(' : ')'));
        }
        return sb.toString();
    }
}
