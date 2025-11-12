package main;

import java.util.Arrays;

public class DPFolding {
    private static final int MIN_LOOP = 3;
    private static final double INF = 1e9;
    private static final double EPS = 1e-6;

    private final EnergyModel model;

    public DPFolding(EnergyModel model) {
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

        for (int i = 0; i < n; i++) {
            W[i][i] = 0.0;
        }

        for (int len = 2; len <= n; len++) {
            for (int i = 0; i + len - 1 < n; i++) {
                int j = i + len - 1;
                if (j - i < MIN_LOOP) {
                    W[i][j] = 0.0;
                    continue;
                }

                // V[i][j]
                if (model.canPair(seq.charAt(i), seq.charAt(j))) {
                    double bestV = INF;

                    // hairpin
                    int loopSize = j - i - 1;
                    if (loopSize >= MIN_LOOP) {
                        double hp = model.hairpinEnergy(loopSize);
                        double closePair = model.pairEnergy(seq.charAt(i), seq.charAt(j));
                        double hairpinTotal = hp + closePair;
                        bestV = Math.min(bestV, hairpinTotal);
                    }

                    // stacked
                    if (i + 1 < j - 1 && V[i + 1][j - 1] < INF
                        && model.canPair(seq.charAt(i + 1), seq.charAt(j - 1))) {
                    	double stack = V[i + 1][j - 1]
                    	        + model.stackEnergy(seq.charAt(i), seq.charAt(i + 1), seq.charAt(j - 1), seq.charAt(j));
                        bestV = Math.min(bestV, stack);
                    }

                    // internal loops
                    double bestInternal = INF;
                    for (int i2 = i + 1; i2 <= j - MIN_LOOP - 1; i2++) {
                        for (int j2 = Math.max(i2 + MIN_LOOP + 1, i + 2); j2 < j; j2++) {
                            if (!model.canPair(seq.charAt(i2), seq.charAt(j2))) continue;
                            if (V[i2][j2] >= INF) continue;
                            double eL = model.internalLoopPenalty(i - j);
                            bestInternal = Math.min(bestInternal, V[i2][j2] + eL);
                        }
                    }
                    bestV = Math.min(bestV, bestInternal);

                    // multi-loop closed
                    if (i + 1 <= j - 1 && WM[i + 1][j - 1] < INF) {
                        bestV = Math.min(bestV, WM[i + 1][j - 1] + model.multiInit());
                    }

                    V[i][j] = bestV;
                }

                // WM[i][j]
                double bestWM = INF;
                // i unpaired inside multi-loop
                if (i + 1 <= j && WM[i + 1][j] < INF) {
                    bestWM = Math.min(bestWM, WM[i + 1][j] + model.multiUnpaired());
                }
                // j unpaired inside multi-loop
                if (i <= j - 1 && WM[i][j - 1] < INF) {
                    bestWM = Math.min(bestWM, WM[i][j - 1] + model.multiUnpaired());
                }
                // single branch: closing pair (i,j) forms a branch
                if (V[i][j] < INF) {
                    bestWM = Math.min(bestWM, V[i][j] + model.multiBranch());
                }
                // split multi-loop into two WM regions
                for (int k = i + 1; k < j; k++) {
                    if (WM[i][k] < INF && WM[k + 1][j] < INF) {
                        bestWM = Math.min(bestWM, WM[i][k] + WM[k + 1][j]);
                    }
                }
                
                WM[i][j] = bestWM;

                double bestW = W[i][j - 1];
                for (int k = i; k <= j - MIN_LOOP - 1; k++) {
                    if (V[k][j] < INF) {
                        double left = (k - 1 >= i) ? W[i][k - 1] : 0.0;
                        bestW = Math.min(bestW, left + V[k][j]);
                    }
                }
                W[i][j] = bestW;
            }
        }

        // Traceback
        int[] pairTo = new int[n];
        Arrays.fill(pairTo, -1);
        tracebackW(0, n - 1, W, V, WM, seq, pairTo);

        String dotBracket = toDotBracket(pairTo);
        double mfe = W[0][n - 1];
        return new Result(dotBracket, mfe, pairTo);
    }

    // Traceback for W
    private void tracebackW(int i, int j, double[][] W, double[][] V, double[][] WM,
                            String seq, int[] pairTo) {
        if (i > j) return;
        if (i == j) return;
        double cur = W[i][j];
        
        if (Math.abs(cur - W[i][j - 1]) < EPS) {
            tracebackW(i, j - 1, W, V, WM, seq, pairTo);
            return;
        }
        for (int k = i; k <= j - MIN_LOOP - 1; k++) {
            if (V[k][j] >= INF) continue;
            double left = (k - 1 >= i) ? W[i][k - 1] : 0.0;
            if (Math.abs(cur - (left + V[k][j])) < 1e-3) {
                tracebackV(k, j, V, WM, seq, pairTo);
                if (k - 1 >= i) tracebackW(i, k - 1, W, V, WM, seq, pairTo);
                return;
            }
        }
        tracebackW(i, j - 1, W, V, WM, seq, pairTo);
    }

    private void tracebackV(int i, int j, double[][] V, double[][] WM,
                            String seq, int[] pairTo) {
        if (i >= j) return;
        if (V[i][j] >= INF) return;
        pairTo[i] = j;
        pairTo[j] = i;

        double cur = V[i][j];

        int loopSize = j - i - 1;
        if (loopSize >= MIN_LOOP) {
            double hairpin = model.hairpinEnergy(loopSize);
            if (Math.abs(cur - hairpin) < 1e-3) {
                return;
            }
        }

        if (i + 1 < j - 1 && V[i + 1][j - 1] < INF) {
            double stack = V[i + 1][j - 1]
                    + model.stackEnergy(seq.charAt(i), seq.charAt(i + 1),
                                        seq.charAt(j - 1), seq.charAt(j));
            if (Math.abs(cur - stack) < 1e-3) {
                tracebackV(i + 1, j - 1, V, WM, seq, pairTo);
                return;
            }
        }

        // internal loop
        for (int i2 = i + 1; i2 <= j - MIN_LOOP - 1; i2++) {
            for (int j2 = Math.max(i2 + MIN_LOOP + 1, i + 2); j2 < j; j2++) {
                if (!model.canPair(seq.charAt(i2), seq.charAt(j2))) continue;
                if (V[i2][j2] >= INF) continue;
                double eL = model.internalLoopPenalty(i2-j2);
                if (Math.abs(cur - (V[i2][j2] + eL)) < 1e-3) {
                    // reconstruct inner pair
                    tracebackV(i2, j2, V, WM, seq, pairTo);
                    return;
                }
            }
        }

        // multi-loop
        if (i + 1 <= j - 1 && WM[i + 1][j - 1] < INF) {
            double multi = WM[i + 1][j - 1] + model.multiInit();
            if (Math.abs(cur - multi) < 1e-3) {
                tracebackWM(i + 1, j - 1, WM, V, seq, pairTo);
                return;
            }
        }
    }

    private void tracebackWM(int i, int j, double[][] WM, double[][] V,
                             String seq, int[] pairTo) {
        if (i > j) return;
        if (WM[i][j] >= INF) return;
        double cur = WM[i][j];

        if (V[i][j] < INF) {
            double cand = V[i][j] + model.multiBranch();
            if (Math.abs(cur - cand) < 1e-3) {
                tracebackV(i, j, V, WM, seq, pairTo);
                return;
            }
        }

        // i unpaired inside multi-loop
        if (i + 1 <= j && WM[i + 1][j] < INF) {
            double cand = WM[i + 1][j] + model.multiUnpaired();
            if (Math.abs(cur - cand) < 1e-3) {
                tracebackWM(i + 1, j, WM, V, seq, pairTo);
                return;
            }
        }

        // j unpaired inside multi-loop
        if (i <= j - 1 && WM[i][j - 1] < INF) {
            double cand = WM[i][j - 1] + model.multiUnpaired();
            if (Math.abs(cur - cand) < 1e-3) {
                tracebackWM(i, j - 1, WM, V, seq, pairTo);
                return;
            }
        }

        // split
        for (int k = i + 1; k < j; k++) {
            if (WM[i][k] < INF && WM[k + 1][j] < INF) {
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
