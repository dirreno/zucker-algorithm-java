package main;

import java.util.*;

public class DPFolding {
    
    private static final double INF = 1e9;
    private int minLoop = 3;
    private EnergyModel model;
    private Result result;

    public DPFolding(EnergyModel model, int minLoop) {
        this.model = model;
        this.minLoop = minLoop;
    }

    public void recurrence(String sequence) {
        String seq = sequence.toUpperCase().replaceAll("\\s+", "");
        final int n = seq.length();

        double[][] W = new double[n][n];
        double[][] V = new double[n][n];
        double[][] WM = new double[n][n];
        

        for (int i = 0; i < n; i++) {
            Arrays.fill(W[i], INF);
            Arrays.fill(V[i], INF);
            Arrays.fill(WM[i], INF);
        }

        for (int len = 2; len <= n; len++) {
            for (int i = 0; i + len - 1 < n; i++) {
                int j = i + len - 1;
                if (j - i < minLoop) {
                    W[i][j] = 0.0;
                    continue;
                }

                // V[i][j]
                if (model.canPair(seq.charAt(i), seq.charAt(j))) {
                    double bestV = INF;

                    // hairpin
                    int loopSize = j - i - 1;
                    if (loopSize >= minLoop) {
                        double hp = model.hairpinEnergy(loopSize);
                        double closePair = model.pairEnergy(seq.charAt(i), seq.charAt(j));
                        double hairpinTotal = hp + closePair;
                        bestV = Math.min(bestV, hairpinTotal);
                    }

                    // stacked
                    if (i + 1 < j - 1 && model.canPair(seq.charAt(i + 1), seq.charAt(j - 1))) {
                    	double stack = V[i + 1][j - 1]
                    	        + model.stackEnergy(seq.charAt(i), seq.charAt(i + 1), seq.charAt(j - 1), seq.charAt(j));
                        bestV = Math.min(bestV, stack);
                    }

                    // internal loops
                    double bestInternal = INF;
                    for (int i2 = i + 1; i2 <= j - minLoop - 1; i2++) {
                        for (int j2 = Math.max(i2 + minLoop + 1, i + 2); j2 < j; j2++) {
                            if (!model.canPair(seq.charAt(i2), seq.charAt(j2))) continue;
                            double eL = model.internalLoopPenalty(i - j);
                            bestInternal = Math.min(bestInternal, V[i2][j2] + eL);
                        }
                    }
                    bestV = Math.min(bestV, bestInternal);

                    // multi-loop closed
                    if (i + 1 <= j - 1) {
                        bestV = Math.min(bestV, WM[i + 1][j - 1] + model.multiInit());
                    }

                    V[i][j] = bestV;
                }

                // WM[i][j]
                double bestWM = INF;
                // i unpaired
                if (i + 1 <= j) {
                    bestWM = Math.min(bestWM, WM[i + 1][j] + model.multiUnpaired());
                }
                // j unpaired
                if (i <= j - 1) {
                    bestWM = Math.min(bestWM, WM[i][j - 1] + model.multiUnpaired());
                }
                // single branch
                
                bestWM = Math.min(bestWM, V[i][j] + model.multiBranch());
                
                // split
                for (int k = i + 1; k < j; k++) {
                	bestWM = Math.min(bestWM, WM[i][k] + WM[k + 1][j]);
                }
                
                WM[i][j] = bestWM;

                double bestW = W[i][j - 1];
                for (int k = i; k <= j - minLoop - 1; k++) {
                	double left = (k - 1 >= i) ? W[i][k - 1] : 0.0;
                	bestW = Math.min(bestW, left + V[k][j]);
                }
                W[i][j] = bestW;
            }
        }

        // Traceback
        int[] pairTo = new int[n];
        Arrays.fill(pairTo, -1);
        tracebackW(0, n - 1, W, V, WM, seq, pairTo);

        String dotBracket = Result.toDotBracket(pairTo);
        double mfe = W[0][n - 1];
        this.result = new Result(dotBracket, mfe, pairTo);
    }

    private void tracebackW(int i, int j, double[][] W, double[][] V, double[][] WM,
                            String seq, int[] pairTo) {
        if (i > j) return;
        if (i == j) return;
        double cur = W[i][j];
        
        if (Math.abs(cur - W[i][j - 1]) < 1e-3) {
            tracebackW(i, j - 1, W, V, WM, seq, pairTo);
            return;
        }
        for (int k = i; k <= j - minLoop - 1; k++) {
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
        pairTo[i] = j;
        pairTo[j] = i;

        double cur = V[i][j];

        int loopSize = j - i - 1;
        if (loopSize >= minLoop) {
            double hairpin = model.hairpinEnergy(loopSize);
            if (Math.abs(cur - hairpin) < 1e-3) {
                return;
            }
        }

        if (i + 1 < j - 1) {
            double stack = V[i + 1][j - 1]
                    + model.stackEnergy(seq.charAt(i), seq.charAt(i + 1),
                                        seq.charAt(j - 1), seq.charAt(j));
            if (Math.abs(cur - stack) < 1e-3) {
                tracebackV(i + 1, j - 1, V, WM, seq, pairTo);
                return;
            }
        }

        // internal loop
        for (int i2 = i + 1; i2 <= j - minLoop - 1; i2++) {
            for (int j2 = Math.max(i2 + minLoop + 1, i + 2); j2 < j; j2++) {
                if (!model.canPair(seq.charAt(i2), seq.charAt(j2))) continue;
                double eL = model.internalLoopPenalty(i2-j2);
                if (Math.abs(cur - (V[i2][j2] + eL)) < 1e-3) {
                    tracebackV(i2, j2, V, WM, seq, pairTo);
                    return;
                }
            }
        }

        // multiloop
        if (i + 1 <= j - 1) {
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
        double cur = WM[i][j];

        double cand = V[i][j] + model.multiBranch();
        if (Math.abs(cur - cand) < 1e-3) {
            tracebackV(i, j, V, WM, seq, pairTo);
            return;
        }

        // i unpaired
        if (i + 1 <= j) {
            cand = WM[i + 1][j] + model.multiUnpaired();
            if (Math.abs(cur - cand) < 1e-3) {
                tracebackWM(i + 1, j, WM, V, seq, pairTo);
                return;
            }
        }

        // j unpaired
        if (i <= j - 1) {
            cand = WM[i][j - 1] + model.multiUnpaired();
            if (Math.abs(cur - cand) < 1e-3) {
                tracebackWM(i, j - 1, WM, V, seq, pairTo);
                return;
            }
        }

        // split
        for (int k = i + 1; k < j; k++) {
            if (Math.abs(cur - (WM[i][k] + WM[k + 1][j])) < 1e-3) {
                tracebackWM(i, k, WM, V, seq, pairTo);
                tracebackWM(k + 1, j, WM, V, seq, pairTo);
                return;
            }
        }
    }

	public Result getResult() {
		return result;
	}

}
