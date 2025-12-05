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
        recurrence(sequence, false);
    }

    public void recurrence(String sequence, boolean verbose) {
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
                            double eL =  model.internalLoopPenalty((i2 - i - 1) + (j - j2 - 1));
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
            if (verbose) {
                System.out.println("After len=" + len + ":");
                printMatrix("W", W, seq, n);
                printMatrix("V", V, seq, n);
                printMatrix("WM", WM, seq, n);
            }
        }

        // Traceback
        int[] pairTo = new int[n];
        Arrays.fill(pairTo, -1);
        tracebackW(0, n - 1, W, V, WM, seq, pairTo, verbose);

        if (verbose) {
            System.out.println("pairTo: " + Arrays.toString(pairTo));
        }

        String dotBracket = Result.toDotBracket(pairTo);
        double mfe = W[0][n - 1];
        this.result = new Result(dotBracket, mfe, pairTo);
    }

    private void tracebackW(int i, int j, double[][] W, double[][] V, double[][] WM,
                            String seq, int[] pairTo, boolean verbose) {
        if (i > j) return;
        if (i == j) return;
        double cur = W[i][j];
        if (verbose) System.out.printf("tracebackW: i=%d j=%d cur=%.2f%n", i, j, cur);
        
        if (Math.abs(cur - W[i][j - 1]) < 1e-3) {
            tracebackW(i, j - 1, W, V, WM, seq, pairTo, verbose);
            return;
        }
        for (int k = i; k <= j - minLoop - 1; k++) {
            double left = (k - 1 >= i) ? W[i][k - 1] : 0.0;
            if (Math.abs(cur - (left + V[k][j])) < 1e-3) {
                if (verbose) System.out.printf("  W-match: split at k=%d (left=%.2f V[%d][%d]=%.2f)%n", k, left, k, j, V[k][j]);
                tracebackV(k, j, V, WM, seq, pairTo, verbose);
                if (k - 1 >= i) tracebackW(i, k - 1, W, V, WM, seq, pairTo, verbose);
                return;
            }
        }
        if (verbose) System.out.printf("  W-match: skip j (W[%d][%d-1]=%.2f)%n", i, j, W[i][j-1]);
        tracebackW(i, j - 1, W, V, WM, seq, pairTo, verbose);
    }

    private void tracebackV(int i, int j, double[][] V, double[][] WM,
                            String seq, int[] pairTo, boolean verbose) {
        pairTo[i] = j;
        pairTo[j] = i;

        double cur = V[i][j];
        if (verbose) System.out.printf("tracebackV: i=%d j=%d cur=%.2f%n", i, j, cur);

        int loopSize = j - i - 1;
        if (loopSize >= minLoop) {
            double hairpin = model.hairpinEnergy(loopSize);
            double closePair = model.pairEnergy(seq.charAt(i), seq.charAt(j));
            double hairpinTotal = hairpin + closePair;
            if (Math.abs(cur - hairpinTotal) < 1e-3) {
                if (verbose) System.out.println("  V-match: hairpin");
                return;
            }
        }

        if (i + 1 < j - 1) {
            double stack = V[i + 1][j - 1]
                    + model.stackEnergy(seq.charAt(i), seq.charAt(i + 1),
                                        seq.charAt(j - 1), seq.charAt(j));
            if (Math.abs(cur - stack) < 1e-3) {
                if (verbose) System.out.printf("  V-match: stack -> recurse (%d,%d)%n", i+1, j-1);
                tracebackV(i + 1, j - 1, V, WM, seq, pairTo, verbose);
                return;
            }
        }

        // internal loop
        for (int i2 = i + 1; i2 <= j - minLoop - 1; i2++) {
            for (int j2 = Math.max(i2 + minLoop + 1, i + 2); j2 < j; j2++) {
                if (!model.canPair(seq.charAt(i2), seq.charAt(j2))) continue;
                double eL = model.internalLoopPenalty((i2 - i - 1) + (j - j2 - 1));
                if (Math.abs(cur - (V[i2][j2] + eL)) < 1e-3) {
                    if (verbose) System.out.printf("  V-match: internal loop -> recurse (%d,%d)%n", i2, j2);
                    tracebackV(i2, j2, V, WM, seq, pairTo, verbose);
                    return;
                }
            }
        }

        // multiloop
        if (i + 1 <= j - 1) {
            double multi = WM[i + 1][j - 1] + model.multiInit();
            if (Math.abs(cur - multi) < 1e-3) {
                if (verbose) System.out.println("  V-match: multiloop");
                tracebackWM(i + 1, j - 1, WM, V, seq, pairTo, verbose);
                return;
            }
        }
    }

    private void tracebackWM(int i, int j, double[][] WM, double[][] V,
                             String seq, int[] pairTo, boolean verbose) {
        if (i > j) return;
        double cur = WM[i][j];
        if (verbose) System.out.printf("tracebackWM: i=%d j=%d cur=%.2f%n", i, j, cur);

        double cand = V[i][j] + model.multiBranch();
        if (Math.abs(cur - cand) < 1e-3) {
            if (verbose) System.out.println("  WM-match: single branch (V)");
            tracebackV(i, j, V, WM, seq, pairTo, verbose);
            return;
        }

        // i unpaired
        if (i + 1 <= j) {
            cand = WM[i + 1][j] + model.multiUnpaired();
            if (Math.abs(cur - cand) < 1e-3) {
                if (verbose) System.out.println("  WM-match: i unpaired");
                tracebackWM(i + 1, j, WM, V, seq, pairTo, verbose);
                return;
            }
        }

        // j unpaired
        if (i <= j - 1) {
            cand = WM[i][j - 1] + model.multiUnpaired();
            if (Math.abs(cur - cand) < 1e-3) {
                if (verbose) System.out.println("  WM-match: j unpaired");
                tracebackWM(i, j - 1, WM, V, seq, pairTo, verbose);
                return;
            }
        }

        // split
        for (int k = i + 1; k < j; k++) {
            if (Math.abs(cur - (WM[i][k] + WM[k + 1][j])) < 1e-3) {
                if (verbose) System.out.printf("  WM-match: split at k=%d%n", k);
                tracebackWM(i, k, WM, V, seq, pairTo, verbose);
                tracebackWM(k + 1, j, WM, V, seq, pairTo, verbose);
                return;
            }
        }
    }

    private void printMatrix(String name, double[][] M, String seq, int n) {
        System.out.println(name + ":");
        // horizontal header
        System.out.print("    ");
        for (int j = 0; j < n; j++) {
            System.out.printf("%6c", seq.charAt(j));
        }
        System.out.println();

        for (int i = 0; i < n; i++) {
            // vertical header (nucleotide at start of row)
            System.out.printf(" %3c", seq.charAt(i));
            for (int j = 0; j < n; j++) {
                if (j < i) {
                    System.out.print("      ");
                } else {
                    double v = M[i][j];
                    if (v > 1e8) System.out.printf("%6s", "INF");
                    else System.out.printf("%6.2f", v);
                }
            }
            System.out.println();
        }
    }

	public Result getResult() {
		return result;
	}

}
