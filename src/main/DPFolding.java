package main;

import java.util.Arrays;

public class DPFolding {
    private static final int MIN_LOOP = 3;
    private static final double INF = 1e9;

    private final EnergyModel model;

    public DPFolding(EnergyModel model) {
        this.model = model;
    }

    public Result recurrence(String seq) {
        seq = seq.toUpperCase().replaceAll("\\s+", "");
        final int n = seq.length();

        double[][] dp = new double[n][n];
        int[][] action = new int[n][n];
        int[][] actionIdx = new int[n][n];

        for (double[] row : dp) Arrays.fill(row, INF);
        for (int i = 0; i < n; i++) dp[i][i] = 0.0;

        for (int len = 2; len <= n; len++) {
            for (int i = 0; i + len - 1 < n; i++) {
                int j = i + len - 1;
                double best = INF;
                int bestAction = 0;
                int bestIdx = -1;

                // j unpaired
                double opt1 = (j - 1 >= i) ? dp[i][j - 1] : 0.0;
                if (opt1 < best) { best = opt1; bestAction = 1; }

                // j paired with k
                for (int k = i; k <= j - MIN_LOOP - 1; k++) {
                    char ak = seq.charAt(k), aj = seq.charAt(j);
                    if (!model.canPair(ak, aj)) continue;

                    double ePair = model.pairEnergy(ak, aj);
                    
                    double loopSize = j - k - 1;
                    double total = ((k - 1 >= i) ? dp[i][k - 1] : 0.0)
                                 + ePair
                                 + (loopSize <= MIN_LOOP ? model.hairpinEnergy((int) loopSize) : dp[k + 1][j - 1]);

                    if (total < best) { best = total; bestAction = 2; bestIdx = k; }
                }

                // split
                for (int t = i; t < j; t++) {
                    double total = dp[i][t] + dp[t + 1][j];
                    if (total < best) { best = total; bestAction = 3; bestIdx = t; }
                }

                dp[i][j] = best;
                action[i][j] = bestAction;
                actionIdx[i][j] = bestIdx;
            }
        }

        int[] pairTo = new int[n];
        Arrays.fill(pairTo, -1);
        if (n > 0) traceback(0, n - 1, action, actionIdx, pairTo, seq);

        String dotBracket = toDotBracket(pairTo);
        double mfe = (n > 0) ? dp[0][n - 1] : 0.0;

        return new Result(dotBracket, mfe, pairTo);
    }

    private void traceback(int i, int j, int[][] action, int[][] idx, int[] pairTo, String seq) {
        if (i >= j) return;
        int act = action[i][j];
        switch (act) {
            case 1 -> traceback(i, j - 1, action, idx, pairTo, seq);
            case 2 -> {
                int k = idx[i][j];
                pairTo[k] = j; pairTo[j] = k;
                if (k - 1 >= i) traceback(i, k - 1, action, idx, pairTo, seq);
                if (k + 1 <= j - 1) traceback(k + 1, j - 1, action, idx, pairTo, seq);
            }
            case 3 -> {
                int t = idx[i][j];
                traceback(i, t, action, idx, pairTo, seq);
                traceback(t + 1, j, action, idx, pairTo, seq);
            }
            default -> { if (j - 1 >= i) traceback(i, j - 1, action, idx, pairTo, seq); }
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
