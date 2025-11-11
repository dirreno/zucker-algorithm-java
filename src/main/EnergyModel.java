package main;

import java.util.HashMap;
import java.util.Map;

public class EnergyModel {
    private static final double INF = 1e9;
    
    public double multiInit() { return 2.0; }      // a
    public double multiBranch() { return 0.1; }    // b
    public double multiUnpaired() { return 0.1; }  // c


    private final Map<String, Double> helixTerminal = new HashMap<>();
    private final Map<Integer, Double> hairpinEnergy = new HashMap<>();
    private final Map<Integer, Double> internalLoopEnergy = new HashMap<>();

    public EnergyModel() {

        // --- Hairpin energies ---
        hairpinEnergy.put(3, 5.4);
        hairpinEnergy.put(4, 5.6);
        hairpinEnergy.put(5, 5.7);
        hairpinEnergy.put(6, 5.4);
        hairpinEnergy.put(7, 6.0);
        hairpinEnergy.put(8, 5.5);
        hairpinEnergy.put(9, 6.4);
        hairpinEnergy.put(10, 6.5);
        hairpinEnergy.put(11, 6.6);
        hairpinEnergy.put(12, 6.7);
        hairpinEnergy.put(13, 6.8);
        hairpinEnergy.put(14, 6.9);
        hairpinEnergy.put(15, 6.9);

        // --- Bulge energies ---
        internalLoopEnergy.put(1, 8.8);
        internalLoopEnergy.put(2, 8.8);
        internalLoopEnergy.put(3, 8.2);
        internalLoopEnergy.put(4, 5.6);
        internalLoopEnergy.put(5, 5.7);
        internalLoopEnergy.put(6, 5.4);
        internalLoopEnergy.put(7, 6.0);
        internalLoopEnergy.put(8, 5.5);
        internalLoopEnergy.put(9, 6.4);
        internalLoopEnergy.put(10, 6.5);
        internalLoopEnergy.put(11, 6.6);
        internalLoopEnergy.put(12, 6.7);
        internalLoopEnergy.put(13, 6.8);
        internalLoopEnergy.put(14, 6.9);
        internalLoopEnergy.put(15, 6.9);
    }

    /** Base pair energy (Watson–Crick + wobble) */
    public double pairEnergy(char a, char b) {
        if ((a == 'G' && b == 'C') || (a == 'C' && b == 'G')) return -3.0;
        if ((a == 'A' && b == 'U') || (a == 'U' && b == 'A')) return -2.0;
        if ((a == 'G' && b == 'U') || (a == 'U' && b == 'G')) return -1.0;
        return INF;
    }
    
    /** Hairpin loop energy */
    public double hairpinEnergy(int loopSize) {
        if (loopSize < 3) return INF;
        if (loopSize > 15) loopSize = 15;
        return hairpinEnergy.get(loopSize);
    }
    
    public double stackEnergy(char i, char i1, char j1, char j) {
        double basePair = pairEnergy(i, j);
        String tetramer = "" + i + i1 + j1 + j;
        double terminal = terminalHelixEnergy(tetramer);
        return basePair + terminal;
    }

    /** Bulge loop penalty */
    public double InternalLoopPenalty(int loopSize) {
        if (loopSize < 1) return INF;
        if (loopSize > 15) loopSize = 15;
        return internalLoopEnergy.get(loopSize);
    }
	
    /** Check pair compatibility */
    public boolean canPair(char a, char b) {
        return pairEnergy(a, b) < INF / 2;
    }
}
