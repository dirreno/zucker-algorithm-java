package main;

import java.util.HashMap;
import java.util.Map;

public class EnergyModel {
    private static final double INF = 1e9;
    
    public double multiInit() { return 0.5; }
    public double multiBranch() { return -0.2; }
    public double multiUnpaired() { return -0.05; }


    private Map<Integer, Double> hairpinEnergy = new HashMap<>();
    private Map<Integer, Double> internalLoopEnergy = new HashMap<>();
    private Map<String, Double> stack = new HashMap<>();

    public EnergyModel() {

    	// Stack
    	stack.put("AA/UU", -0.9);  stack.put("AC/UU", -1.1);  stack.put("AG/UU", -0.9);  stack.put("AU/UU", -0.9);
    	stack.put("AA/UC", -0.9);  stack.put("AC/UC", -1.8);  stack.put("AG/UC", -2.2);  stack.put("AU/UC", -0.1);
    	stack.put("AA/UG", -0.9);  stack.put("AC/UG", -1.3);  stack.put("AG/UG", -2.1);  stack.put("AU/UG", -0.6);
    	stack.put("AA/UA", -0.9);  stack.put("AC/UA", -1.1);  stack.put("AG/UA", -1.4);  stack.put("AU/UA", -0.1);

    	stack.put("CA/UG", -2.1);  stack.put("CC/UG", -2.1);  stack.put("CG/UG", -3.3);  stack.put("CU/UG", -1.4);
    	stack.put("CA/UC", -2.1);  stack.put("CC/UC", -2.4);  stack.put("CG/UC", -2.1);  stack.put("CU/UC", -2.1);
    	stack.put("CA/UA", -2.1);  stack.put("CC/UA", -3.3);  stack.put("CG/UA", -2.4);  stack.put("CU/UA", -1.4);
    	stack.put("CA/UU", -2.1);  stack.put("CC/UU", -2.1);  stack.put("CG/UU", -2.1);  stack.put("CU/UU", -2.1);

    	stack.put("GA/UC", -3.4);  stack.put("GC/UC", -3.3);  stack.put("GG/UC", -2.5);  stack.put("GU/UC", -1.5);
    	stack.put("GA/UG", -3.4);  stack.put("GC/UG", -3.3);  stack.put("GG/UG", -2.5);  stack.put("GU/UG", -1.5);
    	stack.put("GA/UA", -2.2);  stack.put("GC/UA", -2.5);  stack.put("GG/UA", -1.4);  stack.put("GU/UA", +1.3);
    	stack.put("GA/UU", -2.4);  stack.put("GC/UU", -2.5);  stack.put("GG/UU", -2.1);  stack.put("GU/UU", -0.5);

    	stack.put("UA/UA", -1.3);  stack.put("UC/UA", -1.0);  stack.put("UG/UA", -1.0);  stack.put("UU/UA", -1.3);
    	stack.put("UA/UG", -2.4);  stack.put("UC/UG", -1.5);  stack.put("UG/UG", -1.0);  stack.put("UU/UG", -1.3);
    	stack.put("UA/UC", -2.1);  stack.put("UC/UC", -1.4);  stack.put("UG/UC", +0.3);  stack.put("UU/UC", -0.9);
    	stack.put("UA/UU", -0.9);  stack.put("UC/UU", -0.6);  stack.put("UG/UU", -0.5);  stack.put("UU/UU", -0.9);


        // Hairpin
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

        // Internal loop
        internalLoopEnergy.put(1, 1.2);
        internalLoopEnergy.put(2, 1.9);
        internalLoopEnergy.put(3, 5.6);
        internalLoopEnergy.put(4, 4.9);
        internalLoopEnergy.put(5, 4.4);
        internalLoopEnergy.put(6, 4.2);
        internalLoopEnergy.put(7, 4.1);
        internalLoopEnergy.put(8, 4.0);
        internalLoopEnergy.put(9, 3.9);
        internalLoopEnergy.put(10, 3.8);
        internalLoopEnergy.put(11, 3.8);
        internalLoopEnergy.put(12, 3.7);
        internalLoopEnergy.put(13, 3.7);
        internalLoopEnergy.put(14, 3.6);
        internalLoopEnergy.put(15, 3.6);
    }

    public double pairEnergy(char a, char b) {
        if ((a == 'G' && b == 'C') || (a == 'C' && b == 'G')) return -3.0;
        if ((a == 'A' && b == 'U') || (a == 'U' && b == 'A')) return -2.0;
        if ((a == 'G' && b == 'U') || (a == 'U' && b == 'G')) return -1.0;
        return INF;
    }
    
    public double stackEnergy(char i5, char i3, char j5, char j3) {
    	String key = "" + i5 + i3 + "/" + j5 + j3;
        return stack.getOrDefault(key, 0.0);
    }
    
    public double hairpinEnergy(int loopSize) {
        if (loopSize < 3) return INF;
        if (loopSize > 15) loopSize = 15;
        return hairpinEnergy.get(loopSize);
    }
    
    public double internalLoopPenalty(int loopSize) {
        if (loopSize < 1) return INF;
        if (loopSize > 15) loopSize = 15;
        return internalLoopEnergy.get(loopSize);
    }
	
    public boolean canPair(char a, char b) {
        return pairEnergy(a, b) < INF / 2;
    }
}
