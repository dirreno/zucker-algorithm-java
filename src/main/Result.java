package main;

public class Result {
	public String structure;
    public double mfe;
    public int[] pairTo;

    public Result(String structure, double mfe, int[] pairTo) {
        this.structure = structure;
        this.mfe = mfe;
        this.pairTo = pairTo;
    }

    @Override
    public String toString() {
        return String.format("Dot-bracket: %s\nMFE: %.3f kcal/mol", structure, mfe);
    }

    public static String toDotBracket(int[] pairTo) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < pairTo.length; i++) {
            int p = pairTo[i];
            sb.append(p == -1 ? '.' : (p > i ? '(' : ')'));
        }
        return sb.toString();
    }

	public String getStructure() {
		return structure;
	}

	public double getMfe() {
		return mfe;
	}

	public int[] getPairTo() {
		return pairTo;
	}
}
