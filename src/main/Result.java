package main;

public class Result {
	public final String structure;
    public final double mfe;
    public final int[] pairTo;

    public Result(String structure, double mfe, int[] pairTo) {
        this.structure = structure;
        this.mfe = mfe;
        this.pairTo = pairTo;
    }

    @Override
    public String toString() {
        return String.format("Dot-bracket: %s\nMFE: %.3f kcal/mol", structure, mfe);
    }
}
