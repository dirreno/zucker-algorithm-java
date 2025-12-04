package main;

public class TestRun {
    public static void main(String[] args) {
        DPFolding folding = new DPFolding(new EnergyModel(), 3);
        String seq = "GGCGGUGAAAUGCC";
        folding.recurrence(seq);
        Result res = folding.getResult();
        System.out.println(res);
    }
}
