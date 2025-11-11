package main;


public class Main {
    public static void main(String[] args) {
        DPFolding folding = new DPFolding(new EnergyModel());

        String[] sequences = {
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA",
            "GCAAUUGC",
            "GGGAAAUCC"
        };

        for (String seq : sequences) {
            Result result = folding.recurrence(seq);
            System.out.println("Sequence : " + seq);
            System.out.println(result);
            System.out.println();
        }
    }
}
