package main;


public class Main {
    public static void main(String[] args) {
        DPFolding2 folding = new DPFolding2(new EnergyModel());

        String[] sequences = {
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA",
            "GCAAUUGC",
            "GGGAAAUCC"
        };

        for (String seq : sequences) {
            Result result = folding.fold(seq);
            System.out.println("Sequence : " + seq);
            System.out.println(result);
            System.out.println();
        }
    }
}
