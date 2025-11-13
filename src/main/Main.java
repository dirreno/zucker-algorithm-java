package main;


public class Main {
    public static void main(String[] args) {
        DPFolding folding = new DPFolding(new EnergyModel(),3);

        String[] sequences = {
        		"GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA",
        		"UUCUUUUUUAGUGGCAGUAAGCCUGGGAAUGGGGGCGACCCAGGCGUAUGAACAUAGUGUAACGCUCCCC",
        };

        Result previous = null;
        for (String seq : sequences) {
            Result result = folding.recurrence(seq);
            System.out.println("Sequence : " + seq);
            System.out.println(result);
            if (previous != null) {
                int d = Result.distance(previous.structure, result.structure);
                System.out.println("Distance to previous structure: " + d);
            }
            System.out.println();
            previous = result;
        }
    }
}
