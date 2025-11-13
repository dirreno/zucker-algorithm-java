package main;


public class Main {
    public static void main(String[] args) {
        DPFolding folding = new DPFolding(new EnergyModel(),3);

        String[] sequences = {
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA",
            "UUCUUUUUUAGUGGCAGUAAGCCUGGGAAUGGGGGCGACCCAGGCGUAUGAACAUAGUGUAACGCUCCCC"
        };
        for (String seq : sequences) {
        	folding.recurrence(seq);
            Result result = folding.getResult();
            System.out.println("Sequence : " + seq);
            System.out.println(result);
            System.out.println();
        }
    }

    public static int distance(String s1, String s2) {
        if (s1 == null || s2 == null) throw new IllegalArgumentException("Inputs must be non-null");
        if (s1.length() != s2.length()) throw new IllegalArgumentException("Hamming distance requires equal-length strings");
        int diff = 0;
        for (int i = 0; i < s1.length(); i++) {
            if (s1.charAt(i) != s2.charAt(i)) diff++;
        }
        return diff;
    }
}
