package main;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

public class Main {
    public static void main(String[] args) {
        DPFolding folding = new DPFolding(new EnergyModel(),3);
        String[] sequences = fastaToString("data/rnas.fasta");
        for (String seq : sequences) {
            long start = System.currentTimeMillis();
            folding.recurrence(seq);
            Result result = folding.getResult();
            long end = System.currentTimeMillis();
            long ms = end - start;
            System.out.println("Sequence : " + seq);
            System.out.printf("Processing time: %d ms%n", ms);
            System.out.println(result);
            System.out.println();
        }
    }

    public static int distance(String s1, String s2) {
        int diff = 0;
        for (int i = 0; i < s1.length(); i++) {
            if (s1.charAt(i) != s2.charAt(i)) diff++;
        }
        return diff;
    }
    
    public static String[] fastaToString(String filename) {
        File fastaFile = new File(filename);
        List<String> seqs = new ArrayList<>();
        try (FastaSequenceFile fasta = new FastaSequenceFile(fastaFile, true)) {
            ReferenceSequence seq;
            while ((seq = fasta.nextSequence()) != null) {
                seqs.add(new String(seq.getBases()));
            }
        } catch (Exception e) {
            return new String[0];
        }
        return seqs.toArray(new String[0]);
    }
}
