package main;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

public class Main {
    public static void main(String[] args) {
        DPFolding folding = new DPFolding(new EnergyModel(),3);
        String[] sequences = fastaToString("data/rnas.fasta");
        String outputFilename = "data/structures_output.txt";
        java.io.File outputFile = new java.io.File(outputFilename);
        if (outputFile.getParentFile() != null) outputFile.getParentFile().mkdirs();
        try (FileWriter fw = new FileWriter(outputFile, false)) {
        } catch (IOException e) {
            e.printStackTrace();
        }
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
            
            StringBuilder sb = new StringBuilder();
            sb.append(result.getStructure()).append(System.lineSeparator());
            String outStr = sb.toString();
            saveStructure(outputFilename, outStr);
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
    
    public static void saveStructure(String filename, String content) {
        java.io.File outFile = new java.io.File(filename);
        if (outFile.getParentFile() != null) outFile.getParentFile().mkdirs();
        try (BufferedWriter out = new BufferedWriter(new FileWriter(outFile, true))) {
            out.write(content);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
