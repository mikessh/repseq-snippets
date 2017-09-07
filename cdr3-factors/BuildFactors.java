import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Created by mikesh on 9/7/17.
 */
public class BuildFactors {
    private static final int CDR3AA_COL = 3,
            READS_COL = 0;

    public static void main(String[] args) throws IOException {
        String[] sampleFileNames = args[0].split(",");

        int nSamples = sampleFileNames.length;

        // Read the list of public clonotypes

        BufferedReader br;
        String line;

        List<Map<String, Counter>> table = new ArrayList<>();
        List<Counter> sampleTotals = new ArrayList<>();

        for (int i = 0; i < sampleFileNames.length; i++) {
            Map<String, Counter> sampleTable = new HashMap<>();
            Counter sampleTotal = new Counter();

            br = new BufferedReader(new FileReader(sampleFileNames[i]));

            br.readLine(); // skip header

            while ((line = br.readLine()) != null) {
                String[] splitLine = line.split("\t");
                String cdr3aa = splitLine[CDR3AA_COL];
                int len = cdr3aa.length();
                int reads = Integer.parseInt(splitLine[READS_COL]);

                for (int j = 0; j < len; j++) {
                    String signature = len + "\t" + j + "\t" + cdr3aa.charAt(i);
                    Counter counter = sampleTable.computeIfAbsent(signature, k -> new Counter());
                    counter.reads += reads;
                    counter.unique++;
                }

                sampleTotal.unique++;
                sampleTotal.reads += reads;
            }

            table.add(sampleTable);
            sampleTotals.add(sampleTotal);

            if (i % 10 == 0) {
                System.out.println("[" + (new Date()).toString() + "] Scanned " + i + " of " +
                        sampleFileNames.length +
                        " samples.");
            }
        }

        try (PrintWriter pw = new PrintWriter(args[1])) {
            pw.println("sample\tlen\tpos\taa\tunique\treads\tunique.total\treads.total");
            for (int i = 0; i < sampleFileNames.length; i++) {
                Map<String, Counter> sampleTable = table.get(i);
                Counter sampleTotal = sampleTotals.get(i);

                for (Map.Entry<String, Counter> entry : sampleTable.entrySet()) {
                    pw.println(i + "\t" + entry.getKey() + "\t" + entry.getValue() + "\t" + sampleTotal);
                }
            }
        }
    }

    static class Counter {
        int unique;
        long reads;

        @Override
        public String toString() {
            return unique + "\t" + reads;
        }
    }
}
