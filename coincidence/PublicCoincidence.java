import java.io.*;
import java.util.*;
import java.util.concurrent.atomic.AtomicLong;

public class PublicCoincidence {
    private static final int CDR3AA_COL = 3, V_COL = 4, J_COL = 6,
            INCIDENCE_COL = 11; // Column layout for VDJtools format

    public static void main(String[] args) throws IOException {
        // Read parameters

        String[] opts = args[0].split(",");

        double incidenceRatioThreshold = Double.parseDouble(opts[0]),
                pValueThreshold = Double.parseDouble(opts[1]),
                logOddsThreshold = Double.parseDouble(opts[2]);

        String inputPublicListFileName = args[1],
                outputFilePrefix = args[3];

        String[] sampleFileNames = args[2].split(",");

        int nSamples = sampleFileNames.length;

        // Read the list of public clonotypes

        Map<String, ClonotypeInfo> incidenceMap = new HashMap<>();
        Map<Integer, Long> incidenceCounts = new HashMap<>();

        BufferedReader br = new BufferedReader(new FileReader(inputPublicListFileName));
        String line;

        br.readLine(); // skip header

        while ((line = br.readLine()) != null) {
            String[] splitLine = line.split("\t");

            int incidence = Integer.parseInt(splitLine[INCIDENCE_COL]);

            Long count = incidenceCounts.get(incidence);
            if (count == null) {
              count = 0L;
            }
            incidenceCounts.put(incidence, count + 1L);

            // Only add clonotypes with incidence/samples > threshold, e.g. 5%-10% of population
            if (incidence > incidenceRatioThreshold * nSamples) {
                incidenceMap.put(splitLine[CDR3AA_COL],
                        new ClonotypeInfo(splitLine[V_COL], splitLine[J_COL], new BitSet(nSamples)));
            }
        }

        // Write incidence histogram

        try (PrintWriter pw = new PrintWriter(outputFilePrefix + ".incidence.hist.txt")) {
            pw.println("incidence\tcount");
            for (Map.Entry<Integer, Long> entry : incidenceCounts.entrySet()) {
              pw.println(entry.getKey() + "\t" + entry.getValue());
            }
        }

        System.out.println("[" + (new Date()).toString() + "] Loaded " + incidenceMap.size() +
                " public clonotypes");

        // Fill in incidence bit array

        for (int i = 0; i < sampleFileNames.length; i++) {
            br = new BufferedReader(new FileReader(sampleFileNames[i]));

            br.readLine(); // skip header

            while ((line = br.readLine()) != null) {
                String[] splitLine = line.split("\t");
                ClonotypeInfo clonotypeInfo = incidenceMap.get(splitLine[CDR3AA_COL]);
                if (clonotypeInfo != null)
                    clonotypeInfo.setIncidence(i);
            }

            if (i % 10 == 0) {
                System.out.println("[" + (new Date()).toString() + "] Scanned " + i + " of " +
                        sampleFileNames.length +
                        " samples for public incidence.");
            }
        }

        // Run pairwise comparisons

        final Map<String, AtomicLong> volcanoBg = generateVolcanoBg();

        final AtomicLong pairsCounter = new AtomicLong(), goodPairsCounter = new AtomicLong();

        try (PrintWriter pw = new PrintWriter(outputFilePrefix + ".txt")) {
            pw.println("cdr3aa.1\tv.1\tj.1\tcdr3aa.2\tv.2\tj.2\tn12\tn1\tn2\tn.total\tlog.odds\tlog.p.value");
            incidenceMap.entrySet().parallelStream().forEach(entry1 ->
                    incidenceMap.forEach((key2, value2) -> {
                        String key1 = entry1.getKey();
                        ClonotypeInfo value1 = entry1.getValue();
                        if (key1.compareTo(key2) > 0) { // ignore duplicates
                            int n1 = value1.getIncidenceCount(),
                                    n2 = value2.getIncidenceCount();

                            if (n1 > 0 && n2 > 0) {
                                int n12 = value1.getCoincidenceCount(value2);

                                double logOdds = computeLogOdds(n12, n1, n2, nSamples),
                                        pValue = computeP(n12, n1, n2, nSamples),
                                        logPValue = Math.log10(pValue + 1e-100); // can actually hash P-values to speed up

                                AtomicLong bgCounter = volcanoBg.get(getLogPCoord(logPValue) + "\t" +
                                        getLogOddsCoord(logOdds));

                                if (bgCounter != null) {
                                    bgCounter.incrementAndGet();
                                }

                                if (pValue <= pValueThreshold & Math.abs(logOdds) >= logOddsThreshold) {
                                    pw.println(key1 + "\t" + value1.v + "\t" + value1.j + "\t" +
                                            key2 + "\t" + value2.v + "\t" + value2.j + "\t" +
                                            n12 + "\t" +
                                            n1 + "\t" +
                                            n2 + "\t" +
                                            nSamples + "\t" +
                                            logOdds + "\t" +
                                            logPValue
                                    );
                                    goodPairsCounter.incrementAndGet();
                                }
                            }

                            if (pairsCounter.incrementAndGet() % 50_000_000 == 0) {
                                System.out.println("[" + (new Date()).toString() + "] Checked ~" + pairsCounter.get() +
                                        " pairs out of " + (long) incidenceMap.size() * ((incidenceMap.size() - 1) / 2L) +
                                        ", " + goodPairsCounter.get() +
                                        " pairs passing minimal filtering criteria.");
                            }
                        }
                    })
            );
        }

        // Write volcano BG

        try (PrintWriter pw = new PrintWriter(outputFilePrefix + ".volcano.txt")) {
            pw.println("log.p\tlog.odds\tcount");
            for (int logPCoord = 0; logPCoord <= LOG_P_BINS; logPCoord++) {
                for (int logOddsCoord = 0; logOddsCoord <= LOG_ODDS_BINS; logOddsCoord++) {
                    pw.println(getLogPValue(logPCoord) + "\t" +
                            getLogOddsValue(logOddsCoord) + "\t" +
                            volcanoBg.get(logPCoord + "\t" + logOddsCoord).get());
                }
            }
        }

        System.out.println("[" + (new Date()).toString() + "] DONE. Checked " + pairsCounter.get() +
                " pairs out of " + incidenceMap.size() * (incidenceMap.size() - 1) / 2 +
                ", " + goodPairsCounter.get() +
                " pairs passing minimal filtering criteria.");
    }

    private static double computeLogOdds(int n12, int n1, int n2, int nSamples) {
        return Math.max(-20, Math.min(20,
                Math.log10(((double) n12 * nSamples) / n1 / n2)));
    }

    private static double computeP(int n12, int n1, int n2, int nSamples) {
        double p = pHyp(n12, n1, n2, nSamples);

        return Math.min(p, 1 - p);
    }

    private static double pHyp(int n12, int n1, int n2, int nSamples) {
        // http://journals.sagepub.com.sci-hub.cc/doi/pdf/10.2466/pms.1998.87.1.51

        int v = Math.max(0, n1 + n2 - nSamples),
                w = Math.min(n1, n2);

        assert n12 >= v && n12 <= w;

        double pPrev = 1, T = pPrev, S = pPrev;

        for (int i = v + 1; i <= w; i++) {
            pPrev *= ((double) (n1 - i + 1) * (n2 - i + 1)) / i / (nSamples - n1 - n2 + i);
            T += pPrev;

            if (i == n12) {
                S = T - 0.5 * pPrev;
            }
        }

        return S / T;
    }

    private static Map<String, AtomicLong> generateVolcanoBg() {
        Map<String, AtomicLong> volcanoBg = new HashMap<>();
        for (int logPCoord = 0; logPCoord <= LOG_P_BINS; logPCoord++) {
            for (int logOddsCoord = 0; logOddsCoord <= LOG_ODDS_BINS; logOddsCoord++) {
                volcanoBg.put(logPCoord + "\t" + logOddsCoord,
                        new AtomicLong());
            }
        }
        return volcanoBg;
    }

    static final int LOG_ODDS_BINS = 100, LOG_P_BINS = 100;
    static final double MIN_LOG_ODDS = -1, MAX_LOG_ODDS = 1,
            MIN_LOG_P = -10, MAX_LOG_P = 0;

    static int getLogOddsCoord(double logOdds) {
        return (int) ((logOdds - MIN_LOG_ODDS) * LOG_ODDS_BINS / (MAX_LOG_ODDS - MIN_LOG_ODDS));
    }

    static double getLogOddsValue(int coord) {
        return MIN_LOG_ODDS + coord / (double) LOG_ODDS_BINS * (MAX_LOG_ODDS - MIN_LOG_ODDS);
    }

    static int getLogPCoord(double logP) {
        return (int) ((logP - MIN_LOG_P) * LOG_P_BINS / (MAX_LOG_P - MIN_LOG_P));
    }

    static double getLogPValue(int coord) {
        return MIN_LOG_P + coord / (double) LOG_P_BINS * (MAX_LOG_P - MIN_LOG_P);
    }

    private static class ClonotypeInfo {
        final String v, j;
        final BitSet incidence;

        ClonotypeInfo(String v, String j, BitSet incidence) {
            this.v = v;
            this.j = j;
            this.incidence = incidence;
        }

        void setIncidence(int sampleIndex) {
            incidence.set(sampleIndex);
        }

        int getIncidenceCount() {
            return incidence.cardinality();
        }

        int getCoincidenceCount(ClonotypeInfo other) {
            BitSet tmp = (BitSet) incidence.clone();
            tmp.and(other.incidence);
            return tmp.cardinality();
        }
    }
}
