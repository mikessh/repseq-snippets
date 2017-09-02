import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.*;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Stream;

import com.milaboratory.core.alignment.Alignment;
import com.milaboratory.core.mutations.MutationType;
import com.milaboratory.core.mutations.Mutations;
import com.milaboratory.core.sequence.AminoAcidSequence;
import com.milaboratory.core.tree.NeighborhoodIterator;
import com.milaboratory.core.tree.SequenceTreeMap;
import com.milaboratory.core.tree.TreeSearchParameters;

public class AlignCdrAux {
    public static void main(String[] args) throws IOException {
        int maxSubstitutions = Integer.parseInt(args[0]),
                maxIndels = Integer.parseInt(args[1]);

        String inputFileName = args[2], outputFileName = args[3];

        final Map<AminoAcidSequence, Cdr3Info> cdr3AntigenMap = new HashMap<>();
        try (Stream<String> stream = Files.lines(new File(inputFileName).toPath())) {
            stream.forEach(line -> {
                String[] splitString = line.split("\t");
                cdr3AntigenMap.compute(new AminoAcidSequence(splitString[0]),
                        (aminoAcidSequence, cdr3Info) -> {
                            if (cdr3Info == null) {
                                cdr3Info = new Cdr3Info(aminoAcidSequence,
                                        Integer.parseInt(splitString[1]),
                                        Integer.parseInt(splitString[2]),
                                        splitString[3],
                                        splitString[5]);
                            }
                            cdr3Info.addAntigen(splitString[4]);
                            return cdr3Info;
                        });
            });
        }

        System.out.println("Loaded " + cdr3AntigenMap.size() + " cdr3 sequences.");

        final SequenceTreeMap<AminoAcidSequence, Cdr3Info> stm = new SequenceTreeMap<>(AminoAcidSequence.ALPHABET);

        cdr3AntigenMap.entrySet().forEach(kvp -> stm.put(kvp.getKey(), kvp.getValue()));

        final TreeSearchParameters tsp = new TreeSearchParameters(maxSubstitutions, maxIndels, maxIndels);

        String header =
                "unknown.ag\tsame.ag\tcdr3.len\tweight\t" +
                        "align.id\talign.sub.id\t" +
                        "subst\tins\tdel\t" +
                        "mut.type\tmut.pos\tmut.from\tmut.to\t" +
                        "mut.region\tgene\tset";


        try (final PrintWriter pw = new PrintWriter(new File(outputFileName))) {
            pw.println(header);

            final Queue<String> lines = new ConcurrentLinkedQueue<>();

            Thread writeThread = new Thread(() -> {
                while (true) {
                    String value;
                    while ((value = lines.poll()) != null) {
                        if (value.equals("END")) {
                            return;
                        }
                        pw.println(value);
                    }
                }
            });
            writeThread.start();

            final AtomicInteger counter = new AtomicInteger(), mutationCounter = new AtomicInteger(),
                    alignmentIdCounter = new AtomicInteger();

            cdr3AntigenMap.values().parallelStream().forEach(thisCdr3Info -> {
                        AminoAcidSequence thisCdr3 = thisCdr3Info.cdr3;
                        Cdr3Info otherCdr3Info;

                        if (!thisCdr3Info.unknownAntigen()) {
                            NeighborhoodIterator<AminoAcidSequence, Cdr3Info> iter = stm
                                    .getNeighborhoodIterator(thisCdr3, tsp);

                            Map<AminoAcidSequence, List<AlignmentInfo>> alignmentVariants = new HashMap<>();

                            while ((otherCdr3Info = iter.next()) != null) {
                                Alignment<AminoAcidSequence> alignment = iter.getCurrentAlignment();

                                if (thisCdr3Info.eligibleComparison(otherCdr3Info) &&
                                        alignment.getSequence1Range().length() == thisCdr3.size()) { // enforce global (JIC)
                                    final AlignmentInfo alignmentInfo = new AlignmentInfo(
                                            thisCdr3Info.antigensOverlap(otherCdr3Info),
                                            alignment,
                                            thisCdr3Info.unknownAntigen() || otherCdr3Info.unknownAntigen());

                                    alignmentVariants.compute(otherCdr3Info.cdr3,
                                            (cdr3, alignments) -> {
                                                if (alignments == null) {
                                                    alignments = new ArrayList<>();
                                                }
                                                alignments.add(alignmentInfo);
                                                return alignments;
                                            });
                                }
                            }

                            for (List<AlignmentInfo> alignmentInfos : alignmentVariants.values()) {
                                float weight = 1.0f / alignmentInfos.size();
                                int alignmentId = alignmentIdCounter.incrementAndGet(),
                                        alignmentSubId = 0;

                                for (AlignmentInfo alignmentInfo : alignmentInfos) {
                                    Set<Integer> mutatedPositions = new HashSet<>();
                                    Alignment alignment = alignmentInfo.alignment;
                                    AminoAcidSequence reference = (AminoAcidSequence) alignment.getSequence1();

                                    String prefix = (alignmentInfo.unknownAntigen ? "TRUE" : "FALSE") + "\t" +
                                            (alignmentInfo.sameAntigen ? "TRUE" : "FALSE") + "\t" +
                                            reference.size() + "\t" + weight + "\t" +
                                            alignmentId + "\t" + (alignmentSubId++) + "\t";

                                    Mutations mutations = alignment.getAbsoluteMutations();

                                    // Count true number of mismatches

                                    int subst = 0, ins = 0, del = 0;

                                    for (int k = 0; k < mutations.size(); k++) {
                                        switch (mutations.getTypeByIndex(k)) {
                                            case Substitution:
                                                subst++;
                                                break;
                                            case Insertion:
                                                ins++;
                                                break;
                                            case Deletion:
                                                del++;
                                                break;
                                        }
                                    }

                                    prefix += subst + "\t" + ins + "\t" + del + "\t";

                                    for (int k = 0; k < mutations.size(); k++) {
                                        MutationType mutationType = mutations.getTypeByIndex(k);
                                        int pos = mutations.getPositionByIndex(k);
                                        boolean isInsertion = mutationType == MutationType.Insertion;

                                        lines.add(prefix +
                                                shortMutationType(mutationType) + "\t" +
                                                pos + "\t" +
                                                (isInsertion ?
                                                        "-" : mutations.getFromAsSymbolByIndex(k)) + "\t" +
                                                (mutationType == MutationType.Deletion ?
                                                        "-" : mutations.getToAsSymbolByIndex(k)) + "\t" +
                                                thisCdr3Info.getRegion(pos) + "\t" +
                                                thisCdr3Info.gene + "\t" +
                                                thisCdr3Info.set
                                        );

                                        mutationCounter.incrementAndGet();

                                        // Store mismatched positions

                                        if (!isInsertion) {
                                            mutatedPositions.add(pos);
                                        }
                                    }

                                    // Exact matches

                                    for (int i = 0; i < reference.size(); i++) {
                                        char aa = thisCdr3.symbolAt(i);
                                        if (!mutatedPositions.contains(i)) {
                                            lines.add(prefix +
                                                    "E\t" +
                                                    i + "\t" +
                                                    aa + "\t" +
                                                    aa + "\t" +
                                                    thisCdr3Info.getRegion(i) + "\t" +
                                                    thisCdr3Info.gene + "\t" +
                                                    thisCdr3Info.set
                                            );
                                        }
                                    }
                                }
                            }

                            int count = counter.incrementAndGet();

                            if (count % 100 == 0) {
                                System.out.println("[" + (new Date()) + "] " +
                                        "Queried " + count + " of " + cdr3AntigenMap.size() + " cdr3 sequences. " +
                                        "Recorded ~" + mutationCounter.get() + " mutations so far.");
                            }
                        }
                    }
            );

            lines.add("END");

            writeThread.join();

            System.out.println("[" + (new Date()) + "] " +
                    "Done. Queried " + cdr3AntigenMap.size() + " cdr3 sequences. " +
                    "Recorded " + mutationCounter.get() + " mutations.");
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    private static class Cdr3Info {
        final static String NA_ANTIGEN_CHAR = ".";

        final AminoAcidSequence cdr3;
        final Set<String> antigens = new HashSet<>();
        final int vEnd, jStart;
        final String gene, set;

        Cdr3Info(AminoAcidSequence cdr3, int vEnd, int jStart, String gene,
                 String set) {
            this.cdr3 = cdr3;
            this.vEnd = vEnd;
            this.jStart = jStart;
            this.gene = gene;
            this.set = set;
        }

        void addAntigen(String antigen) {
            antigens.add(antigen);
        }

        boolean antigensOverlap(Cdr3Info other) {
            return overlaps(this.antigens, other.antigens);
        }

        boolean unknownAntigen() {
            return antigens.size() == 1 && antigens.contains(NA_ANTIGEN_CHAR);
        }

        boolean eligibleComparison(Cdr3Info otherCdr3Info) {
            return gene.equals(otherCdr3Info.gene) &&
                    // unknown antigens are not searched, so just remove exact match
                    ((otherCdr3Info.unknownAntigen() && cdr3.compareTo(otherCdr3Info.cdr3) != 0) ||
                            // remove exact matches and duplicate comparisons
                            cdr3.compareTo(otherCdr3Info.cdr3) > 0);
        }

        private static boolean overlaps(Set<String> set1, Set<String> set2) {
            if (set1.size() > set2.size()) {
                return overlaps(set2, set1);
            }
            for (String value : set1) {
                if (!value.equals(NA_ANTIGEN_CHAR) && set2.contains(value))
                    return true;
            }
            return false;
        }

        String getRegion(int pos) {
            return pos < vEnd ? "V" : (pos > jStart ? "J" : "N");
        }
    }

    private static class AlignmentInfo {
        final boolean sameAntigen, unknownAntigen;
        final Alignment<AminoAcidSequence> alignment;

        AlignmentInfo(boolean sameAntigen, Alignment<AminoAcidSequence> alignment,
                      boolean unknownAntigen) {
            this.sameAntigen = sameAntigen;
            this.alignment = alignment;
            this.unknownAntigen = unknownAntigen;
        }
    }

    static String shortMutationType(MutationType mutationType) {
        switch (mutationType) {
            case Substitution:
                return "S";
            case Insertion:
                return "I";
            case Deletion:
                return "D";
            default:
                throw new IllegalArgumentException();
        }
    }
}
