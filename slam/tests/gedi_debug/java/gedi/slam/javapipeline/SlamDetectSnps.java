/*
 * GEDI instrumentation patch for debugging SNP detection.
 *
 * This file is based on decompilation of GEDI 1.0.6d's
 * gedi.slam.javapipeline.SlamDetectSnps (CFR 0.152), with minimal edits:
 * - Optional debug logging for a target locus via env vars.
 *
 * Env vars:
 *   GEDI_SNP_DEBUG_LOC   e.g. "1:26429152" or "chr1:26429152"
 *   GEDI_SNP_DEBUG_OUT   output TSV path (default: gedi_snp_debug.tsv)
 *   GEDI_SNP_DEBUG_MAX   max debug lines (default: 5000)
 *
 * Debug log TSV columns:
 *   event, chr, pos, geneRef, readRef, readRegion, distinct, count_c,
 *   mpos, mappedPos, overlap, gBase, rBase, covAtPos, mmAtPos, pval
 *
 * NOTE: This class is intended to be placed first on the Java classpath to
 * override GEDI's original class (jar stays unchanged).
 */
package gedi.slam.javapipeline;

import cern.colt.bitvector.BitVector;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.FunctorUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.functions.EI;
import gedi.util.genomic.CoverageAlgorithm;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Pattern;
import jdistlib.Beta;
import jdistlib.Binomial;

public class SlamDetectSnps extends GediProgram {
    // -------------------------
    // Debug plumbing
    // -------------------------
    private static final Object DEBUG_LOCK = new Object();
    private static volatile LineWriter DEBUG_WRITER = null;
    private static volatile int DEBUG_MAX = 5000;
    private static final AtomicInteger DEBUG_LINES = new AtomicInteger(0);
    private static volatile String DEBUG_CHR = null; // empty string => disabled
    private static volatile int DEBUG_POS = -1; // 1-based, GEDI-style

    private static void initDebugOnce() {
        if (DEBUG_CHR != null) return;
        String loc = System.getenv("GEDI_SNP_DEBUG_LOC");
        if (loc == null || loc.trim().isEmpty()) {
            DEBUG_CHR = "";
            return;
        }
        loc = loc.trim();
        if (loc.startsWith("chr")) loc = loc.substring(3);
        String[] parts = loc.split(":", 2);
        if (parts.length != 2) {
            DEBUG_CHR = "";
            return;
        }
        DEBUG_CHR = parts[0];
        try {
            DEBUG_POS = Integer.parseInt(parts[1]);
        } catch (NumberFormatException e) {
            DEBUG_CHR = "";
            DEBUG_POS = -1;
            return;
        }
        String maxStr = System.getenv("GEDI_SNP_DEBUG_MAX");
        if (maxStr != null && !maxStr.trim().isEmpty()) {
            try {
                DEBUG_MAX = Integer.parseInt(maxStr.trim());
            } catch (NumberFormatException e) {
                // ignore
            }
        }
        String out = System.getenv("GEDI_SNP_DEBUG_OUT");
        if (out == null || out.trim().isEmpty()) out = "gedi_snp_debug.tsv";
        try {
            DEBUG_WRITER = new LineOrientedFile(out).write();
            DEBUG_WRITER.writeLine(
                    "event\tchr\tpos\tgeneRef\treadRef\treadRegion\tdistinct\tcount_c\tmpos\tmappedPos\toverlap\tgBase\trBase\tcovAtPos\tmmAtPos\tpval");
        } catch (Exception e) {
            DEBUG_WRITER = null;
        }
    }

    private static boolean debugEnabled() {
        initDebugOnce();
        return DEBUG_WRITER != null && DEBUG_CHR != null && !DEBUG_CHR.isEmpty() && DEBUG_POS > 0;
    }

    private static boolean isDebugPos(ImmutableReferenceGenomicRegion<String> gene, int pos) {
        if (!debugEnabled()) return false;
        String chr = gene.getReference().toStrandIndependent().getName();
        return DEBUG_CHR.equals(chr) && DEBUG_POS == pos;
    }

    private static void debugLogLine(String line) {
        if (!debugEnabled()) return;
        int n = DEBUG_LINES.incrementAndGet();
        if (n > DEBUG_MAX) return;
        synchronized (DEBUG_LOCK) {
            try {
                DEBUG_WRITER.writeLine(line);
            } catch (Exception e) {
                // ignore
            }
        }
    }

    private static void closeDebug() {
        if (!debugEnabled()) return;
        synchronized (DEBUG_LOCK) {
            try {
                DEBUG_WRITER.close();
            } catch (Exception e) {
                // ignore
            }
        }
    }

    // -------------------------
    // GEDI original code
    // -------------------------
    public SlamDetectSnps(SlamParameterSet params) {
        this.addInput(params.nthreads);
        this.addInput(params.genomic);
        this.addInput(params.reads);
        this.addInput(params.snpConv);
        this.addInput(params.snpPval);
        this.addInput(params.strandness);
        this.addInput(params.no4sUpattern);
        this.addInput(params.newsnp);
        this.addInput(params.prefix);
        this.addOutput(params.snpFile);
        this.addOutput(params.strandnessFile);
    }

    public String execute(GediProgramContext context) throws IOException {
        int nthreads = this.getIntParameter(0);
        Genomic genomic = (Genomic)this.getParameter(1);
        GenomicRegionStorage reads = (GenomicRegionStorage)this.getParameter(2);
        double conv = this.getDoubleParameter(3);
        double pvalCutoff = this.getDoubleParameter(4);
        Strandness strandness = (Strandness)this.getParameter(5);
        String pat = (String)this.getParameter(6);
        boolean newsnp = this.getBooleanParameter(7);

        if (debugEnabled()) {
            context.getLog().info("GEDI SNP debug enabled for " + DEBUG_CHR + ":" + DEBUG_POS + " (max=" + DEBUG_MAX + ", out=" + System.getenv("GEDI_SNP_DEBUG_OUT") + ")");
        }

        Pattern no4sUPattern = Pattern.compile(pat, 2);
        String[] cond = (String[])reads.getMetaDataConditions();
        int[] no4sUIndices = EI.along(cond).filterInt(i -> no4sUPattern.matcher(cond[i]).find()).toIntArray();
        boolean[] no4sU = new boolean[cond.length];
        for (int i : no4sUIndices) {
            no4sU[i] = true;
        }
        if (newsnp) {
            context.getLog().info("Finding SNPs (only using no4sU)...");
        } else {
            context.getLog().info("Finding SNPs...");
        }
        AtomicInteger senseCounter = new AtomicInteger();
        AtomicInteger antisenseCounter = new AtomicInteger();
        genomic.getGenes().ei()
                .progress(context.getProgress(), (int)genomic.getGenes().size(), r -> "In " + (String)r.getData())
                .parallelized(nthreads, 5, ei -> ei.map(l -> newsnp
                        ? SlamDetectSnps.collectNew(pvalCutoff, reads, no4sU, (ImmutableReferenceGenomicRegion<String>)l, senseCounter, antisenseCounter)
                        : SlamDetectSnps.collect(conv, pvalCutoff, reads, (ImmutableReferenceGenomicRegion<String>)l, senseCounter, antisenseCounter)))
                .filter(s -> s.length() > 0)
                .print("Location\tCoverage\tMismatches\tP value", this.getOutputFile(0).getPath());

        context.getLog().info("Auto-Detecting sequencing mode: Sense:" + senseCounter.get() + " Antisense:" + antisenseCounter.get());
        if (senseCounter.get() > antisenseCounter.get() * 2) {
            context.getLog().info("Detected strand-specific sequencing (Sense)");
            if (strandness.equals(Strandness.AutoDetect)) {
                strandness = Strandness.Sense;
            } else {
                context.getLog().info("Overriden by command line: " + strandness.name());
            }
        } else if (antisenseCounter.get() > senseCounter.get() * 2) {
            context.getLog().info("Detected strand-specific sequencing (Antisense)");
            if (strandness.equals(Strandness.AutoDetect)) {
                strandness = Strandness.Antisense;
            } else {
                context.getLog().info("Overriden by command line: " + strandness.name());
            }
        } else {
            context.getLog().info("Detected strand-unspecific sequencing");
            if (strandness.equals(Strandness.AutoDetect)) {
                strandness = Strandness.Unspecific;
            } else {
                context.getLog().info("Overriden by command line: " + strandness.name());
            }
        }
        LineWriter sout = this.getOutputWriter(1);
        sout.writeLine(strandness.name());
        sout.close();

        closeDebug();
        return null;
    }

    public static String collect(double conv, double pvalCutoff, GenomicRegionStorage<AlignedReadsData> reads,
                                 ImmutableReferenceGenomicRegion<String> gene,
                                 AtomicInteger senseCounter, AtomicInteger antisenseCounter) {
        gene = gene.toMutable().transformRegion(r -> r.extendBack(1000).extendFront(1000)).toImmutable();
        CoverageAlgorithm cov = new CoverageAlgorithm(gene).setExample(NumericArray.wrap(0.0));
        int[] co = new int[2];
        FunctorUtils.ChainedIterator rit = reads.ei((ReferenceGenomicRegion)gene).sideEffect(r -> co[0] = co[0] + 1)
                .chain((Iterator)reads.ei((ReferenceGenomicRegion)gene.toMutable().toOppositeStrand()).sideEffect(r -> co[1] = co[1] + 1));

        BitVector counted = new BitVector(300);
        HashMap<Integer, double[]> counter = new HashMap<>();

        for (Object r2o : rit.loop()) {
            @SuppressWarnings("unchecked")
            ImmutableReferenceGenomicRegion<AlignedReadsData> r2 = (ImmutableReferenceGenomicRegion<AlignedReadsData>)r2o;
            if (r2.getRegion().getTotalLength() > counted.size()) {
                counted = new BitVector(r2.getRegion().getTotalLength());
            }
            AlignedReadsData ard = r2.getData();
            for (int d = 0; d < ard.getDistinctSequences(); ++d) {
                counted.clear();
                double c = ard.getTotalCountForDistinct(d, ReadCountMode.All);
                boolean debugThisReadCoversPos = debugEnabled()
                        && DEBUG_CHR.equals(gene.getReference().toStrandIndependent().getName())
                        && r2.getRegion().contains(DEBUG_POS);
                boolean debugSawMismatchAtPos = false;
                boolean debugCountedMismatchAtPos = false;

                for (int v = 0; v < ard.getVariationCount(d); ++v) {
                    if (!ard.isMismatch(d, v)) continue;
                    int mpos = ard.getMismatchPos(d, v);
                    int pos = r2.map(mpos);
                    boolean overlap = ard.isPositionInOverlap(d, mpos);

                    if (debugThisReadCoversPos && pos == DEBUG_POS) {
                        debugSawMismatchAtPos = true;
                    }

                    // Log exact locus mismatches
                    if (isDebugPos(gene, pos)) {
                        String gBase = String.valueOf(ard.getMismatchGenomic(d, v));
                        String rBase = String.valueOf(ard.getMismatchRead(d, v));
                        debugLogLine("mismatch\t" + DEBUG_CHR + "\t" + DEBUG_POS + "\t" +
                                gene.getReference().toString() + "\t" +
                                r2.getReference().toString() + "\t" +
                                r2.getRegion().toString() + "\t" +
                                d + "\t" + String.format("%.3f", c) + "\t" +
                                mpos + "\t" + pos + "\t" + (overlap ? 1 : 0) + "\t" +
                                gBase + "\t" + rBase + "\tNA\tNA\tNA");
                    }
                    // Also log near-by mismatches (helps detect systematic coordinate shifts)
                    if (debugThisReadCoversPos && Math.abs(pos - DEBUG_POS) <= 2) {
                        String gBase = String.valueOf(ard.getMismatchGenomic(d, v));
                        String rBase = String.valueOf(ard.getMismatchRead(d, v));
                        debugLogLine("near\t" + DEBUG_CHR + "\t" + DEBUG_POS + "\t" +
                                gene.getReference().toString() + "\t" +
                                r2.getReference().toString() + "\t" +
                                r2.getRegion().toString() + "\t" +
                                d + "\t" + String.format("%.3f", c) + "\t" +
                                mpos + "\t" + pos + "\t" + (overlap ? 1 : 0) + "\t" +
                                gBase + "\t" + rBase + "\tNA\tNA\tNA");
                    }

                    if (gene.getRegion().contains(pos) && (!overlap || counted.getQuick(mpos))) {
                        double[] arr = counter.computeIfAbsent(pos, x -> new double[2]);
                        arr[1] = arr[1] + c;
                        if (debugThisReadCoversPos && pos == DEBUG_POS) {
                            debugCountedMismatchAtPos = true;
                        }
                        if (isDebugPos(gene, pos)) {
                            debugLogLine("add_mm\t" + DEBUG_CHR + "\t" + DEBUG_POS + "\t" +
                                    gene.getReference().toString() + "\t" +
                                    r2.getReference().toString() + "\t" +
                                    r2.getRegion().toString() + "\t" +
                                    d + "\t" + String.format("%.3f", c) + "\t" +
                                    mpos + "\t" + pos + "\t" + (overlap ? 1 : 0) + "\t" +
                                    "NA\tNA\tNA\t" + String.format("%.3f", arr[1]) + "\tNA");
                        }
                    }
                    counted.putQuick(mpos, true);
                }

                if (debugThisReadCoversPos) {
                    debugLogLine("distinct_cov\t" + DEBUG_CHR + "\t" + DEBUG_POS + "\t" +
                            gene.getReference().toString() + "\t" +
                            r2.getReference().toString() + "\t" +
                            r2.getRegion().toString() + "\t" +
                            d + "\t" + String.format("%.3f", c) + "\t" +
                            "-1\t" + DEBUG_POS + "\t0\t" +
                            (debugSawMismatchAtPos ? "SAW_MISMATCH" : "NO_MISMATCH") + "\t" +
                            (debugCountedMismatchAtPos ? "COUNTED_MISMATCH" : "NOT_COUNTED") + "\tNA\tNA\tNA");
                }
                cov.add(r2.getRegion(), NumericArray.wrap(c));
            }
        }
        senseCounter.addAndGet(co[0]);
        antisenseCounter.addAndGet(co[1]);

        StringBuilder sb = new StringBuilder();
        Iterator it = counter.keySet().iterator();
        while (it.hasNext()) {
            int l = (Integer)it.next();
            double[] count = counter.get(l);
            double covAt = cov.getCoverages(gene.induce(l)).getDouble(0);
            double pval = Beta.cumulative(conv, count[1] + conv, covAt - count[1] + 1.0, true, false);

            if (isDebugPos(gene, l)) {
                debugLogLine("pval\t" + DEBUG_CHR + "\t" + DEBUG_POS + "\t" +
                        gene.getReference().toString() + "\tNA\tNA\t-1\tNA\t-1\t" + l + "\tNA\tNA\tNA\t" +
                        String.format("%.3f", covAt) + "\t" + String.format("%.3f", count[1]) + "\t" + String.format("%.6g", pval));
            }

            if (pval < pvalCutoff) {
                sb.append(gene.getReference().toStrandIndependent()).append(":").append(l)
                        .append("\t").append(String.format("%.1f", count[1]))
                        .append("\t").append(String.format("%.1f", covAt))
                        .append("\t").append(String.format("%.3g", pval)).append("\n");
                continue;
            }
            it.remove();
        }
        if (sb.length() > 0) {
            sb.delete(sb.length() - 1, sb.length());
        }
        return sb.toString();
    }

    public static String collectNew(double pvalCutoff, GenomicRegionStorage<AlignedReadsData> reads, boolean[] no4sU,
                                    ImmutableReferenceGenomicRegion<String> gene,
                                    AtomicInteger senseCounter, AtomicInteger antisenseCounter) {
        gene = gene.toMutable().transformRegion(r -> r.extendBack(1000).extendFront(1000)).toImmutable();
        CoverageAlgorithm cov = new CoverageAlgorithm(gene).setExample(NumericArray.wrap(0.0));
        int[] co = new int[2];
        FunctorUtils.ChainedIterator rit = reads.ei((ReferenceGenomicRegion)gene).sideEffect(r -> co[0] = co[0] + 1)
                .chain((Iterator)reads.ei((ReferenceGenomicRegion)gene.toMutable().toOppositeStrand()).sideEffect(r -> co[1] = co[1] + 1));
        BitVector counted = new BitVector(300);
        HashMap<Integer, double[]> counter = new HashMap<>();
        for (Object r2o : rit.loop()) {
            @SuppressWarnings("unchecked")
            ImmutableReferenceGenomicRegion<AlignedReadsData> r2 = (ImmutableReferenceGenomicRegion<AlignedReadsData>)r2o;
            if (r2.getRegion().getTotalLength() > counted.size()) {
                counted = new BitVector(r2.getRegion().getTotalLength());
            }
            AlignedReadsData ard = r2.getData();
            for (int d = 0; d < ard.getDistinctSequences(); ++d) {
                counted.clear();
                double c = ard.getTotalCountForDistinct(d, ReadCountMode.All, (cond, count) -> no4sU[cond] ? count : 0.0);
                for (int v = 0; v < ard.getVariationCount(d); ++v) {
                    if (!ard.isMismatch(d, v)) continue;
                    int mpos = ard.getMismatchPos(d, v);
                    int pos = r2.map(mpos);
                    boolean overlap = ard.isPositionInOverlap(d, mpos);
                    if (gene.getRegion().contains(pos) && (!overlap || counted.getQuick(mpos))) {
                        double[] arr = counter.computeIfAbsent(pos, x -> new double[2]);
                        arr[1] = arr[1] + c;
                    }
                    counted.putQuick(mpos, true);
                }
                cov.add(r2.getRegion(), NumericArray.wrap(c));
            }
        }
        senseCounter.addAndGet(co[0]);
        antisenseCounter.addAndGet(co[1]);
        StringBuilder sb = new StringBuilder();
        Iterator it = counter.keySet().iterator();
        while (it.hasNext()) {
            int l = (Integer)it.next();
            double[] count2 = counter.get(l);
            int convcount = (int)count2[1];
            int covcount = (int)cov.getCoverages(gene.induce(l)).getDouble(0);
            double pval = Binomial.cumulative((double)(convcount - 1), (double)covcount, 0.001, false, false);
            if (isDebugPos(gene, l)) {
                debugLogLine("pval_new\t" + DEBUG_CHR + "\t" + DEBUG_POS + "\t" +
                        gene.getReference().toString() + "\tNA\tNA\t-1\tNA\t-1\t" + l + "\tNA\tNA\tNA\t" +
                        covcount + "\t" + convcount + "\t" + String.format("%.6g", pval));
            }
            if (pval < pvalCutoff) {
                sb.append(gene.getReference().toStrandIndependent()).append(":").append(l)
                        .append("\t").append(String.format("%.1f", count2[1]))
                        .append("\t").append(String.format("%.1f", covcount))
                        .append("\t").append(String.format("%.3g", pval)).append("\n");
                continue;
            }
            it.remove();
        }
        if (sb.length() > 0) {
            sb.delete(sb.length() - 1, sb.length());
        }
        return sb.toString();
    }

    public static void main(String[] args) {
        // Keep GEDI's original self-test main (minimal).
        ImmutableReferenceGenomicRegion gene = ImmutableReferenceGenomicRegion.parse("1+:1000-2000", "Test");
        MemoryIntervalTreeStorage reads = new MemoryIntervalTreeStorage(DefaultAlignedReadsData.class);
        AlignedReadsDataFactory fac = new AlignedReadsDataFactory(1);
        fac.start().newDistinctSequence().addMismatch(4, 'T', 'C', false).setCount(0, 1);
        reads.add((ReferenceGenomicRegion)ImmutableReferenceGenomicRegion.parse("1+:1100-1110", fac.create()));
        System.out.println("OK");
    }
}

