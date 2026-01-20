/**
 * Unit-level leak-safety checks for CR-compatibility regression components.
 *
 * This is intentionally a single binary so we only link the STAR core once.
 *
 * Compile (recommended via core Makefile target added for this test):
 *   make -C core/legacy/source ASAN=1 compat-unit-tests
 *
 * Run:
 *   ASAN_OPTIONS=detect_leaks=1 LSAN_OPTIONS=detect_leaks=1 \
 *     core/legacy/source/compat_unit_tests
 */

#include "SoloReadBarcode.h"
#include "Transcriptome.h"
#include "Solo.h"
#include "Genome.h"

#include "PackedReadInfo.h"

#include <cstdlib>
#include <iostream>
#include <string>

static int check(bool ok, const std::string &label) {
    if (!ok) {
        std::cerr << "FAIL: " << label << "\n";
        return 1;
    }
    return 0;
}

static void initMinimalSoloParams(Parameters &P, uint64_t cbWLsize) {
    // Enable Solo-ish defaults needed for constructing barcode objects.
    P.pSolo.type = ParametersSolo::CB_UMI_Simple;
    P.pSolo.cbWLyes = true;
    P.pSolo.cbWLsize = cbWLsize;
    P.pSolo.umiL = 12;
}

static int testPackedReadInfo() {
    PackedReadInfo pri;
    pri.init(/*nReads*/ 16, /*wlSize*/ 8, /*umiLength*/ 12);

    pri.set(0, /*cbIdx*/ 1, /*umiPacked*/ 7, /*status*/ 3);
    pri.set(1, /*cbIdx*/ 2, /*umiPacked*/ 11, /*status*/ 0);

    uint32_t cb = pri.getCB(0);
    uint32_t umi = pri.getUMI(0);
    uint8_t st = pri.getStatus(0);
    if (check(cb == 1, "PackedReadInfo: getCB") != 0) return 1;
    if (check(umi == 7, "PackedReadInfo: getUMI") != 0) return 1;
    if (check(st == 3, "PackedReadInfo: getStatus") != 0) return 1;

    pri.setStatus(0, 5);
    if (check(pri.getStatus(0) == 5, "PackedReadInfo: setStatus") != 0) return 1;

    uint32_t cb2 = 0, umi2 = 0;
    uint8_t st2 = 0;
    pri.unpack(pri.pack(3, 12, 1), cb2, umi2, st2);
    if (check(cb2 == 3 && umi2 == 12 && st2 == 1, "PackedReadInfo: pack/unpack") != 0) return 1;

    return 0;
}

static int testSoloReadBarcodeLifecycle() {
    Parameters P;
    initMinimalSoloParams(P, /*cbWLsize*/ 32);

    SoloReadBarcode rb(P);
    // Basic smoke calls: validate paths do not crash and remain bounded.
    rb.addStats(/*cbMatch1*/ -1);

    SoloReadBarcode rb2(P);
    rb2.addStats(1);
    rb.addStats(rb2);

    return 0;
}

static int testTranscriptomeLoadUnload() {
    Parameters P;
    P.quant.yes = true;
    P.quant.gene.yes = true;

    // Use minimal on-disk index shipped with tests.
    P.pGe.gDir = std::string("tests/solo_smoke/ref/star_index");
    P.pGe.sjdbGTFfile = "-";

    Transcriptome tr(P);
    // If constructor succeeded, it loaded at least geneInfo.tab.
    return check(tr.nGe > 0, "Transcriptome: nGe populated");
}

static int testSoloDestructorFreesReadBarSum() {
    Parameters P;
    initMinimalSoloParams(P, /*cbWLsize*/ 3000000);
    // Avoid allocating SoloFeature structures: exercise the CB_samTagOut early-return path
    // which still allocates readBarSum.
    P.pSolo.type = ParametersSolo::CB_samTagOut;

    Transcriptome dummyTr(P);  // quant disabled by default, should return quickly
    Solo solo(/*RAchunk*/ nullptr, P, dummyTr);
    return 0;
}

static int testGenomeChrBinFillLifecycle() {
    Parameters P;
    Genome g(P, P.pGe);

    // Minimal synthetic chromosome layout: chrStart must have nChrReal+1 entries.
    g.nChrReal = 2;
    g.chrStart.clear();
    g.chrStart.push_back(0);
    g.chrStart.push_back(1000);
    g.chrStart.push_back(2000);
    g.genomeChrBinNbases = 100;

    g.chrBinFill();
    int failed = check(g.chrBin != nullptr && g.chrBinN > 0, "Genome: chrBinFill allocates");
    delete[] g.chrBin;
    g.chrBin = nullptr;
    g.chrBinN = 0;
    return failed;
}

int main() {
    int failed = 0;

    failed += testPackedReadInfo();
    failed += testSoloReadBarcodeLifecycle();
    failed += testTranscriptomeLoadUnload();
    failed += testSoloDestructorFreesReadBarSum();
    failed += testGenomeChrBinFillLifecycle();

    if (failed == 0) {
        std::cout << "PASS: CR compat regression component checks\n";
        return 0;
    }
    return 1;
}
