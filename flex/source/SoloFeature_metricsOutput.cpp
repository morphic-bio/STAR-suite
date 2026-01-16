#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "ErrorWarning.h"
#include "streamFuns.h"
#include "GlobalVariables.h"
#include <fstream>
#include <iomanip>

void SoloFeature::writeUMICleanupMetrics() {
    // Only write if UMI correction is enabled
    if (pSolo.umiCorrectionMode == 0 && cellsAllowSet.empty() && assignmentsMap.empty()) {
        return;
    }
    
    // Write metrics even when writing keys.bin (for keys-to-MEX replay workflow)
    // Metrics help validate that clique correction ran correctly
    
    string metricsPath = outputPrefix + "umi_cleanup_metrics.tsv";
    ofstream fp(metricsPath);
    if (!fp.is_open()) {
        P.inOut->logMain << "Warning: Could not open UMI cleanup metrics file for writing: " << metricsPath << endl;
        return;
    }
    
    // Write header (matching bam_to_counts format exactly)
    fp << "phase\tcells_input\tcells_kept\tambiguous_cells\tmolecules_before\tmolecules_after\tumis_before\tumis_after\treads_before\treads_after\tmerges_total\tmax_component\tcomponents_capped\tcomponents_below_threshold\treads_ur_grouped\treads_ur_missing\tcomp_size_1\tcomp_size_2\tcomp_size_3\tcomp_size_4\tcomp_size_gt4\n";
    
    if (pSolo.umiCorrectionMode == 0) {
        // Filter-only phase
        fp << "filter_only\t" << cellsInput << "\t" << cellsKept << "\t" << ambiguousCellsDropped
           << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t"
           << readsURGrouped << "\t" << readsURMissing
           << "\t0\t0\t0\t0\t0\n";
    } else {
        // Filter-only phase (before correction)
        fp << "filter_only\t" << cellsInput << "\t" << cellsKept << "\t" << ambiguousCellsDropped
           << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t"
           << readsURGrouped << "\t" << readsURMissing
           << "\t0\t0\t0\t0\t0\n";
        
        // Correction phase
        fp << "correction_clique\t" << cellsInput << "\t" << cellsKept << "\t" << ambiguousCellsDropped
           << "\t0\t0\t"  // molecules_before, molecules_after (not tracked)
           << umisBeforeTotal << "\t" << umisAfterTotal
           << "\t" << readsBeforeTotal << "\t" << readsAfterTotal
           << "\t" << mergesTotal << "\t" << maxComponentSeen
           << "\t" << componentsCappedTotal << "\t" << componentsBelowThresholdTotal
           << "\t" << readsURGrouped << "\t" << readsURMissing
           << "\t" << componentSizeHist[0] << "\t" << componentSizeHist[1] 
           << "\t" << componentSizeHist[2] << "\t" << componentSizeHist[3] 
           << "\t" << componentSizeHist[4] << "\n";
    }
    
    fp.close();
    P.inOut->logMain << "UMI cleanup metrics written to: " << metricsPath << endl;
}

void SoloFeature::writeQCSummary() {
    // Write aggregated QC summary (simplified version matching bam_to_counts format)
    string qcPath = outputPrefix + "qc_summary.tsv";
    ofstream fp(qcPath);
    if (!fp.is_open()) {
        P.inOut->logMain << "Warning: Could not open QC summary file for writing: " << qcPath << endl;
        return;
    }
    
    // Write header
    fp << "Metric\tValue\tPercentage\tDescription\n";
    
    // Basic counts (using GlobalVariables for total reads, SoloFeature stats for mapped)
    uint64_t totalReads = g_statsAll.readN;
    fp << "Total Reads\t" << totalReads << "\t100.0%\tTotal reads processed\n";
    
    uint64_t totalAlignments = readFeatSum->stats.numMappedToTranscriptome();
    double alignPct = (totalReads > 0) ? 100.0 * totalAlignments / totalReads : 0.0;
    fp << "Total Alignments\t" << totalAlignments << "\t" << std::fixed << std::setprecision(1) 
       << alignPct << "%\tTotal alignments processed\n";
    
    uint64_t alignmentsKept = readFeatSum->stats.V[readFeatSum->stats.yesWLmatch];
    double keptPct = (totalAlignments > 0) ? 100.0 * alignmentsKept / totalAlignments : 0.0;
    fp << "Alignments Kept\t" << alignmentsKept << "\t" << keptPct 
       << "%\tAlignments passing all filters\n";
    
    // UMI correction metrics (if enabled)
    if (pSolo.umiCorrectionMode > 0) {
        fp << "UMIs Before Correction\t" << umisBeforeTotal << "\t100.0%\tUMIs before clique correction\n";
        if (umisBeforeTotal > 0) {
            double afterPct = 100.0 * umisAfterTotal / umisBeforeTotal;
            fp << "UMIs After Correction\t" << umisAfterTotal << "\t" << afterPct 
               << "%\tUMIs after clique correction\n";
            fp << "UMI Merges\t" << mergesTotal << "\t" 
               << std::setprecision(1) << (100.0 * mergesTotal / umisBeforeTotal) 
               << "%\tUMIs merged during correction\n";
        }
        fp << "Components Found\t" << componentsTotal << "\t-\tConnected components found\n";
        fp << "Components Capped\t" << componentsCappedTotal << "\t-\tComponents exceeding max size\n";
        fp << "Components Below Threshold\t" << componentsBelowThresholdTotal 
           << "\t-\tComponents failing ratio threshold\n";
    }
    
    fp.close();
    P.inOut->logMain << "QC summary written to: " << qcPath << endl;
}

