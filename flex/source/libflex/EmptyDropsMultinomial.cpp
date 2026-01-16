#include "EmptyDropsMultinomial.h"
#include "EmptyDropsCRSampler.h"
#include "CRLogProb.h"
#include "SimpleGoodTuring/sgt.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <random>
#include <chrono>
#include <iostream>
#include <fstream>
#include <thread>
#include <atomic>
#include <vector>
#include <cassert>
#include <stdexcept>
using namespace std;

AmbientProfile EmptyDropsMultinomial::computeAmbientProfile(
    const vector<uint32>& ambCount,
    uint32 featuresNumber,
    const vector<uint32>& featDetVec,
    uint32 featDetN
) {
    AmbientProfile profile;
    // Build compact view: only genes with nonzero ambient counts
    vector<uint32> nonzeroIdx;
    nonzeroIdx.reserve(featuresNumber);
    for (uint32 i = 0; i < featuresNumber; i++) {
        if (ambCount[i] > 0) nonzeroIdx.push_back(i);
    }
    profile.featuresDetected = static_cast<uint32>(nonzeroIdx.size());
    
    // Count zero-count genes for Good-Turing zero-bin
    uint32 numZeros = featuresNumber - static_cast<uint32>(nonzeroIdx.size());
    
    if (nonzeroIdx.empty()) {
        profile.ambProfileLogP.resize(featuresNumber, -1e300);
        return profile;
    }
    
    // Count frequencies of ambient gene counts (including zero bin for Good-Turing)
    map<uint32, uint32> ambCountFreq;
    if (numZeros > 0) {
        ambCountFreq[0] = numZeros;  // Add zero-count frequency for Good-Turing
    }
    for (auto idx : nonzeroIdx) {
        ambCountFreq[ambCount[idx]]++;
    }
    uint32 maxFreq = ambCountFreq.rbegin()->first;
    
    // Simple Good-Turing estimation (including zero bin)
    vector<double> ambCountFreqSGT(maxFreq + 1, 0.0);
    double pZeroTotal = 0.0;  // Total probability mass for all zero-count genes combined
    {
        SGT<uint32> sgt;
        // Only add nonzero counts to SGT (don't add the zero bin)
        for (auto& cf : ambCountFreq) {
            if (cf.first > 0) {  // Skip zero bin
                sgt.add(cf.first, cf.second);
            }
        }
        sgt.analyse();
        
        // Estimate probabilities for nonzero counts
        for (uint32 freq = 1; freq <= maxFreq; freq++) {
            sgt.estimate(freq, ambCountFreqSGT[freq]);
        }
        
        // Compute total probability mass for zero-count genes using N1/N formula
        // P0 = N1 / N where N1 = count of genes with freq=1, N = total UMI in ambient
        if (numZeros > 0) {
            uint32 n1 = (ambCountFreq.find(1) != ambCountFreq.end()) ? ambCountFreq[1] : 0;
            uint64_t bigN = 0;
            for (auto& cf : ambCountFreq) {
                if (cf.first > 0) {  // Only count nonzero
                    bigN += static_cast<uint64_t>(cf.first) * static_cast<uint64_t>(cf.second);
                }
            }
            if (bigN > 0 && n1 > 0) {
                pZeroTotal = static_cast<double>(n1) / static_cast<double>(bigN);
            } else {
                pZeroTotal = 1e-10 * numZeros;  // Minimum fallback per gene
            }
        }
    }
    
    // Build ambient profile (probabilities) for all features
    profile.ambProfileLogP.resize(featuresNumber, 0.0);  // Will be set properly below
    double norm1 = 0.0;
    
    // Assign probabilities to nonzero-count genes
    for (auto idx : nonzeroIdx) {
        uint32 c = ambCount[idx];
        double p = (c < ambCountFreqSGT.size()) ? ambCountFreqSGT[c] : 0.0;
        profile.ambProfileLogP[idx] = p; // store probability temporarily
        norm1 += p;
    }
    
    // Add total zero-bin probability to normalization
    norm1 += pZeroTotal;
    
    // Probability per zero-count gene (before normalization)
    double pZeroPerGene = (numZeros > 0) ? (pZeroTotal / static_cast<double>(numZeros)) : 0.0;
    
    // Normalize and convert to log probabilities
    profile.ambProfilePnon0.reserve(nonzeroIdx.size());
    profile.ambProfileLogPnon0.reserve(nonzeroIdx.size());
    if (norm1 > 0.0) {
        // Normalize nonzero genes
        for (auto idx : nonzeroIdx) {
            double cf = profile.ambProfileLogP[idx];
            if (cf > 0) {
                cf /= norm1;
                profile.ambProfilePnon0.push_back(cf);
                double l = log(cf);
                profile.ambProfileLogPnon0.push_back(l);
                profile.ambProfileLogP[idx] = l;
            } else {
                profile.ambProfileLogP[idx] = -1e300;
            }
        }
        
        // Assign normalized probability to zero-count genes
        if (numZeros > 0 && pZeroPerGene > 0) {
            double pZeroNorm = pZeroPerGene / norm1;
            double logPZero = log(pZeroNorm);
            for (uint32 i = 0; i < featuresNumber; i++) {
                if (ambCount[i] == 0) {
                    profile.ambProfileLogP[i] = logPZero;
                }
            }
        } else {
            // Set zero-count genes to very small value
            for (uint32 i = 0; i < featuresNumber; i++) {
                if (ambCount[i] == 0) {
                    profile.ambProfileLogP[i] = -1e300;
                }
            }
        }
    }
    
    return profile;
}

double EmptyDropsMultinomial::logMultinomialPDFsparse(
    const vector<double>& ambProfileLogP,
    const vector<uint32>& countCellGeneUMI,
    uint32 stride,
    uint32 shift,
    int64 start,
    uint32 nGenes,
    const vector<double>& logFactorial
) {
    uint32 sumCount = 0;
    double sumLogFac = 0.0, sumCountLogP = 0.0;
    
    for (uint32 ig = 0; ig < nGenes; ig++) {
        uint64_t idx = start + ig * stride + shift;
        if (idx >= countCellGeneUMI.size()) {
            // Out of bounds access - return invalid value
            return -1e10;
        }
        auto count1 = countCellGeneUMI[idx];
        sumCount += count1;
        
        if (count1 >= logFactorial.size()) {
            // Count exceeds logFactorial size - return invalid value
            return -1e10;
        }
        sumLogFac += logFactorial[count1];
        
        uint64_t geneIdx = start + ig * stride;
        if (geneIdx >= countCellGeneUMI.size()) {
            // Out of bounds access - return invalid value
            return -1e10;
        }
        uint32 geneId = countCellGeneUMI[geneIdx];
        if (geneId >= ambProfileLogP.size()) {
            // Gene index out of bounds for ambient profile - return invalid value
            return -1e10;
        }
        sumCountLogP += ambProfileLogP[geneId] * count1;
    }
    
    if (sumCount >= logFactorial.size()) {
        // Sum count exceeds logFactorial size - return invalid value
        return -1e10;
    }
    
    return logFactorial[sumCount] - sumLogFac + sumCountLogP;
}

vector<double> EmptyDropsMultinomial::computeLogFactorial(uint32 maxCount) {
    vector<double> logFactorial(maxCount + 1);
    logFactorial[0] = 0.0;  // log(0!) = log(1) = 0
    if (maxCount >= 1) {
        logFactorial[1] = 0.0;  // log(1!) = log(1) = 0
        for (uint32 cc = 2; cc <= maxCount; cc++) {
            logFactorial[cc] = logFactorial[cc - 1] + log((double)cc);
        }
    }
    return logFactorial;
}

vector<EmptyDropsResult> EmptyDropsMultinomial::computePValues(
    const AmbientProfile& ambProfile,
    const vector<uint32>& candidateIndices,
    const vector<uint32>& candidateCounts,
    const vector<uint32>& countCellGeneUMI,
    const vector<uint32>& countCellGeneUMIindex,
    const vector<uint32>& nGenePerCB,
    uint32 countMatStride,
    uint32 umiDedupCountIndMain,
    const EmptyDropsParams& params,
    uint32 nSimpleCells,
    const vector<string>& featureNames,
    uint32 nTotalCells,
    const string& debugOutputDir,
    const string& tagName,
    bool enableInvariantChecks
) {
    vector<EmptyDropsResult> results;
    
    if (candidateIndices.empty() || ambProfile.ambProfilePnon0.empty()) {
        return results;
    }
    
    // Use default seed if not specified (matches existing STAR behavior)
    uint64 seed = (params.seed == 0) ? 19760110LLU : params.seed;
    
    // Guard: clamp nSimpleCells to prevent out-of-bounds access
    uint32 nSimple = min(nSimpleCells, static_cast<uint32>(candidateIndices.size()));
    
    // Debug assertion: verify simple cells really sit first (if nSimple > 0)
    #ifndef NDEBUG
    if (nSimple > 0) {
        // In debug builds, we can add additional checks here if needed
        // For now, we rely on the caller to ensure correct ordering
    }
    #endif
    
    // Compute observed log probabilities for each candidate
    // Find maximum count across tail candidates (after simple cells) for logFactorial and simulation sizing
    uint32 maxCount = 0;
    for (uint32 icand = nSimple; icand < candidateCounts.size(); icand++) {
        if (candidateCounts[icand] > maxCount) {
            maxCount = candidateCounts[icand];
        }
    }
    // Also check simple cells for logFactorial sizing (needed for obsLogProb computation)
    for (uint32 icand = 0; icand < nSimple; icand++) {
        if (candidateCounts[icand] > maxCount) {
            maxCount = candidateCounts[icand];
        }
    }
    // Size logFactorial to handle sumCount (sum of gene counts can be up to ~2x UMI count)
    // Add safety margin to prevent out-of-bounds access
    uint32 logFactorialSize = maxCount * 2 + 1000;  // Safety margin for sum of counts
    
    // Timing: logFactorial build
    auto t0 = chrono::high_resolution_clock::now();
    vector<double> logFactorial = computeLogFactorial(logFactorialSize);
    auto t1 = chrono::high_resolution_clock::now();
    auto logFactorialTime = chrono::duration_cast<chrono::milliseconds>(t1 - t0).count();
    
    // The ambient profile and matrix are already compact (built in FlexFilter per compact_copy_plan.md)
    // Gene IDs in countCellGeneUMI are already compact indices
    // ambProfile is already the compact ambient profile
    
    // Ensure ambient consistency: renormalize ambProfilePnon0 at compute time if needed
    vector<double> ambProfilePnon0Norm = ambProfile.ambProfilePnon0;
    vector<double> ambProfileLogPnon0Norm = ambProfile.ambProfileLogPnon0;
    uint32 nNonZeroForSim = static_cast<uint32>(ambProfilePnon0Norm.size());
    double normSum = accumulate(ambProfilePnon0Norm.begin(), ambProfilePnon0Norm.end(), 0.0);
    if (normSum > 0 && fabs(normSum - 1.0) > 1e-10) {
        for (size_t i = 0; i < ambProfilePnon0Norm.size(); i++) {
            ambProfilePnon0Norm[i] /= normSum;
            ambProfileLogPnon0Norm[i] = log(ambProfilePnon0Norm[i]);
        }
    }
    
    cerr << "[EmptyDrops] Ambient profile: " << nNonZeroForSim << " nonzero genes" << endl;
    
    // Compute observed log probabilities (matrix gene IDs are already compact indices)
    vector<double> obsLogProb(candidateIndices.size());
    for (uint32 icand = 0; icand < candidateIndices.size(); icand++) {
        uint32 icell = candidateIndices[icand];
        if (icell >= countCellGeneUMIindex.size()) {
            obsLogProb[icand] = -1e10;
            continue;
        }
        uint64_t startIdx = countCellGeneUMIindex[icell];
        uint32 nGenes = (icell < nGenePerCB.size()) ? nGenePerCB[icell] : 0;
        
        // Extract cell gene IDs (already compact) and counts
        vector<uint32_t> cellGeneIds;
        vector<uint32_t> cellCounts;
        for (uint32 g = 0; g < nGenes; g++) {
            uint64_t geneOffset = startIdx + g * countMatStride;
            if (geneOffset >= countCellGeneUMI.size()) continue;
            uint32 geneId = countCellGeneUMI[geneOffset];
            uint32 count = countCellGeneUMI[geneOffset + umiDedupCountIndMain];
            cellGeneIds.push_back(geneId);
            cellCounts.push_back(count);
        }
        // Compute logprob using compact ambient (ambProfile.ambProfileLogP is compact)
        obsLogProb[icand] = computeDenseMultinomialLogProb(cellCounts, cellGeneIds, ambProfile.ambProfileLogP);
    }
    
    // Initialize results for all candidates
    results.resize(candidateIndices.size());
    for (uint32 icand = 0; icand < candidateIndices.size(); icand++) {
        results[icand].cellIndex = candidateIndices[icand];
        results[icand].obsLogProb = obsLogProb[icand];
    }
    
    // Set simple cells to pass without simulation
    for (uint32 icand = 0; icand < nSimple; icand++) {
        results[icand].monteCarloRank = 0;
        results[icand].pValue = 0.0;
        results[icand].pAdjusted = 0.0;
        results[icand].passesRawP = true;
        results[icand].passesFDR = true;
    }
    
    // Use CR sampler (DropletUtils-style nested multinomial) - always enabled
    if (nSimple < candidateIndices.size() && params.simN > 0) {
        cerr << "[EmptyDrops] Using CR sampler (DropletUtils-style nested multinomial)" << endl;
        
        // Sort tail candidates by total ascending (required for run-length encoding)
        vector<pair<uint32, uint32>> tailCandidates;  // (original index, total)
        for (uint32 icand = nSimple; icand < candidateIndices.size(); icand++) {
            tailCandidates.push_back({icand, candidateCounts[icand]});
        }
        sort(tailCandidates.begin(), tailCandidates.end(), 
             [](const pair<uint32, uint32>& a, const pair<uint32, uint32>& b) {
                 return a.second < b.second;  // Sort by total ascending
             });
        
        // Build run-length encoded data: totalval (distinct totals), totallen (run lengths), prob (logprobs)
        vector<uint32> totalval;
        vector<uint32> totallen;
        vector<double> prob;
        
        // For mapping counts back after sorting by logprob within each total run
        vector<uint32> sortedToOriginal;
        
        if (!tailCandidates.empty()) {
            size_t i = 0;
            while (i < tailCandidates.size()) {
                uint32 currentTotal = tailCandidates[i].second;
                // Compute log(total!) to subtract from observed logprobs
                // This is because CR sampler's simulated logprob doesn't include log(n!) term
                // but computeDenseMultinomialLogProb does. Since log(n!) is constant within a run,
                // we subtract it to make observed and simulated logprobs comparable.
                double logFactTotal = lgamma(static_cast<double>(currentTotal) + 1.0);
                
                // collect run
                vector<pair<double, uint32>> run; // (logprob - logFactTotal, original index)
                size_t j = i;
                while (j < tailCandidates.size() && tailCandidates[j].second == currentTotal) {
                    // Subtract log(n!) from observed logprob for fair comparison with simulated
                    double adjLogProb = obsLogProb[tailCandidates[j].first] - logFactTotal;
                    run.push_back({adjLogProb, tailCandidates[j].first});
                    j++;
                }
                // sort run by logprob ascending (required by montecarlo_pval lower_bound)
                sort(run.begin(), run.end(), [](const pair<double,uint32>& a, const pair<double,uint32>& b) {
                    return a.first < b.first;
                });
                totalval.push_back(currentTotal);
                totallen.push_back(static_cast<uint32>(run.size()));
                for (const auto& r : run) {
                    prob.push_back(r.first);
                    sortedToOriginal.push_back(r.second);
                }
                i = j;
            }
        }
        
        // Build ambient vector for CR sampler from FULL compact ambient profile
        // Must use ambProfileLogP (all compact genes) to match observed logprob calculation
        // ambProfilePnon0 only has nonzero-ambient genes, which creates a dimension mismatch
        vector<double> ambientForSampler;
        ambientForSampler.reserve(ambProfile.ambProfileLogP.size());
        for (size_t i = 0; i < ambProfile.ambProfileLogP.size(); i++) {
            // Convert log prob back to prob; ambProfileLogP includes tiny values for zero-count genes
            double logp = ambProfile.ambProfileLogP[i];
            double p = (logp > -500.0) ? exp(logp) : 0.0;  // Avoid underflow
            ambientForSampler.push_back(p);
        }
        // Renormalize to sum to 1 (should already be close)
        double sumP = accumulate(ambientForSampler.begin(), ambientForSampler.end(), 0.0);
        if (sumP > 0 && fabs(sumP - 1.0) > 1e-10) {
            for (auto& p : ambientForSampler) p /= sumP;
        }
        cerr << "[EmptyDrops] CR sampler using full compact ambient (" << ambientForSampler.size() << " genes, sum=" << sumP << ")" << endl;
        
        // Size check: validate invariants before calling CR sampler (only when explicitly enabled)
        if (enableInvariantChecks) {
            uint32_t sumTotallen = 0;
            for (uint32_t len : totallen) {
                sumTotallen += len;
            }
            if (sumTotallen != prob.size()) {
                string errorMsg = (tagName.empty() ? "" : "tag=" + tagName + " ") + 
                                 "sum(totallen)=" + to_string(sumTotallen) + " != prob.size()=" + to_string(prob.size());
                cerr << "[ED_INVARIANT_FAIL] " << errorMsg << endl;
                throw runtime_error("[ED_INVARIANT_FAIL] " + errorMsg);
            }
            if (ambientForSampler.empty()) {
                string errorMsg = (tagName.empty() ? "" : "tag=" + tagName + " ") + "ambientForSampler.size()=0";
                cerr << "[ED_INVARIANT_FAIL] " << errorMsg << endl;
                throw runtime_error("[ED_INVARIANT_FAIL] " + errorMsg);
            }
        }
        
        // Call CR sampler with ambient vector (compact or dense)
        // Use mcThreads for parallel Monte Carlo (0 = single-threaded)
        vector<uint32_t> above = EmptyDropsCRSampler::montecarloPval(
            totalval, totallen, prob, ambientForSampler, params.simN, seed, params.mcThreads);
        
        // Extract p-values from above vector (which is in sorted order within each total run)
        for (size_t idx = 0; idx < sortedToOriginal.size(); idx++) {
            uint32 origIdx = sortedToOriginal[idx];
            uint32 count = above[idx];
            results[origIdx].monteCarloRank = count;
            results[origIdx].pValue = (double)(count + 1) / (double)(params.simN + 1);
            results[origIdx].passesRawP = (results[origIdx].pValue <= params.rawPvalueThreshold);
        }
        
        // FDR and final processing
        for (auto& res : results) {
            res.pAdjusted = res.pValue;
            res.passesFDR = (res.pAdjusted <= params.FDR);
        }
        
        return results;
    }
    
    // Edge case: no tail candidates or simN is 0 - set all tail candidates to p=1.0
    for (uint32 icand = nSimple; icand < candidateIndices.size(); icand++) {
        results[icand].monteCarloRank = 0;
        results[icand].pValue = 1.0;
        results[icand].pAdjusted = 1.0;
        results[icand].passesRawP = false;
        results[icand].passesFDR = false;
    }
    
    // FDR and final processing
    for (auto& res : results) {
        res.pAdjusted = res.pValue;
        res.passesFDR = (res.pAdjusted <= params.FDR);
    }
    
    return results;
}
