#include "SoloFeature.h"
#include "serviceFuns.cpp"
#include "SimpleGoodTuring/sgt.h"
#include <fstream>
#include <math.h>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <random>
#include <sys/stat.h>

double logMultinomialPDFsparse(const vector<double> &ambProfileLogP, const vector<uint32> &countCellGeneUMI, const uint32 stride, const uint32 shift, const int64 start, const uint32 nGenes, const vector<double> &logFactorial);

namespace {
struct EmptyDropsIndexCount {
    uint32 index;
    uint32 count;
};
}

void SoloFeature::emptyDrops_CR()
{
    if (nCB<=pSolo.cellFilter.eDcr.indMin) {
        P.inOut->logMain << "emptyDrops_CR filtering: total number of cells: nCB=" << nCB <<" is smaller than emptyCellMinIndex="<< pSolo.cellFilter.eDcr.indMin
                         << ", which is the starting index for the *true empty* cells. The additional non-empty cells will not be detected.\n";
        return;
    };
    
    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) <<" ... starting emptyDrops_CR filtering" <<endl;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //find genes that were detected in all cells
    unordered_set<uint32> featDet;   
    for (uint32 icb=0; icb<nCB; icb++) {           
        for (uint32 ig=0; ig<nGenePerCB[icb]; ig++) {
            uint32 irec=countCellGeneUMIindex[icb]+ig*countMatStride;
            if (countCellGeneUMI[irec + pSolo.umiDedup.countInd.main] > 0)
                featDet.insert(countCellGeneUMI[irec]); //gene is present if it's count > 0 for 
        };
    };
    uint32 featDetN=featDet.size(); //total number of detected genes - this should have been done already?
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //indCount - total UMI per cell sorted descending
    vector<EmptyDropsIndexCount> indCount(nCB);
    for (uint32 ii=0; ii<nCB; ii++) {
        indCount[ii].index=ii;
        indCount[ii].count=nUMIperCB[ii];
    };
    std::sort(indCount.begin(), indCount.end(), [](const EmptyDropsIndexCount &ic1, const EmptyDropsIndexCount &ic2) {
                                                    return (ic1.count>ic2.count) || (ic1.count==ic2.count && ic1.index<ic2.index); //descending order by count, ascending by index
                                                });
    
    //////////////////////////////////////////////////////////////////////////////////////////
    //ambient gene counts: sum gene expression over the collection of empty cells
    vector<uint32> ambCount(featuresNumber,0);
    for (auto icb=pSolo.cellFilter.eDcr.indMin; icb<min(nCB,pSolo.cellFilter.eDcr.indMax); icb++) {
        auto icb1 = indCount[icb].index;
        for (uint32 ig=0; ig<nGenePerCB[icb1]; ig++) {
            auto irec = countCellGeneUMIindex[icb1]+ig*countMatStride;
            ambCount[countCellGeneUMI[irec+0]] += countCellGeneUMI[irec + pSolo.umiDedup.countInd.main];
        };
    };    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) <<" ... finished ambient cells counting" <<endl;
    
    ////////////////////////////////////////////////////////////////////
    //frequencies
    map<uint32,uint32> ambCountFreq; //ordered map is not really needed
    for (auto &ac: ambCount) {
        ambCountFreq[ac]++;
    };
    if (ambCountFreq.size()<=1) {//only 0-frequency genes are in the empty cells. This is possible because nCB can contain ome cells with no genes - because of multigene
        P.inOut->logMain << "emptyDrops_CR filtering: empty cells contain no genes\n";
        return;
    };
    ambCountFreq[0] -= (featuresNumber-featDetN); //subtract genes that were not detected in *any* cells
    uint32 maxFreq = ambCountFreq.rbegin()->first;
    
    ///////////////////////////////////////////////////////////////////////
    //SGT
    vector<double> ambCountFreqSGT(maxFreq+1);//up to max frequency
    {//SGT estimate of ambient profile
        SGT<uint32> sgt;
        for (auto &cf: ambCountFreq) {
            if (cf.first != 0)
                sgt.add(cf.first, cf.second);
        };
        sgt.analyse();
        
        for (uint32 freq=0; freq<=maxFreq; freq++) {
            sgt.estimate(freq, ambCountFreqSGT[freq]);
        };
        ambCountFreqSGT[0] /= ambCountFreq[0]; //divide freq=0 probability equally among all undetected genes in ambient profile
    };
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) <<" ... finished SGT"<<endl;
    
    //ambient profile for all features
    vector<double> ambProfileLogP(featuresNumber, 0.0);//logarithm
    vector<double> ambProfilePnon0, ambProfileLogPnon0;//only non-0 genes
    {
        for (uint32 ig=0; ig<featuresNumber; ig++) {
            if (featDet.count(ig)>0) {//this is only needed if normalization below is performed
                ambProfileLogP[ig]=ambCountFreqSGT[ambCount[ig]];
            };
        };
        
        double norm1 = accumulate(ambProfileLogP.begin(), ambProfileLogP.end(), 0.0);
        ambProfileLogPnon0.reserve(ambProfileLogP.size());
        ambProfilePnon0.reserve(ambProfileLogP.size());
        for (auto &cf: ambProfileLogP) {
            if (cf>0) {
                cf /= norm1;//normalization is just in case
                ambProfilePnon0.push_back(cf);
                cf = std::log(cf);
                ambProfileLogPnon0.push_back(cf);
            };
        };
    };
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) <<" ... finished ambient profile"<<endl;
    
    //select candidate cells
    uint32 iCandFirst, iCandLast; //first/last candidate cell in the descending sorted indCount
    uint32 minUMI = 0;
    {
        iCandFirst=filteredCells.nCellsSimple;//candidates start right after the cutoff for the simple filtering
        minUMI = int(pSolo.cellFilter.eDcr.umiMinFracMedian * nUMIperCBsorted[filteredCells.nCellsSimple/2]);//this is not exactly median
        minUMI = max(pSolo.cellFilter.eDcr.umiMin, minUMI);
        for (iCandLast=iCandFirst; iCandLast<iCandFirst+pSolo.cellFilter.eDcr.candMaxN; iCandLast++) {
            if (indCount[iCandLast].count<minUMI)
                break;
        };
        --iCandLast;
        
        time(&rawTime);
        P.inOut->logMain << timeMonthDayTime(rawTime) << " ... candidate cells: minUMI="<< minUMI << "; number of candidate cells=" << iCandLast-iCandFirst+1 <<endl;
        if (iCandLast<iCandFirst)
            return; //no candidate cells to consider
    };
    
    //calculate observed probability for each candidate
    vector<double> obsLogProb(iCandLast-iCandFirst+1);
    {
        vector<double> logFactorial; //tabulate log-factorial
        logFactorial.resize(indCount[iCandFirst].count+1);
        logFactorial[1]=0;
        for (uint32 cc=2; cc<logFactorial.size(); cc++)
            logFactorial[cc]=logFactorial[cc-1]+std::log(cc);
        
        for (uint32 icand=0; icand<obsLogProb.size(); icand++) {
            auto icell=indCount[icand+iCandFirst].index;
            obsLogProb[icand]=logMultinomialPDFsparse(ambProfileLogP, countCellGeneUMI, countMatStride, pSolo.umiDedup.countInd.main, countCellGeneUMIindex[icell], nGenePerCB[icell], logFactorial);
        };
        time(&rawTime);
    }
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... finished observed logProb" <<endl;
    
    //simulate the probabilities for each cell count
    vector<vector<double>> simLogProb(pSolo.cellFilter.eDcr.simN);
    {
        std::discrete_distribution<uint32> distrAmb ( ambProfilePnon0.begin(), ambProfilePnon0.end() );
        auto maxCount=indCount[iCandFirst].count;
        
        //#pragma omp parallel for num_threads(P.runThreadN) //does not increase speed significantly - might be useful for larger number of simulations
        for (uint64 isim=0; isim<simLogProb.size(); isim++) {
            simLogProb[isim].resize(maxCount+1);
            simLogProb[isim][0]=0;

            std::mt19937 rngGen(19760110LLU*(isim+1));

            vector<uint32> currCounts(ambProfilePnon0.size(), 0);
            for (uint32 ic=1; ic<=maxCount; ic++) {
                uint32 ig1 = distrAmb(rngGen);
                currCounts[ig1]++;
                simLogProb[isim][ic] = simLogProb[isim][ic-1] + ambProfileLogPnon0[ig1] + std::log(ic) - std::log(currCounts[ig1]);
            };
        };
    };
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... finished simulations" <<endl;
    
    
    //p-values
    typedef struct{uint32 index; double p; double padj;} IndPPadj;
    vector<IndPPadj> pValues(obsLogProb.size());
    {
        for (uint32 icand=0; icand<obsLogProb.size(); icand++) {
            pValues[icand].index=indCount[icand+iCandFirst].index;
            auto count1=indCount[icand+iCandFirst].count;
            
            //auto funSumLess = [&] (uint32 n, vector<double> sp) { return n + (sp[count1]<obsLogProb[icand]); };
            //uint32 nLowerP = std::accumulate<uint32>(simLogProb.begin(), simLogProb.end(), 0, funSumLess);
            uint32 nLowerP=0;
            //for (uint64 isim=0; isim<simLogProb.size(); isim++) {
            for (auto &sp: simLogProb) {
                nLowerP += ( sp[count1]<obsLogProb[icand] );
            };

            pValues[icand].p=double(1+nLowerP)/(1+simLogProb.size());
        };
        //BH
        std::sort(pValues.begin(), pValues.end(), [](const IndPPadj &ip1, const IndPPadj &ip2) {return (ip1.p < ip2.p);} );
        uint32 rank=0;
        for (auto &ip: pValues) {
            rank++;
            ip.padj=ip.p*pValues.size()/rank;
        };
        for (auto ip=pValues.rbegin()+1; ip!=pValues.rend(); ++ip)
            ip->padj = min(ip->padj, (ip-1)->padj); //make it non-decreasing
    };
    
    uint32 extraCells=0;
    for (auto &ip: pValues) {
        if (ip.padj<=pSolo.cellFilter.eDcr.FDR) {
            ++extraCells;
            filteredCells.filtVecBool[ip.index]=true;
        };
    };
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) <<  " ... finished emptyDrops_CR filtering: number of additional non-ambient cells=" << extraCells <<endl;

    // Emit branch-only audit files so legacy and libscrna can be compared on the same candidate surface.
    {
        string edOutDir = outputPrefixFiltered + "EmptyDrops";
        mkdir(edOutDir.c_str(), P.runDirPerm);

        uint32 ambientStart = min(nCB, pSolo.cellFilter.eDcr.indMin);
        uint32 ambientEnd = min(nCB, pSolo.cellFilter.eDcr.indMax);

        unordered_map<uint32, pair<double, double>> pByIndex;
        pByIndex.reserve(pValues.size() * 2);
        for (const auto &ip : pValues) {
            pByIndex.emplace(ip.index, make_pair(ip.p, ip.padj));
        }

        unordered_map<uint32, double> obsByIndex;
        obsByIndex.reserve(obsLogProb.size() * 2);
        for (uint32 icand = 0; icand < obsLogProb.size(); icand++) {
            obsByIndex.emplace(indCount[iCandFirst + icand].index, obsLogProb[icand]);
        }

        ofstream summaryOut(edOutDir + "/backend_debug_summary.json");
        if (summaryOut.is_open()) {
            summaryOut << "{\n";
            summaryOut << "  \"backend\": \"legacy\",\n";
            summaryOut << "  \"n_total_cells\": " << nCB << ",\n";
            summaryOut << "  \"n_simple_cells\": " << filteredCells.nCellsSimple << ",\n";
            summaryOut << "  \"n_candidates_total\": " << (iCandLast + 1) << ",\n";
            summaryOut << "  \"n_tail_candidates\": " << obsLogProb.size() << ",\n";
            summaryOut << "  \"candidate_first_rank\": " << iCandFirst << ",\n";
            summaryOut << "  \"candidate_last_rank\": " << iCandLast << ",\n";
            summaryOut << "  \"min_umi\": " << minUMI << ",\n";
            summaryOut << "  \"n_ambient_cells\": " << (ambientEnd > ambientStart ? (ambientEnd - ambientStart) : 0) << ",\n";
            summaryOut << "  \"ambient_start_rank\": " << ambientStart << ",\n";
            summaryOut << "  \"ambient_end_rank_exclusive\": " << ambientEnd << ",\n";
            summaryOut << "  \"n_ed_passers\": " << extraCells << ",\n";
            summaryOut << "  \"fdr_threshold\": " << pSolo.cellFilter.eDcr.FDR << ",\n";
            summaryOut << "  \"raw_p_threshold_for_audit\": 0.05,\n";
            summaryOut << "  \"sim_n\": " << pSolo.cellFilter.eDcr.simN << "\n";
            summaryOut << "}\n";
        }

        ofstream candidateOut(edOutDir + "/backend_debug_candidates.tsv");
        if (candidateOut.is_open()) {
            candidateOut << "backend\trank_desc\tbarcode\tcell_index\tumi_count\tis_simple_cell\tis_tail_candidate\tis_ambient_cell\tp_value\tp_adjusted\tpasses_raw_p\tpasses_fdr\tobs_log_prob\n";
            for (uint32 rank = 0; rank <= iCandLast; rank++) {
                uint32 cellIndex = indCount[rank].index;
                bool isSimple = rank < iCandFirst;
                bool isAmbient = rank >= ambientStart && rank < ambientEnd;
                auto pIt = pByIndex.find(cellIndex);
                auto obsIt = obsByIndex.find(cellIndex);
                string pVal = isSimple ? "0" : (pIt != pByIndex.end() ? to_string(pIt->second.first) : "");
                string pAdj = isSimple ? "0" : (pIt != pByIndex.end() ? to_string(pIt->second.second) : "");
                string passRaw = isSimple ? "1" : (pIt != pByIndex.end() && pIt->second.first <= 0.05 ? "1" : "0");
                string passFdr = isSimple ? "1" : (pIt != pByIndex.end() && pIt->second.second <= pSolo.cellFilter.eDcr.FDR ? "1" : "0");
                string obsLog = obsIt != obsByIndex.end() ? to_string(obsIt->second) : "";

                candidateOut << "legacy\t"
                             << rank << '\t'
                             << pSolo.cbWLstr[indCB[cellIndex]] << '\t'
                             << cellIndex << '\t'
                             << indCount[rank].count << '\t'
                             << (isSimple ? 1 : 0) << '\t'
                             << (isSimple ? 0 : 1) << '\t'
                             << (isAmbient ? 1 : 0) << '\t'
                             << pVal << '\t'
                             << pAdj << '\t'
                             << passRaw << '\t'
                             << passFdr << '\t'
                             << obsLog << '\n';
            }
        }

        ofstream ambientOut(edOutDir + "/backend_debug_ambient.tsv");
        if (ambientOut.is_open()) {
            ambientOut << "backend\trank_desc\tbarcode\tcell_index\tumi_count\n";
            for (uint32 rank = ambientStart; rank < ambientEnd; rank++) {
                uint32 cellIndex = indCount[rank].index;
                ambientOut << "legacy\t"
                           << rank << '\t'
                           << pSolo.cbWLstr[indCB[cellIndex]] << '\t'
                           << cellIndex << '\t'
                           << indCount[rank].count << '\n';
            }
        }
    }
    
    return;
};

double logMultinomialPDFsparse(const vector<double> &ambProfileLogP, const vector<uint32> &countCellGeneUMI, const uint32 stride, const uint32 shift, const int64 start, const uint32 nGenes, const vector<double> &logFactorial)
{
    uint32 sumCount=0;
    double sumLogFac=0.0, sumCountLogP=0.0;
    for (uint32 ig=0; ig<nGenes; ig++) {
        auto count1 = countCellGeneUMI[start+ig*stride+shift];
        sumCount += count1;
        sumLogFac += logFactorial[count1];
        sumCountLogP += ambProfileLogP[countCellGeneUMI[start+ig*stride]] * count1;
    };
    
    return logFactorial[sumCount] - sumLogFac + sumCountLogP;
};
