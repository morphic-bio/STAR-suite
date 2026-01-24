# user may define the whole sets of flags
# LDFLAGS
# CPPFLAGS
# CXXFLAGS
# CFLAGS

# or these user-set flags that will be added to standard flags
LDFLAGSextra ?=
CXXFLAGSextra ?=

# user may define the compiler
CXX ?= g++

# Repo layout (STAR-suite)
ROOT_DIR := $(abspath ../../..)
FLEX_DIR := $(ROOT_DIR)/flex
SLAM_DIR := $(ROOT_DIR)/slam
VBEM_DIR := $(ROOT_DIR)/core/features/vbem
YREMOVE_FASTQ_DIR := $(ROOT_DIR)/core/features/yremove_fastq
BAMSORT_DIR := $(ROOT_DIR)/core/features/bamsort
PROCESS_FEATURES_DIR := $(ROOT_DIR)/core/features/process_features

FLEX_SRC_DIR := $(FLEX_DIR)/source
SLAM_SRC_DIR := $(SLAM_DIR)/source
VBEM_SRC_DIR := $(VBEM_DIR)/source
BAMSORT_SRC_DIR := $(BAMSORT_DIR)/source

LIBFLEX_DIR := $(FLEX_SRC_DIR)/libflex
LIBTRIM_DIR := $(VBEM_SRC_DIR)/libtrim
LIBEM_DIR := $(VBEM_SRC_DIR)/libem

FLEXFILTER_DIR := $(FLEX_DIR)/tools/flexfilter
REMOVE_Y_READS_DIR := $(YREMOVE_FASTQ_DIR)/tools/remove_y_reads
TRIMVALIDATE_DIR := $(VBEM_DIR)/tools/trimvalidate

LIBFLEX_LIB := $(LIBFLEX_DIR)/libflex.a
LIBTRIM_LIB := $(LIBTRIM_DIR)/libtrim.a
LIBEM_LIB := $(LIBEM_DIR)/libem.a
PROCESS_FEATURES_DIR := $(ROOT_DIR)/core/features/process_features
PROCESS_FEATURES_LIB := $(PROCESS_FEATURES_DIR)/libprocess_features.a
LIBSCRNA_DIR := $(ROOT_DIR)/core/features/libscrna
LIBSCRNA_LIB := $(LIBSCRNA_DIR)/libscrna.a

VPATH := $(FLEX_SRC_DIR) $(SLAM_SRC_DIR) $(VBEM_SRC_DIR) $(BAMSORT_SRC_DIR)

# pre-defined flags
LDFLAGS_shared := -pthread -Lhtslib -Bstatic -lhts -Bdynamic -lz -lssl -lcrypto
LDFLAGS_static := -static -static-libgcc -pthread -Lhtslib -lhts -lz
LDFLAGS_Mac :=-pthread -lz htslib/libhts.a
LDFLAGS_Mac_static :=-pthread -lz -static-libgcc htslib/libhts.a
LDFLAGS_gdb := $(LDFLAGS_shared)

DATE_FMT = --iso-8601=seconds
ifdef SOURCE_DATE_EPOCH
    BUILD_DATE ?= $(shell date -u -d "@$(SOURCE_DATE_EPOCH)" "$(DATE_FMT)" 2>/dev/null || date -u -r "$(SOURCE_DATE_EPOCH)" "$(DATE_FMT)" 2>/dev/null || date -u "$(DATE_FMT)")
else
    BUILD_DATE ?= $(shell date "$(DATE_FMT)")
endif

BUILD_PLACE ?= $(HOSTNAME):$(shell pwd)

COMPTIMEPLACE := -D'COMPILATION_TIME_PLACE="$(BUILD_DATE) $(BUILD_PLACE)"'


GIT_CHECK := $(shell git status 1> /dev/null 2> /dev/null && echo 0)
ifeq ($(GIT_CHECK),0)
	GIT_DIFF := $(shell git diff --name-only | tr '\n' ' ')
    GIT_BRANCH_COMMIT_DIFF := $(shell git status | head -n 1) ; \
                                                    $(shell git log | head -n 1) ; \
                                                    diff files: $(GIT_DIFF)
endif
GIT_BRANCH_COMMIT_DIFF := -D'GIT_BRANCH_COMMIT_DIFF="$(GIT_BRANCH_COMMIT_DIFF)"'

# Defaults, can be overridden by make arguments or environment
CXXFLAGS ?= -pipe -Wall -Wextra
CFLAGS ?= -pipe -Wall -Wextra -O3
CXXFLAGS_SIMD ?= -mavx2

# Unconditionally set essential flags and optimization options
# Packed ReadInfo is the only supported path (legacy removed)
# Include flex (klib/libflex), vbem (libtrim/libem), slam, and process_features headers
# Include current directory for IncludeDefine.h and other source headers
PROCESS_FEATURES_INC := $(ROOT_DIR)/core/features/process_features/include
CXXFLAGS_common := -std=c++11 -fopenmp -I. -I$(FLEX_SRC_DIR) -I$(LIBFLEX_DIR) -I$(VBEM_SRC_DIR) -I$(LIBEM_DIR) -I$(SLAM_SRC_DIR) -I$(BAMSORT_SRC_DIR) -I$(PROCESS_FEATURES_INC) $(COMPTIMEPLACE) $(GIT_BRANCH_COMMIT_DIFF)
CXXFLAGS_main := -O3 $(CXXFLAGS_common)
CXXFLAGS_gdb := -O0 -g3 $(CXXFLAGS_common)

# AddressSanitizer flags (enabled with ASAN=1)
ifdef ASAN
CFLAGS   += -fsanitize=address -fno-omit-frame-pointer -O1
CXXFLAGS += -fsanitize=address -fno-omit-frame-pointer -O1
LDFLAGS  += -fsanitize=address
endif

# CB/UB Parity validation (enabled with DEBUG_CB_UB_PARITY=1)
ifdef DEBUG_CB_UB_PARITY
CXXFLAGSextra += -DDEBUG_CB_UB_PARITY
endif

##########################################################################################################
OBJECTS = systemFunctions.o funPrimaryAlignMark.o \
	SoloFeature_loadRawMatrix.o SoloFeature_emptyDrops_CR.o soloInputFeatureUMI.o SoloFeature_countSmartSeq.o SoloFeature_redistributeReadsByCB.o \
	SoloFeature_quantTranscript.o SoloFeature_sumThreads.o SoloFeature_countVelocyto.o SoloFeature_countCBgeneUMI.o \
	SoloFeature_collapseUMI_Graph.o SoloFeature_collapseUMIall.o SoloFeature_collapseUMI_fromHash.o SoloFeature_umiCorrection.o SoloFeature_materializeFromHash.o SoloFeature_metricsOutput.o SoloFeature_writeMexFromInlineHashDedup.o SoloFeature_flexfilter.o SoloFeature_addBAMtags.o SoloFeature_writeTagTable.o \
	MexWriter.o MexWriterUtil.o \
	UMICorrector.o \
	ParametersClip_initialize.o ClipMate_clip.o ClipCR4.o opal/opal.o ClipMate_clipChunk.o ClipMate_initialize.o \
	Transcriptome_classifyAlign.o Transcriptome_geneFullAlignOverlap_ExonOverIntron.o Transcriptome_alignExonOverlap.cpp \
	SoloFeature_cellFiltering.o \
	SoloFeature_statsOutput.o bamSortByCoordinate.o SamtoolsSorter.o SoloBarcode.o \
	ParametersSolo.o SoloRead.o SoloRead_record.o \
	SoloReadBarcode.o SoloReadBarcode_getCBandUMI.o SoloBarcode_extractBarcode.o \
	SoloReadFeature.o SoloReadFeature_record.o SoloReadFeature_record_base.o flex/SoloReadFeature_record_flex.o SoloReadFeature_inputRecords.o \
	Solo.o SoloFeature.o SoloFeature_outputResults.o SoloFeature_processRecords.o SoloFeature_prepareReadInfoOnly.o SoloReadInfoLoader.o SoloReadInfoSink.o \
	PackedReadInfo.o InlineCBCorrection.o solo/CbCorrector.o solo/CbBayesianResolver.o GeneResolver.o \
	ReadAlign_transformGenome.o Genome_transformGenome.o Transcript_convertGenomeCigar.o \
	twoPassRunPass1.o samHeaders.o Genome_genomeLoad.o Genome_genomeOutLoad.o Transcript_transformGenome.o ReadAlign_outputSpliceGraphSAM.o \
	ReadAlign_mapOneReadSpliceGraph.o SpliceGraph.o SpliceGraph_swScoreSpliced.o SpliceGraph_swTraceBack.o \
	SpliceGraph_findSuperTr.o sjAlignSplit.o \
	GTF.o GTF_transcriptGeneSJ.o GTF_superTranscript.o SuperTranscriptome.o TranscriptomeFasta.o \
	ReadAlign_outputAlignments.o ReadAlign_resolveAmbiguousCBs.o ReadAlign_slamQuant.o \
		ReadAlign.o STAR.o TranscriptQuantEC.o LibFormatDetection.o TranscriptQuantOutput.o SlamQuant.o SlamSolver.o SlamCompat.o SlamVarianceAnalysis.o SlamReadBuffer.o SlamDump.o SlamQcOutput.o TrimQc.o TrimQcOutput.o SnpMaskBuild.o \
	SharedMemory.o PackedArray.o SuffixArrayFuns.o Parameters.o Parameters_samAttributes.o InOutStreams.o SequenceFuns.o Genome.o ParametersGenome.o Stats.o \
	Transcript.o Transcript_alignScore.o Transcript_generateCigarP.o Chain.o \
	Transcript_variationAdjust.o Variation.o ReadAlign_waspMap.o \
	ReadAlign_storeAligns.o ReadAlign_stitchPieces.o ReadAlign_multMapSelect.o ReadAlign_mapOneRead.o readLoad.o \
	ReadAlignChunk.o ReadAlignChunk_processChunks.o ReadAlignChunk_mapChunk.o \
	OutSJ.o outputSJ.o blocksOverlap.o ThreadControl.o sysRemoveDir.o \
	ReadAlign_maxMappableLength2strands.o binarySearch2.o\
	ReadAlign_outputTranscriptSAM.o ReadAlign_outputTranscriptSJ.o ReadAlign_outputTranscriptCIGARp.o ReadAlign_calcCIGAR.cpp \
	ReadAlign_createExtendWindowsWithAlign.o ReadAlign_assignAlignToWindow.o ReadAlign_oneRead.o \
	ReadAlign_stitchWindowSeeds.o \
	ReadAlign_peOverlapMergeMap.o ReadAlign_mappedFilter.o \
	ParametersChimeric_initialize.o ReadAlign_chimericDetection.o ReadAlign_chimericDetectionOld.o ReadAlign_chimericDetectionOldOutput.o\
	ChimericDetection.o ChimericDetection_chimericDetectionMult.o ReadAlign_chimericDetectionPEmerged.o \
	stitchWindowAligns.o extendAlign.o stitchAlignToTranscript.o \
	ChimericSegment.cpp ChimericAlign.cpp ChimericAlign_chimericJunctionOutput.o ChimericAlign_chimericBAMoutput.o ChimericAlign_chimericStitching.o \
	Genome_genomeGenerate.o genomeParametersWrite.o genomeScanFastaFiles.o genomeSAindex.o \
	Genome_insertSequences.o insertSeqSA.o funCompareUintAndSuffixes.o funCompareUintAndSuffixesMemcmp.o \
	TimeFunctions.o ErrorWarning.o streamFuns.o stringSubstituteAll.o \
	Transcriptome.o Transcriptome_quantAlign.o Transcriptome_geneFullAlignOverlap.o \
	ReadAlign_quantTranscriptome.o Quantifications.o Transcriptome_geneCountsAddAlign.o \
	sjdbLoadFromFiles.o sjdbLoadFromStream.o sjdbPrepare.o sjdbBuildIndex.o sjdbInsertJunctions.o mapThreadsSpawn.o \
	Parameters_readFilesInit.o Parameters_openReadsFiles.cpp Parameters_closeReadsFiles.cpp Parameters_readSAMheader.o \
	bam_cat.o serviceFuns.o GlobalVariables.cpp \
	BAMoutput.o BAMfunctions.o ReadAlign_alignBAM.o BAMbinSortByCoordinate.o signalFromBAM.o bamRemoveDuplicates.o BAMbinSortUnmapped.o \
	ProbeListIndex.o \
	SampleDetector.o \
	SoloBamParsing.o \
	ZGZXTags.o \
	FlexProbeIndex.o \
	CellRangerFormatter.o \
	CrMultiConfig.o CrMultiMexStub.o CrMultiAssign.o CrMultiMerge.o CrMultiProcess.o

# Unit tests for CR-compat regression components (built as a single binary)
OBJECTS_NO_MAIN := $(filter-out STAR.o,$(OBJECTS))
COMPAT_UNIT_TEST_SRC := $(ROOT_DIR)/tests/compat/test_cr_compat_components.cpp

.PHONY: compat-unit-tests
compat-unit-tests : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
compat-unit-tests : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS) $(shell pkg-config --libs glib-2.0 2>/dev/null || echo "")
compat-unit-tests : Depend.list parametersDefault.xxd $(OBJECTS_NO_MAIN) libflex libtrim libem $(COMPAT_UNIT_TEST_SRC)
	$(CXX) -o compat_unit_tests $(CXXFLAGS) $(COMPAT_UNIT_TEST_SRC) $(OBJECTS_NO_MAIN) $(LIBFLEX_LIB) $(LIBTRIM_LIB) $(LIBEM_LIB) $(LDFLAGS)

SOURCES := $(wildcard *.cpp) $(wildcard *.c) $(wildcard flex/*.cpp)

%.o : %.c
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) $<

all: cleanCompileInfo STAR$(SFX)

# Standalone FlexFilter MEX tool (optional, not built by default)
# Requires STAR to be built first (needs libflex.a and MexWriter.o)
.PHONY: flexfilter
flexfilter: STAR$(SFX)
	$(MAKE) -C $(FLEXFILTER_DIR)

# Standalone remove_y_reads FASTQ splitter tool (optional, not built by default)
# Requires htslib to be built first
.PHONY: remove_y_reads
remove_y_reads: htslib/libhts.a
	$(MAKE) -C $(REMOVE_Y_READS_DIR)

# Standalone star_feature_call tool for perturb-seq processing
# Requires libprocess_features.a to be built first
.PHONY: star_feature_call libprocess_features
libprocess_features:
	$(MAKE) -C $(PROCESS_FEATURES_DIR) lib

star_feature_call: CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
star_feature_call: LDFLAGS_PF := -L$(PROCESS_FEATURES_DIR) -lprocess_features -L$(LIBSCRNA_DIR) -lscrna -lm -lpthread -lz -fopenmp -lcairo -Lhtslib -lhts
star_feature_call: libprocess_features htslib/libhts.a $(LIBSCRNA_LIB)
	$(CXX) $(CXXFLAGS) -o star_feature_call star_feature_call.cpp $(LDFLAGS_PF)

# Convenience target to rebuild only the MEX writer components
.PHONY: MexWriterObjs
MexWriterObjs: MexWriter.o MexWriterUtil.o
.PHONY: test_packed_readinfo
test_packed_readinfo:
	$(CXX) -std=c++11 -I./ $(CPPFLAGS) $(CXXFLAGS) ../test/PackedReadInfo_test.cpp PackedReadInfo.cpp -o ../test/pri_test
	../test/pri_test

.PHONY: test_flex_probe_index
test_flex_probe_index: FlexProbeIndex.o
	$(CXX) -std=c++11 -I./ $(CPPFLAGS) $(CXXFLAGS) ../test/test_flex_probe_index.cpp FlexProbeIndex.o -lssl -lcrypto -o ../test/test_flex_probe_index

.PHONY: test_cellranger_format
test_cellranger_format: CellRangerFormatter.o
	$(CXX) -std=c++11 -I./ $(CPPFLAGS) $(CXXFLAGS) ../test/test_cellranger_format.cpp CellRangerFormatter.o -Lhtslib -lhts -lz -o ../test/test_cellranger_format

.PHONY: test_transcriptome_fasta
test_transcriptome_fasta: TranscriptomeFasta.o Genome.o GTF.o GTF_transcriptGeneSJ.o Parameters.o ParametersGenome.o InOutStreams.o SequenceFuns.o ErrorWarning.o streamFuns.o serviceFuns.o
	$(CXX) -std=c++11 -I./ $(CPPFLAGS) $(CXXFLAGS) ../test/test_transcriptome_fasta.cpp TranscriptomeFasta.o Genome.o GTF.o GTF_transcriptGeneSJ.o Parameters.o ParametersGenome.o InOutStreams.o SequenceFuns.o ErrorWarning.o streamFuns.o serviceFuns.o -Lhtslib -lhts -lz -o ../test/test_transcriptome_fasta

.PHONY: test_cellranger_parity
test_cellranger_parity: test_cellranger_format
	../test/run_cellranger_format_parity.sh

.PHONY: test_cellranger_golden
test_cellranger_golden: test_cellranger_format
	../test/test_cellranger_format \
		../test/fixtures/cellranger_format/test_genome.fa \
		../test/fixtures/cellranger_format/test_genes.gtf \
		/tmp/cpp_genome.fa /tmp/cpp_genes.gtf
	diff ../test/fixtures/cellranger_format/expected_genome.fa /tmp/cpp_genome.fa
	diff ../test/fixtures/cellranger_format/expected_genes.gtf /tmp/cpp_genes.gtf

.PHONY: test_cksum
test_cksum: CellRangerFormatter.o
	$(CXX) -std=c++11 -I./ $(CPPFLAGS) $(CXXFLAGS) ../test/test_cksum.cpp CellRangerFormatter.o -Lhtslib -lhts -lz -o ../test/test_cksum
	../test/test_cksum

.PHONY: test_download_untrusted
test_download_untrusted: CellRangerFormatter.o
	$(CXX) -std=c++11 -I./ $(CPPFLAGS) $(CXXFLAGS) ../test/test_download_untrusted.cpp CellRangerFormatter.o -Lhtslib -lhts -lz -o ../test/test_download_untrusted
	../test/test_download_untrusted

.PHONY: test_replace_unverifiable
test_replace_unverifiable: CellRangerFormatter.o
	$(CXX) -std=c++11 -I./ $(CPPFLAGS) $(CXXFLAGS) ../test/test_replace_unverifiable.cpp CellRangerFormatter.o -Lhtslib -lhts -lz -o ../test/test_replace_unverifiable
	../test/test_replace_unverifiable

.PHONY: test_autocksum
test_autocksum: CellRangerFormatter.o
	$(CXX) -std=c++11 -I./ $(CPPFLAGS) $(CXXFLAGS) ../test/test_autocksum.cpp CellRangerFormatter.o -Lhtslib -lhts -lz -o ../test/test_autocksum
	../test/test_autocksum

.PHONY: test_SoloReadInfoBaseline
# NOTE: This target links against most STAR objects (excluding STAR.o) which is heavy for a tiny harness.
# Once the loader abstraction is introduced, consider extracting a smaller test lib to reduce build time.
test_SoloReadInfoBaseline: CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
test_SoloReadInfoBaseline: LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS)
test_SoloReadInfoBaseline: Depend.list parametersDefault.xxd $(filter-out STAR.o,$(OBJECTS))
	$(CXX) $(CXXFLAGS) ../test/SoloReadInfoBaseline_test.cpp $(filter-out STAR.o,$(OBJECTS)) $(LDFLAGS) -o ../test/solo_readinfo_baseline
	../test/solo_readinfo_baseline

.PHONY: test
test: test_packed_readinfo test_SoloReadInfoBaseline test_SoloReadInfoLoader test_SoloReadInfoParity test_MinimalSink test_CbCorrector test_CbBayesianResolver test_CbBayesianResolverPhase2 test_inline_hash_record_capture test_inline_hash_resolver test_inline_hash_clique test_inline_hash_parity test_gene_resolver

.PHONY: test_MinimalSink
test_MinimalSink: CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
test_MinimalSink: LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS)
test_MinimalSink: Depend.list parametersDefault.xxd $(filter-out STAR.o,$(OBJECTS))
ifeq ($(STAR_TEST_SKIP_TXOME),)
	$(CXX) $(CXXFLAGS) ../test/MinimalSink_test.cpp $(filter-out STAR.o,$(OBJECTS)) $(LDFLAGS) -o ../test/minimal_sink_test
	../test/minimal_sink_test
else
	@echo "Skipping MinimalSink_test due to STAR_TEST_SKIP_TXOME=1 or missing STAR_TXOME_DIR"
endif
# Loader unit test
.PHONY: test_SoloReadInfoLoader
test_SoloReadInfoLoader: CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
test_SoloReadInfoLoader: LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS)
test_SoloReadInfoLoader: Depend.list parametersDefault.xxd $(filter-out STAR.o,$(OBJECTS))
	$(CXX) $(CXXFLAGS) ../test/SoloReadInfoLoader_test.cpp $(filter-out STAR.o,$(OBJECTS)) $(LDFLAGS) -o ../test/solo_readinfo_loader
	../test/solo_readinfo_loader

# Parity test (loader minimal vs counting)
.PHONY: test_SoloReadInfoParity
test_SoloReadInfoParity: CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
test_SoloReadInfoParity: LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS)
test_SoloReadInfoParity: Depend.list parametersDefault.xxd $(filter-out STAR.o,$(OBJECTS))
	$(CXX) $(CXXFLAGS) ../test/SoloReadInfoParity_test.cpp $(filter-out STAR.o,$(OBJECTS)) $(LDFLAGS) -o ../test/solo_readinfo_parity
	../test/solo_readinfo_parity

# CbCorrector unit test
.PHONY: test_CbCorrector
test_CbCorrector: CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
test_CbCorrector: solo/CbCorrector.o
	$(CXX) $(CXXFLAGS) -I./ ../test/cb_corrector_test.cpp solo/CbCorrector.o $(LDFLAGS_shared) -o ../test/cb_corrector_test

.PHONY: test_CbBayesianResolver
test_CbBayesianResolver: CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
test_CbBayesianResolver: solo/CbBayesianResolver.o
	$(CXX) $(CXXFLAGS) -I./ ../test/cb_bayesian_resolver_test.cpp solo/CbBayesianResolver.o $(LDFLAGS_shared) -o ../test/cb_bayesian_resolver_test

.PHONY: test_CbBayesianResolverPhase2
test_CbBayesianResolverPhase2: CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
test_CbBayesianResolverPhase2: solo/CbBayesianResolver.o
	$(CXX) $(CXXFLAGS) -I./ ../test/cb_bayesian_resolver_phase2_test.cpp solo/CbBayesianResolver.o $(LDFLAGS_shared) -o ../test/cb_bayesian_resolver_phase2_test
	../test/cb_bayesian_resolver_phase2_test

.PHONY: test_inline_hash_record_capture
test_inline_hash_record_capture: CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS) -DINLINE_HASH_STANDALONE=0
test_inline_hash_record_capture: LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS) $(shell pkg-config --libs glib-2.0 2>/dev/null || echo "")
test_inline_hash_record_capture: Depend.list parametersDefault.xxd $(filter-out STAR.o,$(OBJECTS))
	$(CXX) $(CXXFLAGS) ../test/test_inline_hash_record_capture.cpp $(filter-out STAR.o,$(OBJECTS)) $(LDFLAGS) -o ../test/test_inline_hash_record_capture
	../test/test_inline_hash_record_capture

.PHONY: test_inline_hash_resolver
test_inline_hash_resolver: CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS) -DINLINE_HASH_STANDALONE=0
test_inline_hash_resolver: LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS) $(shell pkg-config --libs glib-2.0 2>/dev/null || echo "")
test_inline_hash_resolver: Depend.list parametersDefault.xxd $(filter-out STAR.o,$(OBJECTS))
	$(CXX) $(CXXFLAGS) ../test/test_inline_hash_resolver.cpp $(filter-out STAR.o,$(OBJECTS)) $(LDFLAGS) -o ../test/test_inline_hash_resolver
	../test/test_inline_hash_resolver

.PHONY: test_inline_hash_clique
test_inline_hash_clique: CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS) -DINLINE_HASH_STANDALONE=0
test_inline_hash_clique: LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS) $(shell pkg-config --libs glib-2.0 2>/dev/null || echo "")
test_inline_hash_clique: Depend.list parametersDefault.xxd $(filter-out STAR.o,$(OBJECTS))
	$(CXX) $(CXXFLAGS) ../test/test_inline_hash_clique.cpp $(filter-out STAR.o,$(OBJECTS)) $(LDFLAGS) -o ../test/test_inline_hash_clique
	../test/test_inline_hash_clique

.PHONY: test_inline_hash_parity
test_inline_hash_parity: CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS) -DINLINE_HASH_STANDALONE=0
test_inline_hash_parity: LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS) $(shell pkg-config --libs glib-2.0 2>/dev/null || echo "")
test_inline_hash_parity: Depend.list parametersDefault.xxd $(filter-out STAR.o,$(OBJECTS))
	$(CXX) $(CXXFLAGS) ../test/test_inline_hash_parity.cpp $(filter-out STAR.o,$(OBJECTS)) $(LDFLAGS) -o ../test/test_inline_hash_parity
	../test/test_inline_hash_parity

.PHONY: test_gene_resolver
test_gene_resolver: CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
test_gene_resolver: GeneResolver.o
	$(CXX) $(CXXFLAGS) -I./ ../test/test_gene_resolver.cpp GeneResolver.o $(LDFLAGS_shared) -o ../test/test_gene_resolver
	../test/test_gene_resolver

.PHONY: test_SamtoolsSorter
test_SamtoolsSorter: CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
test_SamtoolsSorter: LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS)
test_SamtoolsSorter: Depend.list parametersDefault.xxd libflex libtrim libem $(filter-out STAR.o,$(OBJECTS))
	$(CXX) $(CXXFLAGS) ../test/SamtoolsSorter_test.cpp $(filter-out STAR.o,$(OBJECTS)) $(LIBFLEX_LIB) $(LIBTRIM_LIB) $(LIBEM_LIB) $(LDFLAGS) -o ../test/samtools_sorter_test
	../test/samtools_sorter_test

.PHONY: test_trim_parity
test_trim_parity: libtrim
	$(MAKE) -C $(TRIMVALIDATE_DIR) test

.PHONY: test_trim_e2e
test_trim_e2e: STAR
	@echo "Running end-to-end trimming parity test..."
	@bash ../test/integration/trim/run_e2e_trim_parity.sh

.PHONY: ctest
ctest: test

# Pattern rule for root-level .cpp files to ensure CXXFLAGS_common (including third_party) is used
# This applies to files like SoloFeature_flexfilter.cpp that need the third_party include path
# Subdirectories (opal/, solo/, flex/) have their own specific rules that take precedence
%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(CXXFLAGS_common) $< -o $@

flex/%.o: flex/%.cpp
	@mkdir -p flex
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(CXXFLAGS_common) $< -o $@

opal/opal.o : opal/opal.cpp opal/opal.h
	cd opal && \
	$(CXX) -c -I./ -std=c++11 $(CPPFLAGS) $(CXXFLAGS) $(CXXFLAGSextra) $(CXXFLAGS_SIMD) opal.cpp

solo/CbCorrector.o : $(FLEX_SRC_DIR)/solo/CbCorrector.cpp $(FLEX_SRC_DIR)/solo/CbCorrector.h
	@mkdir -p solo
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(CXXFLAGS_common) -I./ $(FLEX_SRC_DIR)/solo/CbCorrector.cpp -o solo/CbCorrector.o

solo/CbBayesianResolver.o : $(FLEX_SRC_DIR)/solo/CbBayesianResolver.cpp $(FLEX_SRC_DIR)/solo/CbBayesianResolver.h
	@mkdir -p solo
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(CXXFLAGS_common) -I./ $(FLEX_SRC_DIR)/solo/CbBayesianResolver.cpp -o solo/CbBayesianResolver.o

.PHONY: clean
clean:
	'rm' -f *.o opal/opal.o solo/*.o STAR STARstatic STARlong star_feature_call Depend.list parametersDefault.xxd
	$(MAKE) -C $(LIBFLEX_DIR) clean
	$(MAKE) -C $(LIBTRIM_DIR) clean
	$(MAKE) -C $(LIBEM_DIR) clean

.PHONY: CLEAN
CLEAN: clean
	$(MAKE) -C htslib clean


.PHONY: clean_solo
clean_solo:
	'rm' -f Solo*.o

.PHONY: cleanCompileInfo
cleanCompileInfo:
	'rm' -f STAR.o Parameters.o

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),CLEAN)
ifneq ($(MAKECMDGOALS),CLEAN)
ifneq ($(MAKECMDGOALS),clean_solo)
ifneq ($(MAKECMDGOALS),STARforMac)
ifneq ($(MAKECMDGOALS),STARforMacGDB)
Depend.list: $(SOURCES) parametersDefault.xxd htslib
	echo $(SOURCES)
	'rm' -f ./Depend.list
	$(CXX) $(CXXFLAGS_common) -MM $(SOURCES) parametersDefault.xxd htslib >> Depend.list
include Depend.list
endif
endif
endif
endif
endif
endif

htslib : htslib/libhts.a

htslib/libhts.a :
	$(MAKE) -C htslib lib-static

parametersDefault.xxd: parametersDefault
	xxd -i parametersDefault > parametersDefault.xxd

# Ensure libflex is always rebuilt (sources may have changed)
.PHONY: libflex
libflex:
	$(MAKE) -C $(LIBFLEX_DIR) libflex.a

# Ensure libtrim is always rebuilt (sources may have changed)
.PHONY: libtrim
libtrim:
	$(MAKE) -C $(LIBTRIM_DIR) libtrim.a

# Ensure libem is always rebuilt (sources may have changed)
.PHONY: libem
libem:
	$(MAKE) -C $(LIBEM_DIR) libem.a

STAR$(SFX) : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
STAR$(SFX) : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS) $(shell pkg-config --libs glib-2.0 2>/dev/null || echo "") -lcairo
STAR$(SFX) : Depend.list parametersDefault.xxd $(OBJECTS) libflex libtrim libem libprocess_features $(LIBSCRNA_LIB)
	$(CXX) -o STAR$(SFX) $(CXXFLAGS) $(OBJECTS) $(LIBFLEX_LIB) $(LIBTRIM_LIB) $(LIBEM_LIB) $(PROCESS_FEATURES_LIB) $(LIBSCRNA_LIB) $(LDFLAGS)

STARstatic$(SFX) : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
STARstatic$(SFX) : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_static) $(LDFLAGS)
STARstatic$(SFX) : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR$(SFX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARlong$(SFX) : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) -D'COMPILE_FOR_LONG_READS' $(CXXFLAGS)
STARlong$(SFX) : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS)
STARlong$(SFX) : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STARlong$(SFX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARlongStatic$(SFX) : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) -D'COMPILE_FOR_LONG_READS' $(CXXFLAGS)
STARlongStatic$(SFX) : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_static) $(LDFLAGS)
STARlongStatic$(SFX) : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STARlong$(SFX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)



POSIXSHARED : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) -DPOSIX_SHARED_MEM $(CXXFLAGS)
POSIXSHARED : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS)
POSIXSHARED : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR$(SFX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

gdb : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_gdb) $(CXXFLAGS)
gdb : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_gdb) $(LDFLAGS) $(shell pkg-config --libs glib-2.0 2>/dev/null || echo "") -lcairo
gdb : Depend.list parametersDefault.xxd $(OBJECTS) libflex libtrim libem libprocess_features $(LIBSCRNA_LIB)
	$(CXX) -o STAR $(CXXFLAGS) $(OBJECTS) $(LIBFLEX_LIB) $(LIBTRIM_LIB) $(LIBEM_LIB) $(PROCESS_FEATURES_LIB) $(LIBSCRNA_LIB) $(LDFLAGS)

gdb-long : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_gdb) -D'COMPILE_FOR_LONG_READS' $(CXXFLAGS)
gdb-long : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_gdb) $(LDFLAGS) $(shell pkg-config --libs glib-2.0 2>/dev/null || echo "")
gdb-long : Depend.list parametersDefault.xxd $(OBJECTS) libflex libtrim libem
	$(CXX) -o STARlong $(CXXFLAGS) $(OBJECTS) $(LIBFLEX_LIB) $(LIBTRIM_LIB) $(LIBEM_LIB) $(LDFLAGS)

STARforMacStatic : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) -D'COMPILE_FOR_MAC' $(CXXFLAGS)
STARforMacStatic : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_Mac_static) $(LDFLAGS)
STARforMacStatic : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STAR $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

STARlongForMacStatic : CXXFLAGS := -D'COMPILE_FOR_LONG_READS' $(CXXFLAGSextra) $(CXXFLAGS_main) -D'COMPILE_FOR_MAC' $(CXXFLAGS)
STARlongForMacStatic : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_Mac_static) $(LDFLAGS)
STARlongForMacStatic : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o STARlong $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)
