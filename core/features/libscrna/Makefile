# libscrna - Shared single-cell RNA analysis library
# Contains: EmptyDrops, OrdMag (Simple EmptyDrops), Occupancy filtering

CXX := g++
AR := ar
CXXFLAGS := -std=c++11 -O3 -Wall -Wextra -fPIC
INCLUDES := -I./include

# Source files
SRCS := src/EmptyDropsCRSampler.cpp \
        src/EmptyDropsMultinomial.cpp \
        src/OrdMagStage.cpp \
        src/OccupancyGuard.cpp \
        src/scrna_api.cpp

# C++ source files (.cc) — separate list so we can keep the .cpp/.cc split
# explicit. AtacEvidenceFromPeaks is consumed both by the standalone CLI in
# tools/ and by STAR's chromap orchestration (so it has to live inside
# libscrna.a, not just in the tool).
CC_SRCS := src/AtacEvidenceFromPeaks.cc

# Object files
OBJS := $(SRCS:.cpp=.o) $(CC_SRCS:.cc=.o)

# Output library
TARGET := libscrna.a

# Public headers (for dependent projects)
PUBLIC_HEADERS := include/scrna_types.h \
                  include/SampleMatrixData.h \
                  include/EmptyDropsCRSampler.h \
                  include/EmptyDropsMultinomial.h \
                  include/OrdMagStage.h \
                  include/OccupancyGuard.h \
                  include/CRLogProb.h \
                  include/scrna_api.h

.PHONY: all clean lib test tools

all: lib

lib: $(TARGET)

# Simple EmptyDrops CLI (wrapper around libscrna C API)
TOOLS_DIR := tools
BIN_DIR := bin
SIMPLEED_TOOL := $(BIN_DIR)/scrna_simpleed
ARC_TABLE_TOOL := $(BIN_DIR)/scrna_arc_barcode_table
MULTIOME_COMPARE_TOOL := $(BIN_DIR)/scrna_compare_multiome_calls
MULTIOME_COMBINE_TOOL := $(BIN_DIR)/scrna_multiome_combine
GEX_EVIDENCE_TOOL := $(BIN_DIR)/scrna_build_gex_evidence
ATAC_EVIDENCE_TOOL := $(BIN_DIR)/scrna_build_atac_evidence

tools: lib $(SIMPLEED_TOOL) $(ARC_TABLE_TOOL) $(MULTIOME_COMPARE_TOOL) $(MULTIOME_COMBINE_TOOL) $(GEX_EVIDENCE_TOOL) $(ATAC_EVIDENCE_TOOL)

$(SIMPLEED_TOOL): $(TOOLS_DIR)/scrna_simpleed.cpp $(TARGET)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(TARGET)
	@echo "Built $(SIMPLEED_TOOL)"

$(ARC_TABLE_TOOL): $(TOOLS_DIR)/scrna_arc_barcode_table.cpp
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $<
	@echo "Built $(ARC_TABLE_TOOL)"

$(MULTIOME_COMPARE_TOOL): $(TOOLS_DIR)/scrna_compare_multiome_calls.cpp
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $<
	@echo "Built $(MULTIOME_COMPARE_TOOL)"

$(MULTIOME_COMBINE_TOOL): $(TOOLS_DIR)/scrna_multiome_combine.cpp
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $<
	@echo "Built $(MULTIOME_COMBINE_TOOL)"

$(GEX_EVIDENCE_TOOL): $(TOOLS_DIR)/scrna_build_gex_evidence.cpp
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $<
	@echo "Built $(GEX_EVIDENCE_TOOL)"

$(ATAC_EVIDENCE_TOOL): $(TOOLS_DIR)/scrna_build_atac_evidence.cpp
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $<
	@echo "Built $(ATAC_EVIDENCE_TOOL)"

$(TARGET): $(OBJS)
	$(AR) rcs $@ $^
	@echo "Built $(TARGET)"

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

src/%.o: src/%.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
	rm -f tests/*.o tests/test_*
	rm -f $(SIMPLEED_TOOL) $(ARC_TABLE_TOOL) $(MULTIOME_COMPARE_TOOL) $(MULTIOME_COMBINE_TOOL) $(GEX_EVIDENCE_TOOL) $(ATAC_EVIDENCE_TOOL)

# Dependency on headers
src/EmptyDropsCRSampler.o: include/EmptyDropsCRSampler.h include/scrna_types.h include/pcg_random.hpp
src/EmptyDropsMultinomial.o: include/EmptyDropsMultinomial.h include/EmptyDropsCRSampler.h include/CRLogProb.h include/scrna_types.h include/SimpleGoodTuring/sgt.h
src/OrdMagStage.o: include/OrdMagStage.h include/scrna_types.h include/pcg_random.hpp
src/OccupancyGuard.o: include/OccupancyGuard.h include/scrna_types.h include/pcg_random.hpp

# Test targets (tests are in ../../tests/emptydrops/)
# Run: make -C ../../tests/emptydrops test
test: lib
	@echo "Tests moved to tests/emptydrops/"
	@echo "Run: make -C ../../tests/emptydrops test"
	$(MAKE) -C ../../tests/emptydrops test

# --- ATAC evidence from peaks (self-contained tool) ----------------------
# Consumes a Chromap fragments TSV plus a MACS3 narrowPeak file and emits
# per-barcode ATAC evidence (atac_peak_region_cutsites,
# atac_peak_region_fragments, atac_fragments, atac_peak_fraction).
#
# Deliberately not part of the aggregate `tools:` target because that target
# also depends on sibling tool sources that may not be tracked. Build this
# tool explicitly:
#   make -C core/features/libscrna scrna_build_atac_evidence_from_peaks
#   make -C core/features/libscrna clean-atac-evidence-from-peaks
ATAC_EVIDENCE_FROM_PEAKS_TOOL := $(BIN_DIR)/scrna_build_atac_evidence_from_peaks

.PHONY: scrna_build_atac_evidence_from_peaks clean-atac-evidence-from-peaks

scrna_build_atac_evidence_from_peaks: $(ATAC_EVIDENCE_FROM_PEAKS_TOOL)

$(ATAC_EVIDENCE_FROM_PEAKS_TOOL): $(TOOLS_DIR)/scrna_build_atac_evidence_from_peaks.cpp $(TARGET)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(TARGET)
	@echo "Built $(ATAC_EVIDENCE_FROM_PEAKS_TOOL)"

clean-atac-evidence-from-peaks:
	rm -f $(ATAC_EVIDENCE_FROM_PEAKS_TOOL)
