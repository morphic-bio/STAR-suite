CC=gcc
CXX=g++
AR=ar

# Directories
SRCDIR = src
INCDIR = include
VPATH = $(SRCDIR)

# libscrna (EmptyDrops library)
LIBSCRNA_DIR = ../libscrna
LIBSCRNA_INC = $(LIBSCRNA_DIR)/include
LIBSCRNA = $(LIBSCRNA_DIR)/libscrna.a

# Library flags (no cairo - heatmaps now use Plotly HTML+JSON)
HTS_LIBS := -lhts

DEFINES :=

# Source files (basenames only)
# Core library sources (shared by all tools)
# Note: plot_histogram.c is needed by assignBarcodes.c for histogram generation
# Note: heatmap.c now uses Plotly (HTML+JSON), no external dependencies
# pf_counts.c and mex_writer.c provide modular counts/output stages
CORE_SRCS = assignBarcodes.c queue.c globals.c utils.c memory.c io.c barcode_match.c pf_api.c call_features.c gmm.c mex_10x.c plot_histogram.c pf_counts.c mex_writer.c emptydrops_bridge.c barcode_filter.c heatmap.c

# CLI-specific sources
MAIN_SRCS = main.c
DEMUX_FASTQ_SRCS = demux_fastq.c
DEMUX_BAM_SRCS = demux_bam.c

# Full source lists for each target
SRCS_BASE = $(MAIN_SRCS) $(CORE_SRCS)
SRCS = $(SRCS_BASE)

# Compiler flags: optimise by default, but if DEBUG=1 build with symbols and no optimisations
ifeq ($(DEBUG),1)
    CFLAGS=-g -ggdb -Wall -O0 -I$(INCDIR) -I$(LIBSCRNA_INC) -fopenmp $(DEFINES)
else
    CFLAGS=-g -Wall -O3 -I$(INCDIR) -I$(LIBSCRNA_INC) -fopenmp $(DEFINES)
endif
LDFLAGS=-lm -lpthread -lz -fopenmp -lstdc++ $(HTS_LIBS)

# Object files
OBJS=$(SRCS:.c=.o)
CORE_OBJS=$(CORE_SRCS:.c=.o)

# Executable names
TARGET=assignBarcodes
DEMUX_TARGET=demux_fastq
DEMUX_BAM_TARGET=demux_bam
CALL_FEATURES_TARGET=call_features

# Library target
LIB_TARGET=libprocess_features.a

# demux_fastq sources/objects
DEMUX_SRCS=demux_fastq.c io.c utils.c globals.c barcode_match.c
DEMUX_OBJS=$(DEMUX_SRCS:.c=.o)

DEMUX_BAM_SRCS_FULL=demux_bam.c utils.c globals.c barcode_match.c io.c memory.c
DEMUX_BAM_OBJS=$(DEMUX_BAM_SRCS_FULL:.c=.o)

# Test programs
TEST_API_TARGET=tests/test_pf_api
TEST_PF_COUNTS_TARGET=tests/test_pf_counts
TEST_MEX_WRITER_TARGET=tests/test_mex_writer
TEST_MEX_PARITY_TARGET=tests/test_mex_writer_parity
TEST_BARCODE_FILTER_TARGET=tests/test_barcode_filter
TEST_OFFSET_DETECTION_TARGET=tests/test_offset_detection

.PHONY: all clean lib tools test

all: lib tools

test: $(TEST_API_TARGET) $(TEST_PF_COUNTS_TARGET) $(TEST_MEX_WRITER_TARGET) $(TEST_MEX_PARITY_TARGET) $(TEST_BARCODE_FILTER_TARGET) $(TEST_OFFSET_DETECTION_TARGET) test-regression

lib: $(LIB_TARGET)

tools: $(TARGET) $(DEMUX_TARGET) $(DEMUX_BAM_TARGET) $(CALL_FEATURES_TARGET)

# Build libscrna first
$(LIBSCRNA):
	$(MAKE) -C $(LIBSCRNA_DIR) lib

$(LIB_TARGET): $(CORE_OBJS)
	$(AR) rcs $@ $^

$(TARGET): main.o $(LIB_TARGET) $(LIBSCRNA)
	$(CXX) -o $@ main.o $(LIB_TARGET) $(LIBSCRNA) $(LDFLAGS)

$(DEMUX_TARGET): $(DEMUX_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(DEMUX_BAM_TARGET): $(DEMUX_BAM_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(CALL_FEATURES_TARGET): call_features_main.o $(LIB_TARGET) $(LIBSCRNA)
	$(CXX) -o $@ call_features_main.o $(LIB_TARGET) $(LIBSCRNA) $(LDFLAGS)

# plot_histogram.o is needed by assignBarcodes but not in core lib
plot_histogram.o: plot_histogram.c
	$(CC) $(CFLAGS) -c $< -o $@

# The implicit rule will now use VPATH to find the .c files in src/
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_API_TARGET): tests/test_pf_api.c $(LIB_TARGET) $(LIBSCRNA)
	$(CXX) $(CFLAGS) -o $@ $< $(LIB_TARGET) $(LIBSCRNA) $(LDFLAGS)

# Note: test_pf_counts.c is pure C, compile with gcc then link with g++ for libscrna
$(TEST_PF_COUNTS_TARGET): tests/test_pf_counts.c $(LIB_TARGET) $(LIBSCRNA)
	$(CC) $(CFLAGS) -c tests/test_pf_counts.c -o tests/test_pf_counts.o
	$(CXX) -o $@ tests/test_pf_counts.o $(LIB_TARGET) $(LIBSCRNA) $(LDFLAGS)

# Note: test_mex_writer.c is pure C, compile with gcc then link with g++ for libscrna
$(TEST_MEX_WRITER_TARGET): tests/test_mex_writer.c $(LIB_TARGET) $(LIBSCRNA)
	$(CC) $(CFLAGS) -c tests/test_mex_writer.c -o tests/test_mex_writer.o
	$(CXX) -o $@ tests/test_mex_writer.o $(LIB_TARGET) $(LIBSCRNA) $(LDFLAGS)

# Parity test: compares mex_write_all() output vs printFeatureCounts() output
$(TEST_MEX_PARITY_TARGET): tests/test_mex_writer_parity.c $(LIB_TARGET) $(LIBSCRNA)
	$(CC) $(CFLAGS) -c tests/test_mex_writer_parity.c -o tests/test_mex_writer_parity.o
	$(CXX) -o $@ tests/test_mex_writer_parity.o $(LIB_TARGET) $(LIBSCRNA) $(LDFLAGS)

# Barcode filter test: tests emptydrops_bridge and barcode_filter modules
$(TEST_BARCODE_FILTER_TARGET): tests/test_barcode_filter.c $(LIB_TARGET) $(LIBSCRNA)
	$(CC) $(CFLAGS) -c tests/test_barcode_filter.c -o tests/test_barcode_filter.o
	$(CXX) -o $@ tests/test_barcode_filter.o $(LIB_TARGET) $(LIBSCRNA) $(LDFLAGS)

# Offset detection test: tests pf_api offset auto-detection and conflict handling
$(TEST_OFFSET_DETECTION_TARGET): tests/test_offset_detection.c $(LIB_TARGET) $(LIBSCRNA)
	$(CC) $(CFLAGS) -c tests/test_offset_detection.c -o tests/test_offset_detection.o
	$(CXX) -o $@ tests/test_offset_detection.o $(LIB_TARGET) $(LIBSCRNA) $(LDFLAGS)

# Regression test: compares assignBarcodes output against saved baseline
.PHONY: test-regression
test-regression: $(TARGET)
	@echo "Running assignBarcodes regression test..."
	@./tests/test_assignbarcodes_regression.sh

clean:
	rm -f *.o $(TARGET) $(DEMUX_TARGET) $(DEMUX_BAM_TARGET) $(CALL_FEATURES_TARGET) $(LIB_TARGET) $(TEST_API_TARGET) $(TEST_PF_COUNTS_TARGET) $(TEST_MEX_WRITER_TARGET) $(TEST_MEX_PARITY_TARGET) $(TEST_BARCODE_FILTER_TARGET) $(TEST_OFFSET_DETECTION_TARGET) tests/*.o
