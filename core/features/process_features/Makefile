CC=gcc
AR=ar

# Directories
SRCDIR = src
INCDIR = include
VPATH = $(SRCDIR)

# Use pkg-config to get the correct flags for cairo
HTS_LIBS := -lhts

CAIRO_CFLAGS :=
CAIRO_LIBS :=
DEFINES :=

# Source files (basenames only)
# Core library sources (shared by all tools)
# Note: plot_histogram.c is needed by assignBarcodes.c for histogram generation
CORE_SRCS = assignBarcodes.c queue.c globals.c utils.c memory.c io.c barcode_match.c pf_api.c call_features.c gmm.c mex_10x.c plot_histogram.c

# CLI-specific sources
MAIN_SRCS = main.c
DEMUX_FASTQ_SRCS = demux_fastq.c
DEMUX_BAM_SRCS = demux_bam.c

# Full source lists for each target
SRCS_BASE = $(MAIN_SRCS) $(CORE_SRCS)
SRCS = $(SRCS_BASE)

# Add cairo flags and define if NO_HEATMAP is not set to 1
ifeq ($(NO_HEATMAP), 1)
	DEFINES += -DNO_HEATMAP
else
	CAIRO_CFLAGS := $(shell pkg-config --cflags cairo)
	CAIRO_LIBS := $(shell pkg-config --libs cairo)
	SRCS += heatmap.c
	CORE_SRCS += heatmap.c
endif

# Compiler flags: optimise by default, but if DEBUG=1 build with symbols and no optimisations
ifeq ($(DEBUG),1)
    CFLAGS=-g -ggdb -Wall -O0 -I$(INCDIR) -fopenmp $(CAIRO_CFLAGS) $(DEFINES)
else
    CFLAGS=-g -Wall -O3 -I$(INCDIR) -fopenmp $(CAIRO_CFLAGS) $(DEFINES)
endif
LDFLAGS=-lm -lpthread -lz -fopenmp $(CAIRO_LIBS) $(HTS_LIBS)

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

# Test program
TEST_API_TARGET=tests/test_pf_api

.PHONY: all clean lib tools test

all: lib tools

test: $(TEST_API_TARGET)

lib: $(LIB_TARGET)

tools: $(TARGET) $(DEMUX_TARGET) $(DEMUX_BAM_TARGET) $(CALL_FEATURES_TARGET)

$(LIB_TARGET): $(CORE_OBJS)
	$(AR) rcs $@ $^

$(TARGET): main.o $(LIB_TARGET)
	$(CC) -o $@ main.o $(LIB_TARGET) $(LDFLAGS)

$(DEMUX_TARGET): $(DEMUX_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(DEMUX_BAM_TARGET): $(DEMUX_BAM_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(CALL_FEATURES_TARGET): call_features_main.o $(LIB_TARGET)
	$(CC) -o $@ call_features_main.o $(LIB_TARGET) $(LDFLAGS)

# plot_histogram.o is needed by assignBarcodes but not in core lib
plot_histogram.o: plot_histogram.c
	$(CC) $(CFLAGS) -c $< -o $@

# The implicit rule will now use VPATH to find the .c files in src/
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_API_TARGET): tests/test_pf_api.c $(LIB_TARGET)
	$(CC) $(CFLAGS) -o $@ $< $(LIB_TARGET) $(LDFLAGS)

clean:
	rm -f *.o $(TARGET) $(DEMUX_TARGET) $(DEMUX_BAM_TARGET) $(CALL_FEATURES_TARGET) $(LIB_TARGET) $(TEST_API_TARGET)
