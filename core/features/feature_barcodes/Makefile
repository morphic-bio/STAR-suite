CC=gcc

# Directories
SRCDIR = src
INCDIR = include
VPATH = $(SRCDIR)

# Use pkg-config to get the correct flags for glib-2.0 and cairo
GLIB_CFLAGS := $(shell pkg-config --cflags glib-2.0)
GLIB_LIBS := $(shell pkg-config --libs glib-2.0)
HTS_LIBS := -lhts

CAIRO_CFLAGS :=
CAIRO_LIBS :=
DEFINES :=

# Source files (basenames only)
SRCS_BASE=main.c assignBarcodes.c queue.c globals.c utils.c memory.c io.c plot_histogram.c barcode_match.c
SRCS = $(SRCS_BASE)

# Add cairo flags and define if NO_HEATMAP is not set to 1
ifeq ($(NO_HEATMAP), 1)
	DEFINES += -DNO_HEATMAP
else
	CAIRO_CFLAGS := $(shell pkg-config --cflags cairo)
	CAIRO_LIBS := $(shell pkg-config --libs cairo)
	SRCS += heatmap.c
endif

# Compiler flags: optimise by default, but if DEBUG=1 build with symbols and no optimisations
ifeq ($(DEBUG),1)
    CFLAGS=-g -ggdb -Wall -O0 -I$(INCDIR) -fopenmp $(GLIB_CFLAGS) $(CAIRO_CFLAGS) $(DEFINES)
else
    CFLAGS=-g -Wall -O3 -I$(INCDIR) -fopenmp $(GLIB_CFLAGS) $(CAIRO_CFLAGS) $(DEFINES)
endif
LDFLAGS=-lm -lpthread -lz -fopenmp $(GLIB_LIBS) $(CAIRO_LIBS) $(HTS_LIBS)

# Object files
OBJS=$(SRCS:.c=.o)

# Executable names
TARGET=assignBarcodes
DEMUX_TARGET=demux_fastq

# demux_fastq sources/objects
DEMUX_SRCS=demux_fastq.c io.c utils.c globals.c barcode_match.c
DEMUX_OBJS=$(DEMUX_SRCS:.c=.o)

DEMUX_BAM_TARGET=demux_bam
DEMUX_BAM_SRCS=demux_bam.c utils.c globals.c barcode_match.c io.c memory.c
DEMUX_BAM_OBJS=$(DEMUX_BAM_SRCS:.c=.o)

.PHONY: all clean

all: $(TARGET) $(DEMUX_TARGET) $(DEMUX_BAM_TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(DEMUX_TARGET): $(DEMUX_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(DEMUX_BAM_TARGET): $(DEMUX_BAM_OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

# The implicit rule will now use VPATH to find the .c files in src/
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET) $(DEMUX_OBJS) $(DEMUX_TARGET) $(DEMUX_BAM_OBJS) $(DEMUX_BAM_TARGET)