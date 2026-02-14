# Makefile for remove_y_reads standalone tool

ROOT_DIR := $(abspath ../../../../..)
HTSLIB_DIR = $(ROOT_DIR)/core/legacy/source/htslib
HTSLIB = $(HTSLIB_DIR)/libhts.a
CC ?= cc
CFLAGS ?= -O3 -Wall -Wextra -pthread
CPPFLAGS ?=
LOCAL_CPPFLAGS = -I$(HTSLIB_DIR) -I$(HTSLIB_DIR)/htslib
LOCAL_LDLIBS = $(HTSLIB) -lz -lpthread

TARGET = remove_y_reads
SRC = remove_y_reads.c

$(TARGET): $(SRC) $(HTSLIB)
	$(CC) $(CPPFLAGS) $(LOCAL_CPPFLAGS) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS) $(LOCAL_LDLIBS)

clean:
	rm -f $(TARGET) $(TARGET).o

.PHONY: clean

$(HTSLIB):
	$(MAKE) -C $(HTSLIB_DIR) lib-static
