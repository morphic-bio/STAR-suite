# Makefile for remove_y_reads standalone tool

ROOT_DIR := $(abspath ../../../../..)
HTSLIB_DIR = $(ROOT_DIR)/core/legacy/source/htslib
HTSLIB = $(HTSLIB_DIR)/libhts.a
CC = gcc
CFLAGS = -O3 -Wall -Wextra -pthread -I$(HTSLIB_DIR) -I$(HTSLIB_DIR)/htslib
LDFLAGS = -L$(HTSLIB_DIR) -Bstatic -lhts -Bdynamic -lz -lpthread

TARGET = remove_y_reads
SRC = remove_y_reads.c

$(TARGET): $(SRC) $(HTSLIB)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

clean:
	rm -f $(TARGET) $(TARGET).o

.PHONY: clean

$(HTSLIB):
	$(MAKE) -C $(HTSLIB_DIR) lib-static
