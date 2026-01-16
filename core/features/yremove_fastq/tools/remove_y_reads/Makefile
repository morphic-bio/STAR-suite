# Makefile for remove_y_reads standalone tool

HTSLIB_DIR = ../../source/htslib
CC = gcc
CFLAGS = -O3 -Wall -Wextra -pthread -I$(HTSLIB_DIR) -I$(HTSLIB_DIR)/htslib
LDFLAGS = -L$(HTSLIB_DIR) -Bstatic -lhts -Bdynamic -lz -lpthread

TARGET = remove_y_reads
SRC = remove_y_reads.c

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

clean:
	rm -f $(TARGET) $(TARGET).o

.PHONY: clean

