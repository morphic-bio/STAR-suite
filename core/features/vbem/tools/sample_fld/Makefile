CC = gcc
CFLAGS = -Wall -O2 -g
LDFLAGS = -lhts -lz -lpthread

# Find htslib
ROOT_DIR := $(abspath ../../../../..)
HTSLIB_DIR = $(ROOT_DIR)/core/legacy/source/htslib
HTSLIB_INC = -I$(HTSLIB_DIR)
HTSLIB_LIB = -L$(HTSLIB_DIR) -lhts
HTSLIB = $(HTSLIB_DIR)/libhts.a

sample_fld: sample_fld.c $(HTSLIB)
	$(CC) $(CFLAGS) $(HTSLIB_INC) -o sample_fld sample_fld.c $(HTSLIB_LIB) $(LDFLAGS)

clean:
	rm -f sample_fld

install: sample_fld
	cp sample_fld $(ROOT_DIR)/bin/

.PHONY: clean install

$(HTSLIB):
	$(MAKE) -C $(HTSLIB_DIR) lib-static
