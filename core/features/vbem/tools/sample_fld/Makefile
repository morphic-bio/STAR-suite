CC = gcc
CFLAGS = -Wall -O2 -g
LDFLAGS = -lhts -lz -lpthread

# Find htslib
HTSLIB_DIR = ../../source/htslib
HTSLIB_INC = -I$(HTSLIB_DIR)
HTSLIB_LIB = -L$(HTSLIB_DIR) -lhts

sample_fld: sample_fld.c
	$(CC) $(CFLAGS) $(HTSLIB_INC) -o sample_fld sample_fld.c $(HTSLIB_LIB) $(LDFLAGS)

clean:
	rm -f sample_fld

install: sample_fld
	cp sample_fld ../../bin/

.PHONY: clean install
