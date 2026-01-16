# Makefile for libtrim - Cutadapt-parity trimming library

CXX ?= g++
CXXFLAGS_common := -std=c++11 -O3 -Wall -Wextra
CXXFLAGS ?= -pipe
AR ?= ar

LIBTRIM_OBJECTS = trim.o quality_trim.o adapter_trim.o

all: libtrim.a

libtrim.a: $(LIBTRIM_OBJECTS)
	$(AR) rcs $@ $^
	ranlib $@

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(CXXFLAGS_common) -I.. $< -o $@

clean:
	rm -f $(LIBTRIM_OBJECTS) libtrim.a

.PHONY: all clean
