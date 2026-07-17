ROOT_DIR := $(abspath ../../..)
FLEX_SRC := $(ROOT_DIR)/flex/source
HTSLIB_DIR := $(ROOT_DIR)/core/legacy/source/htslib
HTSLIB := $(HTSLIB_DIR)/libhts.a

CXX ?= g++
CXXFLAGS ?= -O3 -std=c++11 -Wall -Wextra
VERSION_VALUE := $(shell sed -n 's/^\#define STAR_SUITE_VERSION "\([^"]*\)"/\1/p' $(ROOT_DIR)/core/legacy/source/VERSION)
CPPFLAGS += -I$(FLEX_SRC) -I$(HTSLIB_DIR) -DSTAR_SUITE_VERSION=\"$(VERSION_VALUE)\"

RESOLVER_SRC := $(FLEX_SRC)/solo/MoleculeFirstResolver.cpp
TARGETS := molecule_first_resolver molecule_first_bam_ledger molecule_first_materialize
TEST_TARGET := test_molecule_first_resolver

.DEFAULT_GOAL := all

all: $(TARGETS)

molecule_first_resolver: molecule_first_resolver.cpp $(RESOLVER_SRC) $(FLEX_SRC)/solo/MoleculeFirstResolver.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ molecule_first_resolver.cpp $(RESOLVER_SRC)

molecule_first_bam_ledger: molecule_first_bam_ledger.cpp $(HTSLIB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ molecule_first_bam_ledger.cpp \
		$(HTSLIB) -lz -lssl -lcrypto -ldl -lpthread

molecule_first_materialize: molecule_first_materialize.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ molecule_first_materialize.cpp

$(HTSLIB):
	$(MAKE) -C $(HTSLIB_DIR) lib-static

$(TEST_TARGET): test_molecule_first_resolver.cpp $(RESOLVER_SRC) $(FLEX_SRC)/solo/MoleculeFirstResolver.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ test_molecule_first_resolver.cpp $(RESOLVER_SRC)

test: $(TEST_TARGET) $(TARGETS)
	./$(TEST_TARGET)
	./test_cli.sh ./molecule_first_resolver

clean:
	rm -f $(TARGETS) $(TEST_TARGET)

.PHONY: all test clean
