ROOT_DIR := $(abspath ../../..)
FLEX_SRC := $(ROOT_DIR)/flex/source

CXX ?= g++
CXXFLAGS ?= -O3 -std=c++11 -Wall -Wextra
VERSION_VALUE := $(shell sed -n 's/^\#define STAR_SUITE_VERSION "\([^"]*\)"/\1/p' $(ROOT_DIR)/core/legacy/source/VERSION)
CPPFLAGS += -I$(FLEX_SRC) -DSTAR_SUITE_VERSION=\"$(VERSION_VALUE)\"

RESOLVER_SRC := $(FLEX_SRC)/solo/MoleculeFirstResolver.cpp
TARGET := molecule_first_resolver
TEST_TARGET := test_molecule_first_resolver

.DEFAULT_GOAL := all

all: $(TARGET)

$(TARGET): molecule_first_resolver.cpp $(RESOLVER_SRC) $(FLEX_SRC)/solo/MoleculeFirstResolver.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ molecule_first_resolver.cpp $(RESOLVER_SRC)

$(TEST_TARGET): test_molecule_first_resolver.cpp $(RESOLVER_SRC) $(FLEX_SRC)/solo/MoleculeFirstResolver.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ test_molecule_first_resolver.cpp $(RESOLVER_SRC)

test: $(TEST_TARGET) $(TARGET)
	./$(TEST_TARGET)
	./test_cli.sh ./$(TARGET)

clean:
	rm -f $(TARGET) $(TEST_TARGET)

.PHONY: all test clean
