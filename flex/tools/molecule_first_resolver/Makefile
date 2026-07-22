ROOT_DIR := $(abspath ../../..)
FLEX_SRC := $(ROOT_DIR)/flex/source
CORE_SRC := $(ROOT_DIR)/core/legacy/source
HTSLIB_DIR := $(ROOT_DIR)/core/legacy/source/htslib
HTSLIB := $(HTSLIB_DIR)/libhts.a

CXX ?= g++
CXXFLAGS ?= -O3 -std=c++11 -Wall -Wextra
VERSION_VALUE := $(shell sed -n 's/^\#define STAR_SUITE_VERSION "\([^"]*\)"/\1/p' $(ROOT_DIR)/core/legacy/source/VERSION)
CPPFLAGS += -I$(FLEX_SRC) -I$(CORE_SRC) -I$(HTSLIB_DIR) -DSTAR_SUITE_VERSION=\"$(VERSION_VALUE)\"

RESOLVER_SRC := $(FLEX_SRC)/solo/MoleculeFirstResolver.cpp
MULTIGENE_SRC := $(CORE_SRC)/MultiGeneUmiCr.cpp
SIDECAR_SRC := $(CORE_SRC)/SpatialFeatureSidecar.cpp
TARGETS := molecule_first_resolver molecule_first_bam_ledger molecule_first_materialize \
	molecule_first_gex_reconcile spatial_feature_sidecar_dump spatial_feature_sidecar_join
TEST_TARGET := test_molecule_first_resolver
SIDECAR_JOIN_TEST := test_spatial_feature_sidecar_join
GEX_RECONCILE_CLI_TEST := test_gex_reconcile_cli.sh

.DEFAULT_GOAL := all

all: $(TARGETS)

molecule_first_resolver: molecule_first_resolver.cpp $(RESOLVER_SRC) $(MULTIGENE_SRC) $(FLEX_SRC)/solo/MoleculeFirstResolver.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ molecule_first_resolver.cpp $(RESOLVER_SRC) $(MULTIGENE_SRC)

molecule_first_bam_ledger: molecule_first_bam_ledger.cpp $(HTSLIB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ molecule_first_bam_ledger.cpp \
		$(HTSLIB) -lz -lssl -lcrypto -ldl -lpthread

molecule_first_materialize: molecule_first_materialize.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ molecule_first_materialize.cpp

molecule_first_gex_reconcile: molecule_first_gex_reconcile.cpp $(MULTIGENE_SRC)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ molecule_first_gex_reconcile.cpp $(MULTIGENE_SRC)

spatial_feature_sidecar_dump: spatial_feature_sidecar_dump.cpp $(SIDECAR_SRC) $(CORE_SRC)/SpatialFeatureSidecar.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ spatial_feature_sidecar_dump.cpp $(SIDECAR_SRC) -lcrypto

spatial_feature_sidecar_join: spatial_feature_sidecar_join.cpp $(SIDECAR_SRC) $(CORE_SRC)/SpatialFeatureSidecar.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ spatial_feature_sidecar_join.cpp $(SIDECAR_SRC) -lz -lcrypto

$(HTSLIB):
	$(MAKE) -C $(HTSLIB_DIR) lib-static

$(TEST_TARGET): test_molecule_first_resolver.cpp $(RESOLVER_SRC) $(MULTIGENE_SRC) $(FLEX_SRC)/solo/MoleculeFirstResolver.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ test_molecule_first_resolver.cpp $(RESOLVER_SRC) $(MULTIGENE_SRC)

$(SIDECAR_JOIN_TEST): test_spatial_feature_sidecar_join.cpp $(SIDECAR_SRC) spatial_feature_sidecar_join
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ test_spatial_feature_sidecar_join.cpp $(SIDECAR_SRC) -lz -lcrypto

test: $(TEST_TARGET) $(SIDECAR_JOIN_TEST) $(TARGETS)
	./$(TEST_TARGET)
	./$(SIDECAR_JOIN_TEST) ./spatial_feature_sidecar_join
	./test_cli.sh ./molecule_first_resolver
	./$(GEX_RECONCILE_CLI_TEST) ./molecule_first_gex_reconcile

clean:
	rm -f $(TARGETS) $(TEST_TARGET) $(SIDECAR_JOIN_TEST)

.PHONY: all test clean
