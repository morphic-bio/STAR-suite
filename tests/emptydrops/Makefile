# Makefile for EmptyDrops unit tests
# Per emptydrops_refactor_plan.md

CXX := g++
CXXFLAGS := -std=c++11 -O3 -Wall -Wextra

# libscrna paths
LIBSCRNA_DIR := ../../core/features/libscrna
LIBSCRNA_INC := $(LIBSCRNA_DIR)/include
LIBSCRNA := $(LIBSCRNA_DIR)/libscrna.a

INCLUDES := -I$(LIBSCRNA_INC)
LDFLAGS := -L$(LIBSCRNA_DIR) -lscrna -lpthread

# Test executables
TESTS := test_sampler_deterministic \
         test_ambient_profile \
         test_simple_ed_threshold \
         test_ordmag_ties \
         test_ordmag_floor \
         test_ordmag_parallel \
         test_ordmag_estimator_parity \
         test_thread_permits \
         test_no_tag_mode \
         test_observed_occupancy

.PHONY: all clean test libscrna

all: $(TESTS)

# Build libscrna if needed
$(LIBSCRNA):
	$(MAKE) -C $(LIBSCRNA_DIR) lib

libscrna: $(LIBSCRNA)

# Build test executables
test_sampler_deterministic: test_sampler_deterministic.cpp $(LIBSCRNA)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

test_ambient_profile: test_ambient_profile.cpp $(LIBSCRNA)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

test_simple_ed_threshold: test_simple_ed_threshold.cpp $(LIBSCRNA)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

test_ordmag_ties: test_ordmag_ties.cpp $(LIBSCRNA_INC)/OrdMagRank.h $(LIBSCRNA)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

test_ordmag_floor: test_ordmag_floor.cpp $(LIBSCRNA)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

test_ordmag_parallel: test_ordmag_parallel.cpp $(LIBSCRNA_INC)/ParallelTasks.h $(LIBSCRNA)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

test_ordmag_estimator_parity: test_ordmag_estimator_parity.cpp $(LIBSCRNA)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

test_thread_permits: test_thread_permits.cpp $(LIBSCRNA_INC)/ParallelTasks.h $(LIBSCRNA_INC)/ThreadPermits.h $(LIBSCRNA)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

test_no_tag_mode: test_no_tag_mode.cpp $(LIBSCRNA)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

test_observed_occupancy: test_observed_occupancy.cpp $(LIBSCRNA_INC)/ObservedTagOccupancy.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $<

# Run all tests
test: $(TESTS)
	@echo "=== EmptyDrops Unit Tests ==="
	@failed=0; \
	for t in $(TESTS); do \
		echo ""; \
		echo "Running $$t..."; \
		./$$t || { echo "FAILED: $$t"; failed=1; }; \
	done; \
	echo ""; \
	if [ $$failed -eq 0 ]; then \
		echo "All tests PASSED"; \
	else \
		echo "Some tests FAILED"; \
		exit 1; \
	fi

clean:
	rm -f $(TESTS) test_flex_tag_aware

# Integration fixture: run explicitly with a fresh output directory, then compare
# its saved internal results against scrna_simpleed using check_flex_tag_aware.py.
test_flex_tag_aware: test_flex_tag_aware.cpp $(LIBSCRNA) ../../flex/source/libflex/libflex.a ../../core/legacy/source/MexWriter.o
	$(CXX) $(CXXFLAGS) -fopenmp $(INCLUDES) -I../../flex/source -I../../core/legacy/source -o $@ $< ../../flex/source/libflex/libflex.a ../../core/legacy/source/MexWriter.o $(LDFLAGS)
