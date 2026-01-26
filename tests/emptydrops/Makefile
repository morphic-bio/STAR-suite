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
         test_no_tag_mode

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

test_no_tag_mode: test_no_tag_mode.cpp $(LIBSCRNA)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

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
	rm -f $(TESTS)
