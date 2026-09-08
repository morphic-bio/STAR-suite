CXX ?= g++
CXXFLAGS ?= -O3 -DNDEBUG -std=c++17 -Wall -Wextra -pedantic

.DEFAULT_GOAL := all

all: flex_bam_read_compare flex_mex_cell_stats

flex_bam_read_compare: flex_bam_read_compare.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

flex_mex_cell_stats: flex_mex_cell_stats.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -lz

test: flex_bam_read_compare flex_mex_cell_stats
	bash test_flex_bam_read_compare.sh ./flex_bam_read_compare
	bash test_flex_mex_cell_stats.sh ./flex_mex_cell_stats

test-flex-mex-cell-stats: flex_mex_cell_stats
	bash test_flex_mex_cell_stats.sh ./flex_mex_cell_stats

clean:
	rm -f flex_bam_read_compare flex_mex_cell_stats

.PHONY: all test test-flex-mex-cell-stats clean
