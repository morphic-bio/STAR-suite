# Makefile for EC Filter Test Suite

CXX ?= g++
CXXFLAGS := -std=c++11 -O3 -Wall -Wextra -pthread
LDFLAGS :=

SCRIPT_DIR := $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST))))
LIBEM := $(SCRIPT_DIR)/../../source/libem/libem.a
HTSLIB_DIR := $(SCRIPT_DIR)/../../source/htslib

all: ec_filter_cli

ec_filter_cli: ec_filter_cli.cpp $(LIBEM) $(HTSLIB_DIR)/libhts.a
	$(CXX) $(CXXFLAGS) -I$(SCRIPT_DIR)/../../source -I$(HTSLIB_DIR) ec_filter_cli.cpp $(LIBEM) -L$(HTSLIB_DIR) -lhts -lz -lpthread $(LDFLAGS) -o ec_filter_cli

$(HTSLIB_DIR)/libhts.a:
	$(MAKE) -C $(HTSLIB_DIR) lib-static

$(LIBEM):
	$(MAKE) -C $(SCRIPT_DIR)/../../source/libem

test: ec_filter_cli
	@echo "Running unit tests..."
	@python3 generate_fixtures.py > fixtures.json
	@echo "Fixtures generated. Run ec_filter_cli with fixtures.json to generate results."
	@echo "Then run: python3 verify_filtering.py results.json"

parity-test: ec_filter_cli
	@echo "Running Salmon parity test (synthetic data)..."
	@./run_salmon_parity.sh /tmp/ec_parity_test
	@if grep -q "STATUS: PASS" /tmp/ec_parity_test/parity_report.txt 2>/dev/null; then \
		echo "Parity test PASSED"; \
	elif grep -q "STATUS: WARN" /tmp/ec_parity_test/parity_report.txt 2>/dev/null; then \
		echo "Parity test WARNED"; \
	else \
		echo "Parity test FAILED or report not found"; \
		exit 1; \
	fi

nfcore-parity-test: ec_filter_cli
	@echo "Running Salmon parity test (nf-core/test-datasets)..."
	@./run_nfcore_parity.sh /tmp/nfcore_ec_parity_test
	@if grep -q "STATUS: PASS" /tmp/nfcore_ec_parity_test/parity_report.txt 2>/dev/null; then \
		echo "nf-core parity test PASSED"; \
	elif grep -q "STATUS: WARN" /tmp/nfcore_ec_parity_test/parity_report.txt 2>/dev/null; then \
		echo "nf-core parity test WARNED"; \
	else \
		echo "nf-core parity test FAILED or report not found (CLI may not be implemented yet)"; \
		exit 1; \
	fi

clean:
	rm -f ec_filter_cli fixtures.json results.json

.PHONY: all test parity-test nfcore-parity-test clean
