CXX ?= g++
AR ?= ar

CHROMAP_DIR ?= /mnt/pikachu/Chromap-suite
CHROMAP_LIB ?= $(CHROMAP_DIR)/libchromap.a
CHROMAP_SRC_DIR := $(CHROMAP_DIR)/src
CHROMAP_HTSLIB_DIR := $(CHROMAP_DIR)/third_party/htslib
LIBMACS3_INCLUDE_DIR ?= $(CHROMAP_DIR)/third_party/libMACS3/include
LIBMACS3_LIB ?= $(CHROMAP_DIR)/third_party/libMACS3/lib/libmacs3.a

# libscrna is a contract dep: the contract folds the per-barcode ATAC
# evidence-from-peaks step in-memory after libchromap's MACS3 FRAG step,
# using libscrna's compact row-stream evidence API.
LIBSCRNA_DIR ?= $(abspath $(CURDIR)/../libscrna)
LIBSCRNA_INCLUDE_DIR ?= $(LIBSCRNA_DIR)/include
LIBSCRNA_LIB ?= $(LIBSCRNA_DIR)/libscrna.a

SRC_DIR := src
INC_DIR := include
TOOLS_DIR := tools
OBJ_DIR := objs

LIB := libstar_chromap_contract.a
RUNNER := star_libchromap_contract_runner

CXXFLAGS ?= -std=c++11 -Wall -Wextra -O3 -fopenmp
CXXFLAGS += -I$(INC_DIR) -I$(CHROMAP_SRC_DIR) -I$(CHROMAP_HTSLIB_DIR) -I$(LIBMACS3_INCLUDE_DIR) -I$(LIBSCRNA_INCLUDE_DIR)
DEPFLAGS ?= -MMD -MP
LDFLAGS ?=
CHROMAP_LDFLAGS := -L$(CHROMAP_HTSLIB_DIR) -lhts -lm -lz -lpthread -lcurl -lcrypto -lbz2 -llzma -ldeflate

LIB_OBJS := $(OBJ_DIR)/star_chromap_contract.o
RUNNER_OBJS := $(OBJ_DIR)/star_libchromap_contract_runner.o
DEPS := $(LIB_OBJS:.o=.d) $(RUNNER_OBJS:.o=.d)

.PHONY: all clean chromap-lib

all: $(LIB) $(RUNNER)

chromap-lib:
	$(MAKE) -C $(CHROMAP_DIR) libchromap.a

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/star_chromap_contract.o: $(SRC_DIR)/star_chromap_contract.cpp $(INC_DIR)/star_chromap_contract.h | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

$(OBJ_DIR)/star_libchromap_contract_runner.o: $(TOOLS_DIR)/star_libchromap_contract_runner.cpp $(INC_DIR)/star_chromap_contract.h | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

$(LIB): $(LIB_OBJS)
	$(AR) rcs $@ $^

$(RUNNER): chromap-lib $(LIB) $(RUNNER_OBJS)
	$(CXX) $(CXXFLAGS) $(RUNNER_OBJS) $(LIB) $(CHROMAP_LIB) $(LIBMACS3_LIB) $(LIBSCRNA_LIB) -o $@ $(LDFLAGS) $(CHROMAP_LDFLAGS)

clean:
	rm -rf $(OBJ_DIR) $(LIB) $(RUNNER)

-include $(DEPS)
