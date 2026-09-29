ROOT_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/..)

CORE_DIR := $(ROOT_DIR)/core
LEGACY_DIR := $(CORE_DIR)/legacy
LEGACY_SRC_DIR := $(LEGACY_DIR)/source

FLEX_DIR := $(ROOT_DIR)/flex
SLAM_DIR := $(ROOT_DIR)/slam

VBEM_DIR := $(CORE_DIR)/features/vbem
YREMOVE_FASTQ_DIR := $(CORE_DIR)/features/yremove_fastq
BAMSORT_DIR := $(CORE_DIR)/features/bamsort
FEATURE_BARCODES_DIR := $(CORE_DIR)/features/feature_barcodes
PROCESS_FEATURES_DIR := $(CORE_DIR)/features/process_features
