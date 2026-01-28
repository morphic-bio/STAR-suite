include build/common.mk
include build/core.mk
include build/flex.mk
include build/slam.mk
include build/tools.mk

.PHONY: all default build help clean

.DEFAULT_GOAL := default

# Build selector knobs (apply to default build only):
#   make
#   make default
#   make default INCLUDE="core flex-tools"
#   make default EXCLUDE="slam-tools yremove-tools"
DEFAULT_TARGETS := core flex-tools slam-tools vbem-tools yremove-tools feature-barcodes-tools
ALL_TARGETS := core flex-tools slam-tools vbem-tools yremove-tools feature-barcodes-tools process-features-tools star-feature-call

DEFAULT_BUILD := $(DEFAULT_TARGETS)
ifneq ($(strip $(INCLUDE)),)
DEFAULT_BUILD := $(INCLUDE)
endif
ifneq ($(strip $(EXCLUDE)),)
DEFAULT_BUILD := $(filter-out $(EXCLUDE),$(DEFAULT_BUILD))
endif

default: $(DEFAULT_BUILD)

# Build everything
all: $(ALL_TARGETS)

# Alias for convenience
build: default

help:
	@echo "STAR-suite build targets:"
	@echo "  make core            Build STAR core binary"
	@echo "  make flex            Build core + Flex tools"
	@echo "  make slam            Build core + SLAM tools"
	@echo "  make tools           Build all external tools"
	@echo "  make feature-barcodes-tools  Build feature barcode tools (assignBarcodes/demux)"
	@echo "  make                 Default build (core + common tools)"
	@echo "  make default         Same as make"
	@echo "  make default EXCLUDE=\"slam-tools\""
	@echo "  make default INCLUDE=\"core flex-tools\""
	@echo "  make all             Build everything"
	@echo "  make build           Alias for make default"
	@echo "  make clean           Clean core and tool outputs"

clean: core-clean tools-clean
