include build/common.mk
include build/core.mk
include build/flex.mk
include build/slam.mk
include build/tools.mk

.PHONY: all help clean

all: core flex-tools slam-tools vbem-tools yremove-tools feature-barcodes-tools

help:
	@echo "STAR-suite build targets:"
	@echo "  make core            Build STAR core binary"
	@echo "  make flex            Build core + Flex tools"
	@echo "  make slam            Build core + SLAM tools"
	@echo "  make tools           Build all external tools"
	@echo "  make feature-barcodes-tools  Build feature barcode tools (assignBarcodes/demux)"
	@echo "  make all             Build core + all tools"
	@echo "  make clean           Clean core and tool outputs"

clean: core-clean tools-clean
