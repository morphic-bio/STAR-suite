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

# Release branch helper knobs. These targets are intentionally explicit and do
# not push; pushing release branches/tags stays a human-controlled step.
RELEASE_VERSION ?=
RC ?= 1
DEV_RELEASE_BRANCH ?= dev-release
DEV_RELEASE_BASE ?= master
DEV_RELEASE_VERSION_BRANCH = $(DEV_RELEASE_BRANCH)-v$(RELEASE_VERSION)
RC_TAG = v$(RELEASE_VERSION)-rc$(RC)
DEV_RELEASE_CHECK_SCRIPTS := tests/run_production_module_regression_suite.sh tests/run_star_chromap_macs3_lowmem_smoke_100k.sh

default: $(DEFAULT_BUILD)

# Build everything
all: $(ALL_TARGETS)

# Alias for convenience
build: default

help:
	@echo "STAR-suite build targets:"
	@echo "  make core            Build Chromap-enabled STAR core binary"
	@echo "  make core-portable   Build STAR core without libchromap"
	@echo "  make flex            Build core + Flex tools"
	@echo "  make slam            Build core + SLAM tools"
	@echo "  make tools           Build all external tools"
	@echo "  make feature-barcodes-tools  Build feature barcode tools (assignBarcodes/demux)"
	@echo "  make star-libchromap-contract Build STAR-owned libchromap contract runner"
	@echo "  make                 Default build (core + common tools)"
	@echo "  make default         Same as make"
	@echo "  make default EXCLUDE=\"slam-tools\""
	@echo "  make default INCLUDE=\"core flex-tools\""
	@echo "  make all             Build everything"
	@echo "  make build           Alias for make default"
	@echo "  make clean           Clean core and tool outputs"
	@echo "  make dev-release-help"
	@echo "                       Show dev-release branch and RC helper commands"
	@echo "  make dev-release-branch RELEASE_VERSION=X.Y.Z"
	@echo "                       Create local dev-release and dev-release-vX.Y.Z branches"
	@echo "  make dev-release-check"
	@echo "                       Run local dev-release build/preflight checks"
	@echo "  make dev-release-tag RELEASE_VERSION=X.Y.Z RC=N"
	@echo "                       Create local annotated vX.Y.Z-rcN tag"

clean: core-clean tools-clean

.PHONY: require-release-version dev-release-help dev-release-branch dev-release-check dev-release-tag

require-release-version:
	@if [ -z "$(strip $(RELEASE_VERSION))" ]; then \
		echo "ERROR: set RELEASE_VERSION, for example: make $@ RELEASE_VERSION=1.1.0" >&2; \
		exit 2; \
	fi

dev-release-help:
	@echo "Dev-release branch workflow:"
	@echo "  make dev-release-branch RELEASE_VERSION=1.1.0"
	@echo "  # Optional: add DEV_RELEASE_BASE=<branch-or-sha> to seed dev-release"
	@echo "  git switch dev-release-v1.1.0"
	@echo "  make dev-release-check"
	@echo "  git tag -a v1.1.0-rc1 -m 'STAR Suite v1.1.0-rc1'"
	@echo "  git push origin dev-release dev-release-v1.1.0"
	@echo "  git push origin v1.1.0-rc1"
	@echo
	@echo "Stable promotion after RC validation:"
	@echo "  git switch master && git merge --no-ff dev-release-v1.1.0"
	@echo "  git tag -a v1.1.0 -m 'STAR Suite v1.1.0'"

dev-release-branch: require-release-version
	@if git rev-parse --verify "$(DEV_RELEASE_BRANCH)" >/dev/null 2>&1; then \
		echo "Branch exists: $(DEV_RELEASE_BRANCH)"; \
	else \
		git branch "$(DEV_RELEASE_BRANCH)" "$(DEV_RELEASE_BASE)" && \
		echo "Created branch: $(DEV_RELEASE_BRANCH) from $(DEV_RELEASE_BASE)"; \
	fi
	@if git rev-parse --verify "$(DEV_RELEASE_VERSION_BRANCH)" >/dev/null 2>&1; then \
		echo "Branch exists: $(DEV_RELEASE_VERSION_BRANCH)"; \
	else \
		git branch "$(DEV_RELEASE_VERSION_BRANCH)" HEAD && \
		echo "Created branch: $(DEV_RELEASE_VERSION_BRANCH) from HEAD"; \
	fi
	@echo "Push when ready: git push origin $(DEV_RELEASE_BRANCH) $(DEV_RELEASE_VERSION_BRANCH)"

dev-release-check:
	$(MAKE) core
	@for script in $(DEV_RELEASE_CHECK_SCRIPTS); do \
		echo "Syntax check $$script"; \
		bash -n "$$script"; \
	done
	tests/run_production_module_regression_suite.sh --preflight --scale production-100k

dev-release-tag: require-release-version
	@if ! git diff --quiet || ! git diff --cached --quiet; then \
		echo "ERROR: working tree has uncommitted changes; commit before tagging $(RC_TAG)." >&2; \
		exit 2; \
	fi
	@case "$$(git branch --show-current)" in \
		$(DEV_RELEASE_BRANCH)|$(DEV_RELEASE_BRANCH)-v*) ;; \
		*) echo "ERROR: create RC tags from $(DEV_RELEASE_BRANCH) or $(DEV_RELEASE_BRANCH)-vX.Y.Z." >&2; exit 2 ;; \
	esac
	@if git rev-parse --verify "refs/tags/$(RC_TAG)" >/dev/null 2>&1; then \
		echo "ERROR: tag already exists: $(RC_TAG)" >&2; \
		exit 2; \
	fi
	git tag -a "$(RC_TAG)" -m "STAR Suite $(RC_TAG)"
	@echo "Created local tag $(RC_TAG). Push when ready: git push origin $(RC_TAG)"
