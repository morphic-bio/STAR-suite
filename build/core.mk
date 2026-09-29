.PHONY: core core-star core-portable core-long core-static core-htslib core-clean star-host-lib host-api-tests

core: core-star

core-star:
	$(MAKE) -C $(LEGACY_SRC_DIR) STAR

# Compatibility alias: since 1.10.0 `make core` is already the portable build
# (STAR's bundled HTSlib, no other suite on any path).
core-portable: core-star

core-long:
	$(MAKE) -C $(LEGACY_SRC_DIR) STARlong

core-static:
	$(MAKE) -C $(LEGACY_SRC_DIR) STARstatic

core-htslib:
	$(MAKE) -C $(LEGACY_SRC_DIR)/htslib lib-static

# Host library libstar_suite.a + libstar_suite.link (docs/HOST_API.md).
# Pass HTSLIB=external when the host links other HTSlib users.
star-host-lib:
	$(MAKE) -C $(LEGACY_SRC_DIR) libstar_suite

# Host interface tests (tests/host_api). One submake builds STAR and the
# host library so the two never race in the same tree.
host-api-tests:
	$(MAKE) -C $(LEGACY_SRC_DIR) STAR libstar_suite
	bash $(ROOT_DIR)/tests/host_api/run_host_api_tests.sh

core-clean:
	$(MAKE) -C $(LEGACY_SRC_DIR) clean
