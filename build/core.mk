.PHONY: core core-star core-portable core-long core-static core-htslib core-clean

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

core-clean:
	$(MAKE) -C $(LEGACY_SRC_DIR) clean
