.PHONY: core core-star core-long core-static core-htslib core-clean

core: core-star

core-star:
	$(MAKE) -C $(LEGACY_SRC_DIR) STAR

core-long:
	$(MAKE) -C $(LEGACY_SRC_DIR) STARlong

core-static:
	$(MAKE) -C $(LEGACY_SRC_DIR) STARstatic

core-htslib:
	$(MAKE) -C $(LEGACY_SRC_DIR)/htslib lib-static

core-clean:
	$(MAKE) -C $(LEGACY_SRC_DIR) clean
