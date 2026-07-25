.PHONY: flex flex-tools flexfilter molecule-first-resolver

flex: flex-tools

flex-tools: flexfilter molecule-first-resolver

flexfilter: core
	$(MAKE) -C $(FLEX_DIR)/tools/flexfilter

molecule-first-resolver:
	$(MAKE) -C $(FLEX_DIR)/tools/molecule_first_resolver
