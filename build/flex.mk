.PHONY: flex flex-tools flexfilter

flex: flex-tools

flex-tools: flexfilter

flexfilter: core
	$(MAKE) -C $(FLEX_DIR)/tools/flexfilter
