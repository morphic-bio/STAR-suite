.PHONY: flex flex-tools flexfilter

flex: core flex-tools

flex-tools: core flexfilter

flexfilter:
	$(MAKE) -C $(FLEX_DIR)/tools/flexfilter
