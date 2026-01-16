.PHONY: tools tools-clean vbem-tools yremove-tools
.PHONY: vbem-compute-expected-gc vbem-sample-fld vbem-compute-gc-bias vbem-em-quant
.PHONY: vbem-ec-filter-test vbem-tximport-compat vbem-trimvalidate

tools: flex-tools slam-tools vbem-tools yremove-tools

vbem-tools: vbem-compute-expected-gc vbem-sample-fld vbem-compute-gc-bias vbem-em-quant vbem-ec-filter-test vbem-tximport-compat vbem-trimvalidate

vbem-compute-expected-gc:
	$(MAKE) -C $(VBEM_DIR)/tools/compute_expected_gc

vbem-sample-fld:
	$(MAKE) -C $(VBEM_DIR)/tools/sample_fld

vbem-compute-gc-bias:
	$(MAKE) -C $(VBEM_DIR)/tools/compute_gc_bias

vbem-em-quant:
	$(MAKE) -C $(VBEM_DIR)/tools/em_quant

vbem-ec-filter-test:
	$(MAKE) -C $(VBEM_DIR)/tools/ec_filter_test

vbem-tximport-compat:
	$(MAKE) -C $(VBEM_DIR)/tools/tximport_compat

vbem-trimvalidate:
	$(MAKE) -C $(VBEM_DIR)/tools/trimvalidate

yremove-tools:
	$(MAKE) -C $(YREMOVE_FASTQ_DIR)/tools/remove_y_reads

tools-clean:
	$(MAKE) -C $(VBEM_DIR)/tools/compute_expected_gc clean
	$(MAKE) -C $(VBEM_DIR)/tools/sample_fld clean
	$(MAKE) -C $(VBEM_DIR)/tools/compute_gc_bias clean
	$(MAKE) -C $(VBEM_DIR)/tools/em_quant clean
	$(MAKE) -C $(VBEM_DIR)/tools/ec_filter_test clean
	$(MAKE) -C $(VBEM_DIR)/tools/tximport_compat clean
	$(MAKE) -C $(VBEM_DIR)/tools/trimvalidate clean
	$(MAKE) -C $(SLAM_DIR)/tools/slam_requant clean
	$(MAKE) -C $(SLAM_DIR)/tools/pileup_snp clean
	$(MAKE) -C $(FLEX_DIR)/tools/flexfilter clean
	$(MAKE) -C $(YREMOVE_FASTQ_DIR)/tools/remove_y_reads clean
