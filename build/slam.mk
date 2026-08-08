.PHONY: slam slam-tools slam-requant pileup-snp trim-qc-tools

slam: slam-tools

slam-tools: slam-requant pileup-snp trim-qc-tools

slam-requant: core
	$(MAKE) -C $(SLAM_DIR)/tools/slam_requant

pileup-snp: core
	$(MAKE) -C $(SLAM_DIR)/tools/pileup_snp

trim-qc-tools: core
	$(MAKE) -C $(LEGACY_SRC_DIR) trim_qc_fastq trim_qc_merge
