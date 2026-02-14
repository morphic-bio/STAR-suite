.PHONY: slam slam-tools slam-requant pileup-snp

slam: slam-tools

slam-tools: slam-requant pileup-snp

slam-requant: core
	$(MAKE) -C $(SLAM_DIR)/tools/slam_requant

pileup-snp: core
	$(MAKE) -C $(SLAM_DIR)/tools/pileup_snp
