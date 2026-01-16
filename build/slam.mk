.PHONY: slam slam-tools slam-requant pileup-snp

slam: core slam-tools

slam-tools: slam-requant pileup-snp

slam-requant:
	$(MAKE) -C $(SLAM_DIR)/tools/slam_requant

pileup-snp:
	$(MAKE) -C $(SLAM_DIR)/tools/pileup_snp
