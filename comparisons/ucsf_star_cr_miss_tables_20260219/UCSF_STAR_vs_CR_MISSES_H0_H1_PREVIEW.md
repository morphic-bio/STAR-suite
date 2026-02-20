# UCSF STAR vs CR Missing Rows (H0 and H1)

## H0 (Exact-match misses; STAR m=0, CR raw=0)

| feature | barcode_tru | star_exact_count | in_cr_call_list |
|---|---|---:|---:|
| negative_control_5_P1P2_A | CGGGTCAAGGCATCTT | 200 | 0 |
| LMO2_P1P2_B | ATCCATTTCTGAGCAT | 189 | 0 |
| VDR_P1P2_A | AGGATAAAGGAGAGTA | 183 | 0 |
| TLR4_P1P2_B | CCACTTGGTCGTTGGC | 169 | 0 |
| negative_control_5_P1P2_A | CAACAACCAACAACAA | 162 | 0 |
| NEUROG1_P1P2_A | CCAAGCGAGAAACTCA | 138 | 0 |
| RUNX1_P1P2_B | GCATCGGAGATGCTGG | 113 | 0 |
| DLX5_P1P2_A | GATGGAGCAGACATCT | 102 | 0 |
| ESRRG_P1_A | GTATTTCGTGATACAA | 58 | 0 |
| TBX20_P1P2_A | AAGGAATTCATGTCTT | 52 | 0 |
| HOXA9_P1P2_A | CTCAGTCCATTGAGGG | 44 | 0 |
| PITX2_P2_B | TGGTGATCAACTGATC | 42 | 0 |
| DNMT3B_P1_A | TTTGGTTAGACTGGGT | 37 | 0 |
| ESRRG_P1_A | TCGATTTCATGCAGGA | 34 | 0 |
| MYB_P1P2_B | TAACGACTCTCTCCGA | 28 | 0 |
| DNMT3B_P1_A | GACACGCCACCTAAAC | 27 | 0 |
| EZH1_P1P2_A | CCACAAAAGGGACCAT | 21 | 0 |
| LMX1A_P1P2_B | AGTTCGAGTGTGGACA | 15 | 0 |
| negative_control_5_P1P2_A | CTTGAGATCTGTCCGT | 15 | 0 |
| POU3F2_P1P2_A | TTTACCATCCGATAGT | 13 | 0 |
| TAL1_P2_A | AGACAAAAGTTGGAGC | 10 | 0 |
| HESX1_P1_A | CTCAACCCACGTACTA | 9 | 0 |
| RUNX2_P2_A | ACCTGAACAGCCGTTG | 9 | 0 |
| DNMT3A_P1P2_A | TGTTGAGAGTCTTCCC | 8 | 0 |
| SOX1_P1P2_A | ACGTAACGTCTACTGA | 7 | 0 |

## H1 (1-Hamming rescue misses; STAR delta, CR raw=0)

| feature | barcode_tru | star_m1_minus_m0 | classification |
|---|---|---:|---|
| BCL6_P1P2_A | ATCCACCTCTACTTCA | 136 | shift_to_partner_feature |
| BCL6_P1P2_A | TACACCCAGTTCACTG | 127 | shift_to_partner_feature |
| ZNF239_P1P2_A | ACTATCTTCAACCTTT | 117 | no_partner_signal_in_cr |
| BCL6_P1P2_A | CGAAGGATCTCCTGAC | 108 | shift_to_partner_feature |
| FOXD3_P1_B | AGGTAGGAGAAGATCT | 101 | no_partner_signal_in_cr |
| negative_control_5_P1P2_A | CGGGTCAAGGCATCTT | 87 | no_partner_signal_in_cr |
| negative_control_5_P1P2_A | CAACAACCAACAACAA | 73 | no_partner_signal_in_cr |
| BCL6_P1P2_A | CCTACGTTCCGCTTAC | 70 | shift_to_partner_feature |
| MDM2_P1P2_A | TTTCACATCGCTTGAA | 17 | no_partner_signal_in_cr |
| DNMT3A_P1P2_A | TGTTGAGAGTCTTCCC | 12 | no_partner_signal_in_cr |
| SNAI2_P1P2_B | TCCATGCTCTTCCGTG | 12 | no_partner_signal_in_cr |
| TLR4_P1P2_B | CCACTTGGTCGTTGGC | 10 | shift_to_partner_feature |
| negative_control_5_P1P2_A | CTTGAGATCTGTCCGT | 9 | no_partner_signal_in_cr |
| HOXA9_P1P2_A | CTCAGTCCATTGAGGG | 8 | shift_to_partner_feature |
| ESRRG_P1_A | GTATTTCGTGATACAA | 4 | shift_to_partner_feature |
| LMO2_P1P2_B | ATCCATTTCTGAGCAT | 4 | no_partner_signal_in_cr |
| ESRRG_P1_A | TCGATTTCATGCAGGA | 3 | shift_to_partner_feature |
| MYB_P1P2_B | TAACGACTCTCTCCGA | 3 | shift_to_partner_feature |
| RUNX1_P1P2_B | GCATCGGAGATGCTGG | 3 | shift_to_partner_feature |
| VDR_P1P2_A | AGGATAAAGGAGAGTA | 3 | no_partner_signal_in_cr |
| DNMT3B_P1_A | GACACGCCACCTAAAC | 2 | shift_to_partner_feature |
| NEUROG1_P1P2_A | CCAAGCGAGAAACTCA | 2 | shift_to_partner_feature |
| PITX2_P2_B | TGGTGATCAACTGATC | 2 | shift_to_partner_feature |
| ATF2_P1P2_A | GGAGGATGTTCCGTTC | 1 | no_partner_signal_in_cr |
| BCL6_P1P2_A | CCTCCAATCGCGTGCA | 1 | no_partner_signal_in_cr |
