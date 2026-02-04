# # Will count all as unstranded for fairness of comparison (-s 0)
# 
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_GSE93202_VAT_Macrophages_aging_counts.txt \
#       GSM2447101_VAT_macrophages_3m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2447103_VAT_macrophages_3m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2447105_VAT_macrophages_3m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2447107_VAT_macrophages_24m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2447109_VAT_macrophages_24m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2447111_VAT_macrophages_24m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2447113_VAT_macrophages_24m_male_STAR_Aligned.sortedByCoord.out.bam 
#       
#       
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_GSE93202_Spleen_Macrophages_aging_counts.txt \
#       GSM2447115_Spleen_macrophages_3m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2447116_Spleen_macrophages_3m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2447119_Spleen_macrophages_24m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2447120_Spleen_macrophages_24m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2447122_Spleen_macrophages_24m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2447124_Spleen_macrophages_24m_male_STAR_Aligned.sortedByCoord.out.bam 
#       
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_GSE98249_BAM_macrophages_aging_counts.txt \
#       GSM2589849_BAM_macrophages_3m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2589851_BAM_macrophages_3m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2589853_BAM_macrophages_3m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2589843_BAM_macrophages_24m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2589845_BAM_macrophages_24m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2589847_BAM_macrophages_24m_male_STAR_Aligned.sortedByCoord.out.bam 
#       
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_GSE98249_BMM_macrophages_aging_counts.txt \
#       GSM2589850_BMM_macrophages_3m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2589852_BMM_macrophages_3m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2589854_BMM_macrophages_3m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2589844_BMM_macrophages_24m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2589846_BMM_macrophages_24m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2589848_BMM_macrophages_24m_male_STAR_Aligned.sortedByCoord.out.bam 
#       
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_GSE98401_Brain_microglia_aging_counts.txt \
#       GSM2593467_Microglia_brain_3m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2593468_Microglia_brain_3m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2593469_Microglia_brain_3m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2593464_Microglia_brain_22m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2593465_Microglia_brain_22m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2593466_Microglia_brain_22m_male_STAR_Aligned.sortedByCoord.out.bam 
#       
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_GSE132882_Nerve_Macrophage_aging_counts.txt \
#       GSM2857711_Nerve_Macrophage_3m_female_Rep_1_STAR_Aligned.sortedByCoord.out.bam \
#       GSM2857712_Nerve_Macrophage_3m_female_Rep_2_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3895604_Nerve_Macrophage_15m_female_Rep_1_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3895605_Nerve_Macrophage_15m_female_Rep_2_STAR_Aligned.sortedByCoord.out.bam
#       
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_GSE142580_SkM_Macrophages_aging_counts.txt \
#       GSM4232357_SkM_Macrophages_4m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4232358_SkM_Macrophages_4m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4232359_SkM_Macrophages_4m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4232360_SkM_Macrophages_4m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4232361_SkM_Macrophages_4m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4232352_SkM_Macrophages_18m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4232353_SkM_Macrophages_18m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4232354_SkM_Macrophages_18m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4232355_SkM_Macrophages_18m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4232356_SkM_Macrophages_18m_NA_STAR_Aligned.sortedByCoord.out.bam 
#       
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_GSE154832_eWAT_Phagocytic_SVF_aging_counts.txt \
#       GSM4681350_eWAT_phagocytic_SVF_2m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4681354_eWAT_phagocytic_SVF_2m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4681352_eWAT_phagocytic_SVF_20m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4681356_eWAT_phagocytic_SVF_20m_male_STAR_Aligned.sortedByCoord.out.bam 
#       
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_PRJNA682234_Callus_Macrophages_aging_counts.txt \
#       SAMN16984882_Callus_Macrophages_3m_males_Y_Cont6_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984883_Callus_Macrophages_3m_males_Y_Cont7_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984884_Callus_Macrophages_3m_males_Y_Cont8_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984885_Callus_Macrophages_3m_males_Y_Cont9_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984886_Callus_Macrophages_3m_males_Y_Cont10_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984887_Callus_Macrophages_3m_males_Y_Cont11_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984888_Callus_Macrophages_3m_males_Y_Cont12_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984889_Callus_Macrophages_3m_males_Y_Cont13_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984890_Callus_Macrophages_3m_males_Y_Cont14_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984891_Callus_Macrophages_3m_males_Y_Cont15_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984892_Callus_Macrophages_3m_males_Y_Cont16_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984872_Callus_Macrophages_24m_males_O_Cont6_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984873_Callus_Macrophages_24m_males_O_Cont7_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984874_Callus_Macrophages_24m_males_O_Cont8_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984875_Callus_Macrophages_24m_males_O_Cont9_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984876_Callus_Macrophages_24m_males_O_Cont10_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984877_Callus_Macrophages_24m_males_O_Cont11_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984878_Callus_Macrophages_24m_males_O_Cont12_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984879_Callus_Macrophages_24m_males_O_Cont13_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984880_Callus_Macrophages_24m_males_O_Cont14_STAR_Aligned.sortedByCoord.out.bam \
#       SAMN16984881_Callus_Macrophages_24m_males_O_Cont15_STAR_Aligned.sortedByCoord.out.bam 
#       
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_PRJNA451159_BMDM_Ames_Dwarf_counts.txt \
#       SRX3995441_BMDM_6m_males_Ames_CTL_Ndf_STAR_Aligned.sortedByCoord.out.bam \
#       SRX3995448_BMDM_6m_males_Ames_CTL_Ndf_STAR_Aligned.sortedByCoord.out.bam \
#       SRX3995451_BMDM_6m_males_Ames_CTL_Ndf_STAR_Aligned.sortedByCoord.out.bam \
#       SRX3995437_BMDM_6m_males_Ames_Dwarf_dfdf_STAR_Aligned.sortedByCoord.out.bam \
#       SRX3995438_BMDM_6m_males_Ames_Dwarf_dfdf_STAR_Aligned.sortedByCoord.out.bam \
#       SRX3995439_BMDM_6m_males_Ames_Dwarf_dfdf_STAR_Aligned.sortedByCoord.out.bam 
#        
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_GSE134397_Alveolar_Macrophages_aging_counts.txt \
#       GSM3945726_Alveolar_Macrophages_4m_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945727_Alveolar_Macrophages_4m_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945728_Alveolar_Macrophages_4m_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945734_Alveolar_Macrophages_18m_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945735_Alveolar_Macrophages_18m_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945736_Alveolar_Macrophages_18m_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945737_Alveolar_Macrophages_18m_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945738_Alveolar_Macrophages_18m_Males_STAR_Aligned.sortedByCoord.out.bam 
#       
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_GSE134397_Alveolar_Macrophages_CTL_aging_counts.txt \
#       GSM3945865_Alveolar_Macrophages_AM_4m_control_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945867_Alveolar_Macrophages_AM_4m_control_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945869_Alveolar_Macrophages_AM_4m_control_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945871_Alveolar_Macrophages_AM_4m_control_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945883_Alveolar_Macrophages_AM_18m_control_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945885_Alveolar_Macrophages_AM_18m_control_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945887_Alveolar_Macrophages_AM_18m_control_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945889_Alveolar_Macrophages_AM_18m_control_Males_STAR_Aligned.sortedByCoord.out.bam 
#       
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_GSE134397_Alveolar_Macrophages_Metformin_and_aging_counts.txt \
#       GSM3945865_Alveolar_Macrophages_AM_4m_control_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945867_Alveolar_Macrophages_AM_4m_control_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945869_Alveolar_Macrophages_AM_4m_control_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945871_Alveolar_Macrophages_AM_4m_control_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945883_Alveolar_Macrophages_AM_18m_control_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945885_Alveolar_Macrophages_AM_18m_control_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945887_Alveolar_Macrophages_AM_18m_control_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945889_Alveolar_Macrophages_AM_18m_control_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945873_Alveolar_Macrophages_AM_4m_metformin_60d_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945875_Alveolar_Macrophages_AM_4m_metformin_60d_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945877_Alveolar_Macrophages_AM_4m_metformin_60d_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945879_Alveolar_Macrophages_AM_4m_metformin_60d_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945881_Alveolar_Macrophages_AM_4m_metformin_60d_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945891_Alveolar_Macrophages_AM_18m_metformin_60d_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945893_Alveolar_Macrophages_AM_18m_metformin_60d_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3945895_Alveolar_Macrophages_AM_18m_metformin_60d_Males_STAR_Aligned.sortedByCoord.out.bam 
#        
#        
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_GSE124829_ImmGen_Peritoneal_Macrophages_aging_counts.txt \
#       GSM3555595_ImmGen_Macrophages_Peritoneal_Cavity_female_2_months_old_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3555596_ImmGen_Macrophages_Peritoneal_Cavity_male_2_months_old_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3555597_ImmGen_Macrophages_Peritoneal_Cavity_female_6_months_old_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3555598_ImmGen_Macrophages_Peritoneal_Cavity_male_6_months_old_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3555593_ImmGen_Macrophages_Peritoneal_Cavity_female_17_months_old_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3555594_ImmGen_Macrophages_Peritoneal_Cavity_male_20_months_old_STAR_Aligned.sortedByCoord.out.bam
#       
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_GSE124872_Alveolar_Macrophages_aging_counts.txt \
#       GSM3557702_Alveolar_Macrophages_3m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3557705_Alveolar_Macrophages_3m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3557707_Alveolar_Macrophages_3m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3557710_Alveolar_Macrophages_3m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3557713_Alveolar_Macrophages_3m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3557703_Alveolar_Macrophages_22m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3557704_Alveolar_Macrophages_22m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3557706_Alveolar_Macrophages_22m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3557708_Alveolar_Macrophages_22m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3557709_Alveolar_Macrophages_22m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3557711_Alveolar_Macrophages_22m_NA_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3557715_Alveolar_Macrophages_22m_NA_STAR_Aligned.sortedByCoord.out.bam 
#       
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_GSE131869_Microglia_aging_counts.txt \
#       GSM3823683_Microglia_2m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823684_Microglia_2m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823685_Microglia_2m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823686_Microglia_2m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823687_Microglia_2m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823688_Microglia_2m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823689_Microglia_2m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823690_Microglia_2m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823691_Microglia_2m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823692_Microglia_2m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823693_Microglia_2m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823721_Microglia_12m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823722_Microglia_12m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823723_Microglia_12m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823714_Microglia_12m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823715_Microglia_12m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823716_Microglia_12m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823717_Microglia_12m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823718_Microglia_12m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823719_Microglia_12m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823720_Microglia_12m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823750_Microglia_24m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823755_Microglia_24m_male_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823756_Microglia_24m_male_STAR_Aligned.sortedByCoord.out.bam  \
#       GSM3823747_Microglia_24m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823748_Microglia_24m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823749_Microglia_24m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823751_Microglia_24m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823752_Microglia_24m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823753_Microglia_24m_female_STAR_Aligned.sortedByCoord.out.bam \
#       GSM3823754_Microglia_24m_female_STAR_Aligned.sortedByCoord.out.bam 
#       
#       
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_GSE137028_Microglia_aging_counts.txt \
#       GSM4065891_Microglia_WT_2m1_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4065892_Microglia_WT_2m2_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4065893_Microglia_WT_2m3_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4065894_Microglia_WT_4m1_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4065895_Microglia_WT_4m2_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4065896_Microglia_WT_4m3_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4065897_Microglia_WT_6m1_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4065898_Microglia_WT_6m2_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4065899_Microglia_WT_6m3_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4065900_Microglia_WT_9m1_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4065901_Microglia_WT_9m2_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4065902_Microglia_WT_9m3_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4065903_Microglia_WT_12m1_Males_STAR_Aligned.sortedByCoord.out.bam \
#       GSM4065904_Microglia_WT_12m2_Males_STAR_Aligned.sortedByCoord.out.bam
#       
#       
# featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2021-06-04_GSE145295_Alveolar_Macrophages_aging_counts.txt \
#        GSM4314353_lung_macrophages_2.5m_STAR_Aligned.sortedByCoord.out.bam \
#        GSM4314354_lung_macrophages_2.5m_STAR_Aligned.sortedByCoord.out.bam \
#        GSM4314355_lung_macrophages_2.5m_STAR_Aligned.sortedByCoord.out.bam \
#        GSM4314356_lung_macrophages_2.5m_STAR_Aligned.sortedByCoord.out.bam \
#        GSM4314357_lung_macrophages_18m_STAR_Aligned.sortedByCoord.out.bam \
#        GSM4314358_lung_macrophages_18m_STAR_Aligned.sortedByCoord.out.bam \
#        GSM4314359_lung_macrophages_18m_STAR_Aligned.sortedByCoord.out.bam \
#        GSM4314360_lung_macrophages_18m_STAR_Aligned.sortedByCoord.out.bam 
       
##########################################################################################################################################################################
       

featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2025-03-07_GSE134397_Alveolar_Macrophages_aging_Mixed_Sex_counts.txt \
      GSM3945726_04_months_MixedSex_Alveolar_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3945727_04_months_MixedSex_Alveolar_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3945728_04_months_MixedSex_Alveolar_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3945729_04_months_MixedSex_Alveolar_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3945734_18_months_MixedSex_Alveolar_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3945735_18_months_MixedSex_Alveolar_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3945736_18_months_MixedSex_Alveolar_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3945737_18_months_MixedSex_Alveolar_Macrophages_STAR_Aligned.sortedByCoord.out.bam


featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2025-03-07_GSE190689_Alveolar_Macrophages_aging_UnkSex_counts.txt \
      GSM5729052_06months_NA_Alveolar_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM5729061_06months_NA_Alveolar_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM5729064_06months_NA_Alveolar_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM5729066_06months_NA_Alveolar_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM5729053_18months_NA_Alveolar_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM5729054_18months_NA_Alveolar_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM5729055_18months_NA_Alveolar_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM5729067_18months_NA_Alveolar_Macrophages_STAR_Aligned.sortedByCoord.out.bam


featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2025-03-07_GSE199763_SkinWound_Macrophages_aging_Female_counts.txt \
      GSM5984047_03months_Female_Skin_Wound_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM5984046_03months_Female_Skin_Wound_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM5984045_03months_Female_Skin_Wound_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM5984050_24months_Female_Skin_Wound_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM5984049_24months_Female_Skin_Wound_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM5984048_24months_Female_Skin_Wound_Macrophages_STAR_Aligned.sortedByCoord.out.bam


featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2025-03-07_GSE199879_Spleen_Red_Pulp_Macrophages_aging_Female_counts.txt \
      GSM5989996_02months_Female_Spleen_Red_Pulp_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM5989997_02months_Female_Spleen_Red_Pulp_Macrophages_STAR_Aligned.out.bam \
      GSM5989998_02months_Female_Spleen_Red_Pulp_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM5989999_02months_Female_Spleen_Red_Pulp_Macrophages_STAR_Aligned.out.bam \
      GSM5990003_10months_Female_Spleen_Red_Pulp_Macrophages_STAR_Aligned.out.bam \
      GSM5990002_10months_Female_Spleen_Red_Pulp_Macrophages_STAR_Aligned.out.bam \
      GSM5990001_10months_Female_Spleen_Red_Pulp_Macrophages_STAR_Aligned.out.bam \
      GSM5990000_10months_Female_Spleen_Red_Pulp_Macrophages_STAR_Aligned.out.bam


featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2025-03-07_GSE205395_Skeletal_Muscle_Macrophages_aging_Male_counts.txt \
      GSM6211369_03months_Males_Skeletal_Muscle_Macrophages_STAR_Aligned.out.bam \
      GSM6211370_03months_Males_Skeletal_Muscle_Macrophages_STAR_Aligned.out.bam \
      GSM6211371_26months_Males_Skeletal_Muscle_Macrophages_STAR_Aligned.out.bam \
      GSM6211372_26months_Males_Skeletal_Muscle_Macrophages_STAR_Aligned.out.bam

featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm10_genes.gtf -o 2025-03-07_GSE128830_Peritoneal_Macrophages_aging_counts.txt \
      GSM3686814_03months_ZT0_1_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686815_03months_ZT0_2_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686816_03months_ZT0_3_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686817_03months_ZT4_1_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686818_03months_ZT4_2_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686819_03months_ZT4_3_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686820_03months_ZT8_1_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686821_03months_ZT8_2_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686822_03months_ZT8_3_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686823_03months_ZT12_1_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686824_03months_ZT12_2_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686825_03months_ZT12_3_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686826_03months_ZT16_1_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686827_03months_ZT16_2_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686828_03months_ZT16_3_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686829_03months_ZT20_1_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686830_03months_ZT20_2_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686831_03months_ZT20_3_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686832_03months_ZT24_1_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686833_03months_ZT24_2_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686834_03months_ZT24_3_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686835_20months_ZT0_1_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686836_20months_ZT0_2_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686837_20months_ZT0_3_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686838_20months_ZT4_1_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686839_20months_ZT4_2_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686840_20months_ZT4_3_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686841_20months_ZT8_1_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686842_20months_ZT8_2_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686843_20months_ZT8_3_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686844_20months_ZT12_1_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686845_20months_ZT12_2_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686846_20months_ZT12_3_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686847_20months_ZT16_1_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686848_20months_ZT16_2_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686849_20months_ZT16_3_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686850_20months_ZT20_1_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686851_20months_ZT20_2_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686852_20months_ZT20_3_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686853_20months_ZT24_1_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686854_20months_ZT24_2_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam \
      GSM3686855_20months_ZT24_3_Males_Peritoneal_Macrophages_STAR_Aligned.sortedByCoord.out.bam

