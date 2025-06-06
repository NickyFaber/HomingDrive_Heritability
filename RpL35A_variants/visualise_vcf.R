rm(list = ls())
setwd("/Users/s2018147/Library/CloudStorage/OneDrive-WageningenUniversity&Research/PhD/03_Genetic_diversity/DrosEU pop gen")

library(tidyverse)
library(vcfR)
library(patchwork)
library(pheatmap)
library(viridisLite)
library(ggplotify)

# Figure aesthetics
PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + 
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        plot.title=element_text(size=12, hjust=0.5), 
        legend.title=element_text(size=12),
        legend.position = "bottom", 
        legend.justification = "center",
        axis.title=element_text(size=12))

# Filtering the DrosEU vcf file from Kapun et al. (2020) for Rpl35A region only on command line
# cat DrosEU-PoolSnp_full-filtered_ann.vcf | grep '^3R\t' > 3R_only.vcf
# cat 3R_only.vcf | awk '{ if ($2 >= 5465630 && $2 <= 5467101) print }' > Rpl35A.vcf

# Check annotation positions
# cat GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff | grep 'RpL35A'
# cat GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff | grep 'RpL35A' | awk '{print $3,$4,$5}'

rpl35a_data <- read.vcfR("Rpl35A.vcf")

all_info <- data.frame(rpl35a_data@gt)
vcf_table <- cbind(rpl35a_data@fix, all_info)

vcf_table_long <- pivot_longer(vcf_table, X1_Mauternbach:X51_Valday, names_to = "population", values_to = "info")
vcf_table_long <- separate(vcf_table_long, info, sep = ":",
         into = c("GT","RD","AD","DP","FREQ"))
vcf_table_long$FREQ <- as.numeric(vcf_table_long$FREQ) / 100
vcf_table_long$POS <- as.numeric(vcf_table_long$POS)
vcf_table_long$population <- factor(vcf_table_long$population, 
                                    levels = unique(vcf_table_long$population))

vcf_table_long <- vcf_table_long %>%
  mutate(tested = ifelse(population %in% c("X1_Mauternbach", "X2_Mauternbach", "X31_Munich", "X32_Munich", "X33_Recarei", "X36_Akaa","X37_Akaa", "X39_Karensminde", "X41_Karensminde", "X50_Uman"), 
                         'red', 'grey30'))

whole_region <- data.frame(POS = 5465630:5467101, gene = NA, mRNA = NA, exon = NA, CDS = NA, gRNA = NA)
#whole_region$gene[whole_region$POS >= 5465630 & whole_region$POS <= 5467101] = TRUE
#whole_region$mRNA[whole_region$POS >= 5465630 & whole_region$POS <= 5467101] = TRUE
whole_region$exon[(whole_region$POS >= 5467069 & whole_region$POS <= 5467101) |
                  (whole_region$POS >= 5466980 & whole_region$POS <= 5467008) |
                  (whole_region$POS >= 5466153 & whole_region$POS <= 5466381) |
                  (whole_region$POS >= 5465800 & whole_region$POS <= 5466027) |
                  (whole_region$POS >= 5465630 & whole_region$POS <= 5465731)] = TRUE
whole_region$CDS[(whole_region$POS >= 5466153 & whole_region$POS <= 5466374) |
                 (whole_region$POS >= 5465800 & whole_region$POS <= 5466027) |
                 (whole_region$POS >= 5465708 & whole_region$POS <= 5465731)] = TRUE
whole_region$gRNA[(whole_region$POS >= 5465987 & whole_region$POS <= 5466006) |
                   (whole_region$POS >= 5465961 & whole_region$POS <= 5465980)] = TRUE
whole_region$primer_F[(whole_region$POS >= 5466216 & whole_region$POS <= 5466233)] = TRUE
whole_region$primer_R[(whole_region$POS >= 5465734 & whole_region$POS <= 5465754)] = TRUE
whole_region_long <- pivot_longer(whole_region, cols = gene : primer_R, 
                                  names_to = "component", 
                                  values_to = "present") %>%
  drop_na(present)

whole_region_long$component <- factor(whole_region_long$component, 
                                      levels = c("gene", "mRNA", "exon", "CDS", "gRNA", "primer_F", "primer_R"),
                                      labels = c("gene", "mRNA", "exon", "CDS", "gRNA target", "primer", "primer"))

p1 <- ggplot() +
  geom_hline(yintercept = 3, linewidth = 3) +
  geom_tile(data = whole_region_long[whole_region_long$component=="exon" | whole_region_long$component=="CDS",], aes(x = POS, y = 1, fill = component)) +
  geom_tile(data = whole_region_long[whole_region_long$component=="exon" | whole_region_long$component=="CDS",], aes(x = POS, y = 2, fill = component)) +
  geom_tile(data = whole_region_long[whole_region_long$component=="exon" | whole_region_long$component=="CDS",], aes(x = POS, y = 3, fill = component)) +
  geom_tile(data = whole_region_long[whole_region_long$component=="exon" | whole_region_long$component=="CDS",], aes(x = POS, y = 4, fill = component)) +
  geom_tile(data = whole_region_long[whole_region_long$component=="exon" | whole_region_long$component=="CDS",], aes(x = POS, y = 5, fill = component)) +
  geom_tile(data = whole_region_long[whole_region_long$component=="gRNA target" | whole_region_long$component=="primer",], aes(x = POS, y = 7, fill = component)) +
  scale_fill_manual(values = c("lightblue", "midnightblue", "gold", "hotpink", "hotpink", "hotpink", "hotpink"), name = "Gene components") +
  xlim(c(5467102, 5465629)) +
  xlab("Position on chromosome 3R (bp)") +
  ylab("RpL35A") +
  PaperTheme + theme(legend.position = "top", axis.text.y=element_blank(), axis.ticks.y=element_blank())
p1

p2 <- ggplot() +
  geom_tile(data = vcf_table_long, aes(x = POS, y = 1)) +
  xlim(c(5467102, 5465629)) +
  xlab(element_blank()) +
  ylab("Variant sites") +
  PaperTheme + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
p2

frequencies <- dplyr::select(vcf_table_long, population, POS, FREQ) %>%
  pivot_wider(names_from = POS, values_from = FREQ)
frequencies_matrix <- as.matrix(frequencies[,-1])
row.names(frequencies_matrix) <- frequencies$population

gRNA1 <- as.matrix(rep(0, 48))
colnames(gRNA1) <- "5466006 - 5465987 (gRNA 1)"
gRNA2 <- as.matrix(rep(0, 48))
colnames(gRNA2) <- "5465980 - 5465961 (gRNA 2)"

frequencies_matrix <- cbind(frequencies_matrix[,23:10], gRNA1, gRNA2, frequencies_matrix[,9:1])

p3 <- as.ggplot(pheatmap(mat = frequencies_matrix, 
                         cluster_cols = FALSE,
                         cluster_rows = FALSE,
                         fontsize = 8,
                         color = colorRampPalette(colors=c("white","skyblue"))(100),
                         angle_col = 45,
                         gaps_col = c(14,15,16)))

p <- p1 / p2 / plot_spacer() / p3 +
  plot_layout(heights = c(1, 1, 0, 12)) +
  plot_annotation(tag_levels = 'A'); p

ggsave(plot = p, filename = "RpL35A_variants.pdf", height = 25, width = 20.7, unit = "cm")
ggsave(plot = p, filename = "RpL35A_variants.png", height = 25, width = 20.7, unit = "cm", dpi = 600)





