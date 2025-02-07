

library(rtracklayer)
library(dplyr)
library(here)
library(ggplot2)

TAZ_peaks <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/TAZ_peak/TAZ_peaks.narrowPeak"))
YAP_peaks <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/YAP_peak/YAP_peaks.narrowPeak"))
TEAD4_peaks <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/TEAD4_peak/TEAD4_peaks.narrowPeak"))

YAP_overlap_TAZ_peaks <- subsetByOverlaps(YAP_peaks, TAZ_peaks)

YAP_overlap_TAZ_peaks_overlap_TEAD4 <- subsetByOverlaps(YAP_overlap_TAZ_peaks, TEAD4_peaks)
YAP_overlap_TAZ_peaks_overlap_TEAD4


# Reading the histone modification peaks
H3K4me1 <- import(here("/Users/aayushojha/Documents/chipseq/data/public_data/H3K4me1.bed"))
H3K4me3 <- import(here("/Users/aayushojha/Documents/chipseq/data/public_data/H3K4me3.bed"))
H3K27ac <- import(here("/Users/aayushojha/Documents/chipseq/data/public_data/H3K27ac.bed"))
H3K27ac

active_enhancers <- subsetByOverlaps(H3K4me1, H3K27ac)
inactive_enhancers <- subsetByOverlaps(H3K4me1, H3K27ac, invert = TRUE)
promoters <- subsetByOverlaps(H3K4me3, H3K4me1, invert = TRUE)

active_enhancers  
inactive_enhancers
promoters


# YAP/TAZ/TEAD peaks annotation YAP/TAZ/TEAD peaks were annotated as promoters or enhancers if their summit
# was overlapping with promoter or enhancer regions as defined above.
# Peaks with the summit falling in regions with no H3K4me1 or H3K4me3 peaks, or in NA regions were defined
# as “not assigned” and discarded from subsequent analyses.  

n_active_enhancers <- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4, active_enhancers) %>%
  length()

n_inactive_enhancers <- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4, inactive_enhancers) %>%
  length()

n_promoters <- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4, promoters) %>%
  length()

n_unclassified <- length(YAP_overlap_TAZ_peaks_overlap_TEAD4) - n_active_enhancers - n_inactive_enhancers - n_promoters

n_active_enhancers
n_inactive_enhancers
n_promoters
n_unclassified

# Creating a dataframe that has the values of counts of each category.
annotation_df <- data.frame(category = c("active_enhancers", "inactive_enhancers", "promoters", "unclassified"),
                            peak_number = c(n_active_enhancers, n_inactive_enhancers, n_promoters, n_unclassified)
                            )
annotation_df

# To change the order of categories in pie-chart so that it looks good.  
annotation_df$category <- factor(annotation_df$category, 
                                 levels = c("promoters", "active_enhancers", "inactive_enhancers", "unclassified")
                                 )

# Adding "percentage" and "label" columns to indicate the percentage share of each catrgory in pie-chart. 
annotation_df <- annotation_df %>%
  dplyr::mutate(
    percentage = peak_number/sum(peak_number)*100,
    label = paste0(round(percentage, 1), "%")
  )
annotation_df

# lets plot a pie chart 
library(ggplot2)
library(ggrepel)

colors<- c("#8D1E0F", "#F57D2B", "#FADAC4", "#D4DADA")
ggplot(annotation_df, aes(x = "", y = peak_number, fill = category)) + 
  geom_bar(stat = "identity", width=1, color="white") + 
  coord_polar("y", start = 0) +
  theme_void() + # Remove unnecessary axes
  labs(title = "                   Categorized YAP/TAZ/TEAD4 peaks") +
  scale_fill_manual(values = colors) + 
  geom_text(aes(label = label), position = position_stack(vjust = 0.5))






  