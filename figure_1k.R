
library(rtracklayer)
library(dplyr)
library(here)
library(ggplot2)

TAZ_peaks <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/TAZ_peak/TAZ_peaks.narrowPeak"))
YAP_peaks <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/YAP_peak/YAP_peaks.narrowPeak"))
TEAD4_peaks <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/TEAD4_peak/TEAD4_peaks.narrowPeak"))

YAP_overlap_TAZ_peaks <- subsetByOverlaps(YAP_peaks, TAZ_peaks)
YAP_overlap_TAZ_peaks


library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

# Get the transcripts 
hg38_transcripts <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
head(hg38_transcripts)

# Get the Transcription start site
tss_gr <- promoters(hg38_transcripts, upstream=0, downstream = 1)
head(tss_gr)

distance_to_tss <- distanceToNearest(YAP_overlap_TAZ_peaks, tss_gr)
distance_to_tss

mcols(distance_to_tss)
head(mcols(distance_to_tss)$distance)

YAP_TAZ_dist <- mcols(distance_to_tss)$distance
YAP_TAZ_dist

tss_distance_df <- data.frame(factor = "YAP/TAZ Peaks", distance = YAP_TAZ_dist)

length(tss_distance_df)
nrow(tss_distance_df)  

counts_per_category <- tss_distance_df %>%
  mutate(category = case_when(
    distance < 1000 ~ "<1kb",
    distance >=1000 & distance < 10000 ~ "1-10kb",
    distance >=10000 & distance < 100000 ~ '10-100kb',
    distance > 100000 ~ '>100kb',
  )) %>%
  group_by(factor, category) %>%
  count()

counts_per_category



total_counts<- tss_distance_df %>%
  count(factor, name = 'total')

total_counts  


merged_df<- left_join(counts_per_category, total_counts)

merged_df %>%
  mutate(Percentage = n/total * 100,
         category = factor(category, levels = c("<1kb", "1-10kb", "10-100kb", ">100kb"))
  ) %>%
  ggplot(aes(x= category, y = Percentage, fill = category)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(
    title = "YAP/TAZ peaks associated\nwith regulated genes",
    x = "Distance to TSS (kb)",
    y = "% of peaks"
  ) +
  theme(
    plot.title = element_text(size = 55, angle = 45, hjust = 0.5)
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = c("#EF3E2B", "#F16161", "#F59595", "#FAD1C8")) +
  scale_x_discrete(
    limits = c("<1kb", "1-10kb", "10-100kb", ">100kb"),
    labels = c("<1", "1-10", "10-100", ">100")
  ) +
  theme_classic(base_size = 16, base_line_size = 0.5)










  