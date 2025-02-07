
# We are trying to plot figures from a research paper : 

# The experiment runs ChiP-seq experiment. 

# So far, we have taken the sequencing data from ENA, and we completed our pre-processing. 
# Ultimately, we have bigwig files, peak files and bed files, which are essential formats of ChIP-Seq data analysis. 
# We will now use these data for visualization.  

library(GenomicRanges)
library(rtracklayer)
library(here)
library(dplyr)
library(ggplot2)

TAZ_peaks <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/TAZ_peak/TAZ_peaks.narrowPeak"))
YAP_peaks <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/YAP_peak/YAP_peaks.narrowPeak"))
TEAD4_peaks <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/TEAD4_peak/TEAD4_peaks.narrowPeak"))

YAP_overlap_TAZ_peaks <- subsetByOverlaps(YAP_peaks, TAZ_peaks)
YAP_overlap_TAZ_peaks
YAP_overlap_TAZ_peaks_overlap_TEAD4 <- subsetByOverlaps(YAP_overlap_TAZ_peaks, TEAD4_peaks)
YAP_overlap_TAZ_peaks_overlap_TEAD4

# use rtracklayer to write the GenomicRanges object to file
export(YAP_overlap_TAZ_peaks_overlap_TEAD4,
       con = here("/Users/aayushojha/Documents/chipseq/data/fastq/YAP_TAZ_TEAD4_common.bed"))


library(readr)
counts <- read_tsv(here("/Users/aayushojha/Documents/chipseq/data/fastq/YAP_TAZ_TEAD4_counts.tsv"), col_names = FALSE)
colnames(counts)<- c("chr", "start", "end", "name", "score", "value", "YAP", "TAZ", "TEAD4")
head(counts)

# normalize the counts to CPM (counts per million)
counts <- counts %>%
  mutate(YAP = YAP/23653961 * 10^6,
         TAZ = TAZ/26789648 * 10^6,
         TEAD4 = TEAD4/34332907 * 10^6)

head(counts)


# plot the values TEAD4 vs YAP
ggplot(counts, aes(x=TEAD4, y=YAP)) +
  geom_point()

# Wanting to see where the outlier is, in which chromosome and it bp position
# so that we can know if it is a blacklisted region 
# Know more about blacklisted region from here : 
# https://github.com/Boyle-Lab/Blacklist/blob/master/README.md
counts %>%
  filter(TEAD4 > 40)

#  Blacklisted regions can be removed from our data. But, it is not a blacklisted region.
# Thus, it needs to be included in our figure. 
# To include it, we will scale x and y axis in log2. 
# Additionally, we will also include the r^2 value using ggpmisc package. 


library(ggpmisc)
ggplot(counts, aes(x=TEAD4, y=YAP)) +
  geom_point(color = 'red') +
  geom_smooth(method = "lm", se = FALSE, color = 'black') +
  stat_poly_eq(
    aes(label = ..rr.label..),
    formula = y ~ x,
    parse = TRUE,
    color = 'black'
  ) +
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2') +
  theme_classic(base_size = 14) +
  xlab('TEAD4 signal') +
  ylab('YAP signal')

correlation_coefficent<- cor(log2(counts$TEAD4), log2(counts$YAP))
correlation_coefficent

R_squared<- correlation_coefficent^2
R_squared


# Tead4 vs TAZ signals :
ggplot(counts, aes(x=TEAD4, y=TAZ)) +
  geom_point(color = 'red') +
  geom_smooth(method = "lm", se = FALSE, color = 'black') +
  stat_poly_eq(
    aes(label = ..rr.label..),
    formula = y ~ x,
    parse = TRUE,
    color = 'black'
  ) +
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2') +
  theme_classic(base_size = 14) +
  xlab('TEAD4 signal') +
  ylab('TAZ signal')



# Figure 1f: Absolute distance of YAP peaks (n=7709), TAZ peaks (n=9798), TEAD4 peaks (n=8406) 
# or overlapping YAP/TAZ/TEAD peaks (n=5522) to the nearest TSS.

# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(GenomicFeatures)

# Get the transcripts 
hg38_transcripts <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
head(hg38_transcripts)

# Get the Transcription start site
tss_gr <- promoters(hg38_transcripts, upstream=0, downstream = 1)
head(tss_gr)

distance_to_tss <- distanceToNearest(YAP_peaks, tss_gr)
distance_to_tss

mcols(distance_to_tss)
head(mcols(distance_to_tss)$distance)

YAP_dist <- mcols(distanceToNearest(YAP_peaks, tss_gr))$distance
TAZ_dist <- mcols(distanceToNearest(TAZ_peaks, tss_gr))$distance
TEAD4_dist <- mcols(distanceToNearest(TEAD4_peaks, tss_gr))$distance
YAP_TAZ_TEAD4_dist <- mcols(distanceToNearest(YAP_overlap_TAZ_peaks_overlap_TEAD4, tss_gr))$distance

tss_distance_df <- bind_rows(data.frame(factor ="YAP", distance = YAP_dist),
                             data.frame(factor ="TAZ", distance = TAZ_dist),
                             data.frame(factor ="TEAD4", distance = TEAD4_dist),
                             data.frame(factor ="YAP/TAZ/TEAD4", distance = YAP_TAZ_TEAD4_dist))

head(tss_distance_df)


counts_per_category<- tss_distance_df %>%
  mutate(category = case_when(
    distance < 1000 ~ '<1kb',
    distance >=1000 & distance < 10000 ~ '1-10kb',
    distance >=10000 & distance < 100000 ~ '10-100kb',
    distance > 100000 ~ '>100kb',
  )) %>%
  group_by(factor, category) %>%
  count()
counts_per_category


total_counts<- tss_distance_df %>%
  mutate(category = case_when(
    distance < 1000 ~ '<1kb',
    distance >=1000 & distance < 10000 ~ '1-10kb',
    distance >=10000 & distance < 100000 ~ '10-100kb',
    distance > 100000 ~ '>100kb',
  )) %>%
  count(factor, name = 'total')
total_counts


merged_df<- left_join(counts_per_category, total_counts)
merged_df %>%
  mutate(Percentage = n/total * 100) %>%
  ggplot(aes(x= factor, y = Percentage, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Distance to TSS",
    x = "Group",
    y = "Percentage"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_classic(base_size = 14)


merged_df$category<- factor(merged_df$category, 
                            levels = c("<1kb", "1-10kb", "10-100kb", ">100kb"))
merged_df %>%
  mutate(Percentage = n/total * 100) %>%
  ggplot(aes(x= factor, y = Percentage, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "      Distance to TSS",
    x = "Peaks",
    y = "Percentage"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = c("#EF3E2B", "#F16161", "#F59595", "#FAD1C8")) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Tilts x-axis labels












