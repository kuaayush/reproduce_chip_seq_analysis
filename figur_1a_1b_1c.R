
# We are trying to plot figures from a research paper : 

# The experiment runs ChiP-seq experiment. 

# So far, we have taken the sequencing data from ENA, and we completed our pre-processing. 
# Ultimately, we have bigwig files, peak files and bed files, which are essential formats of ChIP-Seq data analysis. 
# We will now use these data for visualization.  

library(GenomicRanges)
library(rtracklayer)
library(here)

TAZ_peaks <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/TAZ_peak/TAZ_peaks.narrowPeak"))
YAP_peaks <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/YAP_peak/YAP_peaks.narrowPeak"))

TAZ_peaks
YAP_peaks

TAZ_overlap_YAP_peaks <-subsetByOverlaps(TAZ_peaks, YAP_peaks)
length(TAZ_overlap_YAP_peaks)

YAP_overlap_TAZ_peaks <-subsetByOverlaps(YAP_peaks, TAZ_peaks)
length(YAP_overlap_TAZ_peaks)


# install Vennerable from :
#     devtools::install_github("js229/Vennerable")
# It may require you to download additional packages such as:
#     install("RBGL")
library(Vennerable)

n_YAP <- length(YAP_peaks)
n_TAZ <- length(TAZ_peaks)
n_overlap <- length(YAP_overlap_TAZ_peaks)

# Figure 1a :  venn-diagram showing the overlapping peak number of YAP and TAZ 
venn_data <- Venn(SetNames = c("YAP", "TAZ"),
                  Weight = c(
                    "10" = n_YAP - n_overlap , # unique to A n(A) not "n only (A)"
                    "01" = n_TAZ - n_overlap, # unique to B
                    "11" = n_overlap # Intersection
                  ))

plot(venn_data)

# Figure 1b:  venn-diagram showing the overlapping peak number of YAP/TAZ and TEAD4
TEAD4_peaks <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/TEAD4_peak/TEAD4_peaks.narrowPeak"))
TEAD4_peaks

YAP_overlap_TAZ_peaks_overlap_TEAD4 <- subsetByOverlaps(YAP_overlap_TAZ_peaks, TEAD4_peaks)

n_YAP_TAZ <- length(YAP_overlap_TAZ_peaks)
n_TEAD4 <- length(TEAD4_peaks)
n_overlap2 <- length(YAP_overlap_TAZ_peaks_overlap_TEAD4)
n_overlap2

venn_data2 <- Venn(SetNames = c("YAP/TAZ", "TEAD4"),
                   Weight = c(
                     "10" = n_YAP_TAZ-n_overlap2,
                     "01" = n_TEAD4-n_overlap2,
                     "11" = n_overlap2
                   ))

plot(venn_data2)




# Figure 1c: Position of TEAD4 peak summits relative to the summits of the overlapping YAP/TAZ peaks, in a 500 bp window surrounding the summit of YAP/TAZ peaks.

TAZ_summit <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/TAZ_peak/TAZ_summits.bed"))
TAZ_summit<- TAZ_summit[TAZ_summit$name %in% TAZ_overlap_YAP_peaks$name]

TAZ_summit
TAZ_overlap_YAP_peaks

TEAD4_summit <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/TEAD4_peak/TEAD4_summits.bed"))
TEAD4_summit

# moving TAZ peaks in 500bop window 
TAZ_500bp_window <- resize(TAZ_summit,width = 500, fix = "center")
hits <- findOverlaps(TEAD4_summit, TAZ_500bp_window)
hits

summit_distance <- distance(TEAD4_summit[queryHits(hits)], TAZ_summit[subjectHits(hits)])
table(summit_distance)

TEAD4_summit[queryHits(hits)][summit_distance ==0]
TAZ_summit[subjectHits(hits)][summit_distance ==0]


# The built-in distance function returns the pair-wise distance in absolute value.
# Let’s revise it to return negative values when TEAD4 summit precede the TAZ summit and
# positive values when TEAD4 summit follows TAZ summit
signed_distance <- function(A,B) {
  dist <- distance(A, B)
  sign <- ifelse(start(A)< start(B), -1, 1)
  dist * sign
}

library(dplyr)
library(ggplot2)

summit_distance <- signed_distance(TEAD4_summit[queryHits(hits)],
                                   TAZ_summit[subjectHits(hits)])

distance_df <- table(summit_distance) %>%
  tibble::as_tibble()

summit_distance
distance_df

# plotting distance_df
distance_df %>%
  ggplot(aes(x=summit_distance, y=n)) +
  geom_point()

# Summit distance is reordered
distance_df %>%
  mutate(summit_distance = as.numeric(summit_distance)) %>%
  arrange(summit_distance) %>%
  ggplot(aes(x=summit_distance, y=n)) +
  geom_point()

# changing the plot from points to line
distance_df %>%
  mutate(summit_distance = as.numeric(summit_distance)) %>%
  arrange(summit_distance) %>%
  ggplot(aes(x=summit_distance, y=n)) +
  geom_line()

# The plot looks too wigglely. Let’s smooth it by average the number of peaks per 5 bp bin.
df_binned <- distance_df %>%
  mutate(summit_distance = as.numeric(summit_distance)) %>%
  arrange(summit_distance) %>%
  mutate(bin = floor(summit_distance /5)*5) %>%  # Creating bins by grouping every 5bp
  group_by(bin) %>%
  summarize(n = mean(n, na.rm = TRUE)) # Calculate average 'n' for each bin

print(df_binned)

df_binned %>%
  ggplot(aes(x=bin, y = n)) +
  geom_line() +
  scale_x_continuous(breaks = c(-250, 0, 250)) +
  xlab("distance to the summit \nof TAZ peaks (bp)") +
  ylab("peak density") +
  theme_classic(base_size = 14) 





















