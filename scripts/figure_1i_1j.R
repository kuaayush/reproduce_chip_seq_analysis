
# Figure 1i
# Heatmap representing YAP/TAZ/TEAD binding sites located on promoters (top) and enhancers (bottom).
# YAP, TAZ and TEAD4 peaks are ranked from the strongest to weakest signal in TAZ ChIP,
# in a window of ±1kb centered on the summit of TAZ peaks.
# H3K4me1 and H3K4me3 signal in the corresponding genomic regions is shown on the right.

# Here, Heatmaps were generated using a custom R script which considers a 2-kb window centered on peak summits and
# calculates the normalized reads density with a resolution of 50 bp.


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

YAP_summit <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/YAP_peak/YAP_summits.bed"))
YAP_summit

# Reading the histone modification peaks
H3K4me1 <- import(here("/Users/aayushojha/Documents/chipseq/data/public_data/H3K4me1.bed"))
H3K4me3 <- import(here("/Users/aayushojha/Documents/chipseq/data/public_data/H3K4me3.bed"))
H3K27ac <- import(here("/Users/aayushojha/Documents/chipseq/data/public_data/H3K27ac.bed"))

# Define enhancers and promoters
enhancers <- subsetByOverlaps(H3K4me1, H3K4me3, invert = TRUE)
promoters <- subsetByOverlaps(H3K4me3, H3K4me1, invert = TRUE)

# annotate the YAP1/TAZ/TEAD4 peaks:
YAP_enhancers <- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4, enhancers)
YAP_promoters <- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4, promoters)

YAP_enhancers$name %>% head()


YAP_summit_enhancer <- YAP_summit[YAP_summit$name %in% YAP_enhancers$name]
YAP_summit_promoter <- YAP_summit[YAP_summit$name %in% YAP_promoters$name]

# combine these two:
anchors <- c(YAP_summit_enhancer, YAP_summit_promoter)



# import the bi8gwig files 
YAP_bw <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/YAP.bw"))
TAZ_bw <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/TAZ.bw"))
TEAD4_bw <- import(here("/Users/aayushojha/Documents/chipseq/data/fastq/TEAD4.bw"))
YAP_bw

H3K4me1_bw <- import(here("/Users/aayushojha/Documents/chipseq/data/public_data/H3K4me1.bw"))
# H3K4me3_bw <- import(here("/Users/aayushojha/Documents/chipseq/data/public_data/H3K4me3.bw"))
# H3K27ac_bw <- import(here("/Users/aayushojha/Documents/chipseq/data/public_data/H3K27ac.bw"))
# 
# H3K4me1_bw
# H3K4me3_bw
# H3K27ac_bw




# Quantifying the signals in the bins 
# BiocManager::install("EnrichedHeatmap")
library(EnrichedHeatmap)
# extend 1000 bp on each side and use 50bp bin
mat1<- normalizeToMatrix(YAP_bw, anchors, value_column = "score",
                         extend= 1000, mean_mode = "w0", w=50)

mat2<- normalizeToMatrix(TAZ_bw, anchors, value_column = "score",
                         extend= 1000, mean_mode = "w0", w=50)

mat3<- normalizeToMatrix(TEAD4_bw, anchors, value_column = "score",
                         extend= 1000, mean_mode = "w0", w=50)

mat4<- normalizeToMatrix(H3K4me1_bw, anchors, value_column = "score",
                         extend= 1000, mean_mode = "w0", w=50)

# mat5<- normalizeToMatrix(H3K4me3_bw, anchors, value_column = "score",
#                          extend= 5000, mean_mode = "w0", w=100)
# 
# mat6<- normalizeToMatrix(H3K27ac_bw, anchors, value_column = "score",
#                          extend= 5000, mean_mode = "w0", w=100)

dim(mat1)
dim(mat2)
dim(mat3)
dim(mat4)
# dim(mat5)
# dim(mat6)

mat1[1:5, 1:40] 

# before we assign color, we are checking the data range 
quantile(mat1, c(0.1,0.25,0.5,0.9,1))
quantile(mat2, c(0.1,0.25,0.5,0.9,1))
quantile(mat3, c(0.1,0.25,0.5,0.9,1))
quantile(mat4, c(0.1,0.25,0.5, 0.6, 0.7, 0.8, 0.92, 0.93, 0.96, 0.98,1))
# quantile(mat5, c(0.1,0.25,0.5,0.9,1))
# quantile(mat6, c(0.1,0.25,0.5,0.9,1))


col_fun<- circlize::colorRamp2(c(0, 20), c("white", "red"))
col_fun_H3K4me1 <- circlize::colorRamp2(c(0, 5), c("white", "red"))
# col_fun2<- circlize::colorRamp2(c(0, 5), c("white", "red"))
# 
# col_fun_H3K4me3 <- circlize::colorRamp2(c(0, 13), c("white", "red"))
# 


# Separate promoters and enhancers
partition <- c(rep("promoters", length(YAP_promoters)),
               rep("enhancers", length(YAP_enhancers))
               )
partition <- factor(partition, levels = c("promoters", "enhancers"))
partition

partition_hp <-Heatmap(partition, col=structure(2:3, names = c("promoters", "enhancers")),
                       name = "partition", 
                       show_row_names = FALSE, width = unit(3, "mm")
                       )

partition_hp


ht_list <- partition_hp + 
  EnrichedHeatmap(mat1, pos_line = FALSE, column_title = "YAP", name = "YAP", col = col_fun)+ 
  EnrichedHeatmap(mat2, pos_line = FALSE, column_title = "TAZ", name = "TAZ", col = col_fun) + 
  EnrichedHeatmap(mat3, pos_line = FALSE, column_title = "TEAD4", name = "TEAD4", col = col_fun) +
  EnrichedHeatmap(mat4, pos_line = FALSE, column_title = "H3K4me1", name = "H3K4me1", col = col_fun_H3K4me1)
  # EnrichedHeatmap(mat5, pos_line = FALSE, column_title = "H3K4me3", name = "H3K4me3", col = col_fun_H3K4me3) +
  # EnrichedHeatmap(mat6, pos_line = FALSE, column_title = "H3K27ac", name = "H3K4me1", col = col_fun2)

  
library(Cairo)
draw(ht_list, split = partition, main_heatmap = 2)





library(dplyr)
library(tidyr)
library(ggplot2)

# Figure 1j 
# Bimodal distribution of H3K4me1 signal around the summit of YAP/TAZ and TEAD4 peaks.

YAP_mean <- colMeans(mat1)
TAZ_mean <- colMeans(mat2)
TEAD4_mean <- colMeans(mat3)
H3K4me1_mean <- colMeans(mat4)

YAP_mean
H3K4me1_mean

# bind_rows(YAP_mean, TAZ_mean, TEAD4_mean) %>%
#   mutate(factor = factor(c("YAP", "TAZ", "TEAD4"), levels = c("YAP", "TAZ", "TEAD4"))) %>%
#   select(factor, everything()) %>%
#   tidyr::pivot_longer(-factor) %>%
#   ggplot(aes(x=name, y=value))


bind_rows(YAP_mean, TAZ_mean, TEAD4_mean) %>%
  mutate(factor = factor(c("YAP", "TAZ", "TEAD4"), levels = c("YAP", "TAZ", "TEAD4"))) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor) %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  ggtitle("")+
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("red", "blue", "orange")) +
  theme_classic(base_size = 14) +
  ylab("RPKM") +
  xlab("Distance to peak summit")
  
# # if you include hek4me1 as well :
# bind_rows(YAP_mean, TAZ_mean, TEAD4_mean, H3K4me1_mean) %>%
#   mutate(factor = factor(c("YAP", "TAZ", "TEAD4", "H3K4me1"), levels = c("YAP", "TAZ", "TEAD4", "H3K4me1"))) %>%
#   dplyr::select(factor, everything()) %>%
#   tidyr::pivot_longer(-factor) %>%
#   mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
#   ggplot(aes(x=name, y=value)) +
#   geom_line(aes(color = factor, group=factor)) +
#   ggtitle("Bimodal distribution of H3K4me1 signal around the summit of YAP/TAZ and TEAD4 peaks") + 
#   scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
#   scale_color_manual(values = c("#DA9195", "#E07B78", "#605D7D", "orange")) +
#   theme_classic(base_size = 14) +
#   ylab("RPKM") +
#   xlab("Distance to peak summit")
