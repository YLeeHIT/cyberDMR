library(GenomicRanges)
# library(UpSetR)
library(ComplexUpset)
library(ggVennDiagram)

cyber <- read.table("~/project/methDmr/real-data/prostate_cancer/GSE158927/result/cyberDMR.bed", header = FALSE)
metilene <- read.table("~/project/methDmr/real-data/prostate_cancer/GSE158927/result/metilene.bed", header = FALSE)
bsmooth <- read.table("~/project/methDmr/real-data/prostate_cancer/GSE158927/result/bsmooth.bed", header = FALSE)
home <- read.table("~/project/methDmr/real-data/prostate_cancer/GSE158927/result/home.bed", header = FALSE, sep = '\t')
benchmark <- read.table("~/project/methDmr/real-data/prostate_cancer/GSE158927/result/benchmark.bed", header = FALSE, sep = '\t')


gr_cyber <- GRanges(seqnames = cyber$V1, ranges = IRanges(start = cyber$V2, end = cyber$V3))
gr_metilene <- GRanges(seqnames = metilene$V1, ranges = IRanges(start = metilene$V2, end = metilene$V3))
gr_bsmooth <- GRanges(seqnames = bsmooth$V1, ranges = IRanges(start = bsmooth$V2, end = bsmooth$V3))
gr_home <- GRanges(seqnames = home$V1, ranges = IRanges(start = home$V2, end = home$V3))
gr_benchmark <- GRanges(seqnames = benchmark$V1, ranges = IRanges(start = benchmark$V2, end = benchmark$V3))

all_regions <- reduce(c(gr_cyber, gr_metilene, gr_bsmooth, gr_home, gr_benchmark))
# all_regions <- unique(c(gr_cyber, gr_metilene, gr_bsmooth, gr_home, gr_benchmark))

# length(gr_cyber)
sum(df$cyberDMR)
sum(df$Metilene)
sum(df$HOME)
sum(df$BSmooth)
sum(df$Benchmark)


df <- data.frame(
  cyberDMR = countOverlaps(all_regions, gr_cyber) > 0,
  Metilene = countOverlaps(all_regions, gr_metilene) > 0,
  BSmooth = countOverlaps(all_regions, gr_bsmooth) > 0,
  HOME = countOverlaps(all_regions, gr_home) > 0,
  Benchmark = countOverlaps(all_regions, gr_benchmark) > 0
)

upset(df, intersect = c("cyberDMR", "Metilene", "BSmooth", "HOME", "Benchmark"),
      name = "DMR detection overlap")


upset(df, intersect = c("cyberDMR", "Benchmark"))


p_upset <- upset(
  df,
  intersect = c("cyberDMR", "Metilene", "BSmooth", "HOME", "Benchmark"),
  name = "DMR detection overlap",
  base_annotations = list('Intersection size' = intersection_size(text = list(size = 4)) +
      coord_cartesian(ylim = c(0,50000))
  ),
  mode = "intersect",
  width_ratio = 0.4
) +
  theme(  
    axis.text.x = element_blank(), # x轴文本
    axis.text.y = element_text(size = 14), # y轴文本
    strip.text = element_text(size = 16), # 交集组合标签
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 

p_upset

ggsave("/home/user/liyang/project/methDmr/plot/upset_9.6_3.8_test.jpg",p_upset, width =9.6, height = 3.8, dpi = 300,units = "in")





table(df$cyberDMR, df$Benchmark)
table(df$cyberDMR, df$Metilene)
table(df$cyberDMR, df$HOME)
table(df$cyberDMR, df$BSmooth)



# venn --------------------------------------------------------------------

venn_list <- list(
  Benchmark = which(df$Benchmark),
  cyberDMR = which(df$cyberDMR),
  Metilene = which(df$Metilene),
  BSmooth = which(df$BSmooth),
  HOME = which(df$HOME)
)

p_venn <- ggVennDiagram(venn_list, label_alpha = 0, label = "count", label_size = 4, set_size = 4,
                        edge_size = 0.5)  +
  theme(legend.position = "none",) +
  # scale_fill_gradientn(colors = c("#f7fbff", "#6baed6", "#1868b2"))
  scale_fill_gradientn(colors = c("grey", "grey", "grey"))

p_venn

ggsave("/home/user/liyang/project/methDmr/plot/venn_3_2.5_plot3.jpg",p_venn, width =3.5, height = 2.5, dpi = 300,units = "in")

