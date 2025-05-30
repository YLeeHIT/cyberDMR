library(GenomicRanges)
library(ComplexUpset)
library(ggVennDiagram)

cyber <- read.table("../data/real_data/cyberDMR.bed", header = FALSE)
metilene <- read.table("../data/real_data/metilene.bed", header = FALSE)
bsmooth <- read.table("../data/real_data/bsmooth.bed", header = FALSE)
home <- read.table("../data/real_data/home.bed", header = FALSE, sep = '\t')
benchmark <- read.table("../data/real_data/benchmark.bed", header = FALSE, sep = '\t')

gr_cyber <- GRanges(seqnames = cyber$V1, ranges = IRanges(start = cyber$V2, end = cyber$V3))
gr_metilene <- GRanges(seqnames = metilene$V1, ranges = IRanges(start = metilene$V2, end = metilene$V3))
gr_bsmooth <- GRanges(seqnames = bsmooth$V1, ranges = IRanges(start = bsmooth$V2, end = bsmooth$V3))
gr_home <- GRanges(seqnames = home$V1, ranges = IRanges(start = home$V2, end = home$V3))
gr_benchmark <- GRanges(seqnames = benchmark$V1, ranges = IRanges(start = benchmark$V2, end = benchmark$V3))

all_regions <- reduce(c(gr_cyber, gr_metilene, gr_bsmooth, gr_home, gr_benchmark))

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

# upset --------------------------------------------------------------------

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
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size = 14), 
    strip.text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 

p_upset

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
  scale_fill_gradientn(colors = c("grey", "grey", "grey"))

p_venn
