# Setting the working directory -------------------------------------------
# For the laptop
setwd("C:/Users/User/Documents/Lab/crispr/guide_seq2/")
# For the office PC
setwd("S:/Misha/13 - Data Analysis/hbb.htgts/guide.seq")


# Libraries ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggrepel)


# Data import -------------------------------------------------------------
files <- list.files("bed_files/", pattern = ".bed", full.names = TRUE)

# Merging the mappings from different files into a single dataframe via a loop.
maps <- data.frame()
for (i in files) {
  if (file.info(i)$size > 0) {
    maps.tmp <- read.table(i, header = FALSE,
                           sep = "\t",
                           stringsAsFactors = FALSE,
                           fill = TRUE)
    # Due to the unequal column number, some rows "wrap" into two rows. To
    # counteract this behaviour, I filter the wrapped tails of those rows.
    maps.tmp <- maps.tmp %>% filter(str_detect(V1, "chr"))
    # Extracting the sample name from the filename.
    filename <- unlist(str_split(i, pattern = "/"))[2]
    filename <- unlist(str_split(filename, pattern = "\\."))[1]
    dir <- unlist(str_split(filename, pattern = "_"))[1]
    #pcr <- unlist(str_split(filename, pattern = "_"))[2]
    sample <- unlist(str_split(filename, pattern = "_"))[2]
    rep <- unlist(str_split(filename, pattern = "_"))[3]
    # Assigning new columns to store the sample features.
    maps.tmp$dir <- dir; maps.tmp$sample <- sample
    maps.tmp$rep <- rep
    # Appending to the master table.
    maps <- maps %>% bind_rows(maps.tmp)

  }
}
# Removing the entities I no longer need.
rm(maps.tmp, i, filename, sample, rep, dir, files)
# Removing the "single" end analysis.
#maps <- maps %>% filter(mode == "paired")
#maps <- maps %>% mutate(rep = mode) %>% select(-mode)


# I need just the first columns for now.
maps <- maps %>% select(V1, V2, V3, V4, V6, dir, sample, rep)
colnames(maps) <- c("chrom", "start", "end", "read.id",
                    "strand", "dir", "sample", "rep")

# Collapsing reads into the contigs.
contigs <- maps %>% group_by(read.id) %>%
  summarise(chrom.check = n_distinct(chrom),
            chrom = nth(chrom,1),
            start = min(start),
            end = max(end),
            sample.check = n_distinct(sample),
            sample = nth(sample, 1),
            dir = nth(dir, 1),
            rep = nth(rep, 1)) %>%
  ungroup() %>%
  mutate(length = end - start + 1)
# The control values should be 0.
sum(contigs$chrom.check > 1)
sum(contigs$sample.check > 1)
# If they are equal to 0, remove the columns.
if (sum(contigs$chrom.check > 1) == 0) {
  contigs <- contigs %>% select(-chrom.check)
} else {
  contigs.ambig.map <- contigs %>% filter(chrom.check > 1)
  contigs <- contigs %>% filter(chrom.check == 1) %>%
    select(-chrom.check)
  maps.ambig <- maps %>%
    filter(read.id %in% contigs.ambig.map$read.id)
}
if (sum(contigs$sample.check > 1) == 0) {
  contigs <- contigs %>% select(-sample.check)
}
# Checking whether the reads mapped in vicinity to each other. If the distance
# between two reads of the same contig is more than 2000, put the read in a
# separate dataframe.
ggplot(contigs) + geom_histogram(aes(x = length))
contigs.long <- contigs %>% filter(length >= 2500)
maps.long <- maps %>%
  filter(read.id %in% contigs.long$read.id)
maps.ambig$comment <- "Mapped to diff chroms"
maps.long$comment <- "Mapped to distant loci"
maps.suspicious <- bind_rows(maps.ambig, maps.long)
rm(maps.ambig, maps.long, contigs.ambig.map, contigs.long)
maps.suspicious <- maps.suspicious %>%
  mutate(length = end - start + 1)
write.table(maps.suspicious, "ambiguous.mappings.071020.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

contigs <- contigs %>%
  filter(length < 2500, !str_detect(chrom, "_"))
contigs <- contigs %>% unite(sampledirrep, c(sample, dir, rep))
contigs.agg <- contigs %>% group_by(sampledirrep) %>%
  summarise(paired.read.count = n())
contigs <- contigs %>% separate(sampledirrep, c("sample", "dir", "rep"))
#maps$mapping.id <- seq(1, nrow(maps))

# The "on-target" read is defined as a read mapped to a window 1 kb long
# surrounding the integration site in the HBB locus.
# maps$target <- "off"
# maps$target[maps$chrom == "chr11" & maps$start >= 5226340 &
#               maps$end <= 5227340] <- "on"

# Collapsing the mapping into the reads.
# maps.read <- maps %>% group_by(sample,pcr,mode,read.id) %>%
#   summarise(start = min(start), end = max(end), target = nth(target,1),
#             target.conflict = n_distinct(target))
# Target conflict column shows you whether there are reads mapped to the target
# AND to smth else. It should be 1 in all reads, so the output of the following
# should be zero.
# sum(maps.read$target.conflict == 2)

# Summarising the reads for each sample
# maps.read.agg <- maps.read %>% group_by(sample) %>%
#   summarise(reads = n(),
#             reads.on.target = sum(target == "on"))
# Transforming the number of the on-target reads into the percentage.
# maps.read.agg <- maps.read.agg %>%
#   mutate(reads.on.target.perc = reads.on.target*100/reads)
# Summarising the mappings, not the reads, for each sample.
#maps.agg <- maps %>% group_by(sample, pcr, mode) %>%
#  summarise(reads = n_distinct(read.id),
#            mappings = n(),
#            on.target = sum(chrom == "chr11" & end >= 5226340 & end <= 5227340))
# Transforming the number of the on-target mappings into the percentage.
#maps.agg <- maps.agg %>%
#  mutate(on.target.perc = on.target*100/mappings)
# Merging the summary tables.
#maps.agg <- maps.agg %>%
#  left_join(select(maps.read.agg, c("sample", "pcr", "mode", "reads.on.target", "reads.on.target.perc")),
#            by = c("sample", "pcr", "mode"))
# Saving the table.
# write.table(maps.read.agg, file = "ontargetreads_110920.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#
# write.table(maps.agg, file = "ontargetreads.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(maps, file = "all_mappings.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Off-target analysis -----------------------------------------------------
#maps.agg <- maps.agg %>% ungroup() %>%
#  mutate(sample = paste(sample, pcr, mode, sep = "_")) %>%
#  select(-c("pcr", "mode"))
#maps <- maps %>%
#  mutate(sample = paste(sample, pcr, mode, sep = "_")) %>%
#  select(-c("pcr", "mode"))

#maps <- maps %>%
#  filter(!str_detect(sample, "single")) %>%
#  separate(sample, c("sample", "pcr", "mode"), sep = "_")
#maps$sample %>% unique()

#samples <- samples[str_detect(samples, "paired")]

hg38.genes <-
  read.table("C:/Users/User/Documents/Lab/crispr/guide_seq/hg38_genes_NCBI.txt", header = TRUE,
             sep = "\t", stringsAsFactors = FALSE)
hg38.genes <- hg38.genes %>%
  filter(!str_detect(chrom, "_") & str_detect(name, "NM")) %>%
  group_by(name2) %>%
  summarise(chrom = nth(chrom, 1), start = min(txStart),
            end = max(txEnd)) %>%
  mutate(gene.name = name2) %>% select(-name2)

contigs <- contigs %>% mutate(cells = sample) %>%
  unite(sample, c("sample", "dir", "rep"))

samples <- contigs$sample %>% unique()

for(i in samples){
  contigs.tmp <- contigs %>% filter(sample == i)
  contigs.tmp <- contigs.tmp %>% arrange(chrom, end)

  # Assigning the groups based on the distance. The group ID changes when the
  # neighbour is further away than 1000 bp.
  group.id <- 1
  contigs.tmp$group.id <- NA
  contigs.tmp$group.id[1] <- 1
  if(nrow(contigs.tmp) > 1){
    for (j in seq(2, nrow(contigs.tmp))) {
      if (is.na(contigs.tmp$start[j]) | is.na(contigs.tmp$end[j])) {
        next
      }
      if (contigs.tmp$chrom[j] == contigs.tmp$chrom[j-1]){
        if (min(abs(contigs.tmp$start[j] - contigs.tmp$end[j-1]),
                abs(contigs.tmp$start[j-1] - contigs.tmp$end[j]),
                abs(contigs.tmp$start[j] - contigs.tmp$start[j-1]),
                abs(contigs.tmp$end[j] - contigs.tmp$end[j-1])) <= 1000){

        } else {
          group.id <- group.id + 1
        }
      } else {
        group.id <- group.id + 1
      }
      contigs.tmp$group.id[j] <- group.id
    }
    rm(j)
    print(sprintf("%d clusters were identified in this sample.", group.id))
  }


  # Summarising the table based on those clusters.
  contigs.tmp.agg <- contigs.tmp %>% group_by(group.id) %>%
    summarise(mappings = n(),
              chrom = nth(chrom, 1),
              start = min(start),
              end = max(end))
  contigs.tmp.agg <- contigs.tmp.agg %>%
    mutate(mappings.perc = mappings * 100 / sum(mappings))

  contigs.tmp.agg$genes <- NA
  for (k in seq(nrow(contigs.tmp.agg))) {
    target.chrom <- contigs.tmp.agg$chrom[k]
    target.start <- contigs.tmp.agg$start[k]
    target.end <- contigs.tmp.agg$end[k]
    hg38.genes.this.chrom <- hg38.genes %>% filter(chrom == target.chrom)
    hg38.genes.this.region <- hg38.genes.this.chrom %>%
      filter((start >= target.start & start <= target.end) |
               (end >= target.start & end <= target.end) |
               (start <= target.start & end >= target.end))
    if(nrow(hg38.genes.this.region) > 0){
      contigs.tmp.agg$genes[k] <- paste0(hg38.genes.this.region$gene.name, collapse = "|")
    }
    rm(hg38.genes.this.region)
  }
  rm(target.chrom, target.start, target.end, k)
  contigs.tmp.agg <- contigs.tmp.agg %>% arrange(desc(mappings))
  # Naming the non-genic sites.
  contigs.tmp.agg <- contigs.tmp.agg %>%
    mutate(site = paste(chrom,": ", round(start/1000000), "Mb", sep = ""))
  contigs.tmp.agg$genes[is.na(contigs.tmp.agg$genes)] <- contigs.tmp.agg$site[is.na(contigs.tmp.agg$genes)]
  contigs.tmp.agg <- contigs.tmp.agg %>%
    mutate(site = genes) %>% select(-genes)

#   # Filtering the random integration.
#   # Warning: the threshold is completely arbitrary!
#   contigs.tmp.agg$random.integration <- TRUE
#   # A hotspot is deemed random if it occupies less than 0.1% of the reads.
#   contigs.tmp.agg$random.integration[contigs.tmp.agg$mappings.perc >= 0.1] <- FALSE
#
#   # Pie chart ---------------------------------------------------------------
#   # Specify how many top sites should be shown in the pie chart.
#   show.top <- 4 # Hence the pie chart will show 4 sites + others + random.
#   contigs.tmp.agg.random <- contigs.tmp.agg %>% filter(random.integration)
#   contigs.tmp.agg.nonrandom <- contigs.tmp.agg %>% filter(!random.integration)
#   contigs.tmp.agg.nonrandom.top <- contigs.tmp.agg.nonrandom[1:show.top,]
#   contigs.tmp.agg.nonrandom.bottom <-
#     contigs.tmp.agg.nonrandom[(show.top+1):nrow(contigs.tmp.agg.nonrandom),]
#   contigs.tmp.agg.nonrandom.bottom <- contigs.tmp.agg.nonrandom.bottom %>%
#     group_by(1) %>%
#     summarise(group.id = "Other", mappings = sum(mappings),
#               chrom = NA,
#               end = NA, end = NA,
#               mappings.perc = sum(mappings.perc),
#               site = "Other", random.integration = FALSE) %>%
#     select(-1)
#   contigs.tmp.agg.random <- contigs.tmp.agg.random %>%
#     group_by(1) %>%
#     summarise(group.id = "Random", mappings = sum(mappings),
#               chrom = NA,
#               end = NA, end = NA,
#               mappings.perc = sum(mappings.perc),
#               site = "Random", random.integration = TRUE) %>%
#     select(-1)
#   contigs.tmp.agg.nonrandom.top <- contigs.tmp.agg.nonrandom.top %>%
#     mutate(group.id = as.character(group.id)) %>%
#     bind_rows(contigs.tmp.agg.nonrandom.bottom, contigs.tmp.agg.random)
#
#   sites <- contigs.tmp.agg.nonrandom.top$site
#
#   contigs.tmp.agg.nonrandom.top <- contigs.tmp.agg.nonrandom.top %>%
#     mutate(text.pos = cumsum(mappings) - mappings/2,
#            site = factor(site, levels = sites))
#
#   fill.colors <- c(rgb(seq(1,0,-1/(show.top-1)), 0.5, 0.5),
#                    "grey", "white")
#   names(fill.colors) <- sites
# # Generating some ugly pie-charts.
#   ggplot(contigs.tmp.agg.nonrandom.top)+
#     geom_bar(aes(x = 1, y = mappings, fill = site), stat = "identity",
#              position = position_stack(reverse = TRUE),
#              color = "black")+
#     coord_polar(theta = "y", direction = -1)+
#     theme_classic()+
#     geom_label_repel(aes(x = 1, y = text.pos,
#                          label = paste0(site,", ", round(mappings.perc), "%")),
#                      force = 5,
#                      point.padding = 0,
#                      nudge_x = 0.9,
#                      segment.colour = "black",
#                      segment.size = 1)+
#     scale_fill_manual(values = fill.colors)+
#     theme(axis.title = element_blank(),
#           axis.text = element_blank(),
#           axis.line = element_blank(),
#           axis.ticks = element_blank(),
#           legend.position = "none")
#
#   ggsave(paste0(i, "_", format(Sys.time(),"%d%m%y"), ".png"),
#          width = 4, height = 4, device = "png")
  write.table(contigs.tmp.agg, file = paste0(i, "_", format(Sys.time(),"%d%m%y"), ".txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
#  rm(contigs.tmp, contigs.tmp.agg, contigs.tmp.agg.random, contigs.tmp.agg.nonrandom.top,
#     contigs.tmp.agg.nonrandom, contigs.tmp.agg.nonrandom.bottom, fill.colors,
#     show.top)
}


# Overlap with the CRISPR-Gold predictions --------------------------------
predictions <- read.table("C:/Users/User/Documents/Lab/crispr/guide_seq/crispr_gold_predictions.txt",
                          header= TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
predictions <- predictions %>%
  separate(position, into = c("chrom", "start", "end"), sep = "\\:") %>%
  mutate(start = as.numeric(start), end = as.numeric(end))
predictions <- predictions %>% filter(risk != "target")
predictions$off.target.id <- seq(nrow(predictions))
window <- 1000
predictions <- predictions %>%
  mutate(start.precise = start, end.precise = end,
         middle = (start + end)/2,
         start = floor(middle - window/2), end = start + window) %>%
  select(-middle, -top.offtarget.sites, -strand)

files <- list.files("clusters/", pattern = ".txt", full.names = TRUE)

# Importing the clustered data.
clusters <- data.frame()
for (i in files) {
  if (file.info(i)$size > 0) {
    clusters.tmp <- read.table(i, header = TRUE,
                           sep = "\t",
                           stringsAsFactors = FALSE,
                           fill = TRUE)

    # Extracting the sample name from the filename.
    filename <- unlist(str_split(i, pattern = "/"))[2]
    filename <- unlist(str_split(filename, pattern = "\\."))[1]
    sample <- unlist(str_split(filename, pattern = "_"))[1]
    dir <- unlist(str_split(filename, pattern = "_"))[2]
    rep <- unlist(str_split(filename, pattern = "_"))[3]
    #mode <- unlist(str_split(filename, pattern = "_"))[3]
    # Assigning new columns to store the sample features.
    clusters.tmp$sample <- sample; clusters.tmp$dir <- dir; clusters.tmp$rep <- rep
    #clusters.tmp$mode <- mode
    # Appending to the master table.
    clusters <- clusters %>% bind_rows(clusters.tmp)

  }
}
# Removing the entities I no longer need.
rm(clusters.tmp, i, filename, sample, pcr, files, rep, dir)

clusters$predictoff <- NA
for (k in seq(nrow(clusters))) {
  target.chrom <- clusters$chrom[k]
  target.start <- clusters$start[k]
  target.end <- clusters$end[k]
  predictions.this.chrom <- predictions %>% filter(chrom == target.chrom)
  predictions.this.region <- predictions.this.chrom %>%
    filter((start >= target.start & start <= target.end) |
             (end >= target.start & end <= target.end) |
             (start <= target.start & end >= target.end))
  if(nrow(predictions.this.region) > 0){
    clusters$predictoff[k] <- paste0(predictions.this.region$off.target.id, collapse = "|")
  }
  rm(predictions.this.region)
}
rm(target.chrom, target.start, target.end, k)

write.table(clusters,
            file = "clustered_reads_offtargets_071020.txt",
            sep = "\t",
            row.names = FALSE, quote = FALSE)

clusters <- clusters %>% mutate(length = end - start)

# Clusters overlap --------------------------------------------------------
unique(clusters$sample)

#clusters <- lam.table %>% select(-AY, -BY, -AX, -BX, - rel.start, -rel.end)
#clusters <- clusters %>% select(-pcr)
#clusters <- clusters %>% unite(sample, c("sample", "pcr"), sep = "_")
#clusters <- clusters %>% separate(sample,c("sample", "day", "direction"), sep = "-")
clusters <- clusters %>% unite(sample, c("sample", "rep"), sep = "-")

clusters$sample %>% unique()

clusters.new.pooled <- data.frame()
for(i in unique(clusters$sample)){
  clusters.5 <- clusters %>% filter(sample == i, dir == "5GSP")
  clusters.3 <- clusters %>% filter(sample == i, dir == "3GSP")
  clusters.5$start3 <- NA
  clusters.5$end3 <- NA
  clusters.5$mappings3 <- NA
  clusters.5$mappings.perc3 <- NA

  for (k in seq(nrow(clusters.5))) {
    chrom5 <- clusters.5$chrom[k]
    start5 <- clusters.5$start[k]-500
    end5 <- clusters.5$end[k]+500
    clusters.3.this.chrom <- clusters.3 %>% filter(chrom == chrom5)
    clusters.3.this.region <- clusters.3.this.chrom %>%
      filter( (start >= start5 & start <= end5) |
                (end >= start5 & end <= end5) |
                (start <= start5 & end >= end5))
    if(nrow(clusters.3.this.region) > 1){
      warning("Ambiguous overlap!")
      clusters.5$start3[k] <-  min(clusters.3.this.region$start)
      clusters.5$end3[k] <-  max(clusters.3.this.region$end)
      clusters.5$mappings3[k] <-  sum(clusters.3.this.region$mappings)
      clusters.5$mappings.perc3[k] <- mean(clusters.3.this.region$mappings.perc)
    } else if (nrow(clusters.3.this.region) == 1) {
      clusters.5$start3[k] <- clusters.3.this.region$start
      clusters.5$end3[k] <- clusters.3.this.region$end
      clusters.5$mappings3[k] <- clusters.3.this.region$mappings
      clusters.5$mappings.perc3[k] <- clusters.3.this.region$mappings.perc
    }
  }
  clusters.new.pooled <- clusters.new.pooled %>% bind_rows(clusters.5)

}

clusters.new.pooled <- clusters.new.pooled %>% filter(!is.na(start3))
clusters.new.pooled.1kb <- clusters.new.pooled.1kb %>% filter(!is.na(start3))

write.table(clusters.new.pooled, file = "overlap_071020.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(clusters.new.pooled.1kb, file = "overlap.1kb.window.pooled_110920.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

contigs.agg <- contigs.agg %>%
  separate(sample, c("sample", "day", "direction"), sep = "-") %>%
  unite(sample, c("sample", "day"), sep = "-")

clusters.new.pooled.1kb <- clusters.new.pooled.1kb %>%
  left_join(filter(contigs.agg, direction == "5"), by = "sample") %>%
  mutate(paired.read.count.5 = paired.read.count) %>%
  select(-paired.read.count) %>%
  left_join(filter(contigs.agg, direction == "3"), by = "sample") %>%
  mutate(paired.read.count.3 = paired.read.count) %>%
  select(-paired.read.count)
#
# clusters.olap.cas9 <- clusters.olap.cas9 %>%
#   mutate(sample = factor(sample, levels = rev(c("3-1-N60", "4-1-N60", "4-T-N60",
#                                             "3-1-TD", "4-1-TD", "4-T-TD",
#                                             "3-1-TD2", "4-1-TD2", "4-T-TD2"))),
#          site = factor(site, levels = rev(c("HBB", "chr9: 102Mb", "chr12: 124Mb",
#                                         "HBD", "chr3: 93Mb"))))

# ggplot(clusters.olap.cas9%>% filter(str_detect(sample, "1-TD2")) )+
#   geom_bar(aes(x = site, y = mapping.perc.mean, fill = sample), color = "black",
#            stat = "identity", position = position_dodge2(preserve = "single"))+
#   geom_point(aes(x = site, y = mappings.perc3, fill = sample), shape = 21,
#              stat = "identity", position = position_dodge2(width = 1))+
#   geom_point(aes(x = site, y = mappings.perc, fill = sample), shape = 23,
#              stat = "identity", position = position_dodge2(width = 1))+
#
#   theme_classic()+
#   coord_flip()+
#   scale_y_continuous()
#
# ggsave("mappings1.svg", device = "svg", width = 6, height = 4)
#
# write.table(clusters.olap.cas9, file = "clusters.olaped.txt", sep = '\t',
#             quote = FALSE, row.names = FALSE)
#
# clusters <- read.table("clusters.olaped.txt", sep = "\t", header = TRUE,
#                        stringsAsFactors = FALSE)
#
#
# overlap.perc <- function(start1, end1, start2, end2){
#   # Which one is left
#   if(start2 < start1){
#     # Then the second is on the left, swapping the coordinates
#     start.tmp <- start2; start2 <- start1; start1 <- start.tmp
#     end.tmp <- end2; end2 <- end1; end1 <- end.tmp
#   }
#   overlap.length <- end1 - start2 + 1
#   if(end2 < end1){
#     overlap.length <- overlap.length - (end1 - end2)
#   }
#
#   return(overlap.length)
# }
#



chrlen <- read.table("C:/Users/User/Documents/Lab/LongBCRs/bioinformatics/lair/input/GRCh38_chrlength.txt",
                             header = TRUE, stringsAsFactors = FALSE)
# Preparation of the chrlen data
chrlen <- dplyr::select(chrlen, Chromosome, Length.bp)
chrlen$length <- as.numeric(str_replace_all(chrlen$Length.bp, ',', ''))
chrlen$relative.length <- chrlen$length/max(chrlen$length)
chrlen$chrom <- str_c('chr', chrlen$Chromosome)
chrlen <- chrlen %>% select(chrom, length, relative.length)

centromeres <- read.table("C:/Users/User/Documents/Lab/LongBCRs/bioinformatics/lair/input/hg38_centromeres.txt", header = TRUE, stringsAsFactors = FALSE)
centromeres <- centromeres %>% group_by(chrom) %>%
  summarise(start = min(chromStart),
            end = max(chromEnd)) %>%
  ungroup() %>%
  mutate(center = (start+end)/2)

chrlen <- chrlen %>%
  left_join(select(centromeres, chrom, center), by = "chrom")
rm(centromeres)

plot.height <- 1500
plot.width <- 1000
x.margin <- 50
y.margin <- 10
# Chromosomes rectangles thickness
# Triangles parameters
tw <- 6
th <- 10
thick <- 0
x.inner <- plot.width - x.margin
y.inner <- plot.height - y.margin
inner.width <- plot.width - 2*x.margin
inner.height <- plot.height - 2*y.margin
N <- nrow(chrlen)
spacing <- inner.height/(N-1)
chrlen$element.i <- seq(N)

## Points coordinates calculation
chrlen$AX <- x.margin
chrlen$AY <- y.inner - spacing*(chrlen$element.i - 1) - thick/2
chrlen$BX <- x.margin + inner.width*chrlen$relative.length
chrlen$BY <- y.inner - spacing*(chrlen$element.i - 1) + thick/2

## Centromeres coordinates
chrlen$center.relative <- chrlen$center/max(chrlen$length)
chrlen$center.Y <- chrlen$AY
chrlen$center.X <- x.margin + inner.width*chrlen$center.relative
chrlen$length.mbp <- chrlen$length / 1000000


# length.mbp will be used as the text labels on the plot
chrlen$length.mbp <- round(chrlen$length.mbp, 0)
chrlen$length.mbp <- as.character(chrlen$length.mbp)
chrlen$length.mbp <- paste(chrlen$length.mbp, 'Mbp')
color_of_lines <- 'black'

clusters.new.pooled.1kb <- clusters.new.pooled.1kb %>%
  left_join(select(chrlen, AY, chrom), by = "chrom")
clusters.new.pooled.1kb <- clusters.new.pooled.1kb %>%
  mutate(center = (start+end)/2,
         AX = x.margin + center*inner.width/max(chrlen$length))

predictions <- predictions %>%
  left_join(select(chrlen, AY, chrom), by = "chrom") %>%
  mutate(center = (start+end)/2,
         AX = x.margin + center*inner.width/max(chrlen$length))

predictions.copy <- predictions

predictions <- predictions.copy %>%
  #arrange(desc(cfdOfftargetScore)) %>%
  #dplyr::slice(1:top_n) %>%
  filter(chrom != "chrM") %>%
  filter(!str_detect(chrom, "_")) %>%
  mutate(center = (start+end)/2,
         AX = x.margin + inner.width*center/max(chrlen$length),
         BX = AX - 4, CX = AX + 4) %>%
  mutate(BY = AY - 20, CY = BY)

predictions <- predictions %>%
  unite(A, AX, AY, sep = ';', remove = TRUE) %>%
  unite(B, BX, BY, sep = ';', remove = TRUE) %>%
  unite(C, CX, CY, sep = ';', remove = TRUE) %>%
  gather(dot, dot.coord, c("A","B","C")) %>%
  separate(dot.coord, c('dotX', 'dotY'), sep = ';', remove = TRUE) %>%
  mutate_at(c("dotX", "dotY"), as.numeric)

clusters.new.pooled.1kb <- clusters.new.pooled.1kb %>%
  mutate(read.perc.mean = (mappings.perc + mappings.perc3) / 2)

lamplot <- ggplot(data = chrlen %>% filter(chrom != "chrM"))+
  coord_cartesian(xlim = c(0, plot.width), ylim = c(0, plot.height))+
  geom_segment(aes(x = AX, y = AY, xend = BX, yend = BY), alpha = 0.2, color = "black", size = 1)+
  geom_segment(data = clusters.new.pooled.1kb %>% filter(sample %in% c("3-1", "4-1", "4-T")),
               aes(x = AX, y = AY, xend = AX, yend = AY+read.perc.mean.log*30, color = sample),
               alpha = 0.8, size = 3)+
  geom_point(aes(x = center.X, y = center.Y), size = 2, shape = 21, fill = 'deeppink2')+
  geom_text(aes(x = AX-35, y = AY+3, label = chrom), colour = "black")+
  geom_text(aes(x = BX + 50, y = AY+3, label = length.mbp), size = 3, colour = "black")+
  # geom_segment(data = lam %>% filter(reads.count > 1), aes(x = AX-1, y = AY+2.5, xend = AX-1, yend = BY+2.5),
  #              colour = "blue", size = 1)+
  geom_polygon(data = predictions %>% filter(risk != "lim?0"),
               aes(x = dotX, y = dotY, group = off.target.id, fill = risk))+
  scale_fill_manual(values = c("high" = rgb(1,0,0), "low" = rgb(1,0.5,0), "lim?0" = rgb(0.6,0.6,0.6)))+
  theme_void()
lamplot

ggsave(paste('guideseq_',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.svg', sep = ''),
       width = 15, height = 10,
       device = "svg",
       plot = lamplot)
rm(lam.table.temp)