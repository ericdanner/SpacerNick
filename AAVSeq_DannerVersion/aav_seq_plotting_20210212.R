# Setting the working directory -------------------------------------------
setwd("/home/eric/Data/Spaced_Nicking/AAVSeq2_ELANE_PRF1_IL7R_2_2021/bed_files")


# Libraries ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggrepel)






# Eric testing function space -----------


# -----------------------------------------------------

# Merging the mappings from different files into a single dataframe via a loop. Add the name of the file and rep#
maps <- data.frame()
for (i in files) {
  if (file.info(i)$size > 0) {
    maps.tmp <- read.table(i, header = FALSE,
                           sep = "\t",
                           stringsAsFactors = FALSE,
                           fill = TRUE)
    maps.tmp <- maps.tmp %>% mutate(V2 = as.integer(V2))
    # Due to the unequal column number, some rows "wrap" into two rows. To
    # counteract this behaviour, I filter the wrapped tails of those rows.
    maps.tmp <- maps.tmp %>% filter(str_detect(V1, "chr"))
    # Extracting the sample name from the filename.
    filename <- unlist(str_split(i, pattern = "/"))[3]
    sample <- unlist(str_split(filename, pattern = "\\."))[1]
    # Assigning new columns to store the sample features.
    maps.tmp$sample <- sample
    # Appending to the master table.
    maps <- maps %>% bind_rows(maps.tmp)

  }
}
# Removing the entities I no longer need.
rm(maps.tmp, i, filename, sample, files)


# I need just the first columns for now.
maps <- maps %>% select(V1, V2, V3, V4, V6, sample)
colnames(maps) <- c("chrom", "start", "end", "read.id",
                    "strand", "sample")


# Collapsing reads into the continues info from the pair end read
# paired end so checking that both reads map to the same chrom and that they are in the same sample
contigs <- maps %>% group_by(read.id) %>%
  summarise(chrom.check = n_distinct(chrom), #make sure both reads have the same chrom mapping
            chrom = nth(chrom,1),
            start = min(start), #take the minimum of the two read alignments 
            end = max(end),  #takes hte maximum of the two read alignments
            sample.check = n_distinct(sample), #making sure both reads map to the same "sample"
            sample = nth(sample, 1)) %>%
  ungroup() %>% #ungroup them so its a normal chart again
  mutate(length = end - start + 1) #first create the "start and end colums before performing operations
# The control values should be 0. This is checking if any of the reads "contigs" had issues. Each true counts as 1.
sum(contigs$chrom.check > 1) #if each one passes the test, it is 1. So then its false and 0. so sum is 0.
sum(contigs$sample.check > 1)
# If they are equal to 0 and passed the check, remove the columns.
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
  contigs <- contigs %>% select(-sample.check) #if it passes check everything besides the 
}




#  CLEANING UP THE READS BEFORE CROSS TESTING
# Checking whether the reads mapped in vicinity to each other. If the distance
# between two reads of the same contig is more than 2000, put the read in a
# separate dataframe.

ggplot(contigs) + geom_histogram(aes(x = length))
contigs.long <- contigs %>% filter(length >= 2000)
maps.long <- maps %>% filter(read.id %in% contigs.long$read.id)
#maps.ambig$comment <- "Mapped to diff chroms"
maps.long$comment <- "Mapped to distant loci"
maps.suspicious <- maps.long
#maps.suspicious <- bind_rows(maps.ambig, maps.long)
rm(maps.ambig, maps.long, contigs.ambig.map, contigs.long)
maps.suspicious <- maps.suspicious %>%
  mutate(length = end - start + 1)
write.table(maps.suspicious, "global.ambiguous.mappings.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(maps.suspicious, "local.ambiguous.mappings.txt",
 #           sep = "\t", row.names = FALSE, quote = FALSE)

contigs <- contigs %>% filter(length < 2000, !str_detect(chrom, "_"))
contigs.agg <- contigs %>% group_by(sample) %>%
  summarise(paired.read.count = n())




####################################
#
#   this seems to be the major sorting function now that the reads are cleaned up
#
samples <- contigs$sample %>% unique() #make a list of all samples

for(i in samples){
  contigs.tmp <- contigs %>% filter(sample == i) # filter for reads from one sample
  contigs.tmp <- contigs.tmp %>% arrange(chrom, end) 

  # Assigning the groups based on the distance. The group ID changes when the
  # neighbour is further away than 1000 bp.
  group.id <- 1
  contigs.tmp$group.id <- NA #making an empty column
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


  contigs.tmp.agg <- contigs.tmp.agg %>% arrange(desc(mappings))
  # Naming the non-genic sites.
  contigs.tmp.agg <- contigs.tmp.agg %>%
    mutate(site = paste(chrom,": ", round(start/1000000), "Mb", sep = ""))

  write.table(contigs.tmp.agg, file = paste0(i, "_", format(Sys.time(),"%d%m%y"), ".txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
} 
  
  # I put these files into a folder called "clusters_20210210"
  
  
  
  
################################################
#
# Overlap with the CRISPR-Gold predictions 
#  
# First make a table with the off targets in a dataframe
#  

  
predictions <- read.table("ELANE_ILR7_PRF1_offtargets.txt",
                            header= TRUE, sep = "\t",
                            stringsAsFactors = FALSE)
predictions <- predictions %>%
    separate(position, into = c("chrom", "start", "end"), sep = "\\:") %>%
    mutate(start = as.numeric(start), end = as.numeric(end))

  
predictions$offtarget.id <- seq(nrow(predictions))

# window size around Cas9 break  
window <- 2000
# sstart precise is the actual break site. start and end are the window
predictions <- predictions %>%
    mutate(start.precise = start, end.precise = end,
           middle = (start + end)/2,
           start = floor(middle - window/2), end = start + window) %>%
    select(-middle, -strand) #remove middle and strand collumn
  
  

##################################
#
# Second part of adding CRISPR gold data to clusters, Importing the clustered data.

#files <- list.files("clusters_local_20210216/", pattern = ".txt", full.names = TRUE)
files <- list.files("clusters_global_20210216/", pattern = ".txt", full.names = TRUE)

################################3
# 
# Naming the cluster information 
clusters <- data.frame()
for (i in files) {
  if (file.info(i)$size > 0) {
    clusters.tmp <- read.table(i, header = TRUE,
                               sep = "\t",
                               stringsAsFactors = FALSE,
                               fill = TRUE)
    
    filename <- unlist(str_split(i, pattern = "/"))[3]
    
    filename <- unlist(str_split(filename, pattern = "\\."))[1]
    sample <- unlist(str_split(filename, pattern = "_"))
    sample <- sample[1:(length(sample)-1)]
    sample <- paste0(sample, collapse = '_')
    
    clusters.tmp$sample <- sample 
    # Appending to the master table.
    clusters <- clusters %>% bind_rows(clusters.tmp)
  }
}
# Removing the entities I no longer need.
rm(clusters.tmp, i, filename, sample, files)


##################################################
#
# Third part of adding crisr gold data to cluster table:
# This is the part for adding the off-target the the cluster information
#
# make info to add off-target cas9 information to cluster table

clusters$predictoff <- NA
clusters$predictoff.number <- NA
clusters$predictoff.gene <- NA


for (k in seq(nrow(clusters))) {
  target.chrom <- clusters$chrom[k]
  target.start <- clusters$start[k]
  target.end <- clusters$end[k]
  target.sample <- clusters$sample[k]
  predictions.this.chrom <- predictions %>% filter(chrom == target.chrom)
  predictions.this.region <- predictions.this.chrom %>%
    filter((start >= target.start & start <= target.end) |
             (end >= target.start & end <= target.end) |
             (start <= target.start & end >= target.end))
  if(nrow(predictions.this.region) > 0){
    clusters$predictoff[k] <- paste0(predictions.this.region$offtarget.id, collapse = "|")
    clusters$predictoff.number[k] <- paste0(predictions.this.region$name, collapse = "|")
    clusters$predictoff.gene[k] <- paste0(predictions.this.region$gene, collapse = "|")
  }
  rm(predictions.this.region)
}
rm(target.chrom, target.start, target.end, k)

clusters <- clusters[, c(8,2:5,7,9:11)]

##############################################
#
# adding a bit more of the gRNA data and exporting

clusters <- clusters %>% separate(predictoff.gene, into = c("predictoff.gene", NA), sep = "\\|") #removes redundent words

clusters <- clusters %>% #seperating multiple overlaps
  separate(predictoff, into = c("predictoff.1", "predictoff.2"), sep = "\\|")

predictions <- predictions %>% #not sure why we are changing the pridction dataframe
  mutate(offtarget.id = as.character(offtarget.id))

clusters <- clusters %>%
  left_join(predictions, by = c("predictoff.1" = "offtarget.id"))
clusters <- clusters %>%
  left_join(predictions, by = c("predictoff.2" = "offtarget.id"))

clusters <- clusters %>% mutate(length = end - start)



write.table(clusters,
            file = "clustered_reads_global_20210216.txt",
            sep = "\t",
            row.names = FALSE, quote = FALSE,
            dec = ",")


##############################################################################
#
#
#
############################################################################












# Overlap between the samples ---------------------------------------------
clusters.1 <- clusters %>% filter(sample == "N707_N505")
clusters.2 <- clusters %>% filter(sample == "N708_N505")

clusters.1$mappings.2 <- NA
clusters.1$mappings.perc.2 <- NA

for(i in seq(nrow(clusters.1))){
  chrom.tmp <- clusters.1$chrom[i]
  start.tmp <- clusters.1$start[i]
  end.tmp <- clusters.1$end[i]
  clusters.2.tmp <- clusters.2 %>%
    filter(chrom == chrom.tmp)
  if(nrow(clusters.2.tmp) > 0){
    clusters.2.tmp <- clusters.2.tmp %>%
      filter(abs(start - start.tmp) <= 10000 |
               abs(end - end.tmp) <= 10000)
    if(nrow(clusters.2.tmp) > 0){
      clusters.1$mappings.2[i] <- sum(clusters.2.tmp$mappings)
      clusters.1$mappings.perc.2[i] <- sum(clusters.2.tmp$mappings.perc)
    }
  }
}

overlap <- clusters.1 %>% filter(!is.na(mappings.2))
write.table(overlap, file = "20210111_overlap.txt", quote = FALSE, sep = "\t",
            row.names = FALSE)

# Clusters overlap --------------------------------------------------------
clusters.rep <- data.frame()


for(i in unique(clusters$sample)){
  clusters.tmp <- clusters %>% filter(sample == i)
  clusters.tmp <- clusters.tmp %>% arrange(chrom, end)

  # Assigning the groups based on the distance. The group ID changes when the
  # neighbour is further away than 1000 bp.
  site.id <- 1
  clusters.tmp$site.id <- NA
  clusters.tmp$site.id[1] <- 1
  if(nrow(clusters.tmp) > 1){
    for (j in seq(2, nrow(clusters.tmp))) {
      if (is.na(clusters.tmp$start[j]) | is.na(clusters.tmp$end[j])) {
        next
      }
      if (clusters.tmp$chrom[j] == clusters.tmp$chrom[j-1]){
        if (min(abs(clusters.tmp$start[j] - clusters.tmp$end[j-1]),
                abs(clusters.tmp$start[j-1] - clusters.tmp$end[j]),
                abs(clusters.tmp$start[j] - clusters.tmp$start[j-1]),
                abs(clusters.tmp$end[j] - clusters.tmp$end[j-1])) <= 1000){

        } else {
          site.id <- site.id + 1
        }
      } else {
        site.id <- site.id + 1
      }
      clusters.tmp$site.id[j] <- site.id
    }
    rm(j)
    print(sprintf("%d sites were identified in this sample.", site.id))
  }

  # Summarising the table based on those sites.
  clusters.tmp.agg <- clusters.tmp %>% group_by(site.id) %>%
    summarise(replicates = n(),
              mappings.in.reps = str_c(mappings, collapse = " | "),
              mappings = sum(mappings),
              mappings.perc.in.reps = str_c(round(mappings.perc,3), collapse = " | "),
              mappings.perc = mean(mappings.perc),
              chrom = nth(chrom, 1),
              start = min(start),
              end = max(end),
              site = nth(site, 1),
              predictoff = nth(predictoff, 1))

  clusters.tmp.agg <- clusters.tmp.agg %>% arrange(desc(mappings))
  clusters.tmp.agg$sample <- i
  clusters.rep <- clusters.rep %>% bind_rows(clusters.tmp.agg)

}

write.table(clusters.rep, file = "aav_local_031120.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,
            dec = ",")
write.table(clusters, file = "aav_local_separate_replicates_031120.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,
            dec = ",")

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



# 28.11.2020
df <- read.table(file = "aav_totals.txt", header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE)

ggplot(df %>% filter(sample == 2))+
  geom_bar(aes(x = 1, y = mappings.perc, fill = site), position = position_stack(), stat = "identity")+
  coord_polar(theta = "y")+
  theme_classic()

ggsave("pie.sample2.svg", device = "svg", width = 5, height = 5)






####################################################################
###################################################################

#           EXTRA CODE
#
#


# Locus annotation --------------------------------------------------------
# hg38.genes <-
#   read.table("hg38_genes_NCBI.txt", header = TRUE,
#              sep = "\t", stringsAsFactors = FALSE)
# hg38.genes <- hg38.genes %>%
#   filter(!str_detect(chrom, "_") & str_detect(name, "NM")) %>%
#   group_by(name2) %>%
#   summarise(chrom = nth(chrom, 1), start = min(txStart),
#             end = max(txEnd)) %>%
#   mutate(gene.name = name2) %>% select(-name2)
-----------------------------------------------------------------------



################ this is for checking overlaps with a rgion of interst, aka target site #########
# # 28.11.2020
# roi.chr <- "chr11"
# roi.start <- 5225000
# roi.end <- 5228500
# 
# # Sample 3
# maps.3 <- maps %>% filter(sample == "4", rep == "1", chrom == roi.chr,
#                           start >= roi.start, end <= roi.end)
# 
# bin.width <- 10
# bins <- seq(roi.start, roi.end, bin.width)
# 
# bins.starts <- c()
# bins.ends <- c()
# bins.count <- c()
# 
# for(i in seq(2, length(bins))){
#   bin.start <- bins[i-1]
#   bin.end <- bins[i]
#   maps.tmp <- maps.3 %>% filter(start >= bin.start & start <= bin.end |
#                                   end >= bin.start & end <= bin.end |
#                                   start <= bin.start & end >= bin.end)
#   bins.starts <- c(bins.starts, bin.start)
#   bins.ends <- c(bins.ends, bin.end)
#   bins.count <- c(bins.count, nrow(maps.tmp))
# }
# binned.maps <- data.frame(start = bins.starts, end = bins.ends,
#                           reads = bins.count)
# rm(bins.starts, bins.ends, bins.count, bin.width)
# ggplot(binned.maps)+
#   geom_line(aes(x = start, y = reads))+
#   theme_classic()+
#   geom_vline(xintercept = c(5226576, 5227080))+
#   geom_vline(xintercept = c(5225901, 5227808), color = "red")
# 
# ggsave("binnedmaps_sample4_1.svg", device= "svg", width = 6, height = 3)
# 
# 
# ggplot(maps.3)+
#   geom_histogram(aes(x = start), bins = 250)+
#   theme_classic()+
#   geom_vline()
# 

###############################################################3

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

#####################################################################################


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


