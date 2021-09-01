# Setting the working directory -------------------------------------------
# For the laptop
setwd("C:/Users/User/Documents/Lab/crispr/guide_seq2/")
# For the office PC
setwd("S:/Misha/12 - Data Analysis/chu_hbb/")
#for Eric
setwd("/home/eric/Data/Spaced_Nicking/AAVSeq2_ELANE_PRF1_IL7R_2_2021/bed_files")
setwd("/home/eric/Data/Spaced_Nicking/AAVSeq2_ELANE_PRF1_IL7R_2_2021/bed_files")

# Libraries ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggrepel)


# Data import -------------------------------------------------------------
#Misha: files <- list.files("bed_files_ELANE_PRF1_IL7R/", pattern = ".bed", full.names = TRUE)
#local
files <- list.files("local", pattern = ".bed", full.names = TRUE)
files[[1]]

maps.tmp <- read.table(files[[1]], header = FALSE,
                       sep = "\t",
                       stringsAsFactors = FALSE,
                       fill = TRUE)
maps.tmp <- maps.tmp %>% mutate(V2 = as.integer(V2))
                                
sapply(maps.tmp, typeof)

maps.tmp <- maps.tmp %>% filter(str_detect(V1, "chr"))

filename <- unlist(str_split(files[[1]], pattern = "/"))[3]
direction <- unlist(str_split(filename, pattern = "_"))[2]
filename <- unlist(str_split(files[[1]], pattern = "/"))

# Merging the mappings from different files into a single dataframe via a loop.
maps <- data.frame()
for (i in files) {
  if (file.info(i)$size > 0) { #this is checking if the files in the list is bigger than 0 bytes
    maps.tmp <- read.table(1, header = FALSE,
                           sep = "\t",
                           stringsAsFactors = FALSE,
                           fill = TRUE)
    maps.tmp <- maps.tmp %>% mutate(V2 = as.integer(V2)) # %>% is piping from the dplyr package. iris %>% head() %>% summary() is equivalent to summary(head(iris))
    # Due to the unequal column number, some rows "wrap" into two rows. To
    # counteract this behaviour, I filter the wrapped tails of those rows.
    maps.tmp <- maps.tmp %>% filter(str_detect(V1, "chr"))
    # Extracting the sample name from the filename.
    filename <- unlist(str_split(i, pattern = "/"))[2]
    direction <- unlist(str_split(filename, pattern = "_"))[1]
    filename <- unlist(str_split(filename, pattern = "_"))[2]
    rep <- last(unlist(str_split(filename, pattern = "-")))
    
    sample <- unlist(str_split(filename, pattern = "-"))[-1]
    sample <- sample[-length(sample)]
    sample <- paste0(sample, collapse = "_")
    # Assigning new columns to store the sample features.
    maps.tmp$sample <- sample
    maps.tmp$rep <- rep
    maps.tmp$direction <- direction
    # Appending to the master table.
    maps <- maps %>% bind_rows(maps.tmp)

  }
}
# Removing the entities I no longer need.
rm(maps.tmp, i, filename, sample, rep, files, direction)



# Update from 11.01.2021 for Eric -----------------------------------------
#df <- read.table("N707_N505global.sorted.bed", header = FALSE,
#                 sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
#df$sample <- "N707_N505"
#df1 <- read.table("N708_N505global.sorted.bed", header = FALSE,
#                  sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
#df1$sample <- "N708_N505"
#maps <- bind_rows(df, df1)


# I need just the first columns for now.
maps <- maps %>% select(V1, V2, V3, V4, V6, sample, rep, direction)
colnames(maps) <- c("chrom", "start", "end", "read.id",
                    "strand", "sample", "rep", "direction")

# Collapsing reads into the contigs.
contigs <- maps %>% group_by(read.id) %>%
  summarise(chrom.check = n_distinct(chrom),
            chrom = nth(chrom,1),
            start = min(start),
            end = max(end),
            sample.check = n_distinct(sample),
            sample = nth(sample, 1),
            rep = nth(rep, 1),
            direction = nth(direction, 1)) %>%
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
maps.long <- maps %>% filter(read.id %in% contigs.long$read.id)
maps.ambig$comment <- "Mapped to diff chroms"
maps.long$comment <- "Mapped to distant loci"
maps.suspicious <- bind_rows(maps.ambig, maps.long)
rm(maps.ambig, maps.long, contigs.ambig.map, contigs.long)
maps.suspicious <- maps.suspicious %>%
  mutate(length = end - start + 1)
write.table(maps.suspicious, "ambiguous.mappings.20210208.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# if the samples seem good
contigs <- contigs %>% filter(length < 2500, !str_detect(chrom, "_"))
contigs <- contigs %>% unite(samplerep, c(sample, rep))
contigs <- contigs %>% unite(samplerepdir, c(samplerep, direction), remove = FALSE)
contigs.agg <- contigs %>% group_by(samplerepdir) %>%
  summarise(paired.read.count = n())
contigs <- contigs %>%
  separate(samplerep, c("sample1", "sample2", "rep"), extra = "merge", fill = "left") %>%
  unite(sample, c(sample1, sample2)) %>% mutate(sample = str_remove(sample, "NA_"))


# Locus annotation --------------------------------------------------------
hg38.genes <-
  read.table("hg38_genes_NCBI.txt", header = TRUE,
             sep = "\t", stringsAsFactors = FALSE)
hg38.genes <- hg38.genes %>%
  filter(!str_detect(chrom, "_") & str_detect(name, "NM")) %>%
  group_by(name2) %>%
  summarise(chrom = nth(chrom, 1), start = min(txStart),
            end = max(txEnd)) %>%
  mutate(gene.name = name2) %>% select(-name2)

samples <- contigs$samplerepdir %>% unique()

for(i in samples){
  contigs.tmp <- contigs %>% filter(samplerepdir == i)
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

  write.table(contigs.tmp.agg, file = paste0(i, "_", format(Sys.time(),"%d%m%y"), ".txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
# I put these files into a folder called "clusters_20210210"

# Overlap with the CRISPR-Gold predictions --------------------------------
predictions <- read.table("ELANE_ILR7_PRF1_offtargets.txt",
                          header= TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
predictions <- predictions %>%
  separate(position, into = c("chrom", "start", "end"), sep = "\\:") %>%
  mutate(start = as.numeric(start), end = as.numeric(end))
predictions <- predictions %>%
  separate(name, into = c("grna", "onoff", "number"), sep = "_")
predictions <- predictions %>%
  select(-grna, -onoff)
predictions$number[is.na(predictions$number)] <- "0"

predictions$offtarget.id <- seq(nrow(predictions))

window <- 1000
predictions <- predictions %>%
  mutate(start.precise = start, end.precise = end,
         middle = (start + end)/2,
         start = floor(middle - window/2), end = start + window) %>%
  select(-middle, -strand)


# Importing the clustered data.
files <- list.files("clusters_20210215/", pattern = ".txt", full.names = TRUE)

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
    sample <- unlist(str_split(filename, pattern = "_"))
    sample <- sample[-length(sample)]
    sample <- paste0(sample, collapse = "_")
    # Assigning new columns to store the sample features.
    clusters.tmp$sample <- sample #; clusters.tmp$rep <- rep
    #clusters.tmp$mode <- mode
    # Appending to the master table.
    clusters <- clusters %>% bind_rows(clusters.tmp)

  }
}
# Removing the entities I no longer need.
rm(clusters.tmp, i, filename, sample, files, rep)


clusters <- clusters %>%
  mutate(sample = str_replace(sample, "_Cas9", "-Cas9"))

clusters <- clusters %>%
  separate(sample, into = c("sample", "rep", "direction"), sep = "_") %>%
  separate(sample, into = c("sample", "enzyme"), sep = "-")
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
    clusters$predictoff.number[k] <- paste0(predictions.this.region$number, collapse = "|")
    clusters$predictoff.gene[k] <- paste0(predictions.this.region$gene, collapse = "|")
  }
  rm(predictions.this.region)
}
rm(target.chrom, target.start, target.end, k)


clusters <- clusters %>%
  separate(predictoff.gene, into = c("predictoff.gene", NA), sep = "\\|")
clusters <- clusters %>%
  separate(predictoff, into = c("predictoff.1", "predictoff.2"), sep = "\\|")
predictions <- predictions %>%
  mutate(offtarget.id = as.character(offtarget.id))
clusters <- clusters %>%
  left_join(predictions, by = c("predictoff.1" = "offtarget.id"))
clusters <- clusters %>%
  left_join(predictions, by = c("predictoff.2" = "offtarget.id"))

clusters <- clusters %>% mutate(length = end - start)

write.table(clusters,
            file = "clustered_reads_20210210.txt",
            sep = "\t",
            row.names = FALSE, quote = FALSE,
            dec = ",")
