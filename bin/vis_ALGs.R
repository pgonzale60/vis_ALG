#!/usr/bin/env Rscript

library(optparse)
library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scales))
library(ggplot2)
library(gtools)



option_list = list(
  make_option(c("-b", "--busco"), type="character", default=NULL, 
              help="busco full_table.tsv file", metavar="file.tsv"),
  make_option(c("-n", "--nigon"), type="character", default="gene2Nigon_busco20200927.tsv.gz", 
              help="busco id assignment to Nigons [default=%default]", metavar="file.tsv"),
  make_option(c("-w", "--windowSize"), type="integer", default=5e5, 
              help="window size to bin the busco genes [default=%default]. Sequences shorter than twice this integer will not be shown in the plot", metavar="integer"),
  make_option(c("-m", "--minimumGenesPerSequence"), type="integer", default=15, 
              help="sequences (contigs/scaffolds) with less than this number of busco genes will not be shown in the plot [default=%default]", metavar="integer"),
  make_option(c("-o", "--outPlot"), type="character", default="Nigons.jpeg", 
              help="output image [default=%default]. Should include one of the following extensions: eps, ps, tex, pdf, jpeg, tiff, png, bmp or svg", metavar="file"),
  make_option(c("--height"), type="integer", default=6, 
              help="height of plot. Increase this value according to the number of ploted sequences [default=%default]", metavar="integer"),
  make_option(c("--width"), type="integer", default=5, 
              help="width of plot [default=%default]", metavar="integer"),
  make_option(c("-s", "--species"), type="character", default="", 
              help="Title to be italicized in the plot [default=%default]", metavar="Genus_species")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Load data
nigonDict <- read_tsv(opt$nigon,
                      col_types = c(col_character(), col_character()))
busco <- suppressWarnings(read_tsv(opt$busco,
                  col_names = c("Busco_id", "Status", "Sequence",
                                "start", "end", "strand", "Score", "Length",
                                "OrthoDB_url", "Description"),
                  col_types = c("ccciicdicc"),
                  comment = "#"))

windwSize <- opt$windowSize
minimumGenesPerSequence <- opt$minimumGenesPerSequence
spName <- opt$species
if(grepl(".", opt$species)) {
  spName <- paste0("*", sub("_", " ", opt$species), "*")
} 

# Specify Nigon colors
cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b")

# Filter data

fbusco <- filter(busco, !Status %in% c("Missing")) %>%
  left_join(nigonDict, by = c("Busco_id" = "Orthogroup")) %>%
  mutate(nigon = ifelse(is.na(nigon), "-", nigon),
         stPos = start) %>%
  filter(nigon != "-")

consUsco <- group_by(fbusco, Sequence) %>%
  mutate(nGenes = n(),
         mxGpos = max(stPos)) %>%
  ungroup() %>%
  filter(nGenes > minimumGenesPerSequence, mxGpos > windwSize * 2)


# Plot

plNigon <- group_by(consUsco, Sequence) %>%
  mutate(ints = as.numeric(as.character(cut(stPos,
                                            breaks = seq(0, max(stPos), windwSize),
                                            labels = seq(windwSize, max(stPos), windwSize)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + windwSize, ints)) %>%
  count(ints, nigon) %>%
  ungroup() %>%
  mutate(scaffold_f = factor(Sequence,
                             levels = mixedsort(unique(Sequence)))) %>%
  ggplot(aes(fill=nigon, y=n, x=ints-windwSize)) + 
  facet_grid(scaffold_f ~ ., switch = "y") +
  geom_bar(position="stack", stat="identity") +
  ggtitle(spName) +
  theme_minimal() +
  scale_y_continuous(breaks = scales::pretty_breaks(4),
                     position = "right") +
  scale_x_continuous(labels = label_number_si()) +
  scale_fill_manual(values = cols) +
  guides(fill = guide_legend(ncol = 1,
                             title = "Nigon")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        panel.border = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        text = element_text(size=9),
        plot.title = ggtext::element_markdown()
        )

ggsave(opt$out, plNigon, width = opt$width, height = opt$height)

