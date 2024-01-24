# Set working directory
# setwd("/Users/tiagopereira/Dropbox/SK-data/16S-soil-single-nematode/R-codes")

# Install qiime2R
# download.file("https://github.com/jbisanz/qiime2R/archive/master.zip", "source.zip")
# unzip("source.zip")
# install.packages("qiime2R-master", repos = NULL, type="source")

# ggh4x
# install.packages("remotes")
# remotes::install_github("teunbrand/ggh4x")

# library("BiocManager")
# BiocManager::install("microbiome")

# install.packages("remotes")
# remotes::install_github("vmikk/metagMisc")

#Load libraries
library("qiime2R") #Integrating QIIME2 and R for data visualization and analysis using qiime2R.
library("ape") #Analyses of Phylogenetics and Evolution.
library("Biostrings") #Manipulation of biological strings.
library("biomformat") #An interface package for the BIOM file format.
library("phyloseq") #Handling and analysis of high-throughput microbiome census data.
library("Hmisc") #Harrell Miscellaneous. Many functions useful for data analysis, high-level graphics.
library("yaml") #Methods to Convert R Data to YAML and Back.
library("tidyr") #Provides tools to work with tables, dealing with columns and rows as well as individual cells. 
library("stats") #R statistical functions.
library("utils") #R utility functions.
library("ggplot2") #Used for different graphics/plots styles. Very flexible to manage data.
library("readxl") #Used to import data matrices in excel format, including the different sheets.
library("ggthemes") #Provides additional flexibility/features for ggplot2.
library("extrafont") #Allows to use fonts other than basic PostScript fonts.
library("metagMisc") #Miscellaneous functions for metagenomic analysis.
library("vegan") #Community ecology package; it can be use for diverse multivariate analysis.
library("remotes") #R Package Installation from Remote Repositories, Including 'GitHub'
library("mctoolsr") #Microbial community analyses.
library("RColorBrewer") #Provides color schemes for maps/graphics designed by Cynthia Brewer.
library("circlize") #Visualization of huge amounts of information.
library("viridis") #Color Palettes for R.
library("decontam") #Simple statistical identification of contaminated sequence features.
library("plyr") #Provides a set of tools/functions to manipulate datasets/matrices.
library("ALDEx2") #Differential abundance analysis package.
library("tibble") #Simple Data Frames.
library("gplots") #Various R Programming Tools for Plotting Data.
library("dplyr") #A tool for working/manipulating dataframes.
library("ComplexHeatmap") # Package for Complex heatmaps reveal patterns and correlations in multidimensional genomic data.
library("ggpubr") #Provides some easy-to-use functions to customize ggplot2 graphics.
library("ggh4x") #ggh4x package is a ggplot2 extension package.
library("microbiome") #Tools for microbiome analysis.
library("data.table") #Extension of 'data.frame'.
library("scales") #Use for graphical scales map data.
library("gtable") #Provides tools/functions to make easier to work with tables.
library("grid") #Rewrite graphics layout capabilities, also supporting interactions.
library("cowplot") #Simple add-on to ggplot2; provides additional features to improve graphic quality.
library("stringr") #Provides a set of functions to work with strings.
library("ggord") #A package for creating ordination plots with ggplot2.
library("ggfortify") #Unified plotting tools for common statistics and graphics through ggplot2.
library("ggrepel") #Provides geoms for ggplot2 to repel overlapping text labels.
library("factoextra") #Provides functions for multivariate data analysis as well as graphics using ggplot2.
library("ggbiplot") #Implements the biplot (i.e. samples data points and vectors of variables) using ggplot2.
library("devtools") #Allows the installation of R packages.
library("MicEco") #Functions for microbiome analysis based on phyloseq objetcs

# Set plotting theme
theme_set(theme_bw())

# Load sources
source(file = "src/miseqR.R")
source(file = "src/simper.R")
source(file = "src/summary.R")

# Import phyloseq filtered no scale for feeding group analysis 
phy_16_prune_nem_noncontam_prev05_true_all_filt <- readRDS("data/phy_16_prune_nem_noncontam_prev05_true_all_filt.RDS")
class(phy_16_prune_nem_noncontam_prev05_true_all_filt)
smin <- min(sample_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt)) # smin = 2, i.e., minimum number of reads per sample
phy_16_prune_nem_noncontam_prev05_true_all_filt # 3057 taxa and a total of 530 samples.

# Import phyloseq scale
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale <- readRDS("data/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale.RDS")
class(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale)
smin <- min(sample_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale)) # smin = 114, i.e., minimum number of reads per sample
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale # 3055 taxa and a total of 520 samples.

# Aldex on all nematode taxa (i.e., different infraorders) from phyloseq object
# Collapse taxonomic ranks
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank1 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, taxrank = "Phylum") 
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank2 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, taxrank = "Class")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank3 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, taxrank = "Order")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank4 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, taxrank = "Family")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank5 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, taxrank = "Genus")
#################### Aldex2 Non Selected Nematode Groups Phyloseq objects ####################
##### Across Sampling sites
# ASV00011, ASV00004, ASV00001, ASV00005, ASV00006, ASV00191
# Verrucomicrobiota, Firmicutes, Bacteroidota, Proteobacteria, Actinobacteriota, Gemmatimonadota
all_nematodes_rank1_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank1)),
                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank1)$Sample.Site,
                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_rank1_site, "exported_tables/nema/all_nematodes_rank1_site.csv")

# ASV00011, ASV00018, ASV00004, ASV00001, ASV00006, ASV00005
# Verrucomicrobiae, Gammaproteobacteria, Bacilli, Bacteroidia, Actinobacteria, Alphaproteobacteria
all_nematodes_rank2_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank2)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank2)$Sample.Site,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_rank2_site, "exported_tables/nema/all_nematodes_rank2_site.csv")

# ASV00011, ASV00004, ASV00018, ASV00006, ASV00001
# ASV00027, ASV00007, ASV00023, ASV00034, ASV00033
# Opitutales, Staphylococcales, Burkholderiales, Corynebacteriales, Cytophagales
# Pseudomonadales, Sphingomonadales, Xanthomonadales, Caulobacterales, Rhodobacterales
all_nematodes_rank3_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank3)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank3)$Sample.Site,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_rank3_site, "exported_tables/nema/all_nematodes_rank3_site.csv")

# ASV00018, ASV00011, ASV00006, ASV00004, ASV00007, ASV00001
# ASV00082, ASV00033, ASV00034, ASV00027, ASV00023, ASV00025
# Comamonadaceae, Opitutaceae, Mycobacteriaceae, Staphylococcaceae, Sphingomonadaceae, Amoebophilaceae
# Lactobacillaceae, Caulobacteraceae, Rhodobacteraceae, Moraxellaceae, Enterococcaceae, Burkholderiaceae
all_nematodes_rank4_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank4)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank4)$Sample.Site,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_rank4_site, "exported_tables/nema/all_nematodes_rank4_site.csv")


# ASV00011, ASV00006, ASV00004, ASV00018, ASV00001
# ASV00007, ASV00024, ASV00082, ASV00033, ASV00023 
# Lacunisphaera, Mycobacterium, Staphylococcus, Rhizobacter, Candidatus Cardinium
# Novosphingobium, Limnohabitans, Leuconostoc, Enhydrobacter, Tetragenococcus
all_nematodes_rank5_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank5)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank5)$Sample.Site,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_rank5_site, "exported_tables/nema/all_nematodes_rank5_site.csv")

# ASV00011, ASV00006, ASV00018, ASV00004, ASV00001
# ASV00024, ASV00082, ASV00033, ASV00023, ASV00007
# Lacunisphaera, Mycobacterium, Rhizobacter, Staphylococcus, Candidatus Cardinium
# Limnohabitans, Leuconostoc, Novosphingobium
all_nematodes_asv_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale)$Sample.Site,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_asv_site, "exported_tables/nema/all_nematodes_asv_site.csv")

#################### Across Habitats ####################
# ASV00001, ASV00006
# Bacteroidota, Actinobacteriota
all_nematodes_rank1_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank1)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank1)$Habitat,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_rank1_habitat, "exported_tables/nema/all_nematodes_rank1_habitat.csv")

# ASV00001, ASV00018, ASV00006, ASV00251, ASV00004
# Bacteroidia, Gammaproteobacteria, Actinobacteria, Thermoleophilia, Bacilli
all_nematodes_rank2_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank2)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank2)$Habitat,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_rank2_habitat, "exported_tables/nema/all_nematodes_rank2_habitat.csv")

# ASV00004, ASV00018, ASV00001, ASV00023
# ASV00006, ASV00011, ASV00027, ASV00010
# Staphylococcales, Burkholderiales, Cytophagales, Lactobacillales
# Corynebacteriales, Opitutales, Pseudomonadales, Flavobacteriales
all_nematodes_rank3_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank3)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank3)$Habitat,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_rank3_habitat, "exported_tables/nema/all_nematodes_rank3_habitat.csv")

# ASV00011, ASV00004, ASV00018, ASV00006, ASV00001
# Opitutaceae, Staphylococcaceae, Comamonadaceae, Mycobacteriaceae, Amoebophilaceae
all_nematodes_rank4_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank4)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank4)$Habitat,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_rank4_habitat, "exported_tables/nema/all_nematodes_rank4_habitat.csv")

# ASV00011, ASV00006, ASV00018, ASV00004, ASV00001, ASV00007
# Lacunisphaera, Mycobacterium, Rhizobacter, Staphylococcus, Candidatus Cardinium, Novosphingobium
all_nematodes_rank5_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank5)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank5)$Habitat,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_rank5_habitat, "exported_tables/nema/all_nematodes_rank5_habitat.csv")


all_nematodes_asv_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale)),
                                        phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale)$Habitat,
                                        mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_asv_habitat, "exported_tables/nema/all_nematodes_asv_habitat.csv")
######################### Across Soil type ######################### 
# ASV00001, ASV00006
# Verrucomicrobiota, Actinobacteriota
all_nematodes_rank1_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank1)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank1)$Soil.Type,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_rank1_soil_type, "exported_tables/nema/all_nematodes_rank1_soil_type.csv")

# ASV00018, ASV00001, ASV00006
# Gammaproteobacteria, Bacteroidia, Actinobacteria
all_nematodes_rank2_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank2)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank2)$Soil.Type,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_rank2_soil_type, "exported_tables/nema/all_nematodes_rank2_soil_type.csv")

# ASV00018, ASV00010, ASV00006, ASV00034
# Burkholderiales, Flavobacteriales, Corynebacteriales, Rhodobacterales
all_nematodes_rank3_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank3)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank3)$Soil.Type,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_rank3_soil_type, "exported_tables/nema/all_nematodes_rank3_soil_type.csv")

# ASV00018, ASV00006, ASV00010, ASV00034
# Comamonadaceae, Mycobacteriaceae, Weeksellaceae, Rhodobacteraceae
all_nematodes_rank4_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank4)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank4)$Soil.Type,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_rank4_soil_type, "exported_tables/nema/all_nematodes_rank4_soil_type.csv")

# ASV00006, ASV00018, ASV00010, ASV00007
# Mycobacterium, Rhizobacter, Cloacibacterium, Novosphingobium
all_nematodes_rank5_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank5)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_no_sel_group_rank5)$Soil.Type,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_rank5_soil_type, "exported_tables/nema/all_nematodes_rank5_soil_type.csv")

all_nematodes_asv_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale)$Soil.Type,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(all_nematodes_asv_soil_type, "exported_tables/nema/all_nematodes_asv_soil_type.csv")

####################### Selected Nematode Groups Phyloseq Object #######################
# Exclude nematode groups with low representation (e.g., Monhysterida, Trefusiina, and unknown groups)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale)$Suborder)
# Total of six nematode samples
# Nem.4, Nem.274: Phylum(Nematoda),	Class(Chromadorea),	Subclass(Chromadoria), Order(Monhysterida), Suborder(Monhysterina)
# Nem.212: Phylum(Nematoda), Class(Enoplea),	Subclass(Enoplia), Order(Enoplida),	Suborder(Trefusiina)
# Nem.281, Nem.282 Nem.401:	Phylum(Nematoda),	Class(unknown),	Subclass(unknown), Order(unknown), Suborder(unknown)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups <- 
  subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale,
                 Class != "unknown" &
                   Order != "Monhysterida" &
                   Suborder != "Trefusiina") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups)$Suborder)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups # 3022 taxa and 515 samples
# These represent the major nematode groups in our dataset
# Aphelenchids, Cephalobids, Dorylaimids, Plectids, Rhabditids, Tylenchids

# Save sample data of selected groups
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_meta <- sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_meta <- data.frame(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_meta)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_asv <- as.data.frame(otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups, "matrix"))
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_tax <- as.data.frame(tax_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups))

write.csv(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_meta, "exported_tables/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_meta.csv")
write.csv(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_asv, "exported_tables/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_asv.csv")
write.csv(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_tax, "exported_tables/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_tax.csv")

# Import sample data with nematode group information
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_meta_nema_group <- read.csv("exported_tables/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_meta_new.csv", row.names = "Sample")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_meta_nema_group_data <- sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_meta_nema_group)

# Merge new sample_data with phyloseq object of selected nematode groups
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new <- merge_phyloseq(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups, phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_meta_nema_group_data)
head(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new))
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new)$Nema.Group)

############### alpha-diversity among nematode groups
alpha_div_16S_nema_group <- data.frame(
  "Observed" = phyloseq::estimate_richness(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new, measures = "Shannon"),
  "InvSimpson" = phyloseq::estimate_richness(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new, measures = "InvSimpson"),
  "Reads" = phyloseq::sample_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new),
  "Soil.Type" = phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new)$Soil.Type,
  "Habitat" = phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new)$Habitat,
  "Sample.Site" = phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new)$Sample.Site,
  "Nema.Group" = phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new)$Nema.Group)
alpha_div_16S_nema_group$Evenness <- alpha_div_16S_nema_group$Shannon/log(alpha_div_16S_nema_group$Observed)
head(alpha_div_16S_nema_group)

# Rename variable InvSimpson to Simpson
# The function rename & %>% works on dplyr. make sure it is loaded.
alpha_div_16S_nema_group <- alpha_div_16S_nema_group %>%
  rename(Simpson = InvSimpson)
head(alpha_div_16S_nema_group)

# Reorder dataframe, first categorical then numerical variables  
alpha_div_16S_nema_group <- alpha_div_16S_nema_group[, c(7, 6, 5, 8, 4, 1, 2, 3, 9)]
head(alpha_div_16S_nema_group)
write.csv(alpha_div_16S_nema_group, "exported_tables/nema/alpha_div_16S_nema_group.csv")

#Summarize alpha diversity measures by nema group
summary_alpha_div_16S_nema_group <- alpha_div_16S_nema_group %>%
  group_by(Nema.Group) %>%
  summarise(
    count = n(),
    mean_reads = mean(Reads),
    sd_reads = sd(Reads),
    mean_observed = mean(Observed),
    sd_observed = sd(Observed),
    mean_shannon = mean(Shannon),
    sd_shannon = sd(Shannon),
    mean_Simpson = mean(Simpson),
    sd_Simpson = sd(Simpson),
    mean_evenness = mean(Evenness),
    sd_evenness = sd(Evenness))
write.csv(summary_alpha_div_16S_nema_group, "exported_tables/nema/summary_alpha_div_16S_nema_group.csv")

# KW analysis alpha_div_16S_nema_group on all metrics at once
# Remember, numerical variables are from columns 4-8. The test is by soil type, column 3
kw_alpha_div_16S_nema_group_univ <- as.data.frame(sapply(5:9, function(x) kruskal.test(alpha_div_16S_nema_group[,x],
                                                                                  alpha_div_16S_nema_group[,4])))


# Rename columns with the proper variable names
kw_alpha_div_16S_nema_group_univ <- kw_alpha_div_16S_nema_group_univ %>%
  rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)

# Make some adjustments on the table
kw_alpha_div_16S_nema_group_univ <- t(kw_alpha_div_16S_nema_group_univ) # transpose
kw_alpha_div_16S_nema_group_univ <- as_tibble(kw_alpha_div_16S_nema_group_univ, rownames = "Metric") # adding rownames as a column
kw_alpha_div_16S_nema_group_univ <- kw_alpha_div_16S_nema_group_univ[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_div_16S_nema_group_univ) # checking object class

# Save resulting table with fwrite to avoid any issues with characters
data.table::fwrite(kw_alpha_div_16S_nema_group_univ, "exported_tables/nema/kw_alpha_div_16S_nema_group_univ.csv")

# pairwise comparison including p-value adjustment
kw_alpha_div_16S_nema_group_univ_pwt <- as.data.frame(sapply(5:9, function(x) pairwise.t.test(alpha_div_16S_nema_group[,x],
                                                                                              alpha_div_16S_nema_group[,4], p.adjust.method = "BH")))

kw_alpha_div_16S_nema_group_univ_pwt <- kw_alpha_div_16S_nema_group_univ_pwt %>%
  rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)


# Make some adjustments on the table
kw_alpha_div_16S_nema_group_univ_pwt <- t(kw_alpha_div_16S_nema_group_univ_pwt) # transpose
kw_alpha_div_16S_nema_group_univ_pwt <- as_tibble(kw_alpha_div_16S_nema_group_univ_pwt, rownames = "Metric") # adding rownames as a column
class(kw_alpha_div_16S_nema_group_univ_pwt) # checking object class

# Save resulting table with fwrite to avoid and issues with characters
fwrite(kw_alpha_div_16S_nema_group_univ_pwt, "exported_tables/nema/kw_alpha_div_16S_nema_group_univ_pwt.csv")

##### alpha diversity among feeding groups
# update sample data to fix Aporcelaimellus to Predator feeding group
# Import sample data with nematode group information
phy_16_prune_nem_noncontam_prev05_true_all_filt_meta <- read.csv("exported_tables/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_meta_new.csv", row.names = "Sample")
phy_16_prune_nem_noncontam_prev05_true_all_filt_meta_data <- sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_meta)

sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt) <- phy_16_prune_nem_noncontam_prev05_true_all_filt_meta_data
head(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt), n=30)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt)$Feeding.Group)

# Remove unknown feeding group
phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt, Feeding.Group != "unknown") %>%
prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new)$Feeding.Group)
phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new # 527 samples by 38 sample variables

alpha_div_16S_feeding_group <- data.frame(
  "Observed" = phyloseq::estimate_richness(phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new, measures = "Shannon"),
  "InvSimpson" = phyloseq::estimate_richness(phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new, measures = "InvSimpson"),
  "Reads" = phyloseq::sample_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new),
  "Soil.Type" = phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new)$Soil.Type,
  "Habitat" = phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new)$Habitat,
  "Sample.Site" = phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new)$Sample.Site,
  "Feeding.Group" = phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new)$Feeding.Group)
alpha_div_16S_feeding_group$Evenness <- alpha_div_16S_feeding_group$Shannon/log(alpha_div_16S_feeding_group$Observed)
head(alpha_div_16S_feeding_group)

# Rename variable InvSimpson to Simpson
# The function rename & %>% works on dplyr. make sure it is loaded.
alpha_div_16S_feeding_group <- alpha_div_16S_feeding_group %>%
  rename(Simpson = InvSimpson)
head(alpha_div_16S_feeding_group)

# Reorder dataframe, first categorical then numerical variables  
alpha_div_16S_feeding_group <- alpha_div_16S_feeding_group[, c(7, 6, 5, 8, 4, 1, 2, 3, 9)]
head(alpha_div_16S_feeding_group)
write.csv(alpha_div_16S_feeding_group, "exported_tables/nema/alpha_div_16S_feeding_group.csv")

#Summarize alpha diversity measures by nema group
summary_alpha_div_16S_feeding_group <- alpha_div_16S_feeding_group %>%
  group_by(Feeding.Group) %>%
  summarise(
    count = n(),
    mean_reads = mean(Reads),
    sd_reads = sd(Reads),
    mean_observed = mean(Observed),
    sd_observed = sd(Observed),
    mean_shannon = mean(Shannon),
    sd_shannon = sd(Shannon),
    mean_Simpson = mean(Simpson),
    sd_Simpson = sd(Simpson),
    mean_evenness = mean(Evenness),
    sd_evenness = sd(Evenness))
write.csv(summary_alpha_div_16S_feeding_group, "exported_tables/nema/summary_alpha_div_16S_feeding_group.csv")

# KW analysis alpha_div_16S_feeding_group on all metrics at once
# Remember, numerical variables are from columns 4-8. The test is by soil type, column 3
kw_alpha_div_16S_feeding_group_univ <- as.data.frame(sapply(5:9, function(x) kruskal.test(alpha_div_16S_feeding_group[,x],
                                                                                       alpha_div_16S_feeding_group[,4])))


# Rename columns with the proper variable names
kw_alpha_div_16S_feeding_group_univ <- kw_alpha_div_16S_feeding_group_univ %>%
  rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)

# Make some adjustments on the table
kw_alpha_div_16S_feeding_group_univ <- t(kw_alpha_div_16S_feeding_group_univ) # transpose
kw_alpha_div_16S_feeding_group_univ <- as_tibble(kw_alpha_div_16S_feeding_group_univ, rownames = "Metric") # adding rownames as a column
kw_alpha_div_16S_feeding_group_univ <- kw_alpha_div_16S_feeding_group_univ[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_div_16S_feeding_group_univ) # checking object class

# Save resulting table with fwrite to avoid any issues with characters
data.table::fwrite(kw_alpha_div_16S_feeding_group_univ, "exported_tables/nema/kw_alpha_div_16S_feeding_group_univ.csv")

# pairwise comparison including p-value adjustment
kw_alpha_div_16S_feeding_group_univ_pwt <- as.data.frame(sapply(5:9, function(x) pairwise.t.test(alpha_div_16S_feeding_group[,x],
                                                                                              alpha_div_16S_feeding_group[,4], p.adjust.method = "BH")))

kw_alpha_div_16S_feeding_group_univ_pwt <- kw_alpha_div_16S_feeding_group_univ_pwt %>%
  rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)


# Make some adjustments on the table
kw_alpha_div_16S_feeding_group_univ_pwt <- t(kw_alpha_div_16S_feeding_group_univ_pwt) # transpose
kw_alpha_div_16S_feeding_group_univ_pwt <- as_tibble(kw_alpha_div_16S_feeding_group_univ_pwt, rownames = "Metric") # adding rownames as a column
class(kw_alpha_div_16S_feeding_group_univ_pwt) # checking object class

# Save resulting table with fwrite to avoid and issues with characters
fwrite(kw_alpha_div_16S_feeding_group_univ_pwt, "exported_tables/nema/kw_alpha_div_16S_feeding_group_univ_pwt.csv")

# Plot alpha diversity measures
# Change names
soil.labs <- c("Compacted", "Uncompacted")
names(soil.labs) <- c("Compact.Soil", "Not.Compact.Soil")

site.labs <- c("01", "02", "03", "04", "05", "06",
               "07", "08", "09", "10", "11", "12",
               "13", "14", "15", "16", "17", "18")
names(site.labs) <- c("SK.01", "SK.02", "SK.03", "SK.04", "SK.05", "SK.06",
                      "SK.07", "SK.08", "SK.09", "SK.10", "SK.11", "SK.12",
                      "SK.13", "SK.14", "SK.15", "SK.16", "SK.17", "SK.18")

habitat.labs <- c("CHA", "CSS", "NGR", "HLC",  "OWL", "RIP")
names(habitat.labs) <- c("Chaparral", "Coastal.scrub", "Native.grass",
                         "Hollyleaf.cherry", "Oak.wood", "Riparian")

alpha_color_habitat <- c("#52392F", "#83643E","#397A4C", "#C29B6C", "#77C063", "#BEDB92")

alpha_color_nema_group <- c("#46f0f0", "#1F78B4", "#aaffc3", "#000075", "#e6194b", "#bcf60c", "#33A02C")
                           
# Plot alpha-diversity per site within habitat/soil type
alpha_div_16S_nema_group %>%
  gather(key = metric, value = value, c("Reads", "Observed", "Shannon", "Simpson", "Evenness")) %>%
  mutate(metric = factor(metric, levels = c("Reads", "Observed", "Shannon", "Simpson", "Evenness"))) %>%
  ggplot(aes(x = Nema.Group, y = value, fill = Nema.Group)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  geom_jitter(aes(color = Nema.Group), height = 0, width = .2) +
  facet_nested(metric ~ ., scales = "free", labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  scale_fill_manual(values = alpha_color_nema_group) +
  scale_color_manual(values = alpha_color_nema_group) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.1, 0.1)) + # plot as % and removes the internal margins
  scale_x_discrete(
                   expand = c(0.1, 0.1),
                   drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("") + # add the title on y axis
  xlab("Nematode Group") + # add the title on x axis
  theme(legend.position="none")
ggsave("results/nema/16S_nema_alpha_diversity_soil_type_nema_group.pdf", width = 8, height = 6, dpi = 150) # save graphic

# Plot alpha diversity measures
# Change names
soil.labs <- c("Compacted", "Uncompacted")
names(soil.labs) <- c("Compact.Soil", "Not.Compact.Soil")

site.labs <- c("01", "02", "03", "04", "05", "06",
               "07", "08", "09", "10", "11", "12",
               "13", "14", "15", "16", "17", "18")
names(site.labs) <- c("SK.01", "SK.02", "SK.03", "SK.04", "SK.05", "SK.06",
                      "SK.07", "SK.08", "SK.09", "SK.10", "SK.11", "SK.12",
                      "SK.13", "SK.14", "SK.15", "SK.16", "SK.17", "SK.18")

habitat.labs <- c("CHA", "CSS", "NGR", "HLC",  "OWL", "RIP")
names(habitat.labs) <- c("Chaparral", "Coastal.scrub", "Native.grass",
                         "Hollyleaf.cherry", "Oak.wood", "Riparian")

alpha_color_habitat <- c("#52392F", "#83643E","#397A4C", "#C29B6C", "#77C063", "#BEDB92")

alpha_color_nema_group <- c("#46f0f0", "#1F78B4", "#aaffc3", "#000075", "#e6194b", "#bcf60c", "#33A02C")

# Plot alpha-diversity per site within habitat/soil type
alpha_div_16S_feeding_group %>%
  gather(key = metric, value = value, c("Reads", "Observed", "Shannon", "Simpson", "Evenness")) %>%
  mutate(metric = factor(metric, levels = c("Reads", "Observed", "Shannon", "Simpson", "Evenness"))) %>%
  ggplot(aes(x = Feeding.Group, y = value, fill = Feeding.Group)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  geom_jitter(aes(color = Feeding.Group), height = 0, width = .2) +
  facet_nested(metric ~ Soil.Type, scales = "free", labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  scale_fill_manual(values = alpha_color_nema_group) +
  scale_color_manual(values = alpha_color_nema_group) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.1, 0.1)) + # plot as % and removes the internal margins
  scale_x_discrete(
    expand = c(0.1, 0.1),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("") + # add the title on y axis
  xlab("Nematode Group") + # add the title on x axis
  theme(legend.position="none")
ggsave("results/nema/16S_nema_alpha_diversity_soil_type_feeding_group.pdf", width = 8, height = 6, dpi = 150) # save graphic

# Aldex on all nematode taxa from selected nematode groups phyloseq object
# Collapse taxonomic ranks
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank1 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new, taxrank = "Phylum")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank2 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new, taxrank = "Class")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank3 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new, taxrank = "Order")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank4 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new, taxrank = "Family")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank5 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new, taxrank = "Genus")

phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new_rank1 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new, taxrank = "Phylum")
phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new_rank2 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new, taxrank = "Class")
phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new_rank3 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new, taxrank = "Order")
phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new_rank4 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new, taxrank = "Family")
phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new_rank5 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_feeding_group_new, taxrank = "Genus")


phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_rank5_meta <- data.frame(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank5))
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_rank5_asv <- as.data.frame(otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank5, "matrix"))
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_rank5_tax <- as.data.frame(tax_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank5))

write.csv(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_rank5_meta, "exported_tables/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_rank5_meta.csv")
write.csv(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_rank5_asv, "exported_tables/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_rank5_asv.csv")
write.csv(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_rank5_tax, "exported_tables/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_rank5_tax.csv")
######################################## Aldex2 on phyoseq objects######################################## 
#################### Across nematode groups and taxonomic ranks #################### 
# ASV00001, ASV00005, ASV00006
# Bacteroidota, Proteobacteria, Actinobacteriota
nematode_nema_group_rank1 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank1)),
                                      phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank1)$Nema.Group,
                                      mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_nema_group_rank1, "exported_tables/nema/nematode_Nema.Group_rank1.csv")

# ASV00001, ASV00005, ASV00006, ASV00018
# Bacteroidia, Alphaproteobacteria, Actinobacteria, Gammaproteobacteria
nematode_nema_group_rank2 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank2)),
                                      phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank2)$Nema.Group,
                                      mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_nema_group_rank2, "exported_tables/nema/nematode_Nema.Group_rank2.csv")

# ASV00006, ASV00007, ASV00018
# Corynebacteriales, Sphingomonadales, Burkholderiales
nematode_nema_group_rank3 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank3)),
                                      phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank3)$Nema.Group,
                                      mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_nema_group_rank3, "exported_tables/nema/nematode_Nema.Group_rank3.csv")

# ASV00006, ASV00007, ASV00018
# Mycobacteriaceae, Sphingomonadaceae, Comamonadaceae
nematode_nema_group_rank4 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank4)),
                                      phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank4)$Nema.Group,
                                      mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_nema_group_rank4, "exported_tables/nema/nematode_Nema.Group_rank4.csv")

# ASV00006, ASV00018
# Mycobacterium, Rhizobacter
nematode_nema_group_rank5 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank5)),
                                      phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank5)$Nema.Group,
                                      mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_nema_group_rank5, "exported_tables/nema/nematode_Nema.Group_rank5.csv")

# ASV00006, ASV00018
nematode_nema_group_asv <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new)$Nema.Group,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_nema_group_asv, "exported_tables/nema/nematode_Nema.Group_asv.csv")

######################### Across sampling sites and taxonomic ranks on selected nematode groups ######################### 
# ASV00011, ASV00004, ASV00001, ASV00005, ASV00006, ASV00191
# Verrucomicrobiota, Firmicutes, Bacteroidota, Proteobacteria, Proteobacteria, Gemmatimonadota
nematode_site_rank1 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank1)),
                                        phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank1)$Sample.Site,
                                        mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_site_rank1, "exported_tables/nema/nematode_site_rank1.csv")

# ASV00011, ASV00018, ASV00004, ASV00001, ASV00006, ASV00005, ASV00251
# Verrucomicrobiae, Gammaproteobacteria, Bacilli, Bacteroidia, Actinobacteria, Alphaproteobacteria, Actinobacteriota
nematode_site_rank2 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank2)),
                                        phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank2)$Sample.Site,
                                        mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_site_rank2, "exported_tables/nema/nematode_site_rank2.csv")

# ASV00011, ASV00018, ASV00004, ASV00006, ASV00027, ASV00001, ASV00007, ASV00023, ASV00033, ASV00034
# Opitutales, Burkholderiales, Staphylococcales, Corynebacteriales, Pseudomonadales, Cytophagales, Sphingomonadales, Lactobacillales, Caulobacterales, Rhodobacterales
nematode_site_rank3 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank3)),
                                        phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank3)$Sample.Site,
                                        mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_site_rank3, "exported_tables/nema/nematode_site_rank3.csv")

# ASV00011, ASV00018, ASV00006, ASV00004, ASV00001, ASV00007, ASV00033, ASV00082, ASV00023, ASV00034, ASV00027, ASV00025
# Opitutaceae, Comamonadaceae, Mycobacteriaceae, Staphylococcaceae, Amoebophilaceae, Sphingomonadaceae, Caulobacteraceae, Lactobacillaceae, Enterococcaceae, Rhodobacteraceae, Moraxellaceae, Burkholderiaceae
nematode_site_rank4 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank4)),
                                        phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank4)$Sample.Site,
                                        mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_site_rank4, "exported_tables/nema/nematode_site_rank4.csv")

# ASV00011, ASV00006, ASV00018, ASV00004, ASV00001, ASV00024, ASV00007, ASV00082, ASV00033, ASV00023
# Lacunisphaera, Mycobacterium, Rhizobacter, Staphylococcus, Candidatus Cardinium, Limnohabitans, Novosphingobium, Leuconostoc, Enhydrobacter Tetragenococcus
nematode_site_rank5 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank5)),
                                        phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank5)$Sample.Site,
                                        mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_site_rank5, "exported_tables/nema/nematode_site_rank5.csv")

# ASV00011, ASV00006, ASV00018, ASV00004, ASV00001, ASV00024, ASV00007, ASV00082, ASV00023
nematode_site_asv <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new)),
                                     phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new)$Sample.Site,
                                     mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_site_asv, "exported_tables/nema/nematode_site_asv.csv")

######################### Across habitats and taxonomic ranks ######################### 
# ASV00001, ASV00005, ASV00006, - Bacteroidota, Proteobacteria, Actinobacteriota
nematode_habitat_rank1 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank1)),
                                        phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank1)$Habitat,
                                        mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_habitat_rank1, "exported_tables/nema/nematode_habitat_rank1.csv")

# ASV00001, ASV00004, ASV00006, ASV00018, ASV00251 - Bacteroidia, Bacilli, Actinobacteria, Gammaproteobacteria, Thermoleophilia
nematode_habitat_rank2 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank2)),
                                        phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank2)$Habitat,
                                        mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_habitat_rank2, "exported_tables/nema/nematode_habitat_rank2.csv")

# ASV00001, ASV00004, ASV00006, ASV00011, ASV00018, ASV00023, ASV00027
# Cytophagales, Staphylococcales, Corynebacteriales, Opitutales, Burkholderiales, Lactobacillales, Bacillales
nematode_habitat_rank3 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank3)),
                                        phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank3)$Habitat,
                                        mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_habitat_rank3, "exported_tables/nema/nematode_habitat_rank3.csv")

# ASV00001, ASV00004, ASV00006, ASV00011, ASV00018
# Amoebophilaceae, Staphylococcaceae, Mycobacteriaceae, Opitutaceae, Comamonadaceae
nematode_habitat_rank4 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank4)),
                                        phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank4)$Habitat,
                                        mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_habitat_rank4, "exported_tables/nema/nematode_habitat_rank4.csv")

# ASV00001, ASV00004, ASV00006, ASV00011, ASV00018
# Candidatus Cardinium, Staphylococcus, Mycobacterium, Lacunisphaera, Rhizobacter
nematode_habitat_rank5 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank5)),
                                        phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank5)$Habitat,
                                        mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_habitat_rank5, "exported_tables/nema/nematode_habitat_rank5.csv")

nematode_habitat_asv <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new)),
                                        phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new)$Habitat,
                                        mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_habitat_asv, "exported_tables/nema/nematode_habitat_asv.csv")

############################## Across soil type and taxonomic ranks ############################## 
# ASV00001, ASV00006
# Bacteroidota, Actinobacteriota
nematode_soil_type_rank1 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank1)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank1)$Soil.Type,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_soil_type_rank1, "exported_tables/nema/nematode_soil_type_rank1.csv")

# ASV00001, ASV00006, ASV00018
# Bacteroidia, Actinobacteria, Gammaproteobacteria
nematode_soil_type_rank2 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank2)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank2)$Soil.Type,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_soil_type_rank2, "exported_tables/nema/nematode_soil_type_rank2.csv")

# ASV00006, ASV00010, ASV00018, ASV00034
# Corynebacteriales, Flavobacteriales, Burkholderiales, Rhodobacterales
nematode_soil_type_rank3 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank3)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank3)$Soil.Type,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_soil_type_rank3, "exported_tables/nema/nematode_soil_type_rank3.csv")

# ASV00006, ASV00010, ASV00018, ASV00034, ASV00007
# Sphingomonadaceae, Weeksellaceae, Comamonadaceae, Rhodobacteraceae, Opitutaceae
nematode_soil_type_rank4 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank4)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank4)$Soil.Type,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_soil_type_rank4, "exported_tables/nema/nematode_soil_type_rank4.csv")

# ASV00006, ASV00007, ASV00010, ASV00018
# Mycobacterium, Novosphingobium, Cloacibacterium, Rhizobacter
nematode_soil_type_rank5 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank5)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank5)$Soil.Type,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_soil_type_rank5, "exported_tables/nema/nematode_soil_type_rank5.csv")

# ASV00006, ASV00007, ASV00004, ASV00018
nematode_soil_type_asv <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new)),
                                      phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new)$Soil.Type,
                                      mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_soil_type_asv, "exported_tables/nema/nematode_soil_type_asv.csv")

##############################  Across feeding groups and taxononic ranks ##############################
# ASV00004 - Firmicutes
nematode_feeding_group_rank1 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank1)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank1)$Feeding.Group,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_feeding_group_rank1, "exported_tables/nema/nematode_Feeding.Group_rank1.csv")

# None (ASV00004 - Bacilli, p=0.067)
nematode_feeding_group_rank2 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank2)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank2)$Feeding.Group,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_feeding_group_rank2, "exported_tables/nema/nematode_Feeding.Group_rank2.csv")

# ASV00006, ASV00018, ASV00001
# Corynebacteriales, Burkholderiales, Cytophagales (p=0.066)
nematode_feeding_group_rank3 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank3)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank3)$Feeding.Group,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_feeding_group_rank3, "exported_tables/nema/nematode_Feeding.Group_rank3.csv")

# ASV00006, ASV00018, ASV00001
# Mycobacteriaceae, Comamonadaceae, Amoebophilaceae
nematode_feeding_group_rank4 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank4)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank4)$Feeding.Group,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_feeding_group_rank4, "exported_tables/nema/nematode_Feeding.Group_rank4.csv")

# ASV00006, ASV00018, ASV00001
# Mycobacterium, Rhizobacter, Candidatus Cardinium
nematode_feeding_group_rank5 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank5)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rank5)$Feeding.Group,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_feeding_group_rank5, "exported_tables/nema/nematode_Feeding.Group_rank5.csv")


# # ASV00006, ASV00018, ASV00001
nematode_feeding_group_asv <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new)$Feeding.Group,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(nematode_feeding_group_asv, "exported_tables/nema/nematode_Feeding.Group_asv.csv")

########################### Aldex2 on specific nematode groups ###########################
# Plectida taxa from selected phyloseq object
# Collapse taxonomic ranks
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new, Order == "Plectida") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida)$Genus)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida # 297 taxa and 34 samples

# Save sample data of plectida
plectida_meta <- data.frame(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida))
plectida_asv <- as.data.frame(otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida, "matrix"))
plectida_tax <- as.data.frame(tax_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida))
write.csv(plectida_meta, "exported_tables/nema/plectida_meta.csv")
write.csv(plectida_asv, "exported_tables/nema/plectida_asv.csv")
write.csv(plectida_tax, "exported_tables/nema/plectida_tax.csv")

# Import sample data of plectida with clade information
plectida_taxa_clade <- read.csv("exported_tables/nema/plectida_meta_new.csv", row.names = "Sample")
plectida_taxa_clade_data <- sample_data(plectida_taxa_clade)

# Merge new information for Plectida
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new <- merge_phyloseq(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida, plectida_taxa_clade_data)
head(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new))

####### Dorylaimida
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_sel_groups_new, Order == "Dorylaimida") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida)$Order)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida #673 taxa and 71 samples

# Save sample data of dorylaimida
dorylaimida_meta <- data.frame(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida))
dorylaimida_asv <- as.data.frame(otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida, "matrix"))
dorylaimida_tax <- as.data.frame(tax_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida))
write.csv(dorylaimida_meta, "exported_tables/nema/dorylaimida_meta.csv")
write.csv(dorylaimida_asv, "exported_tables/nema/dorylaimida_asv.csv")
write.csv(dorylaimida_tax, "exported_tables/nema/dorylaimida_tax.csv")

# Import sample data of dorylaimida with clade information
dorylaimida_taxa_clade <- read.csv("exported_tables/nema/dorylaimida_meta_new.csv", row.names = "Sample")
dorylaimida_taxa_clade_data <- sample_data(dorylaimida_taxa_clade)

# Merge new information
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new <- merge_phyloseq(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida, dorylaimida_taxa_clade_data)
head(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new))

####### Tylenchomorpha, only tylenchids
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, Infraorder == "Tylenchomorpha") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida)$Family)

# Still has aphelenchids, so let's remove it
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida <-
  subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida,
                 Family != "Aphelenchidae" &
                   Family != "Aphelenchoididae") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida)$Family)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida # 683 taxa and 89 samples

# Save sample data of tylenchida
tylenchida_meta <- data.frame(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida))
tylenchida_asv <- as.data.frame(otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida, "matrix"))
tylenchida_tax <- as.data.frame(tax_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida))
write.csv(tylenchida_meta, "exported_tables/nema/tylenchida_meta.csv")
write.csv(tylenchida_asv, "exported_tables/nema/tylenchida_asv.csv")
write.csv(tylenchida_tax, "exported_tables/nema/tylenchida_tax.csv")

# Import sample data of tylenchida with clade information
tylenchida_taxa_clade <- read.csv("exported_tables/nema/tylenchida_meta_new.csv", row.names = "Sample")
tylenchida_taxa_clade_data <- sample_data(tylenchida_taxa_clade)

# Merge new information
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new <- merge_phyloseq(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida, tylenchida_taxa_clade_data)
head(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new))

####### Aphelenchids
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, Infraorder == "Tylenchomorpha") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida)$Family)

# Still has tylenchids, so let's remove it
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida <-
  subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida,
                 Family != "Tylenchidae" &
                   Family != "Dolichodoridae" &
                   Family != "Anguinidae" &
                   Family != "Meloidogynidae") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida)$Family)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida # 741 taxa and 103 samples

# Save sample data of aphelenchida
aphelenchida_meta <- data.frame(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida))
aphelenchida_asv <- as.data.frame(otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida, "matrix"))
aphelenchida_tax <- as.data.frame(tax_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida))
write.csv(aphelenchida_meta, "exported_tables/nema/aphelenchida_meta.csv")
write.csv(aphelenchida_asv, "exported_tables/nema/aphelenchida_asv.csv")
write.csv(aphelenchida_tax, "exported_tables/nema/aphelenchida_tax.csv")

# Import sample data of aphelenchida with clade information
aphelenchida_taxa_clade <- read.csv("exported_tables/nema/aphelenchida_meta_new.csv", row.names = "Sample")
aphelenchida_taxa_clade_data <- sample_data(aphelenchida_taxa_clade)

# Merge information
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new <- merge_phyloseq(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida, aphelenchida_taxa_clade_data)
head(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new))

####### Cephalobids all genera
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, Infraorder == "Cephalobomorpha") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida)$Genus)
sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida)
smin <- min(sample_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida))

# Save sample data of cephalobida
cephalobida_meta <- data.frame(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida))
cephalobida_asv <- as.data.frame(otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida, "matrix"))
cephalobida_tax <- as.data.frame(tax_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida))
write.csv(cephalobida_meta, "exported_tables/nema/cephalobida_meta.csv")
write.csv(cephalobida_asv, "exported_tables/nema/cephalobida_asv.csv")
write.csv(cephalobida_tax, "exported_tables/nema/cephalobida_tax.csv")

# Import sample data of cephalobida with clade information
cephalobida_taxa_clade <- read.csv("exported_tables/nema/cephalobida_meta_new.csv", row.names = "Sample")
cephalobida_taxa_clade_data <- sample_data(cephalobida_taxa_clade)

# Merge information
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new <- merge_phyloseq(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida, cephalobida_taxa_clade_data)
head(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new))
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new)$Genus)


####### Cephalobids only Acrobeles and Acrobeloides
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new,
                                                                         Genus != "Cervidellus" &
                                                                         Genus != "Chiloplacus" &
                                                                         Genus != "Heterocephalobellus" &
                                                                         Genus != "Stegelletina") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols)$Genus)
smin <- min(sample_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols))
head(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols), n = 20L)

# Save sample data of acrols
acrols_meta <- data.frame(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols))
acrols_asv <- as.data.frame(otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols, "matrix"))
acrols_tax <- as.data.frame(tax_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols))
write.csv(acrols_meta, "exported_tables/nema/acrols_meta.csv")
write.csv(acrols_asv, "exported_tables/nema/acrols_asv.csv")
write.csv(acrols_tax, "exported_tables/nema/acrols_tax.csv")

####### Genus Acrobeles only
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, Genus == "Acrobeles") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles)$Genus)

# Save sample data of acrobeles
acrobeles_meta <- data.frame(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles))
acrobeles_asv <- as.data.frame(otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles, "matrix"))
acrobeles_tax <- as.data.frame(tax_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles))
write.csv(acrobeles_meta, "exported_tables/nema/acrobeles_meta.csv")
write.csv(acrobeles_asv, "exported_tables/nema/acrobeles_asv.csv")
write.csv(acrobeles_tax, "exported_tables/nema/acrobeles_tax.csv")

# Import sample data of aphelenchida with clade information
acrobeles_taxa_clade <- read.csv("exported_tables/nema/acrobeles_meta_new.csv", row.names = "Sample")
acrobeles_taxa_clade_data <- sample_data(acrobeles_taxa_clade)

# Merge information
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new <- merge_phyloseq(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles, acrobeles_taxa_clade_data)
head(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new))

# Exclude unknown acrobeles clades
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new, Clade != "Unknown") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new)$Clade)

####### Genus Acrobeloides only
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new, Genus == "Acrobeloides") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides)$Clade)
head(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides))

# Exclude unknown Acrobeloides clades
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides, Clade != "Unknown") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new)$Clade)

# Save sample data of Acrobeloides
acrobeloides_meta <- data.frame(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides))
acrobeloides_asv <- as.data.frame(otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides, "matrix"))
acrobeloides_tax <- as.data.frame(tax_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides))
write.csv(acrobeloides_meta, "exported_tables/nema/acrobeloides_meta.csv")
write.csv(acrobeloides_asv, "exported_tables/nema/acrobeloides_asv.csv")
write.csv(acrobeloides_tax, "exported_tables/nema/acrobeloides_tax.csv")

####### Panagrolaimids
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, Infraorder == "Panagrolaimomorpha") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida)$Genus)

# Save sample data of Panagrolaimida
panagrolaimida_meta <- data.frame(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida))
panagrolaimida_asv <- as.data.frame(otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida, "matrix"))
panagrolaimida_tax <- as.data.frame(tax_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida))
write.csv(panagrolaimida_meta, "exported_tables/nema/panagrolaimida_meta.csv")
write.csv(panagrolaimida_asv, "exported_tables/nema/panagrolaimida_asv.csv")
write.csv(panagrolaimida_tax, "exported_tables/nema/panagrolaimida_tax.csv")

# Import sample data of panagrolaimida with clade information
panagrolaimida_taxa_clade <- read.csv("exported_tables/nema/panagrolaimida_meta_new.csv", row.names = "Sample")
panagrolaimida_taxa_clade_data <- sample_data(panagrolaimida_taxa_clade)

# Merge information
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new <- merge_phyloseq(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida, panagrolaimida_taxa_clade_data)
head(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new))

# Exclude unknown panagrolaimida clades
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new, Clade != "Unknown") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new)$Clade)

#### Rhabditida
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, Infraorder == "Rhabditomorpha") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida)$Genus)

# Save sample data of rhabditida
rhabditida_meta <- data.frame(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida))
rhabditida_asv <- as.data.frame(otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida, "matrix"))
rhabditida_tax <- as.data.frame(tax_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida))
write.csv(rhabditida_meta, "exported_tables/nema/rhabditida_meta.csv")
write.csv(rhabditida_asv, "exported_tables/nema/rhabditida_asv.csv")
write.csv(rhabditida_tax, "exported_tables/nema/rhabditida_tax.csv")

# Import sample data of rhabditida with clade information
rhabditida_taxa_clade <- read.csv("exported_tables/nema/rhabditida_meta_new.csv", row.names = "Sample")
rhabditida_taxa_clade_data <- sample_data(rhabditida_taxa_clade)

# Merge information
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new <- merge_phyloseq(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida, rhabditida_taxa_clade_data)
head(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new))

# Exclude unknown rhabditida clades
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new, Clade != "unknown") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new)$Clade)

#### Collapse Plectida phyloseq by ranks
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank1 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new, taxrank = "Phylum")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank2 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new, taxrank = "Class")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank3 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new, taxrank = "Order")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank4 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new, taxrank = "Family")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank5 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new, taxrank = "Genus")

#### Collapse Dorylaimida phyloseq by ranks
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank1 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new, taxrank = "Phylum")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank2 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new, taxrank = "Class")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank3 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new, taxrank = "Order")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank4 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new, taxrank = "Family")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank5 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new, taxrank = "Genus")

#### Collapse Tylenchida/tylenchids phyloseq by ranks
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank1 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new, taxrank = "Phylum")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank2 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new, taxrank = "Class")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank3 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new, taxrank = "Order")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank4 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new, taxrank = "Family")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank5 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new, taxrank = "Genus")

#### Collapse Aphelenchida phyloseq by ranks
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank1 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new, taxrank = "Phylum")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank2 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new, taxrank = "Class")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank3 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new, taxrank = "Order")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank4 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new, taxrank = "Family")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank5 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new, taxrank = "Genus")

#### Collapse Cephalobida phyloseq by ranks
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank1 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new, taxrank = "Phylum")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank2 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new, taxrank = "Class")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank3 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new, taxrank = "Order")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank4 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new, taxrank = "Family")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank5 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new, taxrank = "Genus")

#### Collapse Acrols phyloseq by ranks
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank1 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols, taxrank = "Phylum")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank2 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols, taxrank = "Class")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank3 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols, taxrank = "Order")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank4 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols, taxrank = "Family")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank5 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols, taxrank = "Genus")

#### Collapse Acrobeloides phyloseq by ranks
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank1 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new, taxrank = "Phylum")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank2 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new, taxrank = "Class")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank3 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new, taxrank = "Order")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank4 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new, taxrank = "Family")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank5 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new, taxrank = "Genus")

#### Collapse Acrobeles phyloseq by ranks
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank1 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new, taxrank = "Phylum")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank2 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new, taxrank = "Class")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank3 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new, taxrank = "Order")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank4 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new, taxrank = "Family")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank5 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new, taxrank = "Genus")

# Collapse Panagrolaimida phyloseq by ranks
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank1 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new, taxrank = "Phylum")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank2 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new, taxrank = "Class")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank3 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new, taxrank = "Order")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank4 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new, taxrank = "Family")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank5 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new, taxrank = "Genus")

# Collapse Rhabditida phyloseq by ranks
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank1 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new, taxrank = "Phylum")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank2 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new, taxrank = "Class")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank3 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new, taxrank = "Order")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank4 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new, taxrank = "Family")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank5 <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new, taxrank = "Genus")

######################################## Aldex2 on specific nematode groups ######################################## 
######################################## Plectida nematodes and taxonomic ranks ######################################## 
############### Across Plectida clades and taxonomic ranks
# None differentiated abundant
plectida_new_rank1_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank1)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank1)$Clade,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank1_clade, "exported_tables/nema/plectida_new_rank1_clade.csv")

# ASV00027 - Gammaproteobacteria
plectida_new_rank2_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank2)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank2)$Clade,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank2_clade, "exported_tables/nema/plectida_new_rank2_clade.csv")

# None differentiated abundant
plectida_new_rank3_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank3)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank3)$Clade,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank3_clade, "exported_tables/nema/plectida_new_rank3_clade.csv")

# None differentiated abundant
plectida_new_rank4_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank4)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank4)$Clade,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank4_clade, "exported_tables/nema/plectida_new_rank4_clade.csv")

# None differentiated abundant
plectida_new_rank5_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank5)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank5)$Clade,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)

write.csv(plectida_new_rank5_clade, "exported_tables/nema/plectida_new_rank5_clade.csv")

# None differentiated abundant
plectida_new_asv_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new)$Clade,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)

write.csv(plectida_new_asv_clade, "exported_tables/nema/plectida_new_asv_clade.csv")

############### Across Plectida genera and taxonomic ranks
############### Plectus vs. Tylocephalus
# ASV00010 - Bacteroidota
plectida_new_rank1_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank1)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank1)$Genus,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank1_genus, "exported_tables/nema/plectida_new_rank1_genus.csv")

# ASV00010, ASV00027, ASV00085
# Bacteroidia, Gammaproteobacteria, Alphaproteobacteria
plectida_new_rank2_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank2)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank2)$Genus,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank2_genus, "exported_tables/nema/plectida_new_rank2_genus.csv")

# ASV00027, ASV00145
# Pseudomonadales, Burkholderiales
plectida_new_rank3_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank3)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank3)$Genus,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank3_genus, "exported_tables/nema/plectida_new_rank3_genus.csv")

# None differentiated abundant
plectida_new_rank4_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank4)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank4)$Genus,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank4_genus, "exported_tables/nema/plectida_new_rank4_genus.csv")

# None differentiated abundant
plectida_new_rank5_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank5)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank5)$Genus,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank5_genus, "exported_tables/nema/plectida_new_rank5_genus.csv")


# None differentiated abundant
plectida_new_asv_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new)$Genus,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_asv_genus, "exported_tables/nema/plectida_new_asv_genus.csv")

################# Across site and taxonomic ranks
# None differentiated abundant
plectida_new_rank1_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank1)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank1)$Sample.Site,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank1_site, "exported_tables/nema/plectida_new_rank1_site.csv")

# None differentiated abundant
plectida_new_rank2_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank2)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank2)$Sample.Site,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank2_site, "exported_tables/nema/plectida_new_rank2_site.csv")

# None differentiated abundant
plectida_new_rank3_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank3)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank3)$Sample.Site,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank3_site, "exported_tables/nema/plectida_new_rank3_site.csv")

# None differentiated abundant
plectida_new_rank4_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank4)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank4)$Sample.Site,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank4_site, "exported_tables/nema/plectida_new_rank4_site.csv")

# None differentiated abundant
plectida_new_rank5_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank5)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank5)$Sample.Site,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank5_site, "exported_tables/nema/plectida_new_rank5_site.csv")

# None differentiated abundant
plectida_new_asv_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new)),
                                         phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new)$Sample.Site,
                                         mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_asv_site, "exported_tables/nema/plectida_new_asv_site.csv")

################# Across habitat and taxonomic ranks
# None differentiated abundant
plectida_new_rank1_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank1)),
                                         phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank1)$Habitat,
                                         mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank1_habitat, "exported_tables/nema/plectida_new_rank1_habitat.csv")

# None differentiated abundant
plectida_new_rank2_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank2)),
                                         phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank2)$Habitat,
                                         mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank2_habitat, "exported_tables/nema/plectida_new_rank2_habitat.csv")

# None differentiated abundant
plectida_new_rank3_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank3)),
                                         phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank3)$Habitat,
                                         mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank3_habitat, "exported_tables/nema/plectida_new_rank3_habitat.csv")

# None differentiated abundant
plectida_new_rank4_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank4)),
                                         phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank4)$Habitat,
                                         mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank4_habitat, "exported_tables/nema/plectida_new_rank4_habitat.csv")

# None differentiated abundant
plectida_new_rank5_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank5)),
                                         phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank5)$Habitat,
                                         mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank5_habitat, "exported_tables/nema/plectida_new_rank5_habitat.csv")

# None differentiated abundant
plectida_new_asv_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new)$Habitat,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_asv_habitat, "exported_tables/nema/plectida_new_asv_habitat.csv")

############### Between soil type and taxonomic ranks
# ASV00006 - Actinobacteriota
plectida_new_rank1_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank1)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank1)$Soil.Type,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank1_soil_type, "exported_tables/nema/plectida_new_rank1_soil_type.csv")

# ASV00027 - Gammaproteobacteria
plectida_new_rank2_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank2)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank2)$Soil.Type,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank2_soil_type, "exported_tables/nema/plectida_new_rank2_soil_type.csv")

# ASV00145 - Burkholderiales
plectida_new_rank3_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank3)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank3)$Soil.Type,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank3_soil_type, "exported_tables/nema/plectida_new_rank3_soil_type.csv")

# None differentiated abundant
plectida_new_rank4_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank4)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank4)$Soil.Type,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank4_soil_type, "exported_tables/nema/plectida_new_rank4_soil_type.csv")

# None differentiated abundant
plectida_new_rank5_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank5)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new_rank5)$Soil.Type,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_rank5_soil_type, "exported_tables/nema/plectida_new_rank5_soil_type.csv")

# None differentiated abundant
plectida_new_asv_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_plectida_new)$Soil.Type,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(plectida_new_asv_soil_type, "exported_tables/nema/plectida_new_asv_soil_type.csv")

######################################## Cephalobida nematodes and taxonomic ranks ######################################## 
# Across Cephalobida genera and taxonomic ranks
# ASV00004
# Firmicutes
cephalobida_new_rank1_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank1)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank1)$Genus,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank1_genus, "exported_tables/nema/cephalobida_new_rank1_genus.csv")

# ASV00004, ASV00005
# Bacilli, Alphaproteobacteria
cephalobida_new_rank2_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank2)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank2)$Genus,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank2_genus, "exported_tables/nema/cephalobida_new_rank2_genus.csv")

# ASV00006, ASV00007
# Corynebacteriales, Sphingomonadales
cephalobida_new_rank3_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank3)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank3)$Genus,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank3_genus, "exported_tables/nema/cephalobida_new_rank3_genus.csv")

# ASV00006, ASV00007
# Mycobacteriaceae, Sphingomonadaceae
cephalobida_new_rank4_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank4)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank4)$Genus,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank4_genus, "exported_tables/nema/cephalobida_new_rank4_genus.csv")

# ASV00006
# Mycobacterium
cephalobida_new_rank5_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank5)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank5)$Genus,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank5_genus, "exported_tables/nema/cephalobida_new_rank5_genus.csv")

cephalobida_new_asv_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new)$Genus,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_asv_genus, "exported_tables/nema/cephalobida_new_asv_genus.csv")

######################### Across site and taxonomic ranks
# ASV00011, ASV00005
# Verrucomicrobiota, Proteobacteria
cephalobida_new_rank1_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank1)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank1)$Sample.Site,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank1_site, "exported_tables/nema/cephalobida_new_rank1_site.csv")

# ASV00005, ASV00011, ASV00018
# Alphaproteobacteria, Verrucomicrobiae, Gammaproteobacteria
cephalobida_new_rank2_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank2)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank2)$Sample.Site,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank2_site, "exported_tables/nema/cephalobida_new_rank2_site.csv")

# ASV00011, ASV00018, ASV00007
# Opitutales, Burkholderiales, Sphingomonadales
cephalobida_new_rank3_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank3)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank3)$Sample.Site,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank3_site, "exported_tables/nema/cephalobida_new_rank3_site.csv")

# ASV00018, ASV00011, ASV00007, ASV00006, ASV00023
# Comamonadaceae, Opitutaceae, Sphingomonadaceae, Mycobacteriaceae, Enterococcaceae
cephalobida_new_rank4_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank4)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank4)$Sample.Site,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank4_site, "exported_tables/nema/cephalobida_new_rank4_site.csv")

# ASV00018, ASV00011, ASV00006, ASV00007, ASV00024, ASV00023
# Rhizobacter, Lacunisphaera Mycobacterium, Novosphingobium, Limnohabitans, Tetragenococcus
cephalobida_new_rank5_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank5)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank5)$Sample.Site,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank5_site, "exported_tables/nema/cephalobida_new_rank5_site.csv")

# ASV00018, ASV00011, ASV00006, ASV00007
cephalobida_new_asv_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new)$Sample.Site,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_asv_site, "exported_tables/nema/cephalobida_new_asv_site.csv")

######################### Across habitat and taxonomic ranks
# ASV00011 - Verrucomicrobiota
cephalobida_new_rank1_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank1)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank1)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank1_habitat, "exported_tables/nema/cephalobida_new_rank1_habitat.csv")

# ASV00006, ASV00011
# Actinobacteria, Verrucomicrobiae
cephalobida_new_rank2_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank2)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank2)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank2_habitat, "exported_tables/nema/cephalobida_new_rank2_habitat.csv")

# ASV00011, ASV00007, ASV00006, ASV00018
# Opitutales, Sphingomonadales, Corynebacteriales, Burkholderiales
cephalobida_new_rank3_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank3)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank3)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank3_habitat, "exported_tables/nema/cephalobida_new_rank3_habitat.csv")

# ASV00011, ASV00007, ASV00006
# Opitutaceae, Sphingomonadaceae, Mycobacteriaceae
cephalobida_new_rank4_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank4)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank4)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank4_habitat, "exported_tables/nema/cephalobida_new_rank4_habitat.csv")

# ASV00006 - Mycobacterium
# ASV00011 - Lacunisphaera
cephalobida_new_rank5_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank5)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank5)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank5_habitat, "exported_tables/nema/cephalobida_new_rank5_habitat.csv")

# ASV00011, ASV00006, ASV00007
cephalobida_new_asv_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new)$Habitat,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_asv_habitat, "exported_tables/nema/cephalobida_new_asv_habitat.csv")


##################### Across soil type and taxonomic ranks
# None
cephalobida_new_rank1_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank1)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank1)$Soil.Type,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank1_soil_type, "exported_tables/nema/cephalobida_new_rank1_soil_type.csv")

# ASV00005 - Alphaproteobacteria
cephalobida_new_rank2_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank2)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank2)$Soil.Type,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank2_soil_type, "exported_tables/nema/cephalobida_new_rank2_soil_type.csv")

# ASV00007, ASV00004, ASV00018
# Sphingomonadales, Staphylococcales, Burkholderiales
cephalobida_new_rank3_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank3)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank3)$Soil.Type,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank3_soil_type, "exported_tables/nema/cephalobida_new_rank3_soil_type.csv")

# ASV00004, ASV00006, ASV00007, ASV00018
# Staphylococcaceae, Mycobacteriaceae, Sphingomonadaceae, Comamonadaceae
cephalobida_new_rank4_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank4)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank4)$Soil.Type,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank4_soil_type, "exported_tables/nema/cephalobida_new_rank4_soil_type.csv")

# ASV00006, ASV00007, ASV00018
# Mycobacterium, Novosphingobium, Rhizobacter
cephalobida_new_rank5_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank5)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new_rank5)$Soil.Type,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_rank5_soil_type, "exported_tables/nema/cephalobida_new_rank5_soil_type.csv")

# ASV00006, ASV00007, ASV00018
cephalobida_new_asv_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new)),
                                                 phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_cephalobida_new)$Soil.Type,
                                                 mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(cephalobida_new_asv_soil_type, "exported_tables/nema/cephalobida_new_asv_soil_type.csv")

######################################## Acrobeles vs. Acrobeloides nematodes and taxonomic ranks ######################################## 
############# Across Acrols genera and taxonomic ranks
# ASV00010 - Bacteroidota
acrols_rank1_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank1)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank1)$Genus,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank1_genus, "exported_tables/nema/acrols_rank1_genus.csv")

# ASV00005, ASV00010
# Alphaproteobacteria, Bacteroidia
acrols_rank2_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank2)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank2)$Genus,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank2_genus, "exported_tables/nema/acrols_rank2_genus.csv")

# None
acrols_rank3_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank3)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank3)$Genus,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank3_genus, "exported_tables/nema/acrols_rank3_genus.csv")

# None
acrols_rank4_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank4)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank4)$Genus,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank4_genus, "exported_tables/nema/acrols_rank4_genus.csv")

# None
acrols_rank5_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank5)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank5)$Genus,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank5_genus, "exported_tables/nema/acrols_rank5_genus.csv")


# None
acrols_asv_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols)),
                                    phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols)$Genus,
                                    mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_asv_genus, "exported_tables/nema/acrols_asv_genus.csv")

##################### Across site and taxonomic ranks
# ASV00005, ASV00011
# Proteobacteria, Verrucomicrobiota
acrols_rank1_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank1)),
                                   phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank1)$Sample.Site,
                                   mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank1_site, "exported_tables/nema/acrols_rank1_site.csv")

# ASV00005, ASV00018, ASV00011
# Alphaproteobacteria, Gammaproteobacteria, Verrucomicrobiae
acrols_rank2_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank2)),
                                   phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank2)$Sample.Site,
                                   mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank2_site, "exported_tables/nema/acrols_rank2_site.csv")

# ASV00011, ASV00018, ASV00007
# Opitutales, Burkholderiales, Sphingomonadales
acrols_rank3_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank3)),
                                   phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank3)$Sample.Site,
                                   mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank3_site, "exported_tables/nema/acrols_rank3_site.csv")

# ASV00018, ASV00011, ASV00007
# Comamonadaceae, Opitutaceae, Sphingomonadaceae
acrols_rank4_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank4)),
                                   phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank4)$Sample.Site,
                                   mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank4_site, "exported_tables/nema/acrols_rank4_site.csv")

# ASV00018, ASV00011, ASV00024
# Rhizobacter, Lacunisphaera, Limnohabitans
acrols_rank5_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank5)),
                                   phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank5)$Sample.Site,
                                   mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank5_site, "exported_tables/nema/acrols_rank5_site.csv")

# ASV00018, ASV00011, # ASV00007
acrols_asv_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols)),
                                 phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols)$Sample.Site,
                                 mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_asv_site, "exported_tables/nema/acrols_asv_site.csv")

##################### Across habitat and taxonomic ranks
# ASV00011
# Verrucomicrobiota
acrols_rank1_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank1)),
                                      phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank1)$Habitat,
                                      mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank1_habitat, "exported_tables/nema/acrols_rank1_habitat.csv")

# ASV00006, ASV00011
# Actinobacteria, Verrucomicrobiae
acrols_rank2_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank2)),
                                      phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank2)$Habitat,
                                      mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank2_habitat, "exported_tables/nema/acrols_rank2_habitat.csv")

# ASV00007, ASV00010, ASV00011
# Sphingomonadales, Flavobacteriales, Opitutales
acrols_rank3_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank3)),
                                      phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank3)$Habitat,
                                      mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank3_habitat, "exported_tables/nema/acrols_rank3_habitat.csv")

# ASV00007, ASV00010, ASV00011
# Sphingomonadaceae, Weeksellaceae, Opitutaceae
acrols_rank4_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank4)),
                                      phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank4)$Habitat,
                                      mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank4_habitat, "exported_tables/nema/acrols_rank4_habitat.csv")

# ASV00011 - Lacunisphaera
acrols_rank5_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank5)),
                                      phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank5)$Habitat,
                                      mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank5_habitat, "exported_tables/nema/acrols_rank5_habitat.csv")

# ASV00011
acrols_asv_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols)),
                                      phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols)$Habitat,
                                      mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_asv_habitat, "exported_tables/nema/acrols_asv_habitat.csv")

################### Across soil type and taxonomic ranks
# None
acrols_rank1_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank1)),
                                                 phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank1)$Soil.Type,
                                                 mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank1_soil_type, "exported_tables/nema/acrols_rank1_soil_type.csv")

# ASV00005 - Alphaproteobacteria
acrols_rank2_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank2)),
                                                 phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank2)$Soil.Type,
                                                 mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank2_soil_type, "exported_tables/nema/acrols_rank2_soil_type.csv")

# ASV00004, ASV00007
# Staphylococcales, Sphingomonadales
acrols_rank3_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank3)),
                                                 phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank3)$Soil.Type,
                                                 mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank3_soil_type, "exported_tables/nema/acrols_rank3_soil_type.csv")

# ASV00004, ASV00006, ASV00007, ASV00018
# Staphylococcaceae, Mycobacteriaceae, Sphingomonadaceae, Comamonadaceae
acrols_rank4_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank4)),
                                                 phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank4)$Soil.Type,
                                                 mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank4_soil_type, "exported_tables/nema/acrols_rank4_soil_type.csv")

# ASV00004, ASV00006, ASV00007, ASV00018
# Staphylococcus, Mycobacterium, Novosphingobium, Rhizobacter
acrols_rank5_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank5)),
                                                 phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols_rank5)$Soil.Type,
                                                 mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_rank5_soil_type, "exported_tables/nema/acrols_rank5_soil_type.csv")

# ASV00004, ASV00006, ASV00007, ASV00018
acrols_asv_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols)),
                                        phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrols)$Soil.Type,
                                        mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrols_asv_soil_type, "exported_tables/nema/acrols_asv_soil_type.csv")

######################################## Acrobeles nematodes and taxonomic ranks ######################################## 
# Across Acrobeles clades and taxonomic ranks
# None differentiated abundant
acrobeles_new_rank1_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank1)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank1)$Clade,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank1_clade, "exported_tables/nema/acrobeles_new_rank1_clade.csv")

# None differentiated abundant
acrobeles_new_rank2_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank2)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank2)$Clade,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank2_clade, "exported_tables/nema/acrobeles_new_rank2_clade.csv")

# None differentiated abundant
acrobeles_new_rank3_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank3)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank3)$Clade,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank3_clade, "exported_tables/nema/acrobeles_new_rank3_clade.csv")

# None differentiated abundant
acrobeles_new_rank4_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank4)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank4)$Clade,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank4_clade, "exported_tables/nema/acrobeles_new_rank4_clade.csv")

# None differentiated abundant
acrobeles_new_rank5_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank5)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank5)$Clade,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank5_clade, "exported_tables/nema/acrobeles_new_rank5_clade.csv")

# None differentiated abundant
acrobeles_new_asv_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new)$Clade,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_asv_clade, "exported_tables/nema/acrobeles_new_asv_clade.csv")

# Across sites and taxonomic ranks
# None differentiated abundant
acrobeles_new_rank1_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank1)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank1)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank1_site, "exported_tables/nema/acrobeles_new_rank1_site.csv")

# None differentiated abundant
acrobeles_new_rank2_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank2)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank2)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank2_site, "exported_tables/nema/acrobeles_new_rank2_site.csv")

# ASV00007, ASV00004 - Sphingomonadales, Staphylococcales
acrobeles_new_rank3_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank3)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank3)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank3_site, "exported_tables/nema/acrobeles_new_rank3_site.csv")

# ASV00007, ASV00004 - Sphingomonadaceae, Staphylococcaceae
acrobeles_new_rank4_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank4)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank4)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank4_site, "exported_tables/nema/acrobeles_new_rank4_site.csv")

# ASV00025, ASV00004 - Ralstonia, Staphylococcus
acrobeles_new_rank5_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank5)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank5)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank5_site, "exported_tables/nema/acrobeles_new_rank5_site.csv")

# ASV00025, ASV00004
acrobeles_new_asv_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_asv_site, "exported_tables/nema/acrobeles_new_asv_site.csv")

# Across habitats and taxonomic ranks
# None differentiated abundant
acrobeles_new_rank1_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank1)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank1)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank1_habitat, "exported_tables/nema/acrobeles_new_rank1_habitat.csv")

# None differentiated abundant
acrobeles_new_rank2_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank2)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank2)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank2_habitat, "exported_tables/nema/acrobeles_new_rank2_habitat.csv")

# ASV00007, ASV00004 - Sphingomonadales, Staphylococcales
acrobeles_new_rank3_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank3)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank3)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank3_habitat, "exported_tables/nema/acrobeles_new_rank3_habitat.csv")

# ASV00007, ASV00004 - Sphingomonadaceae, Staphylococcaceae
acrobeles_new_rank4_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank4)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank4)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank4_habitat, "exported_tables/nema/acrobeles_new_rank4_habitat.csv")

# ASV00025, ASV00004 - Ralstonia, Staphylococcus
acrobeles_new_rank5_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank5)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank5)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank5_habitat, "exported_tables/nema/acrobeles_new_rank5_habitat.csv")

# ASV00025, ASV00004
acrobeles_new_asv_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new)$Habitat,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_asv_habitat, "exported_tables/nema/acrobeles_new_asv_habitat.csv")

# Between soil type and taxonomic ranks
# ASV00004 - Firmicutes
acrobeles_new_rank1_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank1)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank1)$Soil.Type,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank1_soil_type, "exported_tables/nema/acrobeles_new_rank1_soil_type.csv")

# ASV00004 - Bacilli
acrobeles_new_rank2_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank2)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank2)$Soil.Type,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank2_soil_type, "exported_tables/nema/acrobeles_new_rank2_soil_type.csv")

# ASV00004, ASV00023, ASV00007
# Staphylococcales, Lactobacillales, Sphingomonadales
acrobeles_new_rank3_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank3)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank3)$Soil.Type,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank3_soil_type, "exported_tables/nema/acrobeles_new_rank3_soil_type.csv")

# ASV00004, ASV00006, ASV00007
# Staphylococcaceae, Mycobacteriaceae, Sphingomonadaceae
acrobeles_new_rank4_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank4)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank4)$Soil.Type,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank4_soil_type, "exported_tables/nema/acrobeles_new_rank4_soil_type.csv")

# ASV00004, ASV00006, ASV00025 - Staphylococcus, Mycobacterium, Ralstonia
acrobeles_new_rank5_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank5)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new_rank5)$Soil.Type,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_rank5_soil_type, "exported_tables/nema/acrobeles_new_rank5_soil_type.csv")

#
acrobeles_new_asv_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeles_new)$Soil.Type,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeles_new_asv_soil_type, "exported_tables/nema/acrobeles_new_asv_soil_type.csv")
######################################## Acrobeloides nematodes and taxonomic ranks ######################################## 
# Across Acrobeloides clades and taxonomic ranks
# None differentiated abundant
acrobeloides_new_rank1_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank1)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank1)$Clade,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank1_clade, "exported_tables/nema/acrobeloides_new_rank1_clade.csv")

# ASV00011 -Verrucomicrobiota
acrobeloides_new_rank2_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank2)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank2)$Clade,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank2_clade, "exported_tables/nema/acrobeloides_new_rank2_clade.csv")

# None differentiated abundant
acrobeloides_new_rank3_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank3)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank3)$Clade,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank3_clade, "exported_tables/nema/acrobeloides_new_rank3_clade.csv")

# None differentiated abundant
acrobeloides_new_rank4_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank4)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank4)$Clade,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank4_clade, "exported_tables/nema/acrobeloides_new_rank4_clade.csv")

# None differentiated abundant
acrobeloides_new_rank5_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank5)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank5)$Clade,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank5_clade, "exported_tables/nema/acrobeloides_new_rank5_clade.csv")

# None differentiated abundant
acrobeloides_new_asv_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new)$Clade,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_asv_clade, "exported_tables/nema/acrobeloides_new_asv_clade.csv")

# Across sites and taxonomic ranks
# None differentiated abundant
acrobeloides_new_rank1_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank1)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank1)$Sample.Site,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank1_site, "exported_tables/nema/acrobeloides_new_rank1_site.csv")

# None differentiated abundant
acrobeloides_new_rank2_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank2)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank2)$Sample.Site,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank2_site, "exported_tables/nema/acrobeloides_new_rank2_site.csv")

# ASV00006 - Corynebacteriales
acrobeloides_new_rank3_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank3)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank3)$Sample.Site,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank3_site, "exported_tables/nema/acrobeloides_new_rank3_site.csv")

# None
acrobeloides_new_rank4_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank4)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank4)$Sample.Site,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank4_site, "exported_tables/nema/acrobeloides_new_rank4_site.csv")

# ASV00007 - Novosphingobium
acrobeloides_new_rank5_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank5)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank5)$Sample.Site,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank5_site, "exported_tables/nema/acrobeloides_new_rank5_site.csv")

# ASV00007 - Novosphingobium
acrobeloides_new_asv_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new)$Sample.Site,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_asv_site, "exported_tables/nema/acrobeloides_new_asv_site.csv")

# Across habitats and taxonomic ranks
# None differentiated abundant
acrobeloides_new_rank1_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank1)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank1)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank1_habitat, "exported_tables/nema/acrobeloides_new_rank1_habitat.csv")

# None differentiated abundant
acrobeloides_new_rank2_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank2)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank2)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank2_habitat, "exported_tables/nema/acrobeloides_new_rank2_habitat.csv")

# ASV00006 - Corynebacteriales
acrobeloides_new_rank3_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank3)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank3)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank3_habitat, "exported_tables/nema/acrobeloides_new_rank3_habitat.csv")

# None
acrobeloides_new_rank4_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank4)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank4)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank4_habitat, "exported_tables/nema/acrobeloides_new_rank4_habitat.csv")

# ASV00007 - Novosphingobium
acrobeloides_new_rank5_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank5)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank5)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank5_habitat, "exported_tables/nema/acrobeloides_new_rank5_habitat.csv")

# ASV00007 - Novosphingobium
acrobeloides_new_asv_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new)$Habitat,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_asv_habitat, "exported_tables/nema/acrobeloides_new_asv_habitat.csv")

# Between soil type and taxonomic ranks
# None differentiated abundant
acrobeloides_new_rank1_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank1)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank1)$Soil.Type,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank1_soil_type, "exported_tables/nema/acrobeloides_new_rank1_soil_type.csv")

# ASV00086 - Alphaproteobacteria
acrobeloides_new_rank2_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank2)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank2)$Soil.Type,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank2_soil_type, "exported_tables/nema/acrobeloides_new_rank2_soil_type.csv")

# ASV00007 - Sphingomonadales
acrobeloides_new_rank3_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank3)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank3)$Soil.Type,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank3_soil_type, "exported_tables/nema/acrobeloides_new_rank3_soil_type.csv")

# ASV00007 - Sphingomonadaceae
acrobeloides_new_rank4_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank4)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank4)$Soil.Type,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank4_soil_type, "exported_tables/nema/acrobeloides_new_rank4_soil_type.csv")

# ASV00007 - Novosphingobium
acrobeloides_new_rank5_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank5)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new_rank5)$Soil.Type,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_rank5_soil_type, "exported_tables/nema/acrobeloides_new_rank5_soil_type.csv")

# ASV00007 - Novosphingobium
acrobeloides_new_asv_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new)),
                                                  phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_acrobeloides_new)$Soil.Type,
                                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(acrobeloides_new_asv_soil_type, "exported_tables/nema/acrobeloides_new_asv_soil_type.csv")

######################################## Aphelenchida nematodes and taxonomic ranks ######################################## 
# Aphelenchida genera
# ASV00004 - Firmicutes
aphelenchida_new_rank1_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank1)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank1)$Genus,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank1_genus, "exported_tables/nema/aphelenchida_new_rank1_genus.csv")

# ASV00004, ASV00021 - Bacilli, Alphaproteobacteria
aphelenchida_new_rank2_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank2)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank2)$Genus,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank2_genus, "exported_tables/nema/aphelenchida_new_rank2_genus.csv")

# None differentiated abundant
aphelenchida_new_rank3_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank3)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank3)$Genus,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank3_genus, "exported_tables/nema/aphelenchida_new_rank3_genus.csv")

# None differentiated abundant
aphelenchida_new_rank4_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank4)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank4)$Genus,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank4_genus, "exported_tables/nema/aphelenchida_new_rank4_genus.csv")

# ASV00159 - Sphingomonas (p-value: 0.0545059)
aphelenchida_new_rank5_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank5)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank5)$Genus,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank5_genus, "exported_tables/nema/aphelenchida_new_rank5_genus.csv")


# ASV00159 - Sphingomonas (p-value: 0.0545059)
aphelenchida_new_asv_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new)$Genus,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_asv_genus, "exported_tables/nema/aphelenchida_new_asv_genus.csv")

# Across sites and taxonomic ranks
# ASV00004 - Firmicutes
aphelenchida_new_rank1_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank1)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank1)$Sample.Site,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank1_site, "exported_tables/nema/aphelenchida_new_rank1_site.csv")

# ASV00004 - Bacilli
aphelenchida_new_rank2_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank2)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank2)$Sample.Site,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank2_site, "exported_tables/nema/aphelenchida_new_rank2_site.csv")

# ASV00004 - Staphylococcales
aphelenchida_new_rank3_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank3)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank3)$Sample.Site,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank3_site, "exported_tables/nema/aphelenchida_new_rank3_site.csv")

# ASV00004 - Staphylococcaceae
aphelenchida_new_rank4_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank4)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank4)$Sample.Site,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank4_site, "exported_tables/nema/aphelenchida_new_rank4_site.csv")

# ASV00004 - Staphylococcus
aphelenchida_new_rank5_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank5)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank5)$Sample.Site,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank5_site, "exported_tables/nema/aphelenchida_new_rank5_site.csv")


# ASV00004 - Staphylococcus
aphelenchida_new_asv_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new)$Sample.Site,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_asv_site, "exported_tables/nema/aphelenchida_new_asv_site.csv")

# Across habitats and taxonomic ranks
# ASV00004 - Firmicutes
aphelenchida_new_rank1_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank1)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank1)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank1_habitat, "exported_tables/nema/aphelenchida_new_rank1_habitat.csv")

# ASV00004 - Bacilli
aphelenchida_new_rank2_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank2)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank2)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank2_habitat, "exported_tables/nema/aphelenchida_new_rank2_habitat.csv")

# ASV00004 - Staphylococcales
aphelenchida_new_rank3_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank3)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank3)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank3_habitat, "exported_tables/nema/aphelenchida_new_rank3_habitat.csv")

# ASV00004 - Staphylococcaceae
aphelenchida_new_rank4_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank4)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank4)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank4_habitat, "exported_tables/nema/aphelenchida_new_rank4_habitat.csv")

# ASV00004 - Staphylococcus
aphelenchida_new_rank5_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank5)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank5)$Habitat,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank5_habitat, "exported_tables/nema/aphelenchida_new_rank5_habitat.csv")


# ASV00004 - Staphylococcus
aphelenchida_new_asv_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new)$Habitat,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_asv_habitat, "exported_tables/nema/aphelenchida_new_asv_habitat.csv")

# Between soil type and taxonomic ranks
# None differentiated abundant
aphelenchida_new_rank1_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank1)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank1)$Soil.Type,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank1_soil_type, "exported_tables/nema/aphelenchida_new_rank1_soil_type.csv")

# ASV00004 - Bacilli
aphelenchida_new_rank2_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank2)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank2)$Soil.Type,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank2_soil_type, "exported_tables/nema/aphelenchida_new_rank2_soil_type.csv")

# ASV00004 - Staphylococcales
aphelenchida_new_rank3_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank3)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank3)$Soil.Type,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank3_soil_type, "exported_tables/nema/aphelenchida_new_rank3_soil_type.csv")

# ASV00004, ASV00006 - Staphylococcaceae, Mycobacteriaceae
aphelenchida_new_rank4_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank4)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank4)$Soil.Type,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank4_soil_type, "exported_tables/nema/aphelenchida_new_rank4_soil_type.csv")

# ASV00004, ASV00006, ASV00025 - Staphylococcus, Mycobacterium, Ralstonia
aphelenchida_new_rank5_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank5)),
                                               phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new_rank5)$Soil.Type,
                                               mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_rank5_soil_type, "exported_tables/nema/aphelenchida_new_rank5_soil_type.csv")


# ASV00004, ASV00006, ASV00025 - Staphylococcus, Mycobacterium, Ralstonia
aphelenchida_new_asv_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new)),
                                                  phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_aphelenchida_new)$Soil.Type,
                                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(aphelenchida_new_asv_soil_type, "exported_tables/nema/aphelenchida_new_asv_soil_type.csv")

######################################## Rhabditida nematodes and taxonomic ranks ######################################## 
# Across rhabditida clades and taxonomic ranks
# None differentiated abundant
rhabditida_new_rank1_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank1)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank1)$Clade,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank1_clade, "exported_tables/nema/rhabditida_new_rank1_clade.csv")

# None differentiated abundant
rhabditida_new_rank2_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank2)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank2)$Clade,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank2_clade, "exported_tables/nema/rhabditida_new_rank2_clade.csv")

# None differentiated abundant
rhabditida_new_rank3_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank3)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank3)$Clade,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank3_clade, "exported_tables/nema/rhabditida_new_rank3_clade.csv")

# None differentiated abundant
rhabditida_new_rank4_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank4)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank4)$Clade,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank4_clade, "exported_tables/nema/rhabditida_new_rank4_clade.csv")

# None differentiated abundant
rhabditida_new_rank5_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank5)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank5)$Clade,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank5_clade, "exported_tables/nema/rhabditida_new_rank5_clade.csv")


# None differentiated abundant
rhabditida_new_asv_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new)$Clade,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_asv_clade, "exported_tables/nema/rhabditida_new_asv_clade.csv")

# Across habitats and taxonomic ranks
# ASV00081 - Firmicutes (p=0.058612493)
rhabditida_new_rank1_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank1)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank1)$Sample.Site,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank1_site, "exported_tables/nema/rhabditida_new_rank1_site.csv")

# None
rhabditida_new_rank2_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank2)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank2)$Sample.Site,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank2_site, "exported_tables/nema/rhabditida_new_rank2_site.csv")

# None
rhabditida_new_rank3_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank3)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank3)$Sample.Site,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank3_site, "exported_tables/nema/rhabditida_new_rank3_site.csv")

# None
rhabditida_new_rank4_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank4)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank4)$Sample.Site,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank4_site, "exported_tables/nema/rhabditida_new_rank4_site.csv")

# None
rhabditida_new_rank5_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank5)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank5)$Sample.Site,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank5_site, "exported_tables/nema/rhabditida_new_rank5_site.csv")


# None
rhabditida_new_asv_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new)$Sample.Site,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_asv_site, "exported_tables/nema/rhabditida_new_asv_site.csv")


# Across habitats and taxonomic ranks
# ASV00081 - Firmicutes (p=0.058612493)
rhabditida_new_rank1_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank1)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank1)$Habitat,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank1_habitat, "exported_tables/nema/rhabditida_new_rank1_habitat.csv")

# None
rhabditida_new_rank2_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank2)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank2)$Habitat,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank2_habitat, "exported_tables/nema/rhabditida_new_rank2_habitat.csv")

# None
rhabditida_new_rank3_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank3)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank3)$Habitat,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank3_habitat, "exported_tables/nema/rhabditida_new_rank3_habitat.csv")

# None
rhabditida_new_rank4_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank4)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank4)$Habitat,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank4_habitat, "exported_tables/nema/rhabditida_new_rank4_habitat.csv")

# None
rhabditida_new_rank5_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank5)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank5)$Habitat,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank5_habitat, "exported_tables/nema/rhabditida_new_rank5_habitat.csv")


# None
rhabditida_new_asv_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new)$Habitat,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_asv_habitat, "exported_tables/nema/rhabditida_new_asv_habitat.csv")

# Between soil type and taxonomic ranks
# None differentiated abundant
rhabditida_new_rank1_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank1)),
                                                  phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank1)$Soil.Type,
                                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank1_soil_type, "exported_tables/nema/rhabditida_new_rank1_soil_type.csv")

# ASV00156 - Bacteroidota (p= 0.053166304)
rhabditida_new_rank2_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank2)),
                                                  phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank2)$Soil.Type,
                                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank2_soil_type, "exported_tables/nema/rhabditida_new_rank2_soil_type.csv")

# ASV00004 - Staphylococcales (p = 0.05716674)
rhabditida_new_rank3_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank3)),
                                                  phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank3)$Soil.Type,
                                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank3_soil_type, "exported_tables/nema/rhabditida_new_rank3_soil_type.csv")

# None
rhabditida_new_rank4_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank4)),
                                                  phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank4)$Soil.Type,
                                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank4_soil_type, "exported_tables/nema/rhabditida_new_rank4_soil_type.csv")

# None
rhabditida_new_rank5_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank5)),
                                                  phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new_rank5)$Soil.Type,
                                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_rank5_soil_type, "exported_tables/nema/rhabditida_new_rank5_soil_type.csv")

# None
rhabditida_new_asv_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_rhabditida_new)$Soil.Type,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(rhabditida_new_asv_soil_type, "exported_tables/nema/rhabditida_new_asv_soil_type.csv")

######################################## Panagrolaimida nematodes and taxonomic ranks ######################################## 
# Across Panagrolaimida clades and taxonomic ranks
# None differentiated abundant
panagrolaimida_new_rank1_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank1)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank1)$Clade,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_rank1_clade, "exported_tables/nema/panagrolaimida_new_rank1_clade.csv")

# None differentiated abundant
panagrolaimida_new_rank2_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank2)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank2)$Clade,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_rank2_clade, "exported_tables/nema/panagrolaimida_new_rank2_clade.csv")

# ASV00006 - Corynebacteriales
panagrolaimida_new_rank3_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank3)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank3)$Clade,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_rank3_clade, "exported_tables/nema/panagrolaimida_new_rank3_clade.csv")

# None differentiated abundant
panagrolaimida_new_rank4_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank4)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank4)$Clade,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_rank4_clade, "exported_tables/nema/panagrolaimida_new_rank4_clade.csv")

# None differentiated abundant
panagrolaimida_new_rank5_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank5)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank5)$Clade,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_rank5_clade, "exported_tables/nema/panagrolaimida_new_rank5_clade.csv")


panagrolaimida_new_asv_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new)$Clade,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_asv_clade, "exported_tables/nema/panagrolaimida_new_asv_clade.csv")

# Across sites and taxonomic ranks
# None
panagrolaimida_new_rank1_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank1)),
                                                  phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank1)$Sample.Site,
                                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_rank1_site, "exported_tables/nema/panagrolaimida_new_rank1_site.csv")

# ASV00006 - Actinobacteria
panagrolaimida_new_rank2_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank2)),
                                                  phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank2)$Sample.Site,
                                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_rank2_site, "exported_tables/nema/panagrolaimida_new_rank2_site.csv")

# None
panagrolaimida_new_rank3_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank3)),
                                                  phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank3)$Sample.Site,
                                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_rank3_site, "exported_tables/nema/panagrolaimida_new_rank3_site.csv")

# None
panagrolaimida_new_rank4_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank4)),
                                                  phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank4)$Sample.Site,
                                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_rank4_site, "exported_tables/nema/panagrolaimida_new_rank4_site.csv")

# None
panagrolaimida_new_rank5_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank5)),
                                                  phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank5)$Sample.Site,
                                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_rank5_site, "exported_tables/nema/panagrolaimida_new_rank5_site.csv")

# None
panagrolaimida_new_asv_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new)),
                                                phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new)$Sample.Site,
                                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_asv_site, "exported_tables/nema/panagrolaimida_new_asv_site.csv")


# Across habitats and taxonomic ranks
# None
panagrolaimida_new_rank1_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank1)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank1)$Habitat,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_rank1_habitat, "exported_tables/nema/panagrolaimida_new_rank1_habitat.csv")

# ASV00006 - Actinobacteria
panagrolaimida_new_rank2_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank2)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank2)$Habitat,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_rank2_habitat, "exported_tables/nema/panagrolaimida_new_rank2_habitat.csv")

# None
panagrolaimida_new_rank3_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank3)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank3)$Habitat,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_rank3_habitat, "exported_tables/nema/panagrolaimida_new_rank3_habitat.csv")

# None
panagrolaimida_new_rank4_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank4)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank4)$Habitat,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_rank4_habitat, "exported_tables/nema/panagrolaimida_new_rank4_habitat.csv")

# None
panagrolaimida_new_rank5_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank5)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new_rank5)$Habitat,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_rank5_habitat, "exported_tables/nema/panagrolaimida_new_rank5_habitat.csv")

# None
panagrolaimida_new_asv_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new)),
                                                  phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_panagrolaimida_new)$Habitat,
                                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(panagrolaimida_new_asv_habitat, "exported_tables/nema/panagrolaimida_new_asv_habitat.csv")

######################################## Dorylaimida nematodes and taxonomic ranks ######################################## 
# none taxa differentially abundant
dorylaimida_new_rank1_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank1)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank1)$Clade,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank1_clade, "exported_tables/nema/dorylaimida_new_rank1_clade.csv")

# none taxa differentially abundant
dorylaimida_new_rank2_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank2)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank2)$Clade,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank2_clade, "exported_tables/nema/dorylaimida_new_rank2_clade.csv")

# none taxa differentially abundant
dorylaimida_new_rank3_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank3)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank3)$Clade,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank3_clade, "exported_tables/nema/dorylaimida_new_rank3_clade.csv")

# none taxa differentially abundant
dorylaimida_new_rank4_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank4)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank4)$Clade,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank4_clade, "exported_tables/nema/dorylaimida_new_rank4_clade.csv")

# none taxa differentially abundant
dorylaimida_new_rank5_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank5)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank5)$Clade,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank5_clade, "exported_tables/nema/dorylaimida_new_rank5_clade.csv")

# none taxa differentially abundant
dorylaimida_new_asv_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new)$Clade,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_asv_clade, "exported_tables/nema/dorylaimida_new_asv_clade.csv")

### Across genera
dorylaimida_new_rank1_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank1)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank1)$Genus,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank1_genus, "exported_tables/nema/dorylaimida_new_rank1_genus.csv")

dorylaimida_new_rank2_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank2)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank2)$Genus,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank2_genus, "exported_tables/nema/dorylaimida_new_rank2_genus.csv")

dorylaimida_new_rank3_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank3)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank3)$Genus,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank3_genus, "exported_tables/nema/dorylaimida_new_rank3_genus.csv")

dorylaimida_new_rank4_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank4)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank4)$Genus,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank4_genus, "exported_tables/nema/dorylaimida_new_rank4_genus.csv")

dorylaimida_new_rank5_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank5)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank5)$Genus,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank5_genus, "exported_tables/nema/dorylaimida_new_rank5_genus.csv")

dorylaimida_new_asv_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new)$Genus,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_asv_genus, "exported_tables/nema/dorylaimida_new_asv_genus.csv")


### Across sites and taxonomic ranks
dorylaimida_new_rank1_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank1)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank1)$Sample.Site,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank1_site, "exported_tables/nema/dorylaimida_new_rank1_site.csv")

dorylaimida_new_rank2_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank2)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank2)$Sample.Site,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank2_site, "exported_tables/nema/dorylaimida_new_rank2_site.csv")

dorylaimida_new_rank3_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank3)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank3)$Sample.Site,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank3_site, "exported_tables/nema/dorylaimida_new_rank3_site.csv")

dorylaimida_new_rank4_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank4)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank4)$Sample.Site,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank4_site, "exported_tables/nema/dorylaimida_new_rank4_site.csv")

dorylaimida_new_rank5_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank5)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank5)$Sample.Site,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank5_site, "exported_tables/nema/dorylaimida_new_rank5_site.csv")

dorylaimida_new_asv_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new)),
                                           phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new)$Sample.Site,
                                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_asv_site, "exported_tables/nema/dorylaimida_new_asv_site.csv")

### Across habitats and taxonomic ranks
dorylaimida_new_rank1_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank1)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank1)$Habitat,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank1_habitat, "exported_tables/nema/dorylaimida_new_rank1_habitat.csv")

dorylaimida_new_rank2_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank2)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank2)$Habitat,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank2_habitat, "exported_tables/nema/dorylaimida_new_rank2_habitat.csv")

dorylaimida_new_rank3_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank3)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank3)$Habitat,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank3_habitat, "exported_tables/nema/dorylaimida_new_rank3_habitat.csv")

dorylaimida_new_rank4_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank4)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank4)$Habitat,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank4_habitat, "exported_tables/nema/dorylaimida_new_rank4_habitat.csv")

dorylaimida_new_rank5_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank5)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank5)$Habitat,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank5_habitat, "exported_tables/nema/dorylaimida_new_rank5_habitat.csv")

dorylaimida_new_asv_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new)$Habitat,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_asv_habitat, "exported_tables/nema/dorylaimida_new_asv_habitat.csv")

### Across soil type and taxonomic ranks
dorylaimida_new_rank1_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank1)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank1)$Soil.Type,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank1_soil_type, "exported_tables/nema/dorylaimida_new_rank1_soil_type.csv")

dorylaimida_new_rank2_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank2)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank2)$Soil.Type,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank2_soil_type, "exported_tables/nema/dorylaimida_new_rank2_soil_type.csv")

dorylaimida_new_rank3_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank3)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank3)$Soil.Type,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank3_soil_type, "exported_tables/nema/dorylaimida_new_rank3_soil_type.csv")

dorylaimida_new_rank4_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank4)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank4)$Soil.Type,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank4_soil_type, "exported_tables/nema/dorylaimida_new_rank4_soil_type.csv")

dorylaimida_new_rank5_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank5)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new_rank5)$Soil.Type,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_rank5_soil_type, "exported_tables/nema/dorylaimida_new_rank5_soil_type.csv")

dorylaimida_new_asv_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_dorylaimida_new)$Soil.Type,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(dorylaimida_new_asv_soil_type, "exported_tables/nema/dorylaimida_new_asv_soil_type.csv")

######################################## Tylenchida nematodes and taxonomic ranks ######################################## 
# Clades and taxonomic ranks
# none taxa differentially abundant
tylenchida_new_rank1_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank1)),
                                             phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank1)$Clade,
                                             mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank1_clade, "exported_tables/nema/tylenchida_new_rank1_clade.csv")

# none taxa differentially abundant
tylenchida_new_rank2_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank2)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank2)$Clade,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank2_clade, "exported_tables/nema/tylenchida_new_rank2_clade.csv")

# none taxa differentially abundant
tylenchida_new_rank3_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank3)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank3)$Clade,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank3_clade, "exported_tables/nema/tylenchida_new_rank3_clade.csv")

# none taxa differentially abundant
tylenchida_new_rank4_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank4)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank4)$Clade,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank4_clade, "exported_tables/nema/tylenchida_new_rank4_clade.csv")

# none taxa differentially abundant
tylenchida_new_rank5_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank5)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank5)$Clade,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank5_clade, "exported_tables/nema/tylenchida_new_rank5_clade.csv")

# none taxa differentially abundant
tylenchida_new_asv_clade <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new)$Clade,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_asv_clade, "exported_tables/nema/tylenchida_new_asv_clade.csv")

# Genera and taxonomic ranks
# none taxa differentially abundant
tylenchida_new_rank1_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank1)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank1)$Genus,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank1_genus, "exported_tables/nema/tylenchida_new_rank1_genus.csv")

# none taxa differentially abundant
tylenchida_new_rank2_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank2)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank2)$Genus,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank2_genus, "exported_tables/nema/tylenchida_new_rank2_genus.csv")

# none taxa differentially abundant
tylenchida_new_rank3_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank3)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank3)$Genus,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank3_genus, "exported_tables/nema/tylenchida_new_rank3_genus.csv")

# none taxa differentially abundant
tylenchida_new_rank4_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank4)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank4)$Genus,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank4_genus, "exported_tables/nema/tylenchida_new_rank4_genus.csv")

# none taxa differentially abundant
tylenchida_new_rank5_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank5)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank5)$Genus,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank5_genus, "exported_tables/nema/tylenchida_new_rank5_genus.csv")

# none taxa differentially abundant
tylenchida_new_asv_genus <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new)$Genus,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_asv_genus, "exported_tables/nema/tylenchida_new_asv_genus.csv")

# Sites and taxonomic ranks
# none taxa differentially abundant
tylenchida_new_rank1_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank1)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank1)$Sample.Site,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank1_site, "exported_tables/nema/tylenchida_new_rank1_site.csv")

# none taxa differentially abundant
tylenchida_new_rank2_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank2)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank2)$Sample.Site,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank2_site, "exported_tables/nema/tylenchida_new_rank2_site.csv")

# none taxa differentially abundant
tylenchida_new_rank3_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank3)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank3)$Sample.Site,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank3_site, "exported_tables/nema/tylenchida_new_rank3_site.csv")

# none taxa differentially abundant
tylenchida_new_rank4_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank4)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank4)$Sample.Site,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank4_site, "exported_tables/nema/tylenchida_new_rank4_site.csv")

# none taxa differentially abundant
tylenchida_new_rank5_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank5)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank5)$Sample.Site,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank5_site, "exported_tables/nema/tylenchida_new_rank5_site.csv")

# none taxa differentially abundant
tylenchida_new_asv_site <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new)$Sample.Site,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_asv_site, "exported_tables/nema/tylenchida_new_asv_site.csv")

# Habitats and taxonomic ranks
# none taxa differentially abundant
tylenchida_new_rank1_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank1)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank1)$Habitat,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank1_habitat, "exported_tables/nema/tylenchida_new_rank1_habitat.csv")

# none taxa differentially abundant
tylenchida_new_rank2_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank2)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank2)$Habitat,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank2_habitat, "exported_tables/nema/tylenchida_new_rank2_habitat.csv")

# none taxa differentially abundant
tylenchida_new_rank3_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank3)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank3)$Habitat,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank3_habitat, "exported_tables/nema/tylenchida_new_rank3_habitat.csv")

# none taxa differentially abundant
tylenchida_new_rank4_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank4)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank4)$Habitat,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank4_habitat, "exported_tables/nema/tylenchida_new_rank4_habitat.csv")

# none taxa differentially abundant
tylenchida_new_rank5_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank5)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank5)$Habitat,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank5_habitat, "exported_tables/nema/tylenchida_new_rank5_habitat.csv")

# none taxa differentially abundant
tylenchida_new_asv_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new)),
                                          phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new)$Habitat,
                                          mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_asv_habitat, "exported_tables/nema/tylenchida_new_asv_habitat.csv")

# Soil type and taxonomic ranks
# none taxa differentially abundant
tylenchida_new_rank1_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank1)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank1)$Soil.Type,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank1_soil_type, "exported_tables/nema/tylenchida_new_rank1_soil_type.csv")

# none taxa differentially abundant
tylenchida_new_rank2_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank2)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank2)$Soil.Type,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank2_soil_type, "exported_tables/nema/tylenchida_new_rank2_soil_type.csv")

# none taxa differentially abundant
tylenchida_new_rank3_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank3)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank3)$Soil.Type,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank3_soil_type, "exported_tables/nema/tylenchida_new_rank3_soil_type.csv")

# none taxa differentially abundant
tylenchida_new_rank4_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank4)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank4)$Soil.Type,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank4_soil_type, "exported_tables/nema/tylenchida_new_rank4_soil_type.csv")

# none taxa differentially abundant
tylenchida_new_rank5_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank5)),
                                              phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new_rank5)$Soil.Type,
                                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_rank5_soil_type, "exported_tables/nema/tylenchida_new_rank5_soil_type.csv")

# none taxa differentially abundant
tylenchida_new_asv_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new)),
                                            phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_tylenchida_new)$Soil.Type,
                                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
write.csv(tylenchida_new_asv_soil_type, "exported_tables/nema/tylenchida_new_asv_soil_type.csv")

######################################## Ordinations per nematode group ######################################## 
# Plot PCA based on clr transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_aphelenchids_clr,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_aphelenchids_clr_pca,
  type="samples",
  color="Sample.Site",
  shape = "Habitat") +
  facet_nested(~ Soil.Type, scales = "free") +
  #geom_point(aes(shape = Sample.Site, fill = Sample.Site), size = 3, alpha = 0.7) +
  ggrepel::geom_text_repel(aes(label = Genus.Abbr), size = 2, color = "black", max.overlaps = 30) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A",
                                         "#984EA3", "#FF7F00", "lightblue",
                                         "#A65628", "#F781BF", "#999999",
                                         "#1B9E77", "#D95F02", "#7570B3",
                                         "#E7298A", "#66A61E", "#E6AB02",
                                         "#A6761D", "#666666", "lightblue")) +
                                           scale_shape_manual(values = c(15, 19, 3, 8, 4, 17,
                                                                         1, 2, 3, 5, 25, 20),
                                                              name = "Site",
                                                              breaks = c("SK.01", "SK.02", "SK.03",
                                                                         "SK.05", "SK.06", "SK.07",
                                                                         "SK.09", "SK.10", "SK.11",
                                                                         "SK.16", "SK.17", "SK.18"),
                                                              labels = c("01", "02", "03",
                                                                         "05", "06", "07",
                                                                         "09", "10", "11",
                                                                         "16", "17", "18"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #geom_point(colour = "grey90", size = 1.5) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend right   
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5, legend.position = "right") # adjusts the title of the legend
#ylab("NMDS2") + # add the title on y axis
#xlab("NMDS1") + # add the title on x axis
#ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
#theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
#theme(strip.background =element_blank()) + # remove the background of titles
#annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
#annotate("text", x = 3, y = 2, label ="2D Stress: 0.07") +
#coord_fixed(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca_clr2 / phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca_clr1) +
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_aphelenchida_clr_pca.pdf", width = 6, height = 6, dpi = 200)



phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp_sel <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp,
                                                                                 sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.216" & # Aphelenchus sp1, SK03
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.112" & # Aporcelaimellus.sp2, SK02
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.73" & # Tylenchorhynchus, SK12
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.42" & # unknown Aphelenchida, SK08
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.57" & # Aphelenchus, SK08
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.501" & # Plectus.sp1, SK11
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.98" & # Chiloplacus.sp2, SK02
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.118" & # Chiloplacus.sp1, SK02
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.28" & # Aphelenchoides.sp1, SK16
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.551" & # Aphelenchus.sp1, SK15
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.524" & # Mesorhabditis, SK17
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.288" & # Tylenchus, SK13
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.474" & # Tylenchorhynchus, SK09
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.289" & # Leptonchus, SK13
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.277")  # Tylocephalus, SK10

# Ordinate NMDS
set.seed(1)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp_sel_nmds <- ordinate(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp_sel, 
  method = "NMDS", 
  distance = "bray"
)

theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp_sel,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp_sel_nmds,
  label = "Description",
  type = "samples",
  color = "Genus",
  shape = "Feeding.Group") +
  facet_nested(~ Soil.Type+Habitat, scales = "free") +
  geom_point(aes(color = Genus), alpha = 0.7, size = 3) +
  scale_color_manual(values = c("#52392F", "#C29B6C", "#83643E",
                                         "#397A4C", "#77C063", "#BEDB92"),
                                         name = "Habitat",
                     breaks = c("Chaparral", "Coastal.scrub", "Native.grass", "Hollyleaf.cherry", "Oak.wood", "Riparian"),
                     labels = c("Chaparral", "Coastal Sage Scrub", "Native Grass", "Hollyleaf Cherry", "Oak Woodland", "Riparian"))+
  scale_shape_manual(values = c(19,17),
                     name = "Soil Type",
                     breaks = c("Compact.Soil", "Not.Compact.Soil"),
                     labels = c("Compacted", "Uncompacted"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  #geom_point(colour = "grey90", size = 1.5) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("NMDS2") + # add the title on y axis
  #xlab("NMDS1") + # add the title on x axis
  #ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  #annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
  annotate("text", x = 1.4, y = 0.8, label ="2D Stress: 0.07") +
  ggsave("16S_nem_ShSk_ASV_NMDS_samples_greater_500_prev05_no-mito-chloro.tiff", width = 6, height = 4, dpi = 150)


# CAP analysis
set.seed(1)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp_sel_cap <- ordinate(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp_sel, 
  method = "CAP", 
  distance = "bray",
  formula = ~ Soil.Type + Habitat + Sample.Site
)

theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp_sel,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp_sel_cap,
  #label = "Description",
  type = "samples",
  color = "Genus",
  shape = "Feeding.Group") +
  facet_nested(~ Soil.Type+Habitat+Sample.Site, scales = "free") +
  geom_point(aes(color = Genus), alpha = 0.7, size = 3) +
  theme(legend.position = "bottom")
scale_shape_manual(values = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                   name = "Nematode Order") +
  scale_color_manual(values = c("#52392F", "#C29B6C", "#83643E",
                                         "#397A4C", "#77C063", "#BEDB92"),
                                         name = "Habitat",
                     breaks = c("Chaparral", "Coastal.scrub", "Native.grass", "Hollyleaf.cherry", "Oak.wood", "Riparian"),
                     labels = c("Chaparral", "Coastal Sage Scrub", "Native Grass", "Hollyleaf Cherry", "Oak Woodland", "Riparian"))+
  scale_shape_manual(values = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                     name = "Nematode Order") +
  #breaks = c("Compact.Soil", "Not.Compact.Soil"),
  #labels = c("Compacted", "Uncompacted"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  #geom_point(colour = "grey90", size = 1.5) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("NMDS2") + # add the title on y axis
  #xlab("NMDS1") + # add the title on x axis
  #ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  #annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
  annotate("text", x = 1.4, y = 0.8, label ="2D Stress: 0.07") +
  ggsave("16S_nem_ShSk_ASV_NMDS_samples_greater_500_prev05_no-mito-chloro.tiff", width = 6, height = 4, dpi = 150)

# Prune phyloseq object according to soil type, compact vs. uncompact
phy_16S_prune_nem_clr_compact <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr,
                                                Soil.Type == "Compact.Soil" &
                                                  SSU.nema.qual != "Bad") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_16S_prune_nem_clr_compact
unique(sample_data(phy_16S_prune_nem_clr_compact)$SSU.nema.qual)
any(taxa_sums(phy_16S_prune_nem_clr_compact) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16S_prune_nem_clr_compact) == 0) # gives the number of cases
any(sample_sums(phy_16S_prune_nem_clr_compact) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(sample_sums(phy_16S_prune_nem_clr_compact) == 0) # gives the number of cases
smin <- min(sample_sums(phy_16S_prune_nem_clr_compact))

phy_16S_prune_nem_clr_compact_pca <- phyloseq::ordinate(
  physeq = phy_16S_prune_nem_clr_compact,
  method = "RDA",
  distance = "euclidean"
)

# Plot scree plot to see PCA contributions
phyloseq::plot_scree(phy_16S_prune_nem_clr_compact_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/16S_Nema_ShSk_ASV_noncontam_prev05_true_all_filt_clr_pca.tiff", width = 12, height = 6, dpi = 150)

head(phy_16S_prune_nem_clr_compact_pca$CA$eig)
sapply(phy_16S_prune_nem_clr_compact_pca$CA$eig[1:5], function(x) x / sum(phy_16S_prune_nem_clr_compact_pca$CA$eig))

# Scale axes and plot ordination
phy_16S_prune_nem_clr_compact_pca_clr1 <- phy_16S_prune_nem_clr_compact_pca$CA$eig[1] / sum(phy_16S_prune_nem_clr_compact_pca$CA$eig)
phy_16S_prune_nem_clr_compact_pca_clr2 <- phy_16S_prune_nem_clr_compact_pca$CA$eig[2] / sum(phy_16S_prune_nem_clr_compact_pca$CA$eig)

# Plot PCA based on clr transformation
phyloseq::plot_ordination(
  physeq = phy_16S_prune_nem_clr_compact,
  ordination = phy_16S_prune_nem_clr_compact_pca,
  type="samples",
  color="Sample.Site",
  shape = "Sample.Site") +
  facet_nested(~ Habitat, scales = "free") +
  geom_point(aes(color = Sample.Site), alpha = 0.7, size = 3) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A",
                                         "#984EA3", "#FF7F00", "#FFFF33",
                                         "#A65628", "#F781BF", "#999999"),
                                         name = "Site",
                     breaks = c("SK.01", "SK.02", "SK.03",
                                "SK.13", "SK.14", "SK.15",
                                "SK.16", "SK.17", "SK.18"),
                     labels = c("01", "02", "03",
                                "13", "14", "15",
                                "16", "17", "18"))+
  theme(legend.position = "bottom") +
  scale_shape_manual(values = c(15, 24, 3, 4, 17, 6, 25, 0, 13)) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #geom_point(colour = "grey90", size = 1.5) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("NMDS2") + # add the title on y axis
  #xlab("NMDS1") + # add the title on x axis
  #ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14))
#theme(strip.background =element_blank()) + # remove the background of titles
#annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
#annotate("text", x = 3, y = 2, label ="2D Stress: 0.07") +
#coord_fixed(phy_16_prune_soil_noncontam_prev05_true_all_filt_clr_pca_clr2 / phy_16_prune_soil_noncontam_prev05_true_all_filt_clr_pca_clr1) +
ggsave("results/16S_nema_ShSk_ASV_noncontam_prev05_true_all_filt_clr_no-mito-chloro_pca_plot_genus_sample_site_compacted.tiff", width = 12, height = 6, dpi = 150)

# Prune phyloseq object according to soil type, compact vs. uncompact
phy_16S_prune_nem_clr_uncompact <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr,
                                                  Soil.Type == "Not.Compact.Soil" &
                                                    SSU.nema.qual != "Bad") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_16S_prune_nem_clr_uncompact
unique(sample_data(phy_16S_prune_nem_clr_uncompact)$SSU.nema.qual)

any(taxa_sums(phy_16S_prune_nem_clr_uncompact) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16S_prune_nem_clr_uncompact) == 0) # gives the number of cases
any(sample_sums(phy_16S_prune_nem_clr_uncompact) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(sample_sums(phy_16S_prune_nem_clr_uncompact) == 0) # gives the number of cases
smin <- min(sample_sums(phy_16S_prune_nem_clr_uncompact))

phy_16S_prune_nem_clr_uncompact_pca <- phyloseq::ordinate(
  physeq = phy_16S_prune_nem_clr_uncompact,
  method = "RDA",
  distance = "euclidean"
)

# Plot scree plot to see PCA contributions
phyloseq::plot_scree(phy_16S_prune_nem_clr_uncompact_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/16S_Nema_ShSk_ASV_noncontam_prev05_true_all_filt_clr_pca_uncompact.tiff", width = 12, height = 6, dpi = 150)

head(phy_16S_prune_nem_clr_uncompact_pca$CA$eig)
sapply(phy_16S_prune_nem_clr_uncompact_pca$CA$eig[1:5], function(x) x / sum(phy_16S_prune_nem_clr_uncompact_pca$CA$eig))

# Scale axes and plot ordination
phy_16S_prune_nem_clr_uncompact_pca_clr1 <- phy_16S_prune_nem_clr_uncompact_pca$CA$eig[1] / sum(phy_16S_prune_nem_clr_uncompact_pca$CA$eig)
phy_16S_prune_nem_clr_uncompact_pca_clr2 <- phy_16S_prune_nem_clr_uncompact_pca$CA$eig[2] / sum(phy_16S_prune_nem_clr_uncompact_pca$CA$eig)

# Plot PCA based on clr transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16S_prune_nem_clr_uncompact,
  ordination = phy_16S_prune_nem_clr_uncompact_pca,
  type="samples",
  color="Sample.Site",
  shape = "Sample.Site") +
  facet_nested(~ Habitat, scales = "free") +
  geom_point(aes(color = Sample.Site), alpha = 0.7, size = 3) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3",
                                         "#E7298A", "#66A61E", "#E6AB02",
                                         "#A6761D", "#666666", "lightblue"),
                                         name = "Site",
                     breaks = c("SK.04", "SK.05", "SK.06",
                                "SK.07", "SK.08", "SK.09",
                                "SK.10", "SK.11", "SK.12"),
                     labels = c("04", "05", "06",
                                "07", "08", "09",
                                "10", "11", "12"))+
  scale_shape_manual(values = c(15, 24, 3, 4, 17, 6, 25, 0, 13)) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #geom_point(colour = "grey90", size = 1.5) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("NMDS2") + # add the title on y axis
  #xlab("NMDS1") + # add the title on x axis
  #ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14))
#theme(strip.background =element_blank()) + # remove the background of titles
#annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
#annotate("text", x = 3, y = 2, label ="2D Stress: 0.07") +
#coord_fixed(phy_16_prune_soil_noncontam_prev05_true_all_filt_clr_pca_clr2 / phy_16_prune_soil_noncontam_prev05_true_all_filt_clr_pca_clr1) +
ggsave("results/16S_nema_ShSk_ASV_noncontam_prev05_true_all_filt_clr_no-mito-chloro_pca_plot_genus_sample_site_uncompacted.tiff", width = 12, height = 6, dpi = 150)

phy_16S_prune_nem_uncompact_sel <- subset_samples(phy_16S_prune_nem_uncompact,
                                                  sample_names(phy_16S_prune_nem_uncompact) != "Nem.216" & # Aphelenchus sp1, SK03
                                                    sample_names(phy_16S_prune_nem_uncompact) != "Nem.112" & # Aporcelaimellus.sp2, SK02
                                                    sample_names(phy_16S_prune_nem_uncompact) != "Nem.98" & # Chiloplacus.sp2, SK02
                                                    sample_names(phy_16S_prune_nem_uncompact) != "Nem.28" & # Aphelenchoides.sp1, SK16
                                                    sample_names(phy_16S_prune_nem_uncompact) != "Nem.118" & # Chiloplacus.sp1, SK02
                                                    sample_names(phy_16S_prune_nem_uncompact) != "Nem.289" & # Leptonchus, SK13
                                                    sample_names(phy_16S_prune_nem_uncompact) != "Nem.206") # Acrobeloides.sp3, SK03

# Ordinate NMDS
set.seed(1)
phy_16S_prune_nem_uncompact_sel_nmds <- ordinate(
  physeq = phy_16S_prune_nem_uncompact_sel, 
  method = "NMDS", 
  distance = "bray"
)

theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16S_prune_nem_uncompact_sel,
  ordination = phy_16S_prune_nem_uncompact_sel_nmds,
  #label = "Description",
  type = "samples",
  color = "Genus",
  shape = "Feeding.Group") +
  facet_nested(~ Habitat, scales = "free") +
  #scale_color_manual(values = c("#52392F", "#C29B6C", "#83643E",
  #    "#397A4C", "#77C063", "#BEDB92"),
  #                name = "Habitat",
  #                 breaks = c("Chaparral", "Coastal.scrub", "Native.grass", "Hollyleaf.cherry", "Oak.wood", "Riparian"),
  #                labels = c("Chaparral", "Coastal Sage Scrub", "Native Grass", "Hollyleaf Cherry", "Oak Woodland", "Riparian"))+
  #scale_shape_manual(values = c(19,17),
  #                  name = "Soil Type",
  #                 breaks = c("Compact.Soil", "Not.Compact.Soil"),
  #                labels = c("Compacted", "Uncompacted"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Site), alpha = 0.7, size = 4) +
  #geom_point(colour = "grey90", size = 1.5) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("NMDS2") + # add the title on y axis
  #xlab("NMDS1") + # add the title on x axis
  #ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  #annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
  #annotate("text", x = 1.4, y = 0.8, label ="2D Stress: 0.07") +
  ggsave("16S_nem_ShSk_ASV_NMDS_samples_greater_500_prev05_no-mito-chloro.tiff", width = 6, height = 4, dpi = 150)
phy_16S_prune_nem_uncompact <- subset_samples(phy_16S_prune, Soil.Type == "Not.Compact.Soil")
phy_16S_prune_nem_uncompact


phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp_sel <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp,
                                                                                 sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.216" & # Aphelenchus sp1, SK03
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.112" & # Aporcelaimellus.sp2, SK02
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.73" & # Tylenchorhynchus, SK12
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.277" & # Tylocephalus, SK10
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.57" & # Aphelenchus, SK08
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.42" & # Aphelenchida unknown, SK08
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.98" & # Chiloplacus.sp2, SK02
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.118" & # Chiloplacus.sp1, SK02
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.28" & # Aphelenchoides.sp1, SK16
                                                                                   sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp) != "Nem.501") # Plectus.sp1, SK11


# Ordinate NMDS
set.seed(1)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp_sel_nmds <- ordinate(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp_sel, 
  method = "NMDS", 
  distance = "bray"
)

# Plot NMDS
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp_sel,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp_sel_nmds,
  label = "Description",
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#a65628", "red", "#ffae19",
                                         "#4daf4a", "#1919ff", "darkorchid3"),
                                         name = "Habitat",
                     breaks = c("Chaparral", "Coastal.scrub", "Hollyleaf.cherry", "Native.grass", "Oak.wood", "Riparian"),
                     labels = c("Chaparral", "Coastal Scrub", "Hollyleaf Cherry", "Native Grass", "Oak Wood", "Riparian"))+
  scale_shape_manual(values = c(19,17),
                     name = "Soil Type",
                     breaks = c("Compact.Soil", "Not.Compact.Soil"),
                     labels = c("Compacted", "Uncompacted"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("NMDS2") + # add the title on y axis
  #xlab("NMDS1") + # add the title on x axis
  #ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
  theme(plot.title = element_text(hjust = 1, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  annotate("text", x = 0.004, y = 0.005, label ="2D Stress: 0.28") +
  ggsave("16S_worm_ShSk_ASV_NMDS_greater_500_filt_comb05_no-mito-chloro_selected_samples_habitat.tiff", width = 8, height = 6, dpi = 150)



# import aldex2 results and rename the X variable by OTU
aldex2_phy_hab_result <- read.csv("aldex2_phy_hab.csv")
colnames(aldex2_phy_hab_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

aldex2_class_hab_result <- read.csv("aldex2_class_hab.csv")
colnames(aldex2_class_hab_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

aldex2_order_hab_result <- read.csv("aldex2_order_hab.csv")
colnames(aldex2_order_hab_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

aldex2_fam_hab_result <- read.csv("aldex2_fam_hab.csv")
colnames(aldex2_fam_hab_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

aldex2_genus_hab_result <- read.csv("aldex2_genus_hab.csv")
colnames(aldex2_genus_hab_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

aldex2_fam_soil_result <- read.csv("aldex2_fam_soil.csv")
colnames(aldex2_fam_soil_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

aldex2_genus_soil_result <- read.csv("aldex2_genus_soil.csv")
colnames(aldex2_genus_soil_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")


#Clean up presentation
taxa_info <- data.frame(tax_table(phy_16_prune_soil_noncontam_prev05_true_all_filt_scale))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

otu_count <- data.frame(otu_table(phy_16_prune_soil_noncontam_prev05_true_all_filt_scale))
otu_count <- otu_count %>% rownames_to_column(var = "OTU")

sample_tab_soil <- data.frame(sample_data(phy_16_prune_soil_noncontam_prev05_true_all_filt_scale))
sample_soil_factors <- data.frame(cbind(sample_tab_soil[, c(7, 10, 11, 25)]))

# filter aldex2 results by sig kw.eBH and join the taxanomic information
sig_aldex2_phy_hab_result <- aldex2_phy_hab_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_phy_hab_result <- left_join(sig_aldex2_phy_hab_result, taxa_info)
write.csv(sig_aldex2_phy_hab_result, "sig_aldex2_phy_hab_result.csv") # 8 taxa with significant differential abundance

sig_aldex2_class_hab_result <- aldex2_class_hab_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_class_hab_result <- left_join(sig_aldex2_class_hab_result, taxa_info)
write.csv(sig_aldex2_class_hab_result, "sig_aldex2_class_hab_result.csv") # 28 taxa with significant differential abundance

sig_aldex2_order_hab_result <- aldex2_order_hab_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_order_hab_result <- left_join(sig_aldex2_order_hab_result, taxa_info)
write.csv(sig_aldex2_order_hab_result, "sig_aldex2_order_hab_result.csv") # 75 taxa with significant differential abundance

sig_aldex2_order_hab_result_top20 <- sig_aldex2_order_hab_result %>% top_n(-20, kw.ep)
write.csv(sig_aldex2_order_hab_result_top20, "sig_aldex2_order_hab_result_top20.csv")

sig_aldex2_fam_hab_result <- aldex2_fam_hab_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_fam_hab_result <- left_join(sig_aldex2_fam_hab_result, taxa_info)
write.csv(sig_aldex2_fam_hab_result, "sig_aldex2_fam_hab_result.csv") # 74 taxa with significant differential abundance

sig_aldex2_fam_hab_result_top20 <- sig_aldex2_fam_hab_result %>% top_n(-20, kw.ep)
write.csv(sig_aldex2_fam_hab_result_top20, "sig_aldex2_fam_hab_result_top20.csv")

sig_aldex2_genus_hab_result <- aldex2_genus_hab_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_genus_hab_result <- left_join(sig_aldex2_genus_hab_result, taxa_info)
write.csv(sig_aldex2_genus_hab_result, "sig_aldex2_genus_hab_result.csv") # 196 taxa with significant differential abundance

sig_aldex2_genus_hab_result_top20 <- sig_aldex2_genus_hab_result %>% top_n(-20, kw.ep)
write.csv(sig_aldex2_genus_hab_result_top20, "sig_aldex2_genus_hab_result_top20.csv")

### Comparison among soil types only
sig_aldex2_fam_soil_result <- aldex2_fam_soil_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_fam_soil_result <- left_join(sig_aldex2_fam_soil_result, taxa_info)
write.csv(sig_aldex2_fam_soil_result, "sig_aldex2_fam_soil_result.csv") # 67 taxa with significant differential abundance

sig_aldex2_genus_soil_result <- aldex2_genus_soil_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_genus_soil_result <- left_join(sig_aldex2_genus_soil_result, taxa_info)
write.csv(sig_aldex2_genus_soil_result, "sig_aldex2_genus_soil_result.csv") # 177 taxa with significant differential abundance

# Create clr objects by using OTU ids and original otu table, all significant OTUs and only top 20
sig_aldex2_phy_hab_result_count <- left_join(sig_aldex2_phy_hab_result, otu_count)
clr_phy <- sig_aldex2_phy_hab_result_count[, -(2:10)]
rownames(clr_phy) <- clr_phy$OTU
clr_phy <- clr_phy[, -1]

sig_aldex2_class_hab_result_count <- left_join(sig_aldex2_class_hab_result, otu_count)
clr_class <- sig_aldex2_class_hab_result_count[, -(2:10)]
rownames(clr_class) <- clr_class$OTU
clr_class <- clr_class[, -1]

sig_aldex2_order_hab_result_count <- left_join(sig_aldex2_order_hab_result, otu_count)
clr_order <- sig_aldex2_order_hab_result_count[, -(2:10)]
rownames(clr_order) <- clr_order$OTU
clr_order <- clr_order[, -1]

sig_aldex2_order_hab_result_count_top20 <- left_join(sig_aldex2_order_hab_result_top20, otu_count)
clr_order_top20 <- sig_aldex2_order_hab_result_count_top20[, -(2:10)]
rownames(clr_order_top20) <- clr_order_top20$OTU
clr_order_top20 <- clr_order_top20[, -1]

sig_aldex2_fam_hab_result_count <- left_join(sig_aldex2_fam_hab_result, otu_count)
clr_fam <- sig_aldex2_fam_hab_result_count[, -(2:10)]
rownames(clr_fam) <- clr_fam$OTU
clr_fam <- clr_fam[, -1]

sig_aldex2_fam_hab_result_count_top20 <- left_join(sig_aldex2_fam_hab_result_top20, otu_count)
clr_fam_top20 <- sig_aldex2_fam_hab_result_count_top20[, -(2:10)]
rownames(clr_fam_top20) <- clr_fam_top20$OTU
clr_fam_top20 <- clr_fam_top20[, -1]

sig_aldex2_genus_hab_result_count <- left_join(sig_aldex2_genus_hab_result, otu_count)
clr_genus <- sig_aldex2_genus_hab_result_count[, -(2:10)]
rownames(clr_genus) <- clr_genus$OTU
clr_genus <- clr_genus[, -1]

sig_aldex2_genus_hab_result_count_top20 <- left_join(sig_aldex2_genus_hab_result_top20, otu_count)
clr_genus_top20 <- sig_aldex2_genus_hab_result_count_top20[, -(2:10)]
rownames(clr_genus_top20) <- clr_genus_top20$OTU
clr_genus_top20 <- clr_genus_top20[, -1]

sig_aldex2_genus_soil_result_count <- left_join(sig_aldex2_genus_soil_result, otu_count)
clr_genus_soil <- sig_aldex2_genus_soil_result_count[, -(2:10)]
rownames(clr_genus_soil) <- clr_genus_soil$OTU
clr_genus_soil <- clr_genus_soil[, -1]

# Adjusting zeros on the matrix
shsk_phylum_czm <- cmultRepl(t(clr_phy),  label=0, method="CZM")
shsk_phylum_clr <- t(apply(shsk_phylum_czm, 1, function(x){log(x) - mean(log(x))}))

shsk_class_czm <- cmultRepl(t(clr_class),  label=0, method="CZM")
shsk_class_clr <- t(apply(shsk_class_czm, 1, function(x){log(x) - mean(log(x))}))

shsk_order_czm <- cmultRepl(t(clr_order),  label=0, method="CZM")
shsk_order_clr <- t(apply(shsk_order_czm, 1, function(x){log(x) - mean(log(x))}))

shsk_order_czm_top20 <- cmultRepl(t(clr_order_top20),  label=0, method="CZM")
shsk_order_clr_top20 <- t(apply(shsk_order_czm_top20, 1, function(x){log(x) - mean(log(x))}))

shsk_fam_czm <- cmultRepl(t(clr_fam),  label=0, method="CZM")
shsk_fam_clr <- t(apply(shsk_fam_czm, 1, function(x){log(x) - mean(log(x))}))

shsk_fam_czm_top20 <- cmultRepl(t(clr_fam_top20),  label=0, method="CZM")
shsk_fam_clr_top20 <- t(apply(shsk_fam_czm_top20, 1, function(x){log(x) - mean(log(x))}))

shsk_genus_czm <- cmultRepl(t(clr_genus),  label=0, method="CZM")
shsk_genus_clr <- t(apply(shsk_genus_czm, 1, function(x){log(x) - mean(log(x))}))

shsk_genus_czm_top20 <- cmultRepl(t(clr_genus_top20),  label=0, method="CZM")
shsk_genus_clr_top20 <- t(apply(shsk_genus_czm_top20, 1, function(x){log(x) - mean(log(x))}))

shsk_genus_czm_soil <- cmultRepl(t(clr_genus_soil),  label=0, method="CZM")
shsk_genus_clr_soil <- t(apply(shsk_genus_czm_soil, 1, function(x){log(x) - mean(log(x))}))

# Checking first lines of object, transpose it, and then create a heatmap according to the tax rank
head(shsk_phylum_clr)
shsk_phylum_clr_trav <- t(shsk_phylum_clr)
shsk_phylum_clr[, order(colnames(shsk_phylum_clr))]
heatmap(shsk_phylum_clr_trav, scale = "none", col = bluered(100))

shsk_class_clr_trav <- t(shsk_class_clr)
shsk_class_clr[, order(colnames(shsk_class_clr))]
heatmap(shsk_class_clr_trav, scale = "none", col = bluered(100))

shsk_order_clr_trav <- t(shsk_order_clr)
shsk_order_clr[, order(colnames(shsk_order_clr))]
heatmap(shsk_order_clr_trav, scale = "none", col = bluered(100))

shsk_order_clr_trav_top20 <- t(shsk_order_clr_top20)
shsk_order_clr_top20[, order(colnames(shsk_order_clr_top20))]
heatmap(shsk_order_clr_trav_top20, scale = "none", col = bluered(100))

shsk_fam_clr_trav <- t(shsk_fam_clr)
shsk_fam_clr[, order(colnames(shsk_fam_clr))]
heatmap(shsk_fam_clr_trav, scale = "none", col = bluered(100))

shsk_fam_clr_trav_top20 <- t(shsk_fam_clr_top20)
shsk_fam_clr_top20[, order(colnames(shsk_fam_clr_top20))]
heatmap(shsk_fam_clr_trav_top20, scale = "none", col = bluered(100))

shsk_genus_clr_trav <- t(shsk_genus_clr)
shsk_genus_clr[, order(colnames(shsk_genus_clr))]
heatmap(shsk_genus_clr_trav, scale = "none", col = bluered(100))

shsk_genus_clr_trav_top20 <- t(shsk_genus_clr_top20)
shsk_genus_clr_top20[, order(colnames(shsk_genus_clr_top20))]
heatmap(shsk_genus_clr_trav_top20, scale = "none", col = bluered(100))

shsk_genus_clr_trav_soil <- t(shsk_genus_clr_soil)
shsk_genus_clr_soil[, order(colnames(shsk_genus_clr_soil))]
heatmap(shsk_genus_clr_trav_soil, scale = "none", col = bluered(100))

# Define palette color
col_matrix <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
col_matrix2 <- colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_matrix3 <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Combine the heatmap and the annotation
Z.Score.phylum <- scale(shsk_phylum_clr_trav)

Z.Score.class <- scale(shsk_class_clr_trav)

Z.Score.order <- scale(shsk_order_clr_trav)
Z.Score.order.top20 <- scale(shsk_order_clr_trav_top20)

Z.Score.fam <- scale(shsk_fam_clr_trav)
Z.Score.fam.top20 <- scale(shsk_fam_clr_trav_top20)

Z.Score.genus <- scale(shsk_genus_clr_trav)
Z.Score.genus.top20 <- scale(shsk_genus_clr_trav_top20)

# Define colors for each level of qualitative variables, i.e. soil type and habitat
# Create the heatmap annotation
ha_top = HeatmapAnnotation(
  Soil = as.vector(sample_soil_factors$Soil.Type),
  Habitat = as.vector(sample_soil_factors$Habitat),
  col = list(
    Soil = c("Compact.Soil" = "#674422", "Not.Compact.Soil" = "#004B00"),
    Habitat = c("Chaparral" = "#52392F", "Coastal.scrub" = "#83643E", "Native.grass" = "#C29B6C",
                "Hollyleaf.cherry" = "#397A4C", "Oak.wood" = "#77C063", "Riparian" = "#BEDB92")
  ),
  annotation_legend_param = list(
    Soil = list(
      title = "Soil Type",
      at = c("Compact.Soil", "Not.Compact.Soil"),
      labels = c ("Compacted", "Uncompacted")
    ),
    Habitat = list(
      title = "Habitat",
      at = c("Chaparral", "Coastal.scrub", "Native.grass", "Hollyleaf.cherry", "Oak.wood", "Riparian"),
      labels = c ("Chaparral", "Coastal Sage Scrub", "Native Grass", "Hollyleaf Cherry", "Oak Woodland", "Riparian")
    )
  ))

# Organize total abundance and taxa name
Z.Score.phylum_name <- as.data.frame(Z.Score.phylum)
str(Z.Score.phylum_name)
Z.Score.phylum_name <- rownames_to_column(Z.Score.phylum_name, var = "OTU")
Z.Score.phylum_name <- left_join(Z.Score.phylum_name, taxa_info)
Z.Score.phylum_name <- Z.Score.phylum_name[, -(2:51)]

otu_count_total <- as.data.frame(otu_count)
otu_count_total$Total <- rowSums(otu_count_total[, -1])
head(otu_count_total)
otu_count_total <- otu_count_total[, -(2:51)]

Z.Score.phylum_count_total <- left_join(Z.Score.phylum_name, otu_count_total)
head(Z.Score.phylum_count_total)

ha_right_phy = rowAnnotation(
  Abundance = Z.Score.phylum_count_total$Total)

row_labels_phy = Z.Score.phylum_count_total$Phylum

# Plot heatmap at the phylum level
hm_phylum <- Heatmap(Z.Score.phylum, name = "Z-score, CLR", col = col_matrix3,
                     column_title = "16S rRNA soil microbiome", 
                     column_title_gp = gpar(fontface = "bold", fontsize = 14),
                     top_annotation = ha_top,
                     right_annotation = ha_right_phy,
                     row_title = "Phylum",
                     row_labels = row_labels_phy,
                     row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                     row_names_gp = gpar(fontsize = 10),
                     row_order = order(row_labels_phy),
                     rect_gp = gpar(col = "white", lwd = 1),
                     show_column_names = TRUE,
                     show_heatmap_legend = TRUE)
hm_phylum

Z.Score.class_name <- as.data.frame(Z.Score.class)
str(Z.Score.class_name)
Z.Score.class_name <- rownames_to_column(Z.Score.class_name, var = "OTU")
Z.Score.class_name <- left_join(Z.Score.class_name, taxa_info)
Z.Score.class_name <- Z.Score.class_name[, -(2:51)]

Z.Score.class_count_total <- left_join(Z.Score.class_name, otu_count_total)
head(Z.Score.class_count_total)

ha_right = rowAnnotation(
  Abundance = Z.Score.class_count_total$Total)

row_labels_class = Z.Score.class_count_total$Class

hm_class <- Heatmap(Z.Score.class, name = "Z-score, CLR", col = col_matrix3,
                    column_title = "16S rRNA soil microbiome", 
                    column_title_gp = gpar(fontface = "bold", fontsize = 14),
                    top_annotation = ha_top,
                    right_annotation = ha_right,
                    row_title = "Class",
                    row_labels = row_labels_class,
                    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                    row_names_gp = gpar(fontsize = 10),
                    row_order = order(row_labels_class),
                    rect_gp = gpar(col = "white", lwd = 1),
                    show_column_names = TRUE,
                    show_heatmap_legend = TRUE)
hm_class


Z.Score.order_name <- as.data.frame(Z.Score.order)
str(Z.Score.order_name)
Z.Score.order_name <- rownames_to_column(Z.Score.order_name, var = "OTU")
Z.Score.order_name <- left_join(Z.Score.order_name, taxa_info)
Z.Score.order_name <- Z.Score.order_name[, -(2:51)]

Z.Score.order_count_total <- left_join(Z.Score.order_name, otu_count_total)
head(Z.Score.order_count_total)

ha_right = rowAnnotation(
  Abundance = Z.Score.order_count_total$Total)

row_labels_order = Z.Score.order_count_total$Order

hm_order <- Heatmap(Z.Score.order, name = "Z-score, CLR", col = col_matrix3,
                    column_title = "16S rRNA soil microbiome", 
                    column_title_gp = gpar(fontface = "bold", fontsize = 14),
                    top_annotation = ha_top,
                    right_annotation = ha_right,
                    row_title = "Order",
                    row_labels = row_labels_order,
                    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                    row_names_gp = gpar(fontsize = 10),
                    row_order = order(row_labels_order),
                    rect_gp = gpar(col = "white", lwd = 1),
                    show_column_names = TRUE,
                    show_heatmap_legend = TRUE)
hm_order

Z.Score.fam_name <- as.data.frame(Z.Score.fam)
str(Z.Score.fam_name)
Z.Score.fam_name <- rownames_to_column(Z.Score.fam_name, var = "OTU")
Z.Score.fam_name <- left_join(Z.Score.fam_name, taxa_info)
Z.Score.fam_name <- Z.Score.fam_name[, -(2:51)]

Z.Score.fam_count_total <- left_join(Z.Score.fam_name, otu_count_total)
head(Z.Score.fam_count_total)

ha_right = rowAnnotation(
  Abundance = Z.Score.fam_count_total$Total)

row_labels_fam = Z.Score.fam_count_total$Family

hm_fam <- Heatmap(Z.Score.fam, name = "Z-score, CLR", col = col_matrix3,
                  column_title = "16S rRNA soil microbiome", 
                  column_title_gp = gpar(fontface = "bold", fontsize = 14),
                  top_annotation = ha_top,
                  right_annotation = ha_right,
                  row_title = "Family",
                  row_labels = row_labels_fam,
                  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                  row_names_gp = gpar(fontsize = 10),
                  row_order = order(row_labels_fam),
                  rect_gp = gpar(col = "white", lwd = 1),
                  show_column_names = TRUE,
                  show_heatmap_legend = TRUE)
hm_fam


Z.Score.genus_name <- as.data.frame(Z.Score.genus)
str(Z.Score.genus_name)
Z.Score.genus_name <- rownames_to_column(Z.Score.genus_name, var = "OTU")
Z.Score.genus_name <- left_join(Z.Score.genus_name, taxa_info)
Z.Score.genus_name <- Z.Score.genus_name[, -(2:51)]

otu_count_total <- as.data.frame(otu_count)
otu_count_total$Total <- rowSums(otu_count_total[, -1])
head(otu_count_total)
otu_count_total <- otu_count_total[, -(2:51)]

Z.Score.genus_count_total <- left_join(Z.Score.genus_name, otu_count_total)
head(Z.Score.genus_count_total)

color_abundance = colorRamp2(c(0, 5000, 10000, 15000), c("#EFF3FF", "#BDD7E7", "#6BAED6", "#2171B5"))

ha_right_genus = rowAnnotation(
  Abundance = Z.Score.genus_count_total$Total)

row_labels_genus = Z.Score.genus_count_total$Genus

hm_genus <- Heatmap(Z.Score.genus, name = "Z-score, CLR", col = col_matrix3,
                    #clustering_distance_rows = "pearson",
                    column_title = "16S rRNA soil microbiome", 
                    column_title_gp = gpar(fontface = "bold", fontsize = 14),
                    column_km = 2,
                    border = TRUE,
                    top_annotation = ha_top,
                    right_annotation = ha_right_genus,
                    row_title = "Genus",
                    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
                    row_names_gp = gpar(fontsize = 6),
                    column_names_gp = gpar(fontsize = 10),
                    column_names_rot = 60,
                    row_dend_reorder = TRUE,
                    row_order = order(row_labels_genus),
                    #column_order = order(label)
                    column_dend_height = unit(2, "cm"), 
                    row_dend_width = unit(2, "cm"),
                    rect_gp = gpar(col = "white", lwd = 1),
                    show_column_names = TRUE,
                    row_labels = row_labels_genus,
                    show_row_names = TRUE,
                    show_heatmap_legend = TRUE)
hm_genus
#################### plot heatmap with ggplot2 ####################
# change Z.score objects from matrix to data.frame class
class(Z.Score.phylum)
Z.Score.phylum.df <- as.data.frame(t(Z.Score.phylum))
class(Z.Score.phylum.df)

Z.Score.class.df <- as.data.frame(t(Z.Score.class))

Z.Score.order.df <- as.data.frame(t(Z.Score.order))
Z.Score.order.df.top20 <- as.data.frame(t(Z.Score.order.top20))

Z.Score.fam.df <- as.data.frame(t(Z.Score.fam))
Z.Score.fam.df.top20 <- as.data.frame(t(Z.Score.fam.top20))

Z.Score.genus.df <- as.data.frame(t(Z.Score.genus))
Z.Score.genus.df.top20 <- as.data.frame(t(Z.Score.genus.top20))

# Use rownames to create another variable called sample
Z.Score.phylum.df <- tibble::rownames_to_column(Z.Score.phylum.df, "Sample")
Z.Score.class.df <- tibble::rownames_to_column(Z.Score.class.df, "Sample")
Z.Score.order.df <- tibble::rownames_to_column(Z.Score.order.df, "Sample")
Z.Score.order.df.top20 <- tibble::rownames_to_column(Z.Score.order.df.top20, "Sample")
Z.Score.fam.df <- tibble::rownames_to_column(Z.Score.fam.df, "Sample")
Z.Score.fam.df.top20 <- tibble::rownames_to_column(Z.Score.fam.df.top20, "Sample")
Z.Score.genus.df <- tibble::rownames_to_column(Z.Score.genus.df, "Sample")
Z.Score.genus.df.top20 <- tibble::rownames_to_column(Z.Score.genus.df.top20, "Sample")

#Add sample factors (Sample.Site, Habitat, Soil.Type) to the Z.score object
Z.Score.phylum.df.var <- data.frame(cbind(Z.Score.phylum.df, sample_soil_factors))
Z.Score.class.df.var <- data.frame(cbind(Z.Score.class.df, sample_soil_factors))
Z.Score.order.df.var <- data.frame(cbind(Z.Score.order.df, sample_soil_factors))
Z.Score.order.df.var.top20 <- data.frame(cbind(Z.Score.order.df.top20, sample_soil_factors))
Z.Score.fam.df.var <- data.frame(cbind(Z.Score.fam.df, sample_soil_factors))
Z.Score.fam.df.var.top20 <- data.frame(cbind(Z.Score.fam.df.top20, sample_soil_factors))
Z.Score.genus.df.var <- data.frame(cbind(Z.Score.genus.df, sample_soil_factors))
Z.Score.genus.df.var.top20 <- data.frame(cbind(Z.Score.genus.df.top20, sample_soil_factors))

#Add sample factors (Sample.Site, Habitat, Soil.Type) to the Z.score object
Z.Score.phylum.df.long <- pivot_longer(data = Z.Score.phylum.df.var,
                                       cols = -c(Sample, Sample.Type, Sample.Site, Habitat, Soil.Type),
                                       names_to = "ASV",
                                       values_to = "Z")
class(Z.Score.phylum.df.long)
Z.Score.phylum.df.long <- as.data.frame(Z.Score.phylum.df.long)
class(Z.Score.phylum.df.long)
str(Z.Score.phylum.df.long)
Z.Score.phylum.df.long_merged <- merge(Z.Score.phylum.df.long, sig_aldex2_phy_hab_result, by.x = 6, by.y = 1, all.x = TRUE)

Z.Score.class.df.long <- pivot_longer(data = Z.Score.class.df.var,
                                      cols = -c(Sample, Sample.Type, Sample.Site, Habitat, Soil.Type),
                                      names_to = "ASV",
                                      values_to = "Z")
Z.Score.class.df.long <- as.data.frame(Z.Score.class.df.long)
Z.Score.class.df.long_merged <- merge(Z.Score.class.df.long, sig_aldex2_class_hab_result, by.x = 6, by.y = 1, all.x = TRUE)

Z.Score.order.df.long <- pivot_longer(data = Z.Score.order.df.var,
                                      cols = -c(Sample, Sample.Type, Sample.Site, Habitat, Soil.Type),
                                      names_to = "ASV",
                                      values_to = "Z")
Z.Score.order.df.long <- as.data.frame(Z.Score.order.df.long)
Z.Score.order.df.long_merged <- merge(Z.Score.order.df.long, sig_aldex2_order_hab_result, by.x = 6, by.y = 1, all.x = TRUE)

Z.Score.order.df.long.top20 <- pivot_longer(data = Z.Score.order.df.var.top20,
                                            cols = -c(Sample, Sample.Type, Sample.Site, Habitat, Soil.Type),
                                            names_to = "ASV",
                                            values_to = "Z")
Z.Score.order.df.long.top20 <- as.data.frame(Z.Score.order.df.long.top20)
Z.Score.order.df.long_merged.top20 <- merge(Z.Score.order.df.long.top20, sig_aldex2_order_hab_result_top20, by.x = 6, by.y = 1, all.x = TRUE)

Z.Score.fam.df.long <- pivot_longer(data = Z.Score.fam.df.var,
                                    cols = -c(Sample, Sample.Type, Sample.Site, Habitat, Soil.Type),
                                    names_to = "ASV",
                                    values_to = "Z")
Z.Score.fam.df.long <- as.data.frame(Z.Score.fam.df.long)
Z.Score.fam.df.long_merged <- merge(Z.Score.fam.df.long, sig_aldex2_fam_hab_result, by.x = 6, by.y = 1, all.x = TRUE)

Z.Score.fam.df.long.top20 <- pivot_longer(data = Z.Score.fam.df.var.top20,
                                          cols = -c(Sample, Sample.Type, Sample.Site, Habitat, Soil.Type),
                                          names_to = "ASV",
                                          values_to = "Z")
Z.Score.fam.df.long.top20 <- as.data.frame(Z.Score.fam.df.long.top20)
Z.Score.fam.df.long_merged.top20 <- merge(Z.Score.fam.df.long.top20, sig_aldex2_fam_hab_result_top20, by.x = 6, by.y = 1, all.x = TRUE)

Z.Score.genus.df.long <- pivot_longer(data = Z.Score.genus.df.var,
                                      cols = -c(Sample, Sample.Type, Sample.Site, Habitat, Soil.Type),
                                      names_to = "ASV",
                                      values_to = "Z")
Z.Score.genus.df.long <- as.data.frame(Z.Score.genus.df.long)
Z.Score.genus.df.long_merged <- merge(Z.Score.genus.df.long, sig_aldex2_genus_hab_result, by.x = 6, by.y = 1, all.x = TRUE)

Z.Score.genus.df.long.top20 <- pivot_longer(data = Z.Score.genus.df.var.top20,
                                            cols = -c(Sample, Sample.Type, Sample.Site, Habitat, Soil.Type),
                                            names_to = "ASV",
                                            values_to = "Z")
Z.Score.genus.df.long.top20 <- as.data.frame(Z.Score.genus.df.long.top20)
Z.Score.genus.df.long_merged.top20 <- merge(Z.Score.genus.df.long.top20, sig_aldex2_genus_hab_result_top20, by.x = 6, by.y = 1, all.x = TRUE)

#Change facet habitat labels
soil.type.labs <- c("Compacted Soil","Uncompacted Soil")
names(soil.type.labs) <- c("Compact.Soil","Not.Compact.Soil")

hab.type.labs <- c("Coastal Scrub", "Native Grass", "Chaparral", "Riparian", "Hollyleaf Cherry", "Oak Wood")
names(hab.type.labs) <- c("Coastal.scrub", "Native.grass", "Chaparral", "Riparian", "Hollyleaf.cherry", "Oak.wood")

Z.Score.phylum.df.long_merged$Habitat <- factor(Z.Score.phylum.df.long_merged$Habitat,
                                                levels = c("Coastal.scrub", "Native.grass", "Chaparral", "Riparian", "Hollyleaf.cherry", "Oak.wood"))

heatmap_phy <- ggplot(Z.Score.phylum.df.long_merged, aes(x = Sample, y = Phylum, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ Soil.Type+Habitat, scales = "free", labeller = labeller(Soil.Type = soil.type.labs, Habitat = hab.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 10, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0, size = 10, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_phy
ggsave("16S_heatmap_phylum_soil_habitat_prune_noncontam_prev05_true_filt_greater_500.tiff", width = 10, height = 2, dpi = 150) # save graphic

Z.Score.class.df.long_merged$Habitat <- factor(Z.Score.class.df.long_merged$Habitat,
                                               levels = c("Coastal.scrub", "Native.grass", "Chaparral", "Riparian", "Hollyleaf.cherry", "Oak.wood"))
heatmap_class <- ggplot(Z.Score.class.df.long_merged, aes(x = Sample, y = Class, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ Soil.Type+Habitat, scales = "free", labeller = labeller(Soil.Type = soil.type.labs, Habitat = hab.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  #theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 10, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0, size = 10, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_class
ggsave("16S_heatmap_class_soil_habitat_prune_noncontam_prev05_true_filt_greater_500.tiff", width = 10, height = 4, dpi = 150) # save graphic

Z.Score.order.df.long_merged.top20$Habitat <- factor(Z.Score.order.df.long_merged.top20$Habitat,
                                                     levels = c("Coastal.scrub", "Native.grass", "Chaparral", "Riparian", "Hollyleaf.cherry", "Oak.wood"))
heatmap_order <- ggplot(Z.Score.order.df.long_merged.top20, aes(x = Sample, y = Order, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ Soil.Type+Habitat, scales = "free", labeller = labeller(Soil.Type = soil.type.labs, Habitat = hab.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 10, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0, size = 10, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_order
ggsave("16S_heatmap_order_soil_habitat_prune_noncontam_prev05_true_filt_greater_500.tiff", width = 10, height = 6, dpi = 150) # save graphic


Z.Score.fam.df.long_merged.top20$Habitat <- factor(Z.Score.fam.df.long_merged.top20$Habitat,
                                                   levels = c("Coastal.scrub", "Native.grass", "Chaparral", "Riparian", "Hollyleaf.cherry", "Oak.wood"))
heatmap_fam <- ggplot(Z.Score.fam.df.long_merged.top20, aes(x = Sample, y = Family, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "NULL") + # Add a nicer x-axis title
  facet_nested(~ Soil.Type+Habitat, scales = "free", labeller = labeller(Soil.Type = soil.type.labs, Habitat = hab.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 10, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0, size = 10, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_fam
ggsave("16S_heatmap_fam_soil_habitat_prune_noncontam_prev05_true_filt_greater_500.tiff", width = 10, height = 6, dpi = 150) # save graphic

Z.Score.genus.df.long_merged.top20$Habitat <- factor(Z.Score.genus.df.long_merged.top20$Habitat,
                                                     levels = c("Coastal.scrub", "Native.grass", "Chaparral", "Riparian", "Hollyleaf.cherry", "Oak.wood"))
heatmap_genus <- ggplot(Z.Score.genus.df.long_merged.top20, aes(x = Sample, y = Genus, fill = Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  labs(x= "Sample") + # Add a nicer x-axis title
  facet_nested(~ Soil.Type+Habitat, scales = "free", labeller = labeller(Soil.Type = soil.type.labs, Habitat = hab.type.labs)) + # # facet_grid makes two panels, compacted vs. not compacted soil, and separates per habitat
  theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 10, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 9)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 9)) + # adjusts the title of x axis
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 8, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1.1, size = 8, color = "black")) + # adjusts text of x axis
  scale_y_discrete(limits=rev, expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_genus
ggsave("16S_heatmap_genus_soil_habitat_prune_noncontam_prev05_true_filt_greater_500.tiff", width = 12, height = 4, dpi = 150) # save graphic

Figure_heatmap_taxa <- ggarrange(heatmap_phy,  heatmap_class, heatmap_order, heatmap_fam, heatmap_genus, heights = unit(c(6, 17, 20, 20, 16), "lines"), ncol = 1, nrow = 5, align = "v", common.legend = TRUE, legend = "right") 
Figure_heatmap_taxa
ggsave("16S_heatmap_all_levels_soil_habitat_prune_noncontam_prev05_true_filt_greater_500.tiff", width = 12, height = 20, dpi = 150) # save graphic

# write ASV_Soil matrix to file as csv
ASV_SHSK_soil = as(otu_table(phy_soil_prune_noncontam_comb05_true_scale_filt), "matrix")
ASV_SHSK_soil_rank1 = as(otu_table(phy_soil_prune_noncontam_comb05_true_scale_filt_rank1), "matrix")
ASV_SHSK_soil_rank2 = as(otu_table(phy_soil_prune_noncontam_comb05_true_scale_filt_rank2), "matrix")
ASV_SHSK_soil_rank3 = as(otu_table(phy_soil_prune_noncontam_comb05_true_scale_filt_rank3), "matrix")
ASV_SHSK_soil_rank4 = as(otu_table(phy_soil_prune_noncontam_comb05_true_scale_filt_rank4), "matrix")
ASV_SHSK_soil_rank5 = as(otu_table(phy_soil_prune_noncontam_comb05_true_scale_filt_rank5), "matrix")
ASV_SHSK_soil_rank6 = as(otu_table(phy_soil_prune_noncontam_comb05_true_scale_filt_rank6), "matrix")

# Transpose matrix
if(taxa_are_rows(phy_soil_prune_noncontam_comb05_true_scale_filt)) {ASV_SHSK_soil <- t(ASV_SHSK_soil)}
if(taxa_are_rows(phy_soil_prune_noncontam_comb05_true_scale_filt_rank1)) {ASV_SHSK_soil_rank1 <- t(ASV_SHSK_soil_rank1)}
if(taxa_are_rows(phy_soil_prune_noncontam_comb05_true_scale_filt_rank2)) {ASV_SHSK_soil_rank2 <- t(ASV_SHSK_soil_rank2)}
if(taxa_are_rows(phy_soil_prune_noncontam_comb05_true_scale_filt_rank3)) {ASV_SHSK_soil_rank3 <- t(ASV_SHSK_soil_rank3)}
if(taxa_are_rows(phy_soil_prune_noncontam_comb05_true_scale_filt_rank4)) {ASV_SHSK_soil_rank4 <- t(ASV_SHSK_soil_rank4)}
if(taxa_are_rows(phy_soil_prune_noncontam_comb05_true_scale_filt_rank5)) {ASV_SHSK_soil_rank5 <- t(ASV_SHSK_soil_rank5)}
if(taxa_are_rows(phy_soil_prune_noncontam_comb05_true_scale_filt_rank6)) {ASV_SHSK_soil_rank6 <- t(ASV_SHSK_soil_rank6)}

# Coerce to a data frame
ASV_SHSK_soil_df = as.data.frame(ASV_SHSK_soil)
class(ASV_SHSK_soil_df)

ASV_SHSK_soil_df_rank1 = as.data.frame(ASV_SHSK_soil_rank1)
class(ASV_SHSK_soil_df_rank1)

ASV_SHSK_soil_df_rank2 = as.data.frame(ASV_SHSK_soil_rank2)
ASV_SHSK_soil_df_rank3 = as.data.frame(ASV_SHSK_soil_rank3)
ASV_SHSK_soil_df_rank4 = as.data.frame(ASV_SHSK_soil_rank4)
ASV_SHSK_soil_df_rank5 = as.data.frame(ASV_SHSK_soil_rank5)
ASV_SHSK_soil_df_rank6 = as.data.frame(ASV_SHSK_soil_rank6)

# write ASV_df matrix to file as csv
write.csv(ASV_SHSK_soil_df, "/Users/tiagopereira/Dropbox/SK-data/16S-soil-single-nematode/R-codes/ASV_SHSK_soil_df.csv")
write.csv(ASV_SHSK_soil_df_rank1, "/Users/tiagopereira/Dropbox/SK-data/16S-soil-single-nematode/R-codes/ASV_SHSK_soil_df_rank1.csv")
write.csv(ASV_SHSK_soil_df_rank2, "/Users/tiagopereira/Dropbox/SK-data/16S-soil-single-nematode/R-codes/ASV_SHSK_soil_df_rank2.csv")
write.csv(ASV_SHSK_soil_df_rank3, "/Users/tiagopereira/Dropbox/SK-data/16S-soil-single-nematode/R-codes/ASV_SHSK_soil_df_rank3.csv")
write.csv(ASV_SHSK_soil_df_rank4, "/Users/tiagopereira/Dropbox/SK-data/16S-soil-single-nematode/R-codes/ASV_SHSK_soil_df_rank4.csv")
write.csv(ASV_SHSK_soil_df_rank5, "/Users/tiagopereira/Dropbox/SK-data/16S-soil-single-nematode/R-codes/ASV_SHSK_soil_df_rank5.csv")
write.csv(ASV_SHSK_soil_df_rank6, "/Users/tiagopereira/Dropbox/SK-data/16S-soil-single-nematode/R-codes/ASV_SHSK_soil_df_rank6.csv")

# write TAX_Soil matrix to file as csv
TAX_SHSK_soil = as(tax_table(phy_soil_prune_noncontam_comb05_true_scale_filt), "matrix")
TAX_SHSK_soil_df = as.data.frame(TAX_SHSK_soil)
class(TAX_SHSK_soil_df)

TAX_SHSK_soil_rank1 = as(tax_table(phy_soil_prune_noncontam_comb05_true_scale_filt_rank1), "matrix")
TAX_SHSK_soil_rank1_df = as.data.frame(TAX_SHSK_soil_rank1)
class(TAX_SHSK_soil_rank1_df)

TAX_SHSK_soil_rank2 = as(tax_table(phy_soil_prune_noncontam_comb05_true_scale_filt_rank2), "matrix")
TAX_SHSK_soil_rank2_df = as.data.frame(TAX_SHSK_soil_rank2)
TAX_SHSK_soil_rank3 = as(tax_table(phy_soil_prune_noncontam_comb05_true_scale_filt_rank3), "matrix")
TAX_SHSK_soil_rank3_df = as.data.frame(TAX_SHSK_soil_rank3)
TAX_SHSK_soil_rank4 = as(tax_table(phy_soil_prune_noncontam_comb05_true_scale_filt_rank4), "matrix")
TAX_SHSK_soil_rank4_df = as.data.frame(TAX_SHSK_soil_rank4)
TAX_SHSK_soil_rank5 = as(tax_table(phy_soil_prune_noncontam_comb05_true_scale_filt_rank5), "matrix")
TAX_SHSK_soil_rank5_df = as.data.frame(TAX_SHSK_soil_rank5)
TAX_SHSK_soil_rank6 = as(tax_table(phy_soil_prune_noncontam_comb05_true_scale_filt_rank6), "matrix")
TAX_SHSK_soil_rank6_df = as.data.frame(TAX_SHSK_soil_rank6)

# write TAX_df matrix to file as csv
write.csv(TAX_SHSK_soil_df, "/Users/tiagopereira/Dropbox/SK-data/16S-soil-single-nematode/R-codes/TAX_SHSK_soil_df.csv")
write.csv(TAX_SHSK_soil_rank1_df, "/Users/tiagopereira/Dropbox/SK-data/16S-soil-single-nematode/R-codes/TAX_SHSK_soil_rank1_df.csv")
write.csv(TAX_SHSK_soil_rank2_df, "/Users/tiagopereira/Dropbox/SK-data/16S-soil-single-nematode/R-codes/TAX_SHSK_soil_rank2_df.csv")
write.csv(TAX_SHSK_soil_rank3_df, "/Users/tiagopereira/Dropbox/SK-data/16S-soil-single-nematode/R-codes/TAX_SHSK_soil_rank3_df.csv")
write.csv(TAX_SHSK_soil_rank4_df, "/Users/tiagopereira/Dropbox/SK-data/16S-soil-single-nematode/R-codes/TAX_SHSK_soil_rank4_df.csv")
write.csv(TAX_SHSK_soil_rank5_df, "/Users/tiagopereira/Dropbox/SK-data/16S-soil-single-nematode/R-codes/TAX_SHSK_soil_rank5_df.csv")
write.csv(TAX_SHSK_soil_rank6_df, "/Users/tiagopereira/Dropbox/SK-data/16S-soil-single-nematode/R-codes/TAX_SHSK_soil_rank6_df.csv")

#####################################
phy_soil_prune_noncontam_comb05_true_scale_filt_phylum_20 = tax_glom(phy_soil_prune_noncontam_comb05_true_scale_filt, taxrank = "Phylum") # agglomerate at phylum level
phy_soil_prune_noncontam_comb05_true_scale_filt_phylum_20_sorted <- names(sort(taxa_sums(phy_soil_prune_noncontam_comb05_true_scale_filt_phylum_20), decreasing=TRUE)[1:20]) # sort the top 20 Phyla
phy_soil_prune_noncontam_comb05_true_scale_filt_phylum_20_trans = transform_sample_counts(phy_soil_prune_noncontam_comb05_true_scale_filt_phylum_20, function(x) x/sum(x)) # Transform to rel. abundance
phy_soil_prune_noncontam_comb05_true_scale_filt_phylum_20_prune = prune_taxa(phy_soil_prune_noncontam_comb05_true_scale_filt_phylum_20_sorted, phy_soil_prune_noncontam_comb05_true_scale_filt_phylum_20_trans) # subtract the top20 taxa from the entire phyloseq object
phy_soil_prune_noncontam_comb05_true_scale_filt_phylum_20_prune_df = psmelt(phy_soil_prune_noncontam_comb05_true_scale_filt_phylum_20_prune) # Melt to long format
phy_soil_prune_noncontam_comb05_true_scale_filt_phylum_20_prune_agr = aggregate(Abundance~Soil.Type+Habitat+Phylum, data=phy_soil_prune_noncontam_comb05_true_scale_filt_phylum_20_prune_df, FUN=mean) # add common factors to use for plotting

phy_soil_prune_noncontam_comb05_true_scale_filt_class_20 = tax_glom(phy_soil_prune_noncontam_comb05_true_scale_filt, taxrank = "Class") # agglomerate at Class level
phy_soil_prune_noncontam_comb05_true_scale_filt_class_20_sorted <- names(sort(taxa_sums(phy_soil_prune_noncontam_comb05_true_scale_filt_class_20), decreasing=TRUE)[1:20]) # sort the top 20 Class
phy_soil_prune_noncontam_comb05_true_scale_filt_class_20_trans = transform_sample_counts(phy_soil_prune_noncontam_comb05_true_scale_filt_class_20, function(x) x/sum(x)) # Transform to rel. abundance
phy_soil_prune_noncontam_comb05_true_scale_filt_class_20_prune = prune_taxa(phy_soil_prune_noncontam_comb05_true_scale_filt_class_20_sorted, phy_soil_prune_noncontam_comb05_true_scale_filt_class_20_trans) # subtract the top20 taxa from the entire phyloseq object
phy_soil_prune_noncontam_comb05_true_scale_filt_class_20_prune_df = psmelt(phy_soil_prune_noncontam_comb05_true_scale_filt_class_20_prune) # Melt to long format
phy_soil_prune_noncontam_comb05_true_scale_filt_class_20_prune_agr = aggregate(Abundance~Soil.Type+Habitat+Class, data=phy_soil_prune_noncontam_comb05_true_scale_filt_class_20_prune_df, FUN=mean) # add common factors to use for plotting

phy_soil_prune_noncontam_comb05_true_scale_filt_order_20 = tax_glom(phy_soil_prune_noncontam_comb05_true_scale_filt, taxrank = "Order") # agglomerate at Order level
phy_soil_prune_noncontam_comb05_true_scale_filt_order_20_sorted <- names(sort(taxa_sums(phy_soil_prune_noncontam_comb05_true_scale_filt_order_20), decreasing=TRUE)[1:20]) # sort the top 20 Phyla
phy_soil_prune_noncontam_comb05_true_scale_filt_order_20_trans = transform_sample_counts(phy_soil_prune_noncontam_comb05_true_scale_filt_order_20, function(x) x/sum(x)) # Transform to rel. abundance
phy_soil_prune_noncontam_comb05_true_scale_filt_order_20_prune = prune_taxa(phy_soil_prune_noncontam_comb05_true_scale_filt_order_20_sorted, phy_soil_prune_noncontam_comb05_true_scale_filt_order_20_trans) # subtract the top20 taxa from the entire phyloseq object
phy_soil_prune_noncontam_comb05_true_scale_filt_order_20_prune_df = psmelt(phy_soil_prune_noncontam_comb05_true_scale_filt_order_20_prune) # Melt to long format
phy_soil_prune_noncontam_comb05_true_scale_filt_order_20_prune_agr = aggregate(Abundance~Soil.Type+Habitat+Order, data=phy_soil_prune_noncontam_comb05_true_scale_filt_order_20_prune_df, FUN=mean) # add common factors to use for plotting

# Melt to long format (for ggplot2); prune out phyla below 1% in each sample
phy_soil_prune_noncontam_comb05_true_scale_filt_phylum <- phy_soil_prune_noncontam_comb05_true_scale_filt %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum
phy_soil_prune_noncontam_comb05_true_scale_filt_phylum
str(phy_soil_prune_noncontam_comb05_true_scale_filt_phylum)

# Melt to long format (for ggplot2); prune out class below 1% in each sample
phy_soil_prune_noncontam_comb05_true_scale_filt_class <- phy_soil_prune_noncontam_comb05_true_scale_filt %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at Class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Class)                                      # Sort data frame alphabetically by phylum
phy_soil_prune_noncontam_comb05_true_scale_filt_class

# Melt to long format (for ggplot2); prune out order below 1% in each sample
phy_soil_prune_noncontam_comb05_true_scale_filt_order <- phy_soil_prune_noncontam_comb05_true_scale_filt %>%
  tax_glom(taxrank = "Order") %>%                     # agglomerate at Class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Order)                                      # Sort data frame alphabetically by phylum
phy_soil_prune_noncontam_comb05_true_scale_filt_order

# Melt to long format (for ggplot2); prune out family below 1% in each sample
phy_soil_prune_noncontam_comb05_true_scale_filt_fam <- phy_soil_prune_noncontam_comb05_true_scale_filt %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at Class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Family)                                      # Sort data frame alphabetically by phylum
phy_soil_prune_noncontam_comb05_true_scale_filt_fam

# Melt to long format (for ggplot2); prune out genus below 1% in each sample
phy_soil_prune_noncontam_comb05_true_scale_filt_genus <- phy_soil_prune_noncontam_comb05_true_scale_filt %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Genus)                                      # Sort data frame alphabetically by phylum
phy_soil_prune_noncontam_comb05_true_scale_filt_genus

# Melt to long format (for ggplot2); prune out species below 1% in each sample
phy_soil_prune_noncontam_comb05_true_scale_filt_spp <- phy_soil_prune_noncontam_comb05_true_scale_filt %>%
  tax_glom(taxrank = "Species") %>%                     # agglomerate at Class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Species)                                      # Sort data frame alphabetically by phylum
phy_soil_prune_noncontam_comb05_true_scale_filt_spp


#Get unique taxa for a specific taxonomic rank
unique(phy_soil_prune_noncontam_comb05_true_scale_filt_king$Kingdom)
unique(phy_soil_prune_noncontam_comb05_true_scale_filt_phylum$Phylum)
unique(phy_soil_prune_noncontam_comb05_true_scale_filt_class$Class)
unique(phy_soil_prune_noncontam_comb05_true_scale_filt_order$Order)
unique(phy_soil_prune_noncontam_comb05_true_scale_filt_fam$Family)
unique(phy_soil_prune_noncontam_comb05_true_scale_filt_genus$Genus)
unique(phy_soil_prune_noncontam_comb05_true_scale_filt_spp$Species)

# Define the number of colors you want for each taxonomic level
nb.cols.phylum <- 20
soil_colors_phylum <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.phylum)

nb.cols.class <- 20
soil_colors_class <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.class)

nb.cols.order <- 20
soil_colors_order <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.order)

nb.cols.family <- 20
soil_colors_family <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.family)

nb.cols.genus <- 20
soil_colors_genus <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.family)

#Change facet habitat labels
soil.type.labs <- c("Compacted Soil","Uncompacted Soil")
names(soil.type.labs) <- c("Compact.Soil","Not.Compact.Soil")

# Plot bar chart 16S Soil samples per habitat - Phylum level
ggplot(phy_soil_prune_noncontam_comb05_true_scale_filt_phylum_20_prune_agr, aes(x = Habitat, y = Abundance, fill = Phylum)) + #plotting by sample
  facet_wrap(. ~ Soil.Type, scales = "free_x", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = soil_colors_phylum) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.0, 0.0)) + # removes the internal margins
  scale_x_discrete(
    breaks = c("Chaparral", "Coastal.scrub", "Hollyleaf.cherry", "Native.grass", "Oak.wood", "Riparian"),
    labels = c("Chaparral", "Coastal Scrub", "Hollyleaf Cherry", "Native Grass", "Oak Wood", "Riparian"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  #scale_x_discrete(
  #  breaks = c("SingleWorm, Soil"),
  # labels = c("Single Worm, Soil"),
  # expand = c(0.0, 0.0),
  # drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample Type") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
ggsave("16S_Phylum_Soil_Habitat_ShSk_barplot_greater_500_filt_comb.tiff", width = 6, height = 6, dpi = 150) # save graphic

# Plot bar chart 16S Soil samples per habitat - Class level
ggplot(phy_soil_prune_noncontam_comb05_true_scale_filt_class_20_prune_agr, aes(x = Habitat, y = Abundance, fill = Class)) + #plotting by sample
  facet_wrap(. ~ Soil.Type, scales = "free_x", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = soil_colors_class) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.0, 0.0)) + # removes the internal margins
  scale_x_discrete(
    breaks = c("Chaparral", "Coastal.scrub", "Hollyleaf.cherry", "Native.grass", "Oak.wood", "Riparian"),
    labels = c("Chaparral", "Coastal Scrub", "Hollyleaf Cherry", "Native Grass", "Oak Wood", "Riparian"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  #scale_x_discrete(
  #  breaks = c("SingleWorm, Soil"),
  # labels = c("Single Worm, Soil"),
  # expand = c(0.0, 0.0),
  # drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample Type") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
ggsave("16S_class_Soil_Habitat_ShSk_barplot_greater_500_filt_comb.tiff", width = 6, height = 6, dpi = 150) # save graphic

# Plot bar chart 16S Soil samples per habitat - Order level
ggplot(phy_soil_prune_noncontam_comb05_true_scale_filt_order_20_prune_agr, aes(x = Habitat, y = Abundance, fill = Order)) + #plotting by sample
  facet_wrap(. ~ Soil.Type, scales = "free_x", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = soil_colors_order) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.0, 0.0)) + # removes the internal margins
  scale_x_discrete(
    breaks = c("Chaparral", "Coastal.scrub", "Hollyleaf.cherry", "Native.grass", "Oak.wood", "Riparian"),
    labels = c("Chaparral", "Coastal Scrub", "Hollyleaf Cherry", "Native Grass", "Oak Wood", "Riparian"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  #scale_x_discrete(
  #  breaks = c("SingleWorm, Soil"),
  # labels = c("Single Worm, Soil"),
  # expand = c(0.0, 0.0),
  # drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Sample Type") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
ggsave("16S_order_Soil_Habitat_ShSk_barplot_greater_500_filt_comb.tiff", width = 6, height = 6, dpi = 150) # save graphic

# Plot bar chart 16S Soil samples per habitat - Family level
ggplot(phy_soil_prune_noncontam_comb05_true_scale_filt_fam, aes(x = Habitat, y = Abundance, fill = Family)) + #plotting by sample
  facet_wrap(. ~ Soil.Type, scales = "free_x", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = soil_colors_family) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.0, 0.0)) + # removes the internal margins
  scale_x_discrete(
    breaks = c("Chaparral", "Coastal.scrub", "Hollyleaf.cherry", "Native.grass", "Oak.wood", "Riparian"),
    labels = c("Chaparral", "Coastal Scrub", "Hollyleaf Cherry", "Native Grass", "Oak Wood", "Riparian"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  #scale_x_discrete(
  #  breaks = c("SingleWorm, Soil"),
  # labels = c("Single Worm, Soil"),
  # expand = c(0.0, 0.0),
  # drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (Family > 1%)") + # add the title on y axis
  xlab("Sample Type") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
ggsave("16S_Family_Soil_Habitat_ShSk_barplot_greater_500_filt_comb.tiff", width = 12, height = 6, dpi = 150) # save graphic  

# Plot bar chart 16S Soil samples per habitat - Genus level
ggplot(phy_soil_prune_noncontam_comb05_true_scale_filt_genus, aes(x = Habitat, y = Abundance, fill = Genus)) + #plotting by sample
  facet_wrap(. ~ Soil.Type, scales = "free_x", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = soil_colors_genus) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.0, 0.0)) + # removes the internal margins
  scale_x_discrete(
    breaks = c("Chaparral", "Coastal.scrub", "Hollyleaf.cherry", "Native.grass", "Oak.wood", "Riparian"),
    labels = c("Chaparral", "Coastal Scrub", "Hollyleaf Cherry", "Native Grass", "Oak Wood", "Riparian"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  #scale_x_discrete(
  #  breaks = c("SingleWorm, Soil"),
  # labels = c("Single Worm, Soil"),
  # expand = c(0.0, 0.0),
  # drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (Genus > 1%)") + # add the title on y axis
  xlab("Sample Type") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
ggsave("16S_Genus_Soil_Habitat_ShSk_barplot_greater_500_filt_comb.tiff", width = 12, height = 6, dpi = 150) # save graphic  


########## Plot only the top 20 taxa for the different taxonomic ranks

#### Step 1
# Use psmelt of the phyloseq object, so it can be taking by gglplot  
phy_soil_prune_noncontam_comb05_true_scale_filt_df <- psmelt(phy_soil_prune_noncontam_comb05_true_scale_filt)

# Transform to a dataframe and save as csv
phy_soil_prune_noncontam_comb05_true_scale_filt_df <- as.data.frame(phy_soil_prune_noncontam_comb05_true_scale_filt_df)
write.csv(phy_soil_prune_noncontam_comb05_true_scale_filt_df, "phy_soil_prune_noncontam_comb05_true_scale_filt_df.csv")

#### Step 2 Class Level
df_phy_soil_prune_noncontam_comb05_true_scale_filt <- read.csv("phy_soil_prune_noncontam_comb05_true_scale_filt_df.csv")

split_frame_df_phy_soil_prune_noncontam_comb05_true_scale_filt <- split(df_phy_soil_prune_noncontam_comb05_true_scale_filt, list(df_phy_soil_prune_noncontam_comb05_true_scale_filt$Sample))

rel_ab_df_phy_soil_prune_noncontam_comb05_true_scale_filt <- vector("list", length(split_frame_df_phy_soil_prune_noncontam_comb05_true_scale_filt))
for (i in 1:length(split_frame_df_phy_soil_prune_noncontam_comb05_true_scale_filt)){
  rel_ab_df_phy_soil_prune_noncontam_comb05_true_scale_filt[[i]] <- group_by(split_frame_df_phy_soil_prune_noncontam_comb05_true_scale_filt[[i]]) %>% mutate(Abundance = 100 * Abundance / sum(Abundance))
}

append_df_phy_soil_prune_noncontam_comb05_true_scale_filt <- ldply(rel_ab_df_phy_soil_prune_noncontam_comb05_true_scale_filt, data.frame)

grouped_taxa_df_phy_soil_prune_noncontam_comb05_true_scale_filt <- ddply(append_df_phy_soil_prune_noncontam_comb05_true_scale_filt, .(Class), summarize, Abundance = sum(Abundance))

topx_df_phy_soil_prune_noncontam_comb05_true_scale_filt <-  grouped_taxa_df_phy_soil_prune_noncontam_comb05_true_scale_filt %>% top_n(n = 20, wt = Abundance)

other_group_df_phy_soil_prune_noncontam_comb05_true_scale_filt <- subset(grouped_taxa_df_phy_soil_prune_noncontam_comb05_true_scale_filt, !(Class %in% topx_df_phy_soil_prune_noncontam_comb05_true_scale_filt$Class))

other_group_filter_df_phy_soil_prune_noncontam_comb05_true_scale_filt <- append_df_phy_soil_prune_noncontam_comb05_true_scale_filt %>% filter(Class %in% other_group_df_phy_soil_prune_noncontam_comb05_true_scale_filt$Class)

other_group_filter_df_phy_soil_prune_noncontam_comb05_true_scale_filt$Class <- "Other"

append.topx_df_phy_soil_prune_noncontam_comb05_true_scale_filt <- append_df_phy_soil_prune_noncontam_comb05_true_scale_filt %>% filter(Class %in% topx_df_phy_soil_prune_noncontam_comb05_true_scale_filt$Class)

merged_df_phy_soil_prune_noncontam_comb05_true_scale_filt<- rbind(append.topx_df_phy_soil_prune_noncontam_comb05_true_scale_filt, other_group_filter_df_phy_soil_prune_noncontam_comb05_true_scale_filt)

final_df_phy_soil_prune_noncontam_comb05_true_scale_filt <- as.data.frame(merged_df_phy_soil_prune_noncontam_comb05_true_scale_filt)

final_df_phy_soil_prune_noncontam_comb05_true_scale_filt[is.na(final_df_phy_soil_prune_noncontam_comb05_true_scale_filt)] <- "Other"

# Step 3 plotting it
unique(final_df_phy_soil_prune_noncontam_comb05_true_scale_filt$Class)

soil_colors_unique <- c("#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277", "#52AF43", "#6F9E4C", "#DD9A88", "#F16667", "#E42022",
                                 "#F06C45", "#FDBB69", "#FE982C", "#F78620", "#D9A295", "#B294C7", "#7D54A5", "#9E8099", "#F0EB99", "#DBB466", "gray")
                                 
final_df_phy_soil_prune_noncontam_comb05_true_scale_filt$Class <- factor(final_df_phy_soil_prune_noncontam_comb05_true_scale_filt$Class,
                                                                         levels = c("Acidobacteriae", "Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidia", "Blastocatellia", "Clostridia", "Gammaproteobacteria",
                                                                                    "Gemmatimonadetes", "Holophagae", "KD4-96", "Myxococcia", "Nitrososphaeria", "Oligoflexia", "Planctomycetes",
                                                                                    "Polyangia", "Rubrobacteria", "Thermoleophilia", "Verrucomicrobiae", "Vicinamibacteria", "Other"))

# Plot bar chart 16S Soil samples per habitat - Class level
ggplot(final_df_phy_soil_prune_noncontam_comb05_true_scale_filt, aes(x = Habitat, y = Abundance, fill = Class)) + #plotting by sample
  facet_wrap(. ~ Soil.Type, scales = "free_x", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = soil_colors_unique) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.0, 0.0)) + # removes the internal margins
  scale_x_discrete(
    breaks = c("Chaparral", "Coastal.scrub", "Hollyleaf.cherry", "Native.grass", "Oak.wood", "Riparian"),
    labels = c("Chaparral", "Coastal Scrub", "Hollyleaf Cherry", "Native Grass", "Oak Wood", "Riparian"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  #scale_x_discrete(
  #  breaks = c("SingleWorm, Soil"),
  # labels = c("Single Worm, Soil"),
  # expand = c(0.0, 0.0),
  # drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (Class > 1%)") + # add the title on y axis
  xlab("Sample Type") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
ggsave("16S_Class_Top20_Soil_Habitat_ShSk_barplot_greater_500_filt_comb.tiff", width = 7, height = 4, dpi = 150) # save graphic  

#### Step 1
# Use psmelt of the phyloseq object, so it can be taking by gglplot  
phy_soil_prune_noncontam_comb05_true_scale_filt_df_fam <- psmelt(phy_soil_prune_noncontam_comb05_true_scale_filt)

# Transform to a dataframe and save as csv
phy_soil_prune_noncontam_comb05_true_scale_filt_df_fam <- as.data.frame(phy_soil_prune_noncontam_comb05_true_scale_filt_df_fam)
write.csv(phy_soil_prune_noncontam_comb05_true_scale_filt_df_fam, "phy_soil_prune_noncontam_comb05_true_scale_filt_df_fam.csv")

#### Step 2 Family Level
df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam <- read.csv("phy_soil_prune_noncontam_comb05_true_scale_filt_df_fam.csv")

split_frame_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam <- split(df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam, list(df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam$Sample))

rel_ab_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam <- vector("list", length(split_frame_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam))
for (i in 1:length(split_frame_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam)){
  rel_ab_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam[[i]] <- group_by(split_frame_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam[[i]]) %>% mutate(Abundance = 100 * Abundance / sum(Abundance))
}

append_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam <- ldply(rel_ab_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam, data.frame)

grouped_taxa_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam <- ddply(append_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam, .(Family), summarize, Abundance = sum(Abundance))

topx_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam <-  grouped_taxa_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam %>% top_n(n = 20, wt = Abundance)

other_group_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam <- subset(grouped_taxa_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam, !(Family %in% topx_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam$Family))

other_group_filter_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam <- append_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam %>% filter(Family %in% other_group_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam$Family)

other_group_filter_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam$Family <- "Other"

append.topx_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam <- append_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam %>% filter(Family %in% topx_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam$Family)

merged_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam<- rbind(append.topx_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam, other_group_filter_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam)

final_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam <- as.data.frame(merged_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam)

final_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam[is.na(final_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam)] <- "Other"

# Step 3 plotting it

unique(final_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam$Family)

soil_colors_unique <- c("#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277", "#52AF43", "#6F9E4C", "#DD9A88", "#F16667", "#E42022",
                                 "#F06C45", "#FDBB69", "#FE982C", "#F78620", "#D9A295", "#B294C7", "#7D54A5", "#9E8099", "#F0EB99", "#DBB466", "gray")
                                 
final_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam$Family <- factor(final_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam$Family,
                                                                              levels = c("67-14","Bacillaceae","Beijerinckiaceae","Chitinophagaceae","Chthoniobacteraceae","Comamonadaceae","Gemmatimonadaceae","Geodermatophilaceae",
                                                                                         "Micromonosporaceae","Nitrososphaeraceae","Oxalobacteraceae","Pyrinomonadaceae","Rubrobacteriaceae","Solirubrobacteraceae","Sphingomonadaceae",
                                                                                         "Unassigned","uncultured bacterium","uncultured","Vicinamibacteraceae","Xanthobacteraceae","Other"))
# Plot bar chart 16S Soil samples per habitat - Family level
ggplot(final_df_phy_soil_prune_noncontam_comb05_true_scale_filt_fam, aes(x = Habitat, y = Abundance, fill = Family)) + #plotting by sample
  facet_wrap(. ~ Soil.Type, scales = "free_x", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = soil_colors_unique) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.0, 0.0)) + # removes the internal margins
  scale_x_discrete(
    breaks = c("Chaparral", "Coastal.scrub", "Hollyleaf.cherry", "Native.grass", "Oak.wood", "Riparian"),
    labels = c("Chaparral", "Coastal Scrub", "Hollyleaf Cherry", "Native Grass", "Oak Wood", "Riparian"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  #scale_x_discrete(
  #  breaks = c("SingleWorm, Soil"),
  # labels = c("Single Worm, Soil"),
  # expand = c(0.0, 0.0),
  # drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (Family > 1%)") + # add the title on y axis
  xlab("Sample Type") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
ggsave("16S_Class_Top20_Soil_Habitat_ShSk_barplot_greater_500_filt_comb.tiff", width = 7, height = 4, dpi = 150) # save graphic  

####################################### Soil and Nematode samples together #######################################

# Remove low abundance samples
smin <- min(sample_sums(phy))
phy_prune <- prune_samples(sample_sums(phy)>0, phy)
smin <- min(sample_sums(phy_prune))
phy_prune

# Check for ASVs that have no counted reads
any(taxa_sums(phy_prune) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_prune) == 0) # gives the number of cases
phy_prune <-prune_taxa(taxa_sums(phy_prune) > 0, phy_prune)

#Inspect Library sizes
df_phy_prune <- as.data.frame(sample_data(phy_prune))
df_phy_prune$LibrarySize <- sample_sums(phy_prune)
df_phy_prune <- df_phy_prune[order(df_phy_prune$LibrarySize),]
df_phy_prune$Index <- seq(nrow(df_phy_prune))
ggplot(data = df_phy_prune, aes(x= Index, y= LibrarySize, color = Sample.Type)) +
  geom_point()

# Identify Contaminants - Frequency
contamdf_freq <- isContaminant(phy_prune, method="frequency", conc="QubitConc")
head(contamdf_freq)
tail(contamdf_freq)
table(contamdf_freq$contaminant)
head(which(contamdf_freq$contaminant))

# Explore AVSs that are potentially contaminants
plot_frequency(phy_prune, taxa_names(phy_prune)[c(1, 58)], conc="QubitConc") + 
  xlab("DNA Concentration (ng/ul)")

# Explore AVSs that are potentially contaminants, randomly choosing 4 AVSs.
set.seed(100)
plot_frequency(phy_prune, taxa_names(phy_prune)[sample(which(contamdf_freq$contaminant),4)], conc="QubitConc") + 
  xlab("DNA Concentration (ng/ul)")

# Identify Contaminants - Prevalence
sample_data(phy_prune)$is.neg <- sample_data(phy_prune)$Sample.Control == "Control.Sample"
contamdf_prev <- isContaminant(phy_prune, method="prevalence", neg="is.neg")
table(contamdf_prev$contaminant)
head(which(contamdf_prev$contaminant))

contamdf_prev05 <- isContaminant(phy_prune, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf_prev05$contaminant)

# Make phyloseq object of presence-absence in negative controls and true samples
phy_prune_pa <- transform_sample_counts(phy_prune, function(abund) 1*(abund>0))
phy_prune_pa_neg <- prune_samples(sample_data(phy_prune_pa)$Sample.Control == "Control.Sample", phy_prune_pa)
phy_prune_pa_pos <- prune_samples(sample_data(phy_prune_pa)$Sample.Control == "True.Sample", phy_prune_pa)

# Make data.frame of prevalence in positive and negative samples
df_phy_prune_pa <- data.frame(pa.pos=taxa_sums(phy_prune_pa_pos), pa.neg=taxa_sums(phy_prune_pa_neg),
                              contaminant=contamdf_prev$contaminant)
ggplot(data=df_phy_prune_pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Remove contaminants from phyloseq object - prevalence based
phy_prune_noncontam_prev <- prune_taxa(!contamdf_prev$contaminant, phy_prune)
phy_prune_noncontam_prev
smin <-min(sample_sums(phy_prune_noncontam_prev))
any(sample_sums(phy_prune_noncontam_prev) == 0) # if TRUE, there are samples with no reads, otherwise FALSE
sum(sample_sums(phy_prune_noncontam_prev) == 0) # gives the number of cases

# Scale reads to even depth 
phy_prune_noncontam_prev <- prune_samples(sample_sums(phy_prune_noncontam_prev)>500, phy_prune_noncontam_prev)
smin <- min(sample_sums(phy_prune_noncontam_prev))
phy_prune_noncontam_prev

phy_prune_noncontam_freq <- prune_samples(sample_sums(phy_prune_noncontam_freq)>500, phy_prune_noncontam_freq)
smin <- min(sample_sums(phy_prune_noncontam_freq))
phy_prune_noncontam_freq

# Ordinate PCoA
phy_prune_noncontam_prev_pcoa <- ordinate(
  physeq = phy_prune_noncontam_prev, 
  method = "PCoA", 
  distance = "bray"
)

phy_prune_noncontam_freq_pcoa <- ordinate(
  physeq = phy_prune_noncontam_freq, 
  method = "PCoA", 
  distance = "bray"
)

# Plot PCoA 
phyloseq::plot_ordination(
  physeq = phy_prune_noncontam_prev,
  ordination = phy_prune_noncontam_prev_pcoa,
  color = "Habitat",
  shape = "Sample.Type",
  title = "PCoA of Shipley Skinner bacterial Communities") + 
  scale_color_manual(values = phylum_colors <- c(
    "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
             "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
             "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861")) +
               geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 

phyloseq::plot_ordination(
  physeq = phy_prune_noncontam_freq,
  ordination = phy_prune_noncontam_freq_pcoa,
  color = "Habitat",
  shape = "Sample.Type",
  title = "Filtered by frequency - ASV level") + 
  scale_color_manual(values = phylum_colors <- c(
    "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
             "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
             "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861")) +
               scale_shape_manual(values = c(15, 17, 19, 23, 7, 8, 3)) +
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("NMDS2") + # add the title on y axis
  #xlab("NMDS1") + # add the title on x axis
  ggtitle("Filtered by frequency - ASV level") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  #annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
  #annotate("text", x = 2.0, y = 2.0, label ="2D Stress: 0.11") +
  ggsave("16S_Soil_ShSk_ASV_PCoA_filt_freq.tiff", width = 8, height = 6, dpi = 300)

#subset by sample types, e.g., only soil samples
phy_prune_noncontam_freq_true <- subset_samples(phy_prune_noncontam_freq, Sample.Control =="True.Sample") %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_prune_noncontam_freq_true

# Check for ASVs that have no counted reads
any(taxa_sums(phy_prune_noncontam_freq_true) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_prune_noncontam_freq_true) == 0) # gives the number of cases

# Ordinate PCoA
phy_prune_noncontam_freq_true_pcoa <- ordinate(
  physeq = phy_prune_noncontam_freq_true, 
  method = "PCoA", 
  distance = "bray"
)

# Plot PCoA
phyloseq::plot_ordination(
  physeq = phy_prune_noncontam_freq_true,
  ordination = phy_prune_noncontam_freq_true_pcoa,
  color = "Soil.Type",
  shape = "Sample.Type",
  title = "PCoA of Shipley Skinner bacterial Communities") + 
  geom_point(aes(color = Soil.Type), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) +
  scale_color_manual(values = phylum_colors <- c(
    "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
             "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
             "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861")) +
               theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("NMDS2") + # add the title on y axis
  #xlab("NMDS1") + # add the title on x axis
  ggtitle("Filtered by frequency - ASV level") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  #annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
  #annotate("text", x = 2.0, y = 2.0, label ="2D Stress: 0.11") +
  ggsave("16S_Soil_ShSk_ASV_PCoA_filt_freq_soil-worm_soil-type.tiff", width = 8, height = 6, dpi = 300)

################################################
#subset by sample types, e.g., only soil samples
phy_prune_noncontam_freq_true_soil <-subset_samples(phy_prune_noncontam_freq_true, Sample.Type =="Soil") %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_prune_noncontam_freq_true_soil
any(taxa_sums(phy_prune_noncontam_freq_true_soil) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_prune_noncontam_freq_true_soil) == 0) # gives the number of cases
any(sample_sums(phy_prune_noncontam_freq_true_soil) == 0) # if TRUE, there are samples with no reads, otherwise FALSE
sum(sample_sums(phy_prune_noncontam_freq_true_soil) == 0) # gives the number of cases

# Filter out Eukaryotes, Archaea, chloroplasts and mitochondria,
phy_prune_noncontam_freq_true_soil_filt <- phy_prune_noncontam_freq_true_soil %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Order   != "Chloroplast"
  )
phy_prune_noncontam_freq_true_soil_filt

################################################
phy_prune_nema <-subset_samples(phy_prune, Sample.Type ==c("SingleWorm")) %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_prune_nema

any(taxa_sums(phy_prune_nema) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_prune_nema) == 0) # gives the number of cases
any(sample_sums(phy_prune_nema) == 0) # if TRUE, there are samples with no reads, otherwise FALSE
sum(sample_sums(phy_prune_nema) == 0) # gives the number of cases

# Filter out Eukaryotes, Archaea, chloroplasts and mitochondria,
phy_prune_nema_filt <- phy_prune_nema %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Order   != "Chloroplast"
  )
phy_prune_nema_filt

# mean, max and min of sample read counts
smin <- min(sample_sums(phy_prune_nema_filt))
smean <- mean(sample_sums(phy_prune_nema_filt))
smax <- max(sample_sums(phy_prune_nema_filt))
smedian <- median(sample_sums(phy_prune_nema_filt))

# Scale reads to even depth 
phy_prune_nema_filt_scale <- prune_samples(sample_sums(phy_prune_nema_filt)>500, phy_prune_nema_filt)
smin <- min(sample_sums(phy_prune_nema_filt_scale))
phy_prune_nema_filt_scale

# Transform number of reads to relative abundance
phy_prune_nema_filt_scale_trans = transform_sample_counts(phy_prune_nema_filt_scale, function(x) 100 * x/sum(x))

samples_to_remove <- c("Nem.57.SingleWorm.SK.08", "Nem.42.SingleWorm.SK.08", "Nem.118.SingleWorm.SK.02",
                       "Nem.114.SingleWorm.SK.02", "Nem.503.SingleWorm.SK.11")
subset_samples(phy_prune_nema_filt_scale_trans, Description %in% samples_to_remove)

phy_prune_nema_filt_scale_trans_sample_removed = subset_samples(
  phy_prune_nema_filt_scale_trans, !(Description %in% samples_to_remove))


# Ordinate NMDS
set.seed(1)
phy_prune_nema_filt_scale_trans_sample_removed_nmds <- ordinate(
  physeq = phy_prune_nema_filt_scale_trans_sample_removed, 
  method = "NMDS", 
  distance = "bray"
)

nb.cols <- 9
class_colors4 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

# Plot NMDS
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_prune_nema_filt_scale_trans_sample_removed,
  ordination = phy_prune_nema_filt_scale_trans_sample_removed_nmds,
  type = "samples",
  color = "Nematode.Order",
  shape = "Habitat") +
  scale_color_manual(values = class_colors4,
                     name = "Nematode Order") +
  #scale_shape_manual(values = c(19,17),
  #                  name = "Soil Type",
  #                 breaks = c("Compact.Soil", "Not.Compact.Soil"),
  #                labels = c("Compacted", "Uncompacted"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Nematode.Order), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("NMDS2") + # add the title on y axis
  #xlab("NMDS1") + # add the title on x axis
  ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  #annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
  annotate("text", x = 2.0, y = 2.0, label ="2D Stress: 0.11") +
  ggsave("16S_Nema_ShSk_ASV_NMDS_greater_500.tiff", width = 12, height = 8, dpi = 300)

################################################
phy_prune_soil_nema <-subset_samples(phy_prune, Sample.Type ==c("Soil", "SingleWorm")) %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_prune_soil_nema

any(taxa_sums(phy_prune_soil_nema) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_prune_soil_nema) == 0) # gives the number of cases
any(sample_sums(phy_prune_soil_nema) == 0) # if TRUE, there are samples with no reads, otherwise FALSE
sum(sample_sums(phy_prune_soil_nema) == 0) # gives the number of cases

# Filter out Eukaryotes, Archaea, chloroplasts and mitochondria,
phy_prune_soil_nema_filt <- phy_prune_soil_nema %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Order   != "Chloroplast"
  )
phy_prune_soil_nema_filt

# Melt to long format (for ggploting); prune out phyla below 1% in each sample
phy_prune_soil_nema_filt_phylum <- phy_prune_soil_nema_filt %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum
phy_prune_soil_nema_filt_phylum

phy_prune_soil_nema_filt_class <- phy_prune_soil_nema_filt %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum
phy_prune_soil_nema_filt_class


# Define the number of colors you want
nb.cols <- 22
class_colors2 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)


# Plot bar chart 16S Soil samples per habitat - Phylum level
ggplot(phy_prune_soil_nema_filt_phylum, aes(x = Sample.Type, y = Abundance, fill = Phylum)) + #plotting by sample
  facet_grid(. ~ Habitat) + # By Habitat
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = class_colors2) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.0, 0.0)) + # removes the internal margins
  #scale_x_discrete(
  #  breaks = c("SingleWorm, Soil"),
  # labels = c("Single Worm, Soil"),
  # expand = c(0.0, 0.0),
  # drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (Phyla > 1%)") + # add the title on y axis
  xlab("Sample Type") + # add the title on x axis
  ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  ggsave("16S_Phylum_Soil_Nema_ShSk_barplot.tiff", width = 10, height = 4, dpi = 300) # save graphic

# Define the number of colors you want
nb.cols <- 45
class_colors3 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

ggplot(phy_prune_soil_nema_filt_class, aes(x = Sample.Type, y = Abundance, fill = Class)) + #plotting by sample
  facet_grid(. ~ Habitat) + # By Habitat
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = class_colors3) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.0, 0.0)) + # removes the internal margins
  #scale_x_discrete(
  #  breaks = c("SingleWorm, Soil"),
  # labels = c("Single Worm, Soil"),
  # expand = c(0.0, 0.0),
  # drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (Phyla > 1%)") + # add the title on y axis
  xlab("Sample Type") + # add the title on x axis
  ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  ggsave("16S_Class_Soil_Nema_ShSk_barplot.tiff", width = 12, height = 4, dpi = 300) # save graphic

################################################
phy_soil <- subset_samples(phy, Sample.Type =="Soil") %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_soil

# Check for ASVs that have no counted reads
any(taxa_sums(phy_soil) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_soil) == 0) # gives the number of cases
any(sample_sums(phy_soil) == 0) # if TRUE, there are samples with no reads, otherwise FALSE
sum(sample_sums(phy_soil) == 0) # gives the number of cases

# Filter out Eukaryotes, Archaea, chloroplasts and mitochondria,
phy_soil_filt <- phy_soil %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Order   != "Chloroplast"
  )
phy_soil_filt

################################################
phy_prune_filt <- phy_prune_filt %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Order   != "Chloroplast"
  )
phy_prune_filt

# mean, max and min of sample read counts
smin <- min(sample_sums(phy_soil_filt))
smean <- mean(sample_sums(phy_soil_filt))
smax <- max(sample_sums(phy_soil_filt))
smedian <- median(sample_sums(phy_soil_filt))

# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(phy_soil_filt))

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
################################################
# Melt to long format (for ggploting); prune out phyla below 1% in each sample
phy_soil_filt_phylum <- phy_soil_filt %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# Set colors for plotting
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
           "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
           "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)

# Plot bar chart 16S Soil samples per habitat - Phylum level
ggplot(phy_soil_filt_phylum, aes(x = Habitat, y = Abundance, fill = Phylum)) + #plotting by sample
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = phylum_colors) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.0, 0.0)) + # removes the internal margins
  scale_x_discrete(
    breaks = c("Chaparral", "Coastal.scrub", "Hollyleaf.cherry", "Native.grass", "Oak.wood", "Riparian"),
    labels = c("Chaparral", "Coastal Scrub", "Hollyleaf Cherry", "Native Grass", "Oak Wood", "Riparian"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (Phyla > 1%)") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  ggsave("16S_Phylum_Soil_ShSk.tiff", width = 8, height = 4, dpi = 300) # save graphic

# Plot bar chart 16S Soil samples per Soil Type - Phylum level
ggplot(phy_soil_filt_phylum, aes(x = Soil.Type, y = Abundance, fill = Phylum)) + #plotting by sample
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = phylum_colors) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.0, 0.0)) + # removes the internal margins
  scale_x_discrete(
    breaks = c("Compact.Soil", "Not.Compact.Soil"),
    labels = c("Compacted", "Uncompacted"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (Phyla > 1%)") + # add the title on y axis
  xlab("Soil Type") + # add the title on x axis
  ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  ggsave("16S_Phylum_Soil-Type_ShSk.tiff", width = 8, height = 4, dpi = 300) # save graphic

################################################
# Melt to long format (for ggploting); prune out class below 1% in each sample
phy_soil_filt_class <- phy_soil_filt %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Class)                                      # Sort data frame alphabetically by Class

get_taxa_unique(phy_soil_filt_class)
# Set colors for plotting
class_colors <- viridis(28, option = "D")

# Define the number of colors you want
nb.cols <- 28
class_colors2 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

# Plot bar chart 16S Soil samples per habitat - Class level
ggplot(phy_soil_filt_class, aes(x = Habitat, y = Abundance, fill = Class)) + #plotting by sample
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  #scale_fill_viridis(discrete = FALSE, option = "D") +
  scale_fill_manual(values = class_colors2) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.0, 0.0)) + # removes the internal margins
  scale_x_discrete(
    breaks = c("Chaparral", "Coastal.scrub", "Hollyleaf.cherry", "Native.grass", "Oak.wood", "Riparian"),
    labels = c("Chaparral", "Coastal Scrub", "Hollyleaf Cherry", "Native Grass", "Oak Wood", "Riparian"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (Class > 1%)") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  ggsave("16S_Class_Soil-Habitat_ShSk.tiff", width = 8, height = 4, dpi = 300) # save graphic

# Plot bar chart 16S Soil samples per Soil Type - Class level
ggplot(phy_soil_filt_class, aes(x = Soil.Type, y = Abundance, fill = Class)) + #plotting by sample
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = class_colors2) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.0, 0.0)) + # removes the internal margins
  scale_x_discrete(
    breaks = c("Compact.Soil", "Not.Compact.Soil"),
    labels = c("Compacted", "Uncompacted"),
    expand = c(0.0, 0.0),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance (Class > 1%)") + # add the title on y axis
  xlab("Soil Type") + # add the title on x axis
  ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  ggsave("16S_Class_Soil-Type_ShSk.tiff", width = 8, height = 4, dpi = 300) # save graphic

########################### Ordinations ###########################

# Scale reads to even depth 
phy_soil_filt_prune <- prune_samples(sample_sums(phy_soil_filt)>500, phy_soil_filt)
smin <- min(sample_sums(phy_soil_filt_prune))
phy_soil_filt_prune

# Ordinate PCoA
phy_soil_filt_prune_pcoa <- ordinate(
  physeq = phy_soil_filt_prune, 
  method = "PCoA", 
  distance = "bray"
)

# Ordinate NMDS
set.seed(1)
phy_soil_filt_prune_nmds <- ordinate(
  physeq = phy_soil_filt_prune, 
  method = "NMDS", 
  distance = "bray"
)
phy_soil_filt_prune_nmds

set.seed(1)
phy_prune_noncontam_freq_true_soil_filt_nmds <- ordinate(
  physeq = phy_prune_noncontam_freq_true_soil_filt, 
  method = "NMDS", 
  distance = "bray"
)
phy_prune_noncontam_freq_true_soil_filt_nmds

# Plot PCoA 
plot_ordination(
  physeq = phy_soil_filt_prune,
  ordination = phy_soil_filt_prune_pcoa,
  color = "Habitat",
  shape = "Soil.Type",
  title = "PCoA of Shipley Skinner bacterial Communities") + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19",
                                         "#4daf4a", "#1919ff", "darkorchid3", "magenta")) +
                                           geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 

# Plot NMDS
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_soil_filt_prune,
  ordination = phy_soil_filt_prune_nmds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#a65628", "red", "#ffae19",
                                         "#4daf4a", "#1919ff", "darkorchid3"),
                                         name = "Habitat",
                     breaks = c("Chaparral", "Coastal.scrub", "Hollyleaf.cherry", "Native.grass", "Oak.wood", "Riparian"),
                     labels = c("Chaparral", "Coastal Scrub", "Hollyleaf Cherry", "Native Grass", "Oak Wood", "Riparian"))+
  scale_shape_manual(values = c(19,17),
                     name = "Soil Type",
                     breaks = c("Compact.Soil", "Not.Compact.Soil"),
                     labels = c("Compacted", "Uncompacted"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("NMDS2") + # add the title on y axis
  #xlab("NMDS1") + # add the title on x axis
  ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  #annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
  annotate("text", x = 2.0, y = 2.0, label ="2D Stress: 0.11") +
  ggsave("16S_Soil_ShSk_ASV_NMDS_greater_500.tiff", width = 6, height = 4, dpi = 300)

# Plot NMDS
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_prune_noncontam_freq_true_soil_filt,
  ordination = phy_prune_noncontam_freq_true_soil_filt_nmds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#a65628", "red", "#ffae19",
                                         "#4daf4a", "#1919ff", "darkorchid3"),
                                         name = "Habitat",
                     breaks = c("Chaparral", "Coastal.scrub", "Hollyleaf.cherry", "Native.grass", "Oak.wood", "Riparian"),
                     labels = c("Chaparral", "Coastal Scrub", "Hollyleaf Cherry", "Native Grass", "Oak Wood", "Riparian"))+
  scale_shape_manual(values = c(19,17),
                     name = "Soil Type",
                     breaks = c("Compact.Soil", "Not.Compact.Soil"),
                     labels = c("Compacted", "Uncompacted"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("NMDS2") + # add the title on y axis
  #xlab("NMDS1") + # add the title on x axis
  ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  #annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
  annotate("text", x = 2.0, y = 2.0, label ="2D Stress: 0.09") +
  ggsave("16S_Soil_ShSk_ASV_NMDS_filt_freq.tiff", width = 6, height = 4, dpi = 300)

set.seed(1)
# Calculate bray curtis distance matrix
phy_soil_filt_prune_bray <- phyloseq::distance(phy_soil_filt_prune, method = "bray")

# make a data frame from the sample_data
sampledf_soil <- data.frame(sample_data(phy_soil_filt_prune))

# Adonis test
adonis(phy_soil_filt_prune_bray ~ Habitat, data = sampledf_soil)
calc_pairwise_permanovas(phy_soil_filt_prune_bray, sampledf_soil, "Habitat")
adonis(phy_soil_filt_prune_bray ~ Soil.Type, data = sampledf_soil)
adonis(phy_soil_filt_prune_bray ~ Sample.Site, data = sampledf_soil)

# Homogeneity of dispersion test
beta_habitat <- betadisper(phy_soil_filt_prune_bray, sampledf_soil$Habitat)
permutest(beta_habitat)
beta_soil_type <- betadisper(phy_soil_filt_prune_bray, sampledf_soil$Soil.Type)
permutest(beta_soil_type)
beta_sample_site <- betadisper(phy_soil_filt_prune_bray, sampledf_soil$Sample.Site)
permutest(beta_sample_site)

# Check again for ASVs that have no counted reads
any(taxa_sums(phy_soil_filt) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_soil_prune) == 0) # gives the number of cases

# Transform number of reads to relative abundance and remove NAs
phy_soil_prune_trans = transform_sample_counts(phy_soil_prune, function(x) 100 * x/sum(x))
phy_soil_prune_trans2 <- phyloseq_rm_na_tax(phy_soil_prune)
any(taxa_sums(phy_soil_prune_trans2) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_soil_prune_trans2) == 0) # gives the number of cases

#plot simple barplot
plot_bar(phy_soil, "Habitat", fill = "Kingdom") +
  geom_bar(aes(color=Kingdom, fill=Kingdom), stat="identity", position="stack")


#Create NMDS
phy_soil_prune_ord <- ordinate(phy_soil_prune, "NMDS", "bray")
p1 = plot_ordination(GP1, GP.ord, type="taxa", color="Phylum", title="taxa")
print(p1)

####################### Additional code #######################
set.seed(28132)
phy_soil_prune_noncontam_comb05_true_raref = rarefy_even_depth(phy_soil_prune_noncontam_comb05_true, sample.size = 4223)

phy_soil_prune_noncontam_comb05_true_prop = transform_sample_counts(phy_soil_prune_noncontam_comb05_true, function(x) 500 * x/sum(x))

par(mfrow = c(1, 2))
title = "Sum of reads for each sample, Raref"
plot(sort(sample_sums(phy_soil_prune_noncontam_comb05_true_raref), TRUE), type = "h", main = title, ylab = "reads", 
     ylim = c(0, 1000))
title = "Sum of reads for each sample, Prop"
plot(sort(sample_sums(phy_soil_prune_noncontam_comb05_true_prop), TRUE), type = "h", main = title, ylab = "reads", 
     ylim = c(0, 1000))

sample_variables(phy_soil_prune_noncontam_comb05_true_raref)
top25otus = names(sort(taxa_sums(phy_soil_prune_noncontam_comb05_true_raref), TRUE)[1:25])
taxtab25 = cbind(tax_table(phy_soil_prune_noncontam_comb05_true_raref), genus_25 = NA)
taxtab25[top25otus, "genus_25"] <- as(tax_table(phy_soil_prune_noncontam_comb05_true_raref)[top25otus, "Genus"],
                                      "character")
tax_table(phy_soil_prune_noncontam_comb05_true_raref) <- tax_table(taxtab25)

phy_soil_prune_noncontam_comb05_true_raref <- transform_sample_counts(phy_soil_prune_noncontam_comb05_true_raref, function(x) 100 * x/sum(x))

title = "Relative Abundance"
genus_plot2 <- plot_bar(phy_soil_prune_noncontam_comb05_true_raref, "Habitat", fill = "genus_25", title = title)+
  facet_wrap(~Soil.Type, scales="free_x")
print(genus_plot2)


top19asvs = names(sort(taxa_sums(phy_soil_prune_noncontam_comb05_true_raref), TRUE)[1:19])
taxtab19 = cbind(tax_table(phy_soil_prune_noncontam_comb05_true_raref), family19 = NA)
taxtab19[top19asvs, "family19"] <- as(tax_table(phy_soil_prune_noncontam_comb05_true_raref)[top19asvs, "Family"], 
                                      "character")
tax_table(phy_soil_prune_noncontam_comb05_true_raref) <- tax_table(taxtab19)

title = "Figure 1 Part A (remake), attempt 1"
plot_bar(phy_soil_prune_noncontam_comb05_true_raref, "Habitat", fill = "family19", title = title) + coord_flip()

phy_soil_prune_noncontam_comb05_true_rarefm = merge_samples(phy_soil_prune_noncontam_comb05_true_raref, "Habitat")
sample_data(phy_soil_prune_noncontam_comb05_true_rarefm)$Habitat <- levels(sample_data(phy_soil_prune_noncontam_comb05_true_raref)$Habitat)
phy_soil_prune_noncontam_comb05_true_rarefm = transform_sample_counts(phy_soil_prune_noncontam_comb05_true_rarefm, function(x) 100 * x/sum(x))

title = "Figure 1 Part A (remake), attempt 2"
plot_bar(phy_soil_prune_noncontam_comb05_true_rarefm, "Sample", fill = "family19", title = title) + coord_flip() + 
  ylab("Percentage of Sequences")

phy_soil_prune_noncontam_comb05_true_rarefm19 = prune_taxa(top19asvs, phy_soil_prune_noncontam_comb05_true_rarefm)
title = "Figure 1 Part A (remake), attempt 3"
plot_bar(phy_soil_prune_noncontam_comb05_true_rarefm19, "Sample", fill = "family19", title = title) + coord_flip() + 
  ylab("Percentage of Sequences") + ylim(0, 100)


p <- plot_richness(phy_soil_prune_noncontam_comb05_true, "Habitat", "Soil.Type")
(p <- p + geom_boxplot(data = p$data, aes(x = Habitat, y = value, color = NULL), 
                       alpha = 0.1))
####################### Additional code #######################