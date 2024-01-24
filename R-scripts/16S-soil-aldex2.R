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

# Import phyloseq object scaled (samples > 500 reads)
phy_16_prune_soil_noncontam_prev05_true_filt_scale <- readRDS("data/phy_16_prune_soil_noncontam_prev05_true_filt_scale.RDS")
class(phy_16_prune_soil_noncontam_prev05_true_filt_scale)
phy_16_prune_soil_noncontam_prev05_true_filt_scale # 9622 taxa and 50 samples

######## Aldex2 analysis on phyoseq object ASV level ############
aldex2_phy_soil_asv <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale)),
                                phyloseq::sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale)$Soil.Type,
                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2_phy_soil_asv_table <- data.frame(phyloseq::otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale))
aldex2_phy_soil_asv_table <- rownames_to_column(aldex2_phy_soil_asv_table, var = "OTU")
write.csv(aldex2_phy_soil_asv_table, "aldex2/aldex2_phy_soil_asv_table.csv")

# Write aldex2 results to file as csv
write.csv(aldex2_phy_soil_asv, "aldex2/aldex2_phy_soil_asv.csv")

# Import aldex2 results and rename the X variable by OTU
aldex2_phy_soil_asv_result <- read.csv("aldex2/aldex2_phy_soil_asv.csv")
colnames(aldex2_phy_soil_asv_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

#Clean up presentation
taxa_info <- data.frame(tax_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale))
taxa_info <- taxa_info %>%
  rownames_to_column(var = "OTU")
sample_tab_soil <- data.frame(sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale))
sample_soil_factors <- data.frame(cbind(sample_tab_soil[, c(21, 22, 36, 37)]))

# filter aldex2 results by sig kw.ep and join the taxanomic information
sig_aldex2_phy_soil_asv_result <- aldex2_phy_soil_asv_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_phy_soil_asv_result <- left_join(sig_aldex2_phy_soil_asv_result, taxa_info)
write.csv(sig_aldex2_phy_soil_asv_result, "aldex2/sig_aldex2_phy_soil_asv_result.csv") # 432 taxa with significant differential abundance

# Create clr objects by using OTU ids and original otu table.
# All significant OTUs
sig_aldex2_phy_soil_asv_result_count <- left_join(sig_aldex2_phy_soil_asv_result, aldex2_phy_soil_asv_table)
write.csv(sig_aldex2_phy_soil_asv_result_count, "aldex2/sig_aldex2_phy_soil_asv_result_count.csv")

clr_soil_asv <- sig_aldex2_phy_soil_asv_result_count[, -(2:9)]
rownames(clr_soil_asv) <- clr_soil_asv$OTU
clr_soil_asv <- clr_soil_asv[, -1]

# Adjusting zeros on the matrix and applying log transformation
clr_soil_asv_czm <- cmultRepl(t(clr_soil_asv),  label=0, method="CZM")
clr_soil_asv_czm_tv <- t(apply(clr_soil_asv_czm, 1, function(x){log(x) - mean(log(x))}))
clr_soil_asv_czm <- (apply(clr_soil_asv, 1, function(x){log(x+1) - mean(log(x+1))}))

# Define palette color
col_matrix <- colorRampPalette(brewer.pal(10, "BrBG"))(256)
col_matrix2 <- colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_matrix3 <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_matrix4 <- colorRamp2(c(-2, 0, 2), c("blue", "yellow", "red"))

# Combine the heatmap and the annotation
Z.Score.clr_soil_asv <- scale(t(clr_soil_asv_czm))
heatmap(Z.Score.clr_soil_asv) # simple heatmap

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
      labels = c ("CHA", "CSS", "NGR", "HLC", "OWL", "RIP")
    )
  ))

# Organize total abundance and taxa name
# Taxa name
Z.Score.clr_soil_asv_name <- as.data.frame(Z.Score.clr_soil_asv)
str(Z.Score.clr_soil_asv)
Z.Score.clr_soil_asv_name <- rownames_to_column(Z.Score.clr_soil_asv_name, var = "OTU")
Z.Score.clr_soil_asv_name <- left_join(Z.Score.clr_soil_asv_name, taxa_info)
head(Z.Score.clr_soil_asv_name)
Z.Score.clr_soil_asv_name <- Z.Score.clr_soil_asv_name[, -(2:51)]

# Taxa abundance
aldex2_phy_soil_asv_table_total <- as.data.frame(aldex2_phy_soil_asv_table)
aldex2_phy_soil_asv_table_total$Total <- rowSums(aldex2_phy_soil_asv_table[, -1])
head(aldex2_phy_soil_asv_table_total)
aldex2_phy_soil_asv_table_total <- aldex2_phy_soil_asv_table_total[, -(2:51)]

Z.Score.clr_soil_asv_count_total <- left_join(Z.Score.clr_soil_asv_name, aldex2_phy_soil_asv_table_total)
head(Z.Score.clr_soil_asv_count_total)

# Defining color scheme for row annotations
abundance_col_fun = colorRamp2(c(0, 5000, 10000, 15000, 20000),
                               c("white", "#c7eae5", "#80cdc1","#35978f","#01665e"))

ha_right_phy = rowAnnotation(
  Abundance = Z.Score.clr_soil_asv_count_total$Total, border = FALSE, col = list(Abundance = abundance_col_fun))
row_labels_phy = Z.Score.clr_soil_asv_count_total$OTU

row_labels_phy = Z.Score.clr_soil_asv_count_total$OTU

# Plot heatmap at the phylum level
hm_soil_asv <- Heatmap(Z.Score.clr_soil_asv, name = "Z-score, CLR", col = col_matrix,
                         column_title = "16S rRNA soil microbiome", 
                         column_title_gp = gpar(fontface = "bold", fontsize = 14),
                         column_split = as.vector(as.vector(sample_soil_factors$Soil.Type)),
                         #column_order = order(as.numeric(gsub("column", "", colnames(Z.Score.phylum_hab)))),
                         border = TRUE,
                         top_annotation = ha_top,
                         right_annotation = ha_right_phy,
                         row_title = "ASV",
                         row_labels = row_labels_phy,
                         row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                         row_names_gp = gpar(fontsize = 6),
                         column_names_gp = gpar(fontsize = 6),
                         row_order = order(row_labels_phy),
                         rect_gp = gpar(col = "white", lwd = 1),
                         show_column_names = TRUE,
                         show_heatmap_legend = TRUE)
hm_soil_asv

# Collapse tax table according to taxonomic rank for analysis with ALDEx2
phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank1 <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale, taxrank = "Phylum")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank2 <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale, taxrank = "Class")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank3 <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale, taxrank = "Order")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank4 <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale, taxrank = "Family")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank5 <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale, taxrank = "Genus")

######## Aldex2 analysis on collapsed phyoseq objects############
aldex2_phy_hab <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank1)),
                                phyloseq::sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank1)$Habitat,
                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2_phy_soil <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank1)),
                                 phyloseq::sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank1)$Soil.Type,
                                 mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)

rank1_otu_table <- data.frame(phyloseq::otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank1))
rank1_otu_table <- rownames_to_column(rank1_otu_table, var = "OTU")
write.csv(rank1_otu_table, "aldex2/rank1_otu_table.csv")

aldex2_class_hab <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank2)),
                                  phyloseq::sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank2)$Habitat,
                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2_class_soil <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank2)),
                                   phyloseq::sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank2)$Soil.Type,
                                   mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
rank2_otu_table <- data.frame(phyloseq::otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank2))
rank2_otu_table <- rownames_to_column(rank2_otu_table, var = "OTU")
write.csv(rank2_otu_table, "aldex2/rank2_otu_table.csv")

aldex2_order_hab <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank3)),
                                  phyloseq::sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank3)$Habitat,
                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2_order_soil <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank3)),
                                   phyloseq::sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank3)$Soil.Type,
                                   mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
rank3_otu_table <- data.frame(phyloseq::otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank3))
rank3_otu_table <- rownames_to_column(rank3_otu_table, var = "OTU")
write.csv(rank3_otu_table, "aldex2/rank3_otu_table.csv")

aldex2_fam_hab <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank4)),
                                phyloseq::sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank4)$Habitat,
                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2_fam_soil <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank4)),
                                 phyloseq::sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank4)$Soil.Type,
                                 mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
rank4_otu_table <- data.frame(phyloseq::otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank4))
rank4_otu_table <- rownames_to_column(rank4_otu_table, var = "OTU")
write.csv(rank4_otu_table, "aldex2/rank4_otu_table.csv")

aldex2_genus_hab <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank5)),
                                  phyloseq::sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank5)$Habitat,
                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex2_genus_soil <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank5)),
                                   phyloseq::sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank5)$Soil.Type,
                                   mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)

rank5_otu_table <- data.frame(phyloseq::otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_rank5))
rank5_otu_table <- rownames_to_column(rank5_otu_table, var = "OTU")
write.csv(rank5_otu_table, "aldex2/rank5_otu_table.csv")

# Write aldex2 results to file as csv
# Habitat comparison results
write.csv(aldex2_phy_hab, "aldex2/aldex2_phy_hab.csv")
write.csv(aldex2_class_hab, "aldex2/aldex2_class_hab.csv")
write.csv(aldex2_order_hab, "aldex2/aldex2_order_hab.csv")
write.csv(aldex2_fam_hab, "aldex2/aldex2_fam_hab.csv")
write.csv(aldex2_genus_hab, "aldex2/aldex2_genus_hab.csv")

# Soil comparison results
write.csv(aldex2_phy_soil, "aldex2/aldex2_phy_soil.csv")
write.csv(aldex2_class_soil, "aldex2/aldex2_class_soil.csv")
write.csv(aldex2_order_soil, "aldex2/aldex2_order_soil.csv")
write.csv(aldex2_fam_soil, "aldex2/aldex2_fam_soil.csv")
write.csv(aldex2_genus_soil, "aldex2/aldex2_genus_soil.csv")

# Import aldex2 results and rename the X variable by OTU
aldex2_phy_hab_result <- read.csv("aldex2/aldex2_phy_hab.csv")
colnames(aldex2_phy_hab_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")
aldex2_phy_soil_result <- read.csv("aldex2/aldex2_phy_soil.csv")
colnames(aldex2_phy_soil_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

aldex2_class_hab_result <- read.csv("aldex2/aldex2_class_hab.csv")
colnames(aldex2_class_hab_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")
aldex2_class_soil_result <- read.csv("aldex2/aldex2_class_soil.csv")
colnames(aldex2_class_soil_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

aldex2_order_hab_result <- read.csv("aldex2/aldex2_order_hab.csv")
colnames(aldex2_order_hab_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")
aldex2_order_soil_result <- read.csv("aldex2/aldex2_order_soil.csv")
colnames(aldex2_order_soil_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

aldex2_fam_hab_result <- read.csv("aldex2/aldex2_fam_hab.csv")
colnames(aldex2_fam_hab_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")
aldex2_fam_soil_result <- read.csv("aldex2/aldex2_fam_soil.csv")
colnames(aldex2_fam_soil_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

aldex2_genus_hab_result <- read.csv("aldex2/aldex2_genus_hab.csv")
colnames(aldex2_genus_hab_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")
aldex2_genus_soil_result <- read.csv("aldex2/aldex2_genus_soil.csv")
colnames(aldex2_genus_soil_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

#Clean up presentation
taxa_info <- data.frame(tax_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale))
taxa_info <- taxa_info %>%
  rownames_to_column(var = "OTU")
sample_tab_soil <- data.frame(sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale))
sample_soil_factors <- data.frame(cbind(sample_tab_soil[, c(21, 22, 36, 37)]))

# filter aldex2 results by sig kw.ep and join the taxanomic information
sig_aldex2_phy_hab_result <- aldex2_phy_hab_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_phy_hab_result <- left_join(sig_aldex2_phy_hab_result, taxa_info)
write.csv(sig_aldex2_phy_hab_result, "aldex2/sig_aldex2_phy_hab_result.csv") # 8 taxa with significant differential abundance

sig_aldex2_phy_soil_result <- aldex2_phy_soil_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_phy_soil_result <- left_join(sig_aldex2_phy_soil_result, taxa_info)
write.csv(sig_aldex2_phy_soil_result, "aldex2/sig_aldex2_phy_soil_result.csv") # 8 taxa with significant differential abundance

sig_aldex2_class_hab_result <- aldex2_class_hab_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_class_hab_result <- left_join(sig_aldex2_class_hab_result, taxa_info)
write.csv(sig_aldex2_class_hab_result, "aldex2/sig_aldex2_class_hab_result.csv") # 29 taxa with significant differential abundance

sig_aldex2_class_soil_result <- aldex2_class_soil_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_class_soil_result <- left_join(sig_aldex2_class_soil_result, taxa_info)
write.csv(sig_aldex2_class_soil_result, "aldex2/sig_aldex2_class_soil_result.csv") # 29 taxa with significant differential abundance

sig_aldex2_order_hab_result <- aldex2_order_hab_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_order_hab_result <- left_join(sig_aldex2_order_hab_result, taxa_info)
write.csv(sig_aldex2_order_hab_result, "aldex2/sig_aldex2_order_hab_result.csv") # 73 taxa with significant differential abundance

sig_aldex2_order_soil_result <- aldex2_order_soil_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_order_soil_result <- left_join(sig_aldex2_order_soil_result, taxa_info)
write.csv(sig_aldex2_order_soil_result, "aldex2/sig_aldex2_order_soil_result.csv") # 66 taxa with significant differential abundance

sig_aldex2_fam_hab_result <- aldex2_fam_hab_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_fam_hab_result <- left_join(sig_aldex2_fam_hab_result, taxa_info)
write.csv(sig_aldex2_fam_hab_result, "aldex2/sig_aldex2_fam_hab_result.csv") # 110 taxa with significant differential abundance

sig_aldex2_fam_soil_result <- aldex2_fam_soil_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_fam_soil_result <- left_join(sig_aldex2_fam_soil_result, taxa_info)
write.csv(sig_aldex2_fam_soil_result, "aldex2/sig_aldex2_fam_soil_result.csv") # 100 taxa with significant differential abundance

sig_aldex2_genus_hab_result <- aldex2_genus_hab_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_genus_hab_result <- left_join(sig_aldex2_genus_hab_result, taxa_info)
write.csv(sig_aldex2_genus_hab_result, "aldex2/sig_aldex2_genus_hab_result.csv") # 205 taxa with significant differential abundance

sig_aldex2_genus_soil_result <- aldex2_genus_soil_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_genus_soil_result <- left_join(sig_aldex2_genus_soil_result, taxa_info)
write.csv(sig_aldex2_genus_soil_result, "aldex2/sig_aldex2_genus_soil_result.csv") # 185 taxa with significant differential abundance

# Create clr objects by using OTU ids and original otu table, all significant OTUs and only top 30
sig_aldex2_phy_hab_result_count <- left_join(sig_aldex2_phy_hab_result, rank1_otu_table)
write.csv(sig_aldex2_phy_hab_result_count, "aldex2/sig_aldex2_phy_hab_result_count.csv")
sig_aldex2_phy_soil_result_count <- left_join(sig_aldex2_phy_soil_result, rank1_otu_table)
write.csv(sig_aldex2_phy_soil_result_count, "aldex2/sig_aldex2_phy_soil_result_count.csv")

clr_phy_hab <- sig_aldex2_phy_hab_result_count[, -(2:9)]
rownames(clr_phy_hab) <- clr_phy_hab$OTU
clr_phy_hab <- clr_phy_hab[, -1]

clr_phy_soil <- sig_aldex2_phy_soil_result_count[, -(2:9)]
rownames(clr_phy_soil) <- clr_phy_soil$OTU
clr_phy_soil <- clr_phy_soil[, -1]

sig_aldex2_class_hab_result_count <- left_join(sig_aldex2_class_hab_result, rank2_otu_table)
write.csv(sig_aldex2_class_hab_result_count, "aldex2/sig_aldex2_class_hab_result_count.csv")
sig_aldex2_class_soil_result_count <- left_join(sig_aldex2_class_soil_result, rank2_otu_table)
write.csv(sig_aldex2_class_soil_result_count, "aldex2/sig_aldex2_class_soil_result_count.csv")

clr_class_hab <- sig_aldex2_class_hab_result_count[, -(2:9)]
rownames(clr_class_hab) <- clr_class_hab$OTU
clr_class_hab <- clr_class_hab[, -1]
clr_class_soil <- sig_aldex2_class_soil_result_count[, -(2:9)]
rownames(clr_class_soil) <- clr_class_soil$OTU
clr_class_soil <- clr_class_soil[, -1]

sig_aldex2_order_hab_result_count <- left_join(sig_aldex2_order_hab_result, rank3_otu_table)
write.csv(sig_aldex2_order_hab_result_count, "aldex2/sig_aldex2_order_hab_result_count.csv")
sig_aldex2_order_soil_result_count <- left_join(sig_aldex2_order_soil_result, rank3_otu_table)
write.csv(sig_aldex2_order_soil_result_count, "aldex2/sig_aldex2_order_soil_result_count.csv")

clr_order_hab <- sig_aldex2_order_hab_result_count[, -(2:9)]
rownames(clr_order_hab) <- clr_order_hab$OTU
clr_order_hab <- clr_order_hab[, -1]
clr_order_soil <- sig_aldex2_order_soil_result_count[, -(2:9)]
rownames(clr_order_soil) <- clr_order_soil$OTU
clr_order_soil <- clr_order_soil[, -1]

sig_aldex2_fam_hab_result_count <- left_join(sig_aldex2_fam_hab_result, rank4_otu_table)
write.csv(sig_aldex2_fam_hab_result_count, "aldex2/sig_aldex2_fam_hab_result_count.csv")
sig_aldex2_fam_soil_result_count <- left_join(sig_aldex2_fam_soil_result, rank4_otu_table)
write.csv(sig_aldex2_fam_soil_result_count, "aldex2/sig_aldex2_fam_soil_result_count.csv")

clr_fam_hab <- sig_aldex2_fam_hab_result_count[, -(2:9)]
rownames(clr_fam_hab) <- clr_fam_hab$OTU
clr_fam_hab <- clr_fam_hab[, -1]
clr_fam_soil <- sig_aldex2_fam_soil_result_count[, -(2:9)]
rownames(clr_fam_soil) <- clr_fam_soil$OTU
clr_fam_soil <- clr_fam_soil[, -1]

sig_aldex2_genus_hab_result_count <- left_join(sig_aldex2_genus_hab_result, rank5_otu_table)
write.csv(sig_aldex2_genus_hab_result_count, "aldex2/sig_aldex2_genus_hab_result_count.csv")
sig_aldex2_genus_soil_result_count <- left_join(sig_aldex2_genus_soil_result, rank5_otu_table)
write.csv(sig_aldex2_genus_soil_result_count, "aldex2/sig_aldex2_genus_soil_result_count.csv")

clr_genus_hab <- sig_aldex2_genus_hab_result_count[, -(2:9)]
rownames(clr_genus_hab) <- clr_genus_hab$OTU
clr_genus_hab <- clr_genus_hab[, -1]
clr_genus_soil <- sig_aldex2_genus_soil_result_count[, -(2:9)]
rownames(clr_genus_soil) <- clr_genus_soil$OTU
clr_genus_soil <- clr_genus_soil[, -1]

# Adjusting zeros on the matrix and applying log transformation
shsk_phylum_czm_hab <- cmultRepl(t(clr_phy_hab),  label=0, method="CZM")
shsk_phylum_clr_hab <- t(apply(shsk_phylum_czm_hab, 1, function(x){log(x) - mean(log(x))}))
shsk_phylum_czm_soil <- cmultRepl(t(clr_phy_soil),  label=0, method="CZM")
shsk_phylum_clr_soil <- t(apply(shsk_phylum_czm_soil, 1, function(x){log(x) - mean(log(x))}))

shsk_class_czm_hab <- cmultRepl(t(clr_class_hab),  label=0, method="CZM")
shsk_class_clr_hab <- t(apply(shsk_class_czm_hab, 1, function(x){log(x) - mean(log(x))}))
shsk_class_czm_soil <- cmultRepl(t(clr_class_soil),  label=0, method="CZM")
shsk_class_clr_soil <- t(apply(shsk_class_czm_soil, 1, function(x){log(x) - mean(log(x))}))

shsk_order_czm_hab <- cmultRepl(t(clr_order_hab),  label=0, method="CZM")
shsk_order_clr_hab <- t(apply(shsk_order_czm_hab, 1, function(x){log(x) - mean(log(x))}))
shsk_order_czm_soil <- cmultRepl(t(clr_order_soil),  label=0, method="CZM")
shsk_order_clr_soil <- t(apply(shsk_order_czm_soil, 1, function(x){log(x) - mean(log(x))}))

shsk_fam_czm_hab <- cmultRepl(t(clr_fam_hab),  label=0, method="CZM")
shsk_fam_clr_hab <- t(apply(shsk_fam_czm_hab, 1, function(x){log(x) - mean(log(x))}))
shsk_fam_czm_soil <- cmultRepl(t(clr_fam_soil),  label=0, method="CZM")
shsk_fam_clr_soil <- t(apply(shsk_fam_czm_soil, 1, function(x){log(x) - mean(log(x))}))

shsk_genus_czm_hab <- cmultRepl(t(clr_genus_hab),  label=0, method="CZM")
shsk_genus_clr_hab <- t(apply(shsk_genus_czm_hab, 1, function(x){log(x) - mean(log(x))}))
shsk_genus_czm_soil <- cmultRepl(t(clr_genus_soil),  label=0, method="CZM")
shsk_genus_clr_soil <- t(apply(shsk_genus_czm_soil, 1, function(x){log(x) - mean(log(x))}))

# Define palette color
col_matrix <- colorRampPalette(brewer.pal(10, "BrBG"))(256)
col_matrix2 <- colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_matrix3 <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_matrix4 <- colorRamp2(c(-2, 0, 2), c("blue", "yellow", "red"))

# Combine the heatmap and the annotation
Z.Score.phylum_hab <- scale(t(shsk_phylum_clr_hab))
Z.Score.phylum_soil <- scale(t(shsk_phylum_clr_soil))

Z.Score.class_hab <- scale(t(shsk_class_clr_hab))
Z.Score.class_soil <- scale(t(shsk_class_clr_soil))

Z.Score.order_hab <- scale(t(shsk_order_clr_hab))
Z.Score.order_soil <- scale(t(shsk_order_clr_soil))

Z.Score.fam_hab <- scale(t(shsk_fam_clr_hab))
Z.Score.fam_soil <- scale(t(shsk_fam_clr_soil))

Z.Score.genus_hab <- scale(t(shsk_genus_clr_hab))
Z.Score.genus_soil <- scale(t(shsk_genus_clr_soil))

# Define colors for each level of qualitative variables, i.e. soil type and habitats
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
      labels = c ("CHA", "CSS", "NGR", "HLC", "OWL", "RIP")
    )
  ))

# Organize total abundance and taxa name
Z.Score.phylum_hab_name <- as.data.frame(Z.Score.phylum_hab)
str(Z.Score.phylum_hab_name)
Z.Score.phylum_hab_name <- rownames_to_column(Z.Score.phylum_hab_name, var = "OTU")
Z.Score.phylum_hab_name <- left_join(Z.Score.phylum_hab_name, taxa_info)
head(Z.Score.phylum_hab_name)
Z.Score.phylum_hab_name <- Z.Score.phylum_hab_name[, -(2:51)]

rank1_otu_table_total <- as.data.frame(rank1_otu_table)
rank1_otu_table_total$Total <- rowSums(rank1_otu_table[, -1])
head(rank1_otu_table_total)
rank1_otu_table_total <- rank1_otu_table_total[, -(2:51)]

Z.Score.phylum_hab_count_total <- left_join(Z.Score.phylum_hab_name, rank1_otu_table_total)
head(Z.Score.phylum_hab_count_total)

# Defining color scheme for row annotations
abundance_col_fun = colorRamp2(c(0, 50000, 100000, 150000),
                               c("#c7eae5",
                                          "#80cdc1",
                                          "#35978f",
                                          "#01665e"))
                                          
ha_right_phy = rowAnnotation(
  Abundance = Z.Score.phylum_hab_count_total$Total, border = FALSE, col = list(Abundance = abundance_col_fun))
row_labels_phy = Z.Score.phylum_hab_count_total$Phylum

row_labels_phy = Z.Score.phylum_hab_count_total$Phylum

# Plot heatmap at the phylum level
hm_phylum_hab <- Heatmap(Z.Score.phylum_hab, name = "Z-score, CLR", col = col_matrix,
                         column_title = "16S rRNA soil microbiome", 
                         column_title_gp = gpar(fontface = "bold", fontsize = 14),
                         column_split = as.vector(as.vector(sample_soil_factors$Soil.Type)),
                         #column_order = order(as.numeric(gsub("column", "", colnames(Z.Score.phylum_hab)))),
                         border = TRUE,
                         top_annotation = ha_top,
                         right_annotation = ha_right_phy,
                         row_title = "Phylum",
                         row_labels = row_labels_phy,
                         row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                         row_names_gp = gpar(fontsize = 6),
                         column_names_gp = gpar(fontsize = 6),
                         row_order = order(row_labels_phy),
                         rect_gp = gpar(col = "white", lwd = 1),
                         show_column_names = TRUE,
                         show_heatmap_legend = TRUE)
hm_phylum_hab

# Plot heatmap at the class level
Z.Score.class_hab_name <- as.data.frame(Z.Score.class_hab)
str(Z.Score.class_hab_name)
Z.Score.class_hab_name <- rownames_to_column(Z.Score.class_hab_name, var = "OTU")
Z.Score.class_hab_name <- left_join(Z.Score.class_hab_name, taxa_info)
head(Z.Score.class_hab_name)
Z.Score.class_hab_name <- Z.Score.class_hab_name[, -(2:51)]

rank2_otu_table_total <- as.data.frame(rank2_otu_table)
rank2_otu_table_total$Total <- rowSums(rank2_otu_table[, -1])
head(rank2_otu_table_total)
rank2_otu_table_total <- rank2_otu_table_total[, -(2:51)]

Z.Score.class_hab_count_total <- left_join(Z.Score.class_hab_name, rank2_otu_table_total)
head(Z.Score.class_hab_count_total)

ha_right_class = rowAnnotation(
  Abundance = Z.Score.class_hab_count_total$Total, border = FALSE, col = list(Abundance = abundance_col_fun))
row_labels_class = Z.Score.class_hab_count_total$Class

hm_class <- Heatmap(Z.Score.class_hab, name = "Z-score, CLR", col = col_matrix,
                    column_title = "16S rRNA soil microbiome", 
                    column_title_gp = gpar(fontface = "bold", fontsize = 14),
                    column_split = as.vector(as.vector(sample_soil_factors$Soil.Type)),
                    border = TRUE,
                    top_annotation = ha_top,
                    right_annotation = ha_right_class,
                    row_title = "Class",
                    row_labels = row_labels_class,
                    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                    row_names_gp = gpar(fontsize = 6),
                    column_names_gp = gpar(fontsize = 6),
                    row_order = order(row_labels_class),
                    rect_gp = gpar(col = "white", lwd = 1),
                    show_column_names = FALSE,
                    show_heatmap_legend = TRUE)
hm_class

# Plot heatmap at the order level
Z.Score.order_hab_name <- as.data.frame(Z.Score.order_hab)
str(Z.Score.order_hab_name)
Z.Score.order_hab_name <- rownames_to_column(Z.Score.order_hab_name, var = "OTU")
Z.Score.order_hab_name <- left_join(Z.Score.order_hab_name, taxa_info)
head(Z.Score.order_hab_name)
Z.Score.order_hab_name <- Z.Score.order_hab_name[, -(2:51)]

rank3_otu_table_total <- as.data.frame(rank3_otu_table)
rank3_otu_table_total$Total <- rowSums(rank3_otu_table[, -1])
head(rank3_otu_table_total)
rank3_otu_table_total <- rank3_otu_table_total[, -(2:51)]

Z.Score.order_hab_count_total <- left_join(Z.Score.order_hab_name, rank3_otu_table_total)
head(Z.Score.order_hab_count_total)


ha_right_order = rowAnnotation(
  Abundance = Z.Score.order_hab_count_total$Total, border = FALSE, col = list(Abundance = abundance_col_fun))
row_labels_order = Z.Score.order_hab_count_total$Order

hm_order <- Heatmap(Z.Score.order_hab, name = "Z-score, CLR", col = col_matrix,
                    column_title = "16S rRNA soil microbiome", 
                    column_title_gp = gpar(fontface = "bold", fontsize = 14),
                    column_split = as.vector(as.vector(sample_soil_factors$Soil.Type)),
                    border = TRUE,
                    top_annotation = ha_top,
                    right_annotation = ha_right_order,
                    row_title = "Order",
                    row_labels = row_labels_order,
                    row_title_gp = gpar(fontsize = 8, fontface = "bold"),
                    row_names_gp = gpar(fontsize = 6),
                    column_names_gp = gpar(fontsize = 6),
                    row_order = order(row_labels_order),
                    rect_gp = gpar(col = "white", lwd = 1),
                    show_column_names = TRUE,
                    show_heatmap_legend = TRUE)
hm_order

# Plot heatmap at the family level
Z.Score.fam_hab_name <- as.data.frame(Z.Score.fam_hab)
str(Z.Score.fam_hab_name)
Z.Score.fam_hab_name <- rownames_to_column(Z.Score.fam_hab_name, var = "OTU")
Z.Score.fam_hab_name <- left_join(Z.Score.fam_hab_name, taxa_info)
head(Z.Score.fam_hab_name)
Z.Score.fam_hab_name <- Z.Score.fam_hab_name[, -(2:51)]

rank4_otu_table_total <- as.data.frame(rank4_otu_table)
rank4_otu_table_total$Total <- rowSums(rank4_otu_table[, -1])
head(rank4_otu_table_total)
rank4_otu_table_total <- rank4_otu_table_total[, -(2:51)]

Z.Score.fam_hab_count_total <- left_join(Z.Score.fam_hab_name, rank4_otu_table_total)
head(Z.Score.fam_hab_count_total)


ha_right_fam = rowAnnotation(
  Abundance = Z.Score.fam_hab_count_total$Total, border = FALSE, col = list(Abundance = abundance_col_fun))
row_labels_fam = Z.Score.fam_hab_count_total$Family

hm_fam <- Heatmap(Z.Score.fam_hab, name = "Z-score, CLR", col = col_matrix,
                  column_title = "16S rRNA soil microbiome", 
                  column_title_gp = gpar(fontface = "bold", fontsize = 14),
                  column_split = as.vector(as.vector(sample_soil_factors$Soil.Type)),
                  border = TRUE,
                  top_annotation = ha_top,
                  right_annotation = ha_right_fam,
                  row_title = "Family",
                  row_labels = row_labels_fam,
                  row_title_gp = gpar(fontsize = 8, fontface = "bold"),
                  row_names_gp = gpar(fontsize = 6),
                  column_names_gp = gpar(fontsize = 6),
                  row_order = order(row_labels_fam),
                  rect_gp = gpar(col = "white", lwd = 1),
                  show_column_names = TRUE,
                  show_heatmap_legend = TRUE)
hm_fam


# Plot heatmap at the genus level
Z.Score.genus_hab_name <- as.data.frame(Z.Score.genus_hab)
str(Z.Score.genus_hab_name)
Z.Score.genus_hab_name <- rownames_to_column(Z.Score.genus_hab_name, var = "OTU")
Z.Score.genus_hab_name <- left_join(Z.Score.genus_hab_name, taxa_info)
head(Z.Score.genus_hab_name)
Z.Score.genus_hab_name <- Z.Score.genus_hab_name[, -(2:51)]

rank5_otu_table_total <- as.data.frame(rank5_otu_table)
rank5_otu_table_total$Total <- rowSums(rank5_otu_table[, -1])
head(rank5_otu_table_total)
rank5_otu_table_total <- rank5_otu_table_total[, -(2:51)]

Z.Score.genus_hab_count_total <- left_join(Z.Score.genus_hab_name, rank5_otu_table_total)
head(Z.Score.genus_hab_count_total)

ha_right_genus = rowAnnotation(
  Abundance = Z.Score.genus_hab_count_total$Total, border = FALSE, col = list(Abundance = abundance_col_fun))
row_labels_genus = Z.Score.genus_hab_count_total$Genus

hm_genus <- Heatmap(Z.Score.genus_hab, name = "Z-score, CLR", col = col_matrix,
                    column_title = "16S rRNA soil microbiome", 
                    column_title_gp = gpar(fontface = "bold", fontsize = 14),
                    column_split = as.vector(as.vector(sample_soil_factors$Soil.Type)),
                    border = TRUE,
                    top_annotation = ha_top,
                    right_annotation = ha_right_genus,
                    row_title = "Genus",
                    row_labels = row_labels_genus,
                    row_title_gp = gpar(fontsize = 8, fontface = "bold"),
                    row_names_gp = gpar(fontsize = 6),
                    column_names_gp = gpar(fontsize = 6),
                    row_order = order(row_labels_genus),
                    rect_gp = gpar(col = "white", lwd = 1),
                    show_column_names = TRUE,
                    show_heatmap_legend = TRUE)
hm_genus

#################################  Venn diagrams ################################## 
# Create a dataframe with an ASV variable using original tax_tab
phy_16_prune_soil_noncontam_prev05_true_filt_scale_tab_asv <- as.data.frame(tax_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale)) %>%
  rownames_to_column(var ="ASV")

# Generate Venn diagram using filtered datasets
# Between soil type
v_soil <- ps_venn(phy_16_prune_soil_noncontam_prev05_true_filt_scale, group = "Soil.Type", type = "counts", plot = TRUE) # creates a plot
v_soil_list <- ps_venn(phy_16_prune_soil_noncontam_prev05_true_filt_scale, group = "Soil.Type", type = "counts", plot = FALSE) # creates a list
class(v_soil_list)
v_soil_list_df <- ldply(v_soil_list, data.frame) # it creates a dataframe with the lists
head(v_soil_list_df)
names(v_soil_list_df) <- c("Soil.Type", "ASV") # rename column names
v_soil_list_df_tax <- left_join(v_soil_list_df, phy_16_prune_soil_noncontam_prev05_true_filt_scale_tab_asv)
write.csv(v_soil_list_df_tax, "exported_tables/soil/v_soil_list_df_tax.csv")

# Among habitats uncompacted soil type
phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_uncompact <- subset_samples(phy_16_prune_soil_noncontam_prev05_true_filt_scale, Soil.Type == "Not.Compact.Soil") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_uncompact # 4556 taxa and 23 samples

v_soil_uncompact <- ps_venn(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_uncompact, group = "Habitat", type = "counts", plot = TRUE) # creates a plot
v_soil_uncompact_list <- ps_venn(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_uncompact, group = "Habitat", type = "counts", plot = FALSE) # creates a list
class(v_soil_uncompact_list)
v_soil_uncompact_list_df <- ldply(v_soil_uncompact_list, data.frame) # it creates a dataframe with the lists
head(v_soil_uncompact_list_df)
names(v_soil_uncompact_list_df) <- c("Habitat", "ASV") # rename column names
v_soil_uncompact_list_df_tax <- left_join(v_soil_uncompact_list_df, phy_16_prune_soil_noncontam_prev05_true_filt_scale_tab_asv)
write.csv(v_soil_uncompact_list_df_tax, "exported_tables/soil/v_soil_uncompact_list_df_tax.csv")

# Among habitats compacted soil type
phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_compact <- subset_samples(phy_16_prune_soil_noncontam_prev05_true_filt_scale, Soil.Type == "Compact.Soil") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_compact # 6006 taxa and 27 samples

v_soil_compact <- ps_venn(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_compact, group = "Habitat", type = "counts", plot = TRUE) # creates a plot
v_soil_compact_list <- ps_venn(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_compact, group = "Habitat", type = "counts", plot = FALSE) # creates a list
class(v_soil_compact_list)
v_soil_compact_list_df <- ldply(v_soil_compact_list, data.frame) # it creates a dataframe with the lists
head(v_soil_compact_list_df)
names(v_soil_compact_list_df) <- c("Habitat", "ASV") # rename column names
v_soil_compact_list_df_tax <- left_join(v_soil_compact_list_df, phy_16_prune_soil_noncontam_prev05_true_filt_scale_tab_asv)
write.csv(v_soil_compact_list_df_tax, "exported_tables/soil/v_soil_compact_list_df_tax.csv")
