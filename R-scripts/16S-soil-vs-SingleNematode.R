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
library("dplyr") #A tool for working/manipulating dataframes. It conflicts with plyr
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

# Compare soil vs. worm samples using clean datasets
shsk_16S_asv_soil_clean <- as.data.frame(read.csv("exported_tables/ASV_SHSK_phy_16_prune_soil_noncontam_prev05_true_all_filt.csv"))
shsk_16S_asv_nema_clean <- as.data.frame(read.csv("exported_tables/ASV_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt.csv"))

shsk_16S_asv_soil_nema_clean_merged <- merge(shsk_16S_asv_soil_clean, shsk_16S_asv_nema_clean,
                                             by = "ASV", all = TRUE,
                                             sort = TRUE, no.dups = TRUE)
shsk_16S_asv_soil_nema_clean_merged[is.na(shsk_16S_asv_soil_nema_clean_merged)] <- 0

write.csv(shsk_16S_asv_soil_nema_clean_merged, "exported_tables/shsk_16S_asv_soil_nema_clean_merged.csv")

shsk_16S_soil_sample_clean <- as.data.frame(read.csv("exported_tables/TABLE_SHSK_phy_16_prune_soil_noncontam_prev05_true_all_filt_df.csv"))
shsk_16S_nema_sample_clean <- as.data.frame(read.csv("exported_tables/TABLE_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_df.csv"))

shsk_16S_soil_nema_sample_clean_merged <- bind_rows(shsk_16S_soil_sample_clean, shsk_16S_nema_sample_clean)
write.csv(shsk_16S_soil_nema_sample_clean_merged, "exported_tables/shsk_16S_soil_nema_sample_clean_merged.csv")

shsk_16S_soil_tax_clean <- as.data.frame(read.csv("exported_tables/TAX_SHSK_phy_16_prune_soil_noncontam_prev05_true_all_filt.csv"))
shsk_16S_nema_tax_clean <- as.data.frame(read.csv("exported_tables/TAX_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt.csv"))

shsk_16S_soil_nema_tax_clean_merged <- full_join(shsk_16S_soil_tax_clean, shsk_16S_nema_tax_clean)
write.csv(shsk_16S_soil_nema_tax_clean_merged, "exported_tables/shsk_16S_soil_nema_tax_clean_merged.csv")


# Import asv table, taxonomy table, and sample table in excel format. It can also be imported as csv files.
asv_tab_16S_com <- read_excel("data/shsk_16S_asv_soil_nema_clean_merged.xlsx")
tax_tab_16S_com <- read_excel("data/shsk_16S_soil_nema_tax_clean_merged.xlsx")
sample_tab_16S_com <- read_excel("data/shsk_16S_soil_nema_sample_clean_merged.xlsx")

# Convert asv and taxonomy tables to matrix format; sample table as dataframe
asv_tab_16S_com <- as.matrix(asv_tab_16S_com)
class(asv_tab_16S_com) # check the object class
tax_tab_16S_com <- as.matrix(tax_tab_16S_com)
sample_tab_16S_com <- as.data.frame(sample_tab_16S_com)
class(sample_tab_16S_com) # check the object class

# Rename row names using ASV column for asv table and taxonomy table; rename using Sample column for sample table
rownames(asv_tab_16S_com) <- asv_tab_16S_com[,1]
rownames(tax_tab_16S_com) <- tax_tab_16S_com[,1]
rownames(sample_tab_16S_com) <- sample_tab_16S_com[,1]

# Exclude first column from all three tables; set asv table as numeric
asv_tab_16S_com <- asv_tab_16S_com[,-1]
class(asv_tab_16S_com) <- "numeric"
tax_tab_16S_com <- tax_tab_16S_com[,-1]
sample_tab_16S_com <- sample_tab_16S_com[,-1]

#Transform each table/matrix to phyloseq objects
ASV_16S_com = otu_table(asv_tab_16S_com, taxa_are_rows = TRUE)
TAX_16S_com = tax_table(tax_tab_16S_com)
samples_16S_com = sample_data(sample_tab_16S_com)

# Create a phyloseq object by combine all three previous objects
phy_16S_com <- phyloseq(ASV_16S_com, TAX_16S_com, samples_16S_com)
phy_16S_com
head(otu_table(phy_16S_com)) # check the first 6 lines of your phyloseq object

# Export combined phyloseq object with all true samples filtered mitochondria/chloroplasts
# This can be used as a starting point for further analysis
saveRDS(phy_16S_com, "data/phy_16S_com.RDS")

# Import phyloseq
phy_16S_com <- readRDS("data/phy_16S_com.RDS")
class(phy_16S_com)
smin <- min(sample_sums(phy_16S_com)) # smin = 0 means you have samples with no reads. Total of 530 samples.
unique(sample_data(phy_16S_com)$Experiment)

# Apply different transformations to ASV table, e.g., converting number of reads to relative abundance ()
# Here we use the library microbiome
phy_16S_com_comp <- microbiome::transform(phy_16S_com, "compositional") # transform reads to relative abundance (0.0 - 1.0). Same as # phy_16S_com_comp = transform_sample_counts(phy_16S_com, function(x) 100 * x/sum(x))
phy_16S_com_hel <- microbiome::transform(phy_16S_com, "hellinger")
phy_16S_com_clr <- microbiome::transform(phy_16S_com, "clr") # transform reads using centered log ratio
phy_16S_com_log10 <- microbiome::transform(phy_16S_com, "log10p") # transform log x +1

set.seed(1)
phy_16S_com_comp_nmds <- ordinate(
  physeq = phy_16S_com_comp, 
  method = "NMDS", 
  distance = "bray"
)
############### Ordination using PCA and clr transformed phyloseq objects
phy_16S_com_clr_pca <- phyloseq::ordinate(
  physeq = phy_16S_com_clr,
  method = "RDA",
  distance = "euclidean"
)

# Plot scree plot to see PCA contributions
phyloseq::plot_scree(phy_16S_com_clr_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/all-dataset/16S_Soil_ShSk_ASV_noncontam_prev05_true_all_filt_clr_pca.pdf", width = 12, height = 6, dpi = 150)

head(phy_16S_com_clr_pca$CA$eig)
sapply(phy_16S_com_clr_pca$CA$eig[1:5], function(x) x / sum(phy_16S_com_clr_pca$CA$eig))

# Scale axes and plot ordination
phy_16S_com_clr_pca_clr1 <- phy_16S_com_clr_pca$CA$eig[1] / sum(phy_16S_com_clr_pca$CA$eig)
phy_16S_com_clr_pca_clr2 <- phy_16S_com_clr_pca$CA$eig[2] / sum(phy_16S_com_clr_pca$CA$eig)

# Plot PCA based on clr transformation
phyloseq::plot_ordination(
  physeq = phy_16S_com_clr,
  ordination = phy_16S_com_clr_pca,
  type="samples",
  color="Sample.Type",
  shape = "Soil.Type") + 
  scale_color_manual(values = c("#52392F", "#397A4C"),
                     name = "Sample Type",
                     breaks = c("SingleWorm", "Soil"),
                     labels = c("Nematode", "Soil"))+
  scale_shape_manual(values = c(19,17),
                     name = "Soil Type",
                     breaks = c("Compact.Soil", "Not.Compact.Soil"),
                     labels = c("Compacted", "Uncompacted"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Sample.Type), alpha = 0.7, size = 4) +
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
  theme(strip.background =element_blank())# remove the background of titles
#annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
#annotate("text", x = 3, y = 2, label ="2D Stress: 0.07") +
#coord_fixed(phy_16_prune_soil_noncontam_prev05_true_all_filt_clr_pca_clr2 / phy_16_prune_soil_noncontam_prev05_true_all_filt_clr_pca_clr1) +
ggsave("results/all-dataset/phy_16S_com_clr_pca_plot.pdf", width = 6, height = 4, dpi = 200)

###################################### Barplots ###################################### 
merge_less_than_top_prev05_top20 <- function(phy_16S_com, top=19){
  transformed <- transform_sample_counts(phy_16S_com, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:6] <- "Others"} # 1:7 if there are species level
  }
  return(merged)
}

merge_less_than_top_prev05_top12 <- function(phy_16S_com, top=11){
  transformed <- microbiome::transform(phy_16S_com, "compositional")
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:6] <- "Others"} # 1:7 if there are species level
  }
  return(merged)
}

phy_16S_com_top_phylum <- tax_glom(phy_16S_com, "Phylum")
phy_16S_com_top_phylum_top12 <- merge_less_than_top_prev05_top12(phy_16S_com_top_phylum, top=11)
phy_16S_com_top_phylum_top12_df = psmelt(phy_16S_com_top_phylum_top12) # Melt to long format, it will give a warning messsage because of table headers
phy_16S_com_top_phylum_top12_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Experiment+Phylum, data=phy_16S_com_top_phylum_top12_df, FUN=mean) # add common factors to use for plotting
phy_16S_com_top_phylum_top12_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Experiment+Phylum, data=phy_16S_com_top_phylum_top12_df, FUN=mean) # add common factors to use for plotting
phy_16S_com_top_phylum_top12_agr_soil = aggregate(Abundance~Soil.Type+Experiment+Phylum, data=phy_16S_com_top_phylum_top12_df, FUN=mean) # add common factors to use for plotting
unique(phy_16S_com_top_phylum_top12_agr_site$Phylum)

phy_16S_com_top_class <- tax_glom(phy_16S_com, "Class")
phy_16S_com_top_class_top12 <- merge_less_than_top_prev05_top12(phy_16S_com_top_class, top=11)
phy_16S_com_top_class_top12_df = psmelt(phy_16S_com_top_class_top12) # Melt to long format, it will give a warning messsage because of table headers
phy_16S_com_top_class_top12_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Experiment+Class, data=phy_16S_com_top_class_top12_df, FUN=mean) # add common factors to use for plotting
phy_16S_com_top_class_top12_agr_hab = aggregate(Abundance~Soil.Type+Experiment+Habitat+Class, data=phy_16S_com_top_class_top12_df, FUN=mean) # add common factors to use for plotting
phy_16S_com_top_class_top12_agr_soil = aggregate(Abundance~Soil.Type+Experiment+Class, data=phy_16S_com_top_class_top12_df, FUN=mean) # add common factors to use for plotting
unique(phy_16S_com_top_class_top12_agr_site$Class)

phy_16S_com_top_order <- tax_glom(phy_16S_com, "Order")
phy_16S_com_top_order_top12 <- merge_less_than_top_prev05_top12(phy_16S_com_top_order, top=11)
phy_16S_com_top_order_top12_df = psmelt(phy_16S_com_top_order_top12) # Melt to long format, it will give a warning messsage because of table headers
phy_16S_com_top_order_top12_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Experiment+Order, data=phy_16S_com_top_order_top12_df, FUN=mean) # add common factors to use for plotting
phy_16S_com_top_order_top12_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Experiment+Order, data=phy_16S_com_top_order_top12_df, FUN=mean) # add common factors to use for plotting
phy_16S_com_top_order_top12_agr_soil = aggregate(Abundance~Soil.Type+Experiment+Order, data=phy_16S_com_top_order_top12_df, FUN=mean) # add common factors to use for plotting
unique(phy_16S_com_top_order_top12_agr_site$Order)

phy_16S_com_top_order <- tax_glom(phy_16S_com, "Order")
phy_16S_com_top_order_top20 <- merge_less_than_top_prev05_top20(phy_16S_com_top_order, top=19)
phy_16S_com_top_order_top20_df = psmelt(phy_16S_com_top_order_top20) # Melt to long format, it will give a warning messsage because of table headers
phy_16S_com_top_order_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Experiment+Order, data=phy_16S_com_top_order_top20_df, FUN=mean) # add common factors to use for plotting
phy_16S_com_top_order_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Experiment+Order, data=phy_16S_com_top_order_top20_df, FUN=mean) # add common factors to use for plotting
phy_16S_com_top_order_top20_agr_soil = aggregate(Abundance~Soil.Type+Experiment+Order, data=phy_16S_com_top_order_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16S_com_top_order_top20_agr_site$Order)

# List of 12 distinct colors
colors_top12 <- c("#9a6324", "#46f0f0", "#f032e6","#aaffc3", "#33A02C", "#4363d8", 
                  "#bcf60c", "#f58231", "#6A3D9A", "#e6beff", "#fffac8", "gray")

# List of 20 distinct colors
colors_top20 <- c("#9a6324", "#46f0f0", "#1F78B4", "#aaffc3", "#000075",  "#33A02C", "#4363d8", "#008080", "#FB9A99", "#e6194b",
                  "#f032e6", "#bcf60c", "#fabebe","#f58231", "#ffe119", "#6A3D9A", "#e6beff", "#800000", "#fffac8", "gray")


# Create new labels for important factors
soil.labs <- c("Compacted", "Uncompacted")
names(soil.labs) <- c("Compact.Soil", "Not.Compact.Soil")

habitat.labs <- c("CHA", "CSS", "NGR", "HLC",  "OWL", "RIP")
names(habitat.labs) <- c("Chaparral", "Coastal.scrub", "Native.grass",
                         "Hollyleaf.cherry", "Oak.wood",  "Riparian")

exp.labs <- c("Soil", "Nematodes")
names(exp.labs) <- c("Environment", "Microbiome")

# Put "Others" to the final of the Phylum list - top 12
unique(phy_16S_com_top_phylum_top12_agr_site$Phylum)
phy_16S_com_top_phylum_top12_agr_site$Phylum <- factor(phy_16S_com_top_phylum_top12_agr_site$Phylum,
                                                                                              levels = c("Acidobacteriota", "Actinobacteriota", "Bacteroidota",
                                                                                                         "Chloroflexi", "Cyanobacteria", "Deinococcota",
                                                                                                         "Firmicutes", "Gemmatimonadota", "Planctomycetota",
                                                                                                         "Proteobacteria","Verrucomicrobiota", "Others"))
# Plot by site - Phylum level top 12
ggplot(phy_16S_com_top_phylum_top12_agr_site, aes(x = Sample.Site, y = Abundance, fill = Phylum)) +
  facet_nested(Experiment ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs,
                                   Experiment = exp.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top12) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    breaks = c("SK.01", "SK.02", "SK.03",
               "SK.04", "SK.05", "SK.06",
               "SK.07", "SK.08", "SK.09",
               "SK.10", "SK.11", "SK.12",
               "SK.13", "SK.14", "SK.15",
               "SK.16", "SK.17", "SK.18"),
    labels = c("01", "02", "03",
               "04", "05", "06",
               "07", "08", "09",
               "10", "11", "12",
               "13", "14", "15",
               "16", "17", "18"),
    expand = c(0.0, 0.0),
    drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/all-dataset/phy_16S_com_top_phylum_top12_agr_site_phylum.pdf", width = 9, height = 6, dpi = 200) # save graphic


# Put "Others" to the final of the Class list - top 12
unique(phy_16S_com_top_class_top12_agr_site$Class)
phy_16S_com_top_class_top12_agr_site$Class <- factor(phy_16S_com_top_class_top12_agr_site$Class,
                                                     levels = c("Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidia",
                                                                "Deinococci", "Gammaproteobacteria", "Planctomycetes", "Thermoleophilia",
                                                                "Vampirivibrionia", "Verrucomicrobiae", "Vicinamibacteria", "Others"))
                                                       
# Plot by site - Class level top 12
ggplot(phy_16S_com_top_class_top12_agr_site, aes(x = Sample.Site, y = Abundance, fill = Class)) +
  facet_nested(Experiment ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs,
                                   Experiment = exp.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top12) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    breaks = c("SK.01", "SK.02", "SK.03",
               "SK.04", "SK.05", "SK.06",
               "SK.07", "SK.08", "SK.09",
               "SK.10", "SK.11", "SK.12",
               "SK.13", "SK.14", "SK.15",
               "SK.16", "SK.17", "SK.18"),
    labels = c("01", "02", "03",
               "04", "05", "06",
               "07", "08", "09",
               "10", "11", "12",
               "13", "14", "15",
               "16", "17", "18"),
    expand = c(0.0, 0.0),
    drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/all-dataset/phy_16S_com_top_phylum_top12_agr_site_class.pdf", width = 9, height = 6, dpi = 200) # save graphic

# Put "Others" to the final of the Order list - top 12
unique(phy_16S_com_top_order_top12_agr_site$Order)
phy_16S_com_top_order_top12_agr_site$Order <- factor(phy_16S_com_top_order_top12_agr_site$Order,
                                                     levels = c("Bacillales", "Burkholderiales", "Corynebacteriales", "Cytophagales",
                                                                "Flavobacteriales", "Lactobacillales", "Opitutales", "Pseudomonadales",
                                                                "Rhizobiales", "Sphingomonadales", "Staphylococcales", "Others"))

# Plot by site - order level top 12
ggplot(phy_16S_com_top_order_top12_agr_site, aes(x = Sample.Site, y = Abundance, fill = Order)) +
  facet_nested(Experiment ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs,
                                   Experiment = exp.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top12) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    breaks = c("SK.01", "SK.02", "SK.03",
               "SK.04", "SK.05", "SK.06",
               "SK.07", "SK.08", "SK.09",
               "SK.10", "SK.11", "SK.12",
               "SK.13", "SK.14", "SK.15",
               "SK.16", "SK.17", "SK.18"),
    labels = c("01", "02", "03",
               "04", "05", "06",
               "07", "08", "09",
               "10", "11", "12",
               "13", "14", "15",
               "16", "17", "18"),
    expand = c(0.0, 0.0),
    drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/all-dataset/phy_16S_com_top_phylum_top12_agr_site_order.pdf", width = 9, height = 6, dpi = 200) # save graphic


# Put "Others" to the final of the Order list - top 20
unique(phy_16S_com_top_order_top20_agr_site$Order)
phy_16S_com_top_order_top20_agr_site$Order <- factor(phy_16S_com_top_order_top20_agr_site$Order,
                                                     levels = c("Bacillales", "Burkholderiales", "Caulobacterales", "Corynebacteriales",
                                                                "Cytophagales", "Enterobacterales", "Flavobacteriales", "Lactobacillales",
                                                                "Micrococcales", "Opitutales", "Pseudomonadales", "Rhizobiales",
                                                                "Rhodobacterales", "Rickettsiales", "Solirubrobacterales", "Sphingomonadales",
                                                                "Staphylococcales", "Vicinamibacterales", "Xanthomonadales", "Others"))

# Plot by site - Order level top 20
ggplot(phy_16S_com_top_order_top20_agr_site, aes(x = Sample.Site, y = Abundance, fill = Order)) +
  facet_nested(Experiment ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs,
                                   Experiment = exp.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    breaks = c("SK.01", "SK.02", "SK.03",
               "SK.04", "SK.05", "SK.06",
               "SK.07", "SK.08", "SK.09",
               "SK.10", "SK.11", "SK.12",
               "SK.13", "SK.14", "SK.15",
               "SK.16", "SK.17", "SK.18"),
    labels = c("01", "02", "03",
               "04", "05", "06",
               "07", "08", "09",
               "10", "11", "12",
               "13", "14", "15",
               "16", "17", "18"),
    expand = c(0.0, 0.0),
    drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/all-dataset/phy_16S_com_top_phylum_top20_agr_site_order.pdf", width = 9, height = 6, dpi = 200) # save graphic

head(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_site)
sum(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_site$Abundance)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_site,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_site.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_hab,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_hab.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_soil,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_soil.csv")

###################################### Subset dominant genera ###################################### 
phy_16S_com_mycobac <- phy_16S_com %>%
  subset_taxa(
    Genus  == "Mycobacterium"
  )
get_taxa_unique(phy_16S_com_mycobac, taxonomic.rank = "Genus")
phy_16S_com_mycobac # 85 taxa and 582 samples

phy_16S_com_staphyl <- phy_16S_com %>%
  subset_taxa(
    Genus  == "Staphylococcus"
  )
get_taxa_unique(phy_16S_com_staphyl, taxonomic.rank = "Genus")
phy_16S_com_staphyl # 67 taxa and 582 samples

phy_16S_com_novos <- phy_16S_com %>%
  subset_taxa(
    Genus  == "Novosphingobium"
  )
get_taxa_unique(phy_16S_com_novos, taxonomic.rank = "Genus")
phy_16S_com_novos # 74 taxa and 582 samples

phy_16S_com_lacun <- phy_16S_com %>%
  subset_taxa(
    Genus  == "Lacunisphaera"
  )
get_taxa_unique(phy_16S_com_lacun, taxonomic.rank = "Genus")
phy_16S_com_lacun # 42 taxa and 582 samples

phy_16S_com_ralst <- phy_16S_com %>%
  subset_taxa(
    Genus  == "Ralstonia"
  )
get_taxa_unique(phy_16S_com_ralst, taxonomic.rank = "Genus")
phy_16S_com_ralst # 9 taxa and 582 samples

phy_16S_com_tetrag <- phy_16S_com %>%
  subset_taxa(
    Genus  == "Tetragenococcus"
  )
get_taxa_unique(phy_16S_com_tetrag, taxonomic.rank = "Genus")
phy_16S_com_tetrag # 13 taxa and 582 samples

phy_16S_com_cardi <- phy_16S_com %>%
  subset_taxa(
    Genus  == "Candidatus Cardinium"
  )
get_taxa_unique(phy_16S_com_cardi, taxonomic.rank = "Genus")
phy_16S_com_cardi # 18 taxa and 582 samples

# Aldex on mycobacterium ASVs
phy_16S_com_mycobac_exp <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16S_com_mycobac)),
                                         phyloseq::sample_data(phy_16S_com_mycobac)$Experiment,
                                         mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
phy_16S_com_mycobac_exp_otu_table <- data.frame(phyloseq::otu_table(phy_16S_com_mycobac))
phy_16S_com_mycobac_exp_otu_table <- rownames_to_column(phy_16S_com_mycobac_exp_otu_table, var = "OTU")
write.csv(phy_16S_com_mycobac_exp_otu_table, "exported_tables/all-dataset/phy_16S_com_mycobac_exp_otu_table.csv")


phy_16S_com_staphyl_exp <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16S_com_staphyl)),
                                  phyloseq::sample_data(phy_16S_com_staphyl)$Experiment,
                                  mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
phy_16S_com_staphyl_exp_otu_table <- data.frame(phyloseq::otu_table(phy_16S_com_staphyl))
phy_16S_com_staphyl_exp_otu_table <- rownames_to_column(phy_16S_com_staphyl_exp_otu_table, var = "OTU")
write.csv(phy_16S_com_staphyl_exp_otu_table, "exported_tables/all-dataset/phy_16S_com_staphyl_exp_otu_table.csv")

phy_16S_com_novos_exp <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16S_com_novos)),
                                         phyloseq::sample_data(phy_16S_com_novos)$Experiment,
                                         mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
phy_16S_com_novos_exp_otu_table <- data.frame(phyloseq::otu_table(phy_16S_com_novos))
phy_16S_com_novos_exp_otu_table <- rownames_to_column(phy_16S_com_novos_exp_otu_table, var = "OTU")
write.csv(phy_16S_com_novos_exp_otu_table, "exported_tables/all-dataset/phy_16S_com_novos_exp_otu_table.csv")

phy_16S_com_lacun_exp <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16S_com_lacun)),
                                         phyloseq::sample_data(phy_16S_com_lacun)$Experiment,
                                         mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
phy_16S_com_lacun_exp_otu_table <- data.frame(phyloseq::otu_table(phy_16S_com_lacun))
phy_16S_com_lacun_exp_otu_table <- rownames_to_column(phy_16S_com_lacun_exp_otu_table, var = "OTU")
write.csv(phy_16S_com_lacun_exp_otu_table, "exported_tables/all-dataset/phy_16S_com_lacun_exp_otu_table.csv")

phy_16S_com_ralst_exp <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16S_com_ralst)),
                                       phyloseq::sample_data(phy_16S_com_ralst)$Experiment,
                                       mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
phy_16S_com_ralst_exp_otu_table <- data.frame(phyloseq::otu_table(phy_16S_com_ralst))
phy_16S_com_ralst_exp_otu_table <- rownames_to_column(phy_16S_com_ralst_exp_otu_table, var = "OTU")
write.csv(phy_16S_com_ralst_exp_otu_table, "exported_tables/all-dataset/phy_16S_com_ralst_exp_otu_table.csv")

phy_16S_com_tetrag_exp <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16S_com_tetrag)),
                                       phyloseq::sample_data(phy_16S_com_tetrag)$Experiment,
                                       mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
phy_16S_com_tetrag_exp_otu_table <- data.frame(phyloseq::otu_table(phy_16S_com_tetrag))
phy_16S_com_tetrag_exp_otu_table <- rownames_to_column(phy_16S_com_tetrag_exp_otu_table, var = "OTU")
write.csv(phy_16S_com_tetrag_exp_otu_table, "exported_tables/all-dataset/phy_16S_com_tetrag_exp_otu_table.csv")

phy_16S_com_cardi_exp <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_16S_com_cardi)),
                                       phyloseq::sample_data(phy_16S_com_cardi)$Experiment,
                                       mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
phy_16S_com_cardi_exp_otu_table <- data.frame(phyloseq::otu_table(phy_16S_com_cardi))
phy_16S_com_cardi_exp_otu_table <- rownames_to_column(phy_16S_com_cardi_exp_otu_table, var = "OTU")
write.csv(phy_16S_com_cardi_exp_otu_table, "exported_tables/all-dataset/phy_16S_com_cardi_exp_otu_table.csv")

# Write aldex2 results to file as csv
# Experiment comparison results
write.csv(phy_16S_com_mycobac_exp, "exported_tables/all-dataset/phy_16S_com_mycobac_exp.csv")
write.csv(phy_16S_com_staphyl_exp, "exported_tables/all-dataset/phy_16S_com_staphyl_exp.csv")
write.csv(phy_16S_com_novos_exp, "exported_tables/all-dataset/phy_16S_com_novos_exp.csv")
write.csv(phy_16S_com_lacun_exp, "exported_tables/all-dataset/phy_16S_com_lacun_exp.csv")
write.csv(phy_16S_com_ralst_exp, "exported_tables/all-dataset/phy_16S_com_ralst_exp.csv")
write.csv(phy_16S_com_tetrag_exp, "exported_tables/all-dataset/phy_16S_com_tetrag_exp.csv")
write.csv(phy_16S_com_cardi_exp, "exported_tables/all-dataset/phy_16S_com_cardi_exp.csv")

# Import aldex2 results and rename the X variable by OTU
phy_16S_com_mycobac_exp_result <- read.csv("exported_tables/all-dataset/phy_16S_com_mycobac_exp.csv")
colnames(phy_16S_com_mycobac_exp_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

phy_16S_com_staphyl_exp_result <- read.csv("exported_tables/all-dataset/phy_16S_com_staphyl_exp.csv")
colnames(phy_16S_com_staphyl_exp_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

phy_16S_com_novos_exp_result <- read.csv("exported_tables/all-dataset/phy_16S_com_novos_exp.csv")
colnames(phy_16S_com_novos_exp_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

phy_16S_com_lacun_exp_result <- read.csv("exported_tables/all-dataset/phy_16S_com_lacun_exp.csv")
colnames(phy_16S_com_lacun_exp_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

phy_16S_com_ralst_exp_result <- read.csv("exported_tables/all-dataset/phy_16S_com_ralst_exp.csv")
colnames(phy_16S_com_ralst_exp_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

phy_16S_com_tetrag_exp_result <- read.csv("exported_tables/all-dataset/phy_16S_com_tetrag_exp.csv")
colnames(phy_16S_com_tetrag_exp_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

phy_16S_com_cardi_exp_result <- read.csv("exported_tables/all-dataset/phy_16S_com_cardi_exp.csv")
colnames(phy_16S_com_cardi_exp_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

#Clean up presentation
taxa_info_mycobac <- data.frame(tax_table(phy_16S_com_mycobac))
taxa_info_mycobac <- taxa_info_mycobac %>%
  rownames_to_column(var = "OTU")
sample_tab_mycobac <- data.frame(sample_data(phy_16S_com_mycobac))
sample_tab_mycobac_factors <- data.frame(cbind(sample_tab_mycobac[, c(19, 21, 22, 36, 37)]))

taxa_info_staphyl <- data.frame(tax_table(phy_16S_com_staphyl))
taxa_info_staphyl <- taxa_info_staphyl %>%
  rownames_to_column(var = "OTU")
sample_tab_staphyl <- data.frame(sample_data(phy_16S_com_staphyl))
sample_tab_staphyl_factors <- data.frame(cbind(sample_tab_staphyl[, c(19, 21, 22, 36, 37)]))

taxa_info_novos <- data.frame(tax_table(phy_16S_com_novos))
taxa_info_novos <- taxa_info_novos %>%
  rownames_to_column(var = "OTU")
sample_tab_novos <- data.frame(sample_data(phy_16S_com_novos))
sample_tab_novos_factors <- data.frame(cbind(sample_tab_novos[, c(19, 21, 22, 36, 37)]))

taxa_info_lacun <- data.frame(tax_table(phy_16S_com_lacun))
taxa_info_lacun <- taxa_info_lacun %>%
  rownames_to_column(var = "OTU")
sample_tab_lacun <- data.frame(sample_data(phy_16S_com_lacun))
sample_tab_lacun_factors <- data.frame(cbind(sample_tab_lacun[, c(19, 21, 22, 36, 37)]))

taxa_info_ralst <- data.frame(tax_table(phy_16S_com_ralst))
taxa_info_ralst <- taxa_info_ralst %>%
  rownames_to_column(var = "OTU")
sample_tab_ralst <- data.frame(sample_data(phy_16S_com_ralst))
sample_tab_ralst_factors <- data.frame(cbind(sample_tab_ralst[, c(19, 21, 22, 36, 37)]))

taxa_info_tetrag <- data.frame(tax_table(phy_16S_com_tetrag))
taxa_info_tetrag <- taxa_info_tetrag %>%
  rownames_to_column(var = "OTU")
sample_tab_tetrag <- data.frame(sample_data(phy_16S_com_tetrag))
sample_tab_tetrag_factors <- data.frame(cbind(sample_tab_tetrag[, c(19, 21, 22, 36, 37)]))

taxa_info_cardi <- data.frame(tax_table(phy_16S_com_cardi))
taxa_info_cardi <- taxa_info_cardi %>%
  rownames_to_column(var = "OTU")
sample_tab_cardi <- data.frame(sample_data(phy_16S_com_cardi))
sample_tab_cardi_factors <- data.frame(cbind(sample_tab_cardi[, c(19, 21, 22, 36, 37)]))

# filter aldex2 results by sig kw.ep and join the taxanomic information
sig_phy_16S_com_mycobac_exp_result <- phy_16S_com_mycobac_exp_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_phy_16S_com_mycobac_exp_result <- left_join(sig_phy_16S_com_mycobac_exp_result, taxa_info_mycobac)
write.csv(sig_phy_16S_com_mycobac_exp_result, "exported_tables/all-dataset/sig_phy_16S_com_mycobac_exp_result.csv") # 5 taxa with significant differential abundance

sig_phy_16S_com_staphyl_exp_result <- phy_16S_com_staphyl_exp_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_phy_16S_com_staphyl_exp_result <- left_join(sig_phy_16S_com_staphyl_exp_result, taxa_info_staphyl)
write.csv(sig_phy_16S_com_staphyl_exp_result, "exported_tables/all-dataset/sig_phy_16S_com_staphyl_exp_result.csv") # 5 taxa with significant differential abundance

sig_phy_16S_com_novos_exp_result <- phy_16S_com_novos_exp_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_phy_16S_com_novos_exp_result <- left_join(sig_phy_16S_com_novos_exp_result, taxa_info_novos)
write.csv(sig_phy_16S_com_novos_exp_result, "exported_tables/all-dataset/sig_phy_16S_com_novos_exp_result.csv") # 5 taxa with significant differential abundance

sig_phy_16S_com_lacun_exp_result <- phy_16S_com_lacun_exp_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_phy_16S_com_lacun_exp_result <- left_join(sig_phy_16S_com_lacun_exp_result, taxa_info_lacun)
write.csv(sig_phy_16S_com_lacun_exp_result, "exported_tables/all-dataset/sig_phy_16S_com_lacun_exp_result.csv") # 5 taxa with significant differential abundance

sig_phy_16S_com_ralst_exp_result <- phy_16S_com_ralst_exp_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_phy_16S_com_ralst_exp_result <- left_join(sig_phy_16S_com_ralst_exp_result, taxa_info_ralst)
write.csv(sig_phy_16S_com_ralst_exp_result, "exported_tables/all-dataset/sig_phy_16S_com_ralst_exp_result.csv") # 5 taxa with significant differential abundance

sig_phy_16S_com_tetrag_exp_result <- phy_16S_com_tetrag_exp_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_phy_16S_com_tetrag_exp_result <- left_join(sig_phy_16S_com_tetrag_exp_result, taxa_info_tetrag)
write.csv(sig_phy_16S_com_tetrag_exp_result, "exported_tables/all-dataset/sig_phy_16S_com_tetrag_exp_result.csv") # 5 taxa with significant differential abundance

sig_phy_16S_com_cardi_exp_result <- phy_16S_com_cardi_exp_result %>%
  filter(kw.ep < 0.3) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_phy_16S_com_cardi_exp_result <- left_join(sig_phy_16S_com_cardi_exp_result, taxa_info_cardi)
write.csv(sig_phy_16S_com_cardi_exp_result, "exported_tables/all-dataset/sig_phy_16S_com_cardi_exp_result.csv") # 5 taxa with significant differential abundance

# Create clr objects by using OTU ids and original otu table, all significant OTUs and only top 30
sig_phy_16S_com_mycobac_exp_result_count <- left_join(sig_phy_16S_com_mycobac_exp_result, phy_16S_com_mycobac_exp_otu_table)
write.csv(sig_phy_16S_com_mycobac_exp_result_count, "exported_tables/all-dataset/sig_phy_16S_com_mycobac_exp_result_count.csv")

sig_phy_16S_com_staphyl_exp_result_count <- left_join(sig_phy_16S_com_staphyl_exp_result, phy_16S_com_staphyl_exp_otu_table)
write.csv(sig_phy_16S_com_staphyl_exp_result_count, "exported_tables/all-dataset/sig_phy_16S_com_staphyl_exp_result_count.csv")

sig_phy_16S_com_novos_exp_result_count <- left_join(sig_phy_16S_com_novos_exp_result, phy_16S_com_novos_exp_otu_table)
write.csv(sig_phy_16S_com_novos_exp_result_count, "exported_tables/all-dataset/sig_phy_16S_com_novos_exp_result_count.csv")

sig_phy_16S_com_lacun_exp_result_count <- left_join(sig_phy_16S_com_lacun_exp_result, phy_16S_com_lacun_exp_otu_table)
write.csv(sig_phy_16S_com_lacun_exp_result_count, "exported_tables/all-dataset/sig_phy_16S_com_lacun_exp_result_count.csv")

sig_phy_16S_com_ralst_exp_result_count <- left_join(sig_phy_16S_com_ralst_exp_result, phy_16S_com_ralst_exp_otu_table)
write.csv(sig_phy_16S_com_ralst_exp_result_count, "exported_tables/all-dataset/sig_phy_16S_com_ralst_exp_result_count.csv")

sig_phy_16S_com_tetrag_exp_result_count <- left_join(sig_phy_16S_com_tetrag_exp_result, phy_16S_com_tetrag_exp_otu_table)
write.csv(sig_phy_16S_com_tetrag_exp_result_count, "exported_tables/all-dataset/sig_phy_16S_com_tetrag_exp_result_count.csv")

sig_phy_16S_com_cardi_exp_result_count <- left_join(sig_phy_16S_com_cardi_exp_result, phy_16S_com_cardi_exp_otu_table)
write.csv(sig_phy_16S_com_cardi_exp_result_count, "exported_tables/all-dataset/sig_phy_16S_com_cardi_exp_result_count.csv")
#########################################
# Apply psmelt on phy_16S_com_mycobac object for plot reasons
phy_16S_com_mycobac_df = psmelt(phy_16S_com_mycobac) # Melt to long format, it will give a warning messsage because of table headers
phy_16S_com_mycobac_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Experiment+Genus, data=phy_16S_com_mycobac_df, FUN=mean) # add common factors to use for plotting
phy_16S_com_mycobac_agr_habitat = aggregate(Abundance~Soil.Type+Habitat+Experiment+Genus, data=phy_16S_com_mycobac_df, FUN=mean) # add common factors to use for plotting
phy_16S_com_mycobac_agr_soil = aggregate(Abundance~Soil.Type+Experiment+Genus, data=phy_16S_com_mycobac_df, FUN=mean) # add common factors to use for plotting
unique(phy_16S_com_mycobac_agr_site$Genus)

# Save mycobacterium data
phy_16S_com_mycobac_asv <- as.data.frame(otu_table(phy_16S_com_mycobac, "matrix"))
phy_16S_com_mycobac_tax <- as.data.frame(tax_table(phy_16S_com_mycobac))
phy_16S_com_mycobac_map <- as.data.frame(sample_data(phy_16S_com_mycobac))

write.csv(phy_16S_com_mycobac_asv, "exported_tables/all-dataset/phy_16S_com_mycobac_asv.csv")
write.csv(phy_16S_com_mycobac_tax, "exported_tables/all-dataset/phy_16S_com_mycobac_tax.csv")
write.csv(phy_16S_com_mycobac_map, "exported_tables/all-dataset/phy_16S_com_mycobac_map.csv")

mycobacterium_selec <- read_excel("exported_tables/all-dataset/phy_16S_com_mycobac_asv_selected.xlsx")
class(mycobacterium_selec)

mycobacterium_selec_sum_asv00006 <- summarySE(mycobacterium_selec,
                                              measurevar = c("ASV00006"),
                                              groupvars = c("Soil.Type", "Habitat","Sample.Site", "Experiment"))
write.csv(mycobacterium_selec_sum_asv00006, "exported_tables/all-dataset/mycobacterium_selec_sum_asv00006.csv")

mycobacterium_selec_sum_asv00809 <- summarySE(mycobacterium_selec,
                                              measurevar = c("ASV00809"),
                                              groupvars = c("Soil.Type", "Habitat","Sample.Site", "Experiment"))
write.csv(mycobacterium_selec_sum_asv00809, "exported_tables/all-dataset/mycobacterium_selec_sum_asv00809.csv")

mycobacterium_selec_sum_asv00427 <- summarySE(mycobacterium_selec,
                                              measurevar = c("ASV00427"),
                                              groupvars = c("Soil.Type", "Habitat","Sample.Site", "Experiment"))
write.csv(mycobacterium_selec_sum_asv00427, "exported_tables/all-dataset/mycobacterium_selec_sum_asv00427.csv")

mycobacterium_selec_sum_all <- read.csv("exported_tables/all-dataset/mycobacterium_selec_sum_all.csv")

# Create new labels for important factors
soil.labs <- c("Compacted", "Uncompacted")
names(soil.labs) <- c("Compact.Soil", "Not.Compact.Soil")

habitat.labs <- c("CHA", "CSS", "NGR", "HLC",  "OWL", "RIP")
names(habitat.labs) <- c("Chaparral", "Coastal.scrub", "Native.grass",
                         "Hollyleaf.cherry", "Oak.wood",  "Riparian")

ggplot(mycobacterium_selec_sum_all, aes(x= Sample.Site, y= Read.Count, fill = Experiment)) +
  geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=Read.Count-se, ymax=Read.Count+se), size = 0.5, width = 0.2, position=position_dodge(0.9)) +
  geom_point(size = 2, shape = 21, position=position_dodge(0.9)) +
  scale_fill_manual(values=c("#6A3D9A", "#FF7F00"))+
  facet_nested(ASV ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(
    breaks = c("SK.01", "SK.02", "SK.03",
               "SK.04", "SK.05", "SK.06",
               "SK.07", "SK.08", "SK.09",
               "SK.10", "SK.11", "SK.12",
               "SK.13", "SK.14", "SK.15",
               "SK.16", "SK.17", "SK.18"),
    labels = c("01", "02", "03",
               "04", "05", "06",
               "07", "08", "09",
               "10", "11", "12",
               "13", "14", "15",
               "16", "17", "18"),
    expand = c(0.0, 0.0),
    drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Read Count") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold"))
ggsave("results/all-dataset/mycobacterium_selec_sum_all.pdf", width = 8, height = 6, dpi = 300) # save graphic
#################################  Venn diagrams ################################## 
# Create a dataframe with an ASV variable using original tax_tab
tax_tab_asv <- as.data.frame(tax_table(phy_16S_com)) %>%
  rownames_to_column(var ="ASV")

# Generate Venn diagram using filtered datasets
# Exp1
v_comb <- ps_venn(phy_16S_com, group = "Experiment", type = "counts", plot = TRUE) # creates a plot
ggsave("results/all-dataset/v_comb_16s.pdf", width = 8, height = 6, dpi = 300) # save graphic
v_comb_list <- ps_venn(phy_16S_com, group = "Experiment", type = "counts", plot = FALSE) # creates a list
class(v_comb_list)
v_comb_list_df <- ldply(v_comb_list, data.frame) # it creates a dataframe with the lists
head(v_comb_list_df)
names(v_comb_list_df) <- c("Experiment", "ASV") # rename column names
v_comb_list_df_tax <- left_join(v_comb_list_df, tax_tab_asv)
write.csv(v_comb_list_df_tax, "exported_tables/all-dataset/v_comb_list_df_tax.csv")


###################################### Subset nematodes only ###################################### 
phy_16S_com_nema <- subset_samples(phy_16S_com, Sample.Type == "SingleWorm")
unique(sample_data(phy_16S_com_nema)$Sample.Type)

any(taxa_sums(phy_16S_com_nema) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16S_com_nema) == 0) # gives the number of cases
phy_16S_com_nema_filt <- prune_taxa(taxa_sums(phy_16S_com_nema) > 0, phy_16S_com_nema)
any(taxa_sums(phy_16S_com_nema_filt) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16S_com_nema_filt) == 0) # gives the number of cases

phy_16S_com_nema_filt_clr <- microbiome::transform(phy_16S_com_nema_filt, "clr") # transform reads using centered log ratio


############### Ordination using PCA and clr transformed phyloseq objects
phy_16S_com_nema_filt_clr_pca <- phyloseq::ordinate(
  physeq = phy_16S_com_nema_filt_clr,
  method = "RDA",
  distance = "euclidean"
)

# Plot scree plot to see PCA contributions
phyloseq::plot_scree(phy_16S_com_nema_filt_clr_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/16S_Soil_ShSk_ASV_noncontam_prev05_true_all_filt_clr_pca.tiff", width = 12, height = 6, dpi = 150)

head(phy_16S_com_nema_filt_clr_pca$CA$eig)
sapply(phy_16S_com_nema_filt_clr_pca$CA$eig[1:5], function(x) x / sum(phy_16S_com_nema_filt_clr_pca$CA$eig))

# Scale axes and plot ordination
phy_16S_com_nema_filt_clr_pca_clr1 <- phy_16S_com_nema_filt_clr_pca$CA$eig[1] / sum(phy_16S_com_nema_filt_clr_pca$CA$eig)
phy_16S_com_nema_filt_clr_pca_clr2 <- phy_16S_com_nema_filt_clr_pca$CA$eig[2] / sum(phy_16S_com_nema_filt_clr_pca$CA$eig)

# Plot PCA based on clr transformation
phyloseq::plot_ordination(
  physeq = phy_16S_com_nema_filt_clr,
  ordination = phy_16S_com_nema_filt_clr_pca,
  type="samples",
  color="Habitat",
  shape = "Soil.Type") + 
  scale_color_manual(values = c("#52392F", "#C29B6C", "#83643E",
                                "#397A4C", "#77C063", "#BEDB92"),
                     name = "Habitat",
                     breaks = c("Chaparral", "Coastal.scrub", "Native.grass", "Hollyleaf.cherry", "Oak.wood", "Riparian"),
                     labels = c("CHA", "CSS", "NGR", "HLC", "OWL", "RIP"))+
  scale_shape_manual(values = c(19,17),
                     name = "Soil Type",
                     breaks = c("Compact.Soil", "Not.Compact.Soil"),
                     labels = c("Compacted", "Uncompacted"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Habitat), alpha = 0.7, size = 3) +
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
  #annotate("text", x = 3, y = 2, label ="2D Stress: 0.07") +
  #coord_fixed(phy_16_prune_soil_noncontam_prev05_true_all_filt_clr_pca_clr2 / phy_16_prune_soil_noncontam_prev05_true_all_filt_clr_pca_clr1) +
  stat_ellipse(aes(group = Habitat), linetype = 2)
ggsave("results/phy_16S_com_clr_pca_plot_soiltype_habitat.tiff", width = 6, height = 4, dpi = 200)


########### Define phyloseq objects per taxonomic group
unique(sample_data(phy_16S_com_nema_filt)$Order)

head(sample_data(phy_16S_com_nema_filt))
phy_16S_com_nema_filt_plectida <- subset_samples(phy_16S_com_nema_filt, Order == "Plectida") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16S_com_nema_filt_plectida)$Order)

# Export and Transpose matrix for ggtree analysis
asv_table_plectida = as.data.frame(otu_table(phy_16S_com_nema_filt_plectida))
asv_table_plectida <- rownames_to_column(asv_table_plectida, var = "ASV")
write.csv(asv_table_plectida, "exported_tables/asv_table_plectida.csv")

sample_table_plectida = as.data.frame(sample_data(phy_16S_com_nema_filt_plectida))
write.csv(sample_table_plectida, "exported_tables/sample_table_plectida.csv")

tax_table_plectida <- as.data.frame(tax_table(phy_16S_com_nema_filt_plectida))
tax_table_plectida <- rownames_to_column(tax_table_plectida, var = "ASV")
write.csv(tax_table_plectida, "exported_tables/tax_table_plectida.csv")


phy_16S_com_nema_filt_dorylaimida <- subset_samples(phy_16S_com_nema_filt, Order == "Dorylaimida") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16S_com_nema_filt_dorylaimida)$Order)

phy_16S_com_nema_filt_tylenchids <- subset_samples(phy_16S_com_nema_filt, Infraorder == "Tylenchomorpha") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16S_com_nema_filt_tylenchids)$Infraorder)

phy_16S_com_nema_filt_cephalobids <- subset_samples(phy_16S_com_nema_filt, Infraorder == "Cephalobomorpha") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16S_com_nema_filt_cephalobids)$Infraorder)

phy_16S_com_nema_filt_panagrolaimids <- subset_samples(phy_16S_com_nema_filt, Infraorder == "Panagrolaimomorpha") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16S_com_nema_filt_panagrolaimids)$Infraorder)
