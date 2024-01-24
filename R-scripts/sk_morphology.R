#setwd("/Users/tiagopereira/Dropbox/SK-data/Morphology-data")

library("readxl") #Can be used to import data matrices in excel format, including the different sheets
library("ggplot2") #Used for different graphics/plots styles. Very flexible to manage data
library("ggthemes") #Provides additional flexibility/features for ggplot2
library("extrafont") #Allows to use fonts other than basic PostScript fonts
library("plyr") #Provides a set of tools/functions to manipulate datasets/matrices
library("scales") #Use for graphical scales map data
library("ggpubr") #Provides some easy-to-use functions to customize ggplot2 graphics
library("RColorBrewer") #Provides color schemes for maps/graphics designed by Cynthia Brewer
library("gtable") #Provides tools/functions to make easier to work with tables
library("grid") #Rewrite graphics layout capabilities, also supporting interactions
library("cowplot") #Simple add-on to ggplot2; provides additional features to improve graphic quality.
#library("dplyr") #A tool for working/manipulating dataframes
library("vegan") #Community ecology package; it can be use for diverse multivariate analysis
library("stringr") #Provides a set of functions to work with strings
library("tidyr") #Provides tools to work with tables, dealing with columns and rows as well as individual cells. 
library("ggord") #A package for creating ordination plots with ggplot2
library("ggfortify") #Unified plotting tools for common statistics and graphics through ggplot2
library("ggrepel") #Provides geoms for ggplot2 to repel overlapping text labels
library("factoextra") #Provides functions for multivariate data analysis as well as graphics using ggplot2
library("ggbiplot") #Implements the biplot (i.e. samples data points and vectors of variables) using ggplot2
library("devtools") #Allows the installation of R packages
library("phyloseq") #Provides a series of methods for multivariate analysis from metabarcoding datasets
library("ggh4x") #allow additional ggplot manipulations

source("src/summary.R")

theme_set(theme_bw())

#Import matrices using excel files
sk_mophology_matrix <- as.data.frame(read_excel("data/summary_morph_matrix.xlsx", sheet = "matrix"))
rownames(sk_mophology_matrix) <- sk_mophology_matrix[, 1]
sk_mophology_matrix <- sk_mophology_matrix[, -1]
sk_mophology_matrix <- t(as.data.frame(sk_mophology_matrix))

sk_mophology_tax_feeding <- as.data.frame(read_excel("data/summary_morph_matrix.xlsx", sheet = "taxonomy-feeding"))
rownames(sk_mophology_tax_feeding) <- sk_mophology_tax_feeding[, 1]
sk_mophology_tax_feeding <- sk_mophology_tax_feeding[, -1]
sk_mophology_tax_feeding <- as.matrix(sk_mophology_tax_feeding)

sk_mophology_mapping <- as.data.frame(read_excel("data/summary_morph_matrix.xlsx", sheet = "mapping"))
rownames(sk_mophology_mapping) <- sk_mophology_mapping[, 1]

#Transform to phyloseq objects
specimens = otu_table(sk_mophology_matrix, taxa_are_rows = TRUE)
taxonomy = tax_table(sk_mophology_tax_feeding)
samples = sample_data(sk_mophology_mapping)

# Creating a phyloseq object
morpho_phy <- phyloseq(specimens, taxonomy, samples)
morpho_phy
head(sample_data(morpho_phy))

# Transform abundances to relative abudance
morpho_phy_comp <- microbiome::transform(morpho_phy, "compositional")


# Ordinate NMDS
set.seed(1)
morpho_phy_comp_nmds <- ordinate(
  physeq = morpho_phy_comp, 
  method = "NMDS", 
  distance = "bray",
  k = 3,
  try =50
)

# Plot NMDS
theme_set(theme_bw())

morpho_comp_nmds_plot <- 
  phyloseq::plot_ordination(
  physeq = morpho_phy_comp,
  ordination = morpho_phy_comp_nmds,
  type = "samples",
  color = "Habitat",
  shape = "SoilType") +
  scale_color_manual(values = c("#52392F", "#83643E", "#C29B6C", 
                                "#397A4C", "#77C063", "#BEDB92"),
                     name = "Habitat",
                     breaks = c("Chaparral", "Coastal_scrub", "Native_grass", "Hollyleaf_cherry", "Oak_wood", "Riparian"),
                     labels = c("CHA", "CCS", "NGR", "HLC", "OWL", "RIP"))+
  scale_shape_manual(values = c(19,17),
                     name = "Soil Type",
                     breaks = c("Compacted", "Uncompacted"),
                     labels = c("Compacted", "Uncompacted"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("NMDS2") + # add the title on y axis
  #xlab("NMDS1") + # add the title on x axis
  theme(strip.background =element_blank()) + # remove the background of titles
  annotate("text", x = 0.7, y = 0.9, label ="2D Stress: 0.17") +
  #ggtitle("Nematode Morphology") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =12))
  morpho_comp_nmds_plot
  ggsave("results/morphology/morpho_comp_nmds.pdf", width = 6, height = 4, dpi = 300)
  
######## Combine nMDS plots #########
# Needs to use the morphological R script
Figure_nMDS <- ggarrange(morpho_comp_nmds_plot,
                           phy_18S_prune_baermann_noncontam_prev05_true_scale_filt_comp_nmds_plot,
                           phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nmds_plot,
                           common.legend = TRUE,  legend = "right", ncol = 3, nrow = 1, labels = "AUTO")
Figure_nMDS
ggsave("results/morphology/morpho_16S_soil_18S_baemann_18S_soil_comp_nmds.tiff", width = 14, height = 5, dpi = 200)

# Function to collapse a certain number of taxa into category others
merge_less_than_top_prev05_top20 <- function(morpho_phy, top=19){
  transformed <- transform_sample_counts(morpho_phy, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:13] <- "Others"} # 1:13 for nematode taxonomy including feedig group and c-p value
  }
  return(merged)
}

# Family level top 20
#morpho_phy_family <- tax_glom(morpho_phy, "Family")
#morpho_phy_family_top20 <- merge_less_than_top_prev05_top20(morpho_phy_family, top=19)
#morpho_phy_family_top20_df = psmelt(morpho_phy_family_top20) # Melt to long format
#morpho_phy_family_top20_agr_site = aggregate(Abundance~SoilType+Code+Sample+Family, data=morpho_phy_family_top20_df, FUN=mean) # add common factors to use for plotting
#unique(morpho_phy_family_top20_agr_site$Family)


# Melt to long format (for ggplot2)
  morpho_phy_family <- morpho_phy %>%
    tax_glom(taxrank = "Family") %>%                      # agglomerate at family level
    transform_sample_counts(function(x) {x/sum(x)} ) %>%  # Transform to rel. abundance
    psmelt() %>%                                          # Melt to long format for ggplot2
    #filter(Abundance > 0) %>%                             # Filter out taxa with 0 abundance
    arrange(Family)                                       # Sort data frame alphabetically by family name
  morpho_phy_family
  str(morpho_phy_family)

# Genus level top 20
morpho_phy_genus <- tax_glom(morpho_phy, "Genus")
morpho_phy_genus_top20 <- merge_less_than_top_prev05_top20(morpho_phy_genus, top=19)
morpho_phy_genus_top20_df = psmelt(morpho_phy_genus_top20) # Melt to long format
morpho_phy_genus_top20_agr_site = aggregate(Abundance~SoilType+Code+Sample+Genus, data=morpho_phy_genus_top20_df, FUN=mean) # add common factors to use for plotting
unique(morpho_phy_genus_top20_agr_site$Genus)
  
# Melt to long format (for ggplot2)
#  morpho_phy_genus <- morpho_phy %>%
#    tax_glom(taxrank = "Genus") %>%                      # agglomerate at family level
#    transform_sample_counts(function(x) {x/sum(x)} ) %>%  # Transform to rel. abundance
#    psmelt() %>%                                          # Melt to long format for ggplot2
    #filter(Abundance > 0) %>%                             # Filter out taxa with 0 abundance
#    arrange(Genus)                                       # Sort data frame alphabetically by family name
#  morpho_phy_genus
#  str(morpho_phy_genus)

# Melt to long format (for ggplot2)
morpho_phy_feeding <- morpho_phy %>%
    tax_glom(taxrank = "Feeding.Group") %>%                      # agglomerate at family level
    transform_sample_counts(function(x) {x/sum(x)} ) %>%  # Transform to rel. abundance
    psmelt() %>%                                          # Melt to long format for ggplot2
    #filter(Abundance > 0) %>%                             # Filter out taxa with 0 abundance
    arrange(Feeding.Group)                                       # Sort data frame alphabetically by family name
morpho_phy_feeding

write.csv(morpho_phy_feeding, "exported_tables/morpho_phy_feeding.csv")  
#Get unique names for family
unique(morpho_phy_family$Family)

#Get unique names for genus
unique(morpho_phy_genus$Genus)

#Get unique names for feeding groups
unique(morpho_phy_feeding$Feeding)

#nb.cols.family.morph <- 21
#soil_colors_family <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.family.morph)

# List of 20 distinct colors
colors_top21 <- c("#9a6324", "#46f0f0", "#1F78B4", "#aaffc3", "#e6beff",  "#33A02C", "#4363d8", "#008080", "#FB9A99", "#e6194b",
                         "#f032e6", "#bcf60c", "#fabebe","#f58231", "#ffe119", "#6A3D9A", "#000075", "#fffac8", "#800000", "black", "gray")

colors_top20 <- c("#9a6324", "#46f0f0", "#1F78B4", "#aaffc3", "#e6beff",  "#33A02C", "#4363d8", "#008080", "#FB9A99", "#e6194b",
                  "#f032e6", "#bcf60c", "#fabebe","#f58231", "#ffe119", "#6A3D9A", "#000075", "#fffac8", "#800000","gray")

colors_top12 <- c("#f032e6", "#46f0f0", "#1F78B4", "#aaffc3", "#e6beff",  "#33A02C", "#4363d8", "#008080", "#FB9A99", "#e6194b",
                  "#bcf60c",  "gray", "#fabebe","#f58231", "#ffe119", "#6A3D9A", "#000075", "#fffac8", "#800000", "black")


colors_feeding <- c("#f032e6", "#46f0f0", "#bcf60c", "#aaffc3", "#e6beff", "gray")

#nb.cols.feeding.morph <- 6
#soil_colors_feeding <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.feeding.morph)

habitat.labs <- c("CHA", "CCS", "NGR", "HLC", "OWL", "RIP") # create new names for the variable Habitat
names(habitat.labs) <- c("Chaparral", "Coastal_scrub", "Native_grass", "Hollyleaf_cherry", "Oak_wood", "Riparian") # assign old name to new value
                   
# Plot bar chart morphological data at Family level
morpho_phy_fam_plot <- ggplot(morpho_phy_family, aes(x = Sample, y = Abundance, fill = Family)) + #plotting by sample
    facet_nested(. ~ SoilType+Code, scales = "free_x") + # By Habitat
    geom_bar(stat = "identity", position = "fill") + # adds to 100%
    geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = colors_top21) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9)) + # adjusts text of y axis
    theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
    theme(axis.title.x = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis  
    scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # removes the internal margins
    scale_x_discrete(
     breaks = c("SK01", "SK02", "SK03", "SK04", "SK05", "SK06",
                "SK07", "SK08", "SK09", "SK10", "SK11", "SK12",
                "SK13", "SK14", "SK15", "SK16", "SK17", "SK18"),
     labels = c("01", "02", "03", "04", "05", "06",
                "07", "08", "09", "10", "11", "12",
                "13", "14", "15", "16", "17", "18"),
     expand = c(0.0, 0.0),
     drop = FALSE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
    guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
    theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
    ylab("Relative Abundance") + # add the title on y axis
    xlab("Site") + # add the title on x axis
    theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
    theme(strip.background =element_rect("white")) + # remove the background of titles
    theme(strip.text.x = element_text(
      size = 12, color = "black", face = "bold"))
morpho_phy_fam_plot
ggsave("results/morphology/morpho_phy_family_barplot.pdf", width = 10, height = 5, dpi = 300) # save graphic  

# Plot bar chart morphological data at genus level
morpho_phy_genus_top20_agr_site$Genus <- factor(morpho_phy_genus_top20_agr_site$Genus,
                                                levels = c("Acrobeles", "Acrobeloides", "Aphelenchoides", "Aphelenchus",
                                                           "Aporcelaimellus", "Aporcella", "Cervidellus", "Chiloplacus",
                                                           "Crassolabium", "Ditylenchus", "Filenchus", "Merlinius",
                                                           "Mesorhabditis",  "Panagrolaimus", "Paraphelenchus", "Plectus",
                                                           "Tylenchorhynchus", "Tylenchus", "Tylocephalus", "Others"))

morpho_phy_genus_plot <- ggplot(morpho_phy_genus_top20_agr_site, aes(x = Sample, y = Abundance, fill = Genus)) + #plotting by sample
    facet_nested(. ~ SoilType+Code, scales = "free_x") + # By Habitat
    geom_bar(stat = "identity", position = "fill") + # adds to 100%
    geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = colors_top20) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9)) + # adjusts text of y axis
    theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
    theme(axis.title.x = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
    scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # removes the internal margins
    scale_x_discrete(
      breaks = c("SK01", "SK02", "SK03", "SK04", "SK05", "SK06",
                 "SK07", "SK08", "SK09", "SK10", "SK11", "SK12",
                 "SK13", "SK14", "SK15", "SK16", "SK17", "SK18"),
      labels = c("01", "02", "03", "04", "05", "06",
                 "07", "08", "09", "10", "11", "12",
                 "13", "14", "15", "16", "17", "18"),
      expand = c(0.0, 0.0),
      drop = FALSE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
    guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
    theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
    ylab("Relative Abundance") + # add the title on y axis
    xlab("Site") + # add the title on x axis
    theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
    theme(strip.background =element_rect("white")) + # remove the background of titles
    theme(strip.text.x = element_text(
      size = 12, color = "black", face = "bold"))
  morpho_phy_genus_plot
  ggsave("results/morphology/morpho_phy_genus_barplot.pdf", width = 10, height = 5, dpi = 300) # save graphic  
  

# Plot bar chart morphological data for feeding groups
morpho_phy_feeding_plot <- ggplot(morpho_phy_feeding, aes(x = Sample, y = Abundance, fill = Feeding.Group)) + #plotting by sample
  facet_nested(. ~ SoilType+Code, scales = "free_x") + # By Habitat
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = colors_feeding) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
    theme(axis.text.x = element_text(angle = 0.5, hjust = 1, size = 9)) + # adjusts text of x axis
    theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
    theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
    scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # removes the internal margins
    scale_x_discrete(
      breaks = c("SK01", "SK02", "SK03", "SK04", "SK05", "SK06",
                 "SK07", "SK08", "SK09", "SK10", "SK11", "SK12",
                 "SK13", "SK14", "SK15", "SK16", "SK17", "SK18"),
      labels = c("01", "02", "03", "04", "05", "06",
                 "07", "08", "09", "10", "11", "12",
                 "13", "14", "15", "16", "17", "18"),
      expand = c(0.0, 0.0),
      drop = FALSE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
    guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
    labs(fill = "Feeding group")+
    theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
    ylab("Relative Abundance") + # add the title on y axis
    xlab("Site") + # add the title on x axis
    theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
    theme(strip.background =element_rect("white")) + # remove the background of titles
    theme(strip.text.x = element_text(
      size = 12, color = "black", face = "bold"))
morpho_phy_feeding_plot
ggsave("results/morphology/morpho_phy_feeding_barplot.pdf", width = 10, height = 5, dpi = 300) # save graphic  
  
phy_morphy_combined <- ggarrange(morpho_phy_fam_plot, morpho_phy_genus_plot,
                                   morpho_phy_feeding_plot, legend = "right", ncol = 1, nrow = 3, labels = "AUTO",
                                          align = "v")
phy_morphy_combined
ggsave("results/morphology/Nematoda_phy_morphy_Soil_Habitat_ShSk_barplot_combined.pdf", width = 10, height = 14, dpi = 150)
  
# Calculate summary (i.e., mean and standard error) using summarySE function from summary.R file
# Summaries per habitat and soil type
sk_ph_sum_habitat <- summarySE(sk_ph, measurevar=c("pH"), groupvars=c("SoilType", "Habitat", "Code", "Type"))
sk_ph_sum_habitat

sk_soil_sum_habitat <- summarySE(sk_soil, measurevar=c("Perc"), groupvars=c("SoilType", "Habitat","Code", "Class"))
sk_soil_sum_habitat

sk_organic_sum_habitat <- summarySE(sk_organic, measurevar=c("Perc"), groupvars=c("SoilType", "Habitat","Code", "Type"))
sk_organic_sum_habitat

sk_nutrients_sum_habitat <- summarySE(sk_nutrients, measurevar=c("ppm"), groupvars=c("SoilType", "Habitat","Code", "Type"))
sk_nutrients_sum_habitat

# Summaries per soil type
sk_ph_sum_soiltype <- summarySE(sk_ph, measurevar=c("pH"), groupvars=c("SoilType", "Type"))
sk_ph_sum_soiltype

sk_soil_sum_soiltype <- summarySE(sk_soil, measurevar=c("Perc"), groupvars=c("SoilType", "Class"))
sk_soil_sum_soiltype

sk_organic_sum_soiltype <- summarySE(sk_organic, measurevar=c("Perc"), groupvars=c("SoilType", "Type"))
sk_organic_sum_soiltype

sk_nutrients_sum_soiltype <- summarySE(sk_nutrients, measurevar=c("ppm"), groupvars=c("SoilType", "Type"))
sk_nutrients_sum_soiltype

# Create plots for each variable using ggplot with facet (Soil type x habitat)
sk_soil_sum_habitat_plot <- ggplot(sk_soil_sum_habitat, aes(x= Code, y= Perc, fill = Class)) +
  geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=Perc-se, ymax=Perc+se), size = 0.5, width = 0.2, position=position_dodge(0.9)) +
  geom_point(size = 2, shape = 21, position=position_dodge(0.9)) +
  scale_fill_brewer(palette = "Blues") +
  facet_grid(Class ~ SoilType, scales = "free") +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  theme(strip.text.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.x = element_text(face = "bold", size = 14)) + # adjusts the title of x axis
  theme(axis.title.y = element_blank()) + # adjusts the title of x axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12, color = "black")) + # adjusts text of x axis
  scale_y_continuous(expand = c(0.03, 0.03)) +
  #scale_x_discrete(expand = c(0.16, 0.16)) +
  xlab("Habitat") + # add the title on y axis
  ggtitle("Soil Texture (%)") +
  theme(plot.title = element_text(face = "bold", size = 14, colour = "black", hjust = 0.5)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(1, "lines")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
sk_soil_sum_habitat_plot

sk_ph_sum_habitat_plot <- ggplot(sk_ph_sum_habitat, aes(x= Code, y= pH)) +
  geom_bar(position=position_dodge(), stat = "identity", color = "gray", size = 0.5) +
  geom_errorbar(aes(ymin=pH-se, ymax=pH+se), size = 0.5, width = 0.2, position=position_dodge(width = 0.9)) +
  geom_point(size = 2, shape = 21, fill = "gray") +
  facet_grid(Type ~ SoilType, scales = "free") +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  theme(strip.text.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.x = element_text(face = "bold", size = 14)) + # adjusts the title of x axis
  theme(axis.title.y = element_blank()) + # adjusts the title of x axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12, color = "black")) + # adjusts text of x axis
  scale_y_continuous(expand = c(0.02, 0.02)) +
  xlab("Habitat") + # add the title on y axis
  ggtitle("pH") +
  theme(plot.title = element_text(face = "bold", size = 14, colour = "black", hjust = 0.5)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(1, "lines")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
sk_ph_sum_habitat_plot

sk_nutrients_sum_habitat_plot <- ggplot(sk_nutrients_sum_habitat, aes(x= Code, y= ppm, fill = Type)) +
  geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=ppm-se, ymax=ppm+se), size = 0.5, width = 0.2, position=position_dodge(width = 0.9)) +
  geom_point(size = 2, shape = 21, position=position_dodge(width = 0.9)) +
  scale_fill_brewer(palette = "YlGn") +
  facet_grid(Type ~ SoilType, scales = "free") +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  theme(strip.text.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.x = element_text(face = "bold", size = 14)) + # adjusts the title of x axis
  theme(axis.title.y = element_blank()) + # adjusts the title of x axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12, color = "black")) + # adjusts text of x axis
  scale_y_continuous(expand = c(0.02, 0.02)) +
  xlab("Habitat") + # add the title on y axis
  ggtitle("Nutrients (ppm)") +
  theme(plot.title = element_text(face = "bold", size = 14, colour = "black", hjust = 0.5)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(1, "lines")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
sk_nutrients_sum_habitat_plot

sk_organic_sum_habitat_plot <- ggplot(sk_organic_sum_habitat, aes(x= Code, y= Perc, fill = Type)) +
  geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=Perc-se, ymax=Perc+se), size = 0.5, width = 0.2, position=position_dodge(width = 0.9)) +
  geom_point(size = 2, shape = 21, position=position_dodge(width = 0.9)) +
  scale_fill_brewer(palette = "PuOr") +
  facet_grid(Type ~ SoilType, scales = "free") +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  theme(strip.text.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.y = element_blank()) + # adjusts the title of x axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12, color = "black")) + # adjusts text of x axis
  scale_y_continuous(expand = c(0.02, 0.02)) +
  xlab("Habitat") + # add the title on y axis
  ggtitle("Organics (%) & pH") +
  theme(plot.title = element_text(face = "bold", size = 14, colour = "black", hjust = 0.5)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(1, "lines")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
sk_organic_sum_habitat_plot

environ_var_habitat <- ggarrange(sk_soil_sum_habitat_plot, sk_nutrients_sum_habitat_plot, sk_organic_sum_habitat_plot, 
                                  ncol = 3, nrow = 1, legend = "none",
                                  align = "hv", labels = "AUTO")
environ_var_habitat
ggsave("FigureS1.tiff", width = 10, height =5, dpi = 150) # save graphic

sk_soil_sum_soiltype_plot <- ggplot(sk_soil_sum_soiltype, aes(x= SoilType, y= Perc, fill = Class)) +
  geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=Perc-se, ymax=Perc+se), size = 0.5, width = 0.2, position=position_dodge(0.9)) +
  geom_point(size = 2, shape = 21, position=position_dodge(0.9)) +
  scale_fill_brewer(palette = "Blues") +
  facet_grid(Class ~ SoilType, scales = "free") +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  theme(strip.text.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  scale_y_continuous(expand = c(0.02, 0.02)) +
  xlab("Soil Type") + # add the title on y axis
  ggtitle("Soil Texture (%)") +
  theme(plot.title = element_text(face = "bold", size = 14, colour = "black", hjust = 0.5)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(1, "lines")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
sk_soil_sum_soiltype_plot

sk_ph_sum_soiltype_plot <- ggplot(sk_ph_sum_soiltype, aes(x= SoilType, y= pH)) +
  geom_bar(position=position_dodge(), stat = "identity", color = "gray", size = 0.5) +
  geom_errorbar(aes(ymin=pH-se, ymax=pH+se), size = 0.5, width = 0.2, position=position_dodge(width = 0.9)) +
  geom_point(size = 2, shape = 21, fill = "gray") +
  facet_grid(Type ~ SoilType, scales = "free") +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  theme(strip.text.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  scale_y_continuous(expand = c(0.02, 0.02)) +
  xlab("Soil Type") + # add the title on y axis
  theme(legend.position="right") +
  ggtitle("pH") +
  theme(plot.title = element_text(face = "bold", size = 14, colour = "black", hjust = 0.5)) +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(1, "lines")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
sk_ph_sum_soiltype_plot

sk_nutrients_sum_soiltype_plot <- ggplot(sk_nutrients_sum_soiltype, aes(x= SoilType, y= ppm, fill = Type)) +
  geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=ppm-se, ymax=ppm+se), size = 0.5, width = 0.2, position=position_dodge(width = 0.9)) +
  geom_point(size = 2, shape = 21, position=position_dodge(width = 0.9)) +
  scale_fill_brewer(palette = "YlGn") +
  facet_grid(Type ~ SoilType, scales = "free") +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  theme(strip.text.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  scale_y_continuous(expand = c(0.02, 0.02)) +
  xlab("Soil Type") + # add the title on y axis
  ggtitle("Nutrients (ppm)") +
  theme(plot.title = element_text(face = "bold", size = 14, colour = "black", hjust = 0.5)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(1, "lines")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
sk_nutrients_sum_soiltype_plot

sk_organic_sum_soiltype_plot <- ggplot(sk_organic_sum_soiltype, aes(x= SoilType, y= Perc, fill = Type)) +
  geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=Perc-se, ymax=Perc+se), size = 0.5, width = 0.2, position=position_dodge(width = 0.9)) +
  geom_point(size = 2, shape = 21, position=position_dodge(width = 0.9)) +
  scale_fill_brewer(palette = "PuOr") +
  facet_grid(Type ~ SoilType, scales = "free") +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  theme(strip.text.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  scale_y_continuous(expand = c(0.02, 0.02)) +
  theme(legend.position="right") +
  ggtitle("Organic Material (%) & pH") +
  theme(plot.title = element_text(face = "bold", size = 14, colour = "black", hjust = 0.5)) +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(1, "lines")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
sk_organic_sum_soiltype_plot

environ_var_soiltype <- ggarrange(sk_soil_sum_soiltype_plot, sk_nutrients_sum_soiltype_plot, sk_organic_sum_soiltype_plot, 
                                 ncol = 3, nrow = 1, legend = "none",
                                 align = "hv", labels = "AUTO")
environ_var_soiltype
ggsave("FigureS2.tiff", width = 10, height =5, dpi = 150) # save graphic

# Principal component analysis using Vegan and RDA
# First, organize matrix by removing factors and renaming rows.
sk_pca_matrix <- as.data.frame(sk_abiotic[, -(2:5)])
str(sk_pca_matrix)
rownames(sk_pca_matrix) <- sk_pca_matrix[,1]
sk_pca_matrix  <- sk_pca_matrix[, -1]

# Summarize PCA analysis
sk_pca_1 <- rda(sk_pca_matrix, scale = TRUE)
summary(sk_pca_1)

# Checking the eigenvalues
sk_pca_1$CA$eig

# Expressed as cumulative percentages
round(cumsum(100*sk_pca_1$CA$eig/sum(sk_pca_1$CA$eig)),2)

# Create screeplot
screeplot(sk_pca_1)
abline(a=1,b=0)

# Plot PCA analysis plot
biplot(sk_pca_1)

# Principal component analysis using prcomp function
sk_pca_2 <- prcomp(sk_pca_matrix, scale. = TRUE)
summary(sk_pca_2)

autoplot(sk_pca_2, data = sk_abiotic, colour = "Habitat", shape = "SoilType",
         loadings = TRUE, loadings.colour = 'lightgray', loadings.label = TRUE, loadings.label.size = 4, loadings.label.colour =  "black",
         repel =  TRUE, cex = 0.5) +
  scale_color_manual(values = c("#52392F", "#C29B6C", "#83643E",
                                "#397A4C", "#77C063", "#BEDB92"),
                     name = "Habitat",
                     breaks = c("Chaparral", "Coastal_sage_scrub", "Native_grass", "Hollyleaf_Cherry", "Oak_woodland", "Riparian"),
                     labels = c("Chaparral", "Coastal Sage Scrub", "Native Grass", "Hollyleaf Cherry", "Oak Woodland", "Riparian"))+
  scale_shape_manual(values = c(19,17),
                     name = "Soil Type",
                     breaks = c("Compacted", "Uncompacted"),
                     labels = c("Compacted", "Uncompacted"))+
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(aes(color =  Habitat, shape = SoilType), alpha = 0.7, size = 4) +
  ylab("PC2 (29%)") + # add the title on y axis
  xlab("PC1 (57%)") + # add the title on x axis
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
ggsave("PCA_soiltype_habitat.pdf", width = 6, height =4) # save graphic as pdf for final adjustments using Inkscape/Illustrator
ggsave("PCA_soiltype_habitat.tiff", width = 6, height =4) # save graphic as tiff

#Non-parametric one-way ANOVA - Kruskal-Wallis rank sum test by habitat and soil type
kt_pH_soiltype <- kruskal.test(pH ~ SoilType, data = sk_abiotic) # KW-test by habitat
kt_pH_habitat <- kruskal.test(pH ~ Code, data = sk_abiotic) # KW-test by habitat
parwise_pH_habitat <- pairwise.wilcox.test(sk_abiotic$pH, sk_abiotic$Code,
                                        p.adjust.method = "bonferroni")
write.table(parwise_asv_shannon[["p.value"]], file="parwise_asv_shannon.csv", sep=",")
kruskal.test(H ~ FG, data = kt_diverse_shannon) # KW-test by feeding group


kt_Sand <- kruskal.test(Sand ~ Code, data = sk_abiotic) # KW-test by habitat
parwise_sand <- pairwise.wilcox.test(sk_abiotic$Sand, sk_abiotic$Code,
                                   p.adjust.method = "none")

meio_basin <- read_excel("uni-meio-nema.xlsx",  sheet = "Meio")
nema_basin <- read_excel("uni-meio-nema.xlsx",  sheet = "Nema")
nema_meio_basin <- read_excel("uni-meio-nema.xlsx",  sheet = "nema-meio")
nema_nmds_basin <- read_excel("nema-nmds.xlsx", sheet = "Basin")
nema_nmds_depth <- read_excel("nema-nmds.xlsx", sheet = "Depth")
nema_feeding_basin <- read_excel("feeding.xlsx",  sheet = "Feeding")
nema_trophic_basin <- read_excel("feeding.xlsx",  sheet = "Trophic")
nema_depth <- read_excel("nema-univar.xlsx")
nema_depth_j <- read_excel("nema-univar.xlsx", sheet = "J")
nema_feeding_depth <- read_excel("nema-univar.xlsx", sheet = "Feeding")
nema_trophic_depth <- read_excel("nema-univar.xlsx", sheet = "Trophic")

# Checking imported files
head(meio_basin)
head(nema_basin)
head(nema_nmds_basin)
head(nema_feeding_basin)
head(nema_trophic_basin)
head(nema_depth)

# Summary +/-SE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
# Summary for meiofauna across basins
meio_sum_taxa_basin <- summarySE(meio_basin, measurevar=c("Ntaxa"), groupvars=c("Basin"))
meio_sum_taxa_basin

meio_sum_density_basin <- summarySE(meio_basin, measurevar=c("Density10"), groupvars=c("Basin"))
meio_sum_density_basin
############
# Summary +/-SE for nematodes across basins
nema_sum_gen_basin <- summarySE(nema_basin, measurevar="Ngenera", groupvars=c("Basin"))
nema_sum_gen_basin

nema_sum_density_basin <- summarySE(nema_basin, measurevar="Density10", groupvars=c("Basin"))
nema_sum_density_basin

nema_sum_es_basin <- summarySE(nema_basin, measurevar="ES51", groupvars=c("Basin"))
nema_sum_es_basin

nema_sum_j_basin <- summarySE(nema_basin, measurevar="J", groupvars=c("Basin"))
nema_sum_j_basin

nema_sum_h_basin <- summarySE(nema_basin, measurevar="H", groupvars=c("Basin"))
nema_sum_h_basin

nema_sum_feeding_basin <- summarySE(nema_feeding_basin, measurevar="Abundance", groupvars=c("Basin", "FeedingGroup"))
nema_sum_feeding_basin

nema_suma_trophic_basin <- summarySE(nema_trophic_basin, measurevar="TropDiv", groupvars=c("Basin"))
nema_suma_trophic_basin

#################
# Summary +/-SE for meiofauna and nematodes across basins in the same plot
nema_meio_sum_density_basin <- summarySE(nema_meio_basin, measurevar=c("Density10"), groupvars=c("Basin", "Group"))
nema_meio_sum_density_basin

nema_meio_sum_taxa_basin <- summarySE(nema_meio_basin, measurevar=c("Ntaxa"), groupvars=c("Basin", "Group"))
nema_meio_sum_taxa_basin

#################
nema_sum_trophic_depth <- summarySE(nema_trophic_depth, measurevar="TropDiv", groupvars=c("Basin", "DepthZone"))
nema_sum_trophic_depth

nema_sum_feeding_depth <- summarySE(nema_feeding_depth, measurevar="Abundance", groupvars=c("Basin", "FeedingGroup", "DepthZone"))
nema_sum_feeding_depth

nema_sum_s_depth <- summarySE(nema_depth, measurevar="S", groupvars=c("Basin", "DepthZone"))
nema_sum_s_depth

nema_sum_n_depth <- summarySE(nema_depth, measurevar="N10", groupvars=c("Basin", "DepthZone"))
nema_sum_n_depth

nema_sum_es_depth <- summarySE(nema_depth, measurevar="ES51", groupvars=c("Basin", "DepthZone"))
nema_sum_es_depth

nema_sum_j_depth <- summarySE(nema_depth_j, measurevar="J", groupvars=c("Basin", "DepthZone"))
nema_sum_j_depth

#################
# Barplots of nematodes across depths, Amazon and Ceara only
nema_s_depth_plot <- ggplot(nema_sum_s_depth, aes(x= DepthZone, y= S, fill = Basin)) + 
        geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
        geom_errorbar(aes(ymin=S-se, ymax=S+se), size = 0.5, width = 0.2, position=position_dodge(width = 0.9)) +
        geom_point(size = 2, shape = 21, fill = "#D95F02") +
        scale_fill_manual(values=c("#1F78B4", "#33A02C")) +
        facet_wrap(~ Basin, scales = "free_x") +
        theme(strip.text.x = element_text(size = 14, color = "black", face = "bold")) +
        theme(axis.title.y = element_text(face = "bold", size = 16)) +  # adjusts the title of y axis
        #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black")) + # adjusts text of y axis
        #theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, color = "black")) + # adjusts text of x axis
        scale_y_continuous(expand = c(0, 0)) +
        expand_limits(y=c(0, 22)) +
        scale_x_discrete(limits=c("MB", "LB", "UA")) +
        ylab("Number of genera") + # add the title on y axis
        xlab("Depth") + # add the title on y axis
        theme(legend.position="right") +
        theme(legend.title = element_text(face = "bold")) +
        theme(legend.title.align = 0.5) +
        theme(axis.line = element_line(color = "black"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
nema_s_depth_plot

nema_n_depth_plot <- ggplot(nema_sum_n_depth, aes(x= DepthZone, y= N10, fill = Basin)) + 
        geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
        geom_errorbar(aes(ymin=N10-se, ymax=N10+se), size =  0.5, width = 0.2, position=position_dodge(0.9)) +
        geom_point(size = 2, shape = 21, fill = "#D95F02") +
        ylab(expression(bold(paste(Density~(inds.10~cm^-2))))) + # format title on y axis
        scale_fill_manual(values=c("#1F78B4", "#33A02C")) +
        facet_wrap(~ Basin, scales = "free_x") +
        theme(strip.text.x = element_text(size = 14, color = "black", face = "bold")) +
        theme(axis.title.y = element_text(face = "bold", size = 16)) +  # adjusts the title of y axis
        #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 16, color = "black")) + # adjusts text of y axis
        #theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, color = "black")) + # adjusts text of x axis
        scale_y_continuous(expand = c(0, 0)) +
        expand_limits(y=c(0, 620)) +
        scale_x_discrete(limits=c("MB", "LB", "UA")) +
        #ylab("S") + # add the title on y axis
        xlab("Depth") + # add the title on y axis
        theme(legend.position="right") +
        theme(legend.title = element_text(face = "bold")) +
        theme(legend.title.align = 0.5) +
        theme(axis.line = element_line(color = "black"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
nema_n_depth_plot

nema_es_depth_plot <- ggplot(nema_sum_es_depth, aes(x= DepthZone, y= ES51, fill = Basin)) + 
        geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
        geom_errorbar(aes(ymin=ES51-se, ymax=ES51+se), size =  0.5, width = 0.2, position=position_dodge(0.9)) +
        geom_point(size = 2, shape = 21, fill = "#D95F02") +
        scale_fill_manual(values=c("#1F78B4", "#33A02C")) +
        facet_wrap(~ Basin, scales = "free_x") +
        theme(strip.background = element_blank(),strip.text.x = element_blank()) +
        theme(axis.title.y = element_text(face = "bold", size = 16)) +  # adjusts the title of y axis
        theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
        theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black")) + # adjusts text of y axis
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, color = "black")) + # adjusts text of x axis
        scale_y_continuous(expand = c(0, 0)) +
        expand_limits(y=c(0, 19)) +
        scale_x_discrete(limits=c("MB", "LB", "UA")) +
        ylab("ES51") + # add the title on y axis
        xlab("Depth Zone") + # add the title on y axis
        theme(legend.position="right") +
        theme(legend.title = element_text(face = "bold")) +
        theme(legend.title.align = 0.5) +
        theme(axis.line = element_line(color = "black"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
nema_es_depth_plot 

nema_j_depth_plot <- ggplot(nema_sum_j_depth, aes(x= DepthZone, y= J, fill = Basin)) + 
        geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
        geom_errorbar(aes(ymin=J-se, ymax=J+se), size =  0.5, width = 0.2, position=position_dodge(0.9)) +
        geom_point(size = 2, shape = 21, fill = "#D95F02") +
        scale_fill_manual(values=c("#1F78B4", "#33A02C")) +
        facet_wrap(~ Basin, scales = "free_x") +
        theme(strip.background = element_blank(),strip.text.x = element_blank()) +
        theme(axis.title.y = element_text(face = "bold", size = 16)) +  # adjusts the title of y axis
        theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
        theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black")) + # adjusts text of y axis
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, color = "black")) + # adjusts text of x axis
        scale_y_continuous(expand = c(0, 0)) +
        expand_limits(y=c(0, 1.02)) +
        ylab("J") + # add the title on y axis
        xlab("Depth Zone") + # add the title on y axis
        scale_x_discrete(limits=c("MB", "LB", "UA")) +
        theme(legend.position="right") +
        theme(legend.title = element_text(face = "bold")) +
        theme(legend.title.align = 0.5) +
        theme(axis.line = element_line(color = "black"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
nema_j_depth_plot

nema_feeding_depth_plot <- ggplot(nema_sum_feeding_depth, aes(x= DepthZone, y= Abundance, fill = FeedingGroup)) + 
  geom_bar(stat = "identity", width= 0.8, color = "black", size = 0.5) +
  #geom_point(size = 2, shape = 21, fill = "#D95F02") +
  scale_fill_manual(values = c("#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F")) +
  facet_wrap(~ Basin, scales = "free_x") +
  theme(strip.text.x = element_text(size = 14, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 16)) +  # adjusts the title of y axis
  #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black")) + # adjusts text of y axis
  #theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, color = "black")) + # adjusts text of x axis
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(limits=c("MB", "LB", "UA")) +
  ylab("Feeding Type (%)") + # add the title on y axis
  xlab("Depth zone") + # add the title on y axis
  theme(legend.position="right") +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0.5) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
nema_feeding_depth_plot

nema_trophic_depth_plot <- ggplot(nema_sum_trophic_depth, aes(x= DepthZone, y= TropDiv)) +
  geom_bar(stat = "identity", width= 0.8, fill =  "lightblue", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=TropDiv-se, ymax=TropDiv+se), width= 0.2, size =  0.5) +
  geom_point(size = 2, shape = 21, fill = "#D95F02") +
  facet_wrap(~ Basin, scales = "free_x") +
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 16)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, color = "black")) + # adjusts text of x axis
  scale_y_continuous(position =  "left", expand = c(0, 0)) +
  expand_limits(y=c(0, 0.8)) +
  scale_x_discrete(limits=c("MB", "LB", "UA")) +
  #annotate("text", x = 1, y = 0.7, label = "*", size = 10, fontface = "bold") +
  #annotate("text", x = 4, y = 0.45, label = "*", size = 10, fontface = "bold") +
  ylab("Nematode ITD") + # add the title on y axis
  xlab("Depth zone") + # add the title on y axis
  theme(legend.position="right") +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0.5) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
nema_trophic_depth_plot

Figure4 <- ggarrange(nema_s_depth_plot, nema_n_depth_plot, nema_feeding_depth_plot,
                     nema_es_depth_plot, nema_j_depth_plot, nema_trophic_depth_plot,
                     labels = c("A", "B", "C", "D", "E", "F"), nrow = 2, ncol = 3, legend = FALSE)
Figure4
ggsave("Figure4.tiff", width = 12, height = 6, dpi = 150) # save graphic
######################################

# Barplots of standard error of the mean for meiofauna
meio_taxa_basin_plot <- ggplot(meio_sum_taxa_basin, aes(x= Basin, y= Ntaxa)) + 
        geom_bar(stat = "identity", width= 0.8, fill = "#1F78B4", color = "black", size = 0.5) +
        geom_errorbar(aes(ymin=Ntaxa-se, ymax=Ntaxa+se), width= 0.2, size =  0.5) +
        geom_point(size = 2, shape = 21, fill = "#D95F02") +
        theme(axis.title.y = element_text(face = "bold", size = 16)) +  # adjusts the title of y axis
        theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
        theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black")) + # adjusts text of y axis
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, color = "black")) + # adjusts text of x axis
        scale_y_continuous(expand = c(0, 0)) +
        expand_limits(y=c(0, 4)) +
        scale_x_discrete(limits=c("AMZ", "PAM", "BAR", "CEA")) +
        #annotate("text", x = 1, y = 1.5, label = "*", size = 10, fontface = "bold") +
        #annotate("text", x = 4, y = 3.5, label = "*", size = 10, fontface = "bold") +
        ylab("Number of taxa") + # add the title on y axis
        theme(axis.line = element_line(color = "black"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_blank(), panel.background = element_blank()) # removes the gridlines
meio_taxa_basin_plot

meio_density_basin_plot <- ggplot(meio_sum_density_basin, aes(x= Basin, y= Density10)) + 
        geom_bar(stat = "identity", width= 0.8, fill = "#1F78B4", color = "black", size = 0.5) +
        geom_errorbar(aes(ymin=Density10-se, ymax=Density10+se), width= 0.2, size =  0.5) +
        geom_point(size = 2, shape = 21, fill = "#D95F02") +
        ylab(expression(bold(paste(Density~(inds.10~cm^-2))))) + # format title on y axis
        theme(axis.title.y = element_text(face = "bold", size = 16)) +  # adjusts the title of y axis
        theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
        theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black")) + # adjusts text of y axis
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, color = "black")) + # adjusts text of x axis
        scale_y_continuous(expand = c(0, 0)) +
        expand_limits(y=c(0, 300)) +
        scale_x_discrete(limits=c("AMZ", "PAM", "BAR", "CEA")) +
        #annotate("text", x = 1, y = 240, label = "*", size = 10, fontface = "bold") +
        xlab("Basin") +
        #ylab("Density (inds.10 cm-2)") + # add the title on x axis
        theme(axis.line = element_line(color = "black"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_blank(), panel.background = element_blank()) # removes the gridlines
meio_density_basin_plot

Figure1_meio <- ggarrange(meio_taxa_basin_plot, meio_density_basin_plot, labels = c("A", "B"), nrow = 1, ncol = 2)
Figure1_meio
ggsave("Figure1_meio.tiff", width = 6, height = 4, dpi = 150) # save graphic

######################################

# Barplots of standard error of the mean for meiofauna and nematode together
nema_meio_density_basin_plot <- ggplot(nema_meio_sum_density_basin, aes(x= Basin, y= Density10, fill = Group)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin=Density10-se, ymax=Density10+se), width= 0.2, position = position_dodge(0.9)) +
  geom_point(aes(y=Density10), color = "#FF7F00", size = 2, position = position_dodge(0.9)) +
  scale_fill_manual(values=c("#1F78B4", "#33A02C")) +
  ylab(expression(bold(paste(Density~(inds.10~cm^-2))))) + # format title on y axis
  theme(axis.title.y = element_text(face = "bold", size = 16)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, color = "black")) + # adjusts text of x axis
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y=c(0, 250)) +
  scale_x_discrete(limits=c("AMZ", "PAM", "BAR", "CEA")) +
  #annotate("text", x = 1, y = 240, label = "*", size = 10, fontface = "bold") +
  xlab("Basin") +
  #ylab("Density (inds.10 cm-2)") + # add the title on x axis
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) # removes the gridlines
nema_meio_density_basin_plot

nema_meio_taxa_basin_plot <- ggplot(nema_meio_sum_taxa_basin, aes(x= Basin, y= Ntaxa, fill = Group)) + 
        geom_bar(stat = "identity", position = position_dodge()) +
        geom_errorbar(aes(ymin=Ntaxa-se, ymax=Ntaxa+se), width= 0.2, position = position_dodge(0.9)) +
        geom_point(aes(y=Ntaxa), color = "#FF7F00", size = 2, position = position_dodge(0.9)) +
        scale_fill_manual(values=c("#1F78B4", "#33A02C")) +
        ylab(expression(bold(paste(Density~(inds.10~cm^-2))))) + # format title on y axis
        theme(axis.title.y = element_text(face = "bold", size = 16)) +  # adjusts the title of y axis
        theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
        theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black")) + # adjusts text of y axis
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, color = "black")) + # adjusts text of x axis
        theme(legend.position = "none") +
        scale_y_continuous(expand = c(0, 0), breaks = seq(0, 16, by = 4)) +
        expand_limits(y=c(0, 16)) +
        scale_x_discrete(limits=c("AMZ", "PAM", "BAR", "CEA")) +
        #annotate("text", x = 1, y = 17, label = "*", size = 10, fontface = "bold") +
        xlab("Basin") +
        ylab("Number of Taxa") + # add the title on x axis
        theme(axis.line = element_line(color = "black"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_blank(), panel.background = element_blank()) # removes the gridlines
nema_meio_taxa_basin_plot

Figure1_meio_nema <- ggarrange(nema_meio_taxa_basin_plot, nema_meio_density_basin_plot, labels = c("A", "B"), nrow = 1, ncol = 2)
Figure1_meio_nema
ggsave("Figure1_meio_nema.tiff", width = 6, height = 4, dpi = 150) # save graphic

# Barplots of standard error of the mean for Nematoda
nema_gen_basin_plot <- ggplot(nema_sum_gen_basin, aes(x= Basin, y= Ngenera)) + 
        geom_bar(stat = "identity", width= 0.8, fill = "#33A02C", color = "black", size = 0.5) +
        geom_errorbar(aes(ymin=Ngenera-se, ymax=Ngenera+se), width= 0.2, size =  0.5) +
        geom_point(size = 2, shape = 21, fill = "#D95F02") +
        theme(axis.title.y = element_text(face = "bold", size = 16)) +  # adjusts the title of y axis
        #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
        theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black")) + # adjusts text of y axis
        #theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, color = "black")) + # adjusts text of x axis
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        scale_y_continuous(expand = c(0, 0)) +
        expand_limits(y=c(0, 20)) +
        scale_x_discrete(limits=c("AMZ", "PAM", "BAR", "CEA")) +
        #annotate("text", x = 1, y = 16, label = "*", size = 10, fontface = "bold") +
        ylab("Number of genera") + # add the title on y axis
        #xlab("Basin") +
        theme(axis.line = element_line(color = "black"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_blank(), panel.background = element_blank()) # removes the gridlines
nema_gen_basin_plot

nema_density_basin_plot <- ggplot(nema_sum_density_basin, aes(x= Basin, y= Density10)) + 
        geom_bar(stat = "identity", width= 0.8, fill = "#33A02C", color = "black", size = 0.5) +
        geom_errorbar(aes(ymin=Density10-se, ymax=Density10+se), width= 0.2, size =  0.5) +
        geom_point(size = 2, shape = 21, fill = "#D95F02") +
        ylab(expression(bold(paste(Density~(inds.10~cm^-2))))) + # format title on y axis
        theme(axis.title.y = element_text(face = "bold", size = 16)) +  # adjusts the title of y axis
        #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
        theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black")) + # adjusts text of y axis
        #theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, color = "black")) + # adjusts text of x axis
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        scale_y_continuous(expand = c(0, 0)) +
        expand_limits(y=c(0, 250)) +
        scale_x_discrete(limits=c("AMZ", "PAM", "BAR", "CEA")) +
        #annotate("text", x = 1, y = 240, label = "*", size = 10, fontface = "bold") +
        #xlab("Basin") +
        #ylab("Density (inds.10 cm-2)") + # add the title on x axis
        theme(axis.line = element_line(color = "black"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_blank(), panel.background = element_blank()) # removes the gridlines
nema_density_basin_plot

nema_es_basin_plot <- ggplot(nema_sum_es_basin, aes(x= Basin, y= ES51)) + 
        geom_bar(stat = "identity", width= 0.8, fill = "#33A02C", color = "black", size = 0.5) +
        geom_errorbar(aes(ymin=ES51-se, ymax=ES51+se), width= 0.2, size =  0.5) +
        geom_point(size = 2, shape = 21, fill = "#D95F02") +
        theme(axis.title.y = element_text(face = "bold", size = 16)) +  # adjusts the title of y axis
        theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
        theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black")) + # adjusts text of y axis
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, color = "black")) + # adjusts text of x axis
        scale_y_continuous(expand = c(0, 0)) +
        expand_limits(y=c(0, 20)) +
        scale_x_discrete(limits=c("AMZ", "PAM", "BAR", "CEA")) +
        #annotate("text", x = 1, y = 16, label = "*", size = 10, fontface = "bold") +
        xlab("Basin") +
        ylab("ES(51)") + # add the title on x axis
        theme(axis.line = element_line(color = "black"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_blank(), panel.background = element_blank()) # removes the gridlines
nema_es_basin_plot

nema_h_basin_plot <- ggplot(nema_sum_h_basin, aes(x= Basin, y= H)) + 
        geom_bar(stat = "identity", width= 0.8, fill = "#33A02C", color = "black", size = 0.5) +
        geom_errorbar(aes(ymin=H-se, ymax=H+se), width= 0.2, size =  0.5) +
        geom_point(size = 2, shape = 21, fill = "#D95F02") +
        theme(axis.title.y = element_text(face = "bold", size = 16)) +  # adjusts the title of y axis
        theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
        theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black")) + # adjusts text of y axis
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, color = "black")) + # adjusts text of x axis
        scale_y_continuous(expand = c(0, 0)) +
        expand_limits(y=c(0, 3)) +
        scale_x_discrete(limits=c("AMZ", "PAM", "BAR", "CEA")) +
        #annotate("text", x = 1, y = 2.5, label = "*", size = 10, fontface = "bold") +
        xlab("Basin") +
        ylab("H'") + # add the title on x axis
        theme(axis.line = element_line(color = "black"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_blank(), panel.background = element_blank()) # removes the gridlines
nema_h_basin_plot

nema_j_basin_plot <- ggplot(nema_sum_j_basin, aes(x= Basin, y= J)) + 
        geom_bar(stat = "identity", width= 0.8, fill = "#33A02C", color = "black", size = 0.5) +
        geom_errorbar(aes(ymin=J-se, ymax=J+se), width= 0.2, size =  0.5) +
        geom_point(size = 2, shape = 21, fill = "#D95F02") +
        theme(axis.title.y = element_text(face = "bold", size = 16)) +  # adjusts the title of y axis
        theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
        theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black")) + # adjusts text of y axis
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, color = "black")) + # adjusts text of x axis
        scale_y_continuous(expand = c(0, 0)) +
        expand_limits(y=c(0.25, 1.0)) +
        scale_x_discrete(limits=c("AMZ", "PAM", "BAR", "CEA")) +
        #annotate("text", x = 1, y = 0.97, label = "*", size = 10, fontface = "bold") +
        #annotate("text", x = 4, y = 0.6, label = "*", size = 10, fontface = "bold") +
        xlab("Basin") +
        ylab("J") + # add the title on x axis
        theme(axis.line = element_line(color = "black"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_blank(), panel.background = element_blank()) # removes the gridlines
nema_j_basin_plot

# Barplots of standard error of the mean for meiofauna
nema_feeding_basin_plot <- ggplot(nema_sum_feeding_basin, aes(x= Basin, y= Abundance, fill = FeedingGroup)) + 
        geom_bar(stat = "identity", width= 0.8, color = "black", size = 0.5) +
        #scale_fill_brewer(name = "Feeding Type", palette = "Paired") +
        scale_fill_manual(values = c("#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F")) +
        theme(axis.title.y = element_text(face = "bold", size = 16)) +  # adjusts the title of y axis
        #theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
        theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black")) + # adjusts text of y axis
        #theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, color = "black")) + # adjusts text of x axis
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_discrete(limits=c("AMZ", "PAM", "BAR", "CEA")) +
        ylab("Feeding Type (%)") + # add the title on y axis
        xlab("Basin") + # add the title on y axis
        theme(legend.position="top") +
        theme(legend.title = element_blank()) +
        theme(legend.title.align = 0.5) +
        theme(axis.line = element_line(color = "black"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_blank(), panel.background = element_blank()) # removes the gridlines
nema_feeding_basin_plot

nema_trophic_basin_plot <- ggplot(nema_suma_trophic_basin, aes(x= Basin, y= TropDiv)) +
        geom_bar(stat = "identity", width= 0.8, fill = "lightblue", color = "black", size = 0.5) +
        geom_errorbar(aes(ymin=TropDiv-se, ymax=TropDiv+se), width= 0.2, size =  0.5) +
        geom_point(size = 2, shape = 21, fill = "#D95F02") +
        theme(axis.title.y = element_text(face = "bold", size = 16)) +  # adjusts the title of y axis
        theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
        theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black")) + # adjusts text of y axis
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14, color = "black")) + # adjusts text of x axis
        scale_y_continuous(position =  "left", expand = c(0, 0)) +
        expand_limits(y=c(0, 0.8)) +
        scale_x_discrete(limits=c("AMZ", "PAM", "BAR", "CEA")) +
        #annotate("text", x = 1, y = 0.7, label = "*", size = 10, fontface = "bold") +
        #annotate("text", x = 4, y = 0.45, label = "*", size = 10, fontface = "bold") +
        ylab("Nematode ITD") + # add the title on y axis
        xlab("Basin") + # add the title on y axis
        theme(legend.position="right") +      
        theme(axis.line = element_line(color = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_rect(fill = "transparent", colour = NA)) # removes the gridlines
nema_trophic_basin_plot


Figure2 <- ggarrange(nema_gen_basin_plot, nema_density_basin_plot, nema_feeding_basin_plot,
                     nema_h_basin_plot, nema_j_basin_plot, nema_trophic_basin_plot,
                     labels = c("A", "B", "C", "D", "E", "F"), nrow = 2, ncol = 3,   common.legend = FALSE)
Figure2
ggsave("Figure2.tiff", width = 9, height = 7, dpi = 150) # save graphic

#################### nMDS
nema_nmds_basin$Basin <- factor(nema_nmds_basin$Basin, levels = c("AMZ", "PAM", "BAR", "CEA"))
nema_nmds_basin_plot <- ggplot(nema_nmds_basin) +
        geom_point(aes(x = NMDS1, y = NMDS2, shape=Basin, color=Basin), size=3)+
        #xlim(-2, 2) +
        #ylim(-1.5, 1.5) +
        theme(legend.position = "right") +
        theme(legend.title = element_text(face = "bold", size = 10))+
        theme(legend.title.align = 0.5) +
        annotate("text", x = -1.3, y = 1.5, label = "A", size = 4, fontface = "bold") +
        annotate("text", x = 1.5, y= 1.6, label = "Stress: 0.18", size = 3, fontface = "bold")+
        scale_shape_manual(values=c(15, 17, 18, 19))+
        scale_color_manual(values=c("#6A3D9A", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A"))+
        theme(axis.title.y = element_blank()) +  # adjusts the title of y axis
        theme(axis.title.x = element_blank()) + # adjusts the title of x axis
        theme(axis.text.y = element_blank()) + # adjusts text of y axis
        theme(axis.text.x = element_blank()) + # adjusts text of x axis
        theme(axis.ticks.x = element_blank()) +
        theme(axis.ticks.y = element_blank()) +
        #xlab("nMDS 1") +
        #ylab("nMDS 2") + # add the title on x axis
        coord_fixed(ratio = 1/1.2) +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),  #remove major-grid labels
              panel.grid.minor = element_blank(),  #remove minor-grid labels
              plot.background = element_blank())
nema_nmds_basin_plot
ggsave("Figure3A.tiff", width = 4, height = 2, dpi = 150) # save graphic

# Organize nMDS data to plot depth gradient with different colors (i.e., a gradient color per basin)
nema_nmds_depth
nema_nmds_depth$Depth <-as.factor(nema_nmds_depth$Depth)

nema_nmds_depth1<-as.data.frame(
  nema_nmds_depth %>% 
    group_by(NMDS1, NMDS2, Basin, Depth) %>%                     
    summarize(N= n())) %>%
    mutate(helper = as.character(group_indices(., Basin, Depth)))

# Define the color for control and treatment groups
AMZ <- c("#C6DBEF","#6BAED6", "#08519C")
CEA <- c("#C7E9C0", "#74C476", "#006D2C")
Palette_nmds <-c(AMZ, CEA)

# Specify color in ggplot using "geom_bar (fill=Palette)"
nema_nmds_depth_plot <- ggplot(nema_nmds_depth1, aes(NMDS1, NMDS2, shape = Basin, color = helper))+
  geom_point(size = 3) +
  scale_shape_manual(values=c(15, 19))+
  scale_color_manual(name = "Depth", values = Palette_nmds, labels = c("MB", "LB", "UA", "MB", "LB", "UA")) +
  guides(color=guide_legend(ncol=2)) +
  #xlim(-2, 2) +
  #ylim(-1.5, 1.5) +
  theme(legend.position = "right") +
  theme(legend.title = element_text(face = "bold", size = 10))+
  theme(legend.title.align = 0.5) +
  annotate("text", x = -2, y = 1.2, label = "B", size = 4, fontface = "bold")+
  annotate("text", x = 1.3, y= 1.2, label = "Stress: 0.23", size = 3, fontface = "bold")+
  theme(axis.title.y = element_blank()) +  # adjusts the title of y axis
  theme(axis.title.x = element_blank()) + # adjusts the title of x axis
  theme(axis.text.y = element_blank()) + # adjusts text of y axis
  theme(axis.text.x = element_blank()) + # adjusts text of x axis
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  #xlab("nMDS 1") +
  #ylab("nMDS 2") + # add the title on x axis
  coord_fixed(ratio = 1/1.2)+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
nema_nmds_depth_plot 
ggsave("Figure3B.tiff", width = 5, height = 2, dpi = 150) # save graphic

Figure3 <- ggarrange(nema_nmds_basin_plot, nema_nmds_depth_plot, nrow = 1, ncol = 2)
Figure3

#Non-parametric one-way ANOVA - Kruskal-Wallis rank sum test by Family/Feeding Group - ASV only
kt_diverse_shannon <- read.csv("diverse-rRNA-variants.csv", header = TRUE, sep=",")
kruskal.test(H ~ Family, data = kt_diverse_shannon) # KW-test by family
parwise_asv_shannon <- pairwise.wilcox.test(kt_diverse_shannon$H, kt_diverse_shannon$Family,
                 p.adjust.method = "BH")
write.table(parwise_asv_shannon[["p.value"]], file="parwise_asv_shannon.csv", sep=",")
kruskal.test(H ~ FG, data = kt_diverse_shannon) # KW-test by feeding group

kt_diverse_richness <- read.csv("diverse-rRNA-variants.csv", header = TRUE, sep=",")
kruskal.test(d ~ Family, data = kt_diverse_richness) # KW-test by family
parwise_asv_richness <- pairwise.wilcox.test(kt_diverse_richness$d, kt_diverse_richness$Family,
                 p.adjust.method = "BH")
write.table(parwise_asv_richness[["p.value"]], file="parwise_asv_richness.csv", sep=",")
kruskal.test(d ~ FG, data = kt_diverse_richness) # KW-test by feeding group

kt_diverse_simpson <- read.csv("diverse-rRNA-variants.csv", header = TRUE, sep=",")
kruskal.test(Lambda ~ Family, data = kt_diverse_simpson) # KW-test by family
parwise_asv_simpson <- pairwise.wilcox.test(kt_diverse_simpson$Lambda, kt_diverse_simpson$Family,
                 p.adjust.method = "BH")
write.table(parwise_asv_simpson[["p.value"]], file="parwise_asv_simpson.csv", sep=",")
kruskal.test(Lambda ~ FG, data = kt_diverse_simpson) # KW-test by feeding group

###################

p1 <- ggplot(morpho, aes(fill = Rank, x=Family, y=Count)) +
        geom_bar(position = "dodge", stat="identity") +
        scale_fill_manual(values= c("#00AFBB", "#E7B800", "grey50"))+
        theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
        theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
        theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 11)) + # adjusts text of y axis
        theme(axis.text.x = element_text(angle = 50, hjust = 1.01, vjust = 1.05, size = 11)) + # adjusts text of x axis
        ylab("Number of specimens") + # add the title on y axis
        xlab("Nematode Families") + # add the title on x axis
        #labs(fill = "Levels") + # change legend title
        theme(legend.position = c(0.88, 0.8)) + # add legend within the box with coordinates
        theme(legend.title=element_blank()) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
p1

p1 <- annotate_figure(p1,
                      fig.lab = "A",
                      fig.lab.pos = "top.left",
                      fig.lab.size = 14,
                      fig.lab.face = "bold")
p1			
ggsave("Figure_S1A.tiff", width = 6, height = 4, dpi = 300) # save graphic

p2 <- ggplot() +
        geom_bar(data = feeding, aes(x= TrophicGroup, y= Counts), stat="identity", fill = "#00AFBB", position ="dodge", width= 0.5) +
        geom_bar(data = feeding, aes(x= TrophicGroup, y= RA), stat="identity", fill = "#00AFBB", position ="dodge", width= 0.5) +
        scale_x_discrete(name = "Feeding Groups") +
        scale_y_continuous(name = "Number of specimens",
                           sec.axis = sec_axis(~ . *100/257, name = "Relative abundance",
                                               labels = function(b) {paste0(round(b), "%")})) +
        theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
        theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
        theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12)) + # adjusts text of y axis
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12)) + # adjusts text of x axis
        #ylab("Number of specimens") + # add the title on y axis
        #xlab("Feeding Groups") + # add the title on x axis
        #labs(fill = "Levels") + # change legend title
        theme(legend.position = c(0.85, 0.85)) + # add legend within the box with coordinates
        theme(legend.title=element_blank()) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
p2

