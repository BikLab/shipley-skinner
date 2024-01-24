# Set working directory
#setwd("/Users/tiagopereira/Dropbox/SK-data/18S-soil-single-nematode/R-codes")
#Install qiime2R
#download.file("https://github.com/jbisanz/qiime2R/archive/master.zip", "source.zip")
#unzip("source.zip")
#install.packages("qiime2R-master", repos = NULL, type="source")

# ggh4x
# install.packages("remotes")
#remotes::install_github("teunbrand/ggh4x")

#Load libraries
library("qiime2R")
library("ape")
library("Biostrings")
library("biomformat")
library("phyloseq")
library("Hmisc")
library("yaml")
library("tidyr")
library("stats")
library("utils")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("metagMisc")    # Miscellaneous functions for metagenomic analysis
library("vegan")
library("remotes")
library("mctoolsr")
library("RColorBrewer")
library("viridis")
library("decontam")     # Simple statistical identification of contaminated sequence features
library("plyr")
library("ALDEx2")
library("tibble")
library("gplots")
library("dplyr")        # filter and reformat data frames
library("ComplexHeatmap") # Package for Complex heatmaps reveal patterns and correlations in multidimensional genomic data
library("ggpubr")
library("ggh4x")
library("data.table")
library("scatterplot3d")
library("plotly.microbiome")
library("gg3D")
library("plot3D")
library("plotly")
library("circlize")
library("ggtreeExtra")
library("ggtree")
library("treeio")
library("tidytree")
library("ggstar")
library("ggnewscale")

# Set plotting theme
theme_set(theme_bw())

# Load sources
source("src/miseqR.R")

# Import asv table, taxonomy table, and sample table in excel format
asv_nema_tab_18S <- read_excel("data/SHSK_18S_Nematoda.xlsx", sheet = "ASV")
tax_nema_tab_18S <- read_excel("data/SHSK_18S_Nematoda.xlsx", sheet = "TAX")
sample_nema_tab_18S <- read_excel("data/SHSK_18S_Nematoda.xlsx", sheet = "METADATA")
tree_nema_18S <- read.tree("data/18S-soil-repseq-nematodes-only.tree")

# Convert asv and taxonomy tables to matrix format; sample table as dataframe
asv_nema_tab_18S <- as.matrix(asv_nema_tab_18S)
class(asv_nema_tab_18S)
tax_nema_tab_18S <- as.matrix(tax_nema_tab_18S)
sample_nema_tab_18S <- as.data.frame(sample_nema_tab_18S)

# Rename row names using ASV column for asv table and taxonomy table; rename row names using Sample column for sample table
rownames(asv_nema_tab_18S) <- asv_nema_tab_18S[,1]
rownames(tax_nema_tab_18S) <- tax_nema_tab_18S[,1]
rownames(sample_nema_tab_18S) <- sample_nema_tab_18S[,1]

# Exclude first column from all three tables; set asv table as numeric
asv_nema_tab_18S <- asv_nema_tab_18S[,-1]
class(asv_nema_tab_18S) <- "numeric"
tax_nema_tab_18S <- tax_nema_tab_18S[,-1]
sample_nema_tab_18S <- sample_nema_tab_18S[,-1]

#Transform to phyloseq objects
ASV_18S_nema = otu_table(asv_nema_tab_18S, taxa_are_rows = TRUE)
TAX_18S_nema = tax_table(tax_nema_tab_18S)
samples_18S_nema = sample_data(sample_nema_tab_18S)

# Creating a phyloseq object
phy_18S_nema <- phyloseq(ASV_18S_nema, TAX_18S_nema, samples_18S_nema, tree_nema_18S)
phy_18S_nema # 109 taxa and 45 samples
head(sample_data(phy_18S_nema))

#Visualize some data entries from phyloseq object
sample_names(phy_18S_nema)
rank_names(phy_18S_nema)
sample_variables(phy_18S_nema)

#################### preparing files for ggtree #################### 
# Transform number of reads to relative abundance
phy_18S_nema_comp <- microbiome::transform(phy_18S_nema, transform = "compositional")
phy_18S_nema_comp # 109 taxa and 45 samples

# Export phyloseq object as long format
# Relative abundance will be used for site annotations
phy_18S_nema_comp_df <- psmelt(phy_18S_nema_comp)
write.csv(phy_18S_nema_comp_df, "exported_tables/phy_18S_nema_comp_df.csv")

# Here we used the phyloseq with raw counts
# This will give us the highest abundance for each OTU and in each habitat
# We will plot it as a barplot
# Export phyloseq object as long format
phy_18S_nema_df <- psmelt(phy_18S_nema)

# Summarize data by OTU and Habitat
phy_18S_nema_df_sum <- phy_18S_nema_df %>%
  group_by(OTU, Habitat) %>%
  summarise(
    HigherAbundance = max(Abundance))

nema_sort_higher_abundance <- phy_18S_nema_df_sum[order(phy_18S_nema_df_sum$HigherAbundance,
                                                                                   decreasing = TRUE),]
nema_higher_abundance_unique <- distinct(nema_sort_higher_abundance, OTU, .keep_all = TRUE)

# Sum OTU abundances across all samples
# Save as dataframe and rename variables accordingly
otu_size_all_nema <- as.data.frame(taxa_sums(phy_18S_nema), decreasing = TRUE)
otu_size_all_nema <- rownames_to_column(otu_size_all_nema)
colnames(otu_size_all_nema) [c(1, 2)] <- c("OTU", "Size")

# Join otu_size_all_nema dataframe with melt comp phyloseq
# This will give us the total abundance for each OTU and we will plot as terminal node size
phy_18S_nema_comp_df <- left_join(phy_18S_nema_comp_df, otu_size_all_nema)
nema_sort <- phy_18S_nema_comp_df[order(phy_18S_nema_comp_df$Size,
                                        decreasing = TRUE),]
nema_unique <- distinct(nema_sort, OTU, .keep_all = TRUE)
write.csv(nema_unique, "exported_tables/nema_unique.csv")

nema_unique_high <- nema_unique[, -(2:17)]
nema_unique_high <- left_join(nema_unique_high, nema_higher_abundance_unique)
write.csv(nema_unique_high, "exported_tables/nema_unique_high.csv")
######################################## Plot annotated tree ######################################## 
# the abundance and types (i.e. feeding groups) of nematodes
dat1_nema <- read.csv("data/tree-example/nema_tippoint_attr.csv")

# the abundance of nematodes at different habitats
dat2_nema <- read.csv("data/tree-example/nema_ringheatmap_attr.csv")

# the abundance of nematodes at the habitat of greatest prevalence.
dat3_nema <- read.csv("data/tree-example/nema_barplot_attr.csv")

# Adjust the order of habitats
dat2_nema$Habitat <- factor(dat2_nema$Habitat,
                              levels=c("CHA", "CSS", "NGR",
                                       "HLC", "OWL", "RIP"))

dat3_nema$Habitat <- factor(dat3_nema$Habitat,
                            levels=c("CHA", "CSS", "NGR",
                                     "HLC", "OWL", "RIP"))

# The circular layout tree.
p_nematoda <- ggtree(tree_nema_18S, layout="fan", size=0.15, open.angle=5)

p_nematoda <- p_nematoda %<+% dat1_nema + geom_fruit(geom=geom_star,
                                                         mapping=aes(fill=Genus, starshape=Feeding.Group, size=Size),
                                                         position="identity",starstroke=0.1)+
  scale_fill_manual(values=c("#FFC125","#87CEFA","#7B68EE","#808080","#EE6A50",
                             "#9ACD32","#D15FEE","#FFC0CB","#800080","#8DEEEE",
                             "#006400","#800000","#f58231","#191970","#9a6324",
                             "#46f0f0","#1F78B4","#bcf60c","#e6beff","#4363d8",
                             "#008080","#e6194b","#f032e6","#aaffc3", "lightgreen",
                             "#ffe119","#000075","#fffac8","#9a6324","#FB9A99",
                             "#33A02C","#fabebe","#6A3D9A", "gray", "lightblue", "darkblue", "mediumspringgreen"),
                    guide=guide_legend(keywidth = 0.5, keyheight = 0.5, order=1,
                                       override.aes=list(starshape=15)),
                    na.translate=FALSE)+
  scale_starshape_manual(values=c(15, 1, 3, 5, 7, 2, 10, 9),
                         guide=guide_legend(keywidth = 0.5, keyheight = 0.5, order=2),
                         na.translate=FALSE)+
  scale_size_continuous(range = c(1, 2.5),
                        guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order=3,
                                             override.aes=list(starshape=15)))+
  new_scale_fill()+
  geom_fruit(data=dat2_nema, geom=geom_tile,
             mapping=aes(y=ID, x=Habitat, alpha=Abundance, fill=Habitat),
             color = "grey50", offset = 0.04,size = 0.02)+
  scale_alpha_continuous(range=c(0, 1),
                         guide=guide_legend(keywidth = 0.3, keyheight = 0.3, order=5)) +
  #cladelabels +
  geom_fruit(data=dat3_nema, geom=geom_bar,
             mapping=aes(y=ID, x=HigherAbundance, fill=Habitat),
             pwidth=0.38, orientation="y", stat="identity")+
  scale_fill_manual(values=c("#52392F", "#83643E","#C29B6C", 
                             "#397A4C", "#77C063", "#BEDB92"),
                    guide=guide_legend(keywidth = 0.3, keyheight = 0.3, order=4))+
  geom_treescale(fontsize=1.2, linesize=0.3, x=0.1, y=0.0001) +
  theme(legend.position=c(0.93, 0.5),
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=5),
        legend.text=element_text(size=4),
        legend.spacing.y = unit(0.02, "cm"))
p_nematoda
ggsave("data/tree-example/p_nematoda.pdf", width = 9, height = 8, dpi = 300)

######################################## 
# Check if there are samples with zero reads and remove them
smin <- min(sample_sums(phy_18S_nema)) #smin should be = 54, so need to remove samples

phy_18S_nema_df <- phy_18S_nema %>%
  #transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()
phy_18S_nema_df_agr = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Family, data=phy_18S_nema_df, FUN = "sum") # add common factors to use for plotting
unique(phy_18S_nema_df_agr$Family)
write.csv(phy_18S_nema_df_agr, "exported_tables/phy_18S_nema_df_agr.csv")

phy_18S_nema_df_agr_ra <- read.csv("exported_tables/phy_18S_nema_df_agr.csv")

############### Organize color scales and factors to plot top12 and top20 taxa as barplots
colors_top20 <- c("#9a6324", "#46f0f0", "#1F78B4", "#aaffc3", "#e6beff",  "#33A02C", "#4363d8", "#008080", "#FB9A99", "#e6194b",
                  "#f032e6", "#bcf60c", "#fabebe","#f58231", "#ffe119", "#6A3D9A", "#000075", "#fffac8", "#800000", "gray")


colors_top23 <- c("#9a6324", "#46f0f0", "#1F78B4", "#aaffc3", "#e6beff",  "#33A02C", "#4363d8", "#008080", "#FB9A99", "#e6194b",
                  "#f032e6", "#bcf60c", "#fabebe","#f58231", "#ffe119", "#6A3D9A", "#000075", "#fffac8", "#800000", "gray", "lightblue", "black", "darkgreen")

#Change facet habitat labels
soil.type.labs <- c("Compacted","Uncompacted")
names(soil.type.labs) <- c("Compact.Soil","Not.Compact.Soil")

habitat.type.labs <- c("CHA", "CSS", "HLC", "NGR", "OWL", "RIP")
names(habitat.type.labs) <- c("Chaparral", "Coastal.scrub", "Hollyleaf.cherry", "Native.grass", "Oak.wood", "Riparian")

# Plot bar chart 18S Soil samples per habitat - Nematoda only - family level - per site
phy_18S_nema_df_fam_site <- ggplot(phy_18S_nema_df_agr_ra, aes(x = Sample.Site, y = Abundance, fill = Family)) + #plotting by sample
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.type.labs, Habitat = habitat.type.labs)) +
  geom_bar(stat = "identity", width = 0.95, position = "fill") + # adds to 100%
  #geom_text(aes(label = ifelse(round(RA) >= 5, paste(round(RA, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
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
  ylab("Relative Abundance (%)") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
phy_18S_nema_df_fam_site
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_all_filt_top_genus_top20.pdf", width = 10, height = 5, dpi = 200) # save graphic