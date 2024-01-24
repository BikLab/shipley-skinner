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

# This script deals with treating/cleaning the raw data
# Creating ordinations (nMDS and PCA) based on different transformations
# Creating alpha-diversity metric plots and taxonomic barplots

# Import ASV table, taxonomy table, and sample table in excel format. It can also be imported as csv files, using read.csv if needed.
# Make sure the name of files match in case any updates were made.
# These files have both nematode and soil microbiome data. Later we will split it accordingly.
asv_tab_16S <- read_excel("data/ASV_table_Nov_06_2021.xlsx")
tax_tab_16S <- read_excel("data/Taxonomy_table_Jan_17_2022.xlsx")
sample_tab_16S <- read_excel("data/Sample_table_16S_Mar_31_2023.xlsx")

# Convert ASV and taxonomy tables to matrix format; sample table as dataframe
# These transformations are needed in order to create the phyloseq object later.
asv_tab_16S <- as.matrix(asv_tab_16S)
class(asv_tab_16S) # check the object class
tax_tab_16S <- as.matrix(tax_tab_16S)
sample_tab_16S <- as.data.frame(sample_tab_16S)
class(sample_tab_16S) # check the object class

# Rename row names using ASV column for ASV table and taxonomy table; rename using Sample column for sample table
# Basically our first column of each file will be used as rownames.
rownames(asv_tab_16S) <- asv_tab_16S[,1]
rownames(tax_tab_16S) <- tax_tab_16S[,1]
rownames(sample_tab_16S) <- sample_tab_16S[,1]

# Since we renamed the rows, now we don't need the first column any more, so the matrix has only numeric values (or only taxonomic strings for that matter).
# Exclude first column from all three tables; set ASV table as numeric
asv_tab_16S <- asv_tab_16S[,-1]
class(asv_tab_16S) <- "numeric"
tax_tab_16S <- tax_tab_16S[,-1]
sample_tab_16S <- sample_tab_16S[,-1]

# Transform each table/matrix to a phyloseq component
ASV_16S = otu_table(asv_tab_16S, taxa_are_rows = TRUE)
TAX_16S = tax_table(tax_tab_16S)
samples_16S = sample_data(sample_tab_16S)

# Create a phyloseq object by combine all three previous objects
# If we have a phylogenetic tree and representative seqs, we can also add into this phyloseq object.
phy_16S <- phyloseq(ASV_16S, TAX_16S, samples_16S)
phy_16S # 12818 taxa and 624 samples
head(sample_data(phy_16S)) # check the first 6 lines of your phyloseq object

#Visualize some data entries from the phyloseq object
sample_names(phy_16S) # list sample names
rank_names(phy_16S) # list the taxonomic ranks from your tax_table
sample_variables(phy_16S) # list the variables from your sample_data

# Check if there are samples with zero reads and remove them
smin <- min(sample_sums(phy_16S)) # smin = 0 means (see values on Global environment) you have samples with no reads. Originally sample number was 624.
phy_16S_prune <- prune_samples(sample_sums(phy_16S)>0, phy_16S)
smin <- min(sample_sums(phy_16S_prune)) # smin = 2
phy_16S_prune # 12818 taxa and 620 samples. We remove 4 samples from the original dataset with zero read number.

# Check for ASVs that have no counted reads.
any(taxa_sums(phy_16S_prune) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16S_prune) == 0) # gives the number of cases

# Inspect Library sizes. This gives you an idea of how many reads per sample/sample type.
df_phy_16S_prune <- as.data.frame(sample_data(phy_16S_prune)) # transform sample_data/mapping file from phyloseq object to a dataframe
df_phy_16S_prune$LibrarySize <- sample_sums(phy_16S_prune) # create another variable LibrarySize on the sample_data by summing the read count across all ASVs per sample
df_phy_16S_prune <- df_phy_16S_prune[order(df_phy_16S_prune$LibrarySize),] # reorder sample_data by library size
df_phy_16S_prune$Index <- seq(nrow(df_phy_16S_prune)) # create another variable Index by using the row number of the dataframe

# Plot library sizes using ggplot
# This plot shows the number of reads per sample type (e.g., real samples, controls, blanks, etc.)
# It also gives an idea of potential outliers in the dataset.
ggplot(data = df_phy_16S_prune, aes(x= Index, y= LibrarySize, color = Sample.Type)) +
  geom_point()
ggsave("results/all-dataset/16S_prune_library-size-sample-type-all-samples.pdf", width = 8, height = 6, dpi = 150)

####################################### Only 16S Soil samples #######################################
# Since the dataset has both, nematode and soil microbiomes, we will work with them separately.
# Here we subset our pruned phyloseq object to keep only samples from the soil microbiome dataset.
# Our mapping file has a variable called Experiment with two entries, "Microbiome" for 16S reads from nematodes and "Environment" for 16S reads from soil samples. 
phy_16S_prune_soil <- subset_samples(phy_16S_prune, Experiment == "Environment") # select only soil samples
phy_16S_prune_soil # 12818 taxa and 58 samples. This represents only the soil dataset and its respective controls/blanks

# We can check if we have only soil samples, by exploring the sample_data
unique(sample_data(phy_16S_prune_soil)$Experiment) # shows only Environment with represents soil samples

# Check/remove samples with zero reads if necessary.
# We already remove samples with 0 reads, so these commands have no effect
smin <- min(sample_sums(phy_16S_prune_soil)) #smin = 79
phy_16S_prune_soil <- prune_samples(sample_sums(phy_16S_prune_soil)>0, phy_16S_prune_soil)
smin <- min(sample_sums(phy_16S_prune_soil)) #smin = 79
phy_16S_prune_soil # 12818 taxa and 58 samples, same as previous object since there was no more samples with zero counts.

# Check for ASVs that have zero reads
# We do have many ASVs with count = 0 since we exclude many samples from the nematode microbiome dataset
any(taxa_sums(phy_16S_prune_soil) == 0) # if TRUE, there are ASVs with zero reads, otherwise FALSE
sum(taxa_sums(phy_16S_prune_soil) == 0) # gives the number of cases. Since we select only soil samples, there are lot of ASVs with zero count.

# So, here we prune ASVs with zero counts
phy_16S_prune_soil <-prune_taxa(taxa_sums(phy_16S_prune_soil) > 0, phy_16S_prune_soil) # prune taxa/ASVs with zero reads
phy_16S_prune_soil # 9923 taxa and 58 samples (before was 12818 taxa)
head(sample_data(phy_16S_prune_soil)) # list 6 first rows from sample_data
head(otu_table(phy_16S_prune_soil)) # list 6 first rows from otu_table

# Inspect Library sizes of soil samples only
df_phy_16S_prune_soil <- as.data.frame(sample_data(phy_16S_prune_soil))
df_phy_16S_prune_soil$LibrarySize <- sample_sums(phy_16S_prune_soil)
df_phy_16S_prune_soil <- df_phy_16S_prune_soil[order(df_phy_16S_prune_soil$LibrarySize),]
df_phy_16S_prune_soil$Index <- seq(nrow(df_phy_16S_prune_soil))

# Plot using ggplot for soil samples, all sample types
# Overall, soil blank samples also head a high number of reads, not so the negative control
# Positive controls had a high number of reads as expected
ggplot(data = df_phy_16S_prune_soil, aes(x= Index, y= LibrarySize, color = Sample.Type)) +
  geom_point()
ggsave("results/soil/16S_prune_soil_library-size-sample-type.pdf", width = 8, height = 6, dpi = 150)

# Plot using ggplot for soil samples, control vs. real samples
# Here we simplify by plotting as two categories. The high number of reads in the control samples are the Blank samples
ggplot(data = df_phy_16S_prune_soil, aes(x= Index, y= LibrarySize, color = Sample.Control)) +
  geom_point()
ggsave("results/soil/16S_prune_soil_library-size-sample-control.pdf", width = 8, height = 6, dpi = 150)

# Inspect control samples in the bi-dimensional space.
# NMDS ordination of all soil samples, including controls & blanks.
# This can give us an idea of potential contaminants may affect our dataset
phy_16S_prune_soil_clr <- microbiome::transform(phy_16S_prune_soil, "clr") # Apply centered log ration transformation due to the variability of the data
phy_16S_prune_soil_hel <- microbiome::transform(phy_16S_prune_soil, "hellinger") # Apply hellinger transformation, which standardize samples by totals and takes the sqrt.

# NMDS using clr transformation with all samples
set.seed(1)
phy_16S_prune_soil_clr_nmds <- ordinate(
  physeq = phy_16S_prune_soil_clr, 
  method = "NMDS", 
  distance = "euclidean" # use euclidean distance instead of bray-curtis since the data is clr transformed (i.e., negative values)
)

# NMDS using Hellinger transformation with all samples
set.seed(1)
phy_16S_prune_soil_hel_nmds <- ordinate(
  physeq = phy_16S_prune_soil_hel, 
  method = "NMDS", 
  distance = "bray" # use bray-curtis distance
)

# nMDS plot (clr transformation) with all samples including blanks and controls.
# Positive/negative controls fall in the midle of the plot with two low abundance RIP samples 
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16S_prune_soil_clr,
  ordination = phy_16S_prune_soil_clr_nmds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E","#C29B6C", 
                                "#397A4C", "#77C063", "#BEDB92", "blue", "orange", "purple"),
                     name = "Habitat",
                     breaks = c("Chaparral", "Coastal.scrub", "Native.grass", "Hollyleaf.cherry", "Oak.wood", "Riparian",
                                "NegCtrl", "Soil.Blank", "ZymoMock"),
                     labels = c("CHA", "CSS", "NGR", "HLC", "OWL", "RIP", "NC", "Blank", "PC"))+
  scale_shape_manual(values = c(19, 17, 3),
                     name = "Soil Type",
                     breaks = c("Compact.Soil", "Not.Compact.Soil", "no.data"),
                     labels = c("Compacted", "Uncompacted", "Controls"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  #ggrepel::geom_text_repel(aes(label = Sample.Site, color = Habitat), size = 2, max.overlaps = 30) +
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
  annotate("text", x = 55, y = 45, label ="2D Stress: 0.18")
ggsave("results/soil/phy_16S_prune_soil_clr_nmds_no-filtering.pdf", width = 7, height = 5, dpi = 200)

# nMDS plot (Hellinger transformation) with all samples including blanks and controls.
# Very divergent negative control sample obscure the grouping of samples
# Also, the stress value of the nMDS is extremely low, suggesting insufficient data.
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16S_prune_soil_hel,
  ordination = phy_16S_prune_soil_hel_nmds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E","#C29B6C", 
                                         "#397A4C", "#77C063", "#BEDB92", "blue", "orange", "purple"),
                                         name = "Habitat",
                     breaks = c("Chaparral", "Coastal.scrub", "Native.grass", "Hollyleaf.cherry", "Oak.wood", "Riparian",
                                "NegCtrl", "Soil.Blank", "ZymoMock"),
                     labels = c("CHA", "CSS", "NGR", "HLC", "OWL", "RIP", "NC", "Blank", "PC"))+
  scale_shape_manual(values = c(19, 17, 3),
                     name = "Soil Type",
                     breaks = c("Compact.Soil", "Not.Compact.Soil", "no.data"),
                     labels = c("Compacted", "Uncompacted", "Controls"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  #ggrepel::geom_text_repel(aes(label = Sample.Site, color = Habitat), size = 2, max.overlaps = 30) +
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
  annotate("text", x = 1000, y = 45, label ="2D Stress: 0.00008")
ggsave("results/soil/phy_16S_prune_soil_hel_nmds_no-filtering.pdf", width = 7, height = 5, dpi = 200)

# We can remove NC sample since it is an outlier and it is impacting the nMDS
# Then we can replot the data to see if we have better grouping (i.e., by habitat and soil type)
phy_16S_prune_soil_hel_nc <- subset_samples(phy_16S_prune_soil_hel, Habitat != "NegCtrl") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16S_prune_soil_hel_nc)$Habitat) # now the negative control sample was removed.

# We redo the nMDS with the Hellinger transformation but no NC sample
set.seed(1)
phy_16S_prune_soil_hel_nc_mds <- ordinate(
  physeq = phy_16S_prune_soil_hel_nc, 
  method = "NMDS", 
  distance = "bray" # use bray-curtis distance
) # Stress:     0.09102368

# nMDS plot (Hellinger transformation no NC sample) with all samples including blanks and controls.
# We can see now how the blanks and postive control fall in a different part of the plot.
# We also observe the grouping by soil type/habitats and how two low RIP samples are slightly off at the bottom of the plot.
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16S_prune_soil_hel_nc,
  ordination = phy_16S_prune_soil_hel_nc_mds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E","#C29B6C", 
                                         "#397A4C", "#77C063", "#BEDB92", "orange", "purple"),
                                         name = "Habitat",
                     breaks = c("Chaparral", "Coastal.scrub", "Native.grass", "Hollyleaf.cherry", "Oak.wood", "Riparian",
                                "Soil.Blank", "ZymoMock"),
                     labels = c("CHA", "CSS", "NGR", "HLC", "OWL", "RIP", "Blank", "PC"))+
  scale_shape_manual(values = c(19, 17, 3),
                     name = "Soil Type",
                     breaks = c("Compact.Soil", "Not.Compact.Soil", "no.data"),
                     labels = c("Compacted", "Uncompacted", "Controls"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  #ggrepel::geom_text_repel(aes(label = Sample.Site, color = Habitat), size = 2, max.overlaps = 30) +
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
  annotate("text", x = 2.5, y = 2.5, label ="2D Stress: 0.09")
ggsave("results/soil/phy_16S_prune_soil_hel_NC_nmds_no-filtering.pdf", width = 7, height = 5, dpi = 200)

# We no the blanks and controls are not impacting our data, but we still may have potential contaminants to remove.
# Here we will use R package decontam to identify potential contaminants
# Frequency method: it is based on the DNA concentration (i.e., Qubit value) of each sample in our dataset.
# The assumption is that control/blank samples will have low DNA concentration compared to real samples
phy_16S_prune_soil_contamdf_freq <- isContaminant(phy_16S_prune_soil, method="frequency", conc="QubitConc")
head(phy_16S_prune_soil_contamdf_freq) # list the first 6 ASVs
tail(phy_16S_prune_soil_contamdf_freq) # list the last 6 ASVs
table(phy_16S_prune_soil_contamdf_freq$contaminant) # gives the number of potential contaminants. 146 potential contaminants with this method
head(which(phy_16S_prune_soil_contamdf_freq$contaminant)) # list the first 6 contaminants
phy_16S_prune_soil_contamdf_freq_list <- rownames_to_column(phy_16S_prune_soil_contamdf_freq, var = "ASV") # create a dataframe with all AVSs and the category contaminant (false or true)
head(phy_16S_prune_soil_contamdf_freq_list)
write.csv(phy_16S_prune_soil_contamdf_freq_list, "exported_tables/soil/phy_16S_prune_soil_contamdf_freq_list.csv")

# Explore AVSs that are potentially contaminants, ASV abundance should be inversely correlated with DNA concentration
plot_frequency(phy_16S_prune_soil, taxa_names(phy_16S_prune_soil)[c(41, 62, 85, 122)], conc="QubitConc") + 
  xlab("DNA Concentration (ng/ul)")
ggsave("results/soil/16S_phy_prune_soil_contaminants_freq.pdf", width = 8, height = 6, dpi = 150)

# Explore AVSs that are potentially contaminants, randomly choosing 4 AVSs.
set.seed(100)
plot_frequency(phy_16S_prune_soil, taxa_names(phy_16S_prune_soil)[sample(which(phy_16S_prune_soil_contamdf_freq$contaminant),4)], conc="QubitConc") + 
  xlab("DNA Concentration (ng/ul)")
ggsave("results/soil/16S_phy_prune_soil_contaminants_freq_random.pdf", width = 8, height = 6, dpi = 150)

# Prevalence method: it is based on the frequency of ASVs showing up in the negative controls and blank samples vs. the real samples
sample_data(phy_16S_prune_soil)$is.neg <- sample_data(phy_16S_prune_soil)$Sample.Control == "Control.Sample"
phy_16S_prune_soil_contamdf_prev <- isContaminant(phy_16S_prune_soil, method="prevalence", neg="is.neg")
table(phy_16S_prune_soil_contamdf_prev$contaminant) # it identified 39 potential contaminants
head(which(phy_16S_prune_soil_contamdf_prev$contaminant))
phy_16S_prune_soil_contamdf_prev_list <- rownames_to_column(phy_16S_prune_soil_contamdf_prev, var = "ASV")
write.csv(phy_16S_prune_soil_contamdf_prev_list, "exported_tables/soil/phy_16S_prune_soil_contamdf_prev_list.csv")

# Prevalence method: here we increase threshhold value to 0.5
phy_16S_prune_soil_contamdf_prev05 <- isContaminant(phy_16S_prune_soil, method="prevalence", neg="is.neg", threshold=0.5)
table(phy_16S_prune_soil_contamdf_prev05$contaminant) # it identified 88 potential contaminants
phy_16S_prune_soil_contamdf_prev05_list <- rownames_to_column(phy_16S_prune_soil_contamdf_prev05, var = "ASV")
write.csv(phy_16S_prune_soil_contamdf_prev05_list, "exported_tables/soil/phy_16S_prune_soil_contamdf_prev05_list.csv")

# Make phyloseq object of presence-absence in negative controls and true samples
phy_16S_prune_soil_pa <- transform_sample_counts(phy_16S_prune_soil, function(abund) 1*(abund>0))
phy_16S_prune_soil_pa_neg <- prune_samples(sample_data(phy_16S_prune_soil_pa)$Sample.Control == "Control.Sample", phy_16S_prune_soil_pa)
phy_16S_prune_soil_pa_pos <- prune_samples(sample_data(phy_16S_prune_soil_pa)$Sample.Control == "True.Sample", phy_16S_prune_soil_pa)

# Make data.frame of prevalence in positive and negative samples, threshold of 0.1
df_phy_16S_prune_soil_pa <- data.frame(pa.pos=taxa_sums(phy_16S_prune_soil_pa), pa.neg=taxa_sums(phy_16S_prune_soil_pa_neg),
                              contaminant=phy_16S_prune_soil_contamdf_prev$contaminant)
ggplot(data=df_phy_16S_prune_soil_pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
ggsave("results/soil/phy_16S_prune_soil_contamdf_prev_01.pdf", width = 8, height = 6, dpi = 150)

# Make data.frame of prevalence in positive and negative samples, threshold of 0.5
df_phy_16S_prune_soil_pa05 <- data.frame(pa.pos=taxa_sums(phy_16S_prune_soil_pa), pa.neg=taxa_sums(phy_16S_prune_soil_pa_neg),
                                   contaminant=phy_16S_prune_soil_contamdf_prev05$contaminant)
ggplot(data=df_phy_16S_prune_soil_pa05, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
ggsave("results/soil/phy_16S_prune_soil_contamdf_prev_05.pdf", width = 8, height = 6, dpi = 150)

# Combined method: here we use both options from decontam (i.e., Frequency + Prevalence)
phy_16S_prune_soil_contamdf_comb <- isContaminant(phy_16S_prune_soil, conc="QubitConc", neg="is.neg", threshold=0.1, detailed = TRUE, normalize = TRUE, method="combined")
table(phy_16S_prune_soil_contamdf_comb$contaminant) # it identified 43 potential contaminants
head(which(phy_16S_prune_soil_contamdf_comb$contaminant))
phy_16S_prune_soil_contamdf_comb_list <- rownames_to_column(phy_16S_prune_soil_contamdf_comb, var = "ASV")
write.csv(phy_16S_prune_soil_contamdf_comb_list, "exported_tables/soil/phy_16S_prune_soil_contamdf_comb_list.csv")

# Combined method: here we increase the threshold to 0.5
phy_16S_prune_soil_contamdf_comb05 <- isContaminant(phy_16S_prune_soil, conc="QubitConc", neg="is.neg", threshold=0.5, detailed = TRUE, normalize = TRUE, method="combined")
table(phy_16S_prune_soil_contamdf_comb05$contaminant) # it identified 564 potential contaminants
head(which(phy_16S_prune_soil_contamdf_comb05$contaminant))
phy_16S_prune_soil_contamdf_comb05_list <- rownames_to_column(phy_16S_prune_soil_contamdf_comb05, var = "ASV")
write.csv(phy_16S_prune_soil_contamdf_comb05_list, "exported_tables/soil/phy_16S_prune_soil_contamdf_comb05_list.csv")

# Manual inspection, comparing blanks/controls and list of potential contaminants given by decontam, suggests prev05 is the most appropriated method
# Remove contaminants from prune phyloseq object using prev05.
phy_16S_prune_soil_noncontam_prev05 <- prune_taxa(!phy_16S_prune_soil_contamdf_prev05$contaminant, phy_16S_prune_soil)
phy_16S_prune_soil_noncontam_prev05 # 9835 taxa and 58 samples (originally was 9923 taxa)

# Check for samples with zero reads
smin <-min(sample_sums(phy_16S_prune_soil_noncontam_prev05)) # smin = 79
any(sample_sums(phy_16S_prune_soil_noncontam_prev05) == 0) # if TRUE, there are samples with no reads, otherwise FALSE
sum(sample_sums(phy_16S_prune_soil_noncontam_prev05) == 0) # gives the number of cases

# Check for ASVs that have no counted reads
any(taxa_sums(phy_16S_prune_soil_noncontam_prev05) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16S_prune_soil_noncontam_prev05) == 0) # gives the number of cases

# After removing potential contaminants (i.e. filtering out) with decontam, we can focus on the true samples only
# Subset phyloseq object by sample types, that is, keeping only soil samples and removing controls and blanks.
# Filter ASVs with a count = 0
phy_16S_prune_soil_noncontam_prev05_true <- subset_samples(phy_16S_prune_soil_noncontam_prev05, Sample.Type =="Soil") %>%
  prune_taxa(taxa_sums(.) > 0, .) # eliminates taxa with 0 abundance
phy_16S_prune_soil_noncontam_prev05_true # 9699 taxa and 52 samples

# Check if there is any AVS with count = 0
any(taxa_sums(phy_16S_prune_soil_noncontam_prev05_true) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16S_prune_soil_noncontam_prev05_true) == 0) # gives the number of cases
smin <-min(sample_sums(phy_16S_prune_soil_noncontam_prev05_true)) # smin = 79

# Checking the total number of reads and its distribution across samples after data cleaning
phy_16S_prune_soil_noncontam_prev05_true_readsums_df = data.frame(nreads = sort(taxa_sums(phy_16S_prune_soil_noncontam_prev05_true), TRUE),
                                                                  sorted = 1:ntaxa(phy_16S_prune_soil_noncontam_prev05_true), 
                                                                  type = "ASVs")
phy_16S_prune_soil_noncontam_prev05_true_readsums_df = rbind(phy_16S_prune_soil_noncontam_prev05_true_readsums_df, 
                                                             data.frame(nreads = sort(sample_sums(phy_16S_prune_soil_noncontam_prev05_true),TRUE), 
                                                                        sorted = 1:nsamples(phy_16S_prune_soil_noncontam_prev05_true), type = "Samples"))
title = "Total number of reads"
p = ggplot(phy_16S_prune_soil_noncontam_prev05_true_readsums_df, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
ggsave("results/soil/16S_prune_soil_noncontam_prev05_true-samples.pdf", width = 12, height = 6, dpi = 150)

# Until this point we just remove bad samples/AVSs
# Now we will filter out Chloroplasts and Mitochondria from the data matrix.
phy_16_prune_soil_noncontam_prev05_true_filt <- phy_16S_prune_soil_noncontam_prev05_true %>%
  subset_taxa(
    Family  != "Mitochondria" &
      Order   != "Chloroplast"
  )
phy_16_prune_soil_noncontam_prev05_true_filt # 9626 taxa and 52 samples

# Check if there is any AVS with count = 0
any(taxa_sums(phy_16_prune_soil_noncontam_prev05_true_filt) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16_prune_soil_noncontam_prev05_true_filt) == 0) # gives the number of cases
smin <-min(sample_sums(phy_16_prune_soil_noncontam_prev05_true_filt)) # smin = 79
phy_16_prune_soil_noncontam_prev05_true_filt # 9626 taxa and 52 samples

# Export phyloseq object with all true samples and also filtered mitochondria/chloroplasts
# This can be used as a starting point for further analysis
saveRDS(phy_16_prune_soil_noncontam_prev05_true_filt, "data/phy_16_prune_soil_noncontam_prev05_true_filt.RDS")

# Import phyloseq object
phy_16_prune_soil_noncontam_prev05_true_filt <- readRDS("data/phy_16_prune_soil_noncontam_prev05_true_filt.RDS")
class(phy_16_prune_soil_noncontam_prev05_true_filt) # check it is a phyloseq object
phy_16_prune_soil_noncontam_prev05_true_filt # 9626 taxa and 52 samples
smin <-min(sample_sums(phy_16_prune_soil_noncontam_prev05_true_filt)) # smin = 79

############################### Export and Transpose matrix for Picrust2 analysis ###############################
# Here will keep all samples and ASVs, including low abundance one
# ASV 16S soil table
asv_phy_16_prune_soil_noncontam_prev05_true_filt = as(otu_table(phy_16_prune_soil_noncontam_prev05_true_filt), "matrix")
write.csv(asv_phy_16_prune_soil_noncontam_prev05_true_filt, "picrust/asv_phy_16_prune_soil_noncontam_prev05_true_filt.csv")

# 16S soil sample table (metadata)
table_phy_16_prune_soil_noncontam_prev05_true_filt = as(sample_data(phy_16_prune_soil_noncontam_prev05_true_filt), "matrix")
class(table_phy_16_prune_soil_noncontam_prev05_true_filt)
table_phy_16_prune_soil_noncontam_prev05_true_filt_df <- as.data.frame(table_phy_16_prune_soil_noncontam_prev05_true_filt)
class(table_phy_16_prune_soil_noncontam_prev05_true_filt_df)
write.csv(table_phy_16_prune_soil_noncontam_prev05_true_filt_df, "picrust/table_phy_16_prune_soil_noncontam_prev05_true_filt_df.csv")

# 16S soil taxonomy
tax_phy_16_prune_soil_noncontam_prev05_true_filt = as(tax_table(phy_16_prune_soil_noncontam_prev05_true_filt), "matrix")
tax_phy_16_prune_soil_noncontam_prev05_true_filt_df <- as.data.frame(tax_phy_16_prune_soil_noncontam_prev05_true_filt)
write.csv(tax_phy_16_prune_soil_noncontam_prev05_true_filt_df, "picrust/tax_phy_16_prune_soil_noncontam_prev05_true_filt_df.csv")

############################### Apply different transformations to clean ASV table ############################### 
# Here we use the library microbiome to transform the data matrix
phy_16_prune_soil_noncontam_prev05_true_filt_comp <- microbiome::transform(phy_16_prune_soil_noncontam_prev05_true_filt, "compositional")
# phy_16_prune_soil_noncontam_prev05_true_filt_trans = transform_sample_counts(phy_16_prune_soil_noncontam_prev05_true_filt, function(x) 100 * x/sum(x)) # transform reads to relative abundance (0.0 - 1.0). This command is the equivalent to the microbiom::transform function "compositional"
phy_16_prune_soil_noncontam_prev05_true_filt_clr <- microbiome::transform(phy_16_prune_soil_noncontam_prev05_true_filt, "clr")
phy_16_prune_soil_noncontam_prev05_true_filt_log10 <- microbiome::transform(phy_16_prune_soil_noncontam_prev05_true_filt, "log10p")
phy_16_prune_soil_noncontam_prev05_true_filt_hel <- microbiome::transform(phy_16_prune_soil_noncontam_prev05_true_filt, "hellinger")
# phy_16_prune_soil_noncontam_prev05_true_filt_comp_sqrt <- transform_sample_counts(phy_16_prune_soil_noncontam_prev05_true_filt, function(x) sqrt(x/sum(x))) # same as Hellinger transformation

# NMDS ordination on the different transformed matrices
# All true samples; controls and blank samples and mitochondria/chloroplasts were removed
set.seed(1)
phy_16_prune_soil_noncontam_prev05_true_filt_comp_nmds <- ordinate(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_comp, 
  method = "NMDS", 
  distance = "bray"
) # Stress:     0.06663449

set.seed(1)
phy_16_prune_soil_noncontam_prev05_true_filt_clr_nmds <- ordinate(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_clr, 
  method = "NMDS", 
  distance = "euclidean"
) # Stress: 0.1584519

set.seed(1)
phy_16_prune_soil_noncontam_prev05_true_filt_log10_nmds <- ordinate(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_log10, 
  method = "NMDS", 
  distance = "bray"
) # Stress: 0.07083906

set.seed(1)
phy_16_prune_soil_noncontam_prev05_true_filt_hel_nmds <- ordinate(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_hel, 
  method = "NMDS", 
  distance = "bray"
) # Stress: 0.07480004

# nMDS plot compositional
# The plot looks nice, but the two low abundance RIP samples are still an issue
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_comp,
  ordination = phy_16_prune_soil_noncontam_prev05_true_filt_comp_nmds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E","#C29B6C", 
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
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  #ggrepel::geom_text_repel(aes(label = Sample.Site, color = Habitat), size = 2, max.overlaps = 30) +
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
  annotate("text", x = 3, y = 2, label ="2D Stress: 0.07")
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_comp_nmds.pdf", width = 6, height = 4, dpi = 200)

# nMDS plot clr
# Deals slightly better withe the low abundance RIP samples.
# One RIP sample is still more close to the compacted soil habitats
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_clr,
  ordination = phy_16_prune_soil_noncontam_prev05_true_filt_clr_nmds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E","#C29B6C", 
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
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  #ggrepel::geom_text_repel(aes(label = Sample.Site, color = Habitat), size = 2, max.overlaps = 30) +
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
  annotate("text", x = 45, y = 50, label ="2D Stress: 0.16")
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_clr_nmds.pdf", width = 6, height = 4, dpi = 200)

# nMDS plot log10
# Here RIP samples are too off the plot
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_log10,
  ordination = phy_16_prune_soil_noncontam_prev05_true_filt_log10_nmds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E","#C29B6C", 
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
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  #ggrepel::geom_text_repel(aes(label = Sample.Site, color = Habitat), size = 2, max.overlaps = 30) +
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
  annotate("text", x = -1, y = 2, label ="2D Stress: 0.07")
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_log10_nmds.pdf", width = 6, height = 4, dpi = 200)

# nMDS plot hellinger
# Similar to log10 nMDS plot
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_hel,
  ordination = phy_16_prune_soil_noncontam_prev05_true_filt_hel_nmds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E","#C29B6C", 
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
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  #ggrepel::geom_text_repel(aes(label = Sample.Site, color = Habitat), size = 2, max.overlaps = 30) +
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
  annotate("text", x = -1, y = 2, label ="2D Stress: 0.07")
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_hellinger_nmds.pdf", width = 6, height = 4, dpi = 200)

# Scale reads to even sequencing depth by removing samples with a total number of reads < 500.
smin <-min(sample_sums(phy_16_prune_soil_noncontam_prev05_true_filt)) # smin = 79
phy_16_prune_soil_noncontam_prev05_true_filt_scale <- prune_samples(sample_sums(phy_16_prune_soil_noncontam_prev05_true_filt)>500, phy_16_prune_soil_noncontam_prev05_true_filt)
smin <- min(sample_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale)) # smin = 4207
phy_16_prune_soil_noncontam_prev05_true_filt_scale # 9626 taxa and 50 samples (before we had 52 samples)

# Check for ASVs that have zero reads
any(taxa_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale) == 0) # gives the number of cases

# Remove ASVs that have zero reads and check again
phy_16_prune_soil_noncontam_prev05_true_filt_scale <- prune_taxa(taxa_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale) > 0, phy_16_prune_soil_noncontam_prev05_true_filt_scale)
any(taxa_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale) == 0) # gives the number of cases
phy_16_prune_soil_noncontam_prev05_true_filt_scale # 9622 taxa and 50 samples

############################### Export and Transpose final matrix for future analysis ###############################
# ASV 16S soil table
asv_phy_16_prune_soil_noncontam_prev05_true_filt_scale = as(otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale), "matrix")
write.csv(asv_phy_16_prune_soil_noncontam_prev05_true_filt_scale, "picrust/asv_phy_16_prune_soil_noncontam_prev05_true_filt_scale.csv")

# 16S soil sample table (metadata)
table_phy_16_prune_soil_noncontam_prev05_true_filt_scale = as(sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale), "matrix")
class(table_phy_16_prune_soil_noncontam_prev05_true_filt_scale)
table_phy_16_prune_soil_noncontam_prev05_true_filt_scale_df <- as.data.frame(table_phy_16_prune_soil_noncontam_prev05_true_filt_scale)
class(table_phy_16_prune_soil_noncontam_prev05_true_filt_scale_df)
write.csv(table_phy_16_prune_soil_noncontam_prev05_true_filt_scale_df, "picrust/table_phy_16_prune_soil_noncontam_prev05_true_filt_scale_df.csv")

# 16S soil taxonomy
tax_phy_16_prune_soil_noncontam_prev05_true_filt_scale = as(tax_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale), "matrix")
tax_phy_16_prune_soil_noncontam_prev05_true_filt_scale_df <- as.data.frame(tax_phy_16_prune_soil_noncontam_prev05_true_filt_scale)
write.csv(tax_phy_16_prune_soil_noncontam_prev05_true_filt_scale, "picrust/tax_phy_16_prune_soil_noncontam_prev05_true_filt_scale.csv")

###############################  Export scaled matrix and sample tables for analysis in Primer ############################### 
# Data transformation is done in Primer
asv_phy_16_prune_soil_noncontam_prev05_true_filt_scale_primer = as(otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale), "matrix")
asv_phy_16_prune_soil_noncontam_prev05_true_filt_scale_primer <- t(asv_phy_16_prune_soil_noncontam_prev05_true_filt_scale_primer)
write.csv(asv_phy_16_prune_soil_noncontam_prev05_true_filt_scale_primer, "exported_tables/soil/asv_phy_16_prune_soil_noncontam_prev05_true_filt_scale_primer.csv")

# Metada to be used in Primer
meta_phy_16_prune_soil_noncontam_prev05_true_filt_scale_primer <- sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale)
write.csv(meta_phy_16_prune_soil_noncontam_prev05_true_filt_scale_primer, "exported_tables/soil/meta_phy_16_prune_soil_noncontam_prev05_true_filt_scale_primer.csv")

# Taxonomy to be used in Primer
tax_phy_16_prune_soil_noncontam_prev05_true_filt_scale_primer <- tax_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale)
write.csv(tax_phy_16_prune_soil_noncontam_prev05_true_filt_scale_primer, "exported_tables/soil/tax_phy_16_prune_soil_noncontam_prev05_true_filt_scale_primer.csv")

############################### Apply different transformations to scaled ASV table ############################### 
# Samples with less than 500 reads were previously removed.
phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp <- microbiome::transform(phy_16_prune_soil_noncontam_prev05_true_filt_scale, "compositional")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr <- microbiome::transform(phy_16_prune_soil_noncontam_prev05_true_filt_scale, "clr")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_log10 <- microbiome::transform(phy_16_prune_soil_noncontam_prev05_true_filt_scale, "log10p")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_hel <- microbiome::transform(phy_16_prune_soil_noncontam_prev05_true_filt_scale, "hellinger")

############### Ordination using nMDS on scaled/transformed phyloseq objects #################
# Ordinate NMDS, all true samples > 500 reads
# Compositional transformation
set.seed(1)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_nmds <- ordinate(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp, 
  method = "NMDS", 
  distance = "bray",
  try = 50 #Stress:     0.07094927 
)

# Clr transformation
set.seed(1)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr_nmds <- ordinate(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr, 
  method = "NMDS", 
  distance = "euclidean",
  try = 50 #Stress:     0.1232263 
)

# Log10 transformation
set.seed(1)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_log10_nmds <- ordinate(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_log10, 
  method = "NMDS", 
  distance = "bray",
  try = 50 #Stress:     0.07958938
)

# Hellinger transformation
set.seed(1)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_hel_nmds <- ordinate(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_hel, 
  method = "NMDS", 
  distance = "bray",
  try = 50 #Stress:     0.07578406
)

# nMDS plot scaled comp transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp,
  ordination = phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_nmds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E","#C29B6C", 
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
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  #ggrepel::geom_text_repel(aes(label = Sample.Site, color = Habitat), size = 2, max.overlaps = 30) +
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
  annotate("text", x = 1.4, y = 0.8, label ="2D Stress: 0.07")
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_nmds.pdf", width = 6, height = 4, dpi = 200)

# nMDS plot scaled clr transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr,
  ordination = phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr_nmds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E","#C29B6C", 
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
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  #ggrepel::geom_text_repel(aes(label = Sample.Site, color = Habitat), size = 2, max.overlaps = 30) +
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
  annotate("text", x = 45, y = 55, label ="2D Stress: 0.12")
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr_nmds.pdf", width = 6, height = 4, dpi = 200)

# nMDS plot scaled log10 transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_log10,
  ordination = phy_16_prune_soil_noncontam_prev05_true_filt_scale_log10_nmds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E","#C29B6C", 
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
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  #ggrepel::geom_text_repel(aes(label = Sample.Site, color = Habitat), size = 2, max.overlaps = 30) +
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
  annotate("text", x = 1.4, y = 1.0, label ="2D Stress: 0.08")
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_log10_nmds.pdf", width = 6, height = 4, dpi = 200)

# nMDS plot scaled hellinger transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_hel,
  ordination = phy_16_prune_soil_noncontam_prev05_true_filt_scale_hel_nmds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E","#C29B6C", 
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
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  #ggrepel::geom_text_repel(aes(label = Sample.Site, color = Habitat), size = 2, max.overlaps = 30) +
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
  annotate("text", x = 1.4, y = 1.0, label ="2D Stress: 0.08")
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_hellinger_nmds.pdf", width = 6, height = 4, dpi = 200)

############### Ordination using PCA and clr transformed phyloseq objects
phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr_pca <- phyloseq::ordinate(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr,
  method = "RDA",
  distance = "euclidean"
)

# Plot scree plot to see PCA contributions
phyloseq::plot_scree(phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr_pca.pdf", width = 12, height = 6, dpi = 150)

head(phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr_pca$CA$eig)
sapply(phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr_pca$CA$eig[1:5], function(x) x / sum(phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr_pca$CA$eig))

# Scale axes and plot ordination
phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr_pca_clr1 <- phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr_pca$CA$eig[1] / sum(phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr_pca$CA$eig)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr_pca_clr2 <- phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr_pca$CA$eig[2] / sum(phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr_pca$CA$eig)

# Plot PCA based on clr transformation
phyloseq::plot_ordination(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr,
  ordination = phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr_pca,
  type="samples",
  color="Habitat",
  shape = "Soil.Type") + 
  scale_color_manual(values = c("#52392F", "#83643E","#C29B6C", 
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
  #annotate("text", x = 3, y = 2, label ="2D Stress: 0.07") +
  stat_ellipse(aes(group = Habitat), linetype = 2)
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_clr_pca.pdf", width = 6, height = 4, dpi = 200)

################################  SIMPER analysis using Vegan  ################################ 
# OTU table needs to be transposed for this analysis
# Simper between soil type
phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_soil_type <-
  vegan::simper(t(otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp)),
                sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp)$Soil.Type,
                permutations = 100, trace = FALSE)

# Save results showing cumulative contribution of each ASV and p-value
phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_soil_type_simper_df <- ldply(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_soil_type, data.frame)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_soil_type_simper_df, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_soil_type_simper_df.csv")

# Simper among habitats
phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_habitat <-
  vegan::simper(t(otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp)),
                sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp)$Habitat,
                permutations = 100, trace = FALSE)

# Save results showing cumulative contribution of each ASV and p-value
phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_habitat_simper_df <- ldply(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_habitat, data.frame)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_habitat_simper_df, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_habitat_simper_df.csv")

################################  Prune phyloseq object according to soil type  ################################ 
# Compacted soil
# Compositional transformation
phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_compact <- subset_samples(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp, Soil.Type == "Compact.Soil") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_compact # 6006 taxa and 27 samples

#Check taxa and samples with zeros
any(taxa_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_compact) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_compact) == 0) # gives the number of cases
any(sample_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_compact) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(sample_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_compact) == 0) # gives the number of cases

# Uncompacted soil
# Compositional transformation
phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_uncompact <- subset_samples(phy_16_prune_soil_noncontam_prev05_true_filt_scale, Soil.Type == "Not.Compact.Soil") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_uncompact # 4556 taxa and 23 samples

#Check taxa and samples with zeros
any(taxa_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_uncompact) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_uncompact) == 0) # gives the number of cases
any(sample_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_uncompact) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(sample_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_uncompact) == 0) # gives the number of cases

# Ordination NMDS compacted soil type
set.seed(1)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_compact_nmds <- ordinate(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_compact, 
  method = "NMDS", 
  distance = "bray",
  try = 50 #Stress:     0.1271002 
)

# Ordination NMDS uncompacted soil type
set.seed(1)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_uncompact_nmds <- ordinate(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_uncompact, 
  method = "NMDS", 
  distance = "bray",
  try = 50 #Stress:     0.1087629 
)

# Plot NMDS based on compositional transformation, only samples with >= 500 reads
# Compacted soil
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_compact,
  ordination = phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_compact_nmds,
  type = "samples",
  color = "Sample.Site",
  shape = "Habitat") +
  #facet_nested(~ Soil.Type+Habitat, scales = "free") +
  scale_color_manual(values = c("#52392F", "#52392F", "#52392F",
                                "#83643E", "#83643E", "#83643E",
                                "#C29B6C", "#C29B6C", "#C29B6C"),
                     name = "Site",
                     breaks = c("SK.01", "SK.02", "SK.03",
                                "SK.13", "SK.14", "SK.15",
                                "SK.16", "SK.17", "SK.18"),
                     labels = c("01", "02", "03",
                                "13", "14", "15",
                                "16", "17", "18"))+
  scale_shape_manual(values = c(19, 17, 8),
                     name = "Habitat",
                     breaks = c("Chaparral", "Coastal.scrub", "Native.grass"),
                     labels = c("CHA", "CSS", "NGR")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Sample.Site), alpha = 0.7, size = 4) +
  #geom_point(colour = "grey90", size = 1.5) +
  ggrepel::geom_text_repel(aes(label = Sample.Site), color = "black", size = 2, max.overlaps = 30) +
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
  annotate("text", x = 0.5, y = 0.8, label ="2D Stress: 0.13")
  ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_compact_nmds.pdf", width = 6, height = 4.5, dpi = 200)

# Uncompacted soil
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_uncompact,
  ordination = phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_uncompact_nmds,
  type = "samples",
  color = "Sample.Site",
  shape = "Habitat") +
  #facet_nested(~ Soil.Type+Habitat, scales = "free") +
  scale_color_manual(values = c("#397A4C", "#397A4C", "#397A4C", 
                                "#77C063", "#77C063", "#77C063",
                                "#BEDB92", "#BEDB92", "#BEDB92"),
                     name = "Site",
                     breaks = c("SK.04", "SK.05", "SK.06",
                                "SK.07", "SK.08", "SK.09",
                                "SK.10", "SK.11", "SK.12"),
                     labels = c("04", "05", "06",
                                "07", "08", "09",
                                "10", "11", "12"))+
  scale_shape_manual(values = c(15, 3, 18),
                     name = "Habitat",
                     breaks = c("Hollyleaf.cherry", "Oak.wood", "Riparian"),
                     labels = c("HLC", "OWL", "RIP"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Sample.Site), alpha = 0.7, size = 4) +
  #geom_point(colour = "grey90", size = 1.5) +
  ggrepel::geom_text_repel(aes(label = Sample.Site), color = "black", size = 2, max.overlaps = 30) +
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
  annotate("text", x = 0.8, y = 1.0, label ="2D Stress: 0.1")
  ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_comp_uncompact_nmds.pdf.pdf", width = 6, height = 4.5, dpi = 200)

# Save phyloseq object with samples >= 500 reads (i.e., scaled)
saveRDS(phy_16_prune_soil_noncontam_prev05_true_filt_scale, "data/phy_16_prune_soil_noncontam_prev05_true_filt_scale.RDS")

# Import phyloseq object
phy_16_prune_soil_noncontam_prev05_true_filt_scale <- readRDS("data/phy_16_prune_soil_noncontam_prev05_true_filt_scale.RDS")
class(phy_16_prune_soil_noncontam_prev05_true_filt_scale)
smin <-min(sample_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale)) # smin = 4207

######################### Generate a data.frame with alpha diversity measures #########################
# Calculate alpha-diversity measures (For plot purposes only!)
# This can be done using the different phyloseq alpha diversity measures
# You will get a Warning message for each index since there is no singletons on the dataset
alpha_div_16S_soil <- data.frame(
    "Observed" = phyloseq::estimate_richness(phy_16_prune_soil_noncontam_prev05_true_filt_scale, measures = "Observed"),
    "Shannon" = phyloseq::estimate_richness(phy_16_prune_soil_noncontam_prev05_true_filt_scale, measures = "Shannon"),
    "InvSimpson" = phyloseq::estimate_richness(phy_16_prune_soil_noncontam_prev05_true_filt_scale, measures = "InvSimpson"),
    "Reads" = phyloseq::sample_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale),
    "Soil.Type" = phyloseq::sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale)$Soil.Type,
    "Habitat" = phyloseq::sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale)$Habitat,
    "Sample.Site" = phyloseq::sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale)$Sample.Site)
alpha_div_16S_soil$Evenness <- alpha_div_16S_soil$Shannon/log(alpha_div_16S_soil$Observed)
head(alpha_div_16S_soil)

# Rename variable InvSimpson to Simpson
# The function rename & %>% works on dplyr. make sure it is loaded.
alpha_div_16S_soil <- alpha_div_16S_soil %>%
rename(Simpson = InvSimpson)
head(alpha_div_16S_soil)

# Reorder dataframe, first categorical then numerical variables  
alpha_div_16S_soil <- alpha_div_16S_soil[, c(7, 6, 5, 4, 1, 2, 3, 8)]
head(alpha_div_16S_soil)
write.csv(alpha_div_16S_soil, "exported_tables/soil/alpha_div_16S_soil.csv")

#Summarize alpha diversity measures by site
summary_alpha_16S_soil_site <- alpha_div_16S_soil %>%
    group_by(Soil.Type, Habitat, Sample.Site) %>%
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
write.csv(summary_alpha_16S_soil_site, "exported_tables/soil/summary_alpha_16S_soil_site.csv")

#Summarize alpha diversity measures by habitat
summary_alpha_16S_soil_habitat <- alpha_div_16S_soil %>%
  group_by(Soil.Type, Habitat) %>%
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
write.csv(summary_alpha_16S_soil_habitat, "exported_tables/soil/summary_alpha_16S_soil_habitat.csv")

#Summarize alpha diversity measures by soil type
summary_alpha_16S_soil_type <- alpha_div_16S_soil %>%
  group_by(Soil.Type) %>%
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
write.csv(summary_alpha_16S_soil_type, "exported_tables/soil/summary_alpha_16S_soil_type.csv")

# KW analysis alpha_16S_soil on all metrics at once
# Remember, numerical variables are from columns 4-8. The test is by soil type, column 3
kw_alpha_16S_soil_type_univ <- as.data.frame(sapply(4:8, function(x) kruskal.test(alpha_div_16S_soil[,x],
                                                                             alpha_div_16S_soil[,3])))

# Rename columns with the proper variable names
kw_alpha_16S_soil_type_univ <- kw_alpha_16S_soil_type_univ %>%
  rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)

# Make some adjustments on the table
kw_alpha_16S_soil_type_univ <- t(kw_alpha_16S_soil_type_univ) # transpose
kw_alpha_16S_soil_type_univ <- as_tibble(kw_alpha_16S_soil_type_univ, rownames = "Metric") # adding rownames as a column
kw_alpha_16S_soil_type_univ <- kw_alpha_16S_soil_type_univ[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_16S_soil_type_univ) # checking object class

# Save resulting table with fwrite to avoid any issues with characters
data.table::fwrite(kw_alpha_16S_soil_type_univ, "exported_tables/soil/kw_alpha_16S_soil_type_univ.csv")

# KW analysis alpha_16S_soil on all metrics at once
# Remember, numerical variables are from columns 4-8. The test is by habitat, column 2
kw_alpha_16S_habitat_univ <- as.data.frame(sapply(4:8, function(x) kruskal.test(alpha_div_16S_soil[,x],
                                                                                  alpha_div_16S_soil[,2])))

kw_alpha_16S_habitat_univ <- kw_alpha_16S_habitat_univ %>%
  rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)

# Make some adjustments on the table
kw_alpha_16S_habitat_univ <- t(kw_alpha_16S_habitat_univ) # transpose
kw_alpha_16S_habitat_univ <- as_tibble(kw_alpha_16S_habitat_univ, rownames = "Metric") # adding rownames as a column
kw_alpha_16S_habitat_univ <- kw_alpha_16S_habitat_univ[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_16S_habitat_univ) # checking object class

# Save resulting table with fwrite to avoid and issues with characters
data.table::fwrite(kw_alpha_16S_habitat_univ, "exported_tables/soil/kw_alpha_16S_habitat_univ.csv")

# pairwise comparison including p-value adjustment
kw_alpha_16S_soil_univ_pwt <- as.data.frame(sapply(4:8, function(x) pairwise.t.test(alpha_div_16S_soil[,x],
                                                                                    alpha_div_16S_soil[,2], p.adjust.method = "BH")))
kw_alpha_16S_soil_univ_pwt <- kw_alpha_16S_soil_univ_pwt %>%
  rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)

# Make some adjustments on the table
kw_alpha_16S_soil_univ_pwt <- t(kw_alpha_16S_soil_univ_pwt) # transpose
kw_alpha_16S_soil_univ_pwt <- as_tibble(kw_alpha_16S_soil_univ_pwt, rownames = "Metric") # adding rownames as a column
class(kw_alpha_16S_soil_univ_pwt) # checking object class

# Save resulting table with fwrite to avoid and issues with characters
fwrite(kw_alpha_16S_soil_univ_pwt, "exported_tables/soil/kw_alpha_16S_soil_univ_pwt.csv")

# pairwise comparison including p-value adjustment for each variable individually
# Across habitats, more than two grouping factors
# pwt_alpha_16S_habitat_obs <- pairwise.t.test(alpha_div_16S_soil$Observed, alpha_div_16S_soil$Habitat, p.adjust.method = "BH")
# pwt_alpha_16S_habitat_h <- pairwise.t.test(alpha_div_16S_soil$Shannon, alpha_div_16S_soil$Habitat, p.adjust.method = "BH")
# pwt_alpha_16S_habitat_d <- pairwise.t.test(alpha_div_16S_soil$Simpson, alpha_div_16S_soil$Habitat, p.adjust.method = "BH")
# pwt_alpha_16S_habitat_j <- pairwise.t.test(alpha_div_16S_soil$Evenness, alpha_div_16S_soil$Habitat, p.adjust.method = "BH")
# pwt_alpha_16S_habitat_n <- pairwise.t.test(alpha_div_16S_soil$Reads, alpha_div_16S_soil$Habitat, p.adjust.method = "BH")

# Between soil types
# kw_alpha_16S_soil_obs <- kruskal.test(Observed ~ Soil.Type, data = alpha_div_16S_soil)
# kw_alpha_16S_soil_h <- kruskal.test(Shannon ~ Soil.Type, data = alpha_div_16S_soil)
# kw_alpha_16S_soil_d <- kruskal.test(Simpson ~ Soil.Type, data = alpha_div_16S_soil)
# kw_alpha_16S_soil_j <- kruskal.test(Evenness ~ Soil.Type, data = alpha_div_16S_soil)
# kw_alpha_16S_soil_n <- kruskal.test(Reads ~ Soil.Type, data = alpha_div_16S_soil)

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
  
# Plot alpha-diversity per site within habitat/soil type
alpha_div_16S_soil %>%
    gather(key = metric, value = value, c("Reads", "Observed", "Shannon", "Simpson", "Evenness")) %>%
    mutate(metric = factor(metric, levels = c("Reads", "Observed", "Shannon", "Simpson", "Evenness"))) %>%
    ggplot(aes(x = Sample.Site, y = value)) +
    geom_boxplot(outlier.color = NA, width = 0.5) +
    geom_jitter(aes(color = Habitat), height = 0, width = .2) +
    facet_nested(metric ~ Soil.Type+Habitat, scales = "free", labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
    scale_color_manual(values = alpha_color_habitat) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
    theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
    theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
    scale_y_continuous(expand = c(0.1, 0.1)) + # plot as % and removes the internal margins
    scale_x_discrete(labels = site.labs,
    expand = c(0.1, 0.1),
    drop = FALSE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
    #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
    theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
    ylab("") + # add the title on y axis
    xlab("Site") + # add the title on x axis
    theme(legend.position="none")
ggsave("results/soil/16S_soil_alpha_diversity_phy_soil_16S_prune_noncontam_prev05_true_scale_site.pdf", width = 8, height = 6, dpi = 150) # save graphic

# Plot alpha-diversity per habitat within soil type
alpha_div_16S_soil %>%
  gather(key = metric, value = value, c("Reads", "Observed", "Shannon", "Simpson", "Evenness")) %>%
  mutate(metric = factor(metric, levels = c("Reads", "Observed", "Shannon", "Simpson", "Evenness"))) %>%
  ggplot(aes(x = Habitat, y = value)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  geom_jitter(aes(color = Habitat), height = 0, width = .2) +
  facet_nested(metric ~ Soil.Type, scales = "free", labeller = labeller(Soil.Type = soil.labs)) +
  scale_color_manual(values = alpha_color_habitat) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.1, 0.1)) + # plot as % and removes the internal margins
  scale_x_discrete(labels = habitat.labs,
                   expand = c(0.1, 0.1),
                   drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(legend.position="none")
ggsave("results/soil/16S_soil_alpha_diversity_phy_soil_16S_prune_noncontam_prev05_true_scale_habitat.pdf", width = 6, height = 6, dpi = 150) # save graphic

# Plot alpha-diversity per soil type
alpha_div_16S_soil$Soil.Type <- factor(alpha_div_16S_soil$Soil.Type,
                                         level = c("Compact.Soil", "Not.Compact.Soil"),
                                         labels = c("Compacted", "Uncompacted"))

alpha_color_soil_type <- c("#52392F", "#397A4C")

alpha_div_16S_soil %>%
    gather(key = metric, value = value, c("Reads", "Observed", "Shannon", "Simpson", "Evenness")) %>%
    mutate(metric = factor(metric, levels = c("Reads", "Observed", "Shannon", "Simpson", "Evenness"))) %>%
    ggplot(aes(x = Soil.Type, y = value)) +
    geom_boxplot(outlier.color = NA, width = 0.5) +
    geom_jitter(aes(color = Soil.Type), height = 0, width = .2) +
    facet_grid(metric ~ Soil.Type, scales = "free") +
    scale_color_manual(values = alpha_color_soil_type) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
    theme(axis.text.x = element_blank()) + # adjusts text of x axis
    theme(axis.ticks.x = element_blank())+
    theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
    theme(axis.title.x = element_blank()) + # adjusts the title of x axis
    scale_y_continuous(expand = c(0.1, 0.1)) + # plot as % and removes the internal margins
    #scale_x_discrete(
    # labels = c("CHA", "CSS", "NGR", "HLC",  "OWL", "RIP"),
    #expand = c(0.1, 0.1),
    #drop = FALSE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
    #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
    theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
    ylab("") + # add the title on y axis
    xlab("Soil Type") + # add the title on x axis
    theme(legend.position="none")
ggsave("results/soil/16S_soil_alpha_diversity_phy_soil_16S_prune_noncontam_prev05_true_scale_soil_type.pdf", width = 5, height = 6, dpi = 150) # save graphic

############################# Plotting top 20 most abundant at different taxonomic ranks ###############################
# This function can be used to split the data into the top N and the less abundant taxa
# Use scaled (samples >500 reads) phyloseq object
# Check how many phyla are in the phyloseq object
get_taxa_unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale, taxonomic.rank = "Phylum") # a total of 43 phyla
# Check how many taxonomic levels
head(tax_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale)) # 6 levels

# Function to collapse a certain number of taxa into category others
merge_less_than_top_prev05_top20 <- function(phy_16_prune_soil_noncontam_prev05_true_filt_scale, top=19){
    transformed <- transform_sample_counts(phy_16_prune_soil_noncontam_prev05_true_filt_scale, function(x) x/sum(x))
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
  
merge_less_than_top_prev05_top12 <- function(phy_16_prune_soil_noncontam_prev05_true_filt_scale, top=11){
    transformed <- microbiome::transform(phy_16_prune_soil_noncontam_prev05_true_filt_scale, "compositional")
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
  
# Another option is to collapse phyloseq by taxonomic rank
# phy_16_prune_soil_noncontam_prev05_true_filt_scale_phy <- phy_16_prune_soil_noncontam_prev05_true_filt_scale %>%
# tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
# transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance if not transformed already
# psmelt() %>%                                         # Melt to long format
# filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa if needed
# arrange(Phylum)                                      # Sort data frame alphabetically by phylum
# phy_16_prune_soil_noncontam_prev05_true_filt_scale_phy
# unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_phy$Habitat) # check variables
# sum(phy_16_prune_soil_noncontam_prev05_true_filt_scale_phy$Abundance) # sum values of variable abundance

# We will collapse the phyloseq object per taxonomic rank, from Phylum to Genus, one at time
# Phylum level top 12
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale, "Phylum")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12 <- merge_less_than_top_prev05_top12(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum, top=11)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12) # Melt to long format, it will give a warning messsage because of table headers
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Phylum, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Phylum, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_soil = aggregate(Abundance~Soil.Type+Phylum, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_site$Phylum)

head(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_site)
sum(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_site$Abundance)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_site,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_site.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_hab,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_hab.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_soil,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_soil.csv")

# Phylum level top 20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale, "Phylum")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20 <- merge_less_than_top_prev05_top20(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum, top=19)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20) # Melt to long format
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Phylum, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Phylum, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_soil = aggregate(Abundance~Soil.Type+Phylum, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_site$Phylum)

head(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_site)
sum(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_site$Abundance)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_site,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_site.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_hab,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_hab.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_soil,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_soil.csv")

# Class level top 12
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale, "Class")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12 <- merge_less_than_top_prev05_top12(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class, top=11)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12) # Melt to long format
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Class, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Class, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_soil = aggregate(Abundance~+Soil.Type+Class, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_site$Class)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_site,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_site.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_hab,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_hab.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_soil,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_soil.csv")

# Class level top 20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale, "Class")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20 <- merge_less_than_top_prev05_top20(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class, top=19)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20) # Melt to long format
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Class, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Class, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_agr_soil = aggregate(Abundance~Soil.Type+Class, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_agr_site$Class)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_agr_site,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_agr_site.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_agr_hab,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_agr_hab.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_agr_soil,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_agr_soil.csv")

# Order level top 12
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale, "Order")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12 <- merge_less_than_top_prev05_top12(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order, top=11)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12) # Melt to long format
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Order, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Order, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_agr_soil = aggregate(Abundance~Soil.Type+Order, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_agr_site$Order)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_agr_site,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_agr_site.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_agr_hab,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_agr_hab.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_agr_soil,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_agr_soil.csv")

# Order top 20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale, "Order")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20 <- merge_less_than_top_prev05_top20(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order, top=19)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20) # Melt to long format
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Order, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Order, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_soil = aggregate(Abundance~Soil.Type+Order, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_site$Order)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_site,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_site.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_hab,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_hab.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_soil,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_soil.csv")

# Family top 20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale, "Family")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20 <- merge_less_than_top_prev05_top20(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family, top=19)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20) # Melt to long format
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_agr_soil = aggregate(Abundance~Soil.Type+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_agr_site$Family)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_agr_site,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_agr_site.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_agr_hab,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_agr_hab.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_agr_soil,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_agr_soil.csv")

# Genus top 20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale, "Genus")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20 <- merge_less_than_top_prev05_top20(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus, top=19)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20) # Melt to long format
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Genus, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Genus, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_agr_soil = aggregate(Abundance~Soil.Type+Genus, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_agr_site$Genus)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_agr_site,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_agr_site.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_agr_hab,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_agr_hab.csv")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_agr_soil,
          "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_agr_soil.csv")

############### Organize color scales and factors to plot top12 and top20 taxa as barplots
# List of 20 distinct colors
colors_top20 <- c("#9a6324", "#46f0f0", "#1F78B4", "#aaffc3", "#000075",  "#33A02C", "#4363d8", "#008080", "#FB9A99", "#e6194b",
                         "#f032e6", "#bcf60c", "#fabebe","#f58231", "#ffe119", "#6A3D9A", "#e6beff", "#800000", "#fffac8", "gray")

# Other set of colors
# colors_top20 <- c("#A6CEE3", "#579CC7", "#3688AD", "#8BC395", "#89CB6C", "#40A635", "#919D5F", "#F99392", "#EB494A","#F79C5D",
#                   "#FDA746", "#FE8205", "#E39970", "#BFA5CF", "#8861AC", "#917099", "#E7E099", "#DEB969", "#B15928", "gray")

# List of 12 distinct colors
colors_top12 <- c("#9a6324", "#46f0f0", "#aaffc3", "#33A02C", "#4363d8", "#f032e6",
                  "#bcf60c", "#f58231", "#6A3D9A", "#e6beff", "#fffac8", "gray")


# Create new labels for important factors
soil.labs <- c("Compacted", "Uncompacted")
names(soil.labs) <- c("Compact.Soil", "Not.Compact.Soil")

habitat.labs <- c("CHA", "CSS", "NGR", "HLC",  "OWL", "RIP")
names(habitat.labs) <- c("Chaparral", "Coastal.scrub", "Native.grass",
                         "Hollyleaf.cherry", "Oak.wood",  "Riparian")
       
# Put "Others" to the final of the Phylum list - top 12
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_site$Phylum)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_site$Phylum <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_site$Phylum,
                                                                             levels = c("Acidobacteriota", "Actinobacteriota", "Bacteroidota",
                                                                                        "Chloroflexi", "Crenarchaeota", "Firmicutes",
                                                                                        "Gemmatimonadota", "Myxococcota", "Planctomycetota",
                                                                                        "Proteobacteria","Verrucomicrobiota", "Others"))
# Plot by site - Phylum level top 12
ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_site, aes(x = Sample.Site, y = Abundance, fill = Phylum)) +
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
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
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_site.pdf", width = 9, height = 4, dpi = 200) # save graphic

# Plot by habitat - Phylum level top 12
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_hab$Phylum <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_hab$Phylum,
                                                                                         levels = c("Acidobacteriota", "Actinobacteriota", "Bacteroidota",
                                                                                                    "Chloroflexi", "Crenarchaeota", "Firmicutes",
                                                                                                    "Gemmatimonadota", "Myxococcota", "Planctomycetota",
                                                                                                    "Proteobacteria","Verrucomicrobiota", "Others"))

ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_hab, aes(x = Habitat, y = Abundance, fill = Phylum)) +
  facet_grid(. ~ Soil.Type, scales = "free",
               labeller = labeller(Soil.Type = soil.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top12) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(labels = habitat.labs,
    expand = c(0.0, 0.0),
    drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_habitat.pdf", width = 6, height = 3.5, dpi = 300) # save graphic

# Plot by soil type - Phylum level top 12
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_soil$Phylum <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_soil$Phylum,
                                                                                             levels = c("Acidobacteriota", "Actinobacteriota", "Bacteroidota",
                                                                                                        "Chloroflexi", "Crenarchaeota", "Firmicutes",
                                                                                                        "Gemmatimonadota", "Myxococcota", "Planctomycetota",
                                                                                                        "Proteobacteria","Verrucomicrobiota", "Others"))



ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_agr_soil, aes(x = Soil.Type, y = Abundance, fill = Phylum)) + #plotting by sample
  facet_grid(. ~ Soil.Type, scales = "free",
             labeller = labeller(Soil.Type = soil.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top12) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_blank())+
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(expand = c(0.0, 0.0), drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Soil type") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top12_soil_type.pdf", width = 5, height = 3.5, dpi = 200) # save graphic

# Put "Others" to the final of the Phylum list - top 20
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_site$Phylum)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_site$Phylum <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_site$Phylum,
                                                                                       levels = c("Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota", "Bdellovibrionota",
                                                                                                  "Chloroflexi", "Crenarchaeota", "Cyanobacteria", "Deinococcota", "Dependentiae",
                                                                                                  "Firmicutes", "Gemmatimonadota", "Methylomirabilota", "Myxococcota", "Nitrospirota",
                                                                                                  "Planctomycetota", "Proteobacteria","RCP2-54", "Verrucomicrobiota", "Others"))
# Plot by site - Phylum level top 20
phylum_top20_site <- ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_site, aes(x = Sample.Site, y = Abundance, fill = Phylum)) +
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
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
phylum_top20_site
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_site.pdf", width = 10, height = 5, dpi = 200) # save graphic

# Plot by habitat - Phylum level top 20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_hab$Phylum <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_hab$Phylum,
                                                                                             levels = c("Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota", "Bdellovibrionota",
                                                                                                        "Chloroflexi", "Crenarchaeota", "Cyanobacteria", "Deinococcota", "Dependentiae",
                                                                                                        "Firmicutes", "Gemmatimonadota", "Methylomirabilota", "Myxococcota", "Nitrospirota",
                                                                                                        "Planctomycetota", "Proteobacteria","RCP2-54", "Verrucomicrobiota", "Others"))

ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_hab, aes(x = Habitat, y = Abundance, fill = Phylum)) + #plotting by sample
  facet_nested(. ~ Soil.Type, scales = "free",
               labeller = labeller(Soil.Type = soil.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(labels = habitat.labs,
    expand = c(0.0, 0.0),
    drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_habitat.pdf", width = 6.5, height = 4.5, dpi = 200) # save graphic

# Plot by soil type - Phylum level top 20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_soil$Phylum <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_soil$Phylum,
                                                                                             levels = c("Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota", "Bdellovibrionota",
                                                                                                        "Chloroflexi", "Crenarchaeota", "Cyanobacteria", "Deinococcota", "Dependentiae",
                                                                                                        "Firmicutes", "Gemmatimonadota", "Methylomirabilota", "Myxococcota", "Nitrospirota",
                                                                                                        "Planctomycetota", "Proteobacteria","RCP2-54", "Verrucomicrobiota", "Others"))



ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_agr_soil, aes(x = Soil.Type, y = Abundance, fill = Phylum)) + #plotting by sample
  facet_grid(. ~ Soil.Type, scales = "free",
             labeller = labeller(Soil.Type = soil.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_blank()) + 
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(labels = soil.labs, expand = c(0.0, 0.0), drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Soil type") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_phylum_top20_soil_type.pdf", width = 5, height = 5, dpi = 200) # save graphic

# Put "Others" to the final of the Class list - top 12
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_site$Class)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_site$Class <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_site$Class,
                                                                                     levels = c("Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidia",
                                                                                                "Blastocatellia", "Gammaproteobacteria", "Nitrososphaeria", "Planctomycetes",
                                                                                                "Thermoleophilia", "Verrucomicrobiae", "Vicinamibacteria", "Others"))
# Plot by site - Class level top12
ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_site, aes(x = Sample.Site, y = Abundance, fill = Class)) +
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
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
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_site.pdf", width = 10, height = 5, dpi = 200) # save graphic

# Plot by habitat - Class level top12
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_hab$Class <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_hab$Class,
                                                                                       levels = c("Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidia",
                                                                                                  "Blastocatellia", "Gammaproteobacteria", "Nitrososphaeria", "Planctomycetes",
                                                                                                  "Thermoleophilia", "Verrucomicrobiae", "Vicinamibacteria", "Others"))


ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_hab, aes(x = Habitat, y = Abundance, fill = Class)) + #plotting by sample
  facet_nested(. ~ Soil.Type, scales = "free",
             labeller = labeller(Soil.Type = soil.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top12) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(labels = habitat.labs,
                   expand = c(0.0, 0.0),
                   drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_habitat.pdf", width = 6, height = 4, dpi = 200) # save graphic

# Plot by soil type - Class level top12
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_soil$Class <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_soil$Class,
                                                                                           levels = c("Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidia",
                                                                                                      "Blastocatellia", "Gammaproteobacteria", "Nitrososphaeria", "Planctomycetes",
                                                                                                      "Thermoleophilia", "Verrucomicrobiae", "Vicinamibacteria", "Others"))





ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_agr_soil, aes(x = Soil.Type, y = Abundance, fill = Class)) + #plotting by sample
  facet_grid(. ~ Soil.Type, scales = "free",
             labeller = labeller(Soil.Type = soil.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top12) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_blank())+
  theme(axis.ticks.x = element_blank())+
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(labels = soil.labs, expand = c(0.0, 0.0), drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Soil type") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top12_soil_type.pdf", width = 5, height = 4, dpi = 200) # save graphic

# Put "Others" to the final of the Class list - top 20
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_agr_site$Class)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_agr_site$Class <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_agr_site$Class,
                                                                                       levels = c("Acidobacteriae", "Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidia",
                                                                                                  "Blastocatellia", "Clostridia", "Gammaproteobacteria", "Gemmatimonadetes", "Holophagae",
                                                                                                  "KD4-96", "Myxococcia", "Nitrososphaeria", "Planctomycetes", "Polyangia",
                                                                                                  "Rubrobacteria", "Thermoleophilia", "Verrucomicrobiae", "Vicinamibacteria", "Others"))
# Plot by site - Class level top20
class_top20 <- ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_agr_site, aes(x = Sample.Site, y = Abundance, fill = Class)) +
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
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
class_top20
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20.pdf", width = 10, height = 5, dpi = 200) # save graphic

# Plot by habitat - Class level top20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_agr_hab$Class <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_agr_hab$Class,
                                                                                       levels = c("Acidobacteriae", "Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidia",
                                                                                                  "Blastocatellia", "Clostridia", "Gammaproteobacteria", "Gemmatimonadetes", "Holophagae",
                                                                                                  "KD4-96", "Myxococcia", "Nitrososphaeria", "Planctomycetes", "Polyangia",
                                                                                                  "Rubrobacteria", "Thermoleophilia", "Verrucomicrobiae", "Vicinamibacteria", "Others"))

ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_agr_hab, aes(x = Habitat, y = Abundance, fill = Class)) + #plotting by sample
  facet_nested(. ~ Soil.Type, scales = "free",
               labeller = labeller(Soil.Type = soil.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(labels = habitat.labs,
    expand = c(0.0, 0.0),
    drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_class_top20_habitat.pdf", width = 7, height = 5, dpi = 300) # save graphic

# Put "Others" to the final of the Order list - top 12
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_agr_site$Order)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_agr_site$Order <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_agr_site$Order,
                                                                                     levels = c("Bacillales", "Burkholderiales", "Frankiales", "Gaiellales",
                                                                                                "Gemmatimonadales", "Micrococcales", "Nitrososphaerales", "Rhizobiales",
                                                                                                "Solirubrobacterales", "Sphingomonadales", "Vicinamibacterales", "Others"))
# Plot by site - Order level top12
ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12_agr_site, aes(x = Sample.Site, y = Abundance, fill = Order)) +
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
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
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top12.pdf", width = 10, height = 5, dpi = 200) # save graphic

# Put "Others" to the final of the order list - top 20
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_site$Order)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_site$Order <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_site$Order,
                                                                                     levels = c("Bacillales", "Burkholderiales", "Chitinophagales", "Chthoniobacterales",
                                                                                                "Corynebacteriales", "Cytophagales", "Elsterales", "Frankiales",
                                                                                                "Gaiellales", "Gemmatimonadales", "Micrococcales", "Micromonosporales",
                                                                                                "Nitrososphaerales", "Pseudonocardiales", "Pyrinomonadales", "Rhizobiales",
                                                                                                "Solirubrobacterales", "Sphingomonadales", "Vicinamibacterales", "Others"))
# Plot by site - Order level top20
order_top20 <- ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_site, aes(x = Sample.Site, y = Abundance, fill = Order)) +
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
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
order_top20
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_site.pdf", width = 10, height = 5, dpi = 200) # save graphic

# Plot by habitat - Order level top20
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_hab$Order)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_hab$Order <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_hab$Order,
                                                                                       levels = c("Bacillales", "Burkholderiales", "Chitinophagales", "Chthoniobacterales",
                                                                                                  "Corynebacteriales", "Cytophagales", "Elsterales", "Frankiales",
                                                                                                  "Gaiellales", "Gemmatimonadales", "Micrococcales", "Micromonosporales",
                                                                                                  "Nitrososphaerales", "Pseudonocardiales", "Pyrinomonadales", "Rhizobiales",
                                                                                                  "Solirubrobacterales", "Sphingomonadales", "Vicinamibacterales", "Others"))

order_top20_hab <- ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_hab, aes(x = Habitat, y = Abundance, fill = Order)) + #plotting by sample
  facet_nested(. ~ Soil.Type, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(labels = habitat.labs,
    expand = c(0.0, 0.0),
    drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
order_top20_hab
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_habitat.pdf", width = 7, height = 5, dpi = 200) # save graphic

# Plot by soil type - Order level top20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_soil$Order <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_soil$Order,
                                                                                           levels = c("Bacillales", "Burkholderiales", "Chitinophagales", "Chthoniobacterales",
                                                                                                      "Corynebacteriales", "Cytophagales", "Elsterales", "Frankiales",
                                                                                                      "Gaiellales", "Gemmatimonadales", "Micrococcales", "Micromonosporales",
                                                                                                      "Nitrososphaerales", "Pseudonocardiales", "Pyrinomonadales", "Rhizobiales",
                                                                                                      "Solirubrobacterales", "Sphingomonadales", "Vicinamibacterales", "Others"))

order_top20_soil <- ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_agr_soil, aes(x = Soil.Type, y = Abundance, fill = Order)) + #plotting by sample
facet_nested(. ~ Soil.Type, scales = "free",
              labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(expand = c(0.0, 0.0),
                   drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Soil Type") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
order_top20_soil
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_order_top20_soil.pdf", width = 7, height = 5, dpi = 200) # save graphic

# Put "Others" to the final of the family list - top 20
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_agr_site$Family)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_agr_site$Family <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_agr_site$Family,
                                                                                       levels = c("67-14", "Bacillaceae", "Beijerinckiaceae", "Chitinophagaceae",
                                                                                                  "Chthoniobacteraceae", "Comamonadaceae", "Gaiellales_fam", "Gemmatimonadaceae",
                                                                                                  "Geodermatophilaceae", "Micromonosporaceae", "Nitrososphaeraceae", "Oxalobacteraceae",
                                                                                                  "Pseudonocardiaceae", "Pyrinomonadaceae", "Solirubrobacteraceae", "Sphingomonadaceae",
                                                                                                  "Vicinamibacteraceae", "Vicinamibacterales_fam", "Xanthobacteraceae", "Others"))
# Plot by site - Family level top20
ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_agr_site, aes(x = Sample.Site, y = Abundance, fill = Family)) +
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
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
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_site.pdf", width = 10, height = 5, dpi = 200) # save graphic

# Plot by habitat - Family level top20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_agr_hab$Family <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_agr_hab$Family,
                                                                                         levels = c("67-14", "Bacillaceae", "Beijerinckiaceae", "Chitinophagaceae",
                                                                                                    "Chthoniobacteraceae", "Comamonadaceae", "Gaiellales_fam", "Gemmatimonadaceae",
                                                                                                    "Geodermatophilaceae", "Micromonosporaceae", "Nitrososphaeraceae", "Oxalobacteraceae",
                                                                                                    "Pseudonocardiaceae", "Pyrinomonadaceae", "Solirubrobacteraceae", "Sphingomonadaceae",
                                                                                                    "Vicinamibacteraceae", "Vicinamibacterales_fam", "Xanthobacteraceae", "Others"))

ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_agr_hab, aes(x = Habitat, y = Abundance, fill = Family)) + #plotting by sample
  facet_nested(. ~ Soil.Type, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(labels = habitat.labs,
    expand = c(0.0, 0.0),
    drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_family_top20_habitat.pdf", width = 7, height = 5, dpi = 200) # save graphic

# Put "Others" to the final of the genus list - top 20
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_agr_site$Genus)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_agr_site$Genus <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_agr_site$Genus,
                                                                                       levels = c("67-14_genus", "Bacillus", "Blastococcus", "Candidatus Udaeobacter",   
                                                                                                  "Elsterales_genus","Gaiellales_genus", "Gemmatimonadaceae_genus","Massilia",                
                                                                                                  "Microvirga","Modestobacter", "Mycobacterium","Nitrososphaeraceae_genus", 
                                                                                                  "RB41", "Rubrobacter","Solirubrobacter", "Sphingomonas",
                                                                                                  "Vicinamibacteraceae_genus", "Vicinamibacterales_genus","Xanthobacteraceae_genus","Others"))
# Plot by site - Genus level top20
ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_agr_site, aes(x = Sample.Site, y = Abundance, fill = Genus)) +
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
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
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_site.pdf", width = 10, height = 5, dpi = 200) # save graphic

# Plot by habitat - Genus level top20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_agr_hab$Genus <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_agr_hab$Genus,
                                                                                       levels = c("67-14_genus", "Bacillus", "Blastococcus", "Candidatus Udaeobacter",   
                                                                                                  "Elsterales_genus","Gaiellales_genus", "Gemmatimonadaceae_genus","Massilia",                
                                                                                                  "Microvirga","Modestobacter", "Mycobacterium","Nitrososphaeraceae_genus", 
                                                                                                  "RB41", "Rubrobacter","Solirubrobacter", "Sphingomonas",
                                                                                                  "Vicinamibacteraceae_genus", "Vicinamibacterales_genus","Xanthobacteraceae_genus","Others"))

ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_agr_hab, aes(x = Habitat, y = Abundance, fill = Genus)) + #plotting by sample
  facet_nested(. ~ Soil.Type, scales = "free",
               labeller = labeller(Soil.Type = soil.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(labels = habitat.labs,
    expand = c(0.0, 0.0),
    drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_top_genus_top20_habitat.pdf", width = 7, height = 5, dpi = 200) # save graphic


# Compacted soil only
# No transformation
phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact <- subset_samples(phy_16_prune_soil_noncontam_prev05_true_filt_scale, Soil.Type == "Compact.Soil") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact # 6006 taxa and 27 samples
phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_asv <- otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact, "matrix")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_asv, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_asv.csv")


# Function to collapse a certain number of taxa into category others
merge_less_than_top_prev05_top20_compact <- function(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact, top=19){
  transformed <- transform_sample_counts(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact, function(x) x/sum(x))
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


# Family level top 20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact, "Family")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam_top20 <- merge_less_than_top_prev05_top20_compact(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam, top=19)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam_top20_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam_top20) # Melt to long format
phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam_top20_agr_soil = aggregate(Abundance~Soil.Type+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam_top20_agr_site$Family)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam_top20_agr_site, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam_top20_agr_site.csv")

# Genus Level top 20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact, "Genus")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus_top20 <- merge_less_than_top_prev05_top20_compact(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus, top=19)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus_top20_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus_top20) # Melt to long format
phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Genus, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Genus, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus_top20_agr_soil = aggregate(Abundance~Soil.Type+Genus, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus_top20_agr_site$Genus)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus_top20_agr_site, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus_top20_agr_site.csv")

# Plot by site - Genus level top20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam_top20_agr_site$Family <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam_top20_agr_site$Family,
                                                                                           levels = c("67-14", "Bacillaceae", "Beijerinckiaceae", "Chitinophagaceae",
                                                                                                      "Chthoniobacteraceae", "Comamonadaceae", "Elsterales_fam", "Gaiellales_fam", "Gemmatimonadaceae",
                                                                                                      "Geodermatophilaceae", "Kineosporiaceae", "Micromonosporaceae", "Nitrososphaeraceae",
                                                                                                      "Pyrinomonadaceae", "Rubrobacteriaceae", "Solirubrobacteraceae", "Sphingomonadaceae",
                                                                                                      "Vicinamibacterales_fam", "Xanthobacteraceae", "Others"))

ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam_top20_agr_site, aes(x = Sample.Site, y = Abundance, fill = Family)) + #plotting by sample
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(expand = c(0.0, 0.0),
                   drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_fam_top20_agr_site.pdf", width = 7, height = 5, dpi = 200) # save graphic

# Plot by site - Family level top20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus_top20_agr_site$Genus <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus_top20_agr_site$Genus,
                                                                                               levels = c("67-14_genus", "Bacillus", "Beijerinckiaceae_genus", "Blastococcus",
                                                                                                          "Candidatus Udaeobacter", "Conexibacter", "Elsterales_genus", "Gaiellales_genus",
                                                                                                          "Gemmatimonadaceae_genus", "Geodermatophilus", "Microvirga", "Modestobacter",
                                                                                                          "Nitrososphaeraceae_genus", "RB41", "Rubrobacter", "Solirubrobacter", "Sphingomonas",
                                                                                                          "Vicinamibacteraceae_genus", "Vicinamibacterales_genus", "Others" ))

ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus_top20_agr_site, aes(x = Sample.Site, y = Abundance, fill = Genus)) + #plotting by sample
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(expand = c(0.0, 0.0),
                   drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_compact_genus_top20_agr_site.pdf", width = 7, height = 5, dpi = 200) # save graphic

# Uncompacted soil only
# No transformation
phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact <- subset_samples(phy_16_prune_soil_noncontam_prev05_true_filt_scale, Soil.Type == "Not.Compact.Soil") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact # 4556 taxa and 23 samples
phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_asv <- otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact, "matrix")
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_asv, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_asv.csv")

# Function to collapse a certain number of taxa into category others
merge_less_than_top_prev05_top20_uncompact <- function(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact, top=19){
  transformed <- transform_sample_counts(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact, function(x) x/sum(x))
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


# Family level top 20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact, "Family")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam_top20 <- merge_less_than_top_prev05_top20_uncompact(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam, top=19)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam_top20_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam_top20) # Melt to long format
phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam_top20_agr_soil = aggregate(Abundance~Soil.Type+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam_top20_agr_site$Family)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam_top20_agr_site, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam_top20_agr_site.csv")

# Genus Level top 20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact, "Genus")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus_top20 <- merge_less_than_top_prev05_top20_uncompact(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus, top=19)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus_top20_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus_top20) # Melt to long format
phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Genus, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Genus, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus_top20_agr_soil = aggregate(Abundance~Soil.Type+Genus, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus_top20_agr_site$Genus)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus_top20_agr_site, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus_top20_agr_site.csv")

# Plot by site - Family level top20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam_top20_agr_site$Family <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam_top20_agr_site$Family,
                                                                                               levels = c("67-14", "Bacillaceae", "Beijerinckiaceae", "Chitinophagaceae",
                                                                                                          "Chthoniobacteraceae", "Comamonadaceae", "Gaiellales_fam", "Gemmatimonadaceae",
                                                                                                          "Geodermatophilaceae", "Micrococcaceae", "Nitrosomonadaceae", "Nitrososphaeraceae",
                                                                                                          "Oxalobacteraceae", "Planococcaceae", "Solirubrobacteraceae", "Sphingomonadaceae",
                                                                                                          "Vicinamibacteraceae", "Vicinamibacterales_fam", "Xanthobacteraceae", "Others"))

ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam_top20_agr_site, aes(x = Sample.Site, y = Abundance, fill = Family)) + #plotting by sample
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(expand = c(0.0, 0.0),
                   drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_fam_top20_agr_site.pdf", width = 7, height = 5, dpi = 200) # save graphic

# Plot by site - Genus level top20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus_top20_agr_site$Genus <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus_top20_agr_site$Genus,
                                                                                                 levels = c("67-14_genus", "Bacillus", "Blastococcus", "Candidatus Udaeobacter",
                                                                                                            "Gaiellales_genus", "Gemmatimonadaceae_genus", "KD4-96_genus", "Massilia",
                                                                                                            "Microvirga", "Mycobacterium", "Nitrososphaeraceae_genus", "Planococcus",
                                                                                                            "RB41", "Solirubrobacter", "Sphingomonadaceae_genus", "Sphingomonas",
                                                                                                            "Vicinamibacteraceae_genus", "Vicinamibacterales_genus", "Xanthobacteraceae_genus", "Others" ))

ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus_top20_agr_site, aes(x = Sample.Site, y = Abundance, fill = Genus)) + #plotting by sample
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(expand = c(0.0, 0.0),
                   drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_uncompact_genus_top20_agr_site.pdf", width = 7, height = 5, dpi = 200) # save graphic

################################# Analysis/plots on selected bacteria genera ################################# 
phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_bacillus <- phy_16_prune_soil_noncontam_prev05_true_filt_scale %>%
  subset_taxa(Genus  == "Bacillus")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_bacillus # 49 taxa and 50 samples
head(otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_bacillus))
phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_bacillus_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_bacillus) # Melt to long format
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_bacillus_df, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_bacillus_df.csv")
bacillus_new <- read.csv("exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_bacillus_df_new.csv")

phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_blastococcus <- phy_16_prune_soil_noncontam_prev05_true_filt_scale %>%
  subset_taxa(Genus  == "Blastococcus")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_blastococcus # 18 taxa and 50 samples
head(otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_blastococcus))
phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_blastococcus_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_blastococcus) # Melt to long format
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_blastococcus_df, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_blastococcus_df.csv")
blastococcus_new <- read.csv("exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_blastococcus_df_new.csv")

phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_sphingomonas <- phy_16_prune_soil_noncontam_prev05_true_filt_scale %>%
  subset_taxa(Genus  == "Sphingomonas")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_sphingomonas # 109 taxa and 50 samples
head(otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_sphingomonas))
phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_sphingomonas_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_sphingomonas) # Melt to long format
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_sphingomonas_df, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_sphingomonas_df.csv")
sphingomonas_new <- read.csv("exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_sphingomonas_df_new.csv")

phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_sphingomonadaceae_unknown <- phy_16_prune_soil_noncontam_prev05_true_filt_scale %>%
  subset_taxa(Genus  == "Sphingomonadaceae_genus")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_sphingomonadaceae_unknown # 19 taxa and 50 samples
head(otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_sphingomonadaceae_unknown))
phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_sphingomonadaceae_unknown_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_sphingomonadaceae_unknown) # Melt to long format
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_sphingomonadaceae_unknown_df, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_sphingomonadaceae_unknown_df.csv")
sphingomonadaceae_unknown_new <- read.csv("exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_sphingomonadaceae_unknown_df_new.csv")


phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_nitrososphaeraceae_unknown <- phy_16_prune_soil_noncontam_prev05_true_filt_scale %>%
  subset_taxa(Genus  == "Nitrososphaeraceae_genus")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_nitrososphaeraceae_unknown # 35 taxa and 50 samples
head(otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_nitrososphaeraceae_unknown))
phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_nitrososphaeraceae_unknown_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_nitrososphaeraceae_unknown) # Melt to long format
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_nitrososphaeraceae_unknown_df, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_nitrososphaeraceae_unknown_df.csv")
nitrososphaeraceae_unknown_new <- read.csv("exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_genus_nitrososphaeraceae_unknown_df_new.csv")


# Create new labels for important factors
soil.labs <- c("Compacted", "Uncompacted")
names(soil.labs) <- c("Compact.Soil", "Not.Compact.Soil")

habitat.labs <- c("CHA", "CSS", "NGR", "HLC",  "OWL", "RIP")
names(habitat.labs) <- c("Chaparral", "Coastal.scrub", "Native.grass",
                         "Hollyleaf.cherry", "Oak.wood",  "Riparian")

# Bacillus
ggplot(bacillus_new, aes(x = Sample.Site, y = Abundance, fill = OTU)) + #plotting by sample
  facet_nested(OTU ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  #geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  #geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
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
  ylab("Abundance") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/bacillus_site.pdf", width = 7, height = 4, dpi = 200) # save graphic

# Blastococcus
ggplot(blastococcus_new, aes(x = Sample.Site, y = Abundance, fill = OTU)) + #plotting by sample
  facet_nested(OTU ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  #geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  #geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
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
  ylab("Abundance") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/blastococcus_site.pdf", width = 7, height = 5, dpi = 200) # save graphic

# Sphingomonas
ggplot(sphingomonas_new, aes(x = Sample.Site, y = Abundance, fill = OTU)) + #plotting by sample
  facet_nested(OTU ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  #geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  #geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
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
  ylab("Abundance ") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/sphingomonas_site.pdf", width = 7, height = 3, dpi = 200) # save graphic

# Sphingomonadaceae_unknown
ggplot(sphingomonadaceae_unknown_new, aes(x = Sample.Site, y = Abundance, fill = OTU)) + #plotting by sample
  facet_nested(OTU ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  #geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  #geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
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
  ylab("Abundance") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/sphingomonadaceae_unknown_site.pdf", width = 7, height = 3, dpi = 200) # save graphic

# Nitrososphaeraceae_unknown
ggplot(nitrososphaeraceae_unknown_new, aes(x = Sample.Site, y = Abundance, fill = OTU)) + #plotting by sample
  facet_nested(OTU ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  #geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  #geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
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
  ylab("Abundance") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/nitrososphaeraceae_unknown_site.pdf", width = 7, height = 3, dpi = 200) # save graphic

# selected bacterial genera all together
selected_genera_16s <- read_excel("exported_tables/soil/16S_selected_genera.xlsx")
ggplot(selected_genera_16s, aes(x = Sample.Site, y = Abundance, fill = OTU)) + #plotting by sample
  facet_nested(Genus+OTU ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  #geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  #geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20) +
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
  ylab("Abundance") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/selected_genera_16s_site.pdf", width = 7, height = 13, dpi = 200) # save graphic

########################### Only Major Bacteria Phyla ########################### 
phy_16_prune_soil_noncontam_prev05_true_filt_scale <- readRDS("data/phy_16_prune_soil_noncontam_prev05_true_filt_scale.RDS")
class(phy_16_prune_soil_noncontam_prev05_true_filt_scale)
smin <-min(sample_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale)) # smin = 4207

# Acidobacteriota
phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido <- phy_16_prune_soil_noncontam_prev05_true_filt_scale %>%
  subset_taxa(
    Phylum  == "Acidobacteriota"
  )
phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido # 848 taxa and 50 samples
get_taxa_unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido, taxonomic.rank = "Class")

merge_less_than_top_prev05_top20_acido <- function(phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido, top=19){
  transformed <- transform_sample_counts(phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido, function(x) x/sum(x))
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

# Family level top 20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido, "Family")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20 <- merge_less_than_top_prev05_top20_acido(phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam, top=19)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20) # Melt to long format
phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20_agr_soil = aggregate(Abundance~Soil.Type+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20_agr_site$Family)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20_agr_site, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20_agr_site.csv")

# open fixed dataframe including Others category
phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20_agr_site_new <- read.csv("exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20_agr_site_new.csv")

phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20_agr_site_new$Family <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20_agr_site_new$Family,
                                                                                                 levels = c("Acidobacteriaceae SB1", "Acidobacteriales_fam", "Acidobacteriota_SG22_fam",
                                                                                                 "Acidobacteriota_SG5_fam", "AT-s3-28_fam", "Blastocatellaceae", "Blastocatellia_fam",
                                                                                                 "Bryobacteraceae", "DS-100", "Holophagae SG7", "Koribacteraceae", "Paludibaculum", "Pyrinomonadaceae",
                                                                                                 "Solibacteraceae", "Thermoanaerobaculaceae", "Vicinamibacteraceae", "Vicinamibacterales", "Vicinamibacterales_fam",
                                                                                                 "Vicinamibacteria SG17", "Other"))
ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20_agr_site_new, aes(x = Sample.Site, y = Abundance, fill = Family)) + #plotting by sample
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
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
  xlab("Habitat") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_acido_fam_top20_agr_site_new.pdf", width = 10, height = 5, dpi = 200) # save graphic


# Actinobacteriota
phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino <- phy_16_prune_soil_noncontam_prev05_true_filt_scale %>%
  subset_taxa(
    Phylum  == "Actinobacteriota"
  )
phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino # 1588 taxa and 50 samples
get_taxa_unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino, taxonomic.rank = "Phylum")

merge_less_than_top_prev05_top20_actino <- function(phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino, top=19){
  transformed <- transform_sample_counts(phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino, function(x) x/sum(x))
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

# Family level top 20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino, "Family")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20 <- merge_less_than_top_prev05_top20_actino(phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam, top=19)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20) # Melt to long format
phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20_agr_soil = aggregate(Abundance~Soil.Type+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20_agr_site$Family)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20_agr_site, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20_agr_site.csv")

phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20_agr_site_new <- read.csv("exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20_agr_site_new.csv")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20_agr_site_new$Family <-factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20_agr_site_new$Family,
                                                                                                 levels = c("67-14", "Corynebacteriaceae", "Gaiellaceae",
                                                                                                            "Gaiellales_fam", "Geodermatophilaceae", "Intrasporangiaceae",
                                                                                                            "Kineosporiaceae", "MB-A2-108_fam", "Microbacteriaceae",
                                                                                                            "Micrococcaceae", "Micromonosporaceae", "Mycobacteriaceae",
                                                                                                            "Nocardioidaceae", "Pseudonocardiaceae", "Rubrobacteriaceae",
                                                                                                            "Solirubrobacteraceae", "Streptomycetaceae", "Streptosporangiaceae",
                                                                                                            "Thermomonosporaceae", "Other"))

ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20_agr_site_new, aes(x = Sample.Site, y = Abundance, fill = Family)) + #plotting by sample
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
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
  xlab("Habitat") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_actino_fam_top20_agr_site_new.pdf", width = 10, height = 5, dpi = 200) # save graphic


# Alphaproteobacteria
phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha <- phy_16_prune_soil_noncontam_prev05_true_filt_scale %>%
  subset_taxa(
    Phylum  == "Proteobacteria" &
      Class == "Alphaproteobacteria"
  )
phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha # 1335 taxa and 50 samples
get_taxa_unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha, taxonomic.rank = "Class")

merge_less_than_top_prev05_top20_alpha <- function(phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha, top=19){
  transformed <- transform_sample_counts(phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha, function(x) x/sum(x))
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

# Family level top 20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha, "Family")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20 <- merge_less_than_top_prev05_top20_alpha(phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam, top=19)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20) # Melt to long format
phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20_agr_soil = aggregate(Abundance~Soil.Type+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20_agr_site$Family)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20_agr_site, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20_agr_site.csv")

phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20_agr_site_new <- read.csv("exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20_agr_site_new.csv")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20_agr_site_new$Family <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20_agr_site_new$Family,
                                                                                                 levels = c("Acetobacteraceae", "Azospirillaceae", "Beijerinckiaceae",
                                                                                                            "Caulobacteraceae", "Devosiaceae", "Dongiaceae",
                                                                                                            "Elsterales_fam", "Geminicoccaceae", "Hyphomicrobiaceae",
                                                                                                            "Hyphomonadaceae", "Methyloligellaceae", "Paracaedibacteraceae",
                                                                                                            "Reyranellaceae", "Rhizobiaceae", "Rhizobiales_fam", 
                                                                                                            "Rhodobacteraceae", "Rickettsiaceae", "Sphingomonadaceae",
                                                                                                            "Xanthobacteraceae", "Others"))
ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20_agr_site_new, aes(x = Sample.Site, y = Abundance, fill = Family)) + #plotting by sample
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
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
  xlab("Habitat") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_alpha_fam_top20_agr_site_new.pdf", width = 10, height = 5, dpi = 200) # save graphic

# Gammaproteobacteria
phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma <- phy_16_prune_soil_noncontam_prev05_true_filt_scale %>%
  subset_taxa(
    Phylum  == "Proteobacteria" &
      Class == "Gammaproteobacteria"
  )
phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma # 737 taxa and 50 samples
get_taxa_unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma, taxonomic.rank = "Family")

merge_less_than_top_prev05_top20_gamma <- function(phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma, top=19){
  transformed <- transform_sample_counts(phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma, function(x) x/sum(x))
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

# Family level top 20
phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam <- tax_glom(phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma, "Family")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20 <- merge_less_than_top_prev05_top20_gamma(phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam, top=19)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20) # Melt to long format
phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20_agr_soil = aggregate(Abundance~Soil.Type+Family, data=phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20_agr_site$Family)
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20_agr_site, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20_agr_site.csv")

phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20_agr_site_new <- read.csv("exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20_agr_site_new.csv")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20_agr_site_new$Family <- factor(phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20_agr_site_new$Family,
                                                                                                 levels = c(
                                                                                                   "Alcaligenaceae", "Beggiatoaceae", "Burkholderiaceae",
                                                                                                   "CCD24","Comamonadaceae","Diplorickettsiaceae",
                                                                                                   "Enterobacteriaceae","Gammaproteobacteria_fam","Moraxellaceae",
                                                                                                   "Neisseriaceae","Nitrosomonadaceae","Oxalobacteraceae",
                                                                                                   "Pasteurellaceae","Pseudomonadaceae","Rhodocyclaceae",
                                                                                                   "SC-I-84","Steroidobacteraceae","TRA3-20",
                                                                                                   "Xanthomonadaceae","Others"))
ggplot(phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20_agr_site_new, aes(x = Sample.Site, y = Abundance, fill = Family)) + #plotting by sample
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
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
  xlab("Habitat") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_gamma_fam_top20_agr_site_new.pdf", width = 10, height = 5, dpi = 200) # save graphic

########################### Only Archaea ########################### 
phy_16_prune_soil_noncontam_prev05_true_filt_scale <- readRDS("data/phy_16_prune_soil_noncontam_prev05_true_filt_scale.RDS")
class(phy_16_prune_soil_noncontam_prev05_true_filt_scale)
smin <-min(sample_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale)) # smin = 4207

phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea <- phy_16_prune_soil_noncontam_prev05_true_filt_scale %>%
  subset_taxa(
    Domain  != "Bacteria"
  )
phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea # 75 taxa and 50 samples
get_taxa_unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea, taxonomic.rank = "Domain") # Only Archaea
phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea_df <- as.data.frame(otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea, "Matrix"))
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea_df, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea_df.csv")


phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria <- phy_16_prune_soil_noncontam_prev05_true_filt_scale %>%
  subset_taxa(
    Domain  != "Archaea"
  )
phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria # 9547 taxa and 50 samples
get_taxa_unique(phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria, taxonomic.rank = "Domain") # Only bacteria
phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria_df <- as.data.frame(otu_table(phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria, "Matrix"))
write.csv(phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria_df, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria_df.csv")


# Here we use the library microbiome to transform the data matrix
phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea_comp <- microbiome::transform(phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea, "compositional")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea_hel <- microbiome::transform(phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea, "hellinger")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria_comp <- microbiome::transform(phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria, "compositional")
phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria_hel <- microbiome::transform(phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria, "hellinger")


# NMDS ordination on the different transformed matrices
# All true samples; controls and blank samples and mitochondria/chloroplasts were removed
set.seed(1)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea_comp_nmds <- ordinate(
  phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea_comp, 
  method = "NMDS", 
  distance = "bray"
) # Stress:     0.06281775

set.seed(1)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea_hel_nmds <- ordinate(
  phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea_hel, 
  method = "NMDS", 
  distance = "bray"
) # Stress: 0.09943494

set.seed(1)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria_comp_nmds <- ordinate(
  phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria_comp, 
  method = "NMDS", 
  distance = "bray"
) # Stress:     0.07360296

set.seed(1)
phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria_hel_nmds <- ordinate(
  phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria_hel, 
  method = "NMDS", 
  distance = "bray"
) # Stress: 0.0794533

# nMDS plot compositional Archaea
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea_comp,
  ordination = phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea_comp_nmds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E","#C29B6C", 
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
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  #ggrepel::geom_text_repel(aes(label = Sample.Site, color = Habitat), size = 2, max.overlaps = 30) +
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
  annotate("text", x = 0.5, y = 0.5, label ="2D Stress: 0.06")
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea_comp_nmds.pdf", width = 6, height = 4, dpi = 200)

# nMDS plot hellinger Archaea
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea_hel,
  ordination = phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea_hel_nmds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E","#C29B6C", 
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
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  #ggrepel::geom_text_repel(aes(label = Sample.Site, color = Habitat), size = 2, max.overlaps = 30) +
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
  annotate("text", x = 0.7, y = 0.45, label ="2D Stress: 0.1")
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea_hel_nmds.pdf", width = 6, height = 4, dpi = 200)

# nMDS plot compositional Bacteria
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria_comp,
  ordination = phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria_comp_nmds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E","#C29B6C", 
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
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  #ggrepel::geom_text_repel(aes(label = Sample.Site, color = Habitat), size = 2, max.overlaps = 30) +
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
  annotate("text", x = 1, y = 1, label ="2D Stress: 0.07")
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria_comp_nmds.pdf", width = 6, height = 4, dpi = 200)

# nMDS plot hellinger Bacteria
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria_hel,
  ordination = phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria_hel_nmds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E","#C29B6C", 
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
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  #ggrepel::geom_text_repel(aes(label = Sample.Site, color = Habitat), size = 2, max.overlaps = 30) +
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
  annotate("text", x = 1, y = 1, label ="2D Stress: 0.08")
ggsave("results/soil/phy_16_prune_soil_noncontam_prev05_true_filt_scale_bacteria_hel_nmds.pdf", width = 6, height = 4, dpi = 200)


alpha_div_archaea_soil <- data.frame(
  "Observed" = phyloseq::estimate_richness(phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea, measures = "Shannon"),
  "InvSimpson" = phyloseq::estimate_richness(phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea, measures = "InvSimpson"),
  "Reads" = phyloseq::sample_sums(phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea),
  "Soil.Type" = phyloseq::sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea)$Soil.Type,
  "Habitat" = phyloseq::sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea)$Habitat,
  "Sample.Site" = phyloseq::sample_data(phy_16_prune_soil_noncontam_prev05_true_filt_scale_archaea)$Sample.Site)
alpha_div_archaea_soil$Evenness <- alpha_div_archaea_soil$Shannon/log(alpha_div_archaea_soil$Observed)
head(alpha_div_archaea_soil)

# Rename variable InvSimpson to Simpson
# The function rename & %>% works on dplyr. make sure it is loaded.
alpha_div_archaea_soil <- alpha_div_archaea_soil %>%
  rename(Simpson = InvSimpson)
head(alpha_div_archaea_soil)

# Reorder dataframe, first categorical then numerical variables  
alpha_div_archaea_soil <- alpha_div_archaea_soil[, c(7, 6, 5, 4, 1, 2, 3, 8)]
head(alpha_div_archaea_soil)
write.csv(alpha_div_archaea_soil, "exported_tables/soil/alpha_div_archaea_soil.csv")

#Summarize alpha diversity measures by site
summary_alpha_div_archaea_soil_site <- alpha_div_archaea_soil %>%
  group_by(Soil.Type, Habitat, Sample.Site) %>%
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
write.csv(summary_alpha_div_archaea_soil_site, "exported_tables/soil/summary_alpha_div_archaea_soil_site.csv")

#Summarize alpha diversity measures by habitat
summary_alpha_div_archaea_soil_habitat <- alpha_div_archaea_soil %>%
  group_by(Soil.Type, Habitat) %>%
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
write.csv(summary_alpha_div_archaea_soil_habitat, "exported_tables/soil/summary_alpha_div_archaea_soil_habitat.csv")

#Summarize alpha diversity measures by soil type
summary_alpha_div_archaea_soil_type <- alpha_div_archaea_soil %>%
  group_by(Soil.Type) %>%
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
write.csv(summary_alpha_div_archaea_soil_type, "exported_tables/soil/summary_alpha_div_archaea_soil_type.csv")

# KW analysis alpha_16S_soil on all metrics at once
# Remember, numerical variables are from columns 4-8. The test is by soil type, column 3
kw_alpha_div_archaea_soil_type_univ <- as.data.frame(sapply(4:8, function(x) kruskal.test(alpha_div_archaea_soil[,x],
                                                                                         alpha_div_archaea_soil[,3])))

# Rename columns with the proper variable names
kw_alpha_div_archaea_soil_type_univ <- kw_alpha_div_archaea_soil_type_univ %>%
  rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)

# Make some adjustments on the table
kw_alpha_div_archaea_soil_type_univ <- t(kw_alpha_div_archaea_soil_type_univ) # transpose
kw_alpha_div_archaea_soil_type_univ <- as_tibble(kw_alpha_div_archaea_soil_type_univ, rownames = "Metric") # adding rownames as a column
kw_alpha_div_archaea_soil_type_univ <- kw_alpha_div_archaea_soil_type_univ[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_div_archaea_soil_type_univ) # checking object class

# Save resulting table with fwrite to avoid any issues with characters
data.table::fwrite(kw_alpha_div_archaea_soil_type_univ, "exported_tables/soil/kw_alpha_div_archaea_soil_type_univ.csv")

# KW analysis alpha_16S_soil on all metrics at once
# Remember, numerical variables are from columns 4-8. The test is by habitat, column 2
kw_alpha_div_archaea_habitat_univ <- as.data.frame(sapply(4:8, function(x) kruskal.test(alpha_div_archaea_soil[,x],
                                                                                        alpha_div_archaea_soil[,2])))

kw_alpha_div_archaea_habitat_univ <- kw_alpha_div_archaea_habitat_univ %>%
  rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)

# Make some adjustments on the table
kw_alpha_div_archaea_habitat_univ <- t(kw_alpha_div_archaea_habitat_univ) # transpose
kw_alpha_div_archaea_habitat_univ <- as_tibble(kw_alpha_div_archaea_habitat_univ, rownames = "Metric") # adding rownames as a column
kw_alpha_div_archaea_habitat_univ <- kw_alpha_div_archaea_habitat_univ[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_div_archaea_habitat_univ) # checking object class

# Save resulting table with fwrite to avoid and issues with characters
data.table::fwrite(kw_alpha_div_archaea_habitat_univ, "exported_tables/soil/kw_alpha_div_archaea_habitat_univ.csv")

# pairwise comparison including p-value adjustment
kw_alpha_div_archaea_soil_univ_pwt <- as.data.frame(sapply(4:8, function(x) pairwise.t.test(alpha_div_archaea_soil[,x],
                                                                                            alpha_div_archaea_soil[,2], p.adjust.method = "BH")))
kw_alpha_div_archaea_soil_univ_pwt <- kw_alpha_div_archaea_soil_univ_pwt %>%
  rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)

# Make some adjustments on the table
kw_alpha_div_archaea_soil_univ_pwt <- t(kw_alpha_div_archaea_soil_univ_pwt) # transpose
kw_alpha_div_archaea_soil_univ_pwt <- as_tibble(kw_alpha_div_archaea_soil_univ_pwt, rownames = "Metric") # adding rownames as a column
class(kw_alpha_div_archaea_soil_univ_pwt) # checking object class

# Save resulting table with fwrite to avoid and issues with characters
fwrite(kw_alpha_div_archaea_soil_univ_pwt, "exported_tables/soil/kw_alpha_div_archaea_soil_univ_pwt.csv")

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

# Plot alpha-diversity per site within habitat/soil type
alpha_div_archaea_soil %>%
  gather(key = metric, value = value, c("Reads", "Observed", "Shannon", "Simpson", "Evenness")) %>%
  mutate(metric = factor(metric, levels = c("Reads", "Observed", "Shannon", "Simpson", "Evenness"))) %>%
  ggplot(aes(x = Sample.Site, y = value)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  geom_jitter(aes(color = Habitat), height = 0, width = .2) +
  facet_nested(metric ~ Soil.Type+Habitat, scales = "free", labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  scale_color_manual(values = alpha_color_habitat) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.1, 0.1)) + # plot as % and removes the internal margins
  scale_x_discrete(labels = site.labs,
                   expand = c(0.1, 0.1),
                   drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(legend.position="none")
ggsave("results/soil/16S_archaea_alpha_diversity_phy_soil_16S_prune_noncontam_prev05_true_scale_site.pdf", width = 8, height = 6, dpi = 150) # save graphic

alpha_div_archaea_soil %>%
  gather(key = metric, value = value, c("Reads", "Observed", "Shannon", "Simpson", "Evenness")) %>%
  mutate(metric = factor(metric, levels = c("Reads", "Observed", "Shannon", "Simpson", "Evenness"))) %>%
  ggplot(aes(x = Habitat, y = value)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  geom_jitter(aes(color = Habitat), height = 0, width = .2) +
  facet_nested(metric ~ Soil.Type, scales = "free", labeller = labeller(Soil.Type = soil.labs)) +
  scale_color_manual(values = alpha_color_habitat) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.1, 0.1)) + # plot as % and removes the internal margins
  scale_x_discrete(labels = habitat.labs,
                   expand = c(0.1, 0.1),
                   drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(legend.position="none")
ggsave("results/soil/16S_archaea_alpha_diversity_phy_soil_16S_prune_noncontam_prev05_true_scale_habitat.pdf", width = 8, height = 6, dpi = 150) # save graphic
