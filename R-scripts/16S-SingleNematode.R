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

# Import asv table, taxonomy table, and sample table in excel format. It can also be imported as csv files.
asv_tab_16S <- read_excel("data/ASV_table_Nov_06_2021.xlsx")
tax_tab_16S <- read_excel("data/Taxonomy_table_Jan_17_2022.xlsx")
sample_tab_16S <- read_excel("data/Sample_table_16S_Mar_31_2023.xlsx")

# Convert asv and taxonomy tables to matrix format; sample table as dataframe
asv_tab_16S <- as.matrix(asv_tab_16S)
class(asv_tab_16S) # check the object class
tax_tab_16S <- as.matrix(tax_tab_16S)
sample_tab_16S <- as.data.frame(sample_tab_16S)
class(sample_tab_16S) # check the object class

# Rename row names using ASV column for asv table and taxonomy table; rename using Sample column for sample table
rownames(asv_tab_16S) <- asv_tab_16S[,1]
rownames(tax_tab_16S) <- tax_tab_16S[,1]
rownames(sample_tab_16S) <- sample_tab_16S[,1]

# Exclude first column from all three tables; set asv table as numeric
asv_tab_16S <- asv_tab_16S[,-1]
class(asv_tab_16S) <- "numeric"
tax_tab_16S <- tax_tab_16S[,-1]
sample_tab_16S <- sample_tab_16S[,-1]

#Transform each table/matrix to phyloseq objects
ASV_16S = otu_table(asv_tab_16S, taxa_are_rows = TRUE)
TAX_16S = tax_table(tax_tab_16S)
samples_16S = sample_data(sample_tab_16S)

# Create a phyloseq object by combine all three previous objects
phy_16S <- phyloseq(ASV_16S, TAX_16S, samples_16S)
phy_16S # 12818 taxa and 624 samples
head(sample_data(phy_16S)) # check the first 6 lines of your phyloseq object

#Visualize some data entries from the phyloseq object
sample_names(phy_16S) # list sample names
rank_names(phy_16S) # list the taxonomic ranks from your tax_table
sample_variables(phy_16S) # list the variables from your sample_data

# Check if there are samples with zero reads and remove them
smin <- min(sample_sums(phy_16S)) # smin = 0 means you have samples with no reads. Orignally sample number was 624.
phy_16S_prune <- prune_samples(sample_sums(phy_16S)>0, phy_16S)
smin <- min(sample_sums(phy_16S_prune))
phy_16S_prune # now the phyloseq object has 620 samples, with a smin = 2

# Check for ASVs that have no counted reads
any(taxa_sums(phy_16S_prune) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16S_prune) == 0) # gives the number of cases

# Inspect Library sizes
df_phy_16S_prune <- as.data.frame(sample_data(phy_16S_prune)) # transform sample_data from phyloseq object to a dataframe
df_phy_16S_prune$LibrarySize <- sample_sums(phy_16S_prune) # create another variable LibrarySize by summing the read count across all ASVs per sample
df_phy_16S_prune <- df_phy_16S_prune[order(df_phy_16S_prune$LibrarySize),] # reorder samples by library size
df_phy_16S_prune$Index <- seq(nrow(df_phy_16S_prune)) # create another varible Index by using the row number of the dataframe

# Plot library sizes using ggplot
ggplot(data = df_phy_16S_prune, aes(x= Index, y= LibrarySize, color = Sample.Type)) +
  geom_point()
ggsave("results/all-dataset/16S_prune_library-size-sample-type-all-samples.pdf", width = 8, height = 6, dpi = 150)

####################################### Only nematode samples #######################################
# Keep only samples from the microbiome dataset (i.e., nematodes)
phy_16S_prune_nem <- subset_samples(phy_16S_prune, Experiment == "Microbiome") # select only microbiome samples
phy_16S_prune_nem # 12818 taxa and 562 samples

# Check/remove samples with zero reads if necessary. We already remove samples with 0 reads, so these commands have no effect
smin <- min(sample_sums(phy_16S_prune_nem))
phy_16S_prune_nem <- prune_samples(sample_sums(phy_16S_prune_nem)>0, phy_16S_prune_nem)
smin <- min(sample_sums(phy_16S_prune_nem)) # smin = 2
phy_16S_prune_nem # 12818 taxa and 562 samples

# Check for ASVs that have no counted reads
any(taxa_sums(phy_16S_prune_nem) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16S_prune_nem) == 0) # gives the number of cases
phy_16S_prune_nem <-prune_taxa(taxa_sums(phy_16S_prune_nem) > 0, phy_16S_prune_nem)
phy_16S_prune_nem # 3341 taxa and 562 samples
head(sample_data(phy_16S_prune_nem)) # list 6 first rows from sample_data
head(otu_table(phy_16S_prune_nem)) # list 6 first rows from otu_table

# Inspect Library sizes of nema samples separately
df_phy_16S_prune_nem <- as.data.frame(sample_data(phy_16S_prune_nem))
df_phy_16S_prune_nem$LibrarySize <- sample_sums(phy_16S_prune_nem)
df_phy_16S_prune_nem <- df_phy_16S_prune_nem[order(df_phy_16S_prune_nem$LibrarySize),]
df_phy_16S_prune_nem$Index <- seq(nrow(df_phy_16S_prune_nem))

# Plot using ggplot for nema samples
ggplot(data = df_phy_16S_prune_nem, aes(x= Index, y= LibrarySize, color = Sample.Type)) +
  geom_point()
ggsave("results/nema/16S_nem_library-size-sample-type-nem-samples.pdf", width = 8, height = 6, dpi = 150)

ggplot(data = df_phy_16S_prune_nem, aes(x= Index, y= LibrarySize, color = Sample.Control)) +
  geom_point()
ggsave("results/nema/16S_nem_library-size-sample-control-nem-samples.pdf", width = 8, height = 6, dpi = 150)

# Identify potential contaminants using R package decontam - Frequency method, it is based on the DNA concentration Qubit value from each sample
phy_16S_prune_nem_contamdf_freq <- isContaminant(phy_16S_prune_nem, method="frequency", conc="QubitConc")
head(phy_16S_prune_nem_contamdf_freq)
tail(phy_16S_prune_nem_contamdf_freq)
table(phy_16S_prune_nem_contamdf_freq$contaminant) # 28 TRUE potential contaminants
head(which(phy_16S_prune_nem_contamdf_freq$contaminant))
phy_16S_prune_nem_contamdf_freq_list <- rownames_to_column(phy_16S_prune_nem_contamdf_freq, var = "ASV")
write.csv(phy_16S_prune_nem_contamdf_freq_list, "exported_tables/nema/phy_16S_prune_nem_contamdf_freq_list.csv")

# Explore AVSs that are potentially contaminants, ASV abundance should be inversely correlated with DNA concentration
plot_frequency(phy_16S_prune_nem, taxa_names(phy_16S_prune_nem)[c(65, 135, 142, 184)], conc="QubitConc") + 
  xlab("DNA Concentration (ng/ul)")
ggsave("results/nema/16S_phy_16S_prune_nem_contaminants.pdf", width = 8, height = 6, dpi = 150)

# Explore AVSs that are potentially contaminants, randomly choosing 4 AVSs.
set.seed(100)
plot_frequency(phy_16S_prune_nem, taxa_names(phy_16S_prune_nem)[sample(which(phy_16S_prune_nem_contamdf_freq$contaminant),4)], conc="QubitConc") + 
  xlab("DNA Concentration (ng/ul)")
ggsave("results/nema/16S_phy_16S_prune_nem_contaminants_random.pdf", width = 8, height = 6, dpi = 150)

# Identify potential contaminants using package decontam - Prevalence method, it is based on the negative controls and blank samples included in the experiment
sample_data(phy_16S_prune_nem)$is.neg <- sample_data(phy_16S_prune_nem)$Sample.Control == "Control.Sample"
phy_16S_prune_nem_contamdf_prev <- isContaminant(phy_16S_prune_nem, method="prevalence", neg="is.neg")
table(phy_16S_prune_nem_contamdf_prev$contaminant) # 53 TRUE potential contaminants
head(which(phy_16S_prune_nem_contamdf_prev$contaminant))
phy_16S_prune_nem_contamdf_prev_list <- rownames_to_column(phy_16S_prune_nem_contamdf_prev, var = "ASV")
write.csv(phy_16S_prune_nem_contamdf_prev_list, "exported_tables/nema/phy_16S_prune_nem_contamdf_prev_list.csv")

# Identify Contaminants using package decontam - Prevalence method, increase threshold value to 0.5
phy_16S_prune_nem_contamdf_prev05 <- isContaminant(phy_16S_prune_nem, method="prevalence", neg="is.neg", threshold=0.5)
table(phy_16S_prune_nem_contamdf_prev05$contaminant) # 114 TRUE potential contaminants
phy_16S_prune_nem_contamdf_prev05_list <- rownames_to_column(phy_16S_prune_nem_contamdf_prev05, var = "ASV")
write.csv(phy_16S_prune_nem_contamdf_prev05_list, "exported_tables/nema/phy_16S_prune_nem_contamdf_prev05_list.csv")

# Make phyloseq object of presence-absence in negative controls and true samples
phy_16S_prune_nem_pa <- transform_sample_counts(phy_16S_prune_nem, function(abund) 1*(abund>0))
phy_16S_prune_nem_pa_neg <- prune_samples(sample_data(phy_16S_prune_nem_pa)$Sample.Control == "Control.Sample", phy_16S_prune_nem_pa)
phy_16S_prune_nem_pa_pos <- prune_samples(sample_data(phy_16S_prune_nem_pa)$Sample.Control == "True.Sample", phy_16S_prune_nem_pa)

# Make data.frame of prevalence in positive and negative samples, threshold of 0.1
df_phy_16S_prune_nem_pa <- data.frame(pa.pos=taxa_sums(phy_16S_prune_nem_pa), pa.neg=taxa_sums(phy_16S_prune_nem_pa_neg),
                                      contaminant=phy_16S_prune_nem_contamdf_prev$contaminant)
ggplot(data=df_phy_16S_prune_nem_pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
ggsave("results/nema/phy_16S_prune_nema_contamdf_prev_0.1.pdf", width = 8, height = 6, dpi = 150)

# Make data.frame of prevalence in positive and negative samples, threshold of 0.5
df_phy_16S_prune_nem_pa05 <- data.frame(pa.pos=taxa_sums(phy_16S_prune_nem_pa), pa.neg=taxa_sums(phy_16S_prune_nem_pa_neg),
                                        contaminant=phy_16S_prune_nem_contamdf_prev05$contaminant)
ggplot(data=df_phy_16S_prune_nem_pa05, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
ggsave("results/nema/phy_16S_prune_nema_contamdf_prev_0.5.pdf", width = 8, height = 6, dpi = 150)

# Identify Contaminants using package decontam -  Combined method (i.e., Frequency and Prevalence)
phy_16S_prune_nem_contamdf_comb <- isContaminant(phy_16S_prune_nem, conc="QubitConc", neg="is.neg", threshold=0.1, detailed = TRUE, normalize = TRUE, method="combined")
table(phy_16S_prune_nem_contamdf_comb$contaminant) # 19 TRUE potential contaminants
head(which(phy_16S_prune_nem_contamdf_comb$contaminant))

# Combined method, threshold = 0.5
phy_16S_prune_nem_contamdf_comb05 <- isContaminant(phy_16S_prune_nem, conc="QubitConc", neg="is.neg", threshold=0.5, detailed = TRUE, normalize = TRUE, method="combined")
table(phy_16S_prune_nem_contamdf_comb05$contaminant) # 185 TRUE potential contaminants
head(which(phy_16S_prune_nem_contamdf_comb05$contaminant))

# Remove contaminants from phyloseq object with method prev05. Manual inspection suggests this is the most appropriated method
# Check for samples with zero reads as it seems to be the most appropriate for this dataset.
# 16S nema samples
phy_16S_prune_nem_noncontam_prev05 <- prune_taxa(!phy_16S_prune_nem_contamdf_prev05$contaminant, phy_16S_prune_nem)
phy_16S_prune_nem_noncontam_prev05 # 3227 taxa and 562 samples
smin <-min(sample_sums(phy_16S_prune_nem_noncontam_prev05)) # smin = 0

# remove 16S nema samples with zero reads
phy_16S_prune_nem_noncontam_prev05 <- prune_samples(sample_sums(phy_16S_prune_nem_noncontam_prev05)>0, phy_16S_prune_nem_noncontam_prev05)
any(sample_sums(phy_16S_prune_nem_noncontam_prev05) == 0) # if TRUE, there are samples with no reads, otherwise FALSE
sum(sample_sums(phy_16S_prune_nem_noncontam_prev05) == 0) # gives the number of cases
phy_16S_prune_nem_noncontam_prev05 # 3227 taxa and 559 samples
smin <-min(sample_sums(phy_16S_prune_nem_noncontam_prev05)) # smin = 2

# Check for ASVs that have no counted reads
any(taxa_sums(phy_16S_prune_nem_noncontam_prev05) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16S_prune_nem_noncontam_prev05) == 0) # gives the number of cases

# After removing potential contaminants with decontam, we can focus on the true samples only
# Subset by sample types, that is, only nema samples removing controls and blanks since we already filtered potential contaminants
phy_16S_prune_nem_noncontam_prev05_true <- subset_samples(phy_16S_prune_nem_noncontam_prev05, Sample.Type =="SingleWorm") %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_16S_prune_nem_noncontam_prev05_true # 3130 taxa and 530 samples
any(taxa_sums(phy_16S_prune_nem_noncontam_prev05_true) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16S_prune_nem_noncontam_prev05_true) == 0) # gives the number of cases
smin <-min(sample_sums(phy_16S_prune_nem_noncontam_prev05_true)) # smin = 2

# Checking the total number of reads and its distribution across samples
phy_16S_prune_nem_noncontam_prev05_true_readsums_df = data.frame(nreads = sort(taxa_sums(phy_16S_prune_nem_noncontam_prev05_true), TRUE),
                                                                 sorted = 1:ntaxa(phy_16S_prune_nem_noncontam_prev05_true), 
                                                                 type = "ASVs")
phy_16S_prune_nem_noncontam_prev05_true_readsums_df = rbind(phy_16S_prune_nem_noncontam_prev05_true_readsums_df, 
                                                            data.frame(nreads = sort(sample_sums(phy_16S_prune_nem_noncontam_prev05_true),TRUE), 
                                                                       sorted = 1:nsamples(phy_16S_prune_nem_noncontam_prev05_true), type = "Samples"))
title = "Total number of reads"
p = ggplot(phy_16S_prune_nem_noncontam_prev05_true_readsums_df, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
ggsave("results/nema/16S_nem_ShSk_ASV_Read_Profile_filt_prev05-nem-samples.pdf", width = 12, height = 6, dpi = 150)

# Filter out Chloroplasts and Mitochondria from the data matrix.
phy_16_prune_nem_noncontam_prev05_true_all_filt <- phy_16S_prune_nem_noncontam_prev05_true %>%
  subset_taxa(
    Family  != "Mitochondria" &
      Order   != "Chloroplast"
  )
phy_16_prune_nem_noncontam_prev05_true_all_filt # 3077 taxa and 530 samples
any(taxa_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt) == 0) # gives the number of cases
smin <-min(sample_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt)) # smin = 2

# Export and Transpose matrix for Picrust2 analysis
# ASV 16S nema table
ASV_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt = as(otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt), "matrix")
write.csv(ASV_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt, "exported_tables/nema/ASV_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt.csv")

# Extract metadata from phyloseq object
# 16S nema sample table
TABLE_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt = as(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt), "matrix")
class(TABLE_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt)
TABLE_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_df <- as.data.frame(TABLE_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt)
class(TABLE_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_df)
write.csv(TABLE_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_df, "exported_tables/nema/TABLE_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_df.csv")

# 16S nema taxonomy
TAX_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt  = as(tax_table(phy_16_prune_nem_noncontam_prev05_true_all_filt), "matrix")
TAX_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt <- as.data.frame(TAX_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt)
write.csv(TAX_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt, "exported_tables/nema/TAX_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt.csv")

################# Save/import clean phyloseq object for further analysis ################# 
# Save as RDS
saveRDS(phy_16_prune_nem_noncontam_prev05_true_all_filt, "data/phy_16_prune_nem_noncontam_prev05_true_all_filt.RDS")  

# Import phyloseq
phy_16_prune_nem_noncontam_prev05_true_all_filt <- readRDS("data/phy_16_prune_nem_noncontam_prev05_true_all_filt.RDS")
class(phy_16_prune_nem_noncontam_prev05_true_all_filt)

#################
# Apply different transformations to ASV table, e.g., converting number of reads to relative abundance
# Here we use the library microbiome
phy_16_prune_nem_noncontam_prev05_true_all_filt_comp <- microbiome::transform(phy_16_prune_nem_noncontam_prev05_true_all_filt, "compositional")
# phy_16_prune_nem_noncontam_prev05_true_all_filt_comp = transform_sample_counts(phy_16_prune_nem_noncontam_prev05_true_all_filt, function(x) 100 * x/sum(x)). # This command is the equivalent to the microbiom::transform function "compositional"
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr <- microbiome::transform(phy_16_prune_nem_noncontam_prev05_true_all_filt, "clr")
phy_16_prune_nem_noncontam_prev05_true_all_filt_log10 <- microbiome::transform(phy_16_prune_nem_noncontam_prev05_true_all_filt, "log10p")
phy_16_prune_nem_noncontam_prev05_true_all_filt_hellinger <- microbiome::transform(phy_16_prune_nem_noncontam_prev05_true_all_filt, "hellinger")
# This command is the equivalent to the microbiom::transform function "hellinger"
# phy_16_prune_nem_noncontam_prev05_true_all_filt_hellinger <- transform_sample_counts(phy_16_prune_nem_noncontam_prev05_true_all_filt, function(x) sqrt(x/sum(x)))

########## Save clr data for analysis in Primer
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_asv_primer <- as.data.frame(otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr, "matrix"))
write.csv(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_asv_primer, "exported_tables/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_asv_primer.csv")
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_metada_primer <- as.data.frame(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr))
write.csv(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_metada_primer, "exported_tables/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_metada_primer.csv")
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_tax_primer <- as.data.frame(tax_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr))
write.csv(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_tax_primer, "exported_tables/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_tax_primer.csv")

############### Ordination using PCA and clr transformed phyloseq objects
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca <- phyloseq::ordinate(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr,
  method = "RDA",
  distance = "euclidean"
)

# Plot scree plot to see PCA contributions
phyloseq::plot_scree(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/nema/16S_nem_ShSk_ASV_noncontam_prev05_true_all_filt_clr_pca.pdf", width = 12, height = 6, dpi = 150)

head(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca$CA$eig)
sapply(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca$CA$eig[1:5], function(x) x / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca$CA$eig))

# Scale axes and plot ordination
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca_clr1 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca$CA$eig[1] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca$CA$eig)
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca_clr2 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca$CA$eig[2] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca$CA$eig)

# Plot PCA based on clr transformation
soil.labs <- c("Compacted", "Uncompacted")
names(soil.labs) <- c("Compact.Soil", "Not.Compact.Soil")

habitat.labs <- c("CHA", "CSS", "NGR", "HLC", "OWL", "RIP")
names(habitat.labs) <- c("Chaparral", "Coastal.scrub", "Native.grass", "Hollyleaf.cherry", "Oak.wood", "Riparian")

# change shape variable (Infraorder or Feeding.Group) depending on desired graphic
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca,
  type="samples",
  color="Sample.Site",
  shape = "Infraorder") +
  facet_nested(~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs))+
  geom_point(aes(color = Sample.Site), alpha = 0.7, size = 3) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A",
                                "#984EA3", "#FF7F00", "#FFFF33",
                                "#A65628", "#F781BF", "#999999",
                                "#1B9E77", "#D95F02", "#7570B3",
                                "#E7298A", "#66A61E", "#E6AB02",
                                "#A6761D", "#666666", "lightblue"),
                     name = "Site",
                     breaks = c("SK.01", "SK.02", "SK.03",
                                "SK.13", "SK.14", "SK.15",
                                "SK.16", "SK.17", "SK.18",
                                "SK.04", "SK.05", "SK.06",
                                "SK.07", "SK.08", "SK.09",
                                "SK.10", "SK.11", "SK.12"),
                     labels = c("01", "02", "03",
                                "13", "14", "15",
                                "16", "17", "18",
                                "04", "05", "06",
                                "07", "08", "09",
                                "10", "11", "12"))+
  scale_shape_manual(values = c(15, 2, 3, 4, 8, 16, 17, 18, 9, 10, 1, 12),
                     name = "Infraorder")+
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
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5, legend.position = "bottom") + # adjusts the title of the legend
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold"))
  #ylab("NMDS2") + # add the title on y axis
#xlab("NMDS1") + # add the title on x axis
#ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
#theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
#annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
#annotate("text", x = 3, y = 2, label ="2D Stress: 0.07") +
#coord_fixed(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca_clr2 / phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca_clr1) +
ggsave("results/nema/16S_nem_ShSk_ASV_noncontam_prev05_true_all_filt_clr_no-mito-chloro_pca_plot_sample_site_infraorder.pdf", width = 14, height = 6, dpi = 300)
ggsave("results/nema/16S_nem_ShSk_ASV_noncontam_prev05_true_all_filt_clr_no-mito-chloro_pca_plot_sample_site_feeding.pdf", width = 14, height = 6, dpi = 300)

# Without facet and only infraorder and soil type
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca,
  type="samples",
  color="Infraorder",
  shape = "Soil.Type") +
  #facet_nested(~ Soil.Type+Habitat, scales = "free",
   #            labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs))+
  #geom_point(aes(color = Infraorder), alpha = 0.7, size = 3) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A",
                                "#984EA3", "#FF7F00", "#FFFF33",
                                "#A65628", "#F781BF", "#999999",
                                "#1B9E77", "#D95F02", "#7570B3",
                                "#E7298A", "#66A61E", "#E6AB02",
                                "#A6761D", "#666666", "lightblue"),
                     name = "Nematode Group")+
  scale_shape_manual(values = c(19,17),
                     name = "Soil Type",
                     breaks = c("Compact.Soil", "Not.Compact.Soil"),
                     labels = c("Compacted", "Uncompacted"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Infraorder), alpha = 0.7, size = 3) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5, legend.position = "right") # adjusts the title of the legend
#ylab("NMDS2") + # add the title on y axis
#xlab("NMDS1") + # add the title on x axis
#ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
#theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
#annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
#annotate("text", x = 3, y = 2, label ="2D Stress: 0.07") +
#coord_fixed(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca_clr2 / phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca_clr1) +
ggsave("results/nema/16S_nem_ShSk_ASV_noncontam_prev05_true_all_filt_clr_no-mito-chloro_pca_plot_without_facet_infraorder.pdf", width = 7, height = 5, dpi = 300)


theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca,
  type="samples",
  color="Habitat",
  shape = "Soil.Type") +
  #facet_nested(~ Soil.Type+Habitat, scales = "free",
  #            labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs))+
  #geom_point(aes(color = Infraorder), alpha = 0.7, size = 3) +
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
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5, legend.position = "right") # adjusts the title of the legend
#ylab("NMDS2") + # add the title on y axis
#xlab("NMDS1") + # add the title on x axis
#ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
#theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
#annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
#annotate("text", x = 3, y = 2, label ="2D Stress: 0.07") +
#coord_fixed(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca_clr2 / phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca_clr1) +
ggsave("results/nema/16S_nem_ShSk_ASV_noncontam_prev05_true_all_filt_clr_no-mito-chloro_pca_plot_without_facet_soil_habitat.pdf", width = 7, height = 5, dpi = 300)


########### Define phyloseq objects per habitat
####### Coastal scrub
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr, Habitat == "Coastal.scrub") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs)$Habitat)

phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs,
               sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs) != "Nem.405" & # Plectoidea, SK03
               sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs) != "Nem.391") # Tylenchoidea, SK01

####### Chaparral
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr, Habitat == "Chaparral") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha)$Habitat)


####### Native grass
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr, Habitat == "Native.grass") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr)$Habitat)


####### Riparian
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr, Habitat == "Riparian") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip)$Habitat)

####### Hollyleaf.cherry
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr, Habitat == "Hollyleaf.cherry") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc)$Habitat)

####### Oakwood.land
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr, Habitat == "Oak.wood") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl)$Habitat)

phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl,
                                                                          sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl) != "Nem.479") # Dorylaimidae, SK11


####### PCA analysis for each phyloseq object
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs_pca <- phyloseq::ordinate(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs,
  method = "RDA",
  distance = "euclidean"
)

phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca <- phyloseq::ordinate(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha,
  method = "RDA",
  distance = "euclidean"
)

phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr_pca <- phyloseq::ordinate(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr,
  method = "RDA",
  distance = "euclidean"
)


phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip_pca <- phyloseq::ordinate(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip,
  method = "RDA",
  distance = "euclidean"
)

phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc_pca <- phyloseq::ordinate(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc,
  method = "RDA",
  distance = "euclidean"
)

phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl_pca <- phyloseq::ordinate(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl,
  method = "RDA",
  distance = "euclidean"
)


# Plot scree plot to see PCA contributions coastal scrub
phyloseq::plot_scree(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs_pca.pdf", width = 12, height = 6, dpi = 150)
head(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs_pca$CA$eig)
sapply(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs_pca$CA$eig[1:5], function(x) x / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs_pca$CA$eig))

# Plot scree plot to see PCA contributions chaparral
phyloseq::plot_scree(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca.pdf", width = 12, height = 6, dpi = 150)
head(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca$CA$eig)
sapply(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca$CA$eig[1:5], function(x) x / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca$CA$eig))

# Plot scree plot to see PCA contributions native grass
phyloseq::plot_scree(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr_pca.pdf", width = 12, height = 6, dpi = 150)
head(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr_pca$CA$eig)
sapply(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr_pca$CA$eig[1:5], function(x) x / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr_pca$CA$eig))

# Plot scree plot to see PCA contributions ripirian
phyloseq::plot_scree(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip_pca.pdf", width = 12, height = 6, dpi = 150)
head(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip_pca$CA$eig)
sapply(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip_pca$CA$eig[1:5], function(x) x / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip_pca$CA$eig))

# Plot scree plot to see PCA contributions holy-leaf cherry
phyloseq::plot_scree(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc_pca.pdf", width = 12, height = 6, dpi = 150)
head(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc_pca$CA$eig)
sapply(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc_pca$CA$eig[1:5], function(x) x / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc_pca$CA$eig))

# Plot scree plot to see PCA contributions oak woodland
phyloseq::plot_scree(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl_pca.pdf", width = 12, height = 6, dpi = 150)
head(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl_pca$CA$eig)
sapply(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl_pca$CA$eig[1:5], function(x) x / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl_pca$CA$eig))

# Scale axes and plot ordination coastal scrub
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs_pca_clr1 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs_pca$CA$eig[1] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs_pca$CA$eig)
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs_pca_clr2 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs_pca$CA$eig[2] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs_pca$CA$eig)
# Scale axes and plot ordination coastal chaparral
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca_clr1 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca$CA$eig[1] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca$CA$eig)
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca_clr2 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca$CA$eig[2] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca$CA$eig)
# Scale axes and plot ordination native grass
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr_pca_clr1 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca$CA$eig[1] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr_pca$CA$eig)
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr_pca_clr2 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca$CA$eig[2] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr_pca$CA$eig)
# Scale axes and plot ordination coastal ripirian
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip_pca_clr1 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca$CA$eig[1] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip_pca$CA$eig)
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip_pca_clr2 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca$CA$eig[2] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip_pca$CA$eig)
# Scale axes and plot ordination holy-leaf cherry
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc_pca_clr1 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca$CA$eig[1] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc_pca$CA$eig)
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc_pca_clr2 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca$CA$eig[2] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc_pca$CA$eig)
# Scale axes and plot ordination oak woodland
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl_pca_clr1 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca$CA$eig[1] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl_pca$CA$eig)
phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl_pca_clr2 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca$CA$eig[2] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl_pca$CA$eig)

# Plot PCA based on clr transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs_pca,
  type="samples",
  shape ="Family",
  color = "Sample.Site") +
  #facet_nested(~ Soil.Type+Habitat, scales = "free") +
  geom_point(aes(color = Sample.Site), size = 3, alpha = 0.7) +
  #ggrepel::geom_text_repel(aes(label = Genus.Abbr), size = 2, color = "black", max.overlaps = 30) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"),
                     name = "Site",
                     breaks = c("SK.01", "SK.02", "SK.03"),
                     labels = c("01", "02", "03")) +
  scale_shape_manual(values = c(15, 19, 3, 8, 4, 17,
                                1, 2, 3, 5, 25, 20))+
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
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ccs_pca_pca_plot_sample_site_family.pdf", width = 6, height = 4, dpi = 200)

# Plot PCA based on clr transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca,
  type="samples",
  color="Sample.Site",
  shape = "Family") +
  #facet_nested(~ Soil.Type+Habitat, scales = "free") +
  geom_point(aes(color = Sample.Site), size = 3, alpha = 0.7) +
  #ggrepel::geom_text_repel(aes(label = Genus.Abbr), size = 2, color = "black", max.overlaps = 30) +
  scale_color_manual(values = c("#A65628", "#F781BF", "#999999"),
                     name = "Site",
                     breaks = c("SK.16", "SK.17", "SK.18"),
                     labels = c("16", "17", "18")) +
  scale_shape_manual(values = c(15, 19, 3, 8, 4, 17,
                                1, 2, 3, 5, 25, 20))+
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
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca_pca_plot_sample_site_family.pdf", width = 6, height = 4, dpi = 200)

# Plot PCA based on clr transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr_pca,
  type="samples",
  color="Sample.Site",
  shape = "Family") +
  #facet_nested(~ Soil.Type+Habitat, scales = "free") +
  geom_point(aes(color = Sample.Site), size = 3, alpha = 0.7) +
  #ggrepel::geom_text_repel(aes(label = Genus.Abbr), size = 2, color = "black", max.overlaps = 30) +
  scale_color_manual(values = c("#984EA3", "#FF7F00", "#FFFF33"),
                     name = "Site",
                     breaks = c("SK.13", "SK.14", "SK.15"),
                     labels = c("13", "14", "15")) +
  scale_shape_manual(values = c(15, 19, 3, 8, 4, 17,
                                1, 2, 3, 5, 25, 20))+
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
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr_pca_pca_plot_sample_site_family.pdf", width = 6, height = 4, dpi = 200)

# Plot PCA based on clr transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip_pca,
  type="samples",
  shape = "Family",
  color="Sample.Site") +
  #facet_nested(~ Soil.Type+Habitat, scales = "free") +
  geom_point(aes(color = Sample.Site), size = 3, alpha = 0.7) +
  #ggrepel::geom_text_repel(aes(label = Genus.Abbr), size = 2, color = "black", max.overlaps = 30) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3"),
                     name = "Site",
                     breaks = c("SK.04", "SK.05", "SK.06"),
                     labels = c("04", "05", "06")) +
  scale_shape_manual(values = c(15, 19, 3, 8, 4, 17,
                                1, 2, 3, 5, 25, 20))+
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
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip_pca_pca_plot_sample_site_family.pdf", width = 6, height = 4, dpi = 200)


# Plot PCA based on clr transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc_pca,
  type="samples",
  color="Sample.Site",
  shape = "Family") +
  #facet_nested(~ Soil.Type+Habitat, scales = "free") +
  geom_point(aes(color = Sample.Site), size = 3, alpha = 0.7) +
  #ggrepel::geom_text_repel(aes(label = Genus.Abbr), size = 2, color = "black", max.overlaps = 30) +
  scale_color_manual(values = c("#E7298A", "#66A61E", "#E6AB02"),
                     name = "Site",
                     breaks = c("SK.07", "SK.08", "SK.09"),
                     labels = c("07", "08", "09")) +
  scale_shape_manual(values = c(15, 19, 3, 8, 4, 17,
                                1, 2, 3, 5, 25, 20))+
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
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc_pca_pca_plot_sample_site_family.pdf", width = 6, height = 4, dpi = 200)

# Plot PCA based on clr transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl_pca,
  type="samples",
  color="Sample.Site",
  shape = "Family") +
  #facet_nested(~ Soil.Type+Habitat, scales = "free") +
  geom_point(aes(color = Sample.Site), size = 3, alpha = 0.7) +
  #ggrepel::geom_text_repel(aes(label = Genus.Abbr), size = 2, color = "black", max.overlaps = 30) +
  scale_color_manual(values = c("#A6761D", "#666666", "lightblue"),
                     name = "Site",
                     breaks = c("SK.10", "SK.11", "SK.12"),
                     labels = c("10", "11", "12")) +
  scale_shape_manual(values = c(15, 19, 3, 8, 4, 17,
                                1, 2, 3, 5, 25, 20))+
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
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl_pca_pca_plot_sample_site_family.pdf", width = 6, height = 4, dpi = 200)

######################### scale phyloseq object to a minimum of 100 reads per sample
# Check if there are samples with zero reads and remove them
smin <- min(sample_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt)) # smin = 2.
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale <- prune_samples(sample_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt)>=100, phy_16_prune_nem_noncontam_prev05_true_all_filt)
smin <- min(sample_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale))
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale # 3057 taxa and 520 samples.

any(taxa_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale) == 0) # gives the number of cases
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale <-prune_taxa(taxa_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale) > 0, phy_16_prune_nem_noncontam_prev05_true_all_filt_scale) # two taxa removed
any(taxa_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale) == 0) # gives the number of cases
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale # 3055 taxa and 520 samples.

#########################  Save scaled phyloseq object as RDS #########################  
saveRDS(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, "data/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale.RDS")  

# Import phyloseq
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale <- readRDS("data/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale.RDS")
class(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale)

#########################  
# Export and Transpose matrix for of scaled dataset
# ASV 16S nema table
ASV_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_scale = as(otu_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale), "matrix")
write.csv(ASV_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, "exported_tables/nema/ASV_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_scale.csv")

# Extract metadata from phyloseq object
# 16S nema sample table
TABLE_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_scale = as(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale), "matrix")
class(TABLE_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_scale)
TABLE_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_df <- as.data.frame(TABLE_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_scale)
class(TABLE_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_df)
write.csv(TABLE_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_df, "exported_tables/nema/TABLE_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_df.csv")

# 16S nema taxonomy
TAX_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_scale  = as(tax_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale), "matrix")
TAX_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_scale <- as.data.frame(TAX_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_scale)
write.csv(TAX_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, "exported_tables/nema/TAX_SHSK_phy_16_prune_nem_noncontam_prev05_true_all_filt_scale.csv")

#########################  
# Apply different transformations to ASV table, e.g., converting number of reads to relative abundance
# Here we use the library microbiome
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp <- microbiome::transform(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, "compositional")
# phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_comp = transform_sample_counts(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, function(x) 100 * x/sum(x)). # This command is the equivalent to the microbiom::transform function "compositional"
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr <- microbiome::transform(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, "clr")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_log10 <- microbiome::transform(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, "log10p")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_hellinger <- microbiome::transform(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, "hellinger")
# This command is the equivalent to the microbiom::transform function "hellinger"
# phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_hellinger <- transform_sample_counts(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, function(x) sqrt(x/sum(x)))

############### Ordination using PCA and clr transformed phyloseq objects with scale
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_pca <- phyloseq::ordinate(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr,
  method = "RDA",
  distance = "euclidean"
)

# Plot scree plot to see PCA contributions
phyloseq::plot_scree(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_pca.pdf", width = 12, height = 6, dpi = 150)

head(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_pca$CA$eig)
sapply(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_pca$CA$eig[1:5], function(x) x / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_pca$CA$eig))

# Scale axes and plot ordination
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_pca_clr1 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_pca$CA$eig[1] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_pca$CA$eig)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_pca_clr2 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_pca$CA$eig[2] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_pca$CA$eig)

# Plot PCA based on clr transformation
soil.labs <- c("Compacted", "Uncompacted")
names(soil.labs) <- c("Compact.Soil", "Not.Compact.Soil")

habitat.labs <- c("CHA", "CSS", "NGR", "HLC", "OWL", "RIP")
names(habitat.labs) <- c("Chaparral", "Coastal.scrub", "Native.grass", "Hollyleaf.cherry", "Oak.wood", "Riparian")

# change shape variable (Infraorder or Feeding.Group) depending on desired graphic
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_pca,
  type="samples",
  color="Sample.Site",
  shape = "Infraorder") +
  facet_nested(~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs))+
  geom_point(aes(color = Sample.Site), alpha = 0.7, size = 3) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A",
                                         "#984EA3", "#FF7F00", "#FFFF33",
                                         "#A65628", "#F781BF", "#999999",
                                         "#1B9E77", "#D95F02", "#7570B3",
                                         "#E7298A", "#66A61E", "#E6AB02",
                                         "#A6761D", "#666666", "lightblue"),
                                         name = "Site",
                     breaks = c("SK.01", "SK.02", "SK.03",
                                "SK.13", "SK.14", "SK.15",
                                "SK.16", "SK.17", "SK.18",
                                "SK.04", "SK.05", "SK.06",
                                "SK.07", "SK.08", "SK.09",
                                "SK.10", "SK.11", "SK.12"),
                     labels = c("01", "02", "03",
                                "13", "14", "15",
                                "16", "17", "18",
                                "04", "05", "06",
                                "07", "08", "09",
                                "10", "11", "12"))+
  scale_shape_manual(values = c(15, 2, 3, 4, 8, 16, 17, 18, 9, 10, 1, 12),
                     name = "Family")+
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
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5, legend.position = "bottom") + # adjusts the title of the legend
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold"))
#ylab("NMDS2") + # add the title on y axis
#xlab("NMDS1") + # add the title on x axis
#ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
#theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
#annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
#annotate("text", x = 3, y = 2, label ="2D Stress: 0.07") +
#coord_fixed(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca_clr2 / phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_pca_clr1) +
ggsave("results/nema/16S_nem_ShSk_ASV_noncontam_prev05_true_all_filt_clr_no-mito-chloro_pca_plot_sample_site_infraorder_scale.pdf", width = 14, height = 6, dpi = 300)
ggsave("results/nema/16S_nem_ShSk_ASV_noncontam_prev05_true_all_filt_clr_no-mito-chloro_pca_plot_sample_site_feeding.pdf", width = 14, height = 6, dpi = 300)
#########################

########### Define phyloseq objects per habitat on scaled dataset
####### Coastal scrub

phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr, Habitat == "Coastal.scrub") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs)$Habitat)

phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs,
                                                                          sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs) != "Nem.391") # Tylenchoidea, SK01
                                                                            

####### Chaparral
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr, Habitat == "Chaparral") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha)$Habitat)

####### Native grass
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ngr <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr, Habitat == "Native.grass") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ngr)$Habitat)

####### Riparian
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_rip <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr, Habitat == "Riparian") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_rip)$Habitat)

####### Hollyleaf.cherry
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_hlc <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr, Habitat == "Hollyleaf.cherry") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_hlc)$Habitat)

####### Oakwood.land
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_owl <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr, Habitat == "Oak.wood") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
unique(sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_owl)$Habitat)

phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl <- subset_samples(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl,
                                                                          sample_names(phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl) != "Nem.479") # Dorylaimidae, SK11


####### PCA analysis for each phyloseq object
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs_pca <- phyloseq::ordinate(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs,
  method = "RDA",
  distance = "euclidean"
)

phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca <- phyloseq::ordinate(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha,
  method = "RDA",
  distance = "euclidean"
)

phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ngr_pca <- phyloseq::ordinate(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ngr,
  method = "RDA",
  distance = "euclidean"
)


phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_rip_pca <- phyloseq::ordinate(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_rip,
  method = "RDA",
  distance = "euclidean"
)

phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_hlc_pca <- phyloseq::ordinate(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_hlc,
  method = "RDA",
  distance = "euclidean"
)

phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_owl_pca <- phyloseq::ordinate(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_owl,
  method = "RDA",
  distance = "euclidean"
)


# Plot scree plot to see PCA contributions coastal scrub
phyloseq::plot_scree(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs_pca.pdf", width = 12, height = 6, dpi = 150)
head(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs_pca$CA$eig)
sapply(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs_pca$CA$eig[1:5], function(x) x / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs_pca$CA$eig))

# Plot scree plot to see PCA contributions chaparral
phyloseq::plot_scree(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca.pdf", width = 12, height = 6, dpi = 150)
head(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca$CA$eig)
sapply(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca$CA$eig[1:5], function(x) x / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca$CA$eig))

# Plot scree plot to see PCA contributions native grass
phyloseq::plot_scree(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ngr_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ngr_pca.pdf", width = 12, height = 6, dpi = 150)
head(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ngr_pca$CA$eig)
sapply(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ngr_pca$CA$eig[1:5], function(x) x / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ngr_pca$CA$eig))

# Plot scree plot to see PCA contributions ripirian
phyloseq::plot_scree(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_rip_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_rip_pca.pdf", width = 12, height = 6, dpi = 150)
head(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_rip_pca$CA$eig)
sapply(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_rip_pca$CA$eig[1:5], function(x) x / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_rip_pca$CA$eig))

# Plot scree plot to see PCA contributions holy-leaf cherry
phyloseq::plot_scree(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_hlc_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_hlc_pca.pdf", width = 12, height = 6, dpi = 150)
head(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_hlc_pca$CA$eig)
sapply(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_hlc_pca$CA$eig[1:5], function(x) x / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_hlc_pca$CA$eig))

# Plot scree plot to see PCA contributions oak woodland
phyloseq::plot_scree(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_owl_pca) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_owl_pca.pdf", width = 12, height = 6, dpi = 150)
head(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_owl_pca$CA$eig)
sapply(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_owl_pca$CA$eig[1:5], function(x) x / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_owl_pca$CA$eig))

# Scale axes and plot ordination coastal scrub
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs_pca_clr1 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs_pca$CA$eig[1] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs_pca$CA$eig)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs_pca_clr2 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs_pca$CA$eig[2] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs_pca$CA$eig)
# Scale axes and plot ordination coastal chaparral
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca_clr1 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca$CA$eig[1] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca$CA$eig)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca_clr2 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca$CA$eig[2] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca$CA$eig)
# Scale axes and plot ordination native grass
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ngr_pca_clr1 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca$CA$eig[1] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ngr_pca$CA$eig)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ngr_pca_clr2 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca$CA$eig[2] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ngr_pca$CA$eig)
# Scale axes and plot ordination coastal ripirian
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_rip_pca_clr1 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca$CA$eig[1] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_rip_pca$CA$eig)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_rip_pca_clr2 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca$CA$eig[2] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_rip_pca$CA$eig)
# Scale axes and plot ordination holy-leaf cherry
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_hlc_pca_clr1 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca$CA$eig[1] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_hlc_pca$CA$eig)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_hlc_pca_clr2 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca$CA$eig[2] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_hlc_pca$CA$eig)
# Scale axes and plot ordination oak woodland
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_owl_pca_clr1 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca$CA$eig[1] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_owl_pca$CA$eig)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_owl_pca_clr2 <- phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca$CA$eig[2] / sum(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_owl_pca$CA$eig)

# Plot PCA based on clr transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs_pca,
  type="samples",
  shape ="Family",
  color = "Sample.Site") +
  #facet_nested(~ Soil.Type+Habitat, scales = "free") +
  geom_point(aes(color = Sample.Site), size = 3, alpha = 0.7) +
  #ggrepel::geom_text_repel(aes(label = Genus.Abbr), size = 2, color = "black", max.overlaps = 30) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"),
                     name = "Site",
                     breaks = c("SK.01", "SK.02", "SK.03"),
                     labels = c("01", "02", "03")) +
  scale_shape_manual(values = c(15, 19, 3, 8, 4, 17,
                                1, 2, 3, 5, 25, 20))+
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
#coord_fixed(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_pca_clr2 / phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_pca_clr1) +
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ccs_pca_pca_plot_sample_site_family.pdf", width = 6, height = 4, dpi = 200)

# Plot PCA based on clr transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_cha_pca,
  type="samples",
  color="Sample.Site",
  shape = "Family") +
  #facet_nested(~ Soil.Type+Habitat, scales = "free") +
  geom_point(aes(color = Sample.Site), size = 3, alpha = 0.7) +
  #ggrepel::geom_text_repel(aes(label = Genus.Abbr), size = 2, color = "black", max.overlaps = 30) +
  scale_color_manual(values = c("#A65628", "#F781BF", "#999999"),
                     name = "Site",
                     breaks = c("SK.16", "SK.17", "SK.18"),
                     labels = c("16", "17", "18")) +
  scale_shape_manual(values = c(15, 19, 3, 8, 4, 17,
                                1, 2, 3, 5, 25, 20))+
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
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_cha_pca_pca_plot_sample_site_family.pdf", width = 6, height = 4, dpi = 200)

# Plot PCA based on clr transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_ngr_pca,
  type="samples",
  color="Sample.Site",
  shape = "Family") +
  #facet_nested(~ Soil.Type+Habitat, scales = "free") +
  geom_point(aes(color = Sample.Site), size = 3, alpha = 0.7) +
  #ggrepel::geom_text_repel(aes(label = Genus.Abbr), size = 2, color = "black", max.overlaps = 30) +
  scale_color_manual(values = c("#984EA3", "#FF7F00", "#FFFF33"),
                     name = "Site",
                     breaks = c("SK.13", "SK.14", "SK.15"),
                     labels = c("13", "14", "15")) +
  scale_shape_manual(values = c(15, 19, 3, 8, 4, 17,
                                1, 2, 3, 5, 25, 20))+
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
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_ngr_pca_pca_plot_sample_site_family.pdf", width = 6, height = 4, dpi = 200)

# Plot PCA based on clr transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_rip_pca,
  type="samples",
  shape = "Family",
  color="Sample.Site") +
  #facet_nested(~ Soil.Type+Habitat, scales = "free") +
  geom_point(aes(color = Sample.Site), size = 3, alpha = 0.7) +
  #ggrepel::geom_text_repel(aes(label = Genus.Abbr), size = 2, color = "black", max.overlaps = 30) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3"),
                     name = "Site",
                     breaks = c("SK.04", "SK.05", "SK.06"),
                     labels = c("04", "05", "06")) +
  scale_shape_manual(values = c(15, 19, 3, 8, 4, 17,
                                1, 2, 3, 5, 25, 20))+
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
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_rip_pca_pca_plot_sample_site_family.pdf", width = 6, height = 4, dpi = 200)


# Plot PCA based on clr transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_hlc_pca,
  type="samples",
  color="Sample.Site",
  shape = "Family") +
  #facet_nested(~ Soil.Type+Habitat, scales = "free") +
  geom_point(aes(color = Sample.Site), size = 3, alpha = 0.7) +
  #ggrepel::geom_text_repel(aes(label = Genus.Abbr), size = 2, color = "black", max.overlaps = 30) +
  scale_color_manual(values = c("#E7298A", "#66A61E", "#E6AB02"),
                     name = "Site",
                     breaks = c("SK.07", "SK.08", "SK.09"),
                     labels = c("07", "08", "09")) +
  scale_shape_manual(values = c(15, 19, 3, 8, 4, 17,
                                1, 2, 3, 5, 25, 20))+
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
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_hlc_pca_pca_plot_sample_site_family.pdf", width = 6, height = 4, dpi = 200)

# Plot PCA based on clr transformation
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl,
  ordination = phy_16_prune_nem_noncontam_prev05_true_all_filt_clr_owl_pca,
  type="samples",
  color="Sample.Site",
  shape = "Family") +
  #facet_nested(~ Soil.Type+Habitat, scales = "free") +
  geom_point(aes(color = Sample.Site), size = 3, alpha = 0.7) +
  #ggrepel::geom_text_repel(aes(label = Genus.Abbr), size = 2, color = "black", max.overlaps = 30) +
  scale_color_manual(values = c("#A6761D", "#666666", "lightblue"),
                     name = "Site",
                     breaks = c("SK.10", "SK.11", "SK.12"),
                     labels = c("10", "11", "12")) +
  scale_shape_manual(values = c(15, 19, 3, 8, 4, 17,
                                1, 2, 3, 5, 25, 20))+
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
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_clr_owl_pca_pca_plot_sample_site_family.pdf", width = 6, height = 4, dpi = 200)
#########################
# Plot CAP analyses using Primer coordinates
cap_all <- read_excel("data/nema-16s-cap-coord-pa.xlsx", sheet = "all-cap-pa")
cap_owl <- read_excel("data/nema-16s-cap-coord-pa.xlsx", sheet = "owl-cap-pa")
cap_rip <- read_excel("data/nema-16s-cap-coord-pa.xlsx", sheet = "rip-cap-pa")
cap_hlc <- read_excel("data/nema-16s-cap-coord-pa.xlsx", sheet = "hlc-cap-pa")
cap_ngr <- read_excel("data/nema-16s-cap-coord-pa.xlsx", sheet = "ngr-cap-pa")
cap_cha <- read_excel("data/nema-16s-cap-coord-pa.xlsx", sheet = "cha-cap-pa")
cap_css <- read_excel("data/nema-16s-cap-coord-pa.xlsx", sheet = "css-cap-pa")

# Plot PCA based on clr transformation
soil.labs <- c("Compacted", "Uncompacted")
names(soil.labs) <- c("Compact.Soil", "Not.Compact.Soil")

habitat.labs <- c("CHA", "CSS", "NGR", "HLC", "OWL", "RIP")
names(habitat.labs) <- c("Chaparral", "Coastal.scrub", "Native.grass", "Hollyleaf.cherry", "Oak.wood", "Riparian")

# change shape variable (Infraorder or Feeding.Group) depending on desired graphic
theme_set(theme_bw())
cap_all_tax_plot <- ggplot(cap_all) +
  geom_point(aes(x = CAP1, y = CAP2, color=Sample.Site, shape = Infraorder), size=3) +
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  theme(legend.title = element_text(face = "bold", size = 10))+
  theme(legend.title.align = 0.5) +
  scale_shape_manual(values=c(15, 1, 17, 5, 18, 10, 19, 2, 4, 8, 3, 16))+
  scale_color_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A","#A65628", "#F781BF", "#999999", "#984EA3", "#FF7F00", "#FFFF33",
                              "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "lightblue"))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())+
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5, legend.position = "bottom") + # adjusts the title of the legend
  theme(strip.background =element_rect(color = "black", fill = "white", linewidth = 1, linetype = "solid"),
        strip.text = element_text(size = 12, color = "black", face = "bold")) # Format facet grid title
cap_all_tax_plot
ggsave("results/nema/cap_all_tax_pa_plot.pdf", width = 16, height = 6, dpi = 150) # save graphic

cap_all_feeding_plot <- ggplot(cap_all) +
  geom_point(aes(x = CAP1, y = CAP2, color=Sample.Site, shape = Feeding.Group), size=3) +
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.labs, Habitat = habitat.labs)) +
  theme(legend.title = element_text(face = "bold", size = 10))+
  theme(legend.title.align = 0.5) +
  scale_shape_manual(values=c(15, 1, 17, 5, 18, 10, 19, 2, 4, 8, 3, 16))+
  scale_color_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A","#A65628", "#F781BF", "#999999", "#984EA3", "#FF7F00", "#FFFF33",
                              "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "lightblue"))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())+
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5, legend.position = "bottom") + # adjusts the title of the legend
theme(strip.background =element_rect(color = "black", fill = "white", linewidth = 1, linetype = "solid"),
      strip.text = element_text(size = 12, color = "black", face = "bold")) # Format facet grid title
cap_all_feeding_plot
ggsave("results/nema/cap_all_feeding_pa_plot.pdf", width = 16, height = 6, dpi = 150) # save graphic

cap_owl_plot <- ggplot(cap_owl) +
  geom_point(aes(x = CAP1, y = CAP2, color=Sample.Site, shape = Infraorder), size=3) +
  theme(legend.position = "right") +
  theme(legend.title = element_text(face = "bold", size = 10))+
  theme(legend.title.align = 0.5) +
  scale_shape_manual(values=c(15, 1, 17, 5, 18, 10, 19, 2, 4))+
  scale_color_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A"))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
cap_owl_plot
ggsave("results/nema/cap_owl_plot_infraorder_pa.pdf", width = 7, height = 5, dpi = 150) # save graphic

cap_rip_plot <- ggplot(cap_rip) +
  geom_point(aes(x = CAP1, y = CAP2, color=Sample.Site, shape = Infraorder), size=3) +
  theme(legend.position = "right") +
  theme(legend.title = element_text(face = "bold", size = 10))+
  theme(legend.title.align = 0.5) +
  scale_shape_manual(values=c(15, 1, 17, 5, 18, 10, 19, 2, 4, 11, 3, 6))+
  scale_color_manual(values=c("#A65628", "#F781BF", "#999999"))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
cap_rip_plot
ggsave("results/nema/cap_rip_plot_infraorder_pa.pdf", width = 7, height = 5, dpi = 150) # save graphic

cap_hlc_plot <- ggplot(cap_hlc) +
  geom_point(aes(x = CAP1, y = CAP2, color=Sample.Site, shape = Infraorder), size=3) +
  theme(legend.position = "right") +
  theme(legend.title = element_text(face = "bold", size = 10))+
  theme(legend.title.align = 0.5) +
  scale_shape_manual(values=c(15, 1, 17, 5, 18, 10, 19, 2, 4))+
  scale_color_manual(values=c("#984EA3", "#FF7F00", "#FFFF33"))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
cap_hlc_plot
ggsave("results/nema/cap_hlc_plot_infraorder_pa.pdf", width = 7, height = 5, dpi = 150) # save graphic

cap_ngr_plot <- ggplot(cap_ngr) +
  geom_point(aes(x = CAP1, y = CAP2, color=Sample.Site, shape = Infraorder), size=3) +
  theme(legend.position = "right") +
  theme(legend.title = element_text(face = "bold", size = 10))+
  theme(legend.title.align = 0.5) +
  scale_shape_manual(values=c(15, 1, 17, 5, 18, 10, 19, 2, 4, 11))+
  scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3"))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
cap_ngr_plot
ggsave("results/nema/cap_ngr_plot_infraorder_pa.pdf", width = 7, height = 5, dpi = 150) # save graphic

cap_cha_plot <- ggplot(cap_cha) +
  geom_point(aes(x = CAP1, y = CAP2, color=Sample.Site, shape = Infraorder), size=3) +
  theme(legend.position = "right") +
  theme(legend.title = element_text(face = "bold", size = 10))+
  theme(legend.title.align = 0.5) +
  scale_shape_manual(values=c(15, 1, 17, 5, 18, 10, 19, 2, 4, 11, 3, 14))+
  scale_color_manual(values=c("#E7298A", "#66A61E", "#E6AB02"))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
cap_cha_plot
ggsave("results/nema/cap_cha_plot_infraorder_pa.pdf", width = 7, height = 5, dpi = 150) # save graphic

cap_css_plot <- ggplot(cap_css) +
  geom_point(aes(x = CAP1, y = CAP2, color=Sample.Site, shape = Infraorder), size=3) +
  theme(legend.position = "right") +
  theme(legend.title = element_text(face = "bold", size = 10))+
  theme(legend.title.align = 0.5) +
  scale_shape_manual(values=c(15, 1, 17, 5, 18, 10, 19, 2, 4, 11, 3, 14))+
  scale_color_manual(values=c("#A6761D", "#666666", "lightblue"))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
cap_css_plot
ggsave("results/nema/cap_css_plot_infraorder_pa.pdf", width = 7, height = 5, dpi = 150) # save graphic
######################### Generate a data.frame with alpha diversity measures #########################
# Calculate alpha-diversity measures (For plot purposes only!)
# This can be done using the different phyloseq alpha diversity measures
# You will get a Warning message for each index since there is no singletons on the dataset
alpha_div_16S_nem <- data.frame(
  "Observed" = phyloseq::estimate_richness(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, measures = "Shannon"),
  "InvSimpson" = phyloseq::estimate_richness(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, measures = "InvSimpson"),
  "Reads" = phyloseq::sample_sums(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale),
  "Soil.Type" = phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale)$Soil.Type,
  "Habitat" = phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale)$Habitat,
  "Sample.Site" = phyloseq::sample_data(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale)$Sample.Site)
alpha_div_16S_nem$Evenness <- alpha_div_16S_nem$Shannon/log(alpha_div_16S_nem$Observed)
head(alpha_div_16S_nem)

# Rename variable InvSimpson to Simpson
# The function rename & %>% works on dplyr. make sure it is loaded.
alpha_div_16S_nem <- alpha_div_16S_nem %>%
  rename(Simpson = InvSimpson)
head(alpha_div_16S_nem)

# Reorder dataframe, first categorical then numerical variables  
alpha_div_16S_nem <- alpha_div_16S_nem[, c(7, 6, 5, 4, 1, 2, 3, 8)]
head(alpha_div_16S_nem)
write.csv(alpha_div_16S_nem, "exported_tables/nema/alpha_div_16S_nem.csv")

#Summarize alpha diversity measures by site
summary_alpha_16S_nem_site <- alpha_div_16S_nem %>%
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
write.csv(summary_alpha_16S_nem_site, "exported_tables/nema/summary_alpha_16S_nem_site.csv")

#Summarize alpha diversity measures by habitat
summary_alpha_16S_nem_habitat <- alpha_div_16S_nem %>%
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
write.csv(summary_alpha_16S_nem_habitat, "exported_tables/nema/summary_alpha_16S_nem_habitat.csv")

#Summarize alpha diversity measures by soil type
summary_alpha_16S_nem_soil_type <- alpha_div_16S_nem %>%
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
write.csv(summary_alpha_16S_nem_soil_type, "exported_tables/nema/summary_alpha_16S_nem_soil_type.csv")


# KW analysis alpha_16S_nem on all metrics at once
# Remember, numerical variables are from columns 4-8. The test is by soil type, column 3
kw_alpha_16S_nem_soil_type_univ <- as.data.frame(sapply(4:8, function(x) kruskal.test(alpha_div_16S_nem[,x],
                                                                                  alpha_div_16S_nem[,3])))

# Rename columns with the proper variable names
kw_alpha_16S_nem_soil_type_univ <- kw_alpha_16S_nem_soil_type_univ %>%
  rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)

# Make some adjustments on the table
kw_alpha_16S_nem_soil_type_univ <- t(kw_alpha_16S_nem_soil_type_univ) # transpose
kw_alpha_16S_nem_soil_type_univ <- as_tibble(kw_alpha_16S_nem_soil_type_univ, rownames = "Metric") # adding rownames as a column
kw_alpha_16S_nem_soil_type_univ <- kw_alpha_16S_nem_soil_type_univ[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_16S_nem_soil_type_univ) # checking object class

# Save resulting table with fwrite to avoid any issues with characters
data.table::fwrite(kw_alpha_16S_nem_soil_type_univ, "exported_tables/nema/kw_alpha_16S_nem_soil_type_univ.csv")

# KW analysis alpha_16S_nem on all metrics at once
# Remember, numerical variables are from columns 4-8. The test is by habitat, column 2
kw_alpha_16S_nema_habitat_univ <- as.data.frame(sapply(4:8, function(x) kruskal.test(alpha_div_16S_nem[,x],
                                                                                alpha_div_16S_nem[,2])))

kw_alpha_16S_nema_habitat_univ <- kw_alpha_16S_nema_habitat_univ %>%
  rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)

# Make some adjustments on the table
kw_alpha_16S_nema_habitat_univ <- t(kw_alpha_16S_nema_habitat_univ) # transpose
kw_alpha_16S_nema_habitat_univ <- as_tibble(kw_alpha_16S_nema_habitat_univ, rownames = "Metric") # adding rownames as a column
kw_alpha_16S_nema_habitat_univ <- kw_alpha_16S_nema_habitat_univ[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_16S_nema_habitat_univ) # checking object class

# Save resulting table with fwrite to avoid and issues with characters
data.table::fwrite(kw_alpha_16S_nema_habitat_univ, "exported_tables/nema/kw_alpha_16S_nema_habitat_univ.csv")

# pairwise comparison including p-value adjustment
kw_alpha_16S_nema_habitat_univ_pwt <- as.data.frame(sapply(4:8, function(x) pairwise.t.test(alpha_div_16S_nem[,x],
                                                                                    alpha_div_16S_nem[,2], p.adjust.method = "BH")))
kw_alpha_16S_nema_habitat_univ_pwt <- kw_alpha_16S_nema_habitat_univ_pwt %>%
  rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)

# Make some adjustments on the table
kw_alpha_16S_nema_habitat_univ_pwt <- t(kw_alpha_16S_nema_habitat_univ_pwt) # transpose
kw_alpha_16S_nema_habitat_univ_pwt <- as_tibble(kw_alpha_16S_nema_habitat_univ_pwt, rownames = "Metric") # adding rownames as a column
class(kw_alpha_16S_nema_habitat_univ_pwt) # checking object class

# Save resulting table with fwrite to avoid and issues with characters
fwrite(kw_alpha_16S_nema_habitat_univ_pwt, "exported_tables/nema/kw_alpha_16S_nema_habitat_univ_pwt.csv")

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
alpha_div_16S_nem %>%
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
ggsave("results/nema/16S_nem_alpha_diversity_phy_nem_16S_prune_noncontam_prev05_true_scale_site.pdf", width = 8, height = 6, dpi = 150) # save graphic

# Plot alpha-diversity per habitat within soil type
alpha_div_16S_nem %>%
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
ggsave("results/nema/16S_nem_alpha_diversity_phy_nem_16S_prune_noncontam_prev05_true_scale_habitat.pdf", width = 6, height = 6, dpi = 150) # save graphic

# Plot alpha-diversity per soil type
alpha_div_16S_nem$Soil.Type <- factor(alpha_div_16S_nem$Soil.Type,
                                       level = c("Compact.Soil", "Not.Compact.Soil"),
                                       labels = c("Compacted", "Uncompacted"))

alpha_color_nem_soil_type <- c("#52392F", "#397A4C")

alpha_div_16S_nem %>%
  gather(key = metric, value = value, c("Reads", "Observed", "Shannon", "Simpson", "Evenness")) %>%
  mutate(metric = factor(metric, levels = c("Reads", "Observed", "Shannon", "Simpson", "Evenness"))) %>%
  ggplot(aes(x = Soil.Type, y = value)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  geom_jitter(aes(color = Soil.Type), height = 0, width = .2) +
  facet_grid(metric ~ Soil.Type, scales = "free") +
  scale_color_manual(values = alpha_color_nem_soil_type) +
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
ggsave("results/nema/16S_nem_alpha_diversity_phy_nem_16S_prune_noncontam_prev05_true_scale_nem_soil_type.pdf", width = 5, height = 6, dpi = 150) # save graphic

############################# Plotting top 20 most abundant at different taxonomic ranks ###############################
# This function can be used to split the data into the top N and the less abundant taxa
# Use scaled (samples >100 reads) phyloseq object
# Check how many phyla are in the phyloseq object
get_taxa_unique(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, taxonomic.rank = "Phylum") # a total of 30 phyla
# Check how many taxonomic levels
head(tax_table(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale)) # 6 levels

# Function to collapse a certain number of taxa into category others
merge_less_than_top_prev05_top20 <- function(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, top=19){
  transformed <- transform_sample_counts(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, function(x) x/sum(x))
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

# Phylum level top 20
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_phylum <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, "Phylum")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_phylum_top20 <- merge_less_than_top_prev05_top20(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_phylum, top=19)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_phylum_top20_df = psmelt(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_phylum_top20) # Melt to long format
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_phylum_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Phylum, data=phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_phylum_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_phylum_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Phylum, data=phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_phylum_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_phylum_top20_agr_soil_type = aggregate(Abundance~Soil.Type+Phylum, data=phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_phylum_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_phylum_top20_agr_site$Phylum)

# Class level top 20
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_class <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, "Class")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_class_top20 <- merge_less_than_top_prev05_top20(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_class, top=19)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_class_top20_df = psmelt(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_class_top20) # Melt to long format
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_class_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Class, data=phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_class_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_class_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Class, data=phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_class_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_class_top20_agr_soil_type = aggregate(Abundance~Soil.Type+Class, data=phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_class_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_class_top20_agr_site$Class)

# Order level top 20
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_order <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, "Order")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_order_top20 <- merge_less_than_top_prev05_top20(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_order, top=19)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_order_top20_df = psmelt(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_order_top20) # Melt to long format
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_order_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Order, data=phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_order_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_order_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Order, data=phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_order_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_order_top20_agr_soil_type = aggregate(Abundance~Soil.Type+Order, data=phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_order_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_order_top20_agr_site$Order)

# Family level top 20
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_fam <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, "Family")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_fam_top20 <- merge_less_than_top_prev05_top20(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_fam, top=19)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_fam_top20_df = psmelt(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_fam_top20) # Melt to long format
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_fam_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Family, data=phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_fam_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_fam_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Family, data=phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_fam_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_fam_top20_agr_soil_type = aggregate(Abundance~Soil.Type+Family, data=phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_fam_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_fam_top20_agr_site$Family)

# Genus level top 20
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_genus <- tax_glom(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale, "Genus")
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_genus_top20 <- merge_less_than_top_prev05_top20(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_genus, top=19)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_genus_top20_df = psmelt(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_genus_top20) # Melt to long format
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_genus_top20_agr_site = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+Genus, data=phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_genus_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_genus_top20_agr_hab = aggregate(Abundance~Soil.Type+Habitat+Genus, data=phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_genus_top20_df, FUN=mean) # add common factors to use for plotting
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_genus_top20_agr_soil_type = aggregate(Abundance~Soil.Type+Genus, data=phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_genus_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_genus_top20_agr_site$Genus)

############### Organize color scales and factors to plot top12 and top20 taxa as barplots
# List of 20 distinct colors
colors_top20 <- c("#9a6324", "#46f0f0", "#1F78B4", "#aaffc3", "#000075",  "#33A02C", "#4363d8", "#008080", "#FB9A99", "#e6194b",
                           "#f032e6", "#bcf60c", "#fabebe","#f58231", "#ffe119", "#6A3D9A", "#e6beff", "#800000", "#fffac8", "gray")
                           
# List of 12 distinct colors
colors_top12 <- c("#9a6324", "#46f0f0", "#aaffc3", "#33A02C", "#4363d8", "#f032e6",
                           "#bcf60c", "#f58231", "#6A3D9A", "#e6beff", "#fffac8", "gray")
                           

# Create new labels for important factors
soil.labs <- c("Compacted", "Uncompacted")
names(soil.labs) <- c("Compact.Soil", "Not.Compact.Soil")

habitat.labs <- c("CHA", "CSS", "NGR", "HLC",  "OWL", "RIP")
names(habitat.labs) <- c("Chaparral", "Coastal.scrub", "Native.grass",
                         "Hollyleaf.cherry", "Oak.wood",  "Riparian")

# Put "Others" to the final of the Phylum list - top 20
unique(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_phylum_top20_agr_site$Phylum)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_phylum_top20_agr_site$Phylum <- factor(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_phylum_top20_agr_site$Phylum,
                                                                                              levels = c("Acidobacteriota", "Actinobacteriota", "Bacteroidota",
                                                                                                         "Bdellovibrionota", "Campylobacterota", "Chloroflexi",
                                                                                                         "Crenarchaeota", "Cyanobacteria", "Deinococcota",
                                                                                                         "Dependentiae", "Firmicutes", "Fusobacteriota", 
                                                                                                         "Gemmatimonadota", "Latescibacterota", "MBNT15",
                                                                                                          "Myxococcota", "Planctomycetota", "Proteobacteria",
                                                                                                         "Verrucomicrobiota", "Others"))
# Plot by site - Phylum level top 20
ggplot(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_phylum_top20_agr_site, aes(x = Sample.Site, y = Abundance, fill = Phylum)) +
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
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_phylum_top20_agr_site.pdf", width = 9, height = 5, dpi = 200) # save graphic

# Put "Others" to the final of the Class list - top 20
unique(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_class_top20_agr_site$Class)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_class_top20_agr_site$Class <- factor(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_class_top20_agr_site$Class,
                                                                                            levels = c("Acidobacteriae", "Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidia",
                                                                                                       "Blastocatellia", "Chlamydiae", "Clostridia", "Deinococci", "Gammaproteobacteria",
                                                                                                       "Gemmatimonadetes", "Longimicrobia", "Planctomycetes", "Polyangia",
                                                                                                       "Thermoanaerobacteria", "Thermoleophilia", "Vampirivibrionia",
                                                                                                       "Verrucomicrobiae", "Vicinamibacteria", "Others"))
# Plot by site - Class level top20
ggplot(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_class_top20_agr_site, aes(x = Sample.Site, y = Abundance, fill = Class)) +
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
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_class_top20_agr_site.pdf", width = 10, height = 5, dpi = 200) # save graphic

# Put "Others" to the final of the order list - top 20
unique(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_order_top20_agr_site$Order)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_order_top20_agr_site$Order <- factor(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_order_top20_agr_site$Order,
                                                                                            levels = c("Bacillales", "Burkholderiales", "Caulobacterales", "Corynebacteriales", "Cytophagales", "Enterobacterales",
                                                                                                       "Flavobacteriales",  "Lactobacillales", "Micrococcales", "Obscuribacterales", "Opitutales", "Pseudomonadales",
                                                                                                       "Rhizobiales", "Rhodobacterales", "Rickettsiales","Sphingomonadales",
                                                                                                       "Staphylococcales", "Thermales", "Xanthomonadales", "Others"))
# Plot by site - Order level top20
ggplot(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_order_top20_agr_site, aes(x = Sample.Site, y = Abundance, fill = Order)) +
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
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_order_top20_agr_site.pdf", width = 10, height = 5, dpi = 200) # save graphic

# Put "Others" to the final of the family list - top 20
unique(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_fam_top20_agr_site$Family)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_fam_top20_agr_site$Family <- factor(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_fam_top20_agr_site$Family,
                                                                                              levels = c("Amoebophilaceae", "Bacillaceae", "Burkholderiaceae", "Caulobacteraceae", "Comamonadaceae", 
                                                                                                         "Corynebacteriaceae", "Enterococcaceae", "Moraxellaceae", "Mycobacteriaceae", "Obscuribacteraceae",
                                                                                                         "Opitutaceae", "Oxalobacteraceae", "Pasteurellaceae", "Pseudomonadaceae", "Rhodobacteraceae",
                                                                                                         "Sphingomonadaceae",  "Staphylococcaceae",  "Weeksellaceae", "Xanthobacteraceae", "Others"))
# Plot by site - Family level top20
ggplot(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_fam_top20_agr_site, aes(x = Sample.Site, y = Abundance, fill = Family)) +
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
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_fam_top20_agr_site.pdf", width = 10, height = 5, dpi = 200) # save graphic

# Put "Others" to the final of the genus list - top 20
unique(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_genus_top20_agr_site$Genus)
phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_genus_top20_agr_site$Genus <- factor(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_genus_top20_agr_site$Genus,
                                                                                            levels = c("Bacillus", "Candidatus Cardinium", "Cloacibacterium", "Corynebacterium", "Enhydrobacter",
                                                                                                       "Geobacillus", "Lacunisphaera", "Limnohabitans", "Massilia", "Mycobacterium",
                                                                                                       "Novosphingobium", "Phenylobacterium", "Pseudomonas", "Ralstonia", "Rhizobacter",
                                                                                                       "Sphingobium", "Sphingomonas", "Staphylococcus", "Tetragenococcus","Others"))
# Plot by site - Genus level top20
ggplot(phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_genus_top20_agr_site, aes(x = Sample.Site, y = Abundance, fill = Genus)) +
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
ggsave("results/nema/phy_16_prune_nem_noncontam_prev05_true_all_filt_scale_top_genus_top20_agr_site.pdf", width = 10, height = 5, dpi = 200) # save graphic