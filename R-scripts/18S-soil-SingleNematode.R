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

# Load libraries
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
#library("dplyr") #A tool for working/manipulating dataframes. It may create problems with plyr.
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

# Import asv table, taxonomy table, and sample table in excel format
asv_tab_18S <- read_excel("data/SHSK_18S_ASV_table_Nov_08_2021.xlsx")
tax_tab_18S <- read_excel("data/SHSK_18S_Taxonomy_table_Nov_02_2022.xlsx")
sample_tab_18S <- read_excel("data/SHSK_18S_Sample_table_Nov_08_2021.xlsx")

# Convert asv and taxonomy tables to matrix format; sample table as dataframe
asv_tab_18S <- as.matrix(asv_tab_18S)
class(asv_tab_18S)
tax_tab_18S <- as.matrix(tax_tab_18S)
sample_tab_18S <- as.data.frame(sample_tab_18S)

# Rename row names using ASV column for asv table and taxonomy table; rename row names using Sample column for sample table
rownames(asv_tab_18S) <- asv_tab_18S[,1]
rownames(tax_tab_18S) <- tax_tab_18S[,1]
rownames(sample_tab_18S) <- sample_tab_18S[,1]

# Exclude first column from all three tables; set asv table as numeric
asv_tab_18S <- asv_tab_18S[,-1]
class(asv_tab_18S) <- "numeric"
tax_tab_18S <- tax_tab_18S[,-1]
sample_tab_18S <- sample_tab_18S[,-1]

#Transform to phyloseq objects
ASV_18S = otu_table(asv_tab_18S, taxa_are_rows = TRUE)
TAX_18S = tax_table(tax_tab_18S)
samples_18S = sample_data(sample_tab_18S)

# Creating a phyloseq object
phy_18S <- phyloseq(ASV_18S, TAX_18S, samples_18S)
phy_18S # 5857 taxa and 625 samples
head(sample_data(phy_18S))

#Visualize some data entries from phyloseq object
sample_names(phy_18S)
rank_names(phy_18S)
sample_variables(phy_18S)

# Check if there are samples with zero reads and remove them
smin <- min(sample_sums(phy_18S)) #smin should be = 0
phy_18S_prune <- prune_samples(sample_sums(phy_18S)>0, phy_18S)
smin <- min(sample_sums(phy_18S_prune))
phy_18S_prune #5857 taxa and 599 samples

# Check for ASVs that have no counted reads
any(taxa_sums(phy_18S_prune) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_18S_prune) == 0) # gives the number of cases

#Inspect Library sizes
df_phy_18S_prune <- as.data.frame(sample_data(phy_18S_prune))
df_phy_18S_prune$LibrarySize <- sample_sums(phy_18S_prune)
df_phy_18S_prune <- df_phy_18S_prune[order(df_phy_18S_prune$LibrarySize),]
df_phy_18S_prune$Index <- seq(nrow(df_phy_18S_prune))
ggplot(data = df_phy_18S_prune, aes(x= Index, y= LibrarySize, color = Sample.Type)) +
  geom_point()
ggsave("results/all-dataset/18S_library-size-sample-type-all-samples.pdf", width = 8, height = 6, dpi = 150)

####################################### Only Soil samples #######################################
# Keep only samples from the soil dataset
phy_18S_prune_soil <- subset_samples(phy_18S_prune, Experiment == "Environment") %>%
  prune_samples(sample_sums(.) >0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_18S_prune_soil # 4262 taxa and 53 samples

# Remove samples with zero reads
smin <- min(sample_sums(phy_18S_prune_soil))

# Check for ASVs that have no counted reads
any(taxa_sums(phy_18S_prune_soil) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_18S_prune_soil) == 0) # gives the number of cases
head(sample_data(phy_18S_prune_soil))
head(otu_table(phy_18S_prune_soil))

#Inspect Library sizes
df_phy_18S_prune_soil <- as.data.frame(sample_data(phy_18S_prune_soil))
df_phy_18S_prune_soil$LibrarySize <- sample_sums(phy_18S_prune_soil)
df_phy_18S_prune_soil <- df_phy_18S_prune_soil[order(df_phy_18S_prune_soil$LibrarySize),]
df_phy_18S_prune_soil$Index <- seq(nrow(df_phy_18S_prune_soil))
ggplot(data = df_phy_18S_prune_soil, aes(x= Index, y= LibrarySize, color = Sample.Type)) +
  geom_point()
ggsave("results/soil/18S_Soil_library-size-sample-type-soil-samples.pdf", width = 8, height = 6, dpi = 150)

ggplot(data = df_phy_18S_prune_soil, aes(x= Index, y= LibrarySize, color = Sample.Control)) +
  geom_point()
ggsave("results/soil/18S_Soil_library-size-sample-control-soil-samples.pdf", width = 8, height = 6, dpi = 150)

# Identify Contaminants using package decontam - Frequency method
soil_18S_contamdf_freq <- isContaminant(phy_18S_prune_soil, method="frequency", conc="QubitConc")
head(soil_18S_contamdf_freq)
tail(soil_18S_contamdf_freq)
table(soil_18S_contamdf_freq$contaminant) # 4196 false    66 true
head(which(soil_18S_contamdf_freq$contaminant))
soil_18S_contamdf_freq_list <- rownames_to_column(soil_18S_contamdf_freq, var = "ASV")
write.csv(soil_18S_contamdf_freq_list, "exported_tables/soil/soil_18S_contamdf_freq_list.csv")

# Explore AVSs that are potentially contaminants, ASV abundance should be inversely correlated with DNA contration
plot_frequency(phy_18S_prune_soil, taxa_names(phy_18S_prune_soil)[c(14, 43, 138, 165)], conc="QubitConc") + 
  xlab("DNA Concentration (ng/ul)")

# Explore AVSs that are potentially contaminants, randomly choosing 4 AVSs.
set.seed(100)
plot_frequency(phy_18S_prune_soil, taxa_names(phy_18S_prune_soil)[sample(which(soil_18S_contamdf_freq$contaminant),4)], conc="QubitConc") + 
  xlab("DNA Concentration (ng/ul)")

# Identify Contaminants using package decontam - Prevalence method
sample_data(phy_18S_prune_soil)$is.neg <- sample_data(phy_18S_prune_soil)$Sample.Control == "Control.Sample"
soil_18S_contamdf_prev <- isContaminant(phy_18S_prune_soil, method="prevalence", neg="is.neg")
table(soil_18S_contamdf_prev$contaminant) # 4249 false    13 true
head(which(soil_18S_contamdf_prev$contaminant))
soil_18S_contamdf_prev_list <- rownames_to_column(soil_18S_contamdf_prev, var = "ASV")
write.csv(soil_18S_contamdf_prev, "exported_tables/soil/soil_18S_contamdf_prev.csv")

# Identify Contaminants using package decontam - Prevalence method, increase threshhold value to 0.5
soil_18S_contamdf_prev05 <- isContaminant(phy_18S_prune_soil, method="prevalence", neg="is.neg", threshold=0.5)
table(soil_18S_contamdf_prev05$contaminant) # 4234 false    28 true
head(which(soil_18S_contamdf_prev05$contaminant))
soil_18S_contamdf_prev05_list <- rownames_to_column(soil_18S_contamdf_prev05, var = "ASV")
write.csv(soil_18S_contamdf_prev05_list, "exported_tables/soil/soil_18S_contamdf_prev05_list.csv")

# Make phyloseq object of presence-absence in negative controls and true samples
phy_soil_18S_prune_pa <- transform_sample_counts(phy_18S_prune_soil, function(abund) 1*(abund>0))
phy_soil_18S_prune_pa_neg <- prune_samples(sample_data(phy_soil_18S_prune_pa)$Sample.Control == "Control.Sample", phy_soil_18S_prune_pa)
phy_soil_18S_prune_pa_pos <- prune_samples(sample_data(phy_soil_18S_prune_pa)$Sample.Control == "True.Sample", phy_soil_18S_prune_pa)

# Make data.frame of prevalence in positive and negative samples
df_phy_18S_prune_soil_pa <- data.frame(pa.pos=taxa_sums(phy_soil_18S_prune_pa_pos), pa.neg=taxa_sums(phy_soil_18S_prune_pa_neg),
                              contaminant=soil_18S_contamdf_prev$contaminant)
ggplot(data=df_phy_18S_prune_soil_pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
ggsave("results/soil/soil_18S_contamdf_prev.pdf", width = 8, height = 6, dpi = 150)

# Make data.frame of prevalence in positive and negative samples, threshold of 0.5
df_phy_18S_prune_soil_pa05 <- data.frame(pa.pos=taxa_sums(phy_soil_18S_prune_pa_pos), pa.neg=taxa_sums(phy_soil_18S_prune_pa_neg),
                                   contaminant=soil_18S_contamdf_prev05$contaminant)
ggplot(data=df_phy_18S_prune_soil_pa05, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
ggsave("results/soil/soil_18S_contamdf_prev05.pdf", width = 8, height = 6, dpi = 150)

# Identify Contaminants using package decontam -  Combined method (i.e., Frequency and Prevalence)
soil_18S_contamdf_comb <- isContaminant(phy_18S_prune_soil, conc="QubitConc", neg="is.neg", threshold=0.1, detailed = TRUE, normalize = TRUE, method="combined")
table(soil_18S_contamdf_comb$contaminant) # 4246 false    16 true
head(which(soil_18S_contamdf_comb$contaminant))
soil_18S_contamdf_comb_list <- rownames_to_column(soil_18S_contamdf_comb, var = "ASV")
write.csv(soil_18S_contamdf_comb_list, "exported_tables/soil/soil_18S_contamdf_comb_list.csv")

soil_18S_contamdf_comb05 <- isContaminant(phy_18S_prune_soil, conc="QubitConc", neg="is.neg", threshold=0.5, detailed = TRUE, normalize = TRUE, method="combined")
table(soil_18S_contamdf_comb05$contaminant) # 3986 false    276 true
head(which(soil_18S_contamdf_comb05$contaminant))
soil_18S_contamdf_comb05_list <- rownames_to_column(soil_18S_contamdf_comb05, var = "ASV")
write.csv(soil_18S_contamdf_comb05_list, "exported_tables/soil/soil_18S_contamdf_comb05_list.csv")

# Remove contaminants from phyloseq object - prev05 based - check for samples with zero reads
phy_18S_prune_soil_noncontam_prev05 <- prune_taxa(!soil_18S_contamdf_prev05$contaminant, phy_18S_prune_soil)
phy_18S_prune_soil_noncontam_prev05 # 4234 taxa and 53 samples
smin <-min(sample_sums(phy_18S_prune_soil_noncontam_prev05))
any(sample_sums(phy_18S_prune_soil_noncontam_prev05) == 0) # if TRUE, there are samples with no reads, otherwise FALSE
sum(sample_sums(phy_18S_prune_soil_noncontam_prev05) == 0) # gives the number of cases

# Check for ASVs that have no counted reads
any(taxa_sums(phy_18S_prune_soil_noncontam_prev05) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_18S_prune_soil_noncontam_prev05) == 0) # gives the number of cases
phy_18S_prune_soil_noncontam_prev05 # 4234 taxa and 53 samples

#Subset by sample types, that is, only soil samples removing controls and blanks
phy_18S_prune_soil_noncontam_prev05_true <- subset_samples(phy_18S_prune_soil_noncontam_prev05, Sample.Type =="Soil") %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_18S_prune_soil_noncontam_prev05_true # 4180 taxa and 47 samples
any(taxa_sums(phy_18S_prune_soil_noncontam_prev05_true) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_18S_prune_soil_noncontam_prev05_true) == 0) # gives the number of cases
smin <- min(sample_sums(phy_18S_prune_soil_noncontam_prev05_true)) # smin = 39

#Checking the total number of reads and its distribution across samples
soil_18S_readsums_df = data.frame(nreads = sort(taxa_sums(phy_18S_prune_soil_noncontam_prev05_true), TRUE), sorted = 1:ntaxa(phy_18S_prune_soil_noncontam_prev05_true), 
                              type = "ASVs")
soil_18S_readsums_df = rbind(soil_18S_readsums_df, data.frame(nreads = sort(sample_sums(phy_18S_prune_soil_noncontam_prev05_true), 
                                                                    TRUE), sorted = 1:nsamples(phy_18S_prune_soil_noncontam_prev05_true), type = "Samples"))
title = "Total number of reads"
p = ggplot(soil_18S_readsums_df, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
ggsave("results/soil/18S_Soil_ShSk_ASV_Read_Profile_filt_prev05-soil-samples.pdf", width = 12, height = 6, dpi = 150)

# Scale reads to even sequencing depth by removing samples with a total number of reads < 500
smin <- min(sample_sums(phy_18S_prune_soil_noncontam_prev05_true)) # smin 39
phy_18S_prune_soil_noncontam_prev05_true_scale <- prune_samples(sample_sums(phy_18S_prune_soil_noncontam_prev05_true)>500, phy_18S_prune_soil_noncontam_prev05_true) %>%
  prune_taxa(taxa_sums(.) > 0, .)
smin <- min(sample_sums(phy_18S_prune_soil_noncontam_prev05_true_scale)) # smin 1220
phy_18S_prune_soil_noncontam_prev05_true_scale # 4175 taxa and 45 samples
any(taxa_sums(phy_18S_prune_soil_noncontam_prev05_true_scale) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_18S_prune_soil_noncontam_prev05_true_scale) == 0) # gives the number of cases

# Remove AVSs matching humans
phy_18S_prune_soil_noncontam_prev05_true_scale_filt <- phy_18S_prune_soil_noncontam_prev05_true_scale %>%
  subset_taxa(
    D16  != "Mammalia" &
      D17   != "Homo sapiens (human)"
  )
phy_18S_prune_soil_noncontam_prev05_true_scale_filt # 4157 taxa and 45 samples
any(sample_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(sample_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt) == 0) # gives the number of cases

# write ASV_soil matrix to file as csv
ASV_SHSK_18S_soil = as(otu_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt), "matrix")

# Transpose matrix
if(taxa_are_rows(phy_18S_prune_soil_noncontam_prev05_true_scale_filt)) {ASV_SHSK_18S_soil <- t(ASV_SHSK_18S_soil)}

# Coerce to a data frame
ASV_SHSK_18S_soil_df = as.data.frame(ASV_SHSK_18S_soil)
class(ASV_SHSK_18S_soil_df)
# write ASV_df matrix to file as csv
write.csv(ASV_SHSK_18S_soil_df, "exported_tables/soil/ASV_SHSK_18S_soil_df.csv")

# write TAX_soil matrix to file as csv
TAX_SHSK_18S_soil = as(tax_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt), "matrix")
TAX_SHSK_18S_soil_df = as.data.frame(TAX_SHSK_18S_soil)
class(TAX_SHSK_18S_soil_df)

# write TAX_df matrix to file as csv
write.csv(TAX_SHSK_18S_soil_df, "exported_tables/soil/TAX_SHSK_18S_soil_df.csv")

# write TAX_soil matrix to file as csv
SAMPLE_SHSK_18S_soil = as(sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt), "matrix")
SAMPLE_SHSK_18S_soil_df = as.data.frame(SAMPLE_SHSK_18S_soil)
class(SAMPLE_SHSK_18S_soil_df)

# write Sample_df matrix to file as csv
write.csv(SAMPLE_SHSK_18S_soil_df, "exported_tables/soil/SAMPLE_SHSK_18S_soil_df.csv")

# Save clean phyloseq object with soil data 
saveRDS(phy_18S_prune_soil_noncontam_prev05_true_scale_filt, "data/phy_18S_prune_soil_noncontam_prev05_true_scale_filt.RDS")

# Import phyloseq object for further analysis
phy_18S_prune_soil_noncontam_prev05_true_scale_filt <- readRDS("data/phy_18S_prune_soil_noncontam_prev05_true_scale_filt.RDS")
phy_18S_prune_soil_noncontam_prev05_true_scale_filt # 4157 taxa and 45 samples
smin <- min(sample_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt)) # smin 1220

######################### Generate a data.frame with alpha diversity measures #########################
# Calculate alpha-diversity measures (For plot purposes only!)
alpha_div_18S_soil <- data.frame(
  "Observed" = phyloseq::estimate_richness(phy_18S_prune_soil_noncontam_prev05_true_scale_filt, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phy_18S_prune_soil_noncontam_prev05_true_scale_filt, measures = "Shannon"),
  "InvSimpson" = phyloseq::estimate_richness(phy_18S_prune_soil_noncontam_prev05_true_scale_filt, measures = "InvSimpson"),
  "Reads" = phyloseq::sample_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt),
  "Soil.Type" = phyloseq::sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt)$Soil.Type,
  "Habitat" = phyloseq::sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt)$Habitat,
  "Sample.Site" = phyloseq::sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt)$Sample.Site)
alpha_div_18S_soil$Evenness <- alpha_div_18S_soil$Shannon/log(alpha_div_18S_soil$Observed)
head(alpha_div_18S_soil)

# Rename variable InvSimpson to Simpson and
alpha_div_18S_soil <- alpha_div_18S_soil %>%
  rename(Simpson = InvSimpson)
head(alpha_div_18S_soil)

# Reorder dataframe, first categorical then numerical variables  
alpha_div_18S_soil <- alpha_div_18S_soil[, c(7, 6, 5, 4, 1, 2, 3, 8)]
head(alpha_div_18S_soil)
write.csv(alpha_div_18S_soil, "exported_tables/soil/alpha_div_18S_soil.csv")

#Summarize alpha diversity measures
summary_alpha_18S_soil <- alpha_div_18S_soil %>%
  group_by(Soil.Type, Habitat, Sample.Site) %>%
  summarise(
    count = n(),
    mean_reads = mean(Reads),
    sd_reads = sd(Reads),
    mean_observed = mean(Observed),
    sd_observed = sd(Observed),
    mean_shannon = mean(Shannon),
    sd_shannon = sd(Shannon),
    mean_invsimpson = mean(Simpson),
    sd_invsimpson = sd(Simpson),
    mean_evenness = mean(Evenness),
    sd_evenness = sd(Evenness))
write.csv(summary_alpha_18S_soil, "exported_tables/soil/summary_alpha_18S_soil.csv")

#Summarize alpha diversity measures
summary_alpha_18S_habitat <- alpha_div_18S_soil %>%
  group_by(Soil.Type, Habitat) %>%
  summarise(
    count = n(),
    mean_reads = mean(Reads),
    sd_reads = sd(Reads),
    mean_observed = mean(Observed),
    sd_observed = sd(Observed),
    mean_shannon = mean(Shannon),
    sd_shannon = sd(Shannon),
    mean_invsimpson = mean(Simpson),
    sd_invsimpson = sd(Simpson),
    mean_evenness = mean(Evenness),
    sd_evenness = sd(Evenness))
write.csv(summary_alpha_18S_habitat, "exported_tables/soil/summary_alpha_18S_habitat.csv")

#Summarize alpha diversity measures
summary_alpha_18S_soil_type <- alpha_div_18S_soil %>%
  group_by(Soil.Type) %>%
  summarise(
    count = n(),
    mean_reads = mean(Reads),
    sd_reads = sd(Reads),
    mean_observed = mean(Observed),
    sd_observed = sd(Observed),
    mean_shannon = mean(Shannon),
    sd_shannon = sd(Shannon),
    mean_invsimpson = mean(Simpson),
    sd_invsimpson = sd(Simpson),
    mean_evenness = mean(Evenness),
    sd_evenness = sd(Evenness))
write.csv(summary_alpha_18S_soil_type, "exported_tables/soil/summary_alpha_18S_soil_type.csv")

# KW analysis alpha_18S_soil on all univar at once
# Remember, numerical variables are from columns 4-8. The test is by habitat, column 2
kw_alpha_18S_soil_habitat_univ <- as.data.frame(sapply(4:8, function(x) kruskal.test(alpha_div_18S_soil[,x],
                                                                             alpha_div_18S_soil[,2])))

# Rename columns with the proper variable names
kw_alpha_18S_soil_habitat_univ <- kw_alpha_18S_soil_habitat_univ %>%
  rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)

# Make some adjustments on the table
kw_alpha_18S_soil_habitat_univ <- t(kw_alpha_18S_soil_habitat_univ) # transpose
kw_alpha_18S_soil_habitat_univ <- as_tibble(kw_alpha_18S_soil_habitat_univ, rownames = "Metric") # adding rownames as a column
kw_alpha_18S_soil_habitat_univ <- kw_alpha_18S_soil_habitat_univ[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_18S_soil_habitat_univ) # checking object class

# Save resulting table with fwrite to avoid and issues with characters
fwrite(kw_alpha_18S_soil_habitat_univ, "exported_tables/soil/kw_alpha_18S_soil_habitat_univ.csv")

# KW analysis alpha_18S_soil one by one
kw_alpha_18S_soil_obs <- kruskal.test(Observed ~ Habitat, data = alpha_div_18S_soil)
#kw_alpha_18S_soil_h <- kruskal.test(Shannon ~ Habitat, data = alpha_div_18S_soil)
#kw_alpha_18S_soil_d <- kruskal.test(InvSimpson ~ Habitat, data = alpha_div_18S_soil)
#kw_alpha_18S_soil_j <- kruskal.test(Evenness ~ Habitat, data = alpha_div_18S_soil)
#kw_alpha_18S_soil_n <- kruskal.test(Reads ~ Habitat, data = alpha_div_18S_soil)

# pairwise comparison including p-value adjustment
kw_alpha_18S_soil_habitat_univ_pwt <- as.data.frame(sapply(4:8, function(x) pairwise.t.test(alpha_div_18S_soil[,x],
                                                                                        alpha_div_18S_soil[,2], p.adjust.method = "BH")))
kw_alpha_18S_soil_habitat_univ_pwt <- kw_alpha_18S_soil_habitat_univ_pwt %>%
  rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)

# Make some adjustments on the table
kw_alpha_18S_soil_habitat_univ_pwt <- t(kw_alpha_18S_soil_habitat_univ_pwt) # transpose
kw_alpha_18S_soil_habitat_univ_pwt <- as_tibble(kw_alpha_18S_soil_habitat_univ_pwt, rownames = "Metric") # adding rownames as a column
class(kw_alpha_18S_soil_habitat_univ_pwt) # checking object class

# Save resulting table with fwrite to avoid and issues with characters
fwrite(kw_alpha_18S_soil_habitat_univ_pwt, "exported_tables/soil/kw_alpha_18S_soil_habitat_univ_pwt.csv")

# pairwise comparison including p-value adjustment
# pwt_alpha_18S_soil_obs <- pairwise.t.test(alpha_div_18S_soil$Observed, alpha_div_18S_soil$Habitat, p.adjust.method = "BH")
#pwt_alpha_18S_soil_h <- pairwise.t.test(alpha_div_18S_soil$Shannon, alpha_div_18S_soil$Habitat, p.adjust.method = "BH")
#pwt_alpha_18S_soil_d <- pairwise.t.test(alpha_div_18S_soil$InvSimpson, alpha_div_18S_soil$Habitat, p.adjust.method = "BH")
#pwt_alpha_18S_soil_j <- pairwise.t.test(alpha_div_18S_soil$Evenness, alpha_div_18S_soil$Habitat, p.adjust.method = "BH")
#pwt_alpha_18S_soil_n <- pairwise.t.test(alpha_div_18S_soil$Reads, alpha_div_18S_soil$Habitat, p.adjust.method = "BH")

# Remember, numerical variables are from columns 4-8. The test is by soil type, column 3
kw_alpha_18S_soil_type_univ <- as.data.frame(sapply(4:8, function(x) kruskal.test(alpha_div_18S_soil[,x],
                                                                             alpha_div_18S_soil[,3])))

# Rename columns with the proper variable names
kw_alpha_18S_soil_type_univ <- kw_alpha_18S_soil_type_univ %>%
  rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)

# Make some adjustments on the table
kw_alpha_18S_soil_type_univ <- t(kw_alpha_18S_soil_type_univ) # transpose
kw_alpha_18S_soil_type_univ <- as_tibble(kw_alpha_18S_soil_type_univ, rownames = "Metric") # adding rownames as a column
kw_alpha_18S_soil_type_univ <- kw_alpha_18S_soil_type_univ[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_18S_soil_type_univ) # checking object class

# Save resulting table with fwrite to avoid and issues with characters
fwrite(kw_alpha_18S_soil_type_univ, "exported_tables/soil/kw_alpha_18S_soil_type_univ.csv")

# Plot alpha diversity measures for all three experiments
alpha_div_18S_soil$Habitat <- factor(alpha_div_18S_soil$Habitat,
                                         level = c("Chaparral", "Coastal.scrub", "Native.grass",
                                                   "Hollyleaf.cherry", "Oak.wood", "Riparian"),
                                         labels = c("CHA", "CSS", "NGR", "HLC",  "OWL", "RIP"))

alpha_div_18S_soil$Soil.Type <- factor(alpha_div_18S_soil$Soil.Type,
                                           level = c("Compact.Soil", "Not.Compact.Soil"),
                                           labels = c("Compacted", "Uncompacted"))

alpha_color <- c("#52392F","#83643E", "#C29B6C",  "#397A4C", "#77C063", "#BEDB92")
                   
alpha_div_18S_soil %>%
  gather(key = metric, value = value, c("Reads", "Observed", "Shannon", "Simpson", "Evenness")) %>%
  mutate(metric = factor(metric, levels = c("Reads", "Observed", "Shannon", "Simpson", "Evenness"))) %>%
  ggplot(aes(x = Habitat, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Habitat), height = 0, width = .2) +
  facet_grid(metric ~ Soil.Type, scales = "free") +
  scale_color_manual(values = alpha_color) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(expand = c(0.1, 0.1)) + # plot as % and removes the internal margins
  #scale_x_discrete(
   # labels = c("CHA", "CSS", "NGR", "HLC",  "OWL", "RIP"),
    #expand = c(0.1, 0.1),
    #drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  theme(legend.position="none")
ggsave("results/soil/18S_soil_alpha_diversity_phy_soil_18S_prune_noncontam_prev05_true_scale_habitat.pdf", width = 8, height = 7, dpi = 150) # save graphic

alpha_color_soil <- c("#52392F", "#397A4C")

alpha_div_18S_soil %>%
  gather(key = metric, value = value, c("Reads", "Observed", "Shannon", "Simpson", "Evenness")) %>%
  mutate(metric = factor(metric, levels = c("Reads", "Observed", "Shannon", "Simpson", "Evenness"))) %>%
  ggplot(aes(x = Soil.Type, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Soil.Type), height = 0, width = .2) +
  facet_grid(metric ~ Soil.Type, scales = "free") +
  scale_color_manual(values = alpha_color_soil) +
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
ggsave("results/soil/18S_soil_alpha_diversity_phy_soil_18S_prune_noncontam_prev05_true_scale_soil_type.pdf", width = 5, height = 7, dpi = 150) # save graphic

# Transform number of reads to relative abundance
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp <- microbiome::transform(phy_18S_prune_soil_noncontam_prev05_true_scale_filt, "compositional") # transform reads to relative abundance (0.0 - 1.0)
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_hel <- microbiome::transform(phy_18S_prune_soil_noncontam_prev05_true_scale_filt, "hellinger") # transform reads to relative abundance (0.0 - 1.0) and takes sqrt

# Check for ASVs that have no counted reads
any(taxa_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp) == 0) # gives the number of cases

# Ordinate NMDS
set.seed(1)
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nmds <- ordinate(
  physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp, 
  method = "NMDS", 
  distance = "bray",
  k =3,
  try = 50
)

set.seed(1)
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_hel_nmds <- ordinate(
  physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_hel, 
  method = "NMDS", 
  distance = "bray",
  k =2,
  try = 50
)


# Extract scores from the nmds
nmds_coord_soil <- as.data.frame(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nmds$points)
nmds_coord_soil <- as_tibble(nmds_coord_soil, rownames = "Sample") # adding rownames as a column

# Extract sample data from the phyloseq object
factors_nmds_soil <- as.data.frame(sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp))
factors_nmds_soil <- as_tibble(factors_nmds_soil, rownames = "Sample") # adding rownames as a column

# Combine both dataframes into one
nmds_3d_df_soil <- as.data.frame(left_join(nmds_coord_soil, factors_nmds_soil))
str(nmds_3d_df_soil)

# Create 3d nMDS plot
nmds_3d_df_soil$Habitat <- factor(nmds_3d_df_soil$Habitat,
                             levels = c("Chaparral", "Coastal.scrub", "Native.grass", "Hollyleaf.cherry", "Oak.wood", "Riparian"),
                             labels = c("CHA", "CSS", "NGR", "HLC", "OWL", "RIP"))

nmds_3d_df_soil$Soil.Type <- factor(nmds_3d_df_soil$Soil.Type,
                               levels = c("Compact.Soil", "Not.Compact.Soil"),
                               labels = c("Compacted", "Uncompacted"))

# Set colors per habitat
cols_habitat_soil <- c("#52392F", "#83643E", "#C29B6C",
                       "#397A4C", "#77C063", "#BEDB92")
cols_habitat_soil <- cols_habitat_soil[as.numeric(nmds_3d_df_soil$Habitat)]

# Set shapes per soil type
shap_type_soil <- c(19,17)
shap_type_soil <- shap_type_soil[as.numeric(nmds_3d_df_soil$Soil.Type)]

# Source function for adjusting grid
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')

par(mar=c(1,1,1,1))
sd3 <- scatterplot3d(x = nmds_3d_df_soil$MDS1, y = nmds_3d_df_soil$MDS2, z = nmds_3d_df_soil$MDS3, # adds the three axes
                     xlab = "NMDS1", ylab = "NMDS2", zlab = "NMDS3", # adds axes labels
                     pch = shap_type_soil, color = cols_habitat_soil,
                     grid = FALSE, box = FALSE,
                     col.grid = "gray", lty.grid=par("lty"),
                     cex.symbols = 2,
                     #angle = 60,
                     addgrids3d(nmds_3d_df_soil[, 2:4], grid = c("xy", "xz", "yz")))
legend("right", legend = levels(nmds_3d_df_soil$Habitat),
       pch = 19,
       col = c("#52392F", "#83643E", "#C29B6C", 
               "#397A4C", "#77C063", "#BEDB92"),
       bty = "n",
       inset = -0.25, xpd = TRUE)
legend("bottom", legend = levels(nmds_3d_df_soil$Soil.Type),
       pch = c(19, 17),
       col = "black",
       bty = "n",
       inset = -0.25, xpd = TRUE, horiz = TRUE)
ggsave("results/soil/18S_soil_ShSk_ASV_NMDS_greater_500_filt_prev05_comp_3d.pdf", width = 8, height = 6, dpi = 300)

# Plot NMDS 2D
theme_set(theme_bw())
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nmds_plot <- 
  phyloseq::plot_ordination(
  physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp,
  ordination = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nmds,
  type = "samples",
  color = "Habitat",
  shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E", "#C29B6C", 
                                "#397A4C", "#77C063", "#BEDB92"),
                  name = "Habitat",
                  breaks = c("Chaparral", "Coastal.scrub", "Native.grass", "Hollyleaf.cherry", "Oak.wood", "Riparian"),
                  labels = c("CHA", "CCS", "NGR", "HLC", "OWL", "RIP"))+
  scale_shape_manual(values = c(19,17),
                     name = "Soil Type",
                     breaks = c("Compact.Soil", "Not.Compact.Soil"),
                     labels = c("Compacted", "Uncompacted"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Habitat), alpha = 0.7, size = 3) +
  #geom_point(colour = "grey90", size = 1.5) +
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("NMDS2") + # add the title on y axis
  #xlab("NMDS1") + # add the title on x axis
  #ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  #annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
  annotate("text", x = 1.5, y = 2, label ="2D Stress: 0.13") +
  #ggtitle("18S rRNA - Eukaryotes - Soil") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =12))
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nmds_plot
ggsave("results/soil/18S_Soil_ShSk_ASV_NMDS_greater_500_filt_prev05_comp.pdf", width = 6, height = 4, dpi = 300)

theme_set(theme_bw())
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_hel_nmds_plot <- 
  phyloseq::plot_ordination(
    physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_hel,
    ordination = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_hel_nmds,
    type = "samples",
    color = "Habitat",
    shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E", "#C29B6C", 
                                         "#397A4C", "#77C063", "#BEDB92"),
                                         name = "Habitat",
                     breaks = c("Chaparral", "Coastal.scrub", "Native.grass", "Hollyleaf.cherry", "Oak.wood", "Riparian"),
                     labels = c("CHA", "CCS", "NGR", "HLC", "OWL", "RIP"))+
  scale_shape_manual(values = c(19,17),
                     name = "Soil Type",
                     breaks = c("Compact.Soil", "Not.Compact.Soil"),
                     labels = c("Compacted", "Uncompacted"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Habitat), alpha = 0.7, size = 3) +
  #geom_point(colour = "grey90", size = 1.5) +
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("NMDS2") + # add the title on y axis
  #xlab("NMDS1") + # add the title on x axis
  #ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  #annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
  annotate("text", x = 1.5, y = 2, label ="2D Stress: 0.12") +
  #ggtitle("18S rRNA - Eukaryotes - Soil") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =12))
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_hel_nmds_plot
ggsave("results/soil/18S_Soil_ShSk_ASV_NMDS_greater_500_filt_prev05_hel.pdf", width = 6, height = 4, dpi = 300)

write.csv(tax_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp), "exported_tables/soil/phy_soil_18S_prune_noncontam_prev05_true_scale_filt_comp.csv")

# Compact soil
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_compact <- subset_samples(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp, Soil.Type == "Compact.Soil") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_compact # 2000 taxa and 18 samples

# Uncompact soil
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_uncompact <- subset_samples(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp, Soil.Type == "Not.Compact.Soil") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_uncompact # 2000 taxa and 18 samples

#Check taxa and samples with zeros
any(taxa_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_compact) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_compact) == 0) # gives the number of cases
any(sample_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_compact) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(sample_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_compact) == 0) # gives the number of cases

any(taxa_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_uncompact) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_uncompact) == 0) # gives the number of cases
any(sample_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_uncompact) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(sample_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_uncompact) == 0) # gives the number of cases

# Ordination NMDS compacted soil type
set.seed(1)
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_compact_nmds <- ordinate(
  physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_compact, 
  method = "NMDS", 
  distance = "bray",
  try = 50 #Stress:     0.1484856 
)

# Ordination NMDS uncompacted soil type
set.seed(1)
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_uncompact_nmds <- ordinate(
  physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_uncompact, 
  method = "NMDS", 
  distance = "bray",
  try = 50 #Stress:     0.112822 
)

# Plot NMDS based on compositional transformation, only samples with >= 500 reads
# Compacted soil
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_compact,
  ordination = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_compact_nmds,
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
  annotate("text", x = 1.5, y = 1, label ="2D Stress: 0.15")
ggsave("results/soil/phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_compact.pdf", width = 6, height = 4.5, dpi = 200)

# Uncompacted soil
theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_uncompact,
  ordination = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_uncompact_nmds,
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
  annotate("text", x = -1.5, y = 1.2, label ="2D Stress: 0.11")
ggsave("results/soil/phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_uncompact.pdf", width = 6, height = 4.5, dpi = 200)


# Subset for Nematoda only for Aldex
get_taxa_unique(phy_18S_prune_soil_noncontam_prev05_true_scale_filt, taxonomic.rank = "D13")
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema <- phy_18S_prune_soil_noncontam_prev05_true_scale_filt %>%
  subset_taxa(
    D13  == "Nematoda") %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema
get_taxa_unique(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema, taxonomic.rank = "D13")
smin <- min(sample_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema))

# Export phyloseq object with nematode taxa/ASV only
# This can be used as a starting point for further analysis
saveRDS(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema, "data/phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema.RDS")

# Import phyloseq object
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema <- readRDS("data/phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema.RDS")
class(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema) # check it is a phyloseq object
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema # 109 taxa and 45 samples

############ save only Nematoda to create a new phyloseq object ############ 
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_asv <- as.data.frame(otu_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema, "matrix"))
write.csv(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_asv, "exported_tables/phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_asv.csv")

phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_tax <- as.data.frame(tax_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema))
write.csv(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_tax, "exported_tables/phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_tax.csv")

phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_metadata <- as.data.frame(sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema))
write.csv(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_metadata, "exported_tables/phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_metadata.csv")
############################################################     
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp <- microbiome::transform(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema, "compositional") # transform reads to relative abundance (0.0 - 1.0)
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_hel <- microbiome::transform(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema, "hellinger") # relative abundance and sqrt
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_log <- microbiome::transform(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema, "log10p") # log
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_clr <- microbiome::transform(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema, "clr") # centered log ratio

############### Ordinate NMDS Nematoda only
set.seed(1)
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp_nmds <- ordinate(
  physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp, 
  method = "NMDS", 
  distance = "bray",
  k =2,
  try = 50
)

set.seed(1)
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_hel_nmds <- ordinate(
  physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_hel, 
  method = "NMDS", 
  distance = "bray",
  k =2,
  try = 50
)

set.seed(1)
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_log_nmds <- ordinate(
  physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_log, 
  method = "NMDS", 
  distance = "bray",
  k =2,
  try = 50
)

set.seed(1)
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_clr_nmds <- ordinate(
  physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_clr, 
  method = "NMDS", 
  distance = "euclidean",
  k =2,
  try = 50
)


#########
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp_sel <- subset_samples(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp,
                                                                                    sample_names(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp) != "SK.15.1.Soil" &
                                                                                      sample_names(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp) != "SK.14.1.Soil") # this sample was removed
sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp_sel)


set.seed(1)
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp_sel_nmds <- ordinate(
  physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp_sel, 
  method = "NMDS", 
  distance = "bray",
  k =2,
  try = 50
)

theme_set(theme_bw())
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp_nmds_plot <- 
  phyloseq::plot_ordination(
    physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp,
    ordination = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp_nmds,
    type = "samples",
    color = "Habitat",
    shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E", "#C29B6C", 
                                "#397A4C", "#77C063", "#BEDB92"),
                     name = "Habitat",
                     breaks = c("Chaparral", "Coastal.scrub", "Native.grass", "Hollyleaf.cherry", "Oak.wood", "Riparian"),
                     labels = c("CHA", "CCS", "NGR", "HLC", "OWL", "RIP"))+
  scale_shape_manual(values = c(19,17),
                     name = "Soil Type",
                     breaks = c("Compact.Soil", "Not.Compact.Soil"),
                     labels = c("Compacted", "Uncompacted"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Habitat), alpha = 0.7, size = 3) +
  #geom_point(colour = "grey90", size = 1.5) +
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("NMDS2") + # add the title on y axis
  #xlab("NMDS1") + # add the title on x axis
  #ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  #annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
  annotate("text", x = 0.5, y = 1, label ="2D Stress: 0.14") +
  #ggtitle("18S rRNA - Eukaryotes - Soil") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =12))
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp_nmds_plot
#ggsave("results/soil/18S_soil_nema_filt_prev05_comp_nematode_only.pdf", width = 6, height = 4, dpi = 300)

theme_set(theme_bw())
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_clr_nmds_plot <- 
  phyloseq::plot_ordination(
    physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_clr,
    ordination = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_clr_nmds,
    type = "samples",
    color = "Habitat",
    shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E", "#C29B6C", 
                                         "#397A4C", "#77C063", "#BEDB92"),
                                         name = "Habitat",
                     breaks = c("Chaparral", "Coastal.scrub", "Native.grass", "Hollyleaf.cherry", "Oak.wood", "Riparian"),
                     labels = c("CHA", "CCS", "NGR", "HLC", "OWL", "RIP"))+
  scale_shape_manual(values = c(19,17),
                     name = "Soil Type",
                     breaks = c("Compact.Soil", "Not.Compact.Soil"),
                     labels = c("Compacted", "Uncompacted"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(aes(color = Habitat), alpha = 0.7, size = 3) +
  #geom_point(colour = "grey90", size = 1.5) +
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("NMDS2") + # add the title on y axis
  #xlab("NMDS1") + # add the title on x axis
  #ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  #annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
  annotate("text", x = 0.5, y = 1, label ="2D Stress: 0.14") +
  #ggtitle("18S rRNA - Eukaryotes - Soil") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =12))
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_clr_nmds_plot
#ggsave("results/soil/18S_soil_nema_filt_prev05_hel_nematode_only.pdf", width = 6, height = 4, dpi = 300)

theme_set(theme_bw())
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp_sel_nmds_plot <- 
  phyloseq::plot_ordination(
    physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp_sel,
    ordination = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp_sel_nmds,
    type = "samples",
    color = "Habitat",
    shape = "Soil.Type") +
  scale_color_manual(values = c("#52392F", "#83643E", "#C29B6C", 
                                         "#397A4C", "#77C063", "#BEDB92"),
                                         name = "Habitat",
                     breaks = c("Chaparral", "Coastal.scrub", "Native.grass", "Hollyleaf.cherry", "Oak.wood", "Riparian"),
                     labels = c("CHA", "CCS", "NGR", "HLC", "OWL", "RIP"))+
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
  #geom_point(colour = "grey90", size = 1.5) +
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #ylab("NMDS2") + # add the title on y axis
  #xlab("NMDS1") + # add the title on x axis
  #ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  #annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
  annotate("text", x = 0.5, y = 1, label ="2D Stress: 0.14") +
  #ggtitle("18S rRNA - Eukaryotes - Soil") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =12))
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_comp_sel_nmds_plot
ggsave("results/soil/18S_soil_nema_filt_prev05_comp_nematode_sel_only.pdf", width = 6, height = 4, dpi = 300)

# write ASV_Soil_nema matrix to file as csv
ASV_SHSK_18S_nema_soil = as(otu_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema), "matrix")

# Transpose matrix
if(taxa_are_rows(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema)) {ASV_SHSK_18S_nema_soil <- t(ASV_SHSK_18S_nema_soil)}

# Coerce to a data frame
ASV_SHSK_18S_nema_soil_df = as.data.frame(ASV_SHSK_18S_nema_soil)
class(ASV_SHSK_18S_nema_soil_df)
# write ASV_df matrix to file as csv
write.csv(ASV_SHSK_18S_nema_soil_df, "exported_tables/soil/ASV_SHSK_18S_nema_soil_df.csv")

# write TAX_soil matrix to file as csv
TAX_SHSK_18S_soil_nema = as(tax_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema), "matrix")
TAX_SHSK_18S_soil_nema_df = as.data.frame(TAX_SHSK_18S_soil_nema)
class(TAX_SHSK_18S_soil_nema_df)
# write TAX_df matrix to file as csv
write.csv(TAX_SHSK_18S_soil_nema_df, "exported_tables/soil/TAX_SHSK_18S_soil_nema_df.csv")

# write TAX_soil matrix to file as csv
SAMPLE_SHSK_18S_soil_nema = as(sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema), "matrix")
SAMPLE_SHSK_18S_soil_nema_df = as.data.frame(SAMPLE_SHSK_18S_soil_nema)
class(SAMPLE_SHSK_18S_soil_nema_df)
# write Sample_df matrix to file as csv
write.csv(SAMPLE_SHSK_18S_soil_nema_df, "exported_tables/soil/SAMPLE_SHSK_18S_soil_nema_df.csv")

# Subset for Nematoda only and create barplots
get_taxa_unique(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp, taxonomic.rank = "D13")
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nema <- phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp %>%
  subset_taxa(
    D13  == "Nematoda") %>%
  prune_taxa(taxa_sums(.) > 0, .)
get_taxa_unique(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nema, taxonomic.rank = "D13")

# Save nematode as a matrix
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nema_df <- as.data.frame(otu_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nema))
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nema_df <- rownames_to_column(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nema_df, var = "ASV")
write.csv(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nema_df,
          "exported_tables/soil/phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nema_df.csv")

# Save as long format
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nema_sum <- phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nema %>%
  psmelt()
write.csv(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nema_sum,
          "exported_tables/soil/phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nema_sum.csv")

# Import nematode matrix from previous lines
phy_soil_18S_nema <- read_excel("exported_tables/soil/phy_soil_18S_prune_noncontam_prev05_true_scale_filt_comp_nema.xlsx")

# Check the unique entries for D20
unique(phy_soil_18S_nema$D20)

# Define color palette for family
nb.cols.nema.fam <- 17
soil_colors_nema_fam <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.nema.fam)

# Check the unique entries for D21
unique(phy_soil_18S_nema$D21)

# Define color palette for genus
nb.cols.nema.genus <- 23
soil_colors_nema_genus <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.nema.genus)

colors_top23 <- c("#9a6324", "#46f0f0", "#1F78B4", "#aaffc3", "#e6beff",  "#33A02C", "#4363d8", "#008080", "#FB9A99", "#e6194b",
                  "#f032e6", "#bcf60c", "#fabebe","#f58231", "#ffe119", "#6A3D9A", "#000075", "#fffac8", "#800000", "gray", "lightblue", "black", "darkgreen")

# Check the unique entries for feeding group
unique(phy_soil_18S_nema$Feeding)

#  Define color palette for feeding group
nb.cols.nema.feeding <- 6
soil_colors_nema_feeding <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.nema.feeding)

#Change facet habitat labels
soil.type.labs <- c("Compacted","Uncompacted")
names(soil.type.labs) <- c("Compact.Soil","Not.Compact.Soil")

habitat.type.labs <- c("CHA", "CSS", "HLC", "NGR", "OWL", "RIP")
names(habitat.type.labs) <- c("Chaparral", "Coastal.scrub", "Hollyleaf.cherry", "Native.grass", "Oak.wood", "Riparian")

site.labs <- c("01", "02", "03", "04", "05", "06",
               "07", "08", "09", "10", "11", "12",
               "13", "14", "15", "16", "17", "18")

names(site.labs) <- c("SK.01", "SK.02", "SK.03", "SK.04", "SK.05", "SK.06",
                      "SK.07", "SK.08", "SK.09", "SK.10", "SK.11", "SK.12",
                      "SK.13", "SK.14", "SK.15", "SK.16", "SK.17", "SK.18")

# Plot bar chart 18S Soil samples per habitat - Nematoda only - family level - per site
phy_soil_18S_nema_fam_site <- ggplot(phy_soil_18S_nema, aes(x = Sample.Site, y = Abundance, fill = D20)) + #plotting by sample
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.type.labs,
                                   Habitat = habitat.type.labs)) + # By Habitat
  geom_bar(stat = "identity", width = 0.95, position = "fill") + # adds to 100%
  #geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c(colors_top23)) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)) + # adjusts text of y axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(expand = c(0.0, 0.0), label = site.labs) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  labs(fill = "Family") +
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background = element_rect("white")) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
phy_soil_18S_nema_fam_site
ggsave("results/soil/18S_Nematoda_Soil_Habitat_ShSk_barplot_filt_comp_family_level_site.pdf", width = 8, height = 4, dpi = 150) # save graphic

# Plot bar chart 18S Soil samples per habitat - Nematoda only - family level - per habitat
phy_soil_18S_nema_fam_habitat <- ggplot(phy_soil_18S_nema, aes(x = Habitat, y = Abundance, fill = D20)) + #plotting by sample
  facet_nested(. ~ Soil.Type, scales = "free",
               labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
  geom_bar(stat = "identity", width = 0.95, position = "fill") + # adds to 100%
  scale_fill_manual(values = c(colors_top23)) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)) + # adjusts text of y axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(expand = c(0.0, 0.0), label = habitat.type.labs) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  labs(fill = "Family") +
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background = element_rect("white")) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
phy_soil_18S_nema_fam_habitat
ggsave("results/soil/18S_Nematoda_Soil_Habitat_ShSk_barplot_filt_comp_family_level.pdf", width = 6, height = 4, dpi = 150) # save graphic

# Plot bar chart 18S Soil samples per habitat - Nematoda only
phy_soil_18S_nema_genus_site <- ggplot(phy_soil_18S_nema, aes(x = Sample.Site, y = Abundance, fill = D21)) + #plotting by sample
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.type.labs,
                                   Habitat = habitat.type.labs)) + # By Habitat
    geom_bar(stat = "identity", width = 0.95, position = "fill") + # adds to 100%
  #geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = colors_top23) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)) + # adjusts text of y axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(expand = c(0.0, 0.0), label = site.labs) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  labs(fill = "Genus") +
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background = element_rect("white")) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
phy_soil_18S_nema_genus_site
ggsave("results/soil/18S_Nematoda_Soil_Habitat_Site_ShSk_barplot_filt_comp_genus_level.pdf", width = 12, height = 4, dpi = 150) # save graphic

# Plot bar chart 18S Soil samples per habitat - Nematoda only
phy_soil_18S_nema_genus_habitat <- ggplot(phy_soil_18S_nema, aes(x = Habitat, y = Abundance, fill = D21)) + #plotting by sample
  facet_nested(. ~ Soil.Type, scales = "free",
               labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
  geom_bar(stat = "identity", width = 0.95, position = "fill") + # adds to 100%
  #geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = colors_top23) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)) + # adjusts text of y axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(expand = c(0.0, 0.0), label = habitat.type.labs) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  labs(fill = "Genus") +
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background = element_rect("white")) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
phy_soil_18S_nema_genus_habitat
ggsave("results/soil/18S_Nematoda_Soil_Habitat_ShSk_barplot_filt_comp_genus_level.pdf", width = 8, height = 4, dpi = 300) # save graphic

# Plot bar chart 18S Soil samples per habitat - Nematoda feeding
phy_soil_18S_nema_feeding <- ggplot(phy_soil_18S_nema, aes(x = Sample.Site, y = Abundance, fill = Feeding)) + #plotting by sample
  facet_nested(. ~ Soil.Type+Habitat, scales = "free_x",
               labeller = labeller(Soil.Type = soil.type.labs,
                                   Habitat = habitat.type.labs)) + # By Habitat
  geom_bar(stat = "identity", width = 0.95, position = "fill") + # adds to 100%
  #geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = c("#A6CEE3", "#98D277", "#F16667",
                               "#FE982C", "#B15928")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)) + # adjusts text of y axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(expand = c(0.0, 0.0), label = site.labs) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  labs(fill = "Feeding Group") +
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background = element_rect("white")) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
phy_soil_18S_nema_feeding
ggsave("results/soil/18S_Nematoda_Soil_Habitat_Site_ShSk_barplot_filt_comp_feeding_group.pdf", width = 8, height = 4, dpi = 150) # save graphic

phy_soil_18S_nema_combined <- ggarrange(phy_soil_18S_nema_fam_site, phy_soil_18S_nema_genus_site,
                                        phy_soil_18S_nema_feeding, legend = "right", ncol = 1, nrow = 3, labels = "AUTO",
                                        align = "v")
phy_soil_18S_nema_combined
ggsave("results/soil/18S_Nematoda_Soil_Habitat_ShSk_barplot_combined.pdf", width = 12, height = 12, dpi = 300)

# Checking the abundance of Fungi versus other groups
merge_less_than_top_prev05_top20_D5 <- function(phy_18S_prune_soil_noncontam_prev05_true_scale_filt, top=19){
  transformed <- transform_sample_counts(phy_18S_prune_soil_noncontam_prev05_true_scale_filt, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:7 if there are species level
  }
  return(merged)
}

merge_less_than_top_prev05_top20_D6 <- function(phy_18S_prune_soil_noncontam_prev05_true_scale_filt, top=19){
  transformed <- transform_sample_counts(phy_18S_prune_soil_noncontam_prev05_true_scale_filt, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:7 if there are species level
  }
  return(merged)
}

phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D5 <- tax_glom(phy_18S_prune_soil_noncontam_prev05_true_scale_filt, "D5")
phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D5_top20 <- merge_less_than_top_prev05_top20_D5(phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D5, top=19)
get_taxa_unique(phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D5_top20, taxonomic.rank = "D5")
phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D5_top20_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D5_top20) # Melt to long format
write.csv(phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D5_top20_df, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D5_top20_df.csv")

phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D5_top20_agr = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+D1+D2+D3+D4+D5, data=phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D5_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D5_top20_agr$D5)

phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D6 <- tax_glom(phy_18S_prune_soil_noncontam_prev05_true_scale_filt, "D6")
phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D6_top20 <- merge_less_than_top_prev05_top20_D6(phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D6, top=19)
get_taxa_unique(phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D6_top20, taxonomic.rank = "D6")
phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D6_top20_df = psmelt(phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D6_top20) # Melt to long format
write.csv(phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D6_top20_df, "exported_tables/soil/phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D6_top20_df.csv")

phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D6_top20_agr = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+D1+D2+D3+D4+D5+D6, data=phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D6_top20_df, FUN=mean) # add common factors to use for plotting
unique(phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D6_top20_agr$D6)

# Reorder taxa on D5 rank
phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D5_top20_agr$D5 <- factor(phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D5_top20_agr$D5,
                                                                                       levels = c("Cercomonadidae", "Cercomonas", "Cercozoa",
                                                                                                  "Chlorophyceae", "Choanozoa", "Chrysophyceae",
                                                                                                  "Cryomonadida", "Eustigmatales", "Fungi", "Glissomonadida",
                                                                                                  "Gregarinasina", "Heteromita", "Phragmoplastophyta",
                                                                                                  "Pythium", "Silicofilosea", "Streptophyta",
                                                                                                  "Thecofilosea", "Trebouxiophyceae", "Others", "Unassigned"))
# Reorder taxa on D6 rank
phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D6_top20_agr$D6 <- factor(phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D6_top20_agr$D6,
                                                                               levels = c("Cercomonadidae", "Cercomonas", "Cercozoa",
                                                                                          "Chlorophyceae", "Chrysophyceae", "Chytridiomycota",
                                                                                          "Dikarya", "Embryophyta", "Euglyphida", "Eugregarinorida",
                                                                                          "Eustigmatales", "Heteromita", "Metazoa", "Mucoromycota",
                                                                                          "Phragmoplastophyta", "Rhizaspididae", "Thecofilosea", "Trebouxiophyceae",
                                                                                          "Unassigned", "Others"))    
phy_soil_18S_major_groups_D5_site <- ggplot(phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D5_top20_agr, aes(x = Sample.Site, y = Abundance, fill = D5)) + #plotting by sample
  facet_nested(. ~ Soil.Type+Habitat, scales = "free", labeller = labeller(Soil.Type = soil.type.labs, Habitat = habitat.type.labs)) + # By Habitat
  geom_bar(stat = "identity", width = 0.95, position = "fill") + # adds to 100%
  scale_fill_manual(values = colors_top23) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)) + # adjusts text of y axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_x_discrete(expand = c(0.0, 0.0), label = site.labs) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + 
  labs(fill = "Groups") +
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background = element_rect("white")) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
phy_soil_18S_major_groups_D5_site
ggsave("results/soil/18S_soil_major_groups_D5_barplot_filt_comp.pdf", width = 10, height = 5, dpi = 150) # save graphic

phy_soil_18S_major_groups_D6_site <- ggplot(phy_16_prune_soil_noncontam_prev05_true_all_filt_top_D6_top20_agr, aes(x = Sample.Site, y = Abundance, fill = D6)) + #plotting by sample
  facet_nested(. ~ Soil.Type+Habitat, scales = "free", labeller = labeller(Soil.Type = soil.type.labs, Habitat = habitat.type.labs)) + # By Habitat
  geom_bar(stat = "identity", width = 0.95, position = "fill") + # adds to 100%
  scale_fill_manual(values = colors_top23) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)) + # adjusts text of y axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_x_discrete(expand = c(0.0, 0.0), label = site.labs) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + 
  labs(fill = "Groups") +
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background = element_rect("white")) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
phy_soil_18S_major_groups_D6_site
ggsave("results/soil/18S_soil_major_groups_D6_barplot_filt_comp.pdf", width = 10, height = 5, dpi = 150) # save graphic

# Subset for Metazoan only and create barplots
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_meta <- phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp %>%
  subset_taxa(
    D6  == "Metazoa") %>%
  psmelt()
unique(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_meta$D16)

write.csv(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_meta,
          "exported_tables/soil/phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_meta.csv")

# Import metazoan matrix from previous lines
phy_soil_18S_meta <- read.csv("exported_tables/soil/phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_meta_tax_fixed.csv")
unique(phy_soil_18S_meta$D14)

soil_colors_meta_phy <- c("#ffe119", "#800000", "#bcf60c", "#1F78B4", "lightblue", "#6A3D9A", "#33A02C", "#000075", "#FB9A99", "#e6194b",
                          "#f032e6", "#bcf60c", "#fabebe","#f58231", "#ffe119", "#6A3D9A", "#000075", "#fffac8", "#800000", "gray")

# Plot bar chart 18S Soil samples per habitat - Metazoan only - per site
phy_soil_18S_meta_site_plot <- ggplot(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_meta, aes(x = Sample.Site, y = Abundance, fill = D14)) + #plotting by sample
  facet_nested(. ~ Soil.Type+Habitat, scales = "free", labeller = labeller(Soil.Type = soil.type.labs, Habitat = habitat.type.labs)) + # By Habitat
  geom_bar(stat = "identity", width = 0.95, position = "fill") + # adds to 100%
  #geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = soil_colors_meta_phy) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)) + # adjusts text of y axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(expand = c(0.0, 0.0), label = site.labs) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + 
  labs(fill = "Groups") +
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background = element_rect("white")) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
phy_soil_18S_meta_site_plot
ggsave("results/soil/18S_metazoan_Soil_Habitat_Site_ShSk_barplot_filt_comp_class_level_site.pdf", width = 8, height = 4, dpi = 200) # save graphic

# Plot bar chart 18S Soil samples per habitat - Metazoan only - per habitat
phy_soil_18S_meta_habitat_plot <- ggplot(phy_soil_18S_meta, aes(x = Habitat, y = Abundance, fill = D14)) + #plotting by sample
  facet_nested(. ~ Soil.Type, scales = "free", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
  geom_bar(stat = "identity", width = 0.95, position = "fill") + # adds to 100%
  #geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = soil_colors_meta_phy) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10)) + # adjusts text of y axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(expand = c(0.0, 0.0), label = habitat.type.labs) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + 
  labs(fill = "Groups") +
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background = element_rect("white")) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
phy_soil_18S_meta_habitat_plot
ggsave("results/soil/18S_metazoan_Soil_Habitat_ShSk_barplot_filt_comp_class_level_habitat.pdf", width = 6, height = 4, dpi = 200) # save graphic

# Subset for Fungi only and create barplots
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi <- phy_18S_prune_soil_noncontam_prev05_true_scale_filt %>%
  subset_taxa(D5  == "Fungi") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)

phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi # 1544 taxa and 45 samples
any(taxa_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi) == 0) # gives the number of cases
smin <-min(sample_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi)) # smin = 501

# Apply different transformations to ASV table, e.g., converting number of reads to relative abundance
# Here we use the library microbiome
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi_comp <- microbiome::transform(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi, "compositional")
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungit_clr <- microbiome::transform(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi, "clr")
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi_log10 <- microbiome::transform(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi, "log10p")
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi_hellinger <- microbiome::transform(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi, "hellinger")

# NMDS ordination
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi_comp_nmds <- phyloseq::ordinate(
  physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi_comp,
  method = "NMDS",
  distance = "bray"
) # Stress:     0.1356364

theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi_comp,
  ordination = phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi_comp_nmds,
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
  geom_point(aes(color = Habitat), alpha = 0.7, size = 4) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 9)) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5, legend.position = "right") + # adjusts the title of the legend
  annotate("text", x = 1.5, y = 1.5, label ="2D Stress: 0.14")
ggsave("results/soil/18S_soil_fungi_only_plot_without_facet_soil_habitat.pdf", width = 7, height = 5, dpi = 300)

# Calculate alpha-diversity measures (For plot purposes only!)
# This can be done using the different phyloseq alpha diversity measures
# You will get a Warning message for each index since there is no singletons on the dataset
alpha_div_18S_fungi <- data.frame(
  "Observed" = phyloseq::estimate_richness(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi, measures = "Shannon"),
  "InvSimpson" = phyloseq::estimate_richness(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi, measures = "InvSimpson"),
  "Reads" = phyloseq::sample_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi),
  "Soil.Type" = phyloseq::sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi)$Soil.Type,
  "Habitat" = phyloseq::sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi)$Habitat,
  "Sample.Site" = phyloseq::sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_fungi)$Sample.Site)
alpha_div_18S_fungi$Evenness <- alpha_div_18S_fungi$Shannon/log(alpha_div_18S_fungi$Observed)
head(alpha_div_18S_fungi)

# Rename variable InvSimpson to Simpson
# The function rename & %>% works on dplyr. make sure it is loaded.
alpha_div_18S_fungi <- alpha_div_18S_fungi %>%
  dplyr::rename(Simpson = InvSimpson)
head(alpha_div_18S_fungi)

# Reorder dataframe, first categorical then numerical variables  
alpha_div_18S_fungi <- alpha_div_18S_fungi[, c(7, 6, 5, 4, 1, 2, 3, 8)]
head(alpha_div_18S_fungi)
write.csv(alpha_div_18S_fungi, "exported_tables/soil/alpha_div_18S_fungi.csv")

#Summarize alpha diversity measures by site
#Load library dplyr
alpha_div_18S_fungi_site <- alpha_div_18S_fungi %>%
  group_by(Soil.Type, Habitat, Sample.Site) %>%
    dplyr::summarise(
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
write.csv(alpha_div_18S_fungi_site, "exported_tables/soil/alpha_div_18S_fungi_site.csv")

#Summarize alpha diversity measures by habitat
alpha_div_18S_fungi_habitat <- alpha_div_18S_fungi %>%
  group_by(Soil.Type, Habitat) %>%
  dplyr::summarise(
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
write.csv(alpha_div_18S_fungi_habitat, "exported_tables/soil/alpha_div_18S_fungi_habitat.csv")

#Summarize alpha diversity measures by soil type
alpha_div_18S_fungi_soil_type <- alpha_div_18S_fungi %>%
  group_by(Soil.Type) %>%
  dplyr::summarise(
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
write.csv(alpha_div_18S_fungi_soil_type, "exported_tables/soil/alpha_div_18S_fungi_soil_type.csv")

# KW analysis alpha_16S_nem on all metrics at once
# Remember, numerical variables are from columns 4-8. The test is by soil type, column 3
kw_alpha_div_18S_fungi_soil_type_univ <- as.data.frame(sapply(4:8, function(x) kruskal.test(alpha_div_18S_fungi[,x],
                                                                                        alpha_div_18S_fungi[,3])))

# Rename columns with the proper variable names
kw_alpha_div_18S_fungi_soil_type_univ <- kw_alpha_div_18S_fungi_soil_type_univ %>%
  dplyr::rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)

# Make some adjustments on the table
kw_alpha_div_18S_fungi_soil_type_univ <- t(kw_alpha_div_18S_fungi_soil_type_univ) # transpose
kw_alpha_div_18S_fungi_soil_type_univ <- as_tibble(kw_alpha_div_18S_fungi_soil_type_univ, rownames = "Metric") # adding rownames as a column
kw_alpha_div_18S_fungi_soil_type_univ <- kw_alpha_div_18S_fungi_soil_type_univ[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_div_18S_fungi_soil_type_univ) # checking object class

# Save resulting table with fwrite to avoid any issues with characters
data.table::fwrite(kw_alpha_div_18S_fungi_soil_type_univ, "exported_tables/soil/kw_alpha_div_18S_fungi_soil_type_univ.csv")

# KW analysis alpha_16S_nem on all metrics at once
# Remember, numerical variables are from columns 4-8. The test is by habitat, column 2
kw_alpha_div_18S_fungi_habitat_univ <- as.data.frame(sapply(4:8, function(x) kruskal.test(alpha_div_18S_fungi[,x],
                                                                                          alpha_div_18S_fungi[,2])))

kw_alpha_div_18S_fungi_habitat_univ <- kw_alpha_div_18S_fungi_habitat_univ %>%
  dplyr::rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)

# Make some adjustments on the table
kw_alpha_div_18S_fungi_habitat_univ <- t(kw_alpha_div_18S_fungi_habitat_univ) # transpose
kw_alpha_div_18S_fungi_habitat_univ <- as_tibble(kw_alpha_div_18S_fungi_habitat_univ, rownames = "Metric") # adding rownames as a column
kw_alpha_div_18S_fungi_habitat_univ <- kw_alpha_div_18S_fungi_habitat_univ[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_div_18S_fungi_habitat_univ) # checking object class

# Save resulting table with fwrite to avoid and issues with characters
data.table::fwrite(kw_alpha_div_18S_fungi_habitat_univ, "exported_tables/soil/kw_alpha_div_18S_fungi_habitat_univ.csv")

# pairwise comparison including p-value adjustment
kw_alpha_div_18S_fungi_habitat_univ_pwt <- as.data.frame(sapply(4:8, function(x) pairwise.t.test(alpha_div_18S_fungi[,x],
                                                                                                 alpha_div_18S_fungi[,2], p.adjust.method = "BH")))
kw_alpha_div_18S_fungi_habitat_univ_pwt <- kw_alpha_div_18S_fungi_habitat_univ_pwt %>%
  rename(Reads = V1,
         Observed = V2,
         Shannon = V3,
         Simpson = V4,
         Evenness = V5)

# Make some adjustments on the table
kw_alpha_div_18S_fungi_habitat_univ_pwt <- t(kw_alpha_div_18S_fungi_habitat_univ_pwt) # transpose
kw_alpha_div_18S_fungi_habitat_univ_pwt <- as_tibble(kw_alpha_div_18S_fungi_habitat_univ_pwt, rownames = "Metric") # adding rownames as a column
class(kw_alpha_div_18S_fungi_habitat_univ_pwt) # checking object class

# Save resulting table with fwrite to avoid and issues with characters
fwrite(kw_alpha_div_18S_fungi_habitat_univ_pwt, "exported_tables/soil/kw_alpha_div_18S_fungi_habitat_univ_pwt.csv")

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
alpha_div_18S_fungi %>%
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
ggsave("results/soil/alpha_div_18S_fungi_site.pdf", width = 8, height = 6, dpi = 150) # save graphic

phy_18S_prune_soil_noncontam_prev05_true_scale_filt_protist <- phy_18S_prune_soil_noncontam_prev05_true_scale_filt %>%
  subset_taxa(D1  == "Amorphea" | D1  == "Archaeplastida" | D1  == "Cryptophyceae" | D1  == "Discoba" |  D1  == "SAR") %>%
  subset_taxa( D4 != "Chlorophyta unk." | D4 != "Chloroplastida unk." | D4 != "Phragmoplastophyta") %>%
  subset_taxa(D5 != "Fungi" & D6 != "Metazoa") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)
get_taxa_unique(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_protist, taxonomic.rank = "D4")

phy_18S_prune_soil_noncontam_prev05_true_scale_filt_protist # 1742 taxa and 45 samples
any(taxa_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_protist) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_protist) == 0) # gives the number of cases
smin <-min(sample_sums(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_protist)) # smin = 105

protist_taxa <- as.data.frame(tax_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_protist))
write.csv(protist_taxa, "exported_tables/soil/protist_taxa.csv")




Eukaryota 

Unassigned
    
    Mitochondria" &
      Order   != "Chloroplast"
  )


# Import fungi matrix from previous lines
#phy_soil_18S_fungi <- read.csv("exported_tables/soil/phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_fungi_tax_fixed.csv")
#unique(phy_soil_18S_fungi$Ordertax)

# Function to collapse a certain number of taxa into category others

# Using all samples here, including those with low abundance.
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9 <- tax_glom(phy_18S_prune_soil_noncontam_prev05_true_scale_filt, "D9")
get_taxa_unique(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9, taxonomic.rank = "D9")
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi <- phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9 %>%
                                                                 subset_taxa(
                                                                   D5 == "Fungi") %>%
                                                                prune_samples(sample_sums(.) > 0, .)
get_taxa_unique(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi, taxonomic.rank = "D9")

# Function to collapse a certain number of taxa into category others
merge_less_than_top_prev05_top20_fungi <- function(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi, top=19){
  transformed <- transform_sample_counts(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:23 for eukaryote data
  }
  return(merged)
}

# The construction of a category others is not working here.
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi_top20 <- merge_less_than_top_prev05_top20_fungi(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi, top=19)
get_taxa_unique(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi_top20, taxonomic.rank = "D9")
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi_top20_df = psmelt(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi_top20) # Melt to long format, it will give a warning messsage because of table headers
unique(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi_top20_df$D9)

write.csv(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi_top20_df,
          "exported_tables/soil/phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi_top20_df.csv")

phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi_top20_df_agr <- read.csv("exported_tables/soil/phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi_top20_df_fixed.csv")
unique(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi_top20_df_agr$D9)
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi_top20_df_agr = aggregate(Abundance~Soil.Type+Habitat+Sample.Site+D0+D1+D2+D3+D4+D5+D6+D7+D8+D9, data=phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi_top20_df_agr, FUN=mean) # add common factors to use for plotting
unique(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi_top20_df_agr$D9)

colors_top23_fungi <- c("#9a6324", "#46f0f0", "#1F78B4", "#aaffc3", "#e6beff",  "#33A02C", "#4363d8", "#008080", "#FB9A99", "#e6194b",
                  "#f032e6", "#bcf60c", "#fabebe","#f58231", "#ffe119", "#6A3D9A", "#fffac8", "#000075", "#800000", "gray", "lightblue", "black", "darkgreen")

phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi_top20_df_agr$D9 <- factor(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi_top20_df_agr$D9,
levels = c("Agaricomycetes", "Blastocladiales", "Chytridiomycetes", "Dothideomycetes",  "Eurotiomycetes",
"Exobasidiomycetes", "Fungi", "Glomerales", "Leotiomycetes", "Mortierellales", "Mucorales",
"Orbiliomycetes", "Pezizomycetes", "Pezizomycotina", "Rhizophydiales", "Saccharomycetes",
"Sordariomycetes", "Spizellomycetales", "Tremellomycetes", "Others"))

phy_18S_prune_soil_fungi_site <- ggplot(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi_top20_df_agr, aes(x = Sample.Site, y = Abundance, fill = D9)) + #plotting by sample
  facet_nested(. ~ Soil.Type+Habitat, scales = "free",
               labeller = labeller(Soil.Type = soil.type.labs, Habitat = habitat.type.labs)) +
  geom_bar(stat = "identity", width = 0.95, position = "fill") + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top23_fungi) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(label = site.labs,
    expand = c(0.0, 0.0),
    drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Site") + # add the title on x axis
  labs(fill = "Class") +
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
phy_18S_prune_soil_fungi_site
ggsave("results/soil/18S_fungi_Soil_Habitat_Site_ShSk_barplot_filt_comp_class_level.pdf", width = 10, height = 5, dpi = 200) # save graphic

phy_18S_prune_soil_fungi_habitat <- ggplot(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D9_fungi_top20_df_agr, aes(x = Habitat, y = Abundance, fill = D9)) + #plotting by sample
  facet_nested(. ~ Soil.Type, scales = "free",
               labeller = labeller(Soil.Type = soil.type.labs, Habitat = habitat.type.labs)) +
  geom_bar(stat = "identity", width = 0.95, position = "fill") + # adds to 100%
  #geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top23_fungi) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  scale_x_discrete(label = habitat.type.labs,
                   expand = c(0.0, 0.0),
                   drop = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  labs(fill = "Class") +
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    size = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) # Format facet grid title
phy_18S_prune_soil_fungi_habitat
ggsave("results/soil/18S_fungi_Soil_Habitat_ShSk_barplot_filt_comp_class_level.pdf", width = 8, height = 5, dpi = 200) # save graphic


phy_soil_18S_major_meta_fungi_combined <- ggarrange(phy_soil_18S_major_groups_site, phy_soil_18S_meta_site_plot,
                                                    phy_18S_prune_soil_fungi_site, legend = "right", ncol = 1, nrow = 3, labels = "AUTO",
                                        align = "v")
phy_soil_18S_major_meta_fungi_combined
ggsave("results/soil/phy_soil_18S_major_meta_fungi_combined.pdf", width = 12, height = 14, dpi = 300)

###################
morph <- read_excel("data/SHSK_18S_morph_Nov_10_2021.xlsx")
summary_morph <- morph %>%
  group_by(Genus, Sample.Site) %>%
  summarise(
    count = n())
write.csv(summary_morph, "exported_tables/morphology/summary_morph.csv")

######## Combine nMDS plots #########
# Needs to use the morphological R script
Figure_nMDS <- ggarrange(morpho_comp_nmds_plot,
                         phy_18S_prune_baermann_noncontam_prev05_true_scale_filt_comp_nmds_plot,
                         phy_18S_prune_soil_noncontam_prev05_true_scale_filt_comp_nmds_plot,
                         common.legend = TRUE,  legend = "right", ncol = 3, nrow = 1, labels = "AUTO")
Figure_nMDS
ggsave("results/morphology/morpho_16S_soil_18S_baemann_18S_soil_comp_nmds.pdf", width = 14, height = 5, dpi = 200)

######################  
# write ASV_Soil matrix to file as csv
ASV_SHSK_18S_soil = as(otu_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt), "matrix")

# Transpose matrix
if(taxa_are_rows(phy_18S_prune_soil_noncontam_prev05_true_scale_filt)) {ASV_SHSK_18S_soil <- t(ASV_SHSK_18S_soil)}

# Coerce to a data frame
ASV_SHSK_18S_soil_df = as.data.frame(ASV_SHSK_18S_soil)
class(ASV_SHSK_18S_soil_df)
# write ASV_df matrix to file as csv
write.csv(ASV_SHSK_18S_soil_df, "exported_tables/soil/ASV_SHSK_18S_soil_df.csv")

# write TAX_soil matrix to file as csv
TAX_SHSK_18S_soil = as(tax_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt), "matrix")
TAX_SHSK_18S_soil_df = as.data.frame(TAX_SHSK_18S_soil)
class(TAX_SHSK_18S_soil_df)
# write TAX_df matrix to file as csv
write.csv(TAX_SHSK_18S_soil_df, "exported_tables/soil/TAX_SHSK_18S_soil_df.csv")

# write TAX_soil matrix to file as csv
SAMPLE_SHSK_18S_soil = as(sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt), "matrix")
SAMPLE_SHSK_18S_soil_df = as.data.frame(SAMPLE_SHSK_18S_soil)
class(SAMPLE_SHSK_18S_soil_df)
# write Sample_df matrix to file as csv
write.csv(SAMPLE_SHSK_18S_soil_df, "exported_tables/soil/SAMPLE_SHSK_18S_soil_df.csv")

################ Export final phyloseq ASV table to a text file for analysis in Primer #####################################
phy_18S_prune_soil_noncontam_prev05_true_scale_filt # phyloseq object to be used for these steps
colnames(tax_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt))
head(tax_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt))

# Collapse tax table according to taxonomic rank
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D22 <- tax_glom(phy_18S_prune_soil_noncontam_prev05_true_scale_filt, taxrank = "D22")
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D21 <- tax_glom(phy_18S_prune_soil_noncontam_prev05_true_scale_filt, taxrank = "D21")
phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_D21 <- tax_glom(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema, taxrank = "D21")

######## Aldex2 on phyoseq objects############
aldex2_all_euk_hab <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt)),
                                phyloseq::sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt)$Habitat,
                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
all_euk_otu_table_hab <- data.frame(phyloseq::otu_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt))
all_euk_otu_table_hab <- rownames_to_column(all_euk_otu_table_hab, var = "OTU")
write.csv(all_euk_otu_table_hab, "aldex2/all_euk_otu_table_hab.csv")

aldex2_D22_hab <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D22)),
                                phyloseq::sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D22)$Habitat,
                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
D22_otu_table_hab <- data.frame(phyloseq::otu_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D22))
D22_otu_table_hab <- rownames_to_column(D22_otu_table_hab, var = "OTU")
write.csv(D22_otu_table_hab, "aldex2/D22_otu_table_hab.csv")

aldex2_D21_hab <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D21)),
                                phyloseq::sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D21)$Habitat,
                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
D21_otu_table_hab <- data.frame(phyloseq::otu_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D21))
D21_otu_table_hab <- rownames_to_column(D21_otu_table_hab, var = "OTU")
write.csv(D21_otu_table_hab, "aldex2/D21_otu_table_hab.csv")

aldex2_D21_soiltype <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D21)),
                                phyloseq::sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D21)$Soil.Type,
                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
D21_otu_table_soiltype <- data.frame(phyloseq::otu_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_D21))
D21_otu_table_soiltype <- rownames_to_column(D21_otu_table_soiltype, var = "OTU")
write.csv(D21_otu_table_soiltype, "aldex2/D21_otu_table_soiltype.csv")

aldex2_nema_hab <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_D21)),
                                phyloseq::sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_D21)$Habitat,
                                mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
nema_otu_table_hab <- data.frame(phyloseq::otu_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_D21))
nema_otu_table_hab <- rownames_to_column(nema_otu_table_hab, var = "OTU")
write.csv(nema_otu_table_hab, "aldex2/nema_otu_table_hab.csv")

aldex2_nema_soil_type <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_D21)),
                                 phyloseq::sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_D21)$Soil.Type,
                                 mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
nema_otu_table_soil_type <- data.frame(phyloseq::otu_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema_D21))
nema_otu_table_soil_type <- rownames_to_column(nema_otu_table_soil_type, var = "OTU")
write.csv(nema_otu_table_soil_type, "aldex2/nema_otu_table_soil_type.csv")

# Write aldex2 results to file as csv
write.csv(aldex2_all_euk_hab, "aldex2/aldex2_all_euk_hab.csv")
write.csv(aldex2_D22_hab, "aldex2/aldex2_D22_hab.csv")
write.csv(aldex2_D21_hab, "aldex2/aldex2_D21_hab.csv")
write.csv(aldex2_D21_soiltype, "aldex2/aldex2_D21_soiltype.csv")
write.csv(aldex2_nema_hab, "aldex2/aldex2_nema_hab.csv")
write.csv(aldex2_nema_soil_type, "aldex2/aldex2_nema_soil_type.csv")

# Import aldex2 results and rename the X variable by OTU
aldex2_all_euk_hab_result <- read.csv("aldex2/aldex2_all_euk_hab.csv")
colnames(aldex2_all_euk_hab_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

aldex2_D22_hab_result <- read.csv("aldex2/aldex2_D22_hab.csv")
colnames(aldex2_D22_hab_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

aldex2_D21_hab_result <- read.csv("aldex2/aldex2_D21_hab.csv")
colnames(aldex2_D21_hab_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")
aldex2_D21_soiltype_result <- read.csv("aldex2/aldex2_D21_soiltype.csv")
colnames(aldex2_D21_soiltype_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

aldex2_nema_hab_result <- read.csv("aldex2/aldex2_nema_hab.csv")
colnames(aldex2_nema_hab_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

#Clean up presentation
taxa_info <- data.frame(tax_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt))
taxa_info <- taxa_info %>%
  rownames_to_column(var = "OTU")
sample_tab_soil <- data.frame(sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt))
sample_soil_factors <- data.frame(cbind(sample_tab_soil[, c(21, 22, 37)]))

taxa_info_nema <- data.frame(tax_table(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema))
taxa_info_nema <- taxa_info_nema %>%
  rownames_to_column(var = "OTU")
sample_tab_soil_nema <- data.frame(sample_data(phy_18S_prune_soil_noncontam_prev05_true_scale_filt_nema))
sample_soil_factors_nema <- data.frame(cbind(sample_tab_soil_nema[, c(21, 22, 37)]))

# filter aldex2 results by sig kw.ep and join the taxanomic information
sig_aldex2_all_euk_hab_result <- aldex2_all_euk_hab_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_all_euk_hab_result <- left_join(sig_aldex2_all_euk_hab_result, taxa_info)
write.csv(sig_aldex2_all_euk_hab_result, "aldex2/sig_aldex2_all_euk_hab_result.csv") # 87 taxa with significant differential abundance

sig_aldex2_D22_hab_result <- aldex2_D22_hab_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_D22_hab_result <- left_join(sig_aldex2_D22_hab_result, taxa_info)
write.csv(sig_aldex2_D22_hab_result, "aldex2/sig_aldex2_D22_hab_result.csv") # 87 taxa with significant differential abundance

sig_aldex2_D21_hab_result <- aldex2_D21_hab_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_D21_hab_result <- left_join(sig_aldex2_D21_hab_result, taxa_info)
write.csv(sig_aldex2_D21_hab_result, "aldex2/sig_aldex2_D21_hab_result.csv") # 87 taxa with significant differential abundance

sig_aldex2_nema_hab_result <- aldex2_nema_hab_result %>%
  #filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)
sig_aldex2_nema_hab_result <- left_join(sig_aldex2_nema_hab_result, taxa_info_nema)
write.csv(sig_aldex2_nema_hab_result, "aldex2/sig_aldex2_nema_hab_result.csv")

# Create clr objects by using OTU ids and original otu table, all significant OTUs and only top 20
sig_aldex2_all_euk_hab_result_count <- left_join(sig_aldex2_all_euk_hab_result, all_euk_otu_table_hab)
write.csv(sig_aldex2_all_euk_hab_result_count, "aldex2/sig_aldex2_all_euk_hab_result_count.csv")
clr_all_euk_hab <- sig_aldex2_all_euk_hab_result_count[, -(2:26)]
rownames(clr_all_euk_hab) <- clr_all_euk_hab$OTU
clr_all_euk_hab <- clr_all_euk_hab[, -1]

sig_aldex2_D22_hab_result_count <- left_join(sig_aldex2_D22_hab_result, D22_otu_table_hab)
write.csv(sig_aldex2_D22_hab_result_count, "aldex2/sig_aldex2_D22_hab_result_count.csv")
clr_D22_hab <- sig_aldex2_D22_hab_result_count[, -(2:26)]
rownames(clr_D22_hab) <- clr_D22_hab$OTU
clr_D22_hab <- clr_D22_hab[, -1]

sig_aldex2_D21_hab_result_count <- left_join(sig_aldex2_D21_hab_result, D21_otu_table_hab)
write.csv(sig_aldex2_D21_hab_result_count, "aldex2/sig_aldex2_D21_hab_result_count.csv")
clr_D21_hab <- sig_aldex2_D21_hab_result_count[, -(2:26)]
rownames(clr_D21_hab) <- clr_D21_hab$OTU
clr_D21_hab <- clr_D21_hab[, -1]

sig_aldex2_nema_hab_result_count <- left_join(sig_aldex2_nema_hab_result, nema_otu_table_hab)
write.csv(sig_aldex2_nema_hab_result_count, "aldex2/sig_aldex2_nema_hab_result_count.csv")
clr_nema_hab <- sig_aldex2_nema_hab_result_count[, -(2:26)]
rownames(clr_nema_hab) <- clr_nema_hab$OTU
clr_nema_hab <- clr_nema_hab[, -1]

# Adjusting zeros on the matrix and applying log transformation
clr_all_euk_hab_czm <- cmultRepl(t(clr_all_euk_hab),  label=0, method="CZM")
clr_all_euk_hab_czm_log <- t(apply(clr_all_euk_hab_czm, 1, function(x){log(x) - mean(log(x))}))

clr_D22_hab_czm <- cmultRepl(t(clr_D22_hab),  label=0, method="CZM")
clr_D22_hab_czm_log <- t(apply(clr_D22_hab_czm, 1, function(x){log(x) - mean(log(x))}))

clr_D21_hab_czm <- cmultRepl(t(clr_D21_hab),  label=0, method="CZM")
clr_D21_hab_czm_log <- t(apply(clr_D21_hab_czm, 1, function(x){log(x) - mean(log(x))}))

clr_nema_hab_czm <- cmultRepl(t(clr_nema_hab),  label=0, method="CZM")
clr_nema_hab_czm_log <- t(apply(clr_nema_hab_czm, 1, function(x){log(x) - mean(log(x))}))

# Define palette color
col_matrix <- colorRampPalette(brewer.pal(10, "BrBG"))(256)
col_matrix2 <- colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_matrix3 <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_matrix4 <- colorRamp2(c(-2, 0, 2), c("blue", "yellow", "red"))

# Combine the heatmap and the annotation
Z.Score.clr_all_euk_hab_czm_log <- scale(t(clr_all_euk_hab_czm_log))
Z.Score.clr_D22_hab_czm_log <- scale(t(clr_D22_hab_czm_log))
Z.Score.clr_D21_hab_czm_log <- scale(t(clr_D21_hab_czm_log))
Z.Score.clr_nema_hab_czm_log <- scale(t(clr_nema_hab_czm_log))

# Define colors for each level of qualitative variables, i.e. soil type and habitats
# Create the heatmap annotation for all eukaryotes plotting D21
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
# all euks - ASV level
Z.Score.clr_all_euk_hab_czm_log_name <- as.data.frame(Z.Score.clr_all_euk_hab_czm_log)
str(Z.Score.clr_all_euk_hab_czm_log_name)
Z.Score.clr_all_euk_hab_czm_log_name <- rownames_to_column(Z.Score.clr_all_euk_hab_czm_log_name, var = "OTU")
Z.Score.clr_all_euk_hab_czm_log_name <- left_join(Z.Score.clr_all_euk_hab_czm_log_name, taxa_info)
head(Z.Score.clr_all_euk_hab_czm_log_name)
Z.Score.clr_all_euk_hab_czm_log_name <- Z.Score.clr_all_euk_hab_czm_log_name[, -(2:46)]

all_euk_otu_table_hab_total <- as.data.frame(all_euk_otu_table_hab)
all_euk_otu_table_hab_total$Total <- rowSums(all_euk_otu_table_hab_total[, -1])
head(all_euk_otu_table_hab_total)
all_euk_otu_table_hab_total <- all_euk_otu_table_hab_total[, -(2:46)]

Z.Score.clr_all_euk_hab_czm_log_total <- left_join(Z.Score.clr_all_euk_hab_czm_log_name, all_euk_otu_table_hab_total)
head(Z.Score.clr_all_euk_hab_czm_log_total)

# D22 level
Z.Score.clr_D22_hab_czm_log_name <- as.data.frame(Z.Score.clr_D22_hab_czm_log)
str(Z.Score.clr_D22_hab_czm_log_name)
Z.Score.clr_D22_hab_czm_log_name <- rownames_to_column(Z.Score.clr_D22_hab_czm_log_name, var = "OTU")
Z.Score.clr_D22_hab_czm_log_name <- left_join(Z.Score.clr_D22_hab_czm_log_name, taxa_info)
head(Z.Score.clr_D22_hab_czm_log_name)
Z.Score.clr_D22_hab_czm_log_name <- Z.Score.clr_D22_hab_czm_log_name[, -(2:46)]

D22_otu_table_hab_total <- as.data.frame(D22_otu_table_hab)
D22_otu_table_hab_total$Total <- rowSums(D22_otu_table_hab_total[, -1])
head(D22_otu_table_hab_total)
D22_otu_table_hab_total <- D22_otu_table_hab_total[, -(2:46)]

Z.Score.clr_D22_hab_czm_log_total <- left_join(Z.Score.clr_D22_hab_czm_log_name, D22_otu_table_hab_total)
head(Z.Score.clr_D22_hab_czm_log_total)

#D21 level
Z.Score.clr_D21_hab_czm_log_name <- as.data.frame(Z.Score.clr_D21_hab_czm_log)
str(Z.Score.clr_D21_hab_czm_log_name)
Z.Score.clr_D21_hab_czm_log_name <- rownames_to_column(Z.Score.clr_D21_hab_czm_log_name, var = "OTU")
Z.Score.clr_D21_hab_czm_log_name <- left_join(Z.Score.clr_D21_hab_czm_log_name, taxa_info)
head(Z.Score.clr_D21_hab_czm_log_name)
Z.Score.clr_D21_hab_czm_log_name <- Z.Score.clr_D21_hab_czm_log_name[, -(2:46)]

D21_otu_table_hab_total <- as.data.frame(D21_otu_table_hab)
D21_otu_table_hab_total$Total <- rowSums(D21_otu_table_hab_total[, -1])
head(D21_otu_table_hab_total)
D21_otu_table_hab_total <- D21_otu_table_hab_total[, -(2:46)]

Z.Score.clr_D21_hab_czm_log_total <- left_join(Z.Score.clr_D21_hab_czm_log_name, D21_otu_table_hab_total)
head(Z.Score.clr_D21_hab_czm_log_total)

# Defining color scheme for row annotations
abundance_col_fun = colorRamp2(c(0, 50000, 100000, 150000),
                               c("#c7eae5",
                                          "#80cdc1",
                                          "#35978f",
                                          "#01665e"))
                                          
ha_right_all_euk = rowAnnotation(
  Abundance = Z.Score.clr_all_euk_hab_czm_log_total$Total, border = FALSE, col = list(Abundance = abundance_col_fun))
row_labels_all_euk = Z.Score.clr_all_euk_hab_czm_log_total$OTU

ha_right_D22 = rowAnnotation(
  Abundance = Z.Score.clr_D22_hab_czm_log_total$Total, border = FALSE, col = list(Abundance = abundance_col_fun))
row_labels_D22 = Z.Score.clr_D22_hab_czm_log_total$D22

ha_right_D21 = rowAnnotation(
  Abundance = Z.Score.clr_D21_hab_czm_log_total$Total, border = FALSE, col = list(Abundance = abundance_col_fun))
row_labels_D21 = Z.Score.clr_D21_hab_czm_log_total$D21

# Plot heatmap at D21 level
hm_D21_hab <- Heatmap(Z.Score.clr_D21_hab_czm_log, name = "Z-score, CLR", col = col_matrix,
                         column_title = "18S rRNA soil microbiome", 
                         column_title_gp = gpar(fontface = "bold", fontsize = 14),
                         column_split = as.vector(as.vector(sample_soil_factors$Soil.Type)),
                         #column_order = order(as.numeric(gsub("column", "", colnames(Z.Score.phylum_hab)))),
                         border = TRUE,
                         top_annotation = ha_top,
                         right_annotation = ha_right_D21,
                         row_title = "Phylum",
                         row_labels = row_labels_D21,
                         row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                         row_names_gp = gpar(fontsize = 6),
                         column_names_gp = gpar(fontsize = 6),
                         row_order = order(row_labels_phy),
                         rect_gp = gpar(col = "white", lwd = 1),
                         show_column_names = TRUE,
                         show_heatmap_legend = TRUE)
hm_D21_hab


# Plot heatmap at D22 level
hm_D22_hab <- Heatmap(Z.Score.clr_D22_hab_czm_log, name = "Z-score, CLR", col = col_matrix,
                      column_title = "18S rRNA soil microbiome", 
                      column_title_gp = gpar(fontface = "bold", fontsize = 14),
                      column_split = as.vector(as.vector(sample_soil_factors$Soil.Type)),
                      #column_order = order(as.numeric(gsub("column", "", colnames(Z.Score.phylum_hab)))),
                      border = TRUE,
                      top_annotation = ha_top,
                      right_annotation = ha_right_D22,
                      row_title = "Phylum",
                      row_labels = row_labels_D22,
                      row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                      row_names_gp = gpar(fontsize = 6),
                      column_names_gp = gpar(fontsize = 6),
                      row_order = order(row_labels_D22),
                      rect_gp = gpar(col = "white", lwd = 1),
                      show_column_names = TRUE,
                      show_heatmap_legend = TRUE)
hm_D22_hab

# all euks - ASV level
hm_all_euk_hab <- Heatmap(Z.Score.clr_all_euk_hab_czm_log, name = "Z-score, CLR", col = col_matrix,
                      column_title = "18S rRNA soil microbiome", 
                      column_title_gp = gpar(fontface = "bold", fontsize = 14),
                      column_split = as.vector(as.vector(sample_soil_factors$Soil.Type)),
                      #column_order = order(as.numeric(gsub("column", "", colnames(Z.Score.phylum_hab)))),
                      border = TRUE,
                      top_annotation = ha_top,
                      right_annotation = ha_right_all_euk,
                      row_title = "Phylum",
                      row_labels = row_labels_all_euk,
                      row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                      row_names_gp = gpar(fontsize = 6),
                      column_names_gp = gpar(fontsize = 6),
                      row_order = order(row_labels_all_euk),
                      rect_gp = gpar(col = "white", lwd = 1),
                      show_column_names = TRUE,
                      show_heatmap_legend = TRUE)
hm_all_euk_hab
# Create the heatmap annotation for all eukaryotes plotting D21
ha_top_nema = HeatmapAnnotation(
  Soil = as.vector(sample_soil_factors_nema$Soil.Type),
  Habitat = as.vector(sample_soil_factors_nema$Habitat),
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

Z.Score.clr_nema_hab_czm_log_name <- as.data.frame(Z.Score.clr_nema_hab_czm_log)
str(Z.Score.clr_nema_hab_czm_log_name)
Z.Score.clr_nema_hab_czm_log_name <- rownames_to_column(Z.Score.clr_nema_hab_czm_log_name, var = "OTU")
Z.Score.clr_nema_hab_czm_log_name <- left_join(Z.Score.clr_nema_hab_czm_log_name, taxa_info_nema)
head(Z.Score.clr_nema_hab_czm_log_name)
Z.Score.clr_nema_hab_czm_log_name <- Z.Score.clr_nema_hab_czm_log_name[, -(2:46)]

nema_otu_table_hab_total <- as.data.frame(nema_otu_table_hab)
nema_otu_table_hab_total$Total <- rowSums(nema_otu_table_hab_total[, -1])
head(nema_otu_table_hab_total)
nema_otu_table_hab_total <- nema_otu_table_hab_total[, -(2:46)]

Z.Score.clr_nema_hab_czm_log_total <- left_join(Z.Score.clr_nema_hab_czm_log_name, nema_otu_table_hab_total)
head(Z.Score.clr_nema_hab_czm_log_total)

# Defining color scheme for row annotations
abundance_col_fun = colorRamp2(c(0, 1000, 2500, 5000),
                               c("#c7eae5",
                                          "#80cdc1",
                                          "#35978f",
                                          "#01665e"))
                                          
ha_right_nema = rowAnnotation(
  Abundance = Z.Score.clr_nema_hab_czm_log_total$Total, border = FALSE, col = list(Abundance = abundance_col_fun))
row_labels_nema = Z.Score.clr_nema_hab_czm_log_total$D22

# Plot heatmap nematode at D22 level
hm_nema_hab <- Heatmap(Z.Score.clr_nema_hab_czm_log, name = "Z-score, CLR", col = col_matrix,
                      column_title = "18S rRNA soil microbiome", 
                      column_title_gp = gpar(fontface = "bold", fontsize = 14),
                      column_split = as.vector(as.vector(sample_soil_factors_nema$Soil.Type)),
                      #column_order = order(as.numeric(gsub("column", "", colnames(Z.Score.phylum_hab)))),
                      border = TRUE,
                      top_annotation = ha_top_nema,
                      right_annotation = ha_right_nema,
                      row_title = "Phylum",
                      row_labels = row_labels_nema,
                      row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                      row_names_gp = gpar(fontsize = 6),
                      column_names_gp = gpar(fontsize = 6),
                      row_order = order(row_labels_nema),
                      rect_gp = gpar(col = "white", lwd = 1),
                      show_column_names = TRUE,
                      show_heatmap_legend = TRUE)
hm_nema_hab

# Import matrix, adjust for zeros, and transform data using clr
shsk_18S_nema_czm <- read.table("aldex2/ASV_SHSK_18S_soil_nema_df.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, comment.char="")
rownames(shsk_18S_nema_czm) <- shsk_18S_nema_czm$ASV
shsk_18S_nema_czm <- shsk_18S_nema_czm[, -1]
shsk_18S_nema_czm <- cmultRepl(t(shsk_18S_nema_czm),  label=0, method="CZM")
shsk_18S_nema_clr <- t(apply(shsk_18S_nema_czm, 1, function(x){log(x) - mean(log(x))}))


shsk_18S_euk_otu_count <- read.table("aldex2/ASV_SHSK_18S_soil_euk_df.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, comment.char="")
shsk_18S_euk_otu_sig <- read.csv("aldex2/shsk_soil_18S_aldex_sig_euk.csv")
shsk_18S_euk_otu_sig_count <- left_join(shsk_18S_euk_otu_sig, shsk_18S_euk_otu_count)
rownames(shsk_18S_euk_otu_sig_count) <- shsk_18S_euk_otu_sig_count$ASV

shsk_18S_euk_czm <- shsk_18S_euk_otu_sig_count[, -(1:27)]
shsk_18S_euk_czm <- cmultRepl(t(shsk_18S_euk_czm),  label=0, method="CZM")
shsk_18S_euk_clr <- t(apply(shsk_18S_euk_czm, 1, function(x){log(x) - mean(log(x))}))

shsk_18S_euk_tax <- read.table("aldex2/TAX_SHSK_18S_soil_euk_df.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, comment.char="")



# Organize total abundance and taxa name
Z.Score.clr_D21_hab_czm_log_name <- as.data.frame(Z.Score.clr_D21_hab_czm_log)
str(Z.Score.clr_D21_hab_czm_log_name)
Z.Score.clr_D21_hab_czm_log_name <- rownames_to_column(Z.Score.clr_D21_hab_czm_log_name, var = "OTU")
Z.Score.clr_D21_hab_czm_log_name <- left_join(Z.Score.clr_D21_hab_czm_log_name, taxa_info)
head(Z.Score.clr_D21_hab_czm_log_name)
Z.Score.clr_D21_hab_czm_log_name <- Z.Score.clr_D21_hab_czm_log_name[, -(2:46)]

D21_otu_table_hab_total <- as.data.frame(D21_otu_table_hab)
D21_otu_table_hab_total$Total <- rowSums(D21_otu_table_hab_total[, -1])
head(D21_otu_table_hab_total)
D21_otu_table_hab_total <- D21_otu_table_hab_total[, -(2:46)]

Z.Score.clr_D21_hab_czm_log_total <- left_join(Z.Score.clr_D21_hab_czm_log_name, D21_otu_table_hab_total)
head(Z.Score.clr_D21_hab_czm_log_total)

# Defining color scheme for row annotations
abundance_col_fun = colorRamp2(c(0, 50000, 100000, 150000),
                               c("#c7eae5",
                                          "#80cdc1",
                                          "#35978f",
                                          "#01665e"))
                                          

# Checking first lines of object, transpose it, and then create a heatmap according to the tax rank
head(shsk_18S_nema_clr)
shsk_18S_nema_clr_trav <- t(shsk_18S_nema_clr)
shsk_18S_nema_clr_trav[, order(colnames(shsk_18S_nema_clr_trav))]
heatmap(shsk_18S_nema_clr_trav, scale = "none", col = bluered(100))

head(shsk_18S_euk_clr)
shsk_18S_euk_clr_trav <- t(shsk_18S_euk_clr)
shsk_18S_euk_clr_trav[, order(colnames(shsk_18S_euk_clr_trav))]
heatmap(shsk_18S_euk_clr_trav, scale = "none", col = bluered(100))


shsk_18S_euk_count_total <- left_join(shsk_18S_euk_otu_sig, shsk_18S_euk_otu_count)
shsk_18S_euk_count_total$Total <- rowSums(shsk_18S_euk_count_total[, -(1:27)])
shsk_18S_euk_count_total <- shsk_18S_euk_count_total[, -(3:72)]
shsk_18S_euk_count_total <- shsk_18S_euk_count_total[, -(1)]

shsk_18S_euk_count_total <- left_join(shsk_18S_euk_tax, shsk_18S_euk_count_total)
head(shsk_18S_euk_count_total)

ha_right_nema = rowAnnotation(
  Abundance = shsk_18S_nema_count_total$Total)
row_labels_nema = shsk_18S_nema_count_total$ASV

ha_right_euk = rowAnnotation(
  Abundance = shsk_18S_euk_count_total$Total)
row_labels_euk = shsk_18S_euk_count_total$ASV

# Combine the heatmap and the annotation
col_matrix3 <- colorRamp2(c(-2, 0, 2), c("blue", "#EEEEEE", "red"), space = "RGB")

Z.Score.shsk_18S_nema <- scale(shsk_18S_nema_clr_trav)

Heatmap(Z.Score.shsk_18S_nema, col = col_matrix3)
hm_shsk_18S_nema <- Heatmap(Z.Score.shsk_18S_nema, name = "Z-score, CLR", col = col_matrix3,
                                         column_title = "18S rRNA soil nematodes", 
                                         column_title_gp = gpar(fontface = "bold", fontsize = 14),
                                         top_annotation = ha_top,
                                         right_annotation = ha_right_nema,
                                         row_title = "ASV",
                                         row_labels = row_labels_phy,
                                         row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                                         row_names_gp = gpar(fontsize = 6),
                                         column_names_gp = gpar(fontsize = 6),
                                         #row_order = order(row_labels_nema),
                                         rect_gp = gpar(col = "white", lwd = 1),
                                         show_column_names = TRUE,
                                         show_heatmap_legend = TRUE)
hm_shsk_18S_nema

Z.Score.shsk_18S_euk <- scale(shsk_18S_euk_clr_trav)
Heatmap(Z.Score.shsk_18S_euk, col = col_matrix3)

hm_shsk_18S_euk <- Heatmap(Z.Score.shsk_18S_euk, name = "Z-score, CLR", col = col_matrix3,
                            column_title = "18S rRNA soil eukaryotes", 
                            column_title_gp = gpar(fontface = "bold", fontsize = 14),
                            top_annotation = ha_top,
                            right_annotation = ha_right_euk,
                            row_title = "ASV",
                            row_labels = row_labels_euk,
                            row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                            row_names_gp = gpar(fontsize = 6),
                            column_names_gp = gpar(fontsize = 6),
                            row_order = order(row_labels_euk),
                            rect_gp = gpar(col = "white", lwd = 1),
                            show_column_names = TRUE,
                            show_heatmap_legend = TRUE)
hm_shsk_18S_euk

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

shsk_ssp_czm <- cmultRepl(t(clr_ssp),  label=0, method="CZM")
shsk_ssp_clr <- t(apply(shsk_ssp_czm, 1, function(x){log(x) - mean(log(x))}))

shsk_ssp_czm_top20 <- cmultRepl(t(clr_ssp_top20),  label=0, method="CZM")
shsk_ssp_clr_top20 <- t(apply(shsk_ssp_czm_top20, 1, function(x){log(x) - mean(log(x))}))

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

shsk_ssp_clr_trav <- t(shsk_ssp_clr)
shsk_ssp_clr[, order(colnames(shsk_ssp_clr))]
heatmap(shsk_ssp_clr_trav, scale = "none", col = bluered(100))

shsk_ssp_clr_trav_top20 <- t(shsk_ssp_clr_top20)
shsk_ssp_clr_top20[, order(colnames(shsk_ssp_clr_top20))]
heatmap(shsk_ssp_clr_trav_top20, scale = "none", col = bluered(100))

# Define palette color
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)

# Define colors for each levels of qualitative variables, i.e. soil type and habitat
col = list(Soil = c("Compact.Soil" = "blue", "Not.Compact.Soil" = "lightblue"),
           Habitat = c("Coastal.scrub" = "#1B9E77", "Riparian" = "#D95F02", "Hollyleaf.cherry" = "#7570B3",
                       "Oak.wood" = "#E7298A", "Native.grass" = "#66A61E", "Chaparral" = "#E6AB02"))

# Create the heatmap annotation
ha <- HeatmapAnnotation(
  Soil = phyloseq::sample_data(phy_soil_prune_noncontam_comb05_true_scale_filt)$Soil.Type,
  Habitat = phyloseq::sample_data(phy_soil_prune_noncontam_comb05_true_scale_filt)$Habitat)

# Combine the heatmap and the annotation
Z.Score.phylum <- scale(shsk_phylum_clr_trav)

Z.Score.class <- scale(shsk_class_clr_trav)

Z.Score.order <- scale(shsk_order_clr_trav)
Z.Score.order.top20 <- scale(shsk_order_clr_trav_top20)

Z.Score.fam <- scale(shsk_fam_clr_trav)
Z.Score.fam.top20 <- scale(shsk_fam_clr_trav_top20)

Z.Score.genus <- scale(shsk_genus_clr_trav)
Z.Score.genus.top20 <- scale(shsk_genus_clr_trav_top20)

Z.Score.ssp <- scale(shsk_ssp_clr_trav)
Z.Score.ssp.top20 <- scale(shsk_ssp_clr_trav_top20)

hm_phylum <- Heatmap(Z.Score.phylum, name = "Z-score, CLR",
                     top_annotation = ha,
                     row_title = "Phylum",
                     row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                     row_names_gp = gpar(fontsize = 10),
                     show_column_names = FALSE,
                     show_heatmap_legend = TRUE)
hm_phylum

hm_class <- Heatmap(Z.Score.class, name = "Z-score, CLR",
                     top_annotation = ha,
                     row_title = "Class",
                     row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                     row_names_gp = gpar(fontsize = 10),
                     show_column_names = FALSE,
                     show_heatmap_legend = TRUE)
hm_class

hm_order <- Heatmap(Z.Score.order, name = "Z-score, CLR",
                    top_annotation = ha,
                    row_title = "Order",
                    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                    row_names_gp = gpar(fontsize = 10),
                    show_column_names = FALSE,
                    show_heatmap_legend = TRUE)
hm_order

hm_fam <- Heatmap(Z.Score.fam, name = "Z-score, CLR",
                     top_annotation = ha,
                     row_title = "Family",
                     row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                     row_names_gp = gpar(fontsize = 10),
                     show_column_names = FALSE,
                     show_heatmap_legend = TRUE)
hm_fam

hm_genus <- Heatmap(Z.Score.genus, name = "Z-score, CLR",
                     top_annotation = ha,
                     row_title = "Genus",
                     row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                     row_names_gp = gpar(fontsize = 10),
                     show_column_names = FALSE,
                     show_heatmap_legend = TRUE)
hm_genus

hm_ssp <- Heatmap(Z.Score.ssp, name = "Z-score, CLR",
                    top_annotation = ha,
                    row_title = "Species",
                    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                    row_names_gp = gpar(fontsize = 5),
                    show_column_names = FALSE,
                    show_heatmap_legend = TRUE)
hm_ssp

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

Z.Score.ssp.df <- as.data.frame(t(Z.Score.ssp))
Z.Score.ssp.df.top20 <- as.data.frame(t(Z.Score.ssp.top20))

# Use rownames to create another variable called sample
Z.Score.phylum.df <- tibble::rownames_to_column(Z.Score.phylum.df, "Sample")
Z.Score.class.df <- tibble::rownames_to_column(Z.Score.class.df, "Sample")
Z.Score.order.df <- tibble::rownames_to_column(Z.Score.order.df, "Sample")
Z.Score.order.df.top20 <- tibble::rownames_to_column(Z.Score.order.df.top20, "Sample")
Z.Score.fam.df <- tibble::rownames_to_column(Z.Score.fam.df, "Sample")
Z.Score.fam.df.top20 <- tibble::rownames_to_column(Z.Score.fam.df.top20, "Sample")
Z.Score.genus.df <- tibble::rownames_to_column(Z.Score.genus.df, "Sample")
Z.Score.genus.df.top20 <- tibble::rownames_to_column(Z.Score.genus.df.top20, "Sample")
Z.Score.ssp.df <- tibble::rownames_to_column(Z.Score.ssp.df, "Sample")
Z.Score.ssp.df.top20 <- tibble::rownames_to_column(Z.Score.ssp.df.top20, "Sample")

#Add sample factors (Sample.Site, Habitat, Soil.Type) to the Z.score object
Z.Score.phylum.df.var <- data.frame(cbind(Z.Score.phylum.df, sample_soil_factors))
Z.Score.class.df.var <- data.frame(cbind(Z.Score.class.df, sample_soil_factors))
Z.Score.order.df.var <- data.frame(cbind(Z.Score.order.df, sample_soil_factors))
Z.Score.order.df.var.top20 <- data.frame(cbind(Z.Score.order.df.top20, sample_soil_factors))
Z.Score.fam.df.var <- data.frame(cbind(Z.Score.fam.df, sample_soil_factors))
Z.Score.fam.df.var.top20 <- data.frame(cbind(Z.Score.fam.df.top20, sample_soil_factors))
Z.Score.genus.df.var <- data.frame(cbind(Z.Score.genus.df, sample_soil_factors))
Z.Score.genus.df.var.top20 <- data.frame(cbind(Z.Score.genus.df.top20, sample_soil_factors))
Z.Score.ssp.df.var <- data.frame(cbind(Z.Score.ssp.df, sample_soil_factors))
Z.Score.ssp.df.var.top20 <- data.frame(cbind(Z.Score.ssp.df.top20, sample_soil_factors))

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


Z.Score.ssp.df.long <- pivot_longer(data = Z.Score.ssp.df.var,
                                       cols = -c(Sample, Sample.Type, Sample.Site, Habitat, Soil.Type),
                                       names_to = "ASV",
                                       values_to = "Z")
Z.Score.ssp.df.long <- as.data.frame(Z.Score.ssp.df.long)
Z.Score.ssp.df.long_merged <- merge(Z.Score.ssp.df.long, sig_aldex2_ssp_hab_result, by.x = 6, by.y = 1, all.x = TRUE)

Z.Score.ssp.df.long.top20 <- pivot_longer(data = Z.Score.ssp.df.var.top20,
                                            cols = -c(Sample, Sample.Type, Sample.Site, Habitat, Soil.Type),
                                            names_to = "ASV",
                                            values_to = "Z")
Z.Score.ssp.df.long.top20 <- as.data.frame(Z.Score.ssp.df.long.top20)
Z.Score.ssp.df.long_merged.top20 <- merge(Z.Score.ssp.df.long.top20, sig_aldex2_ssp_hab_result_top20, by.x = 6, by.y = 1, all.x = TRUE)

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
  scale_y_discrete(expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_phy
ggsave("16S_heatmap_phy_Soil_Habitat_greater_500_filt_comb.pdf", width = 10, height = 1.5, dpi = 150) # save graphic

Z.Score.class.df.long_merged$Habitat <- factor(Z.Score.class.df.long_merged$Habitat,
                                                levels = c("Coastal.scrub", "Native.grass", "Chaparral", "Riparian", "Hollyleaf.cherry", "Oak.wood"))
heatmap_class <- ggplot(Z.Score.class.df.long_merged, aes(x = Sample, y = Class, fill = Z)) + 
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
  scale_y_discrete(expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_class
ggsave("16S_heatmap_class_Soil_Habitat_greater_500_filt_comb.pdf", width = 10, height = 4, dpi = 150) # save graphic
  
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
  scale_y_discrete(expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_order
ggsave("16S_heatmap_order_top20__Soil_Habitat_greater_500_filt_comb.pdf", width = 10, height = 6, dpi = 150) # save graphic


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
  scale_y_discrete(expand = c(0, 0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(panel.spacing.x = unit(0, "lines")) +
  labs(fill = "Z-score, CLR")
heatmap_fam
  
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
    scale_y_discrete(expand = c(0, 0)) +
    theme(legend.position="right") +
    theme(legend.title = element_text(face = "bold")) +
    theme(legend.title.align = 0.5) +
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
    theme(plot.margin = unit(c(0,0,0,0), "lines")) +
    theme(panel.spacing.x = unit(0, "lines")) +
    labs(fill = "Z-score, CLR")
heatmap_genus
ggsave("16S_heatmap_genus_Soil_Habitat_greater_500_filt_comb_top20.pdf", width = 12, height = 4, dpi = 150) # save graphic
  
Figure_heatmap_taxa <- ggarrange(heatmap_phy,  heatmap_class, heatmap_order, heatmap_fam, heatmap_genus, heights = unit(c(6, 17, 20, 20, 16), "lines"), ncol = 1, nrow = 5, align = "v", common.legend = TRUE, legend = "right") 
Figure_heatmap_taxa
ggsave("16S_heatmap_all_levels_Soil_Habitat_greater_500_filt_comb.pdf", width = 12, height = 20, dpi = 150) # save graphic

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

# Melt to long format (for ggplot2); prune out V2 below 1% in each sample
phy_soil_prune_noncontam_comb05_true_scale_v2 <- phy_soil_prune_noncontam_comb05_true_scale %>%
  tax_glom(taxrank = "V2") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(V2)                                      # Sort data frame alphabetically by phylum
phy_soil_prune_noncontam_comb05_true_scale_v2
str(phy_soil_prune_noncontam_comb05_true_scale_v2)
class(phy_soil_prune_noncontam_comb05_true_scale_v2)

# Melt to long format (for ggplot2); prune out V3 below 1% in each sample
phy_soil_prune_noncontam_comb05_true_scale_v3 <- phy_soil_prune_noncontam_comb05_true_scale %>%
  tax_glom(taxrank = "V3") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(V3)                                      # Sort data frame alphabetically by phylum
phy_soil_prune_noncontam_comb05_true_scale_v3
str(phy_soil_prune_noncontam_comb05_true_scale_v3)
class(phy_soil_prune_noncontam_comb05_true_scale_v3)

# Melt to long format (for ggplot2); prune out V4 below 1% in each sample
phy_soil_prune_noncontam_comb05_true_scale_v4 <- phy_soil_prune_noncontam_comb05_true_scale %>%
  tax_glom(taxrank = "V4") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(V4)                                      # Sort data frame alphabetically by phylum
phy_soil_prune_noncontam_comb05_true_scale_v4
str(phy_soil_prune_noncontam_comb05_true_scale_v4)
class(phy_soil_prune_noncontam_comb05_true_scale_v4)

# Melt to long format (for ggplot2); prune out V5 below 1% in each sample
phy_soil_prune_noncontam_comb05_true_scale_v5 <- phy_soil_prune_noncontam_comb05_true_scale %>%
  tax_glom(taxrank = "V5") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(V5)                                      # Sort data frame alphabetically by phylum
phy_soil_prune_noncontam_comb05_true_scale_v5
str(phy_soil_prune_noncontam_comb05_true_scale_v5)
class(phy_soil_prune_noncontam_comb05_true_scale_v5)

# Melt to long format (for ggplot2); prune out V6 below 1% in each sample
phy_soil_prune_noncontam_comb05_true_scale_v6 <- phy_soil_prune_noncontam_comb05_true_scale %>%
  tax_glom(taxrank = "V6") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(V6)                                      # Sort data frame alphabetically by phylum
phy_soil_prune_noncontam_comb05_true_scale_v6
str(phy_soil_prune_noncontam_comb05_true_scale_v6)
class(phy_soil_prune_noncontam_comb05_true_scale_v6)

# Melt to long format (for ggplot2); prune out V7 below 1% in each sample
phy_soil_prune_noncontam_comb05_true_scale_v7 <- phy_soil_prune_noncontam_comb05_true_scale %>%
  tax_glom(taxrank = "V7") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(V7)                                      # Sort data frame alphabetically by phylum
phy_soil_prune_noncontam_comb05_true_scale_v7
str(phy_soil_prune_noncontam_comb05_true_scale_v7)
class(phy_soil_prune_noncontam_comb05_true_scale_v7)

# Melt to long format (for ggplot2); prune out V8 below 1% in each sample
phy_soil_prune_noncontam_comb05_true_scale_v8 <- phy_soil_prune_noncontam_comb05_true_scale %>%
  tax_glom(taxrank = "V8") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(V8)                                      # Sort data frame alphabetically by phylum
phy_soil_prune_noncontam_comb05_true_scale_v8
str(phy_soil_prune_noncontam_comb05_true_scale_v8)
class(phy_soil_prune_noncontam_comb05_true_scale_v8)

# Melt to long format (for ggplot2); prune out V9 below 1% in each sample
phy_soil_prune_noncontam_comb05_true_scale_v9 <- phy_soil_prune_noncontam_comb05_true_scale %>%
  tax_glom(taxrank = "V9") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(V9)                                      # Sort data frame alphabetically by phylum
phy_soil_prune_noncontam_comb05_true_scale_v9
str(phy_soil_prune_noncontam_comb05_true_scale_v9)
class(phy_soil_prune_noncontam_comb05_true_scale_v9)

#Get unique taxa for a specific taxonomic rank
unique(phy_soil_prune_noncontam_comb05_true_scale_v2$V2)
unique(phy_soil_prune_noncontam_comb05_true_scale_v3$V3)
unique(phy_soil_prune_noncontam_comb05_true_scale_v4$V4)
unique(phy_soil_prune_noncontam_comb05_true_scale_v5$V5)
unique(phy_soil_prune_noncontam_comb05_true_scale_v6$V6)
unique(phy_soil_prune_noncontam_comb05_true_scale_v7$V7)
unique(phy_soil_prune_noncontam_comb05_true_scale_v8$V8)
unique(phy_soil_prune_noncontam_comb05_true_scale_v9$V9)

# Define the number of colors you want for each taxonomic level
nb.cols.v2 <- 5
soil_colors_v2 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.v2)

nb.cols.v3 <- 8
soil_colors_v3 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.v3)

nb.cols.v4 <- 11
soil_colors_v4 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.v4)

nb.cols.v5 <- 19
soil_colors_v5 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.v5)

nb.cols.v6 <- 15
soil_colors_v6 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.v6)

nb.cols.v7 <- 16
soil_colors_v7 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.v7)

nb.cols.v8 <- 29
soil_colors_v8 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.v8)

nb.cols.v9 <- 53
soil_colors_v9 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.v9)

#Change facet habitat labels
soil.type.labs <- c("Compacted Soil","Uncompacted Soil")
names(soil.type.labs) <- c("Compact.Soil","Not.Compact.Soil")

# Plot bar chart 18S Soil samples per habitat - V2 level
ggplot(phy_soil_prune_noncontam_comb05_true_scale_v2, aes(x = Habitat, y = Abundance, fill = V2)) + #plotting by sample
  facet_wrap(. ~ Soil.Type, scales = "free_x", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = soil_colors_v2) +
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
  ylab("Relative Abundance (V2 > 1%)") + # add the title on y axis
  xlab("Sample Type") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
ggsave("18S_V2_Soil_Habitat_ShSk_barplot_greater_500_filt_comb.pdf", width = 6, height = 4, dpi = 150) # save graphic

# Plot bar chart 18S Soil samples per habitat - V3 level
ggplot(phy_soil_prune_noncontam_comb05_true_scale_v3, aes(x = Habitat, y = Abundance, fill = V3)) + #plotting by sample
  facet_wrap(. ~ Soil.Type, scales = "free_x", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = soil_colors_v3) +
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
  ylab("Relative Abundance (V3 > 1%)") + # add the title on y axis
  xlab("Sample Type") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  theme(strip.text.x = element_text(
      size = 12, color = "black", face = "bold"))
  ggsave("18S_V3_Soil_Habitat_ShSk_barplot_greater_500_filt_comb.pdf", width = 6, height = 4, dpi = 150) # save graphic

  # Plot bar chart 18S Soil samples per habitat - V4 level
  ggplot(phy_soil_prune_noncontam_comb05_true_scale_v4, aes(x = Habitat, y = Abundance, fill = V4)) + #plotting by sample
    facet_wrap(. ~ Soil.Type, scales = "free_x", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
    geom_bar(stat = "identity", position = "fill") + # adds to 100%
    scale_fill_manual(values = soil_colors_v4) +
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
    ylab("Relative Abundance (V4 > 1%)") + # add the title on y axis
    xlab("Sample Type") + # add the title on x axis
    #ggtitle("Soil Bacterial tCommunities at Shipley Skinner") + # add the title on graphic
    theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
    theme(strip.background =element_blank()) + # remove the background of titles
    theme(strip.text.x = element_text(
      size = 12, color = "black", face = "bold"))
  ggsave("18S_V4_Soil_Habitat_ShSk_barplot_greater_500_filt_comb.pdf", width = 6, height = 4, dpi = 150) # save graphic  
  
  # Plot bar chart 18S Soil samples per habitat - V5 level
  ggplot(phy_soil_prune_noncontam_comb05_true_scale_v5, aes(x = Habitat, y = Abundance, fill = V5)) + #plotting by sample
    facet_wrap(. ~ Soil.Type, scales = "free_x", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
    geom_bar(stat = "identity", position = "fill") + # adds to 100%
    scale_fill_manual(values = soil_colors_v5) +
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
    ylab("Relative Abundance (V5 > 1%)") + # add the title on y axis
    xlab("Sample Type") + # add the title on x axis
    #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
    theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
    theme(strip.background =element_blank()) + # remove the background of titles
    theme(strip.text.x = element_text(
      size = 12, color = "black", face = "bold"))
  ggsave("18S_V5_Soil_Habitat_ShSk_barplot_greater_500_filt_comb.pdf", width = 6, height = 5, dpi = 150) # save graphic  
  
  # Plot bar chart 18S Soil samples per habitat - V6 level
  ggplot(phy_soil_prune_noncontam_comb05_true_scale_v6 , aes(x = Habitat, y = Abundance, fill = V6)) + #plotting by sample
    facet_wrap(. ~ Soil.Type, scales = "free_x", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
    geom_bar(stat = "identity", position = "fill") + # adds to 100%
    scale_fill_manual(values = soil_colors_v6) +
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
    ylab("Relative Abundance (V6 > 1%)") + # add the title on y axis
    xlab("Sample Type") + # add the title on x axis
    #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
    theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
    theme(strip.background =element_blank()) + # remove the background of titles
    theme(strip.text.x = element_text(
      size = 12, color = "black", face = "bold"))
  ggsave("18S_V6_Soil_Habitat_ShSk_barplot_greater_500_filt_comb.pdf", width = 6, height = 5, dpi = 150) # save graphic  

  # Plot bar chart 18S Soil samples per habitat - V7 level
  ggplot(phy_soil_prune_noncontam_comb05_true_scale_v7, aes(x = Habitat, y = Abundance, fill = V7)) + #plotting by sample
    facet_wrap(. ~ Soil.Type, scales = "free_x", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
    geom_bar(stat = "identity", position = "fill") + # adds to 100%
    scale_fill_manual(values = soil_colors_v7) +
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
    ylab("Relative Abundance (V7 > 1%)") + # add the title on y axis
    xlab("Sample Type") + # add the title on x axis
    #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
    theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
    theme(strip.background =element_blank()) + # remove the background of titles
    theme(strip.text.x = element_text(
      size = 12, color = "black", face = "bold"))
  ggsave("18S_V7_Soil_Habitat_ShSk_barplot_greater_500_filt_comb.pdf", width = 6, height = 5, dpi = 150) # save graphic  
  
# Filter out Metazoan only,
# For barplots
phy_soil_prune_noncontam_comb05_true_metazoan <- phy_soil_prune_noncontam_comb05_true %>%
  subset_taxa(
    V1  == "Eukaryota" &
      V2  == "Opisthokonta" &
       V3  == "Holozoa" &
        V4  == "Metazoa (Animalia)") %>%
psmelt()
phy_soil_prune_noncontam_comb05_true_metazoan
unique(phy_soil_prune_noncontam_comb05_true_metazoan$V7)
head(phy_soil_prune_noncontam_comb05_true_metazoan)

nb.cols.metazoan <- 13
soil_colors_metazoan <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.metazoan)

#Change facet habitat labels
soil.type.labs <- c("Compacted Soil","Uncompacted Soil")
names(soil.type.labs) <- c("Compact.Soil","Not.Compact.Soil")

# Plot bar chart 18S Soil samples per habitat - Metazoan only
ggplot(phy_soil_prune_noncontam_comb05_true_metazoan, aes(x = Habitat, y = Abundance, fill = V7)) + #plotting by sample
  facet_wrap(. ~ Soil.Type, scales = "free_x", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = soil_colors_metazoan) +
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
  ylab("Relative Abundance, Metazoan Phyla") + # add the title on y axis
  xlab("Sample Type") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
ggsave("18S_Metazoan_Soil_Habitat_ShSk_barplot_filt_comb.pdf", width = 5, height = 4, dpi = 150) # save graphic


# Filter out Fungi only,
# For barplots
phy_soil_prune_noncontam_comb05_true_fungi <- phy_soil_prune_noncontam_comb05_true %>%
  subset_taxa(
    V1  == "Eukaryota" &
      V2  == "Opisthokonta" &
      V3  == "Nucletmycea" &
      V4  == "Fungi") %>%
  psmelt()
phy_soil_prune_noncontam_comb05_true_fungi
unique(phy_soil_prune_noncontam_comb05_true_fungi$V7)
head(phy_soil_prune_noncontam_comb05_true_fungi)

nb.cols.fungi <- 17
soil_colors_fungi <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.fungi)

#Change facet habitat labels
soil.type.labs <- c("Compacted Soil","Uncompacted Soil")
names(soil.type.labs) <- c("Compact.Soil","Not.Compact.Soil")

# Plot bar chart 18S Soil samples per habitat - Fungi only
ggplot(phy_soil_prune_noncontam_comb05_true_fungi, aes(x = Habitat, y = Abundance, fill = V7)) + #plotting by sample
  facet_wrap(. ~ Soil.Type, scales = "free_x", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = soil_colors_fungi) +
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
  ylab("Relative Abundance, Metazoan Phyla") + # add the title on y axis
  xlab("Sample Type") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
ggsave("18S_Fungi_Soil_Habitat_ShSk_barplot_filt_comb.pdf", width = 6, height = 5, dpi = 150) # save graphic

# Filter out Fungi only,
# For barplots
phy_soil_prune_noncontam_comb05_true_pezi <- phy_soil_prune_noncontam_comb05_true %>%
  subset_taxa(
    V1  == "Eukaryota" &
      V2  == "Opisthokonta" &
      V3  == "Nucletmycea" &
      V4  == "Fungi" &
      V5  == "Dikarya" &
      V6  == "Ascomycota" &
      V7  == "Pezizomycotina") %>%
  psmelt()
unique(phy_soil_prune_noncontam_comb05_true_pezi$V9)

nb.cols.pezi <- 30
soil_colors_pezi <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.pezi)

#Change facet habitat labels
soil.type.labs <- c("Compacted Soil","Uncompacted Soil")
names(soil.type.labs) <- c("Compact.Soil","Not.Compact.Soil")

# Plot bar chart 18S Soil samples per habitat - Fungi only
ggplot(phy_soil_prune_noncontam_comb05_true_pezi, aes(x = Habitat, y = Abundance, fill = V9)) + #plotting by sample
  facet_wrap(. ~ Soil.Type, scales = "free_x", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = soil_colors_pezi) +
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
  ylab("Relative Abundance, Metazoan Phyla") + # add the title on y axis
  xlab("Sample Type") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
ggsave("18S_Pezizomycotina_Soil_Habitat_ShSk_barplot_filt_comb.pdf", width = 7, height = 5, dpi = 150) # save graphic

# Filter out Nematoda only,
# For barplots
phy_soil_prune_noncontam_comb05_true_nema <- phy_soil_prune_noncontam_comb05_true %>%
  subset_taxa(
    V1  == "Eukaryota" &
      V2  == "Opisthokonta" &
      V3  == "Holozoa" &
      V4  == "Metazoa (Animalia)" &
      V5 == "Eumetazoa" &
      V6 == "Bilateria" &
      V7 == "Nematoda") %>%
psmelt()

phy_soil_prune_noncontam_comb05_true_nema
unique(phy_soil_prune_noncontam_comb05_true_nema$V9)

nb.cols.nema <- 7
soil_colors_nema <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols.nema)

#Change facet habitat labels
soil.type.labs <- c("Compacted Soil","Uncompacted Soil")
names(soil.type.labs) <- c("Compact.Soil","Not.Compact.Soil")

# Plot bar chart 18S Soil samples per habitat - Nematoda only
ggplot(phy_soil_prune_noncontam_comb05_true_nema, aes(x = Habitat, y = Abundance, fill = V9)) + #plotting by sample
  facet_wrap(. ~ Soil.Type, scales = "free_x", labeller = labeller(Soil.Type = soil.type.labs)) + # By Habitat
  geom_bar(stat = "identity", position = "fill") + # adds to 100%
  scale_fill_manual(values = soil_colors_nema) +
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
  ylab("Relative Abundance, Nematoda") + # add the title on y axis
  xlab("Sample Type") + # add the title on x axis
  #ggtitle("Soil Bacterial Communities at Shipley Skinner") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold"))
ggsave("18S_Nematoda_Soil_Habitat_ShSk_barplot_filt_comb.pdf", width = 5, height = 4, dpi = 150) # save graphic

phy_soil_prune_noncontam_comb05_true_nema_df <- as.data.frame(phy_soil_prune_noncontam_comb05_true_nema)
write.csv(phy_soil_prune_noncontam_comb05_true_nema, "phy_soil_prune_noncontam_comb05_true_nema.csv")

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
ggsave("16S_Class_Top20_Soil_Habitat_ShSk_barplot_greater_500_filt_comb.pdf", width = 7, height = 4, dpi = 150) # save graphic  

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
ggsave("16S_Class_Top20_Soil_Habitat_ShSk_barplot_greater_500_filt_comb.pdf", width = 7, height = 4, dpi = 150) # save graphic  

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
  ggsave("16S_Soil_ShSk_ASV_PCoA_filt_freq.pdf", width = 8, height = 6, dpi = 300)

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
  ggsave("16S_Soil_ShSk_ASV_PCoA_filt_freq_soil-worm_soil-type.pdf", width = 8, height = 6, dpi = 300)

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
  ggsave("16S_Nema_ShSk_ASV_NMDS_greater_500.pdf", width = 12, height = 8, dpi = 300)

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
  ggsave("16S_Phylum_Soil_Nema_ShSk_barplot.pdf", width = 10, height = 4, dpi = 300) # save graphic

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
  ggsave("16S_Class_Soil_Nema_ShSk_barplot.pdf", width = 12, height = 4, dpi = 300) # save graphic

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
  ggsave("16S_Phylum_Soil_ShSk.pdf", width = 8, height = 4, dpi = 300) # save graphic

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
  ggsave("16S_Phylum_Soil-Type_ShSk.pdf", width = 8, height = 4, dpi = 300) # save graphic

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
  ggsave("16S_Class_Soil-Habitat_ShSk.pdf", width = 8, height = 4, dpi = 300) # save graphic

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
  ggsave("16S_Class_Soil-Type_ShSk.pdf", width = 8, height = 4, dpi = 300) # save graphic

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
  ggsave("16S_Soil_ShSk_ASV_NMDS_greater_500.pdf", width = 6, height = 4, dpi = 300)

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
        ggsave("16S_Soil_ShSk_ASV_NMDS_filt_freq.pdf", width = 6, height = 4, dpi = 300)

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


####################################### Only Worm samples #######################################
# Keep only samples from the single worm dataset
phy_worm <- subset_samples(phy, Experiment == "Microbiome")
phy_worm

# Check if there are samples with zero reads
smin <- min(sample_sums(phy_worm))
phy_worm_prune <- prune_samples(sample_sums(phy_worm)>0, phy_worm)
smin <- min(sample_sums(phy_worm_prune))
phy_worm_prune

# Check for ASVs that have no counted reads
any(taxa_sums(phy_worm_prune) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_worm_prune) == 0) # gives the number of cases
phy_worm_prune <-prune_taxa(taxa_sums(phy_worm_prune) > 0, phy_worm_prune)
phy_worm_prune

#Inspect Library sizes
df_phy_worm_prune <- as.data.frame(sample_data(phy_worm_prune))
df_phy_worm_prune$LibrarySize <- sample_sums(phy_worm_prune)
df_phy_worm_prune <- df_phy_worm_prune[order(df_phy_worm_prune$LibrarySize),]
df_phy_worm_prune$Index <- seq(nrow(df_phy_worm_prune))
ggplot(data = df_phy_worm_prune, aes(x= Index, y= LibrarySize, color = Sample.Type)) +
  geom_point()
ggsave("18S_Soil_library-size-sample-type-worm-samples.pdf", width = 8, height = 6, dpi = 150)

# Identify Contaminants using package decontam - Frequency method
worm_contamdf_freq <- isContaminant(phy_worm_prune, method="frequency", conc="QubitConc")
head(worm_contamdf_freq)
tail(worm_contamdf_freq)
table(worm_contamdf_freq$contaminant)
head(which(worm_contamdf_freq$contaminant))

# Explore AVSs that are potentially contaminants, ASV abundance should be inversely correlated with DNA contration
plot_frequency(phy_worm_prune, taxa_names(phy_worm_prune)[c(12, 117, 279, 283)], conc="QubitConc") + 
  xlab("DNA Concentration (ng/ul)")

# Explore AVSs that are potentially contaminants, randomly choosing 4 AVSs.
set.seed(100)
plot_frequency(phy_worm_prune, taxa_names(phy_worm_prune)[sample(which(worm_contamdf_freq$contaminant),4)], conc="QubitConc") + 
  xlab("DNA Concentration (ng/ul)")

# Identify Contaminants using package decontam - Prevalence method
sample_data(phy_worm_prune)$is.neg <- sample_data(phy_worm_prune)$Sample.Control == "Control.Sample"
worm_contamdf_prev <- isContaminant(phy_worm_prune, method="prevalence", neg="is.neg")
table(worm_contamdf_prev$contaminant)
head(which(worm_contamdf_prev$contaminant))

# Identify Contaminants using package decontam - Prevalence method, increase threshhold value to 0.5
worm_contamdf_prev05 <- isContaminant(phy_worm_prune, method="prevalence", neg="is.neg", threshold=0.5)
table(worm_contamdf_prev05$contaminant)

# Make phyloseq object of presence-absence in negative controls and true samples
phy_worm_prune_pa <- transform_sample_counts(phy_worm_prune, function(abund) 1*(abund>0))
phy_worm_prune_pa_neg <- prune_samples(sample_data(phy_worm_prune_pa)$Sample.Control == "Control.Sample", phy_worm_prune_pa)
phy_worm_prune_pa_pos <- prune_samples(sample_data(phy_worm_prune_pa)$Sample.Control == "True.Sample", phy_worm_prune_pa)

# Make data.frame of prevalence in positive and negative samples
df_phy_worm_prune_pa <- data.frame(pa.pos=taxa_sums(phy_worm_prune_pa_pos), pa.neg=taxa_sums(phy_worm_prune_pa_neg),
                                   contaminant=worm_contamdf_prev$contaminant)
ggplot(data=df_phy_worm_prune_pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Make data.frame of prevalence in positive and negative samples, threshold of 0.5
df_phy_worm_prune_pa05 <- data.frame(pa.pos=taxa_sums(phy_worm_prune_pa_pos), pa.neg=taxa_sums(phy_worm_prune_pa_neg),
                                     contaminant=worm_contamdf_prev05$contaminant)
ggplot(data=df_phy_worm_prune_pa05, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Identify Contaminants using package decontam -  Combined method (i.e., Frequency and Prevalence)
worm_contamdf_comb <- isContaminant(phy_worm_prune, conc="QubitConc", neg="is.neg", threshold=0.1, detailed = TRUE, normalize = TRUE, method="combined")
table(worm_contamdf_comb$contaminant)
head(which(worm_contamdf_comb$contaminant))

worm_contamdf_comb05 <- isContaminant(phy_worm_prune, conc="QubitConc", neg="is.neg", threshold=0.5, detailed = TRUE, normalize = TRUE, method="combined")
table(worm_contamdf_comb05$contaminant)
head(which(worm_contamdf_comb05$contaminant))

# Remove contaminants from phyloseq object - combined05 based - check for samples with zero reads
phy_worm_prune_noncontam_comb05 <- prune_taxa(!worm_contamdf_comb05$contaminant, phy_worm_prune)
phy_worm_prune_noncontam_comb05
smin <-min(sample_sums(phy_worm_prune_noncontam_comb05))
any(sample_sums(phy_worm_prune_noncontam_comb05) == 0) # if TRUE, there are samples with no reads, otherwise FALSE
sum(sample_sums(phy_worm_prune_noncontam_comb05) == 0) # gives the number of cases
phy_worm_prune_noncontam_comb05 <-prune_samples(sample_sums(phy_worm_prune_noncontam_comb05) > 0, phy_worm_prune_noncontam_comb05)

# Check for ASVs that have no counted reads
any(taxa_sums(phy_worm_prune_noncontam_comb05) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_worm_prune_noncontam_comb05) == 0) # gives the number of cases
phy_worm_prune_noncontam_comb05 <-prune_taxa(taxa_sums(phy_worm_prune_noncontam_comb05) > 0, phy_worm_prune_noncontam_comb05)
phy_worm_prune_noncontam_comb05

#Subset by sample types, that is, only worm samples removing controls and blanks
phy_worm_prune_noncontam_comb05_true <- subset_samples(phy_worm_prune_noncontam_comb05, Sample.Type =="SingleWorm") %>%
  prune_taxa(taxa_sums(.) > 0, .)
phy_worm_prune_noncontam_comb05_true
any(taxa_sums(phy_worm_prune_noncontam_comb05_true) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_worm_prune_noncontam_comb05_true) == 0) # gives the number of cases

#Checking the total number of reads and its distribution across samples
worm_readsums_df = data.frame(nreads = sort(taxa_sums(phy_worm_prune_noncontam_comb05_true), TRUE), sorted = 1:ntaxa(phy_worm_prune_noncontam_comb05_true), 
                              type = "ASVs")
worm_readsums_df = rbind(worm_readsums_df, data.frame(nreads = sort(sample_sums(phy_worm_prune_noncontam_comb05_true), 
                                                                    TRUE), sorted = 1:nsamples(phy_worm_prune_noncontam_comb05_true), type = "Samples"))
title = "Total number of reads"
p = ggplot(worm_readsums_df, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
ggsave("18S_worm_ShSk_ASV_Read_Profile_filt_comb05-worm-samples.pdf", width = 12, height = 6, dpi = 150)

# Scale reads to even sequencing depth by removing samples with a total number of reads < 500
phy_worm_prune_noncontam_comb05_true_scale <- prune_samples(sample_sums(phy_worm_prune_noncontam_comb05_true)>500, phy_worm_prune_noncontam_comb05_true)
smin <- min(sample_sums(phy_worm_prune_noncontam_comb05_true_scale))
phy_worm_prune_noncontam_comb05_true_scale
any(taxa_sums(phy_worm_prune_noncontam_comb05_true_scale) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_worm_prune_noncontam_comb05_true_scale) == 0) # gives the number of cases
phy_worm_prune_noncontam_comb05_true_scale <-prune_taxa(taxa_sums(phy_worm_prune_noncontam_comb05_true_scale) > 0, phy_worm_prune_noncontam_comb05_true_scale)
phy_worm_prune_noncontam_comb05_true_scale

# Transform number of reads to relative abundance
phy_worm_prune_noncontam_comb05_true_scale
phy_worm_prune_noncontam_comb05_true_scale_trans = transform_sample_counts(phy_worm_prune_noncontam_comb05_true_scale, function(x) 100 * x/sum(x))

# Check for ASVs that have no counted reads
any(taxa_sums(phy_worm_prune_noncontam_comb05_true_scale_trans) == 0) # if TRUE, there are ASVs with no reads, otherwise FALSE
sum(taxa_sums(phy_worm_prune_noncontam_comb05_true_scale_trans) == 0) # gives the number of cases

# Remove ASVs that have no reads
phy_worm_prune_noncontam_comb05_true_scale_trans <- prune_taxa(taxa_sums(phy_worm_prune_noncontam_comb05_true_scale_trans) > 0, phy_worm_prune_noncontam_comb05_true_scale_trans)
phy_worm_prune_noncontam_comb05_true_scale_trans

# Ordinate NMDS
set.seed(1)
phy_worm_prune_noncontam_comb05_true_scale_trans_nmds <- ordinate(
  physeq = phy_worm_prune_noncontam_comb05_true_scale_trans, 
  method = "NMDS", 
  distance = "bray"
)

# Plot NMDS
unique(sample_data(phy_worm_prune_noncontam_comb05_true_scale_trans)$Habitat)

theme_set(theme_bw())
phyloseq::plot_ordination(
  physeq = phy_worm_prune_noncontam_comb05_true_scale_trans,
  ordination = phy_worm_prune_noncontam_comb05_true_scale_trans_nmds,
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
  #ggtitle("Shipley Skinner bacterial Communities - ASV level") + # add the title on graphic
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +
  theme(strip.background =element_blank()) + # remove the background of titles
  #annotate("text", x = -0.50, y = 0.51, label = "B", size = 9, fontface = 2) +
  annotate("text", x = 120, y = 120, label ="2D Stress: 0.001") +
  ggsave("18S_worm_ShSk_ASV_NMDS_greater_500_filt_comb.pdf", width = 6, height = 4, dpi = 150)

# write ASV_Soil matrix to file as csv
ASV_SHSK_18S_worm = as(otu_table(phy_worm_prune_noncontam_comb05_true_scale), "matrix")

# Transpose matrix
if(taxa_are_rows(phy_worm_prune_noncontam_comb05_true_scale)) {ASV_SHSK_18S_worm <- t(ASV_SHSK_18S_worm)}

# Coerce to a data frame
ASV_SHSK_18S_worm_df = as.data.frame(ASV_SHSK_18S_worm)
class(ASV_SHSK_18S_worm_df)
# write ASV_df matrix to file as csv
write.csv(ASV_SHSK_18S_worm_df, "/Users/tiagopereira/Dropbox/SK-data/18S-soil-single-nematode/R-codes/ASV_SHSK_18S_worm_df.csv")

# write TAX_worm matrix to file as csv
TAX_SHSK_18S_worm = as(tax_table(phy_worm_prune_noncontam_comb05_true_scale), "matrix")
TAX_SHSK_18S_worm_df = as.data.frame(TAX_SHSK_18S_worm)
class(TAX_SHSK_18S_worm_df)
# write TAX_df matrix to file as csv
write.csv(TAX_SHSK_18S_worm_df, "/Users/tiagopereira/Dropbox/SK-data/18S-soil-single-nematode/R-codes/TAX_SHSK_18S_worm_df.csv")

# write TAX_worm matrix to file as csv
SAMPLE_SHSK_18S_worm = as(sample_data(phy_worm_prune_noncontam_comb05_true_scale), "matrix")
SAMPLE_SHSK_18S_worm_df = as.data.frame(SAMPLE_SHSK_18S_worm)
class(SAMPLE_SHSK_18S_worm_df)
# write Sample_df matrix to file as csv
write.csv(SAMPLE_SHSK_18S_worm_df, "/Users/tiagopereira/Dropbox/SK-data/18S-soil-single-nematode/R-codes/SAMPLE_SHSK_18S_worm_df.csv")
