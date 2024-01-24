setwd("/Users/tiagopereira/Documents/shipley-skinner/Abiotic-data")

library("devtools") #Allows the installation of R packages
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
#library("dplyr") #A tool for working/manipulating dataframes. Sometimes it causes an error if using functions having plyr.
library("vegan") #Community ecology package; it can be use for diverse multivariate analysis
library("stringr") #Provides a set of functions to work with strings
library("tidyr") #Provides tools to work with tables, dealing with columns and rows as well as individual cells. 
library("ggord") #A package for creating ordination plots with ggplot2
library("ggfortify") #Unified plotting tools for common statistics and graphics through ggplot2
library("ggrepel") #Provides geoms for ggplot2 to repel overlapping text labels
library("factoextra") #Provides functions for multivariate data analysis as well as graphics using ggplot2
library("ggbiplot") #Implements the biplot (i.e. samples data points and vectors of variables) using ggplot2
library("ggh4x") #ggplot extension

source("summary.R")

theme_set(theme_bw())

# Import matrices using excel files
# Entire matrix (i.e., all variables and factors)
sk_abiotic <- read_excel("abiotic.xlsx", sheet = "abiotic")

#Import matrices per variable
sk_soil <- read_excel("abiotic.xlsx", sheet = "soil") # soil texture
sk_organic <- read_excel("abiotic.xlsx", sheet = "organic") # organic material
sk_nutrients <- read_excel("abiotic.xlsx", sheet = "nutrients") # essential nutrients
sk_ph_BD <- read_excel("abiotic.xlsx", sheet = "pH-BD") # pH and bulk density measures

# Import alpha diversity results from both 16S and 18S soil dataset
alpha_div_comb_soil <- read_excel("alpha_div_16S_18S_soil.xlsx") # alpha diversity matrix

# Calculate summary (i.e., mean and standard error) using summarySE function from summary.R file
# Summaries per habitat and soil type
sk_soil_sum_habitat <- summarySE(sk_soil, measurevar=c("Perc"), groupvars=c("SoilType", "Habitat","Code", "Class"))
sk_soil_sum_habitat
write.csv(sk_soil_sum_habitat, "sk_soil_sum_habitat.csv")

sk_organic_sum_habitat <- summarySE(sk_organic, measurevar=c("Perc"), groupvars=c("SoilType", "Habitat","Code", "Type"))
sk_organic_sum_habitat
write.csv(sk_organic_sum_habitat, "sk_organic_sum_habitat.csv")

sk_nutrients_sum_habitat <- summarySE(sk_nutrients, measurevar=c("ppm"), groupvars=c("SoilType", "Habitat","Code", "Type"))
sk_nutrients_sum_habitat

sk_ph_BD_sum_habitat <- summarySE(sk_ph_BD, measurevar=c("pH_BD"), groupvars=c("SoilType", "Habitat", "Code", "Type"))
sk_ph_BD_sum_habitat

# Summaries per soil type
sk_soil_sum_soiltype <- summarySE(sk_soil, measurevar=c("Perc"), groupvars=c("SoilType", "Class"))
sk_soil_sum_soiltype

sk_organic_sum_soiltype <- summarySE(sk_organic, measurevar=c("Perc"), groupvars=c("SoilType", "Type"))
sk_organic_sum_soiltype

sk_nutrients_sum_soiltype <- summarySE(sk_nutrients, measurevar=c("ppm"), groupvars=c("SoilType", "Type"))
sk_nutrients_sum_soiltype

sk_ph_BD_sum_soiltype <- summarySE(sk_ph_BD, measurevar=c("pH_BD"), groupvars=c("SoilType", "Type"))
sk_ph_BD_sum_soiltype


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

environ_var_habitat_1 <- ggarrange(sk_soil_sum_habitat_plot, sk_nutrients_sum_habitat_plot,
                                 ncol = 2, nrow = 1, legend = "none",
                                 align = "hv", labels = "AUTO")
environ_var_habitat_1
ggsave("FigureS1.pdf", width = 8, height =5, dpi = 150) # save graphic

sk_organic_sum_habitat_plot <- ggplot(sk_organic_sum_habitat, aes(x= Code, y= Perc, fill = Type)) +
  geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=Perc-se, ymax=Perc+se), size = 0.5, width = 0.2, position=position_dodge(width = 0.9)) +
  geom_point(size = 2, shape = 21, position=position_dodge(width = 0.9)) +
  scale_fill_brewer(palette = "Greens") +
  facet_grid(Type ~ SoilType, scales = "free") +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  theme(strip.text.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.x = element_text(face = "bold", size = 16)) + # adjusts the title of x axis
  theme(axis.title.y = element_blank()) + # adjusts the title of x axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12, color = "black")) + # adjusts text of x axis
  scale_y_continuous(expand = c(0.02, 0.02)) +
  xlab("Habitat") + # add the title on y axis
  ggtitle("Organic Material (%)") +
  theme(plot.title = element_text(face = "bold", size = 14, colour = "black", hjust = 0.5)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(1, "lines")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
sk_organic_sum_habitat_plot

sk_ph_BD_sum_habitat_plot <- ggplot(sk_ph_BD_sum_habitat, aes(x= Code, y= pH_BD, fill = Type)) +
  geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=pH_BD-se, ymax=pH_BD+se), size = 0.5, width = 0.2, position=position_dodge(width = 0.9)) +
  geom_point(size = 2, shape = 21, position=position_dodge(width = 0.9)) +
  scale_fill_brewer(palette = "PuOr") +
  facet_grid(Type ~ SoilType, scales = "free") +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  theme(strip.text.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.x = element_text(face = "bold", size = 14)) + # adjusts the title of x axis
  theme(axis.title.y = element_blank()) + # adjusts the title of x axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12, color = "black")) + # adjusts text of x axis
  scale_y_continuous(expand = c(0.02, 0.02)) +
  xlab("Habitat") + # add the title on y axis
  ggtitle(expression(bold("pH & BD"~(g/cm^3)))) +
  theme(plot.title = element_text(face = "bold", size = 14, colour = "black", hjust = 0.5)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(1, "lines")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
sk_ph_BD_sum_habitat_plot

environ_var_habitat_2 <- ggarrange(sk_organic_sum_habitat_plot, sk_ph_BD_sum_habitat_plot,
                                  ncol = 2, nrow = 1, legend = "none",
                                  align = "hv", labels = "AUTO")
environ_var_habitat_2
ggsave("FigureS2.pdf", width = 8, height =5, dpi = 150) # save graphic

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

environ_var_soiltype1 <- ggarrange(sk_soil_sum_soiltype_plot, sk_nutrients_sum_soiltype_plot,
                                  ncol = 2, nrow = 1, legend = "none",
                                  align = "hv", labels = "AUTO")
environ_var_soiltype1
ggsave("FigureS3.pdf", width = 8, height =5, dpi = 150) # save graphic

sk_organic_sum_soiltype_plot <- ggplot(sk_organic_sum_soiltype, aes(x= SoilType, y= Perc, fill = Type)) +
  geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=Perc-se, ymax=Perc+se), size = 0.5, width = 0.2, position=position_dodge(width = 0.9)) +
  geom_point(size = 2, shape = 21, position=position_dodge(width = 0.9)) +
  scale_fill_brewer(palette = "Greens") +
  facet_grid(Type ~ SoilType, scales = "free") +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  theme(strip.text.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  scale_y_continuous(expand = c(0.02, 0.02)) +
  theme(legend.position="right") +
  ggtitle("Organic Material (%)") +
  theme(plot.title = element_text(face = "bold", size = 14, colour = "black", hjust = 0.5)) +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(1, "lines")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
sk_organic_sum_soiltype_plot

sk_ph_BD_sum_soiltype_plot <- ggplot(sk_ph_BD_sum_soiltype, aes(x= SoilType, y= pH_BD, fill = Type)) +
  geom_bar(position=position_dodge(), stat = "identity", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin=pH_BD-se, ymax=pH_BD+se), size = 0.5, width = 0.2, position=position_dodge(width = 0.9)) +
  geom_point(size = 2, shape = 21, position=position_dodge(width = 0.9)) +
  scale_fill_brewer(palette = "PuOr") +
  facet_grid(Type ~ SoilType, scales = "free") +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  theme(strip.text.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  scale_y_continuous(expand = c(0.02, 0.02)) +
  xlab("Soil Type") + # add the title on y axis
  theme(legend.position="right") +
  ggtitle(expression(bold("pH & BD"~(g/cm^3)))) +
  theme(plot.title = element_text(face = "bold", size = 14, colour = "black", hjust = 0.5)) +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(1, "lines")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
sk_ph_BD_sum_soiltype_plot


environ_var_soiltype2 <- ggarrange(sk_organic_sum_soiltype_plot, sk_ph_BD_sum_soiltype_plot,
                                   ncol = 2, nrow = 1, legend = "none",
                                   align = "hv", labels = "AUTO")
environ_var_soiltype2
ggsave("FigureS4.pdf", width = 8, height =5, dpi = 150) # save graphic

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
  scale_color_manual(values = c("#52392F", "#83643E","#C29B6C", 
                                "#397A4C", "#77C063", "#BEDB92"),
                     name = "Habitat",
                     breaks = c("Chaparral", "Coastal_sage_scrub", "Native_grass", "Hollyleaf_Cherry", "Oak_woodland", "Riparian"),
                     labels = c("CHA", "CSS", "NGR", "HLC", "OWL", "RIP"))+
  scale_shape_manual(values = c(19,17),
                     name = "Soil Type",
                     breaks = c("Compacted", "Uncompacted"),
                     labels = c("Compacted", "Uncompacted"))+
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(aes(color =  Habitat, shape = SoilType), alpha = 0.7, size = 4) +
  ylab("PC2 (25.5%)") + # add the title on y axis
  xlab("PC1 (55.7%)") + # add the title on x axis
  theme(legend.position="right") +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.title.align = 0.5) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # removes the gridlines
ggsave("PCA_soiltype_habitat.pdf", width = 6, height =4) # save graphic as pdf for final adjustments using Inkscape/Illustrator

#Non-parametric one-way ANOVA - Kruskal-Wallis rank sum test by habitat and soil type

# KW analysis on all environmental variables at once
# Remember, numerical variables are from columns 6-15. The test is by soiltype, column 4
kw_env_soiltype <- as.data.frame(sapply(6:15, function(x) kruskal.test(sk_abiotic[,x],
                                                                       sk_abiotic[,4])))

kt_pH_soiltype <- kruskal.test(pH ~ SoilType, data = sk_abiotic) # KW-test by soil type
kt_pH_habitat <- kruskal.test(pH ~ Code, data = sk_abiotic) # KW-test by habitat
parwise_pH_habitat <- pairwise.wilcox.test(sk_abiotic$pH, sk_abiotic$Code,
                                        p.adjust.method = "BH")
write.table(parwise_pH_habitat[["p.value"]], file="parwise_pH_habitat.csv", sep=",")

kt_N_soiltype <- kruskal.test(N ~ SoilType, data = sk_abiotic) # KW-test by soil type
kt_N_habitat <- kruskal.test(N ~ Code, data = sk_abiotic) # KW-test by habitat
parwise_N_habitat <- pairwise.wilcox.test(sk_abiotic$N, sk_abiotic$Code,
                                           p.adjust.method = "bonferroni")
write.table(parwise_N_habitat[["p.value"]], file="parwise_N_habitat.csv", sep=",")

kt_P_soiltype <- kruskal.test(P ~ SoilType, data = sk_abiotic) # KW-test by soil type
kt_P_habitat <- kruskal.test(P ~ Code, data = sk_abiotic) # KW-test by habitat
parwise_P_habitat <- pairwise.wilcox.test(sk_abiotic$P, sk_abiotic$Code,
                                          p.adjust.method = "bonferroni")
write.table(parwise_P_habitat[["p.value"]], file="parwise_P_habitat.csv", sep=",")

kt_K_soiltype <- kruskal.test(K ~ SoilType, data = sk_abiotic) # KW-test by soil type
kt_K_habitat <- kruskal.test(K ~ Code, data = sk_abiotic) # KW-test by habitat
parwise_K_habitat <- pairwise.wilcox.test(sk_abiotic$K, sk_abiotic$Code,
                                           p.adjust.method = "bonferroni")
write.table(parwise_K_habitat[["p.value"]], file="parwise_K_habitat.csv", sep=",")

kt_OM_soiltype <- kruskal.test(OM ~ SoilType, data = sk_abiotic) # KW-test by soil type
kt_OM_habitat <- kruskal.test(OM ~ Code, data = sk_abiotic) # KW-test by habitat
parwise_OM_habitat <- pairwise.wilcox.test(sk_abiotic$OM, sk_abiotic$Code,
                                           p.adjust.method = "bonferroni")
write.table(parwise_OM_habitat[["p.value"]], file="parwise_OM_habitat.csv", sep=",")

kt_OC_soiltype <- kruskal.test(OC ~ SoilType, data = sk_abiotic) # KW-test by soil type
kt_OC_habitat <- kruskal.test(OC ~ Code, data = sk_abiotic) # KW-test by habitat
parwise_OC_habitat <- pairwise.wilcox.test(sk_abiotic$OC, sk_abiotic$Code,
                                           p.adjust.method = "bonferroni")
write.table(parwise_OC_habitat[["p.value"]], file="parwise_OC_habitat.csv", sep=",")

kt_Sand_soiltype <- kruskal.test(Sand ~ SoilType, data = sk_abiotic) # KW-test by soil type
kt_Sand_habitat <- kruskal.test(Sand ~ Code, data = sk_abiotic) # KW-test by habitat
parwise_Sand_habitat <- pairwise.wilcox.test(sk_abiotic$Sand, sk_abiotic$Code,
                                           p.adjust.method = "bonferroni")
write.table(parwise_Sand_habitat[["p.value"]], file="parwise_Sand_habitat.csv", sep=",")

kt_Silt_soiltype <- kruskal.test(Silt ~ SoilType, data = sk_abiotic) # KW-test by soil type
kt_Silt_habitat <- kruskal.test(Silt ~ Code, data = sk_abiotic) # KW-test by habitat
parwise_Silt_habitat <- pairwise.wilcox.test(sk_abiotic$Silt, sk_abiotic$Code,
                                           p.adjust.method = "bonferroni")
write.table(parwise_Silt_habitat[["p.value"]], file="parwise_Silt_habitat.csv", sep=",")

kt_Clay_soiltype <- kruskal.test(Clay ~ SoilType, data = sk_abiotic) # KW-test by soil type
kt_Clay_habitat <- kruskal.test(Clay ~ Code, data = sk_abiotic) # KW-test by habitat
parwise_Clay_habitat <- pairwise.wilcox.test(sk_abiotic$Clay, sk_abiotic$Code,
                                           p.adjust.method = "bonferroni")
write.table(parwise_Clay_habitat[["p.value"]], file="parwise_Clay_habitat.csv", sep=",")

kt_BD_soiltype <- kruskal.test(BD ~ SoilType, data = sk_abiotic) # KW-test by soil type
kt_BD_habitat <- kruskal.test(BD ~ Code, data = sk_abiotic) # KW-test by habitat
parwise_BD_habitat <- pairwise.wilcox.test(sk_abiotic$BD, sk_abiotic$Code,
                                           p.adjust.method = "bonferroni")
write.table(parwise_BD_habitat[["p.value"]], file="parwise_BD_habitat.csv", sep=",")

# Change names of factor for alpha diversity plot
alpha_div_comb_soil$Soil.Type <- factor(alpha_div_comb_soil$Soil.Type,
                                        level = c("Compact.Soil", "Not.Compact.Soil"),
                                        labels = c("Compacted", "Uncompacted")) 

alpha_div_comb_soil$Habitat <- factor(alpha_div_comb_soil$Habitat,
                                     level = c("Chaparral", "Coastal.scrub", "Native.grass",
                                               "Hollyleaf.cherry", "Oak.wood", "Riparian"),
                                     labels = c("CHA", "CSS", "NGR", "HLC",  "OWL", "RIP"))

alpha_color_habitat <- c("#52392F", "#83643E","#C29B6C", "#397A4C", "#77C063", "#BEDB92")

alpha_div_comb_soil %>%
  gather(key = metric, value = value, c("Reads", "Observed", "Shannon", "Simpson", "Evenness")) %>%
  mutate(metric = factor(metric, levels = c("Reads", "Observed", "Shannon", "Simpson", "Evenness"))) %>%
  ggplot(aes(x = Habitat, y = value, fill = Gene)) +
  geom_boxplot(outlier.color = NA, width = 0.8) +
  #geom_jitter(aes(color = Habitat), height = 0, width = 0.15) +
  geom_point(aes(colour = Habitat), shape = 21, position=position_jitterdodge()) +
  facet_nested(metric ~ Soil.Type, scales = "free") +
  scale_color_manual(values = alpha_color_habitat) +
  scale_fill_manual(values = c("lightgray", "white")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold")) + # adjusts text of y axis
  #theme(axis.ticks.x = element_blank()) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  #scale_y_continuous(expand = c(0.1, 0.1)) + # plot as % and removes the internal margins
  theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  ylab("") + # add the title on y axis
  xlab("Habitat") + # add the title on x axis
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) # adjusts the title of the legend
ggsave("alpha_div_comb_soil.pdf", width = 8, height = 6, dpi = 300) # save graphic
