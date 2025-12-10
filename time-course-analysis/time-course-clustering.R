# Purpose:
# Identification of dominant temporal expression trajectories in
# longitudinal proteomics data using fuzzy clustering.
#
# Overview:
# This script demonstrates an exploratory time-course analysis aimed at
# capturing coordinated temporal behaviour across proteins rather than
# maximising discrete group separation. Fuzzy clustering (Mfuzz) is used
# to identify dominant expression trajectories across predefined timepoints.
#
# Downstream analyses focus on:
#   - Characterising major temporal patterns
#   - Functional enrichment of high-membership proteins
#   - Network-level validation of cluster coherence
#
# Key principles illustrated:
#   - Time is treated as a biological variable, not a categorical contrast
#   - Cluster number is chosen to reflect dominant qualitative patterns
#   - Functional and network analyses are used to support interpretation,
#     not to infer causality

###### Theme ######
# Set up consistent theme (Helvetica, 6 pt base text)
theme_ajpcell <- function() {
  theme_classic(base_size = 6, base_family = "Helvetica") +
    theme(
      axis.title.x = element_text(size = 10, margin = margin(t = 4)),
      axis.title.y = element_text(size = 10, margin = margin(r = 4)),
      axis.text = element_text(size = 8),
      legend.position = "top",
      legend.title = element_blank()
    )
}
######################################################
library(plyr)
library(dplyr)
library(ggplot2)
library(STRINGdb) #package ‘STRINGdb’ was built under R version 4.3.1 
library(qvalue)


# abd represents data frame of Control and Exp conditions within each animal at each time point.
# experiment was unilateral, 1 leg = stim, the other is internal control
# abundance data quantified across 0, 2, 10, 20, and 30 days of loading

abd$Time<- revalue(abd$Time, c("P0"="D00",
                               "P10"="D10",
                               "P20"="D20",
                               "P30"="D30"))

# Check / set correct format of each column
str(abd)
# make a list of columns that should be factors
factor.cols <- c(3,4,5)
abd[factor.cols] <- lapply(abd[factor.cols], factor)
head(abd)

# Create within animal ID (unite Time and Replicate) to allow for repeated-measures identifier
abd$Animal <- paste(abd$Time,'-',abd$n)
head(abd)


# Run complete 2-way ANOVA approach to extract (main effects and interaction effects)
ME.Time <- ddply(abd, ~Accession, function(x) summary(aov(ug.protein ~ Time*Condition+Error(Animal),data=x))[[1]][[1]][[5]][1])
colnames(ME.Time) <- c("Accession", "ME.Time")
ME.Condition<- ddply(abd, ~Accession, function(x) summary(aov(ug.protein ~ Time*Condition+Error(Animal),data=x))[[2]][[1]][[5]][1])
colnames(ME.Condition) <- c("Accession", "ME.Condition")
Interaction <- ddply(abd, ~Accession, function(x)  summary(aov(ug.protein ~ Time*Condition+Error(Animal),data=x))[[2]][[1]][[5]][2])
colnames(Interaction) <- c("Accession", "Interaction")
p.vals <- merge(ME.Time, ME.Condition, by = c(1))
p.vals <- merge(p.vals, Interaction, by = c(1))
head(p.vals)

output <- ddply(abd, c("Accession", "Time", "Condition"), summarise, mean = mean(ug.protein),
                                                                      sd = sd(ug.protein))
head(output)
output.wide <- spread(output, Condition, mean)

output.wide <- output %>%
  pivot_wider(
    names_from = Condition,
    values_from = c(mean, sd),
    names_sep = "_"
  )

head(output.wide)
stats.output <- merge(output.wide, p.vals, by= "Accession")
stats.output <- stats.output[order(stats.output$Interaction),]

# Use the Qvalue package to calculate FDR and adjusted p-values to correct for multiple comparisons
stats.output$q.val.Interaction <- qvalue(stats.output$Interaction)$qvalues
stats.output$BH.Interaction <- qvalue(stats.output$Interaction, pi0=1)$qvalues	# N.B. BH correction not q values !
stats.output$FC <- stats.output$mean_Stim/stats.output$mean_Ctrl
stats.output$log.FC <- log(stats.output$FC,2) ##Log2
head(stats.output)
stats.output<- stats.output[order(stats.output$BH.Interaction),]


# Extract significatn interaction at BH-adjusted P value < 0.05
BH.05 <- subset(stats.output, BH.Interaction <= 0.05)
head(BH.05)
print(length(unique(BH.05$Accession))) # 187 proteins

head(stats.output)

write.csv(stats.output, "2_Way_ANOVA_Stats_Output.csv", row.names = F)


# Subset proteins that only exhibit a significant interaction across the experiment
sig <- subset(stats.output, Interaction<=0.05)


# collect raw data from the abd dataframe for the proteins that are of statistical interest
bp.data <- subset(abd, Accession %in% sig$Accession)
head(bp.data)


## Quick view of all boxplots of significant responses - to assess if normal ##
ggplot(bp.data, aes(x=Condition, y=ug.protein, fill=Time)) +
  geom_boxplot() +
  facet_wrap(~Accession, scales="free_y")

# -----------------------------------------------------------------------------
####### Post-hoc mFuzz of the interaction effects #####
# -----------------------------------------------------------------------------
# Read in the statistical output file and filtering interactions, reshape df ready for mFUZZ clustering

library(tidyr)
library(plyr)
library(reshape2)
library(Mfuzz)


abd<- read.csv("2_Way_ANOVA_Stats_Output.csv")
head(abd)

abd$Time <- revalue(abd$Time, c("D00"="0",
                                "D10"="1",
                                "D20"="2",
                                "D30"="3"))
head(abd)
sig.int <- subset(abd, BH.Interaction < 0.05)


# reshape by the log2FC (Exp/ Con data) for best assessment of time-specific changes
sig <- dcast(sig.int, Accession ~ Time, value.var = "log.FC")
head(sig)

#reorder columns
matrix <- sig
head(matrix)

# Add accessions as rownames for later indexing
rownames(matrix) <- matrix$Accession
matrix$Accession <- NULL
matrix <- data.matrix(matrix)
head(matrix)

# Add a mnaual time point column for mFUZZ
timepoint <- c(0,1,2,3)

# bind that to the dataframe
test_data <- rbind(timepoint, matrix)
row.names(test_data)[1]<-"time"
head(test_data)


#save it to a temp file so ti doesnt clutter up directory
tmp <- tempfile()
write.table(test_data,file=tmp, sep='\t', quote = F, col.names=NA)


#read it back in as an expression set
data <- table2eset(file=tmp)

head(data)

#standardise the data
data.s <- standardise(data)
head(data.s)

m.1 <- mestimate(data.s)
m.1 #estimate of fuzziness
#3.048

#plot data to find optimum number of clusters
Dmin(data.s, m=m.1, crange = seq(2,10,1), repeats = 3, visu = TRUE)
# inspection of the elbow plot indicates 2 clusters is likely sufficient for most efficient separation of adaptive trajectories

clust = 2
c <- mfuzz(data.s, c=clust, m=m.1)

# Adjust the y-axis margin
par(mar = c(2, y_margin, 3, 2))

# Generate the plot with modified y-axis size
mfuzz.plot2(data.s, cl = c, mfrow = c(3, 2), min.mem = 0.5, time.labels = c("D00","D10", "D20", "D30"),  colo = "fancy", x11 = F)
mfuzzColorBar(col='fancy',horizontal = T, main="Membership value")
dev.off()


#Biological analyiss of the clusters #
# Need to use the clustering output to align with intitial input df to extract members of different clusters and their raw values
acore <- acore(data.s,c,min.acore = 0)
str(acore)
cluster.ids <- do.call(rbind, lapply(seq_along(acore), function(i){data.frame(cluster=i, acore[[i]])}))
head(cluster.ids)

#Rename and convert to factors
colnames(cluster.ids) <- c("cluster", "Accession", "membership.val")
cluster.ids$cluster <- as.factor(cluster.ids$cluster)
cluster.ids$Accession <- as.factor(cluster.ids$Accession)

str(cluster.ids)


cluster.plot <- merge(sig, cluster.ids, by="Accession")
head(cluster.plot)
write.csv(cluster.plot, "MFuzz sig interactions ug.protein BH<0.05 2 clusters.csv", row.names = F)



#-----------------------------------------------------------------------------------------
# Analysis of each cluster 

# Conduct string network analysis for inferences on protein-protein interactions within each trajecotry of response
# combine with gene set enrichment analysis to allow for empirical investiation of protein groups within each response. 
#-----------------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(forcats)

cluster.plot<- read.csv("MFuzz sig interactions ug.protein BH<0.05 2 clusters.csv")
### Now need to subset the data based on the cluster groupings
head(cluster.plot)
colnames(cluster.plot) <- c("Accession", "Ctrl","Day 10", "Day 20", "Day 30", "cluster", "Membership value")
rownames(cluster.plot) <- cluster.plot$Accession
cluster.plot$Accession <- NULL
cluster.plot$`Membership value` <- NULL


#-----------------------------------------------------------------------------------------
# String Network Analysis of each cluster 
#-----------------------------------------------------------------------------------------
# make local STRING database for Mouse
dir.create("STRINGdb")
string_db <- STRINGdb$new(version="12", species=10116, # select rat database
                          score_threshold=400, input_directory="STRINGdb")

#-----------------------------------------------------------------------------------------------
# map data to STRING gene identifiers, ensuring there is a clean background data sheet
abd<- read.csv("2_Way_ANOVA_Stats_Output.csv")
abd <- abd[!duplicated(abd$Accession), ]
background <- string_db$map(abd, "Accession", removeUnmappedRows=TRUE) # provide the dataframe name and the name of the column that contains existing identifiers

# Warning:  we couldn't map to STRING 3 % of your identifiers
head(background)

# Set experimental background for appropriate adjustement for dataset/ identification bias
string_db$set_background(background$STRING_id)
string_db$backgroundV

#-------------------------------------------------------------
# conduct STRING GO enrichment analysis 
#-------------------------------------------------------------
c1 <- subset(cluster.plot, cluster == "1")
c2 <- subset(cluster.plot, cluster == "2")

#-------------------------------------------------------------
#### Cluster 1 ####
#-------------------------------------------------------------
head(c1)
c1$names <- rownames(c1)
length(unique(c1$names)) # 74
x <- subset(background, Accession %in% c1$names)
length(unique(x$Accession)) # 67 mapped to background
x.enrichment <- string_db$get_enrichment(x$STRING_id)
head(x.enrichment) # too big to show here
write.csv(x.enrichment, file="STRING GO Enrichment of Cluster 1.csv", row.names=F)

x.enrichment <- read.csv("STRING GO Enrichment of Cluster 1.csv")

bubble.plot <- subset(x.enrichment, category == 'Component' | category == 'Process'| category == 'Function'|category == 'KEGG')

# Optional: filter or preprocess data
bubble.plot.filtered <- bubble.plot %>%
  filter(fdr < 0.05) %>%
  mutate(
    log_fdr = -log10(fdr),
    description = fct_reorder(description, log_fdr)
  )

ggplot(bubble.plot.filtered, aes(x = log_fdr, y = description)) +
  
  # Lollipop stems
  geom_segment(aes(x = 0, xend = log_fdr, yend = description, colour = category),
               linewidth = 0.5) +
  
  # Lollipop heads
  geom_point(aes(size = number_of_genes, colour = category), alpha = 0.85) +
  
  # Reference FDR threshold
  geom_vline(xintercept = -log10(0.05), linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  
  # Size scale and legend tuning
  scale_size_continuous(
    range = c(0.5, 4),
    breaks = c(15, 30, 45),
    limits = c(0, 50),
    name = "Protein Count"
  ) +
  
  # Colour scale
  scale_colour_brewer(palette = "Set2") +
  
  # Axis and legend labels
  labs(
    x = expression("Enrichment (–log"[10]*" FDR)"),
    y = NULL,
    colour = "Category"
  ) +
  
  # AJP-styled theme
  theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(
    axis.title.x = element_text(size = 8, colour = "black", margin = margin(t = 4)),
    axis.text.x = element_text(size = 6, colour = "black"),
    axis.text.y = element_text(size = 6, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.4),
    axis.ticks = element_line(colour = "black", linewidth = 0.4),
    legend.title = element_text(size = 8, colour = "black"),
    legend.text = element_text(size = 6, colour = "black"),
    legend.key.size = unit(0.5, "lines"),
    panel.grid = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(
    size = guide_legend(override.aes = list(size = c(2.5, 3.5, 4.5)))
  )

ggsave(
  filename = "fig_enrichment_lollipop_C1.pdf",
  width = 110,           # mm: single-column panel
  height = 60,          # mm: short layout
  units = "mm",
  device = cairo_pdf,   # ensures vector text/font rendering
  dpi = 600)

getwd()

#-----------------------------------------------------------------------------------------------
# Visualise STRING network
#-----------------------------------------------------------------------------------------------
pdf("Protein network enriched Cluster 1.pdf")
string_db$plot_network(unique(x$STRING_id))
dev.off()

#-------------------------------------------------------------
#### Cluster 2 ####
#-------------------------------------------------------------
head(c2)
c2$names <- rownames(c2)
length(unique(c2$names)) # 113
x <- subset(background, Accession %in% c2$names)

string_db$set_background(background$STRING_id)
string_db$backgroundV

length(unique(x$Accession)) # 110 mapped to background
x.enrichment <- string_db$get_enrichment(x$STRING_id)
head(x.enrichment) # too big to show here
write.csv(x.enrichment, file="STRING GO Enrichment of Cluster 2.csv", row.names=F)

x.enrichment <- read.csv("STRING GO Enrichment of Cluster 2.csv")

bubble.plot <- subset(x.enrichment, category == 'Component' | category == 'Process'| category == 'Function'|category == 'KEGG')

# Optional: filter or preprocess data
bubble.plot.filtered <- bubble.plot %>%
  filter(fdr < 0.0001) %>%
  mutate(
    log_fdr = -log10(fdr),
    description = fct_reorder(description, log_fdr)
  )

ggplot(bubble.plot.filtered, aes(x = log_fdr, y = description)) +
  
  # Lollipop stems
  geom_segment(aes(x = 0, xend = log_fdr, yend = description, colour = category),
               linewidth = 0.5) +
  
  # Lollipop heads
  geom_point(aes(size = number_of_genes, colour = category), alpha = 0.85) +
  
  # Reference FDR threshold
  geom_vline(xintercept = -log10(0.05), linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  
  # Size scale and legend tuning
  scale_size_continuous(
    range = c(1, 4),
    breaks = c(20, 50, 80),
    limits = c(20, 80),
    name = "Protein Count"
  ) +
  
  # Colour scale
  scale_colour_brewer(palette = "Set2") +
  
  # Axis and legend labels
  labs(
    x = expression("Enrichment (–log"[10]*" FDR)"),
    y = NULL,
    colour = "Category"
  ) +
  
  # AJP-styled theme
  theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(
    axis.title.x = element_text(size = 8, colour = "black", margin = margin(t = 4)),
    axis.text.x = element_text(size = 6, colour = "black"),
    axis.text.y = element_text(size = 6, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.4),
    axis.ticks = element_line(colour = "black", linewidth = 0.4),
    legend.title = element_text(size = 8, colour = "black"),
    legend.text = element_text(size = 6, colour = "black"),
    legend.key.size = unit(0.5, "lines"),
    panel.grid = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(
    size = guide_legend(override.aes = list(size = c(2.5, 3.5, 4.5)))
  )

ggsave(
  filename = "fig_enrichment_lollipop_C2.pdf",
  width = 110,           # mm: single-column panel
  height = 60,          # mm: short layout
  units = "mm",
  device = cairo_pdf,   # ensures vector text/font rendering
  dpi = 600)


#-----------------------------------------------------------------------------------------------
# Visualise STRING network
#-----------------------------------------------------------------------------------------------
pdf("Protein network enriched Cluster 2.pdf")
string_db$plot_network(unique(x$STRING_id))
dev.off()


