# -------------------------------------------------------------------------
# Script: 03_differential_analysis_example.R
#
# Purpose:
# Demonstration of differential protein-level analysis and functional
# interpretation using quantitative proteomics data.
#
# Overview:
# This script illustrates a typical downstream workflow following
# quality control and exploratory analysis, including:
#   - Differential abundance testing between conditions
#   - Multiple-testing correction
#   - Integration of protein abundance and synthesis-rate outputs
#   - Functional interpretation via Gene Ontology enrichment
#
# -------------------------------------------------------------------------

# --- Load libraries ----------------------------------------------------------
library(ggrepel)
library(qvalue)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(patchwork)
library(readxl)
library(stringr)
library(tidyverse)

# Explicitly set tidyverse functions to avoid conflicts
select   <- dplyr::select
filter   <- dplyr::filter
mutate   <- dplyr::mutate
summarise <- dplyr::summarise
arrange  <- dplyr::arrange
rename   <- dplyr::rename

# --- Utility functions -------------------------------------------------------

# Simplify enrichment results (for speed)
simplify_top <- function(ego, top_n = 100, cutoff = 0.7) {
  ego_df <- as.data.frame(ego) %>% arrange(p.adjust) %>% slice_head(n = top_n)
  ego_sub <- ego[ego_df$ID]
  simplify(ego_sub, cutoff = cutoff, by = "p.adjust", select_fun = min, measure = "Wang")
}

# Save helper
save_fig <- function(plot, name, n_rows = 1, width = 180) {
  height <- 70 * n_rows
  ggsave(
    filename = paste0(name, ".pdf"),
    plot = plot, width = width, height = height, units = "mm",
    dpi = 600, bg = "white"
  )
  ggsave(
    filename = paste0(name, ".svg"),
    plot = plot, width = width, height = height, units = "mm",
    dpi = 600, bg = "white"
  )
}


# Global default theme for all figure panels
theme_set(
  theme_classic(base_size = 8) +
    theme(
      aspect.ratio = 1,
      axis.text = element_text(size = 6, colour = "black"),
      axis.title = element_text(size = 8, face = "plain"),
      axis.text.y = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 7),
      legend.text  = element_text(size = 6),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )
)

# Common coordinate systems
coord_volcano   <- coord_cartesian(xlim = c(-5, 5), ylim = c(0, 5), clip = "off")
coord_enrichment <- coord_cartesian(xlim = c(0, 0.5), clip = "off")
# --- Part 1: FSR data processing ---------------------------------------------
fsr <- readRDS("PEA Cell Exp - Proteo_ADPT data.RDS")

head(fsr)

fsr$fsr <- fsr$Kdeg*100
fsr$log_fsr <- log10(fsr$fsr)

# Load and subset
I <- subset(fsr, condition != "P")

# 1. Statistics
# --- 1A. Compute per-accession ANOVA p-values --------------------------------
p.val <- I %>%
  group_by(accession) %>%
  summarise(
    p.val = summary(aov(fsr ~ condition, data = cur_data()))[[1]][["Pr(>F)"]][1],
    .groups = "drop"
  )

# --- 1B. Calculate mean and SD for each accession × condition ----------------
output <- I %>%
  group_by(accession, description, condition) %>%
  summarise(
    mean = mean(fsr, na.rm = TRUE),
    SD   = sd(fsr, na.rm = TRUE),
    .groups = "drop"
  )

# --- 1C. Pivot to wide format ------------------------------------------------
output_wide <- output %>%
  pivot_wider(
    id_cols = c(accession, description),
    names_from  = condition,
    values_from = c(mean, SD),
    names_glue  = "{condition}_{.value}"
  )

# --- 1D. Merge and calculate summary statistics ------------------------------
stats_output <- output_wide %>%
  left_join(p.val, by = "accession") %>%
  mutate(
    q.val     = qvalue(p.val)$qvalues,
    p.adj     = qvalue(p.val, pi0 = 1)$qvalues,
    IvsVC.FC  = log2(I_mean / VC_mean),
    direction = case_when(
      p.val > 0.05 ~ "NS",
      IvsVC.FC > 0 ~ "Up",
      IvsVC.FC < 0 ~ "Down",
      TRUE         ~ "Neutral"
    ),
    negLog10P = -log10(p.val)
  )

length(unique(stats_output$accession)) #1542 protein FSRs.
I_fsr_sig <- subset(stats_output, p.val <= 0.05) # 172 proteins sig P < 0.05
summary(I_fsr_sig$q.val)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1040  0.1332  0.1472  0.1476  0.1620  0.1745 

head(stats_output)

write.csv(stats_output, "FSR_Ibu_vs_VC_stats_output.csv", row.names = F)

I_fsr_sig %>%
  filter(direction == "Up") %>%
  summarise(n_up = n_distinct(accession)) #165 increased in FSR vs VC

I_fsr_sig %>%
  filter(direction == "Down") %>%
  summarise(n_down = n_distinct(accession)) #7 decreased in FSR vs VC

# Largely increases in protien synthesis; anabolic response?

# --- Part 2: Volcano plot (FSR) ---------------------------------------------
ibu_fsr_volcano <- ggplot(stats_output, aes(x = IvsVC.FC, y = negLog10P)) +
  geom_point(aes(color = direction), alpha = 0.4, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = c(Down = "#4575b4", Up = "firebrick", NS = "grey70")) +
  labs(
    x = expression(log[2]~"Fold Change (IBU / VC)"),
    y = expression(-log[10]~"P-value")
  ) +
  theme_classic(base_size = 8) +
  theme(legend.position = "none")+
  coord_cartesian(clip = "off") +          # ensures equal margins
  theme(
    aspect.ratio = 1,                      # square panel ratio
    plot.margin = margin(2, 2, 2, 2),      # uniform padding
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 6),
    legend.position = "none",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 6))+
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 5))

ibu_fsr_volcano

# --- Part 3: GO Enrichment (FSR Upregulated) --------------------------------
# Map UniProt → gene symbols
id_map <- read_excel("idmapping_2025_10_14.xlsx")

# Flexible rename that works even if columns vary slightly
id_map <- id_map %>%
  rename_with(~ "accession", matches("Entry.?Name", ignore.case = TRUE)) %>%
  rename_with(~ "gene.name", matches("Gene.?Names", ignore.case = TRUE)) %>%
  mutate(
    gene.name = str_extract(gene.name, "^[^ ]+")  # take first symbol if multiple
  )

head(id_map)
head(fsr)

fsr <- fsr %>%
  mutate(accession = paste0(accession, "_MOUSE"))

bg <- fsr

head(bg)

bg <- bg[,c(1,2,3)]

write.csv(bg, "FSR_GO_enrichment_background_list.csv", row.names = F)


up_ibu <- subset(stats_output, direction == "Up") %>%
  mutate(accession = paste0(accession, "_MOUSE")) %>%
  left_join(select(id_map, accession, gene.name), by = "accession")

sym_univ_chr <- unique(trimws(as.character(bg$gene.name)))
sym_hits_chr <- unique(trimws(as.character(up_ibu$gene.name)))

ent_univ <- mapIds(org.Mm.eg.db, keys = sym_univ_chr, column = "ENTREZID",
                   keytype = "SYMBOL", multiVals = "first")
ent_hits <- mapIds(org.Mm.eg.db, keys = sym_hits_chr, column = "ENTREZID",
                   keytype = "SYMBOL", multiVals = "first")

universe_entrez <- as.character(na.omit(ent_univ))
hits_entrez <- as.character(na.omit(ent_hits))

ego_ibu_up <- enrichGO(
  gene = hits_entrez,
  universe = universe_entrez,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.99,
  qvalueCutoff = 0.99,
  readable = TRUE
)

head(ego_ibu_up)

# --- Part 4: Dot plot (GO enrichment) ---------------------------------------
plot_df <- as.data.frame(ego_ibu_up) %>%
  arrange(pvalue) %>%
  slice_head(n = 100) %>%
  mutate(
    Description = str_squish(Description),
    Description = str_replace_all(Description, regex("regulation", ignore_case = TRUE), "reg."),
    Description = case_when(
      str_detect(Description, fixed("adaptive immune response based on somatic recombination")) ~
        "adaptive immune reg. (IgSF)",
      TRUE ~ Description
    )
  ) %>%
  arrange(RichFactor) %>%                          # order by RichFactor
  mutate(Description = factor(Description, levels = Description))  # lock factor order

head(plot_df)

write.csv(plot_df, "FSR_Ibu_GO_enrichment_results_list.csv", row.names = F)


# Take top 10 for plotting #
plot_df <- as.data.frame(ego_ibu_up) %>%
  arrange(pvalue) %>%
  slice_head(n = 10) %>%
  mutate(
    Description = str_squish(Description),
    Description = str_replace_all(Description, regex("regulation", ignore_case = TRUE), "reg."),
    Description = case_when(
      str_detect(Description, fixed("adaptive immune response based on somatic recombination")) ~
        "adaptive immune reg. (IgSF)",
      TRUE ~ Description
    )
  ) %>%
  arrange(RichFactor) %>%                          # order by RichFactor
  mutate(Description = factor(Description, levels = Description))  # lock factor order

ibu_fsr_up_dot <- ggplot(plot_df, aes(x = RichFactor, y = Description)) + 
  geom_segment(aes(x = 0, xend = RichFactor, y = Description, yend = Description),
               colour = "grey70", linewidth = 0.5) +
  geom_point(aes(size = Count, fill = p.adjust),
             shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_steps2(
    low = "firebrick", mid = 'white', high = "#eeaf61",
    midpoint = 0.5, limits = c(0, 1),
    n.breaks = 11, name = "adjusted p-value",
    guide = guide_colorsteps(barwidth = 0.6, barheight = 5, ticks = TRUE)
  ) +
  scale_size(range = c(1, 4),
             breaks = seq(floor(min(plot_df$Count)), ceiling(max(plot_df$Count)), by = 10),
             name = "Count") +
  labs(x = "Enrichment factor", y = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.y = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  ) +coord_cartesian(clip = "off") +          # ensures equal margins
  theme(
    aspect.ratio = 1,                      # square panel ratio
    plot.margin = margin(5, 5, 5, 5),      # uniform padding
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 6),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 6))+
  coord_cartesian(xlim = c(0, 0.5))

ibu_fsr_up_dot


# --- Part 5: Combine (FSR volcano + dot) ------------------------------------
fig_ibu_fsr <- ibu_fsr_volcano + ibu_fsr_up_dot +
  plot_annotation(tag_levels = 'A', tag_suffix = ".") &
  theme(plot.margin = margin(2, 2, 2, 2, "mm"))

fig_ibu_fsr

save_fig(fig_ibu_fsr, "fig_ibu_fsr", n_rows = 1)



# --- Part 6: ABD data analysis (same workflow condensed) ---------------------
# --- Load & clean data ----------------------------------------------------
library(tidyverse)
library(qvalue)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(readxl)
library(stringr)
library(patchwork)

abd_wide[abd_wide == 0] <- NA
abd_wide <- abd_wide[rowSums(is.na(abd_wide[, -c(1:2)])) == 0, ]

head(abd_wide)

# --- 2. Reshape & tidy -------------------------------------------------------

abd_df <- abd_wide %>%
  pivot_longer(cols = 3:11, names_to = "sample", values_to = "abd") %>%
  separate(sample, into = c("condition", "time", "replicate"), sep = "_") %>%
  mutate(accession = str_remove(accession, "_MOUSE"))

# Filter relevant conditions
abd_I <- abd_df %>% filter(condition != "P")

# --- 3. Stats: one-way ANOVA per accession ----------------------------------

abd_pval <- abd_I %>%
  group_by(accession) %>%
  summarise(
    p.val = summary(aov(abd ~ condition, data = cur_data()))[[1]][["Pr(>F)"]][1],
    .groups = "drop"
  )

abd_summary <- abd_I %>%
  group_by(accession, description, condition) %>%
  summarise(mean = mean(abd), SD = sd(abd), .groups = "drop") %>%
  pivot_wider(
    id_cols = c(accession, description),
    names_from = condition,
    values_from = c(mean, SD),
    names_glue = "{condition}_{.value}"
  )

abd_stats <- abd_summary %>%
  left_join(abd_pval, by = "accession") %>%
  mutate(
    q.val = qvalue(p.val)$qvalues,
    p.adj = qvalue(p.val, pi0 = 1)$qvalues,
    IvsVC.FC = log2(I_mean / VC_mean),
    direction = case_when(
      p.val > 0.05 ~ "NS",
      IvsVC.FC > 0 ~ "Up",
      IvsVC.FC < 0 ~ "Down",
      TRUE ~ "Neutral"
    ),
    negLog10P = -log10(p.val)
  )

head(abd_stats)

write.csv(abd_stats, "ABD_Ibu_vs_VC_stats_output.csv", row.names = F)

length(unique(abd_stats$accession)) #3085 protein abundances.
I_abd_sig <- subset(abd_stats, p.val <= 0.05) # 561 proteins sig P < 0.05
summary(I_abd_sig$q.val)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.03464 0.08555 0.10413 0.10778 0.13059 0.15754 

I_abd_sig %>%
  filter(direction == "Up") %>%
  summarise(n_up = n_distinct(accession)) #352 increased vs VC

I_abd_sig %>%
  filter(direction == "Down") %>%
  summarise(n_down = n_distinct(accession)) #209 decreased vs VC

# --- 4. Volcano Plot ---------------------------------------------------------

ibu_abd_volcano <- ggplot(abd_stats, aes(x = IvsVC.FC, y = negLog10P)) +
  geom_point(aes(color = direction), alpha = 0.4, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(
    values = c(Down = "#4575b4", Up = "firebrick", NS = "grey70")
  ) +
  labs(
    x = expression(log[2]~"Fold Change (IBU / VC)"),
    y = expression(-log[10]~"P-value")
  ) +
  theme_classic(base_size = 8) +
  theme(legend.position = "none")+
  coord_cartesian(clip = "off") +          # ensures equal margins
  theme(
    aspect.ratio = 1,                      # square panel ratio
    plot.margin = margin(2, 2, 2, 2),      # uniform padding
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 6),
    legend.position = "none",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 6))+
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 5))

ibu_abd_volcano


# --- 5. Gene name mapping ----------------------------------------------------
id_map <- read_excel("idmapping_2025_10_14.xlsx") %>%
  rename_with(~ "accession", matches("Entry.?Name", ignore.case = TRUE)) %>%
  rename_with(~ "gene.name", matches("Gene.?Names", ignore.case = TRUE)) %>%
  mutate(gene.name = str_extract(gene.name, "^[^ ]+"))

abd_df <- abd_df %>%
  mutate(accession = paste0(accession, "_MOUSE")) %>%
  left_join(dplyr::select(id_map, accession, gene.name), by = "accession") %>%
  mutate(gene.name = ifelse(is.na(gene.name), accession, gene.name))

bg_abd <- abd_df %>% select(accession, gene.name)

head(bg_abd)

write.csv(bg_abd, "ABD_GO_enrichment_background_list.csv", row.names = F)


# --- 6. GO enrichment (Up & Down) -------------------------------------------
run_enrichGO_mouse <- function(hit_genes, bg_genes) {
  sym_univ_chr <- unique(trimws(as.character(bg_genes)))
  sym_hits_chr <- unique(trimws(as.character(hit_genes)))
  
  ent_univ <- mapIds(org.Mm.eg.db, keys = sym_univ_chr,
                     column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  ent_hits <- mapIds(org.Mm.eg.db, keys = sym_hits_chr,
                     column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  
  universe_entrez <- as.character(na.omit(ent_univ))
  hits_entrez     <- as.character(na.omit(ent_hits))
  
  enrichGO(
    gene          = hits_entrez,
    universe      = universe_entrez,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.99,
    qvalueCutoff  = 0.99,
    readable      = TRUE
  )
}

# --- 7. Run enrichment for Up & Down ----------------------------------------
up_hits   <- abd_stats %>% filter(direction == "Up") %>% pull(accession)
down_hits <- abd_stats %>% filter(direction == "Down") %>% pull(accession)

up_genes   <- abd_df %>% filter(accession %in% paste0(up_hits, "_MOUSE")) %>% pull(gene.name)
down_genes <- abd_df %>% filter(accession %in% paste0(down_hits, "_MOUSE")) %>% pull(gene.name)
bg_genes   <- bg_abd$gene.name

ego_abd_up   <- run_enrichGO_mouse(up_genes, bg_genes)
ego_abd_down <- run_enrichGO_mouse(down_genes, bg_genes)

up_plot_df <- as.data.frame(ego_abd_up) %>%
  arrange(pvalue) %>%
  slice_head(n = 100) %>%
  mutate(
    Description = str_squish(Description),
    Description = str_replace_all(Description, regex("regulation", ignore_case = TRUE), "reg."),
    Description = case_when(
      str_detect(Description, fixed("adaptive immune response based on somatic recombination")) ~
        "adaptive immune reg. (IgSF)",
      TRUE ~ Description
    )
  ) %>%
  arrange(RichFactor) %>%                          # order by RichFactor
  mutate(Description = factor(Description, levels = Description))  # lock factor order

head(up_plot_df)

write.csv(up_plot_df, "ABD_Ibu_Upreg_GO_enrichment_results_list.csv", row.names = F)


down_plot_df <- as.data.frame(ego_abd_down) %>%
  arrange(pvalue) %>%
  slice_head(n = 100) %>%
  mutate(
    Description = str_squish(Description),
    Description = str_replace_all(Description, regex("regulation", ignore_case = TRUE), "reg."),
    Description = case_when(
      str_detect(Description, fixed("nucleobase-containing small molecule metabolic process")) ~
        "nucleobase-containing-metabolic process",
      TRUE ~ Description
    )
  ) %>%
  arrange(RichFactor) %>%                          # order by RichFactor
  mutate(Description = factor(Description, levels = Description))  # lock factor order

head(down_plot_df)

write.csv(down_plot_df, "ABD_Ibu_Downreg_GO_enrichment_results_list.csv", row.names = F)




# --- 8. Plot enrichment results (top 10 only) ---------------------------------------------
up_plot_df <- as.data.frame(ego_abd_up) %>%
  arrange(pvalue) %>%
  slice_head(n = 100) %>%
  mutate(
    Description = str_squish(Description),
    Description = str_replace_all(Description, regex("regulation", ignore_case = TRUE), "reg."),
    Description = case_when(
      str_detect(Description, fixed("adaptive immune response based on somatic recombination")) ~
        "adaptive immune reg. (IgSF)",
      TRUE ~ Description
    )
  ) %>%
  arrange(RichFactor) %>%                          # order by RichFactor
  mutate(Description = factor(Description, levels = Description))  # lock factor order

down_plot_df <- as.data.frame(ego_abd_down) %>%
  arrange(pvalue) %>%
  slice_head(n = 10) %>%
  mutate(
    Description = str_squish(Description),
    Description = str_replace_all(Description, regex("regulation", ignore_case = TRUE), "reg."),
    Description = case_when(
      str_detect(Description, fixed("nucleobase-containing small molecule metabolic process")) ~
        "nucleobase-containing-metabolic process",
      TRUE ~ Description
    )
  ) %>%
  arrange(RichFactor) %>%                          # order by RichFactor
  mutate(Description = factor(Description, levels = Description))  # lock factor order
plot_go_dot_up <- ggplot(up_plot_df, aes(x = RichFactor, y = Description)) + 
  geom_segment(aes(x = 0, xend = RichFactor, y = Description, yend = Description),
               colour = "grey70", linewidth = 0.5) +
  geom_point(aes(size = Count, fill = p.adjust),
             shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_steps2(
    low = "firebrick", mid = 'white', high = "#eeaf61",
    midpoint = 0.5, limits = c(0, 1),
    n.breaks = 11, name = "adjusted p-value",
    guide = guide_colorsteps(barwidth = 0.6, barheight = 5, ticks = TRUE)
  ) +
  scale_size(range = c(1, 4),
             breaks = seq(floor(8), ceiling(18), by = 2),
             name = "Count") +
  labs(x = "Enrichment factor", y = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.y = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  ) +coord_cartesian(clip = "off") +          # ensures equal margins
  theme(
    aspect.ratio = 1,                      # square panel ratio
    plot.margin = margin(5, 5, 5, 5),      # uniform padding
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 6),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 6))+
  coord_cartesian(xlim = c(0, 0.5))

plot_go_dot_up


plot_go_dot_dwn <- ggplot(down_plot_df, aes(x = RichFactor, y = Description)) + 
  geom_segment(aes(x = 0, xend = RichFactor, y = Description, yend = Description),
               colour = "grey70", linewidth = 0.5) +
  geom_point(aes(size = Count, fill = p.adjust),
             shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_steps2(
    low = "#4575b4", mid = 'white', high = "#eeaf61",
    midpoint = 0.5, limits = c(0, 1),
    n.breaks = 11, name = "adjusted p-value",
    guide = guide_colorsteps(barwidth = 0.6, barheight = 5, ticks = TRUE)
  ) +
  scale_size(range = c(1, 4),
             breaks = seq(floor(0), ceiling(max(plot_df$Count)), by = 5),
             name = "Count") +
  labs(x = "Enrichment factor", y = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.y = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  ) +coord_cartesian(clip = "off") +          # ensures equal margins
  theme(
    aspect.ratio = 1,                      # square panel ratio
    plot.margin = margin(5, 5, 5, 5),      # uniform padding
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 6),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 6))+
  coord_cartesian(xlim = c(0, 0.5))

plot_go_dot_dwn

# --- 9. Combine ABD figure ---------------------------------------------------

ibu_abd_fig <- (ibu_abd_volcano + plot_go_dot_up) &
  plot_annotation(tag_levels = "A", tag_suffix = ".") &
  theme(plot.margin = margin(2, 2, 2, 2, "mm")) # Switch between including up and down plot to maintain size ratio consistencies etc


ibu_abd_fig

save_fig(ibu_abd_fig, "fig_ibu_abd", n_rows = 1)

# ===============================================================
# End of Script
# ===============================================================

