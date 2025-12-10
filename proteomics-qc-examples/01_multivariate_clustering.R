# Script: 02_multivariate_exploration_example.R
#
# Purpose:
# Exploratory and supervised multivariate analysis of quantitative
# proteomics data to assess global structure and condition-related patterns.
#
# Overview:
# This script demonstrates a standard analytical workflow applied after
# initial quality control has confirmed acceptable data completeness and
# inter-sample consistency. Analyses include:
#   - Principal Component Analysis (PCA) for unsupervised structure
#   - Partial Least Squares Discriminant Analysis (PLS-DA) for supervised
#     pattern assessment
#   - Cross-validation to evaluate classification performance

head(abd_clean)

meta_cols <- c(
  "accession", "accession_full", "accession_short",
  "gene_name", "description", "peptide_count", "unique_peptides", "mass"
)

sample_cols <- setdiff(names(abd_clean), meta_cols)

pca_mat <- abd_clean %>%
  drop_na(all_of(sample_cols)) %>%
  select(all_of(sample_cols))

pca_input <- t(pca_mat)

pca_res <- prcomp(pca_input, scale. = TRUE)

scores <- as.data.frame(pca_res$x)
scores$sample <- rownames(scores)

sample_meta <- meta   # whatever your metadata table is called

sample_meta <- sample_meta %>%
  mutate(sample = paste0("M_", sample_code)) # paste M_ to match sample coding

scores <- scores %>% left_join(sample_meta, by = "sample")

scores_annot <- scores %>%
  mutate(
    condition_colour = case_when(
      Condition == "IMMOB" ~ "red",
      Condition == "RT" ~ "blue",
      TRUE ~ "grey60"
    ),
    fill_value = case_when(
      Time == "pre" ~ "white",
      Time == "post" ~ condition_colour,   # filled with the condition colour
      TRUE ~ "white"
    )
  )

p_pca <- ggplot(scores_annot, aes(x = PC1, y = PC2)) +
  
  geom_point(
    aes(
      colour = condition_colour,
      fill = fill_value
    ),
    shape = 21,        # circle with border + fill
    size = 4,
    stroke = 1.1,
    alpha = 0.95
  ) +
  
  scale_colour_identity() +
  scale_fill_identity() +
  
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    subtitle = "Red = IMMOB | Blue = RT | Empty = pre | Filled = post",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "%)")
  )

p_pca


plot_ly(
  scores_annot,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  type = "scatter3d",
  mode = "markers",
  marker = list(
    size = 6,
    color = ~condition_colour,
    symbol = ~ifelse(Time == "pre", "circle-open", "circle"),
    line = list(width = 1, color = ~condition_colour)
  ),
  text = ~sample,
  hoverinfo = "text"
) %>%
  layout(
    title = "3D PCA (PC1–PC2–PC3)",
    scene = list(
      xaxis = list(title = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "%)")),
      yaxis = list(title = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "%)")),
      zaxis = list(title = paste0("PC3 (", round(summary(pca_res)$importance[2,3] * 100, 1), "%)"))
    )
  )


# Pretty clear no separation in 2D or 3D PCA plots. - Groups cannot be separated in unsupervised way. 


# Try PLS-DA, supervised clustering analysis #
mat <- abd_clean %>%
  drop_na(all_of(sample_cols)) %>%
  select(all_of(sample_cols))

X <- t(mat)   # predictor matrix (samples × proteins)

Y <- scores_annot$Condition   # factor vector

stopifnot(all(rownames(X) == scores_annot$sample))

plsda_res <- plsda(X, Y, ncomp = 3)

set.seed(123)
perf_res <- perf(plsda_res, validation = "Mfold", folds = 5, progressBar = TRUE)

perf_res$error.rate

# $overall
# max.dist centroids.dist mahalanobis.dist
# comp1 0.4821429      0.4821429        0.4821429
# comp2 0.5178571      0.4285714        0.5000000
# comp3 0.3750000      0.4642857        0.3750000
# 
# $BER
# max.dist centroids.dist mahalanobis.dist
# comp1 0.4821429      0.4821429        0.4821429
# comp2 0.5178571      0.4285714        0.5000000
# comp3 0.3750000      0.4642857        0.3750000

perf_res$BER

perf_res$error.rate.class

scores_pls <- as.data.frame(plsda_res$variates$X)
scores_pls$sample <- rownames(scores_pls)

scores_pls <- scores_pls %>%
  left_join(scores_annot, by = "sample")


p_plsda <- ggplot(scores_pls, aes(x = comp1, y = comp2)) +
  geom_point(
    aes(colour = condition_colour, fill = fill_value),
    size = 4, shape = 21, stroke = 1.1, alpha = 0.9
  ) +
  scale_colour_identity() +
  scale_fill_identity() +
  theme_bw(base_size = 13) +
  labs(
    x = "PLS-DA Component 1",
    y = "PLS-DA Component 2"
  )

p_plsda

vip <- vip(plsda_res)
top20_vip <- sort(vip[,1], decreasing = TRUE)[1:20]
top20_vip

perf_res$features$stable

p_clustering <- p_pca + p_plsda
p_clustering

ggsave(
  filename = file.path(dir_figs, "abd_myofib_clustering_plots.pdf"),
  plot = p_clustering,
  width = 8,
  height = 5,
  dpi = 600
)
# No meaningful seperation across PCA or PLS-DA. Overall, meaning no dramatic global shift in protein abundance data.

