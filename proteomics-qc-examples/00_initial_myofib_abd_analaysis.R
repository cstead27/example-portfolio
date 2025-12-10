
### 01_qc_proteomics_example.R ###

# Purpose:
# Demonstration of quality control and exploratory analysis
# for quantitative proteomics data.
#
# This script illustrates:
# - Missingness assessment
# - Abundance distribution checks
# - Inter-sample correlation structure
#
# Data are simulated / anonymised and provided for demonstration only.


# Load packaged #
library(tidyverse)
library(pheatmap)
library(patchwork)
library(ggrepel)
library(plotly)
library(mixOmics)
library(limma)

# Load project paths #
source(here::here("R/paths.R"))

# Load progenesis output cleaning helper #
source(here::here("R/abundance_data_cleaning.R"))


#### list files in raw data folder ####
raw_files <- list.files(dir_raw, full.names = T)
raw_files


#### Load and inspect meta data ####
meta <- read_csv(raw_files[3])
glimpse(meta)
head(meta)
#   Sample            Participant sample_code Condition Time    Day
#   <chr>             <chr>             <dbl> <chr>     <chr> <dbl>
# 1 SLRI_01_D-2_IMMOB SLRI_01               1 IMMOB     pre      -2
# 2 SLRI_01_D0_-4h_RT SLRI_01               2 RT        pre       0
# 3 SLRI_01_D10_RT    SLRI_01               3 RT        post     10
# 4 SLRI_01_D10_IMMOB SLRI_01               4 IMMOB     post     10
# 5 SLRI_02_D-2_IMMOB SLRI_02               5 IMMOB     pre      -2
# 6 SLRI_02_D0_-4h_RT SLRI_02               6 RT        pre       0

colnames(meta)
# [1] "Sample"      "Participant" "sample_code" "Condition"   "Time"        "Day"  

dim(meta)

unique(names(meta))


#### Load and inspect abd data ####
raw <- read_csv(raw_files[1])
head(raw)

# apply cleaning source code #
abd_clean <- clean_abundance_data(raw_files[1])
head(abd_clean)


# quantification summary #
total_peptides <- sum(abd_clean$peptide_count)
total_peptides 

unique_peptides <- sum(abd_clean$unique_peptides)
unique_peptides

proteins_with_unique <- abd_clean %>%
  filter(unique_peptides >= 1) %>%
  distinct(accession) %>%
  nrow()

proteins_with_unique

# In total we quantify the abundance of 4031 peptides from the myofibrillar fraction of skeletal muscle with 2902 uniques peptides quantified belonging to 290 proteins

abd_clean <- abd_clean %>%
  filter(unique_peptides >= 1)


### proteins with complete experimental data ###
meta_cols <- c(
  "accession", "accession_full",
  "gene_name", "description", "peptide_count", "unique_peptides", "mass"
)

sample_cols <- setdiff(names(abd_clean), meta_cols)

abd_clean <- abd_clean %>%
  mutate(across(all_of(sample_cols),
                ~ na_if(., 0)))   # convert zeros → NA

n_complete_proteins <- abd_clean %>%
  drop_na(all_of(sample_cols)) %>%   # keep only rows with NO missing in sample cols
  nrow()

n_complete_proteins #273 myofibrillar proteins with complete data. 

## Data missingness profile ##
na_summary <- abd_clean %>%
  summarise(across(all_of(sample_cols), ~ sum(is.na(.)))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "sample",
    values_to = "na_count"
  )

ggplot(na_summary, aes(x = sample, y = na_count)) +
  geom_col(fill = "grey60") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Missing Data per Sample (Protein-level)",
    x = "Sample",
    y = "Number of NA values"
  )

na_prop <- abd_clean %>%
  summarise(across(all_of(sample_cols), ~ mean(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "sample", values_to = "na_prop")

p_bar <-ggplot(na_prop, aes(sample, na_prop)) +
  geom_col(fill = "grey50") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    y = "Missing proportion"
  )

p_bar

# Missingness is relatively balanced, no sample dominating the missingness. 


complete_mat <- abd_clean %>%
  drop_na(all_of(sample_cols)) %>%
  select(all_of(sample_cols))

cor_mat <- cor(complete_mat, use = "pairwise.complete.obs", method = "pearson")

p_heat <- pheatmap(
  cor_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = FALSE,
  fontsize = 6)

summary(cor_mat)

# M_1              M_2              M_3              M_4              M_5        
# Min.   :0.9114   Min.   :0.8771   Min.   :0.8810   Min.   :0.8841   Min.   :0.8984  
# 1st Qu.:0.9863   1st Qu.:0.9867   1st Qu.:0.9877   1st Qu.:0.9850   1st Qu.:0.9850  
# Median :0.9905   Median :0.9911   Median :0.9926   Median :0.9889   Median :0.9876  
# Mean   :0.9876   Mean   :0.9878   Mean   :0.9887   Mean   :0.9853   Mean   :0.9857  
# 3rd Qu.:0.9940   3rd Qu.:0.9937   3rd Qu.:0.9952   3rd Qu.:0.9926   3rd Qu.:0.9914  
# Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
# M_6              M_7              M_8              M_9              M_10       
# Min.   :0.8869   Min.   :0.8208   Min.   :0.8943   Min.   :0.8208   Min.   :0.8823  
# 1st Qu.:0.9856   1st Qu.:0.9403   1st Qu.:0.9885   1st Qu.:0.8850   1st Qu.:0.9876  
# Median :0.9887   Median :0.9540   Median :0.9926   Median :0.8962   Median :0.9912  
# Mean   :0.9863   Mean   :0.9517   Mean   :0.9893   Mean   :0.8954   Mean   :0.9879  
# 3rd Qu.:0.9932   3rd Qu.:0.9678   3rd Qu.:0.9950   3rd Qu.:0.9040   3rd Qu.:0.9949  
# Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
# M_11             M_12             M_13             M_14             M_15       
# Min.   :0.8836   Min.   :0.8881   Min.   :0.8830   Min.   :0.8784   Min.   :0.8808  
# 1st Qu.:0.9881   1st Qu.:0.9837   1st Qu.:0.9856   1st Qu.:0.9878   1st Qu.:0.9900  
# Median :0.9923   Median :0.9899   Median :0.9903   Median :0.9917   Median :0.9932  
# Mean   :0.9890   Mean   :0.9865   Mean   :0.9876   Mean   :0.9883   Mean   :0.9896  
# 3rd Qu.:0.9955   3rd Qu.:0.9935   3rd Qu.:0.9942   3rd Qu.:0.9946   3rd Qu.:0.9955  
# Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
# M_16             M_17             M_18             M_19             M_20       
# Min.   :0.8946   Min.   :0.8987   Min.   :0.9097   Min.   :0.9012   Min.   :0.8996  
# 1st Qu.:0.9923   1st Qu.:0.9916   1st Qu.:0.9821   1st Qu.:0.9898   1st Qu.:0.9908  
# Median :0.9950   Median :0.9944   Median :0.9889   Median :0.9930   Median :0.9935  
# Mean   :0.9914   Mean   :0.9908   Mean   :0.9845   Mean   :0.9893   Mean   :0.9908  
# 3rd Qu.:0.9967   3rd Qu.:0.9959   3rd Qu.:0.9937   3rd Qu.:0.9945   3rd Qu.:0.9959  
# Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
# M_21             M_22             M_23             M_24             M_25       
# Min.   :0.9096   Min.   :0.8936   Min.   :0.8837   Min.   :0.8994   Min.   :0.9063  
# 1st Qu.:0.9837   1st Qu.:0.9798   1st Qu.:0.9877   1st Qu.:0.9910   1st Qu.:0.9840  
# Median :0.9908   Median :0.9856   Median :0.9917   Median :0.9942   Median :0.9883  
# Mean   :0.9859   Mean   :0.9834   Mean   :0.9888   Mean   :0.9910   Mean   :0.9851  
# 3rd Qu.:0.9937   3rd Qu.:0.9899   3rd Qu.:0.9951   3rd Qu.:0.9959   3rd Qu.:0.9921  
# Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
# M_26             M_27             M_28             M_29             M_30       
# Min.   :0.8944   Min.   :0.8861   Min.   :0.8982   Min.   :0.8962   Min.   :0.9011  
# 1st Qu.:0.9849   1st Qu.:0.9874   1st Qu.:0.9806   1st Qu.:0.9806   1st Qu.:0.9897  
# Median :0.9911   Median :0.9915   Median :0.9878   Median :0.9868   Median :0.9929  
# Mean   :0.9863   Mean   :0.9877   Mean   :0.9829   Mean   :0.9834   Mean   :0.9899  
# 3rd Qu.:0.9946   3rd Qu.:0.9947   3rd Qu.:0.9929   3rd Qu.:0.9903   3rd Qu.:0.9956  
# Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
# M_31             M_32             M_33             M_34             M_35       
# Min.   :0.8900   Min.   :0.8852   Min.   :0.9057   Min.   :0.8921   Min.   :0.8823  
# 1st Qu.:0.9865   1st Qu.:0.9888   1st Qu.:0.9879   1st Qu.:0.9884   1st Qu.:0.9834  
# Median :0.9904   Median :0.9923   Median :0.9902   Median :0.9926   Median :0.9884  
# Mean   :0.9877   Mean   :0.9891   Mean   :0.9883   Mean   :0.9891   Mean   :0.9856  
# 3rd Qu.:0.9944   3rd Qu.:0.9957   3rd Qu.:0.9943   3rd Qu.:0.9949   3rd Qu.:0.9933  
# Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
# M_36             M_37             M_38             M_39             M_40       
# Min.   :0.8614   Min.   :0.9046   Min.   :0.9065   Min.   :0.9017   Min.   :0.8963  
# 1st Qu.:0.9760   1st Qu.:0.9925   1st Qu.:0.9890   1st Qu.:0.9732   1st Qu.:0.9890  
# Median :0.9827   Median :0.9946   Median :0.9934   Median :0.9776   Median :0.9931  
# Mean   :0.9795   Mean   :0.9917   Mean   :0.9891   Mean   :0.9757   Mean   :0.9888  
# 3rd Qu.:0.9894   3rd Qu.:0.9965   3rd Qu.:0.9960   3rd Qu.:0.9815   3rd Qu.:0.9953  
# Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
# M_41             M_42             M_43             M_44             M_45       
# Min.   :0.8991   Min.   :0.9038   Min.   :0.8695   Min.   :0.8797   Min.   :0.8933  
# 1st Qu.:0.9897   1st Qu.:0.9802   1st Qu.:0.9611   1st Qu.:0.9819   1st Qu.:0.9836  
# Median :0.9933   Median :0.9885   Median :0.9669   Median :0.9875   Median :0.9875  
# Mean   :0.9901   Mean   :0.9836   Mean   :0.9656   Mean   :0.9852   Mean   :0.9844  
# 3rd Qu.:0.9958   3rd Qu.:0.9935   3rd Qu.:0.9709   3rd Qu.:0.9928   3rd Qu.:0.9898  
# Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
# M_46             M_47             M_48             M_49             M_50       
# Min.   :0.8923   Min.   :0.9156   Min.   :0.9205   Min.   :0.9061   Min.   :0.8853  
# 1st Qu.:0.9896   1st Qu.:0.9864   1st Qu.:0.9826   1st Qu.:0.9915   1st Qu.:0.9831  
# Median :0.9936   Median :0.9917   Median :0.9898   Median :0.9943   Median :0.9887  
# Mean   :0.9901   Mean   :0.9879   Mean   :0.9855   Mean   :0.9907   Mean   :0.9856  
# 3rd Qu.:0.9958   3rd Qu.:0.9944   3rd Qu.:0.9943   3rd Qu.:0.9962   3rd Qu.:0.9917  
# Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
# M_51             M_52             M_53             M_54             M_55       
# Min.   :0.9077   Min.   :0.9097   Min.   :0.9020   Min.   :0.9098   Min.   :0.9018  
# 1st Qu.:0.9853   1st Qu.:0.9872   1st Qu.:0.9877   1st Qu.:0.9843   1st Qu.:0.9908  
# Median :0.9901   Median :0.9914   Median :0.9925   Median :0.9901   Median :0.9937  
# Mean   :0.9871   Mean   :0.9886   Mean   :0.9890   Mean   :0.9869   Mean   :0.9900  
# 3rd Qu.:0.9927   3rd Qu.:0.9954   3rd Qu.:0.9959   3rd Qu.:0.9950   3rd Qu.:0.9959  
# Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
# M_56       
# Min.   :0.9001  
# 1st Qu.:0.9908  
# Median :0.9933  
# Mean   :0.9907  
# 3rd Qu.:0.9958  
# Max.   :1.0000  

cor_summary <- as.data.frame(cor_mat) %>%
  mutate(sample = rownames(.)) %>%
  pivot_longer(cols = -sample, names_to = "other_sample", values_to = "r") %>%
  filter(sample != other_sample) %>%    # remove self-correlations = 1
  group_by(sample) %>%
  summarise(
    min      = min(r),
    q1       = quantile(r, 0.25),
    median   = median(r),
    mean     = mean(r),
    q3       = quantile(r, 0.75),
    max      = max(r),
    .groups = "drop"
  ) %>%
  arrange(mean)

cor_summary %>% arrange(median)

head(cor_summary) # bottom 10 median cor. 
# sample   min    q1 median  mean    q3   max
# <chr>  <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl>
# 1 M_9    0.821 0.885  0.896 0.894 0.903 0.921
# 2 M_7    0.821 0.940  0.954 0.951 0.967 0.983
# 3 M_43   0.869 0.961  0.967 0.965 0.971 0.988
# 4 M_39   0.902 0.973  0.977 0.975 0.981 0.990
# 5 M_36   0.861 0.976  0.983 0.979 0.989 0.997
# 6 M_28   0.898 0.980  0.988 0.983 0.993 0.999

cor_long <- as.data.frame(cor_mat) %>%
  mutate(sample1 = rownames(.)) %>%
  pivot_longer(-sample1, names_to = "sample2", values_to = "r") %>%
  filter(sample1 != sample2)           # remove perfect 1.00 self-corr

ggplot(cor_long, aes(x = sample1, y = r)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  ylim(0.7, 1.00)+
  labs(
    title = "Correlation Distribution per Sample (Complete-Case Proteins)",
    x = "Sample",
    y = "Correlation (Pearson r)"
  )

global_stats <- cor_long %>%
  summarise(
    min_r = min(r),
    q1_r  = quantile(r, 0.25),
    med_r = median(r),
    q3_r  = quantile(r, 0.75),
    max_r = max(r)
  )
global_stats

# # A tibble: 1 × 5
# min_r  q1_r med_r  q3_r max_r
# <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 0.821 0.984 0.990 0.994 0.999


p_box<- ggplot(cor_long, aes(x = sample1, y = r)) +
  geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.4) +
  
  # horizontal lines
  geom_hline(aes(yintercept = global_stats$min_r),   linetype = "dotted", colour = "red") +
  geom_hline(aes(yintercept = global_stats$q1_r),    linetype = "dashed", colour = "red") +
  geom_hline(aes(yintercept = global_stats$med_r),   linetype = "solid",  colour = "red") +
  geom_hline(aes(yintercept = global_stats$q3_r),    linetype = "dashed", colour = "red") +

  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  ylim(0.7,1)+
  labs(
    x = "",
    y = "Correlation"
  )

# abundance distribution plot #

abd_means <- abd_clean %>%
  mutate(mean_abundance = rowMeans(across(all_of(sample_cols)), na.rm = TRUE))

rank_df <- abd_means %>%
  arrange(desc(mean_abundance)) %>%
  mutate(rank = row_number())

rank_df <- rank_df %>%
  mutate(
    label =
      case_when(
        rank <= 10 ~ accession,                     # Top 10
        rank > n() - 10 ~ accession,                # Bottom 10
        TRUE ~ NA_character_                               # No label
      )
  )

p_rank<-  ggplot(rank_df, aes(x = rank, y = mean_abundance)) +
  geom_point(size = 1, alpha = 0.7, colour = "grey40") +
  scale_y_log10() +
  
  # Boxed labels for top & bottom 10
  geom_label_repel(
    aes(label = label),
    na.rm = TRUE,
    size = 3,
    label.size = 0.20,          # border thickness
    fill = "white",
    colour = "black",
    label.r = unit(0.15, "lines"),
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf,         # allow unlimited spreading
    min.segment.length = 0,
    segment.size = 0.25,
    segment.color = "grey50",
    force = 5,
    force_pull = 1,
    seed = 10
  ) +
  
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    x = "Protein rank",
    y = "Mean abundance (log10 scale)"
  )

p_rank

# Make QC figure #
g_heat <- p_heat$gtable

p_heat_pw <- patchwork::wrap_elements(full = g_heat)

p_combined <- (p_bar + p_rank)/ (p_heat_pw + p_box) +
  plot_layout(widths = c(1, 1, 1.2, 1))

p_combined

ggsave(
  filename = file.path(dir_figs, "abd_myofib_qc_plots.pdf"),
  plot = p_combined,
  width = 9,
  height = 8,
  dpi = 600
)


