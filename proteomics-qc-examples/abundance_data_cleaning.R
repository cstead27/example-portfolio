# ---------------------------------------------------------------
# clean_abundance_data()
# Purpose:
#   Cleans Progenesis/Progenesis-like protein abundance exports:
#   - removes unwanted metadata columns
#   - promotes second row to column names
#   - extracts GN= gene names
#   - standardises accession fields
#   - returns fully tidy, analysis-ready tibble
# ---------------------------------------------------------------
clean_abundance_data <- function(path) {
  # ---- Load raw ----
  raw <- readr::read_csv(path, col_names = FALSE, show_col_types = FALSE)
  
  # ---- Identify header rows ----
  header_row1 <- as.character(raw[1, ])  # marker row (contains "Raw abundance")
  header_row3 <- as.character(raw[3, ])  # actual headers
  
  # ---- Find cut point ("Raw abundance") ----
  cut_idx <- which(trimws(header_row1) == "Raw abundance")
  if (length(cut_idx) == 0) {
    stop("'Raw abundance' not found in first row of file.")
  }
  
  # ---- Subset columns before the Raw abundance marker ----
  dat <- raw[, 1:(cut_idx - 1)]
  
  # ---- Promote row 3 to column names ----
  colnames(dat) <- header_row3[1:(cut_idx - 1)]
  
  # ---- Drop first three rows (junk, condition, header) ----
  dat <- dat[-c(1, 2, 3), ]
  
  # ---- Convert to tibble and infer types ----
  dat <- dat %>%
    tibble::as_tibble() %>%
    readr::type_convert()
  
  # ---- Clean column names only for metadata fields ----
  # We'll temporarily clean all, but restore sample names later.
  clean_names_all <- janitor::make_clean_names(names(dat))
  
  # Save a map to restore sample columns exactly
  name_map <- tibble(old = names(dat), new = clean_names_all)
  names(dat) <- clean_names_all
  
  # ---- Remove Progenesis metadata columns ----
  cols_to_remove <- c(
    "confidence_score", "anova_p", "q_value", "max_fold_change", "power",
    "highest_mean_condition", "lowest_mean_condition"
  )
  dat <- dat %>% select(-any_of(cols_to_remove))
  
  # ---- Extract gene_name from GN= ----
  dat <- dat %>%
    mutate(
      gene_name = stringr::str_extract(description, "(?<=GN=)[A-Za-z0-9_.-]+")
    )
  
  # ---- Clean description text ----
  dat <- dat %>%
    mutate(
      description = description %>%
        stringr::str_remove_all(
          "OS=[^ ]+|OX=[0-9]+|GN=[A-Za-z0-9_.-]+|PE=[0-9]+|SV=[0-9]+"
        ) %>%
        stringr::str_squish()
    )
  
  # ---- Standardise accession fields ----
  dat <- dat %>%
    relocate(gene_name, .after = accession) %>%
    mutate(
      accession = purrr::map_chr(stringr::str_split(accession, ";"), ~ .x[1]),
      accession_full  = accession,
      accession = stringr::str_remove(accession, "_HUMAN$")
    ) %>%
    relocate(accession_full, accession, .before = gene_name)
  
  # ---- Restore original sample column names ----
  # Identify columns that are *not* metadata
  meta_cols <- c("accession_full", "accession",
                 "gene_name", "description", "peptide_count", "unique_peptides", "mass")
  sample_cols <- setdiff(names(dat), meta_cols)
  
  # Map cleaned → original names for sample columns only
  restore_map <- name_map %>%
    filter(new %in% sample_cols)
  names(dat)[match(restore_map$new, names(dat))] <- restore_map$old
  
  # ---- Reorder metadata columns to the front ----
  dat <- dat %>% relocate(all_of(meta_cols))
  
  return(dat)
}