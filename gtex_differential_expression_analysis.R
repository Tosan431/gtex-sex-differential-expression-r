# ==============================================================
# GTEx Sex Differential Expression Analysis in R
# Employer-facing portfolio version
# Author: Tosan Akpituren
# MSc Health Data Science for Applied Precision Medicine
# University of Dundee
# ==============================================================

# --------------------------------------------------------------
# 0. Package setup
# --------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "ggplot2",
  "tibble",
  "stringr"
)

installed <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!(pkg %in% installed)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Optional packages used in the original dissertation workflow:
# limma, biomaRt, clusterProfiler, DOSE
# These may require Bioconductor:
# install.packages("BiocManager")
# BiocManager::install(c("limma", "biomaRt", "clusterProfiler", "DOSE"))

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# --------------------------------------------------------------
# 1. Expected input structure
# --------------------------------------------------------------
# This script is a clean portfolio version of the original workflow.
# Large GTEx source files are not included in this GitHub repository.
#
# Expected inputs when reproducing from source:
# - tissue expression matrix: rows = samples, columns = genes
# - phenotype/covariate file containing SUBJID, SEX, AGE, DTHHRDY
# - genomic annotation file containing gene ID, chromosome and position
#
# Example object names used below:
# expression_df      - gene expression matrix plus sample identifier
# phenotype_df       - sample phenotype/covariate table
# annotation_df      - gene annotation table
# merged_data        - merged expression + phenotype table
# results_df         - differential expression results

# --------------------------------------------------------------
# 2. Data cleaning and merge template
# --------------------------------------------------------------

clean_subject_id <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_trim()
}

# Example:
# expression_df <- read_csv("data/aorta_expression.csv")
# phenotype_df  <- read_csv("data/aorta_phenotype.csv")
#
# expression_df <- expression_df %>%
#   mutate(SUBJID = clean_subject_id(SUBJID))
#
# phenotype_df <- phenotype_df %>%
#   mutate(SUBJID = clean_subject_id(SUBJID))
#
# merged_data <- expression_df %>%
#   inner_join(phenotype_df, by = "SUBJID")

# --------------------------------------------------------------
# 3. Differential expression function
# --------------------------------------------------------------

run_sex_differential_expression <- function(merged_data,
                                            gene_columns,
                                            sex_column = "SEX") {

  sex_vector <- merged_data[[sex_column]]

  output <- lapply(gene_columns, function(gene) {

    model_df <- data.frame(
      expression = merged_data[[gene]],
      sex = sex_vector
    )

    model_df <- model_df[complete.cases(model_df), ]

    fit <- lm(expression ~ sex, data = model_df)
    fit_summary <- summary(fit)

    beta <- coef(fit)[["sex"]]
    p_value <- coef(fit_summary)[2, "Pr(>|t|)"]

    data.frame(
      Gene = gene,
      Beta = beta,
      P.Value = p_value
    )
  })

  results <- bind_rows(output)

  results <- results %>%
    mutate(
      BH_Adjusted_P = p.adjust(P.Value, method = "BH"),
      Bonferroni_Adjusted_P = p.adjust(P.Value, method = "bonferroni"),
      logP = -log10(P.Value)
    )

  return(results)
}

# Example:
# gene_columns <- setdiff(names(merged_data), c("SUBJID", "SEX", "AGE", "DTHHRDY"))
# results_df <- run_sex_differential_expression(merged_data, gene_columns)

# --------------------------------------------------------------
# 4. Significant gene filtering
# --------------------------------------------------------------

filter_significant_genes <- function(results_df,
                                     method = c("BH", "Bonferroni"),
                                     threshold = 0.05) {

  method <- match.arg(method)

  if (method == "BH") {
    filtered <- results_df %>%
      filter(BH_Adjusted_P < threshold)
  } else {
    filtered <- results_df %>%
      filter(Bonferroni_Adjusted_P < threshold)
  }

  return(filtered)
}

# Example:
# significant_bh <- filter_significant_genes(results_df, method = "BH")
# significant_bonferroni <- filter_significant_genes(results_df, method = "Bonferroni")

# --------------------------------------------------------------
# 5. Volcano plot
# --------------------------------------------------------------

create_volcano_plot <- function(results_df,
                                title = "Volcano Plot of Differential Expression",
                                output_file = "figures/volcano_plot.png") {

  p <- ggplot(results_df, aes(x = Beta, y = logP)) +
    geom_point(alpha = 0.6, size = 1.4) +
    geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed") +
    geom_hline(yintercept = 1.3, linetype = "dashed") +
    labs(
      title = title,
      x = "Log Fold Change (Beta)",
      y = "-log10(P-value)"
    ) +
    theme_minimal(base_size = 13)

  ggsave(output_file, p, width = 11, height = 6.5, dpi = 300)
  return(p)
}

# Example:
# create_volcano_plot(results_df, "Volcano Plot - Breast Mammary Tissue")

# --------------------------------------------------------------
# 6. Manhattan plot
# --------------------------------------------------------------

create_manhattan_plot <- function(manhattan_df,
                                  title = "Manhattan Plot of Significant Genomic Associations",
                                  output_file = "figures/manhattan_plot.png") {

  p <- ggplot(manhattan_df, aes(x = BP, y = -log10(P), colour = as.factor(CHR))) +
    geom_point(alpha = 0.75, size = 1.2) +
    labs(
      title = title,
      x = "Genomic Position",
      y = "-log10(P-value)",
      colour = "Chromosome"
    ) +
    theme_minimal(base_size = 13)

  ggsave(output_file, p, width = 12, height = 7, dpi = 300)
  return(p)
}

# Expected manhattan_df columns:
# SNP / Gene ID, CHR, BP, P

# --------------------------------------------------------------
# 7. P-value distribution
# --------------------------------------------------------------

create_pvalue_histogram <- function(results_df,
                                    output_file = "figures/pvalue_histogram.png") {

  p <- ggplot(results_df, aes(x = P.Value)) +
    geom_histogram(bins = 60, fill = "steelblue", colour = "white") +
    labs(
      title = "Distribution of Gene Expression P-values",
      x = "P-value",
      y = "Frequency"
    ) +
    theme_minimal(base_size = 13)

  ggsave(output_file, p, width = 10, height = 6, dpi = 300)
  return(p)
}

# --------------------------------------------------------------
# 8. Gene type distribution
# --------------------------------------------------------------

create_gene_type_distribution <- function(annotated_results,
                                          output_file = "figures/gene_type_distribution.png") {

  gene_type_counts <- as.data.frame(table(annotated_results$type))
  names(gene_type_counts) <- c("Gene_Type", "Count")

  gene_type_counts <- gene_type_counts %>%
    arrange(desc(Count)) %>%
    slice_head(n = 12)

  p <- ggplot(gene_type_counts, aes(x = reorder(Gene_Type, Count), y = Count)) +
    geom_col(fill = "darkgreen") +
    coord_flip() +
    labs(
      title = "Top Gene Types in Annotated Results",
      x = "Gene Type",
      y = "Number of Genes"
    ) +
    theme_minimal(base_size = 13)

  ggsave(output_file, p, width = 10, height = 7, dpi = 300)
  return(p)
}

# --------------------------------------------------------------
# 9. KEGG / DAVID pathway visualisation
# --------------------------------------------------------------

create_pathway_barplot <- function(pathway_df,
                                   output_file = "figures/kegg_pathway_enrichment.png") {

  p <- ggplot(pathway_df, aes(x = reorder(Description, -log10(p.adjust)),
                              y = -log10(p.adjust))) +
    geom_col(fill = "darkorange") +
    coord_flip() +
    labs(
      title = "Top Enriched KEGG Pathways",
      x = "KEGG Pathway",
      y = "-log10 Adjusted P-value"
    ) +
    theme_minimal(base_size = 13)

  ggsave(output_file, p, width = 11, height = 7, dpi = 300)
  return(p)
}

# --------------------------------------------------------------
# 10. Fisher exact test template for GWAS overlap
# --------------------------------------------------------------

run_gwas_overlap_test <- function(total_genes,
                                  significant_genes,
                                  gwas_genes) {

  total_genes <- unique(total_genes)
  significant_genes <- unique(significant_genes)
  gwas_genes <- unique(gwas_genes)

  a <- length(intersect(significant_genes, gwas_genes))
  b <- length(setdiff(significant_genes, gwas_genes))
  c <- length(intersect(setdiff(total_genes, significant_genes), gwas_genes))
  d <- length(setdiff(setdiff(total_genes, significant_genes), gwas_genes))

  contingency_table <- matrix(
    c(a, b, c, d),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      DE_Status = c("DEG", "Not_DEG"),
      GWAS_Status = c("In_GWAS", "Not_In_GWAS")
    )
  )

  fisher_result <- fisher.test(contingency_table)

  return(list(
    contingency_table = contingency_table,
    fisher_result = fisher_result
  ))
}

# --------------------------------------------------------------
# 11. Export results
# --------------------------------------------------------------

# Example:
# write_csv(results_df, "results/all_differential_expression_results.csv")
# write_csv(significant_bh, "results/significant_genes_bh.csv")
# write_csv(significant_bonferroni, "results/significant_genes_bonferroni.csv")

# --------------------------------------------------------------
# End of script
# --------------------------------------------------------------
