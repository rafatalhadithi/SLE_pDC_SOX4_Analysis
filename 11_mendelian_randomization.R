# 11_mendelian_randomization.R
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ggplot2)
  library(dplyr)
})

# Authenticate OpenGWAS API Token
user_jwt <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJzYy5yX2FsaGFkaXRoaUB1b2FuYmFyLmVkdS5pcSIsImlhdCI6MTc5MzYxOTg1NiwiZXhwIjoxNzk0ODI5NDU2fQ.HrADPccgVmEnekBt568Wkro1DytkPFrkDZ7AEuW0nMfSyS7FHONsuc-jUJbe8uQziHFp5ZkIqIHnn5_4zC_jtbJXv53HYmIxDhRZtOlnZeleu8FFAFAMlr82H4yL11Z6rYE1Y2TQY2QggLPyRMH1coH-wakHb6sBvUDi_EowzHoa8XAhGa3qRklL3hmfXXIF0MFQJ2-VKW0DrC3Db_6ft2rmyZgSGCNpbLj6svWtN6krc-wDP4QNPk0JYHLr3bb9vNajtaN336l3nvNBqybW2P86YqKrTKAiBmAZc_oTcJ8fFcNLGdWveSVyZJzOleQJ2YFBa91gT7WiGvh61lTUEg"
Sys.setenv(OPENGWAS_JWT = user_jwt)
Sys.setenv(OPENGWAS_TOKEN = user_jwt)

cat("--- Step 1: Initializing Two-Sample MR Pipeline ---\n")
output_dir <- "figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# eQTLGen Ensembl IDs for target network
exposures <- c("eqtl-a-ENSG00000124766", "eqtl-a-ENSG00000121966", "eqtl-a-ENSG00000185507")
gene_map <- c(
  "eqtl-a-ENSG00000124766" = "SOX4", 
  "eqtl-a-ENSG00000121966" = "CXCR4", 
  "eqtl-a-ENSG00000185507" = "IRF7"
)

cat("Fetching genetic instruments (eQTLs) for SOX4, CXCR4, and IRF7...\n")
exp_dat <- extract_instruments(outcomes = exposures, p1 = 5e-6, clump = TRUE)

if (is.null(exp_dat) || nrow(exp_dat) == 0) {
  stop("CRITICAL: No genetic instruments found for these genes on the OpenGWAS server.")
}

exp_dat$exposure <- gene_map[exp_dat$id.exposure]
cat(sprintf("Successfully extracted %d genetic variants across the network.\n", nrow(exp_dat)))

cat("\n--- Step 2: Fetching True SLE GWAS Outcome Data ---\n")
outcome_id <- "ebi-a-GCST003156" 
cat(sprintf("Querying Systemic Lupus Erythematosus GWAS (%s) for overlapping variants...\n", outcome_id))
out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = outcome_id)

if (is.null(out_dat) || nrow(out_dat) == 0) {
  stop("CRITICAL: None of the exposure SNPs were found in the outcome GWAS dataset.")
}

cat("\n--- Step 3: Harmonizing DNA Strands ---\n")
dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat)

cat("\n--- Step 4: Running Causal Inference & Sensitivity Tests ---\n")
res <- mr(dat)
res_or <- generate_odds_ratios(res)

het <- mr_heterogeneity(dat)
pleio <- tryCatch({ mr_pleiotropy_test(dat) }, error = function(e) { NULL })

cat("\n=== PRIMARY MR CAUSAL ESTIMATES ===\n")
print(res_or[, c("exposure", "method", "nsnp", "or", "or_lci95", "or_uci95", "pval")])

if (!is.null(het) && nrow(het) > 0) {
  cat("\n=== HETEROGENEITY TEST (Cochran's Q) ===\n")
  print(het[, c("exposure", "method", "Q", "Q_df", "Q_pval")])
  write.csv(het, file = file.path(output_dir, "Table_MR_Heterogeneity.csv"), row.names = FALSE)
}

if (!is.null(pleio) && nrow(pleio) > 0 && "eg_intercept" %in% colnames(pleio)) {
  cat("\n=== HORIZONTAL PLEIOTROPY TEST (MR-Egger Intercept) ===\n")
  print(pleio[, c("exposure", "eg_intercept", "se", "pval")])
  write.csv(pleio, file = file.path(output_dir, "Table_MR_Pleiotropy.csv"), row.names = FALSE)
} else {
  cat("\n=== HORIZONTAL PLEIOTROPY TEST ===\n")
  cat("Note: Insufficient SNPs (N < 3) on some exposures to fit MR-Egger intercept regression.\n")
}

write.csv(res_or, file = file.path(output_dir, "Table_MR_Causal_Estimates.csv"), row.names = FALSE)

cat("\n--- Step 5: Generating Q1 Publication Plots (Cleaned) ---\n")
# 1. Scatter Plots
p_scatter <- mr_scatter_plot(res, dat)
for (key in names(p_scatter)) {
  matched_gene <- "Network_Gene"
  for (id in names(gene_map)) {
    if (grepl(id, key) || grepl(gene_map[id], key)) {
      matched_gene <- gene_map[id]
      break
    }
  }
  
  # Override default aesthetics
  p <- p_scatter[[key]] + 
    theme_classic(base_size = 14) +
    labs(
      title = sprintf("%s Causality in SLE", matched_gene),
      x = sprintf("Genetic Effect on %s Expression", matched_gene),
      y = "Genetic Effect on SLE Risk",
      color = "MR Method"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90", linetype = "dashed")
    )
    
  ggsave(file.path(output_dir, sprintf("Figure_7A_MR_Scatter_%s.pdf", matched_gene)), plot = p, width = 8, height = 7)
  ggsave(file.path(output_dir, sprintf("Figure_7A_MR_Scatter_%s.png", matched_gene)), plot = p, width = 8, height = 7, dpi = 300)
}

# 2. Forest Plots
res_single <- mr_singlesnp(dat)
p_forest <- mr_forest_plot(res_single)
for (key in names(p_forest)) {
  matched_gene <- "Network_Gene"
  for (id in names(gene_map)) {
    if (grepl(id, key) || grepl(gene_map[id], key)) {
      matched_gene <- gene_map[id]
      break
    }
  }
  
  # Override default aesthetics
  p <- p_forest[[key]] + 
    theme_classic(base_size = 14) +
    labs(
      title = sprintf("Causal Estimate of %s eQTLs on SLE", matched_gene),
      x = "MR Effect Size on SLE",
      y = "Genetic Variant (SNP)"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title = element_text(face = "bold"),
      legend.position = "none",
      panel.grid.major.x = element_line(color = "grey90", linetype = "dashed")
    )
    
  ggsave(file.path(output_dir, sprintf("Figure_7B_MR_Forest_%s.pdf", matched_gene)), plot = p, width = 8, height = 8)
  ggsave(file.path(output_dir, sprintf("Figure_7B_MR_Forest_%s.png", matched_gene)), plot = p, width = 8, height = 8, dpi = 300)
}

cat(sprintf("\nSUCCESS: Mendelian Randomization complete! High-resolution figures saved to: %s\n", output_dir))