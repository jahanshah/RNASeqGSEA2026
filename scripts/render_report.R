#!/usr/bin/env Rscript
# Render the final HTML report

user_lib <- path.expand("~/R/library")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

options(repos = c(CRAN = "https://cloud.r-project.org"))
for (p in c("rmarkdown","knitr","kableExtra","DT","ggridges")) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, lib = user_lib)
}

BASE <- "/Users/jahanshah/Documents/Consultant-NGS/Ahmed/Projects2026/Project-RNAseq/RNASeqGSEA2026"

rmarkdown::render(
  input       = file.path(BASE, "scripts", "GSEA_report.Rmd"),
  output_file = "GSEA_final_report.html",
  output_dir  = BASE,
  envir       = new.env(),
  quiet       = FALSE
)

message("Report saved to: ", file.path(BASE, "GSEA_final_report.html"))
