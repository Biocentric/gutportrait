#!/usr/bin/env Rscript
#
# Harmonise the reference panel's studies with each other (MMUPHin::adjust_batch).
#
# ONLY THE PANEL IS CORRECTED. The personal samples form their own batch, and three samples
# cannot separate their lab from their biology: correcting that batch would shrink them
# toward the panel and erase exactly the deviation the pipeline exists to measure. They are
# passed through unchanged. Correction removes study-to-study shifts WITHIN the panel, so
# the reference distribution reflects people rather than a mixture of lab protocols.
#
# Downstream, every position is computed on both the corrected and the uncorrected panel;
# anything that flips between them is reported as unstable, not as a finding.
#
# usage: batch_correct.R --counts GENUS_COUNTS.tsv --samplesheet SHEET.csv --prefix P
#        [--covariates age_cat,sex]

suppressPackageStartupMessages(library(MMUPHin))

args <- commandArgs(trailingOnly = TRUE)
opt <- list(counts = NULL, samplesheet = NULL, prefix = "out", covariates = "")
i <- 1
while (i <= length(args)) {
    key <- gsub("-", "_", sub("^--", "", args[i]))
    if (!key %in% names(opt)) stop("unknown option ", args[i])
    opt[[key]] <- args[i + 1]
    i <- i + 2
}
if (is.null(opt$counts) || is.null(opt$samplesheet)) stop("--counts and --samplesheet are required")

counts <- as.matrix(read.delim(opt$counts, check.names = FALSE, row.names = 1))
sheet <- read.csv(opt$samplesheet, colClasses = "character", na.strings = character(0))
sheet[is.na(sheet)] <- ""
# A set of bands ("30s|40s", a study's inclusion criteria) says nothing about one person:
# it is unknown for every per-person model here.
for (v in intersect(c("age_cat", "sex", "bmi_cat"), colnames(sheet))) {
    sheet[[v]][grepl("|", sheet[[v]], fixed = TRUE)] <- ""
}
panel <- sheet[sheet$role == "panel", , drop = FALSE]
rownames(panel) <- panel$sample
panel$study <- factor(panel$study)

study_n <- table(panel$study)
wanted <- strsplit(opt$covariates, ",")[[1]]
wanted <- wanted[wanted != ""]
# MMUPHin needs complete covariates; use only those filled for every panel sample.
covs <- wanted[vapply(wanted, function(v) v %in% colnames(panel) && all(panel[[v]] != ""), logical(1))]
dropped <- setdiff(wanted, covs)

corrected <- counts
if (length(study_n) >= 2 && all(study_n >= 2)) {
    fit <- adjust_batch(feature_abd = counts[, panel$sample, drop = FALSE], batch = "study",
                        covariates = if (length(covs)) covs else NULL, data = panel,
                        control = list(verbose = FALSE))
    corrected[, panel$sample] <- fit$feature_abd_adj
    note <- "panel studies harmonised; personal samples unchanged"
} else {
    note <- if (length(study_n) < 2) "single panel study: nothing to correct, table passed through" else
        sprintf("a panel study has fewer than 2 samples (%s): not corrected",
                paste(names(study_n)[study_n < 2], collapse = ", "))
}

out <- data.frame(Genus = rownames(corrected), corrected, check.names = FALSE)
write.table(out, sprintf("%s.genus_counts_corrected.tsv", opt$prefix), sep = "\t", quote = FALSE, row.names = FALSE)
summary <- data.frame(method = "MMUPHin::adjust_batch", batch = "study", panel_studies = length(study_n),
    panel_samples = nrow(panel), covariates = paste(covs, collapse = ","),
    covariates_dropped_incomplete = paste(dropped, collapse = ","), note = note)
write.table(summary, sprintf("%s.correction_summary.tsv", opt$prefix), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("batch correction: %d panel samples in %d studies; %s\n", nrow(panel), length(study_n), note))
