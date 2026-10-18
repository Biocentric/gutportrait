#!/usr/bin/env Rscript
#
# How much of the panel's variation is study, and how much is biology?
#
# PERMANOVA (vegan::adonis2) on Aitchison distance: Euclidean distance between
# centred-log-ratio transformed genus profiles, the distance that respects
# compositionality. Panel samples only, with the same panel-only prevalence filter and
# pseudocount as the percentiles, so both modules look at the same genera.
#
# Models (DECISIONS.md D12):
#   study                 R2(study) for the global banner (orientation only if > threshold)
#   study + <variable>    for every metadata variable with >= 2 levels: marginal R2 for
#                         both terms. A lifestyle contrast counts as a finding only if
#                         R2(variable) > R2(study). Where the variable's samples come from
#                         a single study, study cannot be separated from it; that model is
#                         fitted without study and marked as confounded.
#
# usage: batch_diagnose.R --counts GENUS_COUNTS.tsv --samplesheet SHEET.csv --prefix P
#        [--min-prevalence 0.1] [--pseudocount 0.5] [--permutations 999] [--seed 20261018]
#        [--banner-r2 0.10] [--min-n 10]

suppressPackageStartupMessages(library(vegan))

args <- commandArgs(trailingOnly = TRUE)
opt <- list(counts = NULL, samplesheet = NULL, prefix = "out", min_prevalence = 0.1,
            pseudocount = 0.5, permutations = 999, seed = 20261018, banner_r2 = 0.10, min_n = 10)
i <- 1
while (i <= length(args)) {
    key <- gsub("-", "_", sub("^--", "", args[i]))
    if (!key %in% names(opt)) stop("unknown option ", args[i])
    val <- args[i + 1]
    opt[[key]] <- if (is.numeric(opt[[key]])) as.numeric(val) else val
    i <- i + 2
}
if (is.null(opt$counts) || is.null(opt$samplesheet)) stop("--counts and --samplesheet are required")

counts <- read.delim(opt$counts, check.names = FALSE, row.names = 1)
sheet <- read.csv(opt$samplesheet, colClasses = "character", na.strings = character(0))
sheet[is.na(sheet)] <- ""
# A set of bands ("30s|40s", a study's inclusion criteria) says nothing about one person:
# it is unknown for every per-person model here.
for (v in intersect(c("age_cat", "sex", "bmi_cat"), colnames(sheet))) {
    sheet[[v]][grepl("|", sheet[[v]], fixed = TRUE)] <- ""
}
panel <- sheet[sheet$role == "panel", , drop = FALSE]
missing <- setdiff(panel$sample, colnames(counts))
if (length(missing) > 0) stop("panel samples missing from the counts table: ", paste(head(missing), collapse = ", "))

x <- t(as.matrix(counts[, panel$sample, drop = FALSE]))          # samples x genera
prevalence <- colMeans(x > 0)
x <- x[, prevalence >= opt$min_prevalence, drop = FALSE]
logx <- log(x + opt$pseudocount)
clr <- logx - rowMeans(logx)
d_all <- dist(clr)                                                 # Aitchison distance
rownames(panel) <- panel$sample

rows <- list()
add <- function(model, term, df, r2, f, p, n, note) {
    rows[[length(rows) + 1]] <<- data.frame(model = model, term = term, df = df,
        R2 = round(r2, 4), F = round(f, 3), p = p, n = n, note = note)
}
fit <- function(formula, data, by) {
    set.seed(opt$seed)
    adonis2(formula, data = data, permutations = opt$permutations, by = by)
}
sub_dist <- function(keep) as.dist(as.matrix(d_all)[keep, keep])

n_studies <- length(unique(panel$study))
r2_study <- NA
if (n_studies >= 2) {
    d <- d_all
    a <- fit(d ~ study, panel, "terms")
    r2_study <- a["study", "R2"]
    add("study", "study", a["study", "Df"], r2_study, a["study", "F"], a["study", "Pr(>F)"],
        nrow(panel), "")
} else {
    add("study", "study", NA, NA, NA, NA, nrow(panel),
        "single panel study: between-study variance cannot be estimated")
}

for (v in intersect(c("age_cat", "sex", "bmi_cat", "country"), colnames(panel))) {
    keep <- panel$sample[panel[[v]] != ""]
    sub <- panel[keep, , drop = FALSE]
    if (length(keep) < opt$min_n || length(unique(sub[[v]])) < 2) {
        add(paste("study +", v), v, NA, NA, NA, NA, length(keep),
            sprintf("not testable: %d samples with %s, %d level(s)", length(keep), v, length(unique(sub[[v]]))))
        next
    }
    d <- sub_dist(keep)
    # Collinear with study: every level of the variable belongs to exactly one study and
    # every study to one level (e.g. country when each study is from one country). The
    # marginal test is then undefined (df 0); the variable cannot be separated from study.
    n_pairs <- nrow(unique(sub[, c("study", v)]))
    collinear <- n_pairs == length(unique(sub$study)) && n_pairs == length(unique(sub[[v]]))
    if (collinear) {
        a <- fit(as.formula(paste("d ~", v)), sub, "terms")
        add(v, v, a[v, "Df"], a[v, "R2"], a[v, "F"], a[v, "Pr(>F)"], length(keep),
            sprintf("confounded: %s is collinear with study (one level per study)", v))
    } else if (length(unique(sub$study)) >= 2) {
        a <- fit(as.formula(paste("d ~ study +", v)), sub, "margin")
        for (term in c("study", v)) {
            add(paste("study +", v), term, a[term, "Df"], a[term, "R2"], a[term, "F"], a[term, "Pr(>F)"],
                length(keep), "marginal")
        }
    } else {
        a <- fit(as.formula(paste("d ~", v)), sub, "terms")
        add(v, v, a[v, "Df"], a[v, "R2"], a[v, "F"], a[v, "Pr(>F)"], length(keep),
            sprintf("confounded: all samples with %s come from study %s", v, sub$study[1]))
    }
}

res <- do.call(rbind, rows)
write.table(res, sprintf("%s.permanova.tsv", opt$prefix), sep = "\t", quote = FALSE, row.names = FALSE)

banner <- !is.na(r2_study) && r2_study > opt$banner_r2
summary <- data.frame(panel_samples = nrow(panel), panel_studies = n_studies,
    genera_tested = ncol(x), r2_study = round(r2_study, 4), banner_threshold = opt$banner_r2,
    study_banner = banner, permutations = opt$permutations, seed = opt$seed,
    pseudocount = opt$pseudocount, min_prevalence = opt$min_prevalence)
write.table(summary, sprintf("%s.batch_summary.tsv", opt$prefix), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("PERMANOVA on %d panel samples in %d studies, %d genera: R2(study) = %s%s\n",
    nrow(panel), n_studies, ncol(x), format(round(r2_study, 4)),
    if (banner) "  -> ORIENTATION-ONLY BANNER" else ""))
