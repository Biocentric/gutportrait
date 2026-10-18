#!/usr/bin/env Rscript
#
# Positive control, injected signal (DECISIONS.md D12a): does batch correction preserve an
# effect we planted ourselves? If it erases one we put there, it would erase real ones too,
# and the run STOPS (exit status 1) with the diagnostic instead of producing results.
#
# 1. Pick a high-prevalence (>= --min-prevalence), mid-abundance genus in the panel.
# 2. Multiply its counts by --fold in a random half of each panel study (fixed seed).
# 3. Run the same correction as BATCH_CORRECT on the spiked and on the unspiked table.
# 4. Per spiked sample, the planted effect that survived correction is
#       CLR(spiked, corrected) - CLR(unspiked, corrected)
#    against the expected CLR shift log(fold) * (1 - 1/D) for D genera.
# PASS if the median surviving effect has the planted sign, a paired Wilcoxon signed-rank
# test gives p < --alpha, and the mean retained fraction is >= --retain.
#
# The test is paired on purpose: comparing spiked against unspiked samples would mostly
# measure natural person-to-person variation and could fail even if correction were perfect.
#
# usage: positive_control.R --counts GENUS_COUNTS.tsv --samplesheet SHEET.csv --prefix P
#        [--fold 2] [--alpha 0.01] [--retain 0.5] [--min-prevalence 0.9]
#        [--clr-prevalence 0.1] [--pseudocount 0.5] [--seed 20261018]

suppressPackageStartupMessages(library(MMUPHin))

args <- commandArgs(trailingOnly = TRUE)
opt <- list(counts = NULL, samplesheet = NULL, prefix = "out", fold = 2, alpha = 0.01, retain = 0.5,
            min_prevalence = 0.9, clr_prevalence = 0.1, pseudocount = 0.5, seed = 20261018)
i <- 1
while (i <= length(args)) {
    key <- gsub("-", "_", sub("^--", "", args[i]))
    if (!key %in% names(opt)) stop("unknown option ", args[i])
    opt[[key]] <- if (is.numeric(opt[[key]])) as.numeric(args[i + 1]) else args[i + 1]
    i <- i + 2
}

counts <- as.matrix(read.delim(opt$counts, check.names = FALSE, row.names = 1))
sheet <- read.csv(opt$samplesheet, colClasses = "character", na.strings = character(0))
sheet[is.na(sheet)] <- ""
panel <- sheet[sheet$role == "panel", , drop = FALSE]
rownames(panel) <- panel$sample
panel$study <- factor(panel$study)
x <- counts[, panel$sample, drop = FALSE]

result <- function(verdict, ...) {
    r <- data.frame(verdict = verdict, ..., fold = opt$fold, alpha = opt$alpha, retain_required = opt$retain,
                    seed = opt$seed, check.names = FALSE)
    write.table(r, sprintf("%s.positive_control.tsv", opt$prefix), sep = "\t", quote = FALSE, row.names = FALSE)
    r
}

study_n <- table(panel$study)
if (length(study_n) < 2 || any(study_n < 2)) {
    result("NOT_APPLICABLE", genus = "", n_spiked = 0, expected_clr_shift = NA, median_surviving = NA,
           retained = NA, p = NA, note = "fewer than two panel studies: no correction to test")
    cat("positive control: not applicable (no batch correction with fewer than two panel studies)\n")
    quit(status = 0)
}

rel <- sweep(x, 2, colSums(x), "/")
prev <- rowMeans(x > 0)
cand <- setdiff(names(prev)[prev >= opt$min_prevalence], "Unassigned")
if (length(cand) == 0) stop("no genus reaches the prevalence needed for the positive control")
med <- sort(apply(rel[cand, , drop = FALSE], 1, median))
genus <- names(med)[ceiling(length(med) / 2)]              # mid-abundance among the common genera

set.seed(opt$seed)
spiked <- unlist(lapply(split(panel$sample, panel$study), function(s) sample(s, floor(length(s) / 2))))
xs <- x
xs[genus, spiked] <- round(xs[genus, spiked] * opt$fold)

correct <- function(m) adjust_batch(feature_abd = m, batch = "study", data = panel,
                                    control = list(verbose = FALSE))$feature_abd_adj
keep <- rownames(x)[prev >= opt$clr_prevalence]
clr <- function(m) { l <- log(m[keep, , drop = FALSE] + opt$pseudocount); sweep(l, 2, colMeans(l)) }
c_spiked <- clr(correct(xs))
c_plain <- clr(correct(x))

surviving <- c_spiked[genus, spiked] - c_plain[genus, spiked]
expected <- log(opt$fold) * (1 - 1 / length(keep))
retained <- mean(surviving) / expected
p <- suppressWarnings(wilcox.test(surviving, mu = 0, alternative = "two.sided")$p.value)
right_sign <- sign(median(surviving)) == sign(log(opt$fold))
pass <- right_sign && p < opt$alpha && retained >= opt$retain

r <- result(if (pass) "PASS" else "FAIL", genus = genus, n_spiked = length(spiked),
            expected_clr_shift = round(expected, 4), median_surviving = round(median(surviving), 4),
            retained = round(retained, 3), p = signif(p, 3),
            note = if (pass) "planted effect survives correction" else
                "correction erased or distorted a planted effect: results would not be trustworthy")
msg <- sprintf("positive control %s: %s x%g in %d panel samples; %.0f%% of the planted CLR shift retained, paired Wilcoxon p = %.2g",
               r$verdict, genus, opt$fold, length(spiked), 100 * retained, p)
if (!pass) {
    message("STOP: ", msg)
    message("The batch correction removed an effect planted on purpose, so it would also remove real ones.")
    quit(status = 1)
}
cat(msg, "\n")
