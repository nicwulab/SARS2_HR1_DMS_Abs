library(data.table)
library(dplyr)
esc <- as.data.frame(read.table(file = 'result/S2HR1_bind_scores.tsv', sep = '\t', header = TRUE))
expr <- as.data.frame(read.table(file = 'data/S2HR1_exp_fus_scores.tsv', sep = '\t', header = TRUE))

# Only include mutants that appear in both escape and expression/fusion experiments.
# Append expression score for mutants found in 'esc'. Discard mutants in 'esc' that do not have expression scores.
common_muts <- intersect(esc$mut, expr$mut)

esc <- esc[esc$mut %in% common_muts, ]

# Reorder mutants.
esc <- esc[match(common_muts, esc$mut), ]
expr <- expr[match(common_muts, expr$mut), ]

esc$exp_score <- expr$exp_score
esc$fus_score <- expr$fus_score
esc <- data.table(esc)
bind_score_avg <- rowMeans(select(esc, A107_score, A214_score, A218_score))
residual <- (lm(bind_score_avg~esc$exp_score)$residual)
esc$bind_score_avg <- bind_score_avg
esc$residual <- residual

write.table(esc, file='result/S2HR1_scores_common.tsv', quote=FALSE, sep='\t', col.names = NA)
