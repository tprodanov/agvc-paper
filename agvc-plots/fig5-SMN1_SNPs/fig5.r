#!/usr/bin/env Rscript

pdf(NULL)

suppressMessages(library(ggplot2))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(stringr))
suppressMessages(library(colorspace))

Sys.setenv('VROOM_CONNECTION_SIZE' = 131072 * 10)

add_comma <- function(s) { format(as.numeric(s), big.mark = ',') }

load_vcf <- function(filename) {
    vcf <- suppressMessages(read_delim(filename, '\t', comment = '##'))
    colnames(vcf)[1] <- 'CHROM'
    vcf
}

get_ac_matrix <- function(vcf) {
    gt_ix <- strsplit(vcf$FORMAT, ':') |>
        sapply(function(x) {
            i <- which(x == 'GT')
            if (length(i) == 0) { NA } else { i }
        })
    gq_ix <- strsplit(vcf$FORMAT, ':') |>
        sapply(function(x) {
            i <- which(x == 'GQ')
            if (length(i) == 0) { NA } else { i }
        })
    
    samples <- colnames(vcf)[10:ncol(vcf)]
    ac_matrix <- sapply(samples, function(sample) {
        split_col <- strsplit(vcf[[sample]], ':', fixed = T)
        sapply(seq_along(split_col),
            function(i) {
                j <- gt_ix[i]
                k <- gq_ix[i]
                if (is.na(j) | is.na(k)) {
                    gt <- '.'
                } else {
                    gt <- split_col[[i]][j]
                }
                if (gt == '.') {
                    c(NA, NA, NA)
                } else {
                    gq <- split_col[[i]][k]
                    ac0 <- str_count(gt, '0')
                    ac1 <- str_count(gt, '/') + 1 - ac0
                    gq <- as.numeric(gq)
                    c(ac0, ac1, gq)
                }
            })
    }) |> as.vector()
    ac_matrix2 <- array(ac_matrix, , dim = c(3, nrow(vcf), length(samples))) |>
        aperm(c(2, 3, 1))
    # Vars, Samples, REF/ALT
    dimnames(ac_matrix2) <- list(
        add_comma(vcf$POS),
        samples,
        c('Ref', 'AnyAlt', 'Qual'))
    ac_matrix2
}

vcf <- load_vcf('SMN1-AFR.vcf.gz')
ac_matrix <- get_ac_matrix(vcf)

pos2 <- with(vcf,
    setNames(sub('.*pos2=[^:]+:([0-9]+):.*', '\\1', INFO) |> add_comma(), add_comma(POS)))

TAG_PSV <- '70,951,946'
VARS <- c('70,923,922', '70,951,020', '70,952,074')
QUAL <- 10

sample_cn0 <- ac_matrix[TAG_PSV,,] |> as.data.frame() |>
    rownames_to_column('sample') |>
    filter(Qual >= QUAL) |>
    with(setNames(sprintf('%s,%s', Ref, AnyAlt), sample))

CNS <- c('2,1', '2,2', '4,0', '3,1')
samples <- names(sample_cn0)[sample_cn0 %in% CNS]
sample_cn <- sample_cn0[samples]

LIMIT <- 15
ac_matrix2 <- ac_matrix[VARS, samples,]
ac_long <- as.data.frame(ac_matrix2) |>
    rownames_to_column('var') |>
    pivot_longer(-var, names_to = 'sample_ty', values_to = 'ac') |>
    separate(sample_ty, c('sample', 'allele')) |>
    mutate(allele = recode_factor(allele, 'Ref' = 'ref', 'AnyAlt' = 'alt')) |>
    pivot_wider(names_from = 'allele', values_from = 'ac') |>
    filter(Qual >= QUAL) |>
    mutate(
        cn = factor(sample_cn[sample], levels = CNS),
        agcn = as.numeric(substr(cn, 1, 1)) + as.numeric(substr(cn, 3, 3)),
    ) |>
    group_by(sample) |>
    filter(all(ref + alt == agcn)) |>
    ungroup() |>
    arrange(var) |>
    mutate(
        var = sprintf('%s\n%s', var, pos2[var]),
        var = factor(var, levels = unique(var)),
        dosage = sprintf('%d/%d', ref, ref + alt),
        genotype = paste0(str_dup('A', times = ref), str_dup('a', times = alt)),
    )
panel_sizes <- count(ac_long, var, cn, name = 'panel_size')
bars <- count(ac_long, var, dosage, genotype, ref, cn, name = 'bar_size') |>
    left_join(panel_sizes, join_by(var, cn)) |>
    mutate(perc = 100 * bar_size / panel_size)

colors <- colorspace::sequential_hcl(7, 'Dark mint')[5:1] |>
    setNames(as.character(0:4))
ggplot(bars) +
    geom_bar(aes(dosage, perc, fill = as.character(ref)), stat = 'identity') +
    geom_text(aes(dosage, perc, label = bar_size),
        vjust = ifelse(bars$perc >= LIMIT, 1.3, -0.3),
        color = ifelse(bars$perc >= LIMIT, 'white', colors[5]),
        size = 2, family = 'Carlito') +
    facet_grid(var ~ cn, scales = 'free_x', space = 'free_x') +
    scale_x_discrete('Allele dosage') +
    scale_y_continuous('Percentage of samples', breaks = c(0, 50, 100),
        expand = expansion(add = c(1, 5))) +
    scale_fill_manual(values = colors) +
    ggtitle('*SMN1, SMN2* copy number') +
    theme_bw() +
    theme(
        text = element_text(family = 'Carlito'),
        plot.title = ggtext::element_markdown(size = 10, hjust = 0.5,
            margin = margin(b = 2)),
        panel.border = element_rect(color = 'gray80'),
        panel.grid = element_blank(),
        strip.background = element_rect(color = NA, fill = 'gray90'),
        strip.text.x = element_text(margin = margin(2, 2, 2, 2)),
        strip.text.y = element_text(size = 7, margin = margin(2, 2, 2, 2)),
        legend.position = 'none',
    )
ggsave('fig5.svg', width = 8, height = 5, scale = 0.75)
