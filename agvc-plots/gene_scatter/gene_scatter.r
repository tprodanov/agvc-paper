#!/usr/bin/env Rscript
pdf(NULL)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))

f_score <- function(precision, recall, beta = 1) {
    ((1 + beta * beta) * precision * recall / (beta * beta * precision + recall)) |>
        replace_na(0)
}

add_metrics <- function(df) {
    mutate(df,
        n_vars = tp + fn,
        precision = replace_na(tp / (tp + fp), 0),
        recall = replace_na(tp / (tp + fn), 0),
        # f0.25 = f_score(precision, recall, 0.25),
        # f0.5 = f_score(precision, recall, 0.5),
        f1 = f_score(precision, recall),
    )
}

round_to <- function(x, precision) {
    round(x / precision) * precision
}

#dupl_genes <- readLines('dupl_genes.txt')
filenames <- Sys.glob('evals/*.csv.gz')

full_roc <- lapply(filenames,
    function(filename) read.csv(filename, sep = '\t', comment = '#') |>
        mutate(filename = basename(filename))) %>%
    do.call(rbind, .)

# Splitting filenames
full_roc2 <- full_roc |>
    filter(var_type == 'all') |> select(-var_type) |>
    mutate(filename = sub('.csv.gz$', '', filename) %>%
            sub('parascopy', 'parascopy.NA', .)) |>
    separate(filename, c('sample', 'tool', 'tool_qual', 'cn'))

# Trust Parascopy qual column, for other tools, use quality from the `QX` part of the filename.
full_roc3 <- filter(full_roc2,
        region != '*',
        (tool == 'parascopy' & qual == 'high') | (tool_qual == 'Q10' & qual == 'any')) |>
    select(-c(qual, tool_qual))

# Split region into gene & PSV type
full_roc4 <- full_roc3 |>
    separate(region, c('gene', 'psv_type'), sep = '@') |>
    mutate(
        psv_type = recode_factor(psv_type,
            'nonpsv' = 'Non-PSV', 'nontrivialpsv' = 'PSV', 'trivialpsv' = 'Trivial PSV'),
        #tool = recode_factor(tool, 'freebayes' = 'Freebayes', 'gatk' = 'GATK', 'parascopy' = 'Parascopy'),
        )

# Take all variants together, high quality, duplicated genes
roc <- full_roc4 |>
    select(-cn) |>
    mutate(sample_type = ifelse(startsWith(sample, 'Sim'), 'Simulated', 'GIAB') |>
            factor(levels = c('Simulated', 'GIAB')))

# Sum across different CNs and PSV types for the same gene
roc2 <- filter(roc, psv_type != 'Trivial PSV') |>
    group_by(gene, tool, sample, sample_type) |>
    summarize_at(c('tp', 'fp', 'fn'), sum) |>
    ungroup()

# Take mean across all samples
roc3 <- group_by(roc2, gene, tool, sample_type) |>
    summarize_at(c('tp', 'fp', 'fn'), mean) |>
    ungroup() |>
    add_metrics()

# Subset with genes with enough variants
MIN_VARS <- 10
keep_genes <- filter(roc3, tool == 'parascopy') |>
    mutate(keep = n_vars >= MIN_VARS) |>
    select(gene, sample_type, keep)
filter(keep_genes, keep) %>%
    aggregate(gene ~ sample_type, ., function(x) length(unique(x)))

roc4 <- left_join(roc3, keep_genes, join_by(gene, sample_type)) |>
    filter(keep) |> select(-keep)

roc_long <- select(roc4, -c(tp, fp, fn)) |>
    pivot_longer(cols = c(precision, recall, f1), names_to = 'metric') |>
    mutate(
        tool = recode_factor(tool,
            'freebayes' = 'Freebayes', 'gatk' = 'GATK',
            'deepvariant' = 'DeepVariant', 'parascopy' = 'Parascopy'),
        metric = recode_factor(metric,
            'precision' = 'Precision',
            'recall' = 'Recall',
            'f1' = 'F₁ score'))
roc_lwc <- mutate(roc_long, value_round = round_to(value, 0.01)) |>
    select(-n_vars) |>
    pivot_wider(names_from = 'tool', values_from = c('value', 'value_round')) |>
    rename_with(function(name) sub('value_', '', name) %>% sub('round_', 'rnd.', .), starts_with('value'))

palette <- 'Rocket'
legend_col <- colorspace::sequential_hcl(11, palette = palette)[3]

gene_arrows <- filter(roc_lwc, gene %in% c('CFC1', 'NEB', 'SMN1')) |>
    mutate(x = rnd.GATK, y = rnd.Parascopy) |>
    select(gene, sample_type, metric, x, y) |>
    arrange(x) |>
    group_by(sample_type, metric) |>
    mutate(
        x_rank = rank(x, ties.method = 'first') - (n() + 1) / 2,
        dist_left = c(0.5, diff(x)),
        dist_right = c(diff(x), 0.5),
        ) |>
    ungroup() |>
    mutate(
        # x2 = pmin(x + 0.1 * x_rank + 0.02, x + 0.13),
        x2 = (x + 0.02) - (1 - dist_right) * 0.2 + (1 - dist_left) * 0.2,
        y2 = case_when(metric == 'Recall' ~ 0.73, metric == 'Precision' ~ 0.88, T ~ 0.85),
        curvature = 2 * (x - x2))
roc_lwc2 <- mutate(roc_lwc, n_genes = 1) |> # + str_count(gene, ',')) |>
    group_by(sample_type, metric, rnd.GATK, rnd.Parascopy) |>
    summarize(n_genes = sum(n_genes), .groups = 'keep') |>
    ungroup()

filter(roc_lwc2, rnd.GATK < 1 | rnd.Parascopy < 1) |>
    filter(rnd.GATK == max(rnd.GATK))

FONT <- 'Roboto'
ggplot(filter(roc_lwc2, rnd.GATK < 1 | rnd.Parascopy < 1 | sample_type == 'GIAB')) +
    geom_abline(linewidth = 0.3, color = 'gray30') +
    annotate('rect', xmin = 1, ymin = 1, xmax = Inf, ymax = Inf, fill = 'white') +
    lapply(split(gene_arrows, 1:nrow(gene_arrows)), function(dat) {
        geom_curve(
            aes(x = x, y = y, xend = x2, yend = y2),
            data = dat,
            color = legend_col,
            curvature = dat['curvature'],
            ncp = 3, angle = 90, alpha = 0.5,
            )
    }) +
    geom_label(aes(x2, y2, label = gene), data = gene_arrows,
        color = legend_col, family = FONT, fontface = 'italic', vjust = 1., size = 3,
        linewidth = 0., label.padding = unit(0.01, 'lines')
        ) +
    geom_point(aes(rnd.GATK, rnd.Parascopy,
        size = n_genes, color = rnd.Parascopy - rnd.GATK)) + #, color = '#AF0065') +
    geom_point(aes(x, y), size = 2.5, color = 'gray25',
        data = data.frame(x = 1, y = 1,
            sample_type = factor('Simulated', levels = levels(roc_lwc$sample_type)),
            metric = unique(roc_lwc2$metric))) +
    # annotate('point', x = 1, y = 1, size = 2.5, color = 'gray25') +
    facet_grid(metric ~ sample_type, scales = 'free_x', space = 'free_x') +
    scale_x_continuous('GATK accuracy',
        breaks = seq(0, 1, 0.2), expand = expansion(add = 0.02)) +
    scale_y_continuous('Parascopy accuracy',
        breaks = seq(0, 1, 0.2), limits = c(NA, 1), expand = expansion(add = 0.027)) +
    colorspace::scale_color_continuous_sequential('Rocket', begin = 0.2) +
    scale_size('Number of genes', breaks = c(1, 10, 20), range = c(0.2, 2)) +
    guides(color = 'none',
        size = guide_legend(override.aes = list(color = legend_col, alpha = 0.8))) +
    theme_bw() +
    theme(
        text = element_text(family = FONT),
        panel.grid = element_blank(),
        strip.background = element_rect(color = NA, fill = NA),
        strip.text.x = element_text(margin = margin(2, 2, 2, 2), face = 'bold'),
        strip.text.y = element_text(margin = margin(l = 0.5, 2, 2, 2), face = 'bold'),
        panel.border = element_rect(color = 'gray80', linewidth = 0.5),
        panel.spacing.x = unit(0.7, 'lines'),
        legend.position = 'bottom',
        legend.margin = margin(t = -9),
        legend.title = element_text(size = 10),
        legend.text = element_text(margin = margin(l = -2, r = 2)),
        plot.margin = margin(2, 2, 2, 2)
    )
ggsave('gene_scatter.svg', width = 8, height = 6, scale = 0.8, bg = 'white')
ggsave('gene_scatter.png', width = 8, height = 6, scale = 0.8, bg = 'white', dpi = 500)

group_by(roc_lwc, sample_type, metric) |>
    summarize(n_under = sum(GATK > Parascopy)) |>
    ungroup()
group_by(roc_lwc, sample_type, metric) |>
    summarize(n_less = sum(GATK > Parascopy), n_greater = sum(Parascopy > GATK)) |>
    ungroup()

for (sample_type in c('Simulated', 'GIAB')) {
    total_genes <- roc_lwc[roc_lwc$sample_type == sample_type,] |>
        with(length(unique(gene)))
        #with(sum(1 + str_count(unique(gene), ',')))
    cat(sprintf('%s: %3d genes\n', sample_type, total_genes))
    for (thresh in c(0, 10, 25)) {
        for (metric in levels(roc_lwc$metric)) {
            cat(sprintf('    %-8s, at least %d p.p. better\n', metric, thresh))
            for (tool in c('GATK', 'Freebayes', 'DeepVariant')) {
                n_genes <- roc_lwc[
                    roc_lwc$sample_type == sample_type
                    & roc_lwc$metric == metric
                    & roc_lwc$Parascopy > roc_lwc[tool] + 0.01 * thresh
                    , ] |>
                    filter(!is.na(gene)) |>
                    with(length(unique(gene)))
                    #with(sum(1 + str_count(unique(gene), ',')))
                cat(sprintf('        than %12s: %3d genes\n', tool, n_genes))
            }
        }
    }
}

filter(roc_lwc2, rnd.Parascopy == 1 & rnd.GATK == 1)

filter(roc2, tool == 'parascopy', gene == 'ANKUB1')
filter(roc4, gene == 'NEB')

highly_similar <- read.csv('overlapped_genes.s0.99_m0.csv', sep = '\t', header = F) |>
    filter(V2 >= 3000) |> with(V1)
genes_by_improv <- filter(roc_lwc, sample_type == 'Simulated' & metric == 'Recall') |>
    mutate(d = Parascopy - GATK) |>
    select(gene, d) |>
    separate_rows(gene, sep = ',') |>
    mutate(simil = gene %in% highly_similar)
count(genes_by_improv, d >= 0.25, simil)
fisher.test(matrix(c(341, 119, 42, 175), ncol = 2))
