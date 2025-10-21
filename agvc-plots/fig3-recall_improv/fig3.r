#!/usr/bin/env Rscript
pdf(NULL)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))

add_metrics <- function(df) {
    mutate(df,
        n_vars = tp_base + fn,
        precision = replace_na(tp_call / (tp_call + fp), 0),
        recall = replace_na(tp_base / (tp_base + fn), 0),
        f1 = replace_na(2 * precision * recall / (precision + recall), 0),
    )
}

round_to <- function(x, precision) {
    round(x / precision) * precision
}

dupl_genes <- readLines('dupl_genes.txt')
filenames <- Sys.glob('evals/*.csv.gz')

full_roc <- lapply(filenames,
    function(filename) read.csv(filename, sep = '\t', comment = '#') |>
        mutate(filename = basename(filename))) %>%
    do.call(rbind, .)

# Splitting filenames
full_roc2 <- mutate(full_roc,
    filename = sub('.csv.gz$', '', filename) %>% sub('parascopy', 'parascopy.NA', .)) |>
    separate(filename, c('sample', 'tool', 'tool_qual', 'cn'))

# Trust Parascopy qual column, for other tools, use quality from the `QX` part of the filename.
full_roc3 <- filter(full_roc2,
    region != '*'
    & (tool == 'parascopy' | (tool_qual == 'Q1' & qual == 'any') | (tool_qual == 'Q10' & qual == 'high'))) |>
    select(-tool_qual)

# Split region into gene & PSV type
full_roc4 <- full_roc3 |>
    separate(region, c('gene', 'psv_type'), sep = '@') |>
    mutate(
        psv_type = recode_factor(psv_type,
            'nonpsv' = 'Non-PSV', 'nontrivialpsv' = 'PSV', 'trivialpsv' = 'Trivial PSV'),
        #tool = recode_factor(tool, 'freebayes' = 'Freebayes', 'gatk' = 'GATK', 'parascopy' = 'Parascopy'),
        )

# Take all variants together, high quality, duplicated genes
roc <- filter(full_roc4,
        qual == 'high'
        & var_type == 'all'
        & sample != 'SimFixed'
        & gene %in% dupl_genes) |>
    select(-c(qual, var_type, cn)) |>
    mutate(sample_type = ifelse(startsWith(sample, 'Sim'), 'Simulated', 'GIAB'))

# Sum across different CNs and PSV types for the same gene
roc2 <- filter(roc, psv_type != 'Trivial PSV') |>
    group_by(gene, tool, sample, sample_type) |>
    summarize_at(c('tp_base', 'tp_call', 'fp', 'fn'), sum) |>
    ungroup()

# Take mean across all samples
roc3 <- group_by(roc2, gene, tool, sample_type) |>
    summarize_at(c('tp_base', 'tp_call', 'fp', 'fn'), mean) |>
    ungroup() |>
    add_metrics()

# Subset with genes with enough variants
MIN_VARS <- 10
keep_genes <- filter(roc3, tool == 'parascopy') |>
    mutate(keep = n_vars >= MIN_VARS) |>
    select(gene, sample_type, keep)
roc4 <- left_join(roc3, keep_genes, join_by(gene, sample_type)) |>
    filter(keep) |> select(-keep)

roc_long <- select(roc4, -c(tp_base, tp_call, fp, fn)) |>
    pivot_longer(cols = c(precision, recall, f1), names_to = 'metric') |>
    mutate(
        tool = recode_factor(tool, 'freebayes' = 'Freebayes', 'gatk' = 'GATK', 'parascopy' = 'Parascopy'),
        metric = recode_factor(metric, 'precision' = 'Precision', 'recall' = 'Recall', 'f1' = 'F₁ score'))
roc_lwc <- mutate(roc_long, value = round_to(value, 0.02)) |>
    select(-n_vars) |>
    pivot_wider(names_from = 'tool', values_from = 'value')

palette <- 'Rocket'
legend_col <- colorspace::sequential_hcl(11, palette = palette)[3]

ggplot(filter(roc_lwc, GATK < 1 | Parascopy < 1)) +
    geom_abline(linewidth = 0.3, color = 'gray30') +
    annotate('rect', xmin = 1, ymin = 1, xmax = Inf, ymax = Inf, fill = 'white') +
    geom_count(aes(GATK, Parascopy, color = - GATK + Parascopy)) +
    annotate('point', x = 1, y = 1, size = 2, color = 'gray25') +
    facet_grid(metric ~ sample_type, scales = 'free_x', space = 'free_x') +
    scale_x_continuous('GATK accuracy',
        breaks = seq(0, 1, 0.2), expand = expansion(add = 0.02)) +
    scale_y_continuous('Parascopy accuracy',
        breaks = seq(0, 1, 0.2), limits = c(0.6, 1), expand = expansion(add = 0.027)) +
    colorspace::scale_color_continuous_sequential('Rocket', begin = 0.2) +
    scale_size('Number of genes', breaks = c(1, 10, 20), range = c(0.3, 2.5)) +
    guides(color = 'none',
        size = guide_legend(override.aes = list(color = legend_col, alpha = 0.8))) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(margin = margin(2, 2, 2, 2), face = 'bold'),
        panel.border = element_rect(color = 'gray80', linewidth = 0.5),
        panel.spacing.x = unit(0.7, 'lines'),
        legend.position = 'bottom',
        legend.margin = margin(t = -9),
        legend.title = element_text(size = 10),
        legend.text = element_text(margin = margin(l = -2, r = 2)),
        plot.margin = margin(2, 2, 2, 2)
    )
ggsave('improv_recall.svg', width = 8, height = 6, scale = 0.8, bg = 'white')

