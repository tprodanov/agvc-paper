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
toolx <- c('Freebayes' = 1, 'GATK' = 1.4, 'Parascopy' = 2.15)

colors <- c('#5576A0', '#EF7747', '#FFDD87')
ggplot(roc_long) +
    # geom_hline(yintercept = c(0, 1), linewidth = 0.3, linetype = '15', color = 'gray30') +
    annotate('rect', xmin = -Inf, xmax = Inf, ymin = 0 - 0.005, ymax = 1 + 0.005, fill = 'gray97') +
    geom_violin(aes(toolx[tool], value, group = tool, fill = tool),
        kernel = 'g', bw = 0.015, width = 1, linewidth = 0.3) +
    facet_grid(metric ~ sample_type) +
    scale_x_continuous('Variant caller', breaks = toolx, expand = c(0, 0)) +
    scale_y_continuous('Value', breaks = seq(0, 1, 0.2)) +
    scale_fill_manual(values = colors) +
    coord_cartesian(xlim = c(0.8, 2.65)) +
    theme_bw() +
    theme(
        text = element_text(family = 'Source Sans 3'),
        panel.grid = element_blank(),
        legend.position = 'none',
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(margin = margin(2, 2, 2, 2), face = 'bold'),
        panel.border = element_blank(),
        panel.spacing.x = unit(1, 'lines'),
        # axis.line = element_line(linewidth = 0.25, color = 'gray30'),
    )
ggsave('~/Downloads/1.png', width = 8, height = 6, dpi = 500, scale = 0.8)
ggsave('improv_recall2.png', width = 8, height = 6, dpi = 500, scale = 0.8)

roc_lwc <- mutate(roc_long, value = round_to(value, 0.02)) |>
    select(-n_vars) |>
    pivot_wider(names_from = 'tool', values_from = 'value')

palette <- 'Rocket'
legend_col <- colorspace::sequential_hcl(11, palette = palette)[3]

ggplot(roc_lwc) +
    geom_abline(linewidth = 0.3, color = 'gray30') +
    annotate('rect', xmin = 1, ymin = 1, xmax = Inf, ymax = Inf, fill = 'white') +
    geom_count(aes(GATK, Parascopy,
        color = - GATK + Parascopy,
        )) +
    facet_grid(metric ~ sample_type, scales = 'free_x', space = 'free_x') +
    scale_x_continuous(breaks = seq(0, 1, 0.2), expand = expansion(add = c(0.02, 0.04))) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0.6, 1), expand = expansion(add = 0.05)) +
    colorspace::scale_color_continuous_sequential('Rocket', begin = 0.2) +
    scale_size_area('Number of genes', breaks = c(1, 10, 50, 150)) +
    guides(color = 'none',
        size = guide_legend(override.aes = list(color = legend_col, alpha = 0.8))) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(margin = margin(2, 2, 2, 2), face = 'bold'),
        panel.border = element_rect(color = 'gray80', linewidth = 0.5),
        # panel.border = element_blank(),
        panel.spacing.x = unit(0.7, 'lines'),
        legend.position = 'bottom',
        legend.margin = margin(t = -7),
        legend.text = element_text(margin = margin(l = 0, r = 4)),
    )
ggsave('improv_recall3.png', width = 8, height = 5, dpi = 500, scale = 0.8)


roc_wide <- select(roc3, gene, tool, sample_type, n_vars, precision, recall, f1) |>
    pivot_wider(names_from = 'tool',
        values_from = c('n_vars', 'precision', 'recall', 'f1'), names_sep = '.') |>
    filter(!is.na(n_vars.parascopy) & !is.na(n_vars.freebayes))


ggplot(roc3)

y_range <- with(roc_wide2, range(recall.parascopy - recall.freebayes))
colors <- colorspace::sequential_hcl(7, palette = 'Emrld')[c(2, 4)]




ggplot(roc_wide2) + # Split by sample type and variant type (+ by tool?)
    geom_abline() +
    coord_fixed() +
    geom_count(
        aes(round_to(recall.freebayes, 0.02), round_to(recall.parascopy, 0.02)))

arrange(to_plot, -recall.freebayes + recall.parascopy)
filter(to_plot, gene == 'TNRC18') |> as.data.frame()
filter(roc3, gene == 'TNRC18' & sample_type == 'Simulated' & psv_type != 'Trivial PSV' &
        tool %in% c('freebayes', 'parascopy'))

ggplot(roc_giab) +
    geom_violin(aes(psv_type, recall.improv))

ggplot(roc_giab) +
    geom_histogram(aes(recall.improv, fill = psv_type, color = psv_type,
        y = ifelse(fill == 'Non-PSV', 1, -1) * after_stat(count)),
        position = 'identity', binwidth = 0.01, boundary = 0, linewidth = 0.1) +
    scale_x_continuous('Improvement in recall', breaks = seq(0, 0.8, 0.2),
        expand = expansion(mult = 0.01)) +
    scale_y_continuous('Number of genes', labels = abs, expand = expansion(mult = 0.02)) +
    # coord_cartesian(ylim = c(-20, 20)) +
    scale_fill_manual(NULL, values = colors) +
    scale_color_manual(NULL, values = colors) +
    theme_bw() +
    theme(
        text = element_text(family = 'Carlito'),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'inside',
        legend.justification = c('right', 'top'),
        legend.position.inside = c(0.98, 0.9),
    )
ggsave('improv_recall.svg', width = 10, height = 5, scale = 0.65)

# Calculate stats

roc6 <- filter(roc3, sample_type == 'Simulated' & psv_type != 'Trivial PSV') |>
    group_by(gene, tool) |>
    summarize_at(c('tp_base', 'tp_call', 'fp', 'fn'), sum, na.rm = T) |>
    ungroup() |>
    add_metrics() |>
    pivot_wider(names_from = 'tool',
        values_from = c('tp_base', 'tp_call', 'fp', 'fn', 'n_vars', 'precision', 'recall', 'f1'), names_sep = '.') |>
    filter(!is.na(n_vars.parascopy)) |>
    mutate(recall.std = pmax(recall.freebayes, recall.gatk, na.rm = T)) |> as.data.frame()
MIN_VARS <- 10
with(roc6, sum(n_vars.parascopy >= MIN_VARS))
with(roc6, sum(n_vars.parascopy >= MIN_VARS & recall.parascopy > replace_na(recall.std, 0)))
with(roc6, sum(n_vars.parascopy >= MIN_VARS & recall.parascopy > replace_na(recall.gatk, 0)))
with(roc6, sum(n_vars.parascopy >= MIN_VARS & recall.parascopy > replace_na(recall.freebayes, 0)))
with(roc6, sum(n_vars.parascopy >= MIN_VARS & recall.parascopy >= replace_na(recall.std, 0) + 0.25))
with(roc6, sum(n_vars.parascopy >= MIN_VARS & recall.parascopy >= replace_na(recall.gatk, 0) + 0.25))
with(roc6, sum(n_vars.parascopy >= MIN_VARS & recall.parascopy >= replace_na(recall.freebayes, 0) + 0.25))
with(roc6, sum(n_vars.parascopy >= MIN_VARS & recall.parascopy >= replace_na(recall.std, 0) + 0.1))
with(roc6, sum(n_vars.parascopy >= MIN_VARS & recall.parascopy >= replace_na(recall.gatk, 0) + 0.1))
with(roc6, sum(n_vars.parascopy >= MIN_VARS & recall.parascopy >= replace_na(recall.freebayes, 0) + 0.1))

filter(roc6, gene == 'CFC1')
filter(roc6, grepl('SMN[12]', gene))
filter(roc4, gene == 'SMN1' & psv_type != 'Trivial PSV') |> as.data.frame()
filter(roc4, gene == 'SMN1') |> as.data.frame()

filter(roc6, gene == 'NEB')
