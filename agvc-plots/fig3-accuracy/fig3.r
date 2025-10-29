#!/usr/bin/env Rscript
pdf(NULL)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggh4x))

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
    & (tool == 'parascopy' | (tool_qual == 'Q1' & qual == 'any') | (tool_qual == 'Q10' & qual == 'high')))

full_roc4 <- mutate(full_roc3,
    region = sub('^(non-psv|non-trivial_psv|trivial_psv)_', '\\1@', region)) |>
    separate(region, c('psv_type', 'difficulty'), sep = '@') |>
    mutate(
        psv_type = recode_factor(psv_type, 'trivial_psv' = 'Trivial PSVs',
                'non-trivial_psv' = 'Non-trivial PSVs', 'non-psv' = 'Non-PSVs'),
        difficulty = recode_factor(difficulty,
            'all' = 'All duplications', 'difficult' = 'High sequence similarity'),
        sample = recode_factor(sample, 'SimPoly' = 'Simulated'),
        tool = recode_factor(tool, 'freebayes' = 'Freebayes', 'gatk' = 'GATK', 'parascopy' = 'Parascopy'),
    )
roc <- group_by(full_roc4, sample, tool, qual, psv_type, difficulty, var_type) |>
    summarize_at(vars(tp_base, tp_call, fp, fn), sum) |>
    ungroup() |>
    filter(qual == 'high' & var_type == 'all') |>
    mutate(psv_type = recode_factor(psv_type, 'Non-trivial PSVs' = 'PSVs')) |>
    select(!c(qual, var_type)) |>
    mutate(
        precision = tp_call / (tp_call + fp),
        recall = tp_base / (tp_base + fn),
        f1 = 2 * precision * recall / (precision + recall))
roc_long <- select(roc, sample, tool, psv_type, difficulty, precision, recall, f1) |>
    pivot_longer(c(precision, recall, f1), names_to = 'metric') |>
    mutate(metric = recode_factor(metric,
        'precision' = 'Precision', 'recall' = 'Recall', 'f1' = 'F₁ score'))

roc_long2 <- filter(roc_long, !grepl('^F', metric) & psv_type != 'Trivial PSVs')
nvars <- filter(roc, tool == 'Parascopy' & psv_type != 'Trivial PSVs') |>
    group_by(sample) |>
    mutate(total = tp_base + fn, tool = factor(tool)) |>
    ungroup() |>
    mutate(metric = 'Recall')

FONT <- 'Source Sans 3'

draw_plot <- function(df, sample) {
    colors <- c('#5576A0', '#EF7747', '#FFDD87')
    w <- 0.9; u <- w / 2
    s <- sample
    ggplot(filter(df, sample == s)) +
        geom_bar(aes(as.numeric(tool), value, fill = tool),
            stat = 'identity', position = position_dodge(),
            color = 'gray10', linewidth = 0.1, width = w) +
        geom_text(aes(3.95, 0.02,
            label = sprintf('%s variants', scales::label_comma(precision = 0)(total))),
            data = filter(nvars, sample == s),
            family = FONT, angle = -90, vjust = 1, hjust = 1, size = 3) +
        facet_nested(difficulty ~ psv_type + metric, scales = 'free', space = 'free') +
        scale_fill_manual('Variant caller', values = colors) +
        scale_x_continuous(NULL,
            breaks = 1:3, labels = levels(roc_long2$tool),
            expand = expansion(add = 0.02), limits = c(1 - u - 0.1, NA)) +
        scale_y_continuous('Precision / Recall',
            breaks = seq(0, 1, 0.2), expand = expansion(add = 0.005), limits = c(0, 1),
            guide = guide_axis(minor.ticks = TRUE)) +
        ggtitle(NULL, subtitle = sample) +
        theme_bw() +
        theme(
            text = element_text(family = FONT),
            axis.title.y = element_text(size = 10, margin = margin(r = 2)),
            axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(vjust = c(0, 0.5, 0.5, 0.5, 0.5, 1)),
            axis.minor.ticks.length = rel(0.4),
            panel.background = element_rect(fill = NA),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.border = element_rect(color = NA),
            strip.background = element_rect(color = 'gray90', fill = 'gray90'),
            strip.text = element_text(size = 8, margin = margin(t = 1, b = 1, l = 2, r = 2), color = 'black'),
            legend.position = 'none',
            legend.key.size = unit(0.8, 'lines'),
            legend.margin = margin(b = -4, t = 0),
            legend.title = element_text(size = 9.5, face = 'bold', margin = margin(l = 5, r = 5)),
            legend.text = element_text(size = 9, margin = margin(l = 2, r = 3)),
            strip.placement = 'outside',
            panel.spacing.x = unit(0.5, 'lines'),
            panel.spacing.y = unit(0.3, 'lines'),
            plot.subtitle = element_text(size = 10)
        )
}

cowplot::plot_grid(
    draw_plot(roc_long2, 'HG007'),
    draw_plot(roc_long2, 'Simulated'),
    labels = letters, label_fontfamily = FONT,
    label_x = 0.02, label_y = 1.004
    ) |>
    suppressWarnings()
ggsave('barplot.svg', width = 10, height = 6, scale = 0.7)

# filter(roc_long2, sample == 'Simulated', difficulty != 'All duplications') |>
#     mutate(value = round(value, 2))
# 
# filter(roc_long2, sample == 'Simulated' & tool == 'GATK') |>
#     mutate(value = round(value, 2))
# 
# filter(roc_long2, sample %in% c('Simulated', 'HG007') &
#         difficulty == 'All duplications' & tool == 'Freebayes') |>
#     mutate(value = round(value, 2)) |>
#     arrange(psv_type, metric, desc(sample))
