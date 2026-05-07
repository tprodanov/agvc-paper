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
    filename = sub('.csv.gz$', '', filename) %>%
        sub('parascopy', 'parascopy.NA', .)) |>
    separate(filename, c('sample', 'tool', 'tool_qual', 'cn'))

# Trust Parascopy qual column, for other tools, use quality from the `QX` part of the filename.
full_roc3 <- filter(full_roc2,
    region != '*',
    tool == 'parascopy'
        | (tool_qual == 'Q1' & qual == 'any')
        | (tool_qual == 'Q10' & qual == 'high'))

full_roc4 <- full_roc3 |>
    separate(region, c('psv_type', 'difficulty'), sep = '_') |>
    mutate(
        psv_type = recode_factor(psv_type, 'trivialpsv' = 'Trivial PSVs',
                'nontrivialpsv' = 'Non-trivial PSVs', 'nonpsv' = 'Non-PSVs'),
        difficulty = recode_factor(difficulty,
            'all' = 'All duplications', 'difficult' = 'High sequence similarity'),
        sample = recode_factor(sample, 'SimPoly' = 'Simulated'),
        tool = recode_factor(tool,
            'freebayes' = 'Freebayes', 'gatk' = 'GATK',
            'deepvariant' = 'DeepVariant', 'parascopy' = 'Parascopy'),
    )
roc <- group_by(full_roc4, sample, tool, qual, psv_type, difficulty, var_type) |>
    summarize_at(vars(tp, fp, fn), sum) |>
    ungroup() |>
    filter(qual == 'high' & var_type == 'all') |>
    mutate(psv_type = recode_factor(psv_type, 'Non-trivial PSVs' = 'PSVs')) |>
    select(!c(qual, var_type)) |>
    mutate(
        precision = replace_na(tp / (tp + fp), 0),
        recall = replace_na(tp / (tp + fn), 0),
        f1 = 2 * precision * recall / (precision + recall))
roc_long <- select(roc, sample, tool, psv_type, difficulty, precision, recall, f1) |>
    pivot_longer(c(precision, recall, f1), names_to = 'metric') |>
    mutate(metric = recode_factor(metric,
        'precision' = 'Precision', 'recall' = 'Recall', 'f1' = 'F₁ score'))

roc_long2 <- filter(roc_long, !grepl('^F', metric) & psv_type != 'Trivial PSVs')
nvars <- filter(roc, tool == 'Parascopy' & psv_type != 'Trivial PSVs') |>
    group_by(sample) |>
    mutate(total = tp + fn, tool = factor(tool)) |>
    ungroup() |>
    mutate(metric = 'Recall')

FONT <- 'Roboto'
colors <- c('Freebayes' = '#5576A0', 'GATK' = '#EF7747',
    'DeepVariant' = '#823329', 'Parascopy' = '#FFDC5E')

draw_plot <- function(df, sample) {
    w <- 0.9; u <- w / 2
    s <- sample
    df <- filter(df, sample == s) |> mutate(tool = droplevels(tool))
    nvars <- filter(nvars, sample == s) |> mutate(tool = droplevels(tool))
    tools <- levels(df$tool)
    n_tools <- length(tools)
    
    ggplot(df) +
        geom_bar(aes(as.numeric(tool), value, fill = tool),
            stat = 'identity', position = position_dodge(),
            color = 'gray10', linewidth = 0.1, width = w) +
        geom_text(aes(n_tools + 1.02, 0.02,
            label = sprintf('%s variants', scales::label_comma(precision = 0)(total))),
            data = nvars,
            family = FONT, angle = -90, vjust = 1, hjust = 1, size = 2.7) +
        facet_nested(difficulty ~ psv_type + metric, scales = 'free', space = 'free') +
        scale_fill_manual('Variant caller', values = colors) +
        scale_x_continuous(NULL,
            breaks = 1:n_tools, labels = tools,
            expand = expansion(add = 0.02),
            limits = c(1 - u - 0.1, NA)) +
        scale_y_continuous('Precision / Recall',
            breaks = seq(0, 1, 0.2), expand = expansion(add = 0.005), limits = c(0, 1),
            guide = guide_axis(minor.ticks = TRUE)) +
        ggtitle(NULL, subtitle = sample) +
        theme_bw() +
        theme(
            text = element_text(family = FONT),
            axis.title.y = element_text(size = 10, margin = margin(r = 2)),
            axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(vjust = c(0.1, 0.5, 0.5, 0.5, 0.5, 0.9)),
            axis.minor.ticks.length = rel(0.4),
            panel.background = element_rect(fill = NA),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.border = element_rect(color = NA),
            strip.background = element_rect(color = 'gray90', fill = 'gray90'),
            strip.text = element_text(size = 8, margin = margin(t = 1, b = 1, l = 2, r = 2),
                color = 'black'),
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
    label_x = 0.02, label_y = 1.01,
    rel_widths = c(0.55, 0.45),
    align = 'h'
    ) |>
    suppressWarnings()
ggsave('barplot.svg', width = 10, height = 6, scale = 0.7)

filter(roc_long2, sample == 'HG007' & tool == 'Parascopy')
filter(roc_long2, sample == 'HG007' & tool == 'GATK')
filter(roc_long2, sample == 'Simulated' & tool == 'GATK')
filter(roc_long2, sample %in% c('Simulated', 'HG007')
    & difficulty == 'All duplications' & psv_type == 'Non-PSVs' & metric == 'Recall')

filter(roc_long2, sample %in% c('HG007', 'Simulated')) |>
    pivot_wider(names_from = sample, values_from = value) |>
    mutate(diff = HG007 - Simulated, within_3p = abs(diff) <= 0.03) |>
    filter(tool != 'DeepVariant') |>
    filter(difficulty != 'High sequence similarity') |>
    print(n = 32)
# 
# # Accuracy in predicting non-trivial PSVs
# acc_nontriv <- roc |>
#     mutate(truecall = 0.5 * (tp_base + tp_call), falsecall = 0.5 * (fp + fn)) |>
#     select(-c(precision, recall, f1, tp_base, tp_call, fp, fn)) |>
#     filter(psv_type != 'Non-PSVs') |>
#     mutate(psv_type = ifelse(psv_type == 'PSVs', 'nontriv', 'trivial')) |>
#     pivot_wider(names_from = psv_type, values_from = c(truecall, falsecall)) |>
#     mutate(tp = truecall_nontriv, fp = falsecall_trivial, fn = falsecall_nontriv) |>
#     mutate(
#         precision = tp / (tp + fp),
#         recall = tp / (tp + fn),
#         #f1 = 2 * precision * recall / (precision + recall)
#     )
# acc_nontriv2 <- acc_nontriv |>
#     select(sample, tool, difficulty, precision, recall) |>
#     pivot_longer(c(precision, recall), names_to = 'metric') |>
#     mutate(
#         metric = recode_factor(metric,
#             'precision' = 'Precision', 'recall' = 'Recall') #, 'f1' = 'F₁ score')
#     )
# 
# draw_nontriv_acc <- function(df, sample) {
#     w <- 0.9; u <- w / 2
#     s <- sample
#     df <- filter(df, sample == s) |> mutate(tool = droplevels(tool))
#     
#     ggplot(df) +
#         geom_bar(aes(tool, value, fill = tool),
#             stat = 'identity', position = position_dodge(),
#             color = 'gray10', linewidth = 0.1, width = w) +
#         facet_grid(difficulty ~ metric, scales = 'free', space = 'free') +
#         scale_fill_manual('Variant caller', values = colors) +
#         scale_x_discrete(NULL, expand = expansion(add = u + 0.1)) +
#         scale_y_continuous('Precision / Recall',
#             breaks = seq(0, 1, 0.2), expand = expansion(add = 0.005), limits = c(0, 1),
#             guide = guide_axis(minor.ticks = TRUE)) +
#         ggtitle(NULL, subtitle = sprintf('Non-trivial PSV detection, %s', s) %>%
#             sub('Simulated', 'simulated data', .)) +
#         theme_bw() +
#         theme(
#             text = element_text(family = FONT),
#             axis.title.y = element_text(size = 10, margin = margin(r = 2)),
#             axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
#             axis.text.y = element_text(vjust = c(0.1, 0.5, 0.5, 0.5, 0.5, 0.9)),
#             axis.minor.ticks.length = rel(0.4),
#             panel.background = element_rect(fill = NA),
#             panel.grid.minor = element_blank(),
#             panel.grid.major.x = element_blank(),
#             panel.grid.major.y = element_blank(),
#             panel.border = element_rect(color = NA),
#             strip.background = element_rect(color = 'gray90', fill = 'gray90'),
#             strip.text = element_text(size = 8, margin = margin(t = 1, b = 1, l = 2, r = 2),
#                 color = 'black'),
#             legend.position = 'none',
#             legend.key.size = unit(0.8, 'lines'),
#             legend.margin = margin(b = -4, t = 0),
#             legend.title = element_text(size = 9.5, face = 'bold', margin = margin(l = 5, r = 5)),
#             legend.text = element_text(size = 9, margin = margin(l = 2, r = 3)),
#             strip.placement = 'outside',
#             panel.spacing.x = unit(0.3, 'lines'),
#             panel.spacing.y = unit(0.3, 'lines'),
#             plot.subtitle = element_text(size = 10)
#         )
# }
# 
# library(patchwork)
# 
# # layout <- c(
# #     area(t = 1, b = 1, l = 1, r = 2),
# #     area(t = 1, b = 1, l = 3, r = 4),
# #     area(t = 2, b = 2, l = 2, r = 2),
# #     area(t = 2, b = 2, l = 3, r = 3)
# # )
# 
# draw_plot(roc_long2, 'HG007') +
#     draw_plot(roc_long2, 'Simulated') +
#     draw_nontriv_acc(acc_nontriv2, 'HG007') +
#     draw_nontriv_acc(acc_nontriv2, 'Simulated') +
#     plot_layout(design = 'aabb\n#cd#', widths = c(0.35, 0.9, 0.7, 0.35)) +
#     plot_annotation(tag_levels = 'a') &
#     theme(
#         plot.tag = element_text(family = FONT, face = 'bold', size = 16),
#         plot.tag.position = c(0.05, 0.99),
#         plot.margin = margin(l = 3, t = 3, r = 3, b = 3))
# ggsave('barplot.svg', width = 10, height = 10.5, scale = 0.7)
# 
# filter(roc_long, psv_type == 'Trivial PSVs', sample == 'Simulated') |> print(n = 30)
# filter(roc_long, psv_type == 'PSVs', sample_type == 'GIAB') |> print(n = 30)
