#!/usr/bin/Rscript

pdf(NULL)
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(cowplot))
cowplot::set_null_device('agg')

load_gts <- function(filename) {
    read.csv(filename, sep = ' ') |>
        separate('ad', sep = ',', into = c('ad0', 'ad1', 'ad2'),
            convert = T, fill = 'right') |>
        mutate(
            qual = suppressWarnings(as.numeric(qual)),
            ac0 = str_count(gt, '0'),
            cn = str_count(gt, '/') + 1,
            ad1 = ad1 + replace_na(ad2, 0),
            af0 = ad0 / (ad0 + ad1),
        )
}

gts <- load_gts('NCF1-GTGT.csv.gz')

gts2 <- filter(gts, qual >= 10) |>
    count(pop, ac0) |>
    group_by(pop) |> mutate(perc = 100 * n / sum(n)) |> ungroup() |>
    complete(pop, ac0, fill = list(n = 0, perc = 0))

FONT <- 'Roboto'
fill_colors <- c('#001219', '#005f73', '#ca6702', '#ae2012') |>
    setNames(as.character(1:4))
(zoomed <- ggplot(filter(gts2, ac0 != 2)) +
    geom_bar(aes(pop, perc, fill = factor(ac0), group = ac0),
        stat = 'identity', position = position_dodge(width = 0.9),
        color = NA) +
    scale_y_continuous(NULL,
        breaks = seq(0, 100, 3), minor_breaks = NULL,
        expand = c(0, 0)) +
    scale_x_discrete(NULL, expand = expansion(add = 0.48)) +
    scale_fill_manual(values = fill_colors) +
    theme_bw() +
    theme(
        text = element_text(family = FONT),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none',
    ))

(full <- ggplot(gts2) +
    geom_bar(aes(pop, perc, fill = factor(ac0), group = ac0),
        stat = 'identity', position = position_dodge(width = 0.9),
        color = NA) +
    scale_y_continuous('Percentage of samples',
        breaks = seq(0, 100, 20),
        expand = c(0, 0), limits = c(0, 100)) +
    scale_x_discrete(NULL, expand = expansion(add = 0.48)) +
    scale_fill_manual('Aggregate genotype', values = fill_colors,
        labels = function(x) sprintf('%s⫽K', x)) +
    theme_bw() +
    theme(
        text = element_text(family = FONT),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'top',
        legend.key.size = unit(0.5, 'lines'),
        legend.title = element_text(margin = margin(r = 8), vjust = 0.8),
        legend.text = element_text(family = 'Fira Sans Condensed',
            size = 8.5, margin = margin(l = 3, r = 4))
    ))
legend <- cowplot::get_plot_component(full, 'guide-box', return_all = T)

##############################

gts_wes <- load_gts('WES-NCF1-GTGT.csv.gz')
filt_gts_wes <- filter(gts_wes, ad0 <= 500 & cn == 6)

colors <- c('#001219', '#0a9396', '#ee9b00', '#ae2012') |>
    setNames(as.character(1:4))
(dotplot <- ggplot(filt_gts_wes) +
    geom_abline(slope = 5:2 / 1:4, color = colors, linewidth = 1.4, alpha = 0.3) +
    geom_point(aes(ad0, ad1, color = factor(ac0), shape = qual < 10), size = 1.) +
    scale_x_continuous('Reads supporting GTGT allele') +
    scale_y_continuous('Reads supporting GT allele') +
    scale_color_manual('Aggregate genotype',
        values = colors, labels = function(x) sprintf('%s⫽6', x)) +
    scale_shape_manual('Quality', values = c(19, 4), labels = c('High', 'Low')) +
    guides(
        color = guide_legend(order = 1, override.aes = list(size = 2)),
        shape = guide_legend(order = 2, override.aes = list(size = 2))) +
    theme_bw() +
    theme(
        text = element_text(family = FONT),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'top',
        legend.box = 'vertical',
        legend.spacing.y = unit(-0.25, 'lines'),
        legend.box.margin = margin(t = -4, b = -4),
        legend.margin = margin(r = 25),
        legend.background = element_blank(),
        # legend.title = element_text(margin = margin()),
        legend.text = element_text(
            family = 'Fira Sans Condensed', size = 8.5,
            margin = margin(l = -3, r = -2)),
    ))

ggdraw(xlim = c(0, 3.3), ylim = c(0, 1)) +
    draw_plot(full + theme(legend.position = 'none'),
        x = 0, y = 0.06, width = 1.0, height = 0.93) +
    draw_plot(zoomed, x = 1.1, y = 0.06, width = 0.9, height = 0.93) +
    draw_plot(legend, x = 0.7, y = -0.005, width = 0.4, height = 0.08) +
    draw_plot(dotplot, x = 2, y = 0, width = 1.3, height = 1) +
    annotate('segment', x = 1.0, xend = 1.1, y = 0.121, yend = 0.121,
        color = 'gray40', linetype = '32') +
    annotate('segment', x = 1.0, xend = 1.1, y = 0.248, yend = 0.96,
        color = 'gray40', linetype = '32')
ggsave('ncf1_b-d.svg', width = 8, height = 4, dpi = 600, scale = 1.05)
