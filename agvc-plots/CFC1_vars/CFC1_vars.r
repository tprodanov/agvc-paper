#!/usr/bin/env Rscript
pdf(NULL)

suppressMessages(library(ggplot2))
suppressMessages(library(ggchicklet))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))

source('spaceout_points.r')

exons <- read.csv('exons.bed', sep = '\t', header = F) |>
    setNames(c('chrom', 'start', 'end'))
START <- 130592876 # 0-based left-most position of the end codon
END <- 130599287 # 1-based right-most position of the start codon

aa_to_x <- function(aa, exon) {
    aa + exon * 15
}

exons <- arrange(exons, start) |>
    mutate(,
        exon_id = n():1,
        cds_start = pmax(START, start),
        cds_end = pmin(END, end),
        n_aa = (cds_end - cds_start) / 3,
        end_aa = rev(cumsum(rev(n_aa))),
        start_aa = c(end_aa[2:n()], 0),
        
        xstart = aa_to_x(start_aa, exon_id),
        xend = aa_to_x(end_aa, exon_id),
    )


add_prefix <- function(prefix) {
    function(x) paste0(prefix, '_', x)
}

variants <- bind_rows(
    read.csv('lof_vars.csv', sep = '\t') |> rename_at(vars(-impact), add_prefix('lof')),
    read.csv('common_vars.csv', sep = '\t') |> rename_at(vars(-impact), add_prefix('common')),
    read.csv('rare_vars.csv', sep = '\t') |> rename_at(vars(-impact), add_prefix('rare'))
)

aa_table <- read.csv('amino_acid_table.csv', sep = ',')
aa_sub <- with(aa_table, setNames(X1_letter_code, tolower(X3_letter_code)))

POS_PATTERN <- '^[^0-9]*([0-9]+)'
variants <- mutate(variants,
    aa_pos = as.numeric(sub('^[^0-9]*([0-9]+).*', '\\1', impact)),
    exon = findInterval(aa_pos, rev(exons$start_aa)),
    x = aa_to_x(aa_pos, exon),
    impact_short = stringr::str_replace_all(tolower(impact), aa_sub) %>%
        sub('_', '-\n', .),
    type = case_when(
        grepl('fs', impact) ~ 'Frameshift',
        grepl('*', impact, fixed = T) ~ 'Stop gained',
        T ~ 'Missense'
    ) |> factor(levels = c('Missense', 'Frameshift', 'Stop gained')),
    top = type == 'Missense',
    appear_case = pmax(common_af_case, rare_ac_case, na.rm = T) |> replace_na(0) > 0,
    appear_control = pmax(common_af_control, rare_ac_control, rare_ac_1kgp, na.rm = T) |>
        replace_na(0) > 0,
) # |> group_by(impact_short) |> slice_head(n = 1) |> ungroup()

xlim <- range(variants$x) + c(-10, 70)
SPACING <- 5
vars_top <- filter(variants, top) |>
    mutate(text_x = spaceout_points(x,
        limit = xlim, iterations = 1000, spacing = SPACING, force = 0.1, attraction = 0.001))
vars_btm <- filter(variants, !top) |>
    mutate(text_x = spaceout_points(x,
        limit = xlim, iterations = 1000, spacing = SPACING, force = 1, attraction = 0.1))

xbreaks <- with(exons,
    data.frame(aa = c(floor(end_aa), ceiling(start_aa)), exon_id = exon_id)) |>
    mutate(x = aa_to_x(aa, exon_id))
colors <- c('Missense' = '#5954d5', 'Frameshift' = '#E91E63', 'Stop gained' = '#FF6E00')

CIRCLE_NUDGE <- 1.05
CIRCLE_SIZE <- 5
SEMICIRCLE_L <- '\u25D6' # Left semicircle
SEMICIRCLE_R <- '\u25D7' # Right semicircle

FONT <- 'Roboto'

EXON_Y <- 0.4
LINESTART_Y <- 0.39
LINEEND_Y <- 1.5
CIRCLE_Y <- 1.65
TEXT_Y <- 1.8
LIMIT_Y <- 2.45

ggplot(exons) +
    geom_hline(yintercept = 0, color = 'gray10', linetype = '11') +
    geom_curve(
        data = vars_btm,
        aes(x = x, xend = text_x, y = -LINESTART_Y, yend = -LINEEND_Y),
        curvature = -0.1, ncp = 3, angle = 160,
        ) +
    geom_curve(
        data = vars_top,
        aes(x = x, xend = text_x, y = LINESTART_Y, yend = LINEEND_Y),
        curvature = -0.1, ncp = 3, angle = 160,
        ) +
    geom_point(
        data = vars_btm,
        aes(x = text_x + CIRCLE_NUDGE, y = -CIRCLE_Y, shape = 'Control/1KGP'),
        size = CIRCLE_SIZE,
        ) +
    geom_point(
        data = filter(vars_top, appear_control),
        aes(x = text_x + CIRCLE_NUDGE, y = CIRCLE_Y, shape = 'Control/1KGP'),
        size = CIRCLE_SIZE,
        ) +
    geom_point(
        data = filter(vars_top, appear_case),
        aes(x = text_x - CIRCLE_NUDGE, y = CIRCLE_Y, shape = 'Case'),
        size = CIRCLE_SIZE,
        ) +
    geom_label(
        data = vars_btm,
        aes(x = text_x, y = -TEXT_Y, label = impact_short, color = type),
        angle = 90, hjust = 1, family = FONT, show.legend = F,
        lineheight = 0.6, size = 2.5,
        border.color = NA,
        ) +
    geom_label(
        data = vars_top,
        aes(x = text_x, y = TEXT_Y, label = impact_short, color = type),
        angle = 90, hjust = 0, family = FONT, show.legend = F,
        lineheight = 0.6, size = 2.5,
        linewidth = ifelse(!is.na(vars_top$common_af_1kgp), 0.2, 0),
        label.padding = unit(0.2, 'lines'), label.r = unit(0.2, 'lines'),
        ) +
    geom_rrect(aes(xmin = xstart, xmax = xend, ymin = -EXON_Y, ymax = EXON_Y),
        fill = 'gray10') +
    geom_segment(
        aes(x = x, xend = x, y = -EXON_Y, yend = EXON_Y, color = type),
        data = variants) +
    scale_x_continuous('Amino acid',
        breaks = xbreaks$x, labels = xbreaks$aa, minor_breaks = NULL,
        expand = expansion(add = c(3, 3))) +
    scale_y_continuous(NULL, breaks = NULL,
        limits = c(-LIMIT_Y, LIMIT_Y),
        expand = expansion(add = 0.01)) +
    scale_color_manual('Impact  ',
        breaks = names(colors),
        values = colors,
        guide = guide_legend(order = 1),
        labels = function(x) sprintf(' %s ', x),
        ) +
    scale_shape_manual('   Group',
        values = c(SEMICIRCLE_L, SEMICIRCLE_R),
        guide = guide_legend(order = 2),
        ) +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = FONT),
        legend.position = 'top',
        legend.margin = margin(b = -10),
        legend.box = 'horizontal',
        legend.title = element_text(size = 9.2, margin = margin(l = 5, r = 2), face = 'bold'),
        legend.text = element_text(size = 9, margin = margin(l = 0, r = 0)),
    )
ggsave('CFC1_vars.svg', width = 6, height = 3, dpi = 500, scale = 1.1)
