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

variables <- read.csv('vars.csv', sep = '\t')

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

aa_table <- read.csv('amino_acid_table.csv', sep = ',')
aa_sub <- with(aa_table, setNames(X1_letter_code, tolower(X3_letter_code)))

POS_PATTERN <- '^[^0-9]*([0-9]+)'
variables <- mutate(variables,
    aa_pos = as.numeric(sub('^[^0-9]*([0-9]+).*', '\\1', impact)),
    exon = findInterval(aa_pos, rev(exons$start_aa)),
    x = aa_to_x(aa_pos, exon),
    impact_short = stringr::str_replace_all(tolower(impact), aa_sub) %>%
        sub('_', '-\n', .),
    type2 = recode_factor(type,
        'missense' = 'Missense', 'deletion' = 'Deletion',
        'frameshift' = 'Frameshift', 'stop-gained' = 'Stop gained'),
    top = as.numeric(type2) <= 2,
) |> group_by(impact) |>
    slice_max(replace_na(ac, 0) + replace_na(ac_pcgc, 0), n = 1, with_ties = F) |>
    ungroup()

xlim <- range(variables$x) + c(-10, 70)
SPACING <- 5
vars_top <- filter(variables, top) |>
    mutate(text_x = spaceout_points(x,
        limit = xlim, iterations = 1000, spacing = SPACING, force = 1, attraction = 0.001))
vars_btm <- filter(variables, !top) |>
    mutate(text_x = spaceout_points(x,
        limit = xlim, iterations = 1000, spacing = SPACING, force = 1, attraction = 0.1))

xbreaks <- with(exons,
    data.frame(aa = c(floor(end_aa), ceiling(start_aa)), exon_id = exon_id)) |>
    mutate(x = aa_to_x(aa, exon_id))
colors <- c(
    'Missense' =  '#5954d5',
    'Deletion' = '#ebac23',
    'Frameshift' = '#c7332f',
    'Stop gained' = '#ff9287'
)
EXON_Y1 <- -0.4
EXON_Y2 <- 0.4
Y_SHIFT <- 1.5

FONT <- 'Roboto'
ggplot(exons) +
    geom_hline(yintercept = 0, color = 'gray10', linetype = '11') +
    geom_curve(
        data = vars_btm,
        aes(x = x, xend = text_x, y = EXON_Y1 + 0.01, yend = EXON_Y1 - Y_SHIFT + 0.05),
        curvature = -0.1, ncp = 3, angle = 160,
        ) +
    geom_curve(
        data = vars_top,
        aes(x = x, xend = text_x, y = EXON_Y2 - 0.01, yend = EXON_Y2 + Y_SHIFT - 0.05),
        curvature = -0.1, ncp = 3, angle = 160,
        ) +
    geom_text(aes(x = text_x, y = EXON_Y1 - Y_SHIFT, label = impact_short, color = type2),
        data = vars_btm,
        angle = 90, hjust = 1, family = FONT, show.legend = F,
        lineheight = 0.6, size = 2.5,
        ) +
    geom_text(aes(x = text_x, y = EXON_Y2 + Y_SHIFT, label = impact_short, color = type2),
        data = vars_top,
        angle = 90, hjust = 0, family = FONT, show.legend = F,
        lineheight = 0.6, size = 2.5,
        ) +
    geom_rrect(aes(xmin = xstart, xmax = xend, ymin = EXON_Y1, ymax = EXON_Y2),
        fill = 'gray10') +
    geom_segment(
        aes(x = x, xend = x, y = EXON_Y1, yend = EXON_Y2, color = type2),
        data = variables) +
    scale_x_continuous('Amino acid',
        breaks = xbreaks$x, labels = xbreaks$aa, minor_breaks = NULL,
        expand = expansion(add = c(3, 3))) +
    scale_y_continuous(NULL, breaks = NULL,
        limits = c(-1, 1) * (EXON_Y2 + Y_SHIFT + 0.6),
        expand = expansion(add = 0.01)) +
    scale_color_manual(NULL, values = colors) +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = FONT),
        legend.position = 'top',
        legend.margin = margin(b = -15),
        legend.box = 'vertical',
        legend.title = element_text(size = 9, face = 'bold')
    )
ggsave('CFC1_vars.svg', width = 6, height = 3, dpi = 500, scale = 1.2)
