#!/usr/bin/env Rscript
pdf(NULL)

suppressMessages(library(ggplot2))
suppressMessages(library(ggh4x))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(cowplot))
suppressMessages(library(bit64))

cowplot::set_null_device('agg')
FONT <- 'Source Sans 3'

cmp <- read.csv('cmp_WGS_WES.csv.gz', sep = '\t', comment = '#') |>
    filter(
        count > 50
        & ref_cn > 2
        & exons != '*'
        & pmax(wgs_AF, wes_AF) != 0
        & allele_ix > 0
        & !grepl(',', alts)
        & allele_ix == 1
        # & !psv
    ) |>
    rename(chr_pos = pos) |>
    separate(chr_pos, into = c('chr', 'pos'), sep = ':', convert = T, remove = F) |>
    mutate(rare = wgs_AF < 0.01)

# Select one position from each cluster.
decluster <- function(cmp, dist) {
    last <- as.integer64(-2e9)
    filter(cmp,
        sapply(ext_pos, function(x) {
           if (x - last >= dist) {
               last <<- x;
               T
           } else {
               F
           }
        }))
}

decluster <- function(cmp, n_points, diameter) {
    cmp <- mutate(cmp, ext_pos = as.numeric(factor(chr)) * as.integer64(1e9) + pos) |>
        arrange(ext_pos)
    pos <- cmp$ext_pos
    s <- length(pos)
    neighbors <- rep(0, s)
    i <- 1
    ignore <- sapply(1:s, function(j) {
        x <- pos[j]
        while (x - pos[i] > diameter) {
            i <<- i + 1
        }
        if (i < j) {
            j_1 <- j - 1
            neighbors[i:j_1] <<- neighbors[i:j_1] + 1
        }
        neighbors[j] <<- neighbors[j] + j - i + 1
        return(0)
    })
    cmp$ext_pos <- NULL
    cmp[neighbors < n_points,]
}

# cmp2 <- cmp
cmp2 <- decluster(cmp, 10, 150)

events <- read.csv('cmp_WGS_WES.events.csv.gz', sep = '\t') |>
    inner_join(select(cmp2, chr_pos, rare, ref_cn), join_by(pos == chr_pos)) |>
    mutate(
        wgs_ac = str_count(wgs_gt, '0'),
        wes_ac = str_count(wes_gt, '0'),
        gt_cn = str_count(wgs_gt, '/') + 1,
    )
events_aggr <- filter(events, gt_cn %in% c(4, 6)) %>%
    aggregate(count ~ gt_cn + rare + wgs_ac + wes_ac, ., sum) |>
    mutate(diagonal = wgs_ac == wes_ac)

draw_subplot <- function(events_aggr, cn, rare) {
    curr_events <- events_aggr[events_aggr$gt_cn == cn & events_aggr$rare == rare,]
    left <- cn == 4; right <- !left
    btm <- rare; top <- !btm
    
    minx <- 0
    maxx <- cn
    miny <- if (rare) { cn - 2 } else { 0 }
    maxy <- if (rare) { cn } else { 6 }
    
    ggplot(curr_events, aes(wes_ac, wgs_ac)) +
        geom_tile(aes(fill = diagonal)) +
        geom_text(aes(label = prettyNum(count, big.mark = ','), color = diagonal),
            size = 2, family = 'Fira Sans Condensed') +
        scale_fill_manual(values = chess_colors, guide = 'none') +
        scale_color_manual(values = c('black', chess_colors[1]), guide = 'none') +
        scale_x_continuous(NULL,
            labels = function(x) sprintf('%d/%d', x, cn),
            breaks = 0:cn, limits = c(minx - 0.5, maxx + 0.5), expand = expansion(),
            sec.axis = if (top) {
                dup_axis(name = sprintf('Copy number %d', cn), breaks = NULL)
                } else { waiver() }
            ) +
        scale_y_continuous(NULL,
            labels = function(x) sprintf('%d/%d', x, cn),
            breaks = 0:cn, limits = c(miny - 0.5, maxy + 0.5), expand = expansion(),
            sec.axis = if (right) {
                dup_axis(name = ifelse(rare, 'Rare\n(AF < 1%)', 'Common\n(AF ≥ 1%)'), breaks = NULL)
                } else { waiver() }
            ) +
        theme_bw() +
        theme(
            text = element_text(family = FONT),
            panel.border = element_rect(color = 'black', linewidth = 0.3),
            panel.grid = element_blank(),
            panel.background = element_rect(color = NA),
            axis.text.x = element_text(size = 7, margin = margin(t = 1)),
            axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5, margin = margin(r = 1)),
            axis.ticks = element_blank(),
            axis.title.x.bottom = element_blank(),
            axis.title.x.top = element_text(size = 10),
            axis.title.y.left = element_blank(),
            axis.title.y.right = element_text(size = 10),
            panel.spacing.x = unit(0.4, 'lines'),
            plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
        )
}

chess_colors <- c('#FAFED0', '#006837')
subplots <- list()
for (cn in c(4, 6)) {
    for (rare in 0:1) {
        key <- sprintf('%d%s', cn, ifelse(rare, 'R', 'C'))
        subplots[[key]] <- draw_subplot(events_aggr, cn, rare)
    }
}

x_lab <- ggplot() + 
  annotate(geom = 'text', x = 0, y = 0, label = 'WES allele dosage', family = FONT) +
  coord_cartesian(clip = 'off') +
  theme_void()
y_lab <- ggplot() + 
  annotate(geom = 'text', x = 0, y = 0, label = 'WGS allele dosage', family = FONT, angle = 90) +
  coord_cartesian(clip = 'off') +
  theme_void()

(g_ev <- plot_grid(
    subplots$`4C`, subplots$`6C`, subplots$`4R`, subplots$`6R`,
    ncol = 2,
    rel_widths = c(0.35, 0.65),
    rel_heights = c(0.7, 0.3)) %>%
    plot_grid(y_lab, ., rel_widths = c(0.05, 0.95)) %>%
    plot_grid(., x_lab, rel_heights = c(0.95, 0.05), ncol = 1))

#########################

sum_acs <- group_by(events, pos) |>
    summarise(
        wgs_ac = sum(wgs_ac * count),
        wes_ac = sum(wes_ac * count),
        total_alleles = sum(gt_cn * count),
        ref_cn = unique(ref_cn),
    ) |>
    mutate(
        wgs_af = wgs_ac / total_alleles,
        wes_af = wes_ac / total_alleles,
    )
r2 <- with(sum_acs, cor(wgs_af, wes_af)) ^ 2

# Reduce the number of points in the output SVG file
ROUND_FACTOR <- 0.001
groupped_afs <- mutate(sum_acs,
    wgs_af = round(wgs_af / ROUND_FACTOR) * ROUND_FACTOR,
    wes_af = round(wes_af / ROUND_FACTOR) * ROUND_FACTOR) |>
    count(wgs_af, wes_af)

# RColorBrewer::brewer.pal(9, 'YlGn')[9]
point_color <- '#004529'
(g_af <- ggplot(groupped_afs) +
    annotate('segment', x = 0, y = 0, xend = 1, yend = 1,
        color = point_color, alpha = 0.2) +
    annotate('text', x = 0.05, y = 0.99, label = sprintf('r² = %.4f', r2),
        hjust = 0, vjust = 1, family = FONT, size = 3.5) +
    geom_point(aes(1 - wes_af, 1 - wgs_af, alpha = pmin(0.2 * n, 1)),
        color = point_color, size = 1.5, show.legend = F) +
    coord_fixed(ratio = 1.2) +
    scale_x_continuous('WES alternate allele frequency',
        limits = c(0, 1), breaks = seq(0, 1, 0.2),
        expand = expansion(add = 0.02)) +
    scale_y_continuous('WGS alternate allele frequency',
        limits = c(0, 1), breaks = seq(0, 1, 0.2),
        expand = expansion(add = 0.02)) +
    scale_alpha_identity() +
    theme_bw() +
    theme(
        text = element_text(family = FONT),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = 5, r = 2, b = 1, l = 7),
    ))

cowplot::plot_grid(g_ev, g_af,
    rel_widths = c(0.58, 0.42),
    labels = letters,
    label_fontfamily = FONT,
    label_y = 1.01,
    label_x = c(0, 0)
    )
ggsave('concordance.svg', bg = 'white', width = 10, height = 5, dpi = 500, scale = 0.72)

#################

# group_by(events, gt_cn, rare) |>
#     filter(gt_cn %in% c(4, 6)) |>
#     summarize(
#         n_vars = length(unique(pos)),
#         n_gts = sum(count),
#         concordance = sum((wes_ac == wgs_ac) * count) / n_gts, .groups = 'keep') |>
#     ungroup() |>
#     as.data.frame()
# 
# filter(events, gt_cn %in% c(4, 6) & rare) |>
#     summarize(n_vars = length(unique(pos)), n_gts = sum(count))
