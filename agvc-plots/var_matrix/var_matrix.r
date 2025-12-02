#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tibble))
suppressMessages(library(tidyr))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(require(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(glue))
ht_opt$message <- F

# ------ Constants ------

GENE <- 'CFC1'
# Cluster variants?
CLUSTER_VARS <- T

CN <- 4
# Type of variant: all/psvs
VAR_TYPE <- 'psvs'
# Remove variants with allele fraction < AF or > 1 - AF
AF <- 0.01
# Take samples that have appropriate CN for this fraction of variants
CN_MATCH <- 0.98
# How many variants must variant have
VAR_AVAIL_THRESH <- 0.95
# How many variants must sample have
SAMPLE_AVAIL_THRESH <- 0.95

# Take only these populations
POPULATIONS <- c('EUR', 'AFR')
# Before plotting, subsample to this number of individuals
SUBSAMPLE <- 100
SEED <- 142198621

# ------ Load and process variants ------

print_count <- function(vars) {
    cat(sprintf('%.0fk entries, %s samples, %s variants\n',
        nrow(vars) / 1000, length(unique(vars$sample)), length(unique(vars$pos))))
    vars
}

vars <- read.csv(sprintf('%s/%s.csv.gz', GENE, VAR_TYPE),
        sep = '\t', header = F) |>
    setNames(c('pos', 'pop', 'sample', 'gt')) |>
    mutate(
        ac = ifelse(gt == '.', NA, str_count(gt, '0')),
        cn = ifelse(gt == '.', NA, str_count(gt, '/') + 1),
    ) |>
    print_count()

if (GENE == 'SMN1') {
    start <- -Inf; end <- 70956711
} else {
    start <- -Inf; end <- Inf
}
filt_vars <- filter(vars, between(pos, start, end)) |> print_count() |>
    filter(pop %in% POPULATIONS) |> print_count()
vars_complete <- complete(filt_vars,
        pos, nesting(pop, sample), fill = list(ac = NA, cn = NA)) |>
    print_count()

# Remove samples with copy number variants
cn_stable_samples <- filter(vars_complete, !is.na(cn)) |>
    group_by(sample) |>
    summarize(cn_frac = sum(cn == CN) / n()) |>
    ungroup() |>
    filter(cn_frac >= CN_MATCH) |>
    with(as.vector(sample))
vars_stablecn <- filter(vars_complete, sample %in% cn_stable_samples) |>
    mutate(ac = ifelse(cn == CN, ac, NA)) |>
    print_count()

# Remove variants and samples with many unavailable values
vars_avail <- group_by(vars_stablecn, pos) |>
    filter(mean(!is.na(ac)) >= VAR_AVAIL_THRESH) |>
    ungroup() |> print_count() |>
    group_by(sample) |>
    filter(mean(!is.na(ac)) >= SAMPLE_AVAIL_THRESH) |>
    ungroup() |> print_count()

# Remove variants with low/high AF
vars_common <- group_by(vars_avail, pos) |>
    filter(between(sum(ac, na.rm = T) / sum(cn, na.rm = T), AF, 1 - AF)) |>
    ungroup() |> print_count()

# Subsample samples
set.seed(SEED)
sel_samples <- select(vars_common, pop, sample) |>
    arrange(sample) |> unique() |>
    group_by(pop) |>
    reframe(sample = base::sample(sample, SUBSAMPLE)) |>
    ungroup()

# Convert into matrix
vars_subsample <- filter(vars_common, sample %in% sel_samples$sample)
ac_matrix <- select(vars_subsample, pos, sample, ac) |>
    pivot_wider(names_from = 'sample', values_from = 'ac') |>
    mutate(pos = format(pos, big.mark = ',')) |>
    column_to_rownames('pos') |>
    as.matrix()
n_vars <- nrow(ac_matrix)
n_samples <- ncol(ac_matrix)

# Colors

colors_to_ramp <- function(colors, min_v=NULL, max_v=NULL, data=NULL) {
    if (is.null(min_v)) {
        min_v <- min(data)
    }
    if (is.null(max_v)) {
        max_v <- max(data)
    }
    n <- length(colors)
    colorRamp2(seq(min_v, max_v, length = n), colors)
}

var_color_ramp <- RColorBrewer::brewer.pal(9, 'PuRd') |> colors_to_ramp(-0.3 * n_vars, n_vars)
heatmap_ramp <- RColorBrewer::brewer.pal(9, 'YlGnBu') |> rev() |>
    colors_to_ramp(0, CN + 0.5)
gt_colors <- structure(heatmap_ramp(0:CN), names = as.character(0:CN))

# ------ Cluster variants ------

var_order <- if (CLUSTER_VARS) { hclust(dist(ac_matrix))$order } else { 1:n_vars }
ac_matrix2 <- ac_matrix[var_order,]
var_colors <- var_color_ramp(1:n_vars)[var_order]

# ------ Draw the matrix ------

FILENAME <- 'plots/{gene}_{plot}_{pop}_{vartype}_AF{af}.{ext}'
for (curr_pop in POPULATIONS) {
    out_filename <- glue(FILENAME, gene = GENE, plot = 'matrix',
        pop = curr_pop, vartype = VAR_TYPE, af = AF, ext = 'svg')
    svg(out_filename, width = 5, height = 6, family = 'Roboto', bg = 'transparent')
    
    pop_matrix <- ac_matrix2[,filter(sel_samples, pop == curr_pop)$sample]
    Heatmap(pop_matrix,
        heatmap_legend_param = list(
            title = gt_render('Aggregate genotype'),
            title_gp = gpar(fontface = 'plain'),
            direction = 'horizontal',
            title_position = 'topcenter',
            nrow = 1,
            at = 0:CN,
            labels = sprintf('%d⫽%d', 0:CN, CN)
            ),
        col = gt_colors,
        na_col = 'deeppink2',
        use_raster = F,
        raster_quality = 4,
        row_names_gp = grid::gpar(fontsize = 7,
            col = if (CLUSTER_VARS) { var_colors } else { 'black' }),
        row_names_side = 'right',
        row_dend_side = 'right',
        rect_gp = gpar(col = "white", lwd = 0.2),
        show_column_names = F,
        column_title_side = 'bottom',
        row_title = NULL,
        cluster_rows = F
        ) |>
    draw(heatmap_legend_side = 'top')
    while (dev.cur() > 1) { dev.off() }
    
    # Add proxy column with average genotype AAaa to ground correlations
    pop_matrix2 <- cbind(pop_matrix, rep(CN / 2, nrow(pop_matrix)))
    var_cor <- cor(t(pop_matrix2), use = 'pairwise.complete.obs') |>
        as.data.frame() |>
        rownames_to_column(var = 'var1') |>
        pivot_longer(-var1, names_to = 'var2', values_to = 'r') |>
        mutate(
            var1 = factor(var1, levels = rownames(pop_matrix)),
            var2 = factor(var2, levels = rownames(pop_matrix)),
            val = r * r)
    
    ggplot(filter(var_cor, as.numeric(var1) < as.numeric(var2))) +
        geom_tile(aes(var1, var2, fill = val)) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_fill_gradient(NULL,
            low = 'white', high = '#ff3f00',
            na.value = 'white', limits = c(0, 1),
            breaks = seq(0, 1, 0.2)) +
        coord_equal() +
        theme_minimal() +
        theme(
            text = element_text(family = 'Roboto'),
            plot.background = element_rect(fill = NA, color = NA),
            panel.border = element_rect(fill = NA, color = NA),
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            legend.position = 'inside',
            legend.position.inside = c(0.995, 0.03),
            legend.justification.inside = c('right', 'bottom'),
        )
    out_filename <- glue(FILENAME, gene = GENE, plot = 'cor',
        pop = curr_pop, vartype = VAR_TYPE, af = AF, ext = 'svg')
    ggsave(out_filename, bg = 'transparent', width = 5, height = 5)
}
