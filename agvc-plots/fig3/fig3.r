suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))

dupl_genes <- readLines('dupl_genes.txt')
filenames <- Sys.glob('evals/*.csv.gz')

full_roc <- lapply(filenames,
    function(filename) read.csv(filename, sep = '\t', comment = '#') |>
        mutate(filename = filename)) %>%
    do.call(rbind, .)

# Splitting filenames
full_roc2 <- mutate(full_roc,
    filename = sub('.csv.gz$', '', filename) %>% sub('evals/', '', .) %>% sub('parascopy', 'parascopy.NA', .)) |>
    separate(filename, c('sample', 'tool', 'tool_qual', 'cn'))

# Trust Parascopy qual column, for other tools, use quality from the `QX` part of the filename.
full_roc3 <- filter(full_roc2,
    region != '*'
    & (tool == 'parascopy' | (tool_qual == 'Q1' & qual == 'any') | (tool_qual == 'Q10' & qual == 'high')))

# Split region into gene & PSV type
full_roc4 <- mutate(full_roc3,
    region = sub('-(non-psv|non-trivial_psv|trivial_psv)$', '@\\1', region)) |>
    separate(region, c('gene', 'psv_type'), sep = '@')

# Take all variants together, high quality, duplicated genes
roc <- filter(full_roc4,
        qual == 'high'
        & var_type == 'all'
        & sample != 'SimFixed'
        & gene %in% dupl_genes) |>
    select(-c(qual, tool_qual, var_type)) |>
    mutate(sample_type = ifelse(startsWith(sample, 'Sim'), 'Simulated', 'GIAB'))

# Sum by gene, PSV type, tool and sample (some genes may have multiple entries for different CNs)
roc2 <- group_by(roc, gene, psv_type, tool, sample, sample_type) |>
    summarize_at(c('tp_base', 'tp_call', 'fp', 'fn'), sum) |>
    ungroup()

# Take mean across all samples
roc3 <- group_by(roc2, gene, psv_type, tool, sample_type) |>
    summarize_at(c('tp_base', 'tp_call', 'fp', 'fn'), mean) |>
    ungroup() |>
    mutate(
        n_vars = tp_base + fn,
        precision = replace_na(tp_call / (tp_call + fp), 0),
        recall = replace_na(tp_base / (tp_base + fn), 0),
        f1 = replace_na(2 * precision * recall / (precision + recall), 0),
    )

roc4 <- select(roc3, gene, psv_type, tool, sample_type, n_vars, precision, recall, f1) |>
    pivot_wider(names_from = 'tool',
        values_from = c('n_vars', 'precision', 'recall', 'f1'), names_sep = '.') |>
    filter(!is.na(n_vars.parascopy) & !is.na(n_vars.freebayes))

MIN_VARS <- 10
roc5 <- filter(roc4, pmin(n_vars.parascopy, n_vars.freebayes) >= MIN_VARS) |>
    mutate(
        recall.improv = recall.parascopy - recall.freebayes,
        f1.improv = f1.parascopy - f1.freebayes,
        ) |>
    group_by(psv_type, sample_type) |>
    arrange(-recall.improv) |>
    mutate(ix = 1:n()) |>
    ungroup()

roc_giab <- filter(roc5, sample_type == 'GIAB' & psv_type != 'trivial_psv')
max_ix <- max(roc_giab$ix)
y_range <- range(roc_giab$recall.improv)

roc_giab2 <- mutate(roc_giab, psv_type = recode_factor(psv_type, 'non-psv' = 'Non-PSV',
        'non-trivial_psv' = 'PSV'))
colors <- colorspace::sequential_hcl(7, palette = 'Emrld')[c(2, 4)]

ggplot(roc_giab2) +
    geom_histogram(aes(recall.improv, fill = psv_type, color = psv_type,
        y = ifelse(fill == 'Non-PSV', 1, -1) * after_stat(count)),
        position = 'identity', binwidth = 0.035, center = 0, linewidth = 0.1) +
    scale_x_continuous('Improvement in recall', breaks = seq(0, 0.8, 0.2),
        expand = expansion(mult = 0.01)) +
    scale_y_continuous('Number of genes', labels = abs, expand = expansion(mult = 0.02)) +
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
ggsave('improv_recall.svg', width = 10, height = 5, scale = 0.7)
