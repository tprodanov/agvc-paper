#!/usr/bin/Rscript

pdf(NULL)
suppressMessages(library(ggplot2))
suppressMessages(library(colorspace))

colors <- sequential_hcl(7, 'Dark mint')[5:1]

FONT <- 'Roboto'
for (i in 1:2) {
    df <- read.csv(sprintf('data%d.csv', i), sep = '\t')
    ggplot(df) +
        geom_bar(aes(ref_dosage, perc, fill = factor(ref_dosage)), stat = 'identity') +
        facet_grid(. ~ group) +
        scale_x_continuous('Aggregate genotype',
            expand = expansion(add = 0.2), breaks = 0:4,
            labels = function (x) sprintf('%d⫽4', x)) +
        scale_y_continuous(expand = expansion(add = c(1, 5))) +
        scale_fill_manual(values = colors) +
        theme_bw() +
        theme(
            text = element_text(family = FONT),
            panel.border = element_rect(color = 'black'),
            panel.grid = element_blank(),
            axis.title.x = element_text(size = 9),
            axis.text.x = element_text(family = 'Fira Sans Condensed', size = 7.5),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            strip.background = element_rect(color = NA, fill = 'white'),
            strip.text = element_text(size = 9, margin = margin(2, 2, 2, 2)),
            legend.position = 'none',
            plot.margin = margin(1, 1, 1, 1),
        )
    ggsave(sprintf('subplot%s.svg', i),
        width = 4.5, height = 2, scale = 0.65)
}
