#!/usr/bin/Rscript

pdf(NULL)
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))

freqs <- read.csv('frequencies.csv', sep = '\t') |>
    mutate(ac0 = str_count(haplotype, '0'), ac1 = str_count(haplotype, '1'))

ggplot(freqs) +
    geom_point(aes(ac1, log10_freq))

ggplot(freqs) +
    geom_histogram(aes(freq), binwidth = 0.01)

