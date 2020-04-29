suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2)) 
theme_set(theme_bw())

source('utils-simple.R')

dat.all  = get_all_data_csv()
dir.plot = 'plots/'

pdf(paste0(dir.plot,'plot-data.pdf'), width = 24, height = 6)
dat.all %>%
  plot_all_data() %>%
  plot()
dev.off()


pdf(paste0(dir.plot,'plot-data-canada.pdf'), width = 12, height = 6)
dat.all %>%
  plot_data_report_canada() 
dev.off()


