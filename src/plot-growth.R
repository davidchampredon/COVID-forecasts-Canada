suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2)) 
theme_set(theme_bw())

source('utils-simple.R')
source('estim-R.R')

dir.plot = 'plots/'
dat.all  = get_all_data_csv()
provs    = unique(dat.all$Province)



pdf(paste0(dir.plot,'plot-data-dbl.pdf'), width = 12, height = 6)
plot_dbl_time()
dev.off()

pdf(paste0(dir.plot,'plot-data-dbl-slide.pdf'), width = 12, height = 8)
plot_dbl_sliding()
dev.off()

pdf(paste0(dir.plot,'plot-data-R0.pdf'), width = 12, height = 6)
plot_R0()
dev.off()

pdf(paste0(dir.plot,'plot-data-Rt.pdf'), 
    width = 8, height = 6)
plot_Rt()
dev.off()

pdf(paste0(dir.plot,'plot-data-R0-slide.pdf'), width = 12, height = 8)
plot_R0_sliding()
dev.off()

# pdf(paste0(dir.plot,'plot-data-growth-cone.pdf'), width = 8, height = 8)
# plot_growth_rate_cone()
# dev.off()

for(prov in provs){
  
  message(paste(' --- Plotting growth for', prov,'---'))
  
  pdf(paste0(dir.plot,'plot-data-',prov,'.pdf'), 
      width = 16, height = 8)
  w     = 14
  level = 0.95
  try({
    gridExtra::grid.arrange(
      plot_data_prov_dbl(prov, w),
      plot_dbl_evol(prov, w, level),
      nrow = 1)
  })
  dev.off()
  
  pdf(paste0(dir.plot,'plot-data-Rt-',prov,'.pdf'),
      width = 10, height = 4)
  plot_Rt_prov(prov)
  dev.off()
  
  pdf(paste0(dir.plot,'plot-regr-',prov,'.pdf'))
  plot_r_regression(prov)
  dev.off()
}


