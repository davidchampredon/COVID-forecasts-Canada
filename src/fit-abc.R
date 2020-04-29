###
###  ABC FIT OF THE SIMPLE MODEL
###

args     = commandArgs(trailingOnly = TRUE)
prov     = args[1]
scenario = args[2]  

if(0){ # - - - DEBUG
    message(' ********** WARNING: DEBUG  ON ************')
    prov = 'ON' # 'ON'
    scenario = 'ISO1' # 'BASELINE'
}

source('abc-fcts.R')
source('prm_load.R')

dir.out = 'out/'
ncpumax = parallel::detectCores()

# Behaviour functions w.r.t. scenario:
B.fct = get(paste0("behave_", scenario))
B.prm = get(paste0("B.", scenario))


fcst = forecast_unit_abc(prms, 
                         prms.fit, 
                         prms.prior, 
                         B.fct, 
                         B.prm)

# --- Dump outputs 

save_lambda_post(fcst, prov, scenario)

# ---- Plots 

fname = paste0('plots/abcfit-prm-',prov,'-',scenario,'.pdf')
pdf(fname)
plot_post_abc(fcst$abcfit)
dev.off()

fname = paste0('plots/abcfit-traj-',prov,'-',scenario,'.pdf')
pdf(fname, width = 12, height = 5)
plot_post_traj_pos(fcst, prms.fit, B.prm, positives)
dev.off()

fname = paste0('plots/abcfit-death-',prov,'-',scenario,'.pdf')
pdf(fname, width = 12, height = 5)
plot_post_traj_death(fcst, prms.fit, B.prm, deaths)
dev.off()



rdataname <- paste0('fcst-',prov,'-', scenario, '.RData')
save(list = c('fcst', 'prms', 'dat'), file = rdataname)

