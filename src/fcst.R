source('fcst-fcts.R')

args <- commandArgs(trailingOnly = TRUE)

prov      = args[1]
scenario  = args[2]  

# - - - DEBUG
if(0){
    message(' ********** WARNING: DEBUG  ON ************')
    prov = 'ON'
    scenario = 'BASELINE'
}

# Call to load all parameters.
# Must be after `prov` is defined!
source('prm_load.R')

# ---- RUN ----

message(paste("Forecasting", prov, scenario, "..."))

fcst <- forecast_unit(prms)

rdataname <- paste0('fcst-',prov,'-', scenario, '.RData')

save(list = c('fcst', 'prms', 'dat'), file = rdataname)

