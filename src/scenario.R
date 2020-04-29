library(stringr)


scenario_matcher <- function(scen) {
    
    # Check if the scenario name is known:
    check <- FALSE
    if(grepl('BASELINE', scen)) check <- TRUE
    if(grepl('^INST_RED', scen)) check <- TRUE
    if(grepl('^RED', scen)) check <- TRUE
    if(grepl('^KAPPA', scen)) check <- TRUE
    if(grepl('^PULSE', scen)) check <- TRUE
    stopifnot(check)
    
    x <- list()
        
    # Switch between the scenario names:
    
    if( grepl('^INST_RED_\\d{2}$',scen) ){
        r = as.numeric(stringr::str_extract(scen,'\\d{2}$') )
        x[['R0.factor']] = 1 - r/100
    }
    
    # R0 reduces exponential by XX% in YY days
    if( grepl('^RED_\\d{2}_\\d{2}d$',scen) ){     # scen = "RED_30_07d"
        y     = stringr::str_extract_all(scen,'\\d{2}') 
        reduc = as.numeric(y[[1]][1]) / 100
        b     = as.numeric(y[[1]][2])
        mult  = 1 -reduc
        x[['kappa']]  <- -log(mult)/b
        x[['lambda']] <- mult
    }

    if(scen == "KAPPA_30") {
        x[['kappa']] = log(2)/30
    }
     
    if(scen == "KAPPA_90") {
        x[['kappa']] = log(5/4)/90
    }   
    
    
    if( grepl('^PULSE_\\d{2}d_\\d{2}d_\\d{2}$',scen) ) {
        # Example with  "PULSE_14d_07d_50" 
        # R0 pulsed:
        # 14 days no pulse
        # then after 14 days, pulse of 50% (relative) for 7d
        # then back to no pulse for 14 days, etc.
        y = stringr::str_extract_all(scen,'\\d{2}') 
        x[['pulse.quiet.length']] <- as.numeric(y[[1]][1])
        x[['pulse.length']]       <- as.numeric(y[[1]][2])
        x[['pulse.height']]       <- as.numeric(y[[1]][3])/100
    }
    
    return(x)
}

#' Converts the scenario to more human readable for plots.
convert_scen <- function(x) {
    res = x
    
    if(x=='BASELINE') res <- 'Baseline'
    
    if( grepl('^INST_RED_\\d{2}$',x) ){
        r = stringr::str_extract(x,'\\d{2}$') 
        res <- paste0('Instant ',r,'%')
    }    
    
    if( grepl('^RED_\\d{2}_\\d{2}d$',x) ){     
        y     = stringr::str_extract_all(x,'\\d{2}') 
        reduc = y[[1]][1]
        b     = y[[1]][2]
        res <- paste0(reduc,'%', ' in ',b,' days')
    }
    
    if( grepl('^PULSE_\\d{2}d_\\d{2}d_\\d{2}$',x) ) {
        y = stringr::str_extract_all(x,'\\d{2}') 
        d1 <- as.numeric(y[[1]][1])
        d2 <- as.numeric(y[[1]][2])
        h  <- as.numeric(y[[1]][3])
        res <- paste0('Pulse ',d1,'/',d2,' ',h,'%')
    }
    
    return(res)
}

# DEBUG
# scen = c('BASELINE', 'BASELINE', 
#          'INST_RED_20', 'INST_RED_20',
#          'INST_RED_50',
#          'RED_30_14d', 'RED_30_14d', 'RED_30_14d',
#          'RED_50_28d', 'RED_50_28d',
#          'KAPPA_39', 'TOTO')

scenario_display <- function(scen) {
    res <- sapply(scen, convert_scen)
    return(res)   
}