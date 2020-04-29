suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2)) 
theme_set(theme_bw())

source('utils-simple.R')


add_type <- function(x, L) {
    y = L[[ x ]] 
    y$type = x
    return(y)
}

merge_sims <- function(sims){
    X   <- c('inc', 'testpos', 'hosp','critical', 'death')    
    tmp <- lapply(X, add_type, L=sims)
    df  <- do.call('rbind', tmp)
    return(df)
}

merge_digest <- function(digest) {
    X = c('inc', 'testpos', 'hosp','critical', 'death')    
    
    tmp.ss <- lapply(X, add_type, L=digest[[ 'ss' ]])    
    tmp.pk <- lapply(X, add_type, L=digest[[ 'peaks' ]])    
    tmp.cu <- lapply(X, add_type, L=digest[[ 'cum' ]])    
    
    return(list(
    ss    = do.call('rbind', tmp.ss),
    peaks = do.call('rbind', tmp.pk),
    cum   = do.call('rbind', tmp.cu)    
    ))
}


abcpost_to_df <- function(abc.post, prov, scenario) {
    # prov = 'AB'
    # scenario = 'ISO1'
    # abcpost =  tmp.abc.post[[1]] 
    
    df = data.frame(abc.post) %>%
        select(-idx) %>%
        pivot_longer(cols = everything()) %>%
        mutate(Province = as.character(prov), 
               scenario = as.character(scenario))
    return(df)
}


# ==== RUN ====

rdatas <- system('ls fcst*.RData', intern = TRUE)
n      <- length(rdatas)

prms.all      <- list()
tmp.sims      <- list()
tmp.digest.ss <- list()
tmp.digest.pk <- list()
tmp.digest.cu <- list()
tmp.data      <- list()
tmp.abc.post  <- list()


for(i in 1:n){     # i=1
    message(paste("Digesting", rdatas[i],'...'))
    load(rdatas[i])
    
    prms.all[[i]] <- prms
    tmp.sims[[i]] <- merge_sims(sims = fcst$sims)
    dgst          <- merge_digest(digest = fcst$digest)
    
    tmp.digest.ss[[i]] <- dgst$ss
    tmp.digest.pk[[i]] <- dgst$peaks
    tmp.digest.cu[[i]] <- dgst$cum
    
    tmp.abc.post[[i]] <- fcst$abcfit$posteriors %>%
        abcpost_to_df(prov = prms$prov, scenario = prms$scenario)
    
    # Need this condition, else data duplicated 
    # as many time as the number of scenario:
    if(grepl('BASELINE', rdatas[i])) 
        tmp.data[[i]] <- dat
}

sims         = do.call('rbind', tmp.sims)
digest.ss    = do.call('rbind', tmp.digest.ss)
digest.peaks = do.call('rbind', tmp.digest.pk)
digest.cum   = do.call('rbind', tmp.digest.cu)
abc.post     = do.call('rbind', tmp.abc.post)
dat.all      = do.call('rbind', tmp.data)

# --- Canada ----

# TODO: change the code such that 'Canada' is 
# _automatically_ consistent with the format for 
# individual provinces. 

message('Digesting Canada...')

sims.canada <- sims %>%
    group_by(iter, t, date, scenario, type) %>%
    summarise(value = sum(value)) %>%
    mutate(Province = 'Canada')
sims.canada$value[sims.canada$value==0] <- 0.4

ci <- c(0.50 , 0.95)

# `digest.ss`
q = ci_to_quantiles(ci)
digest.ss.canada <- sims.canada %>%
    group_by(t,date,scenario, type) %>%
    summarize(m = mean(value), 
              md = median(value),
              qvlo = quantile(value, q[1], na.rm = TRUE),
              qlo  = quantile(value, q[2], na.rm = TRUE),
              qhi  = quantile(value, q[3], na.rm = TRUE),
              qvhi = quantile(value, q[4], na.rm = TRUE)
              ) %>%
    mutate(Province = 'Canada')


# `digest.peaks.canada`
anchor.date = digest.peaks$pk.date[1] - digest.peaks$pkt[1]
digest.peaks.canada  <- sims.canada %>%
    group_by(iter, Province,scenario, type) %>%
    summarize(pkv = max(value),
              pkt = t[which.max(value)]) %>%
    mutate(pk.date = anchor.date + pkt)

# digest.cum.canada
digest.cum.canada <- sims.canada %>%
    group_by(iter, Province, scenario, type) %>%
    summarise(s = sum(value))


rdataname <- paste0('digest.RData')

save(list = c('sims', 
              'sims.canada',
              'digest.ss', 
              'digest.peaks', 
              'digest.cum',
              'digest.ss.canada', 
              'digest.peaks.canada', 
              'digest.cum.canada',
              'abc.post',
              'dat.all'), 
     file = rdataname)



