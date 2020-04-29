suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2)) 
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(gridExtra))
library(scales)
theme_set(theme_bw())

source('resude.R')
source('utils-simple.R')
source('estim-R.R')
source('scenario.R')

#' Retrieve data for tested positive:
get_data_testpos <- function(first.date) {
  dat <- get_all_data_csv() %>%
    filter(Province == prov, 
           name=='positive',
           Date >= first.date)
  return(dat)  
}

#' Retrieve data for deaths:
get_data_deaths <- function(first.date) {
  dat <- get_all_data_csv() %>%
    filter(Province == prov, 
           name=='death',
           Date >= first.date)
  return(dat)  
}



#' Re-scale "up" from positive tests 
#' to true incidence:
positives_to_true_incidence <- function(rho, positives, lag = 0) {
  
  np = length(positives)
  n.abc = length(rho)
  
  # Mean true incidences
  # Columns = times
  # Rows    = stochastic iterations
  M = (1.0/rho) %*% t(positives)
  
  # If reported = 0 this does not mean
  # there is no incidence at all. 
  # This is a quick fix (TODO: change that)
  M[M==0] <- 0.5
  
  # Add observation uncertainty:
  inc.init = matrix(rpois(n = n.abc * np, 
                          lambda = M), 
                    nrow=n.abc, byrow=F)
  return(inc.init)
}

#' Converts TRUE incidence to tested positive
convert_true_pos <- function(df, rho, cv, lag = 0) {
  res <- df
  # true incidence:
  inc <- res$value
  n   <- length(inc)
  
  # Retrieve the reporting proportions
  # used at each iteration of the simulations:
  rhos = rho[res$iter]
  
  # --- Positive ~ Gamma(mean = rho * true, CV):
  m = rhos * inc
  # Prevent mean=0 that makes Gamma crash:
  m[m<0.1] <- 0.1  
  # Convert parameterization:
  z = sapply(m, FUN = gamma_prm_conv, cv=cv)
  sh.vec <- z[1,]
  ra.vec <- z[2,]
  res$value <- rgamma(n=n, shape=sh.vec, rate=ra.vec)
  
  return(res)
}

#' Converts tested positive cases to hospitalized 
convert_pos_hosp <- function(df, h) {
  res <- df
  lambda <-  df$value * h
  lambda[lambda < 0.1] <- 0.1
  res$value <- rpois(n=nrow(df), lambda = lambda)
  return(res)
}

#' Converts hospitalized to critical cases:
convert_hosp_critical <- function(df, kappa) {
  res <- df
  lambda <-  df$value * kappa
  lambda[lambda < 0.1] <- 0.1
  res$value <- rpois(n=nrow(df), lambda = lambda)
  return(res)
}

#' Converts critical cases to death:
#' @param df Dataframe of critical cases.
#' @param d Numerical. Proportion of critical case that will die.
#' @param z Numerical. Temporal adjustment.
convert_critical_death <- function(df, d, death.adj = 0) {
  n = nrow(df)
  t = 1:n
  res <- df
  lambda <-  df$value * d * exp(death.adj*t)
  lambda[lambda < 0.1] <- 0.1
  res$value <- rpois(n, lambda = lambda)
  # plot(df$value, typ='l', log='y', ylim=c(0.01,1e5))
  # lines(res$value, typ='o')
  return(res)
}


#' Simulate epidemic forward in time:
sim_incidence <- function(n.mc, 
                          horiz,
                          pop.eff, 
                          R0.m, 
                          R0.sd,
                          scenario,
                          positives,
                          prop.true.pos.beta1,
                          prop.true.pos.beta2) {
  
  np = length(positives)
  sim <- matrix(nrow = n.mc, ncol=horiz+np)
  
  # Re-scale "up" from positive tests 
  # to true incidence:
  
  # Draw reporting proportion:
  rho <- rbeta(n = n.mc, 
               shape1 = prop.true.pos.beta1, 
               shape2 = prop.true.pos.beta2)
  
  # pdf('plot-rho.pdf')
  # hist(rho, breaks = seq(0,1,by=0.025), col='grey', xlim=c(0,1))
  # dev.off()
  
  # Mean true incidences:
  M = (1.0/rho) %*% t(positives)
  
  # Add observation uncertainty:
  inc.init = matrix(rpois(n=n.mc*np, lambda = M), nrow=n.mc, byrow=F)
  
  # Sample R0 that will be used in fwd sims
  # (R0 has been previously fitted on _positives_): 
  R0.smpl <- rnorm(n=n.mc, mean = R0.m, sd = R0.sd)
  
  # Generation interval parameters:
  gi <- get_GI_prms()
  
  # Human behaviour prms:
  alpha <- get_prm('alpha', 'prm-model.csv')
  kappa <- get_prm('kappa', 'prm-model.csv')
  R0.factor <- 1.0 # can be changed by a scenario
  lambda <- 0.0
  
  # Intervention:
  pulse.height = 0
  pulse.length = 0
  pulse.quiet.length = 0
  
  # * * * WARNING * * * 
  # Scenario matching MUST be the last call
  # before simulating to make sure the parameters
  # have their values corresponding to the scenario! 
  
  # Scenario (must be just before simulating):
  ps = scenario_matcher(scenario)
  nm = names(ps)
  for(i in seq_along(ps)) 
    assign(nm[i], ps[[nm[i]]] )
  
  
  # Simulate epidemic fwd in time:
  for(i in 1:n.mc)  # i=1
  {  
    # Remove already infected from effective population size:
    pop_size = pop.eff - sum(inc.init[i,])
    
    tmp <- resude_simulate(pop_size = pop_size, 
                           I.init  = inc.init[i,], 
                           R0      = R0.factor * R0.smpl[i],
                           alpha   = alpha, 
                           kappa   = kappa,
                           GI_span = round(gi$mean + 2*gi$variance), 
                           GI_mean = gi$mean, 
                           GI_var  = gi$variance,
                           horizon = horiz, 
                           distrib.proc = 'rgamma',
                           phi     = gi$gam.disp, 
                           lambda  = lambda,
                           pulse.height       = pulse.height,
                           pulse.length       = pulse.length,
                           pulse.quiet.length = pulse.quiet.length,
                           seed = i
    )
    # This is the _true_ incidence:
    sim[i, ] <- tmp$I
  }
  tmp <- as.data.frame(sim)
  tmp$iter = 1:n.mc
  df <- pivot_longer(tmp, cols = starts_with('V')) %>%
    mutate(t = as.numeric(stringr::str_remove(name, 'V'))) %>%
    select(-name)
  
  return( list(df = df, rho = rho) )
}


# Peaks:
calc_peaks <- function(df, anchor.date) {
  scenario =  df$scenario[1]
  df.pk = NULL
  try({
    df.pk <- df %>%
      group_by(iter) %>%
      summarize(pkv = max(value),
                pkt = t[which.max(value)]) %>%
      mutate(pk.date = as_date(anchor.date + pkt))%>%
      mutate(Province = df$Province[1]) %>%
      mutate(scenario = scenario) 
  })
  return(df.pk)
}

# Total epidemic cumul:
calc_cumul <- function(df) {
  scenario =  df$scenario[1]
  df.cum <- df %>%
    group_by(iter) %>%
    summarise(s = sum(value))%>%
    mutate(Province = df$Province[1]) %>%
    mutate(scenario = scenario) 
  return(df.cum)
}

# Sumary stats trajectory
calc_ss <- function(df, CI, anchor.date) {
  qt <- c(rev(0.5-CI/2), 0.5+CI/2)
  scenario =  df$scenario[1]
  
  # For log-scale cosmetics:
  eps <- 0.1
  
  dfs <- df %>%
    group_by(t) %>%
    summarise(m    = max(eps,mean(value, na.rm = TRUE)),
              md   = max(eps,median(value, na.rm = TRUE)),
              qvlo = max(eps,quantile(value, probs = qt[1], na.rm = TRUE)),
              qlo  = max(eps,quantile(value, probs = qt[2], na.rm = TRUE)),
              qhi  = max(eps,quantile(value, probs = qt[3], na.rm = TRUE)),
              qvhi = max(eps,quantile(value, probs = qt[4], na.rm = TRUE))) %>%
    mutate(Province = df$Province[1]) %>%
    mutate(scenario = scenario) %>%
    mutate(date = lubridate::as_date(anchor.date + t))
  
  return(dfs)
}


calc_duration <- function(df, np) {
  x <- df %>%
    filter(0.1 < value & value < 5 & t > np) %>%
    group_by(iter) %>%
    summarize(last.time = min(t))
  return(x)
}

# Distribution incidence `w` days ahead
distrib_future_inc <- function(df, w, positives){
  
  # w =7
  z <- filter(df, t >np & t <= np+w) %>%
    group_by(iter) %>%
    summarise(cum.future = sum(value))
  return(z)
}


proba_inc_above_past <- function(df, w, positives, mult) {
  
  z <- distrib_future_inc(df, w, positives)
  
  # Calculate past incidence through past window:
  np = length(positives)
  cum.past <- sum(positives[(np-w):np])
  
  x = (z$cum.future > mult * cum.past)
  mean(x)
}


proba_peak_timing <- function(df.pk, 
                              anchor.date, 
                              before.after, 
                              target.date ) {
  
  dt = as.numeric( target.date - anchor.date )
  
  if(before.after=='before') x <- df.pk$pkt <= dt
  if(before.after=='after')  x <- df.pk$pkt >= dt
  p = mean(x)  
  return(p)
}


simulate_fwd <- function(prms) {
  
  # --- Unpack:
  
  prov  <- prms$prov
  n.mc  <- prms$n.mc
  horiz <- prms$horiz
  
  pop.full <- prms$pop.full
  
  prop.susceptible <- prms$prop.susceptible
  
  # Severity Cascade:  
  prop.true.pos.beta1 = prms$prop.true.pos.beta1
  prop.true.pos.beta2 = prms$prop.true.pos.beta2
  prop.pos.hosp       = prms$prop.pos.hosp
  prop.hosp.critical  = prms$prop.hosp.critical
  prop.critical.death = prms$prop.critical.death
  
  positives = prms$positives
  R0.m      = prms$R0.m
  R0.sd     = prms$R0.sd
  
  scenario = prms$scenario
  
  # --- Simulate:
  
  pop.eff = pop.full * prop.susceptible
  
  sims <- sim_incidence(n.mc          = n.mc, 
                        horiz         = horiz,
                        pop.eff       = pop.eff, 
                        R0.m          = R0.m, 
                        R0.sd         = R0.sd,
                        scenario      = scenario,
                        positives     = positives,
                        prop.true.pos.beta1 = prop.true.pos.beta1,
                        prop.true.pos.beta2 = prop.true.pos.beta2) 
  
  df <- sims[['df']] %>%
    mutate(Province = prov) %>%
    mutate(scenario = scenario) %>%
    mutate(date = lubridate::as_date(prms$anchor.date + t))
  
  rho <- sims[['rho']]
  
  df.pos       <- convert_true_pos(df, rho = rho, cv =0.2) #TODO: read from file
  df.hosp      <- convert_pos_hosp(df.pos, h = prop.pos.hosp)
  df.critical  <- convert_hosp_critical(df.hosp, kappa = prop.hosp.critical)
  df.death     <- convert_critical_death(df.critical, 
                                         d = prop.critical.death,
                                         )
  
  return(list(
    inc     = df,
    testpos = df.pos,
    hosp    = df.hosp,
    critical= df.critical,
    death   = df.death
  ))
  
  
}


digest_sims <- function(sims, anchor.date, CI = c(0.50, 0.95)) {
  res <- list(
    ss    = lapply(sims, calc_ss, CI=CI, anchor.date = anchor.date),
    peaks = lapply(sims, calc_peaks, anchor.date = anchor.date),
    cum   = lapply(sims, calc_cumul)
  )
  return(res)
}


forecast_unit <- function(prms) {
  sims <- simulate_fwd(prms)
  
  ds <- digest_sims(sims, 
                    anchor.date = prms$anchor.date, 
                    CI = prms$CI)
  return(list(sims = sims, 
              digest = ds))
}




