

GI_dist <- function(t, GI_span, GI_mean, GI_var){
  tvec <- 0:GI_span
  
  GI_k     <- GI_mean^2/GI_var
  GI_scale <- GI_mean/GI_k
  
  tmp  <- tvec^(GI_k-1) * exp(-tvec/GI_scale)
  tmp2 <- t^(GI_k-1) * exp(-t/GI_scale)
  return(tmp2/sum(tmp))
}


indic_interval <- function(t,a,b) {
  res = 0
  if(a <= t & t <= b) res = 1
  return(res)
}


resude_simulate <- function(pop_size, 
                            I.init,
                            R0, 
                            alpha, 
                            kappa, 
                            GI_span, 
                            GI_mean, 
                            GI_var,
                            horizon,
                            distrib.proc,
                            phi = NULL,
                            t.kappa = 0,   # When kappa kicks in
                            lambda = 0,
                            pulse.height = 0,
                            pulse.length = 0,
                            pulse.quiet.length = 0,
                            seed=123) {
  set.seed(seed)
  
  d.found <- FALSE
  if(distrib.proc %in% c('rpois', 'rgamma')) d.found <- TRUE
  if(!d.found){
    message(paste('ERROR: distribution process',
                  distrib.proc,'not implemented. Aborting!'))
  }
  if(distrib.proc=='rgamma' & is.null(phi)){
    message('ERROR: distribution is Gamma but parameter phi=NULL. Aborting!')
    stop()
  }
  
  # Typically used to generate synthetic data:
  if(length(I.init)==1){
    I <- vector()
    S <- vector()
    I[1] <- I.init
    S[1] <- pop_size - I.init
    numobs <- 1
  }
  
  # Typically used when forecasting 
  # an already-started epidemic:
  if(length(I.init)>1){
    numobs <- length(I.init)
    I <- I.init
    S <- pop_size - cumsum(I.init)
  }
  
  k.pulse <- 0
  if(pulse.length>0) k.pulse = round(horizon / (pulse.length + pulse.quiet.length))+1
  
  # Simulate forward:
  for(t in (numobs+1):(numobs+horizon)){
    
    # GI distribution:
    z <- 0
    for(j in 1:min(GI_span,t-1)){
      GI_j <- GI_dist(j, GI_span, GI_mean, GI_var)
      # print(GI_j)
      z <- z + GI_j * I[t-j]
    }
    
    # pulse multiplier:
    ps = 0
    if(k.pulse>0){
      for(k in 1:k.pulse){
        ps = ps + indic_interval(t, 
                                 (k+1)*pulse.quiet.length + k*pulse.length,
                                 (k+1)*(pulse.quiet.length + pulse.length))
      }
    }
    mult.pulse    = 1.0 + pulse.height * ps
    
    tk = max(t - t.kappa, 0)
    mult.decrease = max( lambda, exp(-kappa*tk) )
    mult          = mult.decrease * mult.pulse
    
    # print(paste("DEBUG: t = ",t,"mult = ",mult))
    
    # Expected (mean) incidence:    
    I.tmp <- (S[t-1]/ pop_size)^(1+alpha) * R0 * mult * z 
    
    # Add uncertainty in transmission process with
    # choice between:
    # - Poisson distribution (~Normal when large, hence symetrical & thin tails)
    # - Gamma (skewed, fat for smaller values but also thicker tail for large values)
    if(distrib.proc=='rpois')
    {
      I[t] <- rpois(n = 1, 
                    lambda =  min(I.tmp, S[t-1]) )
    }
    
    if(distrib.proc=='rgamma')
    {
      m <- min(I.tmp, S[t-1])
      I[t] <- rgamma(n = 1, 
                     shape = m / (1+phi), 
                     rate  = 1/(1+phi)  )
    }
    
    # Update susceptible number:
    S[t] <- max(0, S[t-1] - I[t])
  }
  return(list(S=S, I=I))
}


RESuDe.simulate.wrap <- function(pop_size, 
                                 I.init,
                                 R0, 
                                 alpha, 
                                 kappa, 
                                 GI_span, 
                                 GI_mean, 
                                 GI_var,
                                 horizon,
                                 last.obs,
                                 do.plot = FALSE,
                                 distrib.proc,
                                 seed=123){
  # Generate full epidemic
  syn.data <- resude_simulate(pop_size = pop_size, 
                              I.init = I.init,
                              R0 = R0,
                              alpha = alpha,
                              kappa = kappa, 
                              GI_span = GI_span, 
                              GI_mean = GI_mean, 
                              GI_var = GI_var,
                              horizon = horizon,
                              distrib.proc = distrib.proc,
                              seed = seed)
  syn.inc.full <- syn.data$I
  
  # Just take the start of epidemic
  # (will forecast the rest)
  syn.inc <- syn.inc.full[1:last.obs]
  numobs <- length(syn.inc)
  
  if(do.plot){
    # Generation interval distribution:
    par(mfrow=c(1,1))
    tt <- 0:GI_span
    plot(tt,GI_dist(t = tt,GI_span = GI_span,
                    GI_mean = GI_mean, 
                    GI_var = GI_var),
         typ='l',lwd=6,col='blue',xlab="",ylab="",main="GI distrib")
    grid()
    
    # Incidence:
    par(mfrow=c(1,2))
    plot(syn.inc.full,typ='s', lwd=1, main="Incidence data",
         xlab="time",ylab="")
    lines(syn.inc,typ='s', lwd=6)
    plot(syn.inc.full,typ='o',log="y", lwd=1, xlab="time",ylab="",cex=0.5)
    lines(syn.inc,typ='o', lwd=6) ; grid()
  }
  
  return(list(syn.inc.full = syn.inc.full,
              syn.inc = syn.inc))
}

#' Create aggregated periodic observations
#' @param dat Full data (e.g. daily)
#' @param prd Period of aggregation (in days)
#' @param report.proba Reporting probability for the binomial observation process.
#' @param last.obs Last observation date.
#' 
create.aggreg.prd.obs <- function(dat, 
                                  prd, 
                                  distrib.obs,
                                  report.proba,  
                                  last.obs) {
  
  inc.full <- data.frame(t = 1:length(dat),
                         inc = dat)
  # period of observations:
  t.obs <- numeric()
  sip   <- numeric() # <-- sum incidence during period
  
  # Calculate the sum of incidence
  # during each period:
  cnt <- 1
  for(i in 1:nrow(inc.full)){
    if(inc.full$t[i] %% prd == 0){
      tmp        <- sum(inc.full$inc[(i-prd+1):i])
      sip[cnt]   <- tmp
      t.obs[cnt] <- inc.full$t[i]
      cnt <- cnt + 1
    }
  }
  # Aggregated incidence:
  dat.aggr <- data.frame(t = t.obs, inc = sip)
  
  # Add observation error:
  d.found <- FALSE
  sampled.inc <- NULL
  
  if(distrib.obs=='rbinom'){
    d.found <- TRUE
    sampled.inc <- rbinom(n = length(t.obs), 
                          size = round(sip), 
                          prob = report.proba)
  }
  
  if(distrib.obs=='rnorm'){
    # Normal approx to binomial
    d.found <- TRUE
    sampled.inc <- rnorm(n = length(t.obs), 
                         mean = report.proba * sip, 
                         sd = sqrt(report.proba (1-report.proba)* sip))
    sampled.inc[sampled.inc<0] <- 0
  }
  
  if(!d.found){
    message(paste('ERROR: observation distribution [',
                  distrib.obs,
                  '] not implemented!'))
    stop()
  }
  
  dat.aggr.obs <- data.frame(t   = t.obs, 
                             inc = sampled.inc )
  
  # Crop observation up to last observation time:
  dat.aggr.obs.past <- subset(dat.aggr.obs, t <= last.obs)
  
  return(list(inc.full    = inc.full, 
              dat.aggr = dat.aggr, 
              dat.aggr.obs = dat.aggr.obs,
              dat.aggr.obs.past = dat.aggr.obs.past))
  # dat.prd.obs = dat.prd.obs))
}




resude_simulate_b <- function(pop_size, 
                              I.init,
                              R0, 
                              alpha, 
                              GI_span, 
                              GI_mean, 
                              GI_var,
                              horizon,
                              distrib.proc,
                              phi = NULL,
                              B.fct = NULL,
                              B.prm = NULL,
                              seed=123) {
  set.seed(seed)
  
  d.found <- FALSE
  if(distrib.proc %in% c('rpois', 'rgamma')) d.found <- TRUE
  if(!d.found){
    message(paste('ERROR: distribution process',
                  distrib.proc,'not implemented. Aborting!'))
  }
  if(distrib.proc=='rgamma' & is.null(phi)){
    message('ERROR: distribution is Gamma but parameter phi=NULL. Aborting!')
    stop()
  }
  
  # Typically used to generate synthetic data:
  if(length(I.init)==1){
    I <- vector()
    S <- vector()
    I[1] <- I.init
    S[1] <- pop_size - I.init
    numobs <- 1
  }
  
  # Typically used when forecasting 
  # an already-started epidemic:
  if(length(I.init)>1){
    numobs <- length(I.init)
    I <- I.init
    S <- pop_size - cumsum(I.init)
  }
  
  # Simulate forward:
  for(t in (numobs+1):(numobs+horizon)){
    
    # GI distribution:
    z <- 0
    for(j in 1:min(GI_span,t-1)){
      GI_j <- GI_dist(j, GI_span, GI_mean, GI_var)
      # print(GI_j)
      z <- z + GI_j * I[t-j]
    }
    
    # message("DEBUG:")
    # message(paste(B.prm, collapse = ' ; '))
    
    # Behaviour changes:
    Bt = B.fct(t, B.prm)
    
    # Expected (mean) incidence:    
    I.tmp <- (S[t-1]/ pop_size)^(1+alpha) * R0 * Bt * z 
    
    # Add uncertainty in transmission process with
    # choice between:
    # - Poisson distribution (~Normal when large, hence symetrical & thin tails)
    # - Gamma (skewed, fat for smaller values but also thicker tail for large values)
    if(distrib.proc=='rpois')
    {
      I[t] <- rpois(n = 1, 
                    lambda =  min(I.tmp, S[t-1]) )
    }
    
    if(distrib.proc=='rgamma')
    {
      m <- min(I.tmp, S[t-1])
      I[t] <- rgamma(n = 1, 
                     shape = m / (1+phi), 
                     rate  = 1/(1+phi)  )
    }
    
    # Update susceptible number:
    S[t] <- max(0, S[t-1] - I[t])
  }
  return(list(S=S, I=I))
}

 # ---- debug---- 

if(0){  #DEBUG
  
  
  prmb = list(t1=60, m1=0.50)
  
  
  sim = resude_simulate(pop_size = 1e6, 
                           I.init = c(12,24,55), 
                           R0 = 1.75, 
                           alpha = 1, 
                         kappa = - log(prmb$m1) / 7,
                           GI_span = 14, 
                           GI_mean = 5, 
                           GI_var = 5, 
                           horizon = 180, 
                           t.kappa = prmb$t1,
                        lambda = prmb$m1,
                           distrib.proc = 'rgamma', 
                           phi = 0.7
                          )
  
  simb = resude_simulate_b(pop_size = 1e6, 
                           I.init = c(12,24,55), 
                           R0 = 1.75, 
                           alpha = 1, 
                           GI_span = 14, 
                           GI_mean = 5, 
                           GI_var = 5, 
                           horizon = 180, 
                           distrib.proc = 'rgamma', 
                           phi = 0.7, 
                           B.fct = behave_iso1, 
                           B.prm = prmb)
  
  plot(simb$I, log='y', ylim=c(1,1e5), typ='o')
  lines(sim$I, col='red')
  abline(v=prmb$t1, lty=2)
  grid()
    
}








