# ---- Params ----

first.date  <- lubridate::as_date("2020-03-06")

dat       = get_data_testpos(first.date = first.date)
positives = dat$value
deaths    = get_data_deaths(first.date)$value

# plot(deaths/positives,typ='b')

last.obs.date <- max(dat$Date)

anchor.date <- dat$Date[1] - dat$time[1]
R0estim     <- estim_R0_simple(prov)
pop.full    <- get_prm(prov, csv.name = '../data/populations.csv')

file.prms = 'prm-model.csv'

n.mc                <- get_prm('n.mc', file.prms)
horiz               <- get_prm('horizon', file.prms)
prop.susceptible    <- get_prm('prop.susceptible', file.prms)

prop.true.pos.beta1  <- get_prm('prop.true.pos.beta1', file.prms)
prop.true.pos.beta2  <- get_prm('prop.true.pos.beta2', file.prms)

# ---  Severity cascade:
prop.pos.hosp        <- get_prm('prop.pos.hosp', file.prms)
prop.hosp.critical   <- get_prm('prop.hosp.critical', file.prms)
prop.critical.death  <- get_prm('prop.critical.death', file.prms)

# ---- Lags:
file.lags = 'prm-lags.csv'
lag.inf.pos        = get_prm('lag_infection_to_report', file.lags)
lag.pos.hosp       = get_prm('lag_report_to_hosp', file.lags)
lag.hosp.critical  = get_prm('lag_hosp_to_critical', file.lags)
lag.critical.death = get_prm('lag_critical_to_death', file.lags)

# ----
prms <- list(
    prov                = prov,
    scenario            = scenario,
    n.mc                = n.mc,
    horiz               = horiz,
    pop.full            = pop.full,
    prop.susceptible    = prop.susceptible,
    prop.true.pos.beta1 = prop.true.pos.beta1,
    prop.true.pos.beta2 = prop.true.pos.beta2,
    prop.pos.hosp       = prop.pos.hosp, 
    prop.hosp.critical  = prop.hosp.critical, 
    prop.critical.death = prop.critical.death,
    positives           = positives,
    deaths              = deaths,
    R0.m                = R0estim$mean,
    R0.sd               = R0estim$sd,
    anchor.date         = first.date,
    t.kappa             = 20, # DELETE WHEN SURE
    CI                  = c(0.50, 0.95),
    lag.inf.pos         = lag.inf.pos,
    lag.pos.hosp        = lag.pos.hosp,
    lag.hosp.critical   = lag.hosp.critical,
    lag.critical.death  = lag.critical.death
)

# === ABC priors ====

cpumax = parallel::detectCores()

file.abc   = 'prm-abc.csv'

n.abc      = get_prm('n.abc', file.abc)
R0.mean    = get_prm('prior.R0.mean', file.abc)
R0.sd      = get_prm('prior.R0.sd', file.abc)
alpha.min  = get_prm('prior.alpha.min', file.abc)
alpha.max  = get_prm('prior.alpha.max', file.abc)
lambda.min = get_prm('prior.lambda.min', file.abc)
lambda.max = get_prm('prior.lambda.max', file.abc)
recency    = get_prm('recency', file.abc)
ncpus      = get_prm('ncpus', file.abc)

p.crit.death.min = get_prm('prior.prop.crit.death.min', file.abc)
p.crit.death.max = get_prm('prior.prop.crit.death.max', file.abc)
adj.death.min    = get_prm('adj.death.min', file.abc)
adj.death.max    = get_prm('adj.death.max', file.abc)

weight.death = get_prm('weight.death', file.abc)

prms.prior = list(n.abc = n.abc,
                  prop.true.pos.beta1 = prop.true.pos.beta1,
                  prop.true.pos.beta2 = prop.true.pos.beta2, 
                  R0.mean    = R0.mean,
                  R0.sd      = R0.sd,
                  alpha.min  = alpha.min,
                  alpha.max  = alpha.max,
                  lambda.min = lambda.min,
                  lambda.max = lambda.max, 
                  p.crit.death.min = p.crit.death.min,
                  p.crit.death.max = p.crit.death.max, 
                  adj.death.min = adj.death.min,
                  adj.death.max = adj.death.max)

#TODO: do something smarter with t.fit.start/end
t.fit.start = 12 # When data starts to be not too noisy
t.fit.end   = ifelse(scenario == 'BASELINE', 23, length(positives))

prms.fit = list(recency     = recency,
                t.fit.start = t.fit.start, 
                t.fit.end   = t.fit.end,
                p.abc       = ifelse(n.abc >= 5e3, 0.02, 30/n.abc),
                n.cpus      = ifelse(ncpus>0, ncpus, cpumax-ncpus),
                weight.death = weight.death)





