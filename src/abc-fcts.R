library(snowfall)
source('fcst-fcts.R')
source('resude.R')
source('behaviour.R')

# Draw ABC priors
draw_priors <- function(prms.prior) {
    
    n = prms.prior$n.abc
    
    rho = rbeta(n = n, 
                shape1 = prms.prior$prop.true.pos.beta1, 
                shape2 = prms.prior$prop.true.pos.beta2)
    
    R0 = rnorm(n = n,
               mean = prms.prior$R0.mean,
               sd   = prms.prior$R0.sd)
    R0[R0 <= 0] <- 0.1
    
    alpha = runif(n = n, 
                  min = prms.prior$alpha.min, 
                  max = prms.prior$alpha.max)
    
    lambda = runif(n = n, 
                   min = prms.prior$lambda.min, 
                   max = prms.prior$lambda.max)
    
    prop.crit.death = runif(n = n, 
                            min = prms.prior$p.crit.death.min, 
                            max = prms.prior$p.crit.death.max)

    prop.crit.death.adj = runif(n = n, 
                                min = prms.prior$adj.death.min, 
                                max = prms.prior$adj.death.max)
    
    return(list(rho    = rho, 
                alpha  = alpha, 
                R0     = R0, 
                lambda = lambda,
                prop.crit.death     = prop.crit.death,
                prop.crit.death.adj = prop.crit.death.adj))
}

#' Distance between prior simulations and observations. Used for ABC.
calc_distance <- function(sim, obs, t.fit.start, t.fit.end, 
                          recency) {
    
    # Range upon wich we fit:
    rng = t.fit.start:t.fit.end
    
    # Weights:
    w = 1
    if(recency>0) w = rng^recency
    
    # Distance calculation:
    tmp = (sim[rng] - obs[rng])^2 * w / sum(w)
    d   = sqrt(sum(tmp)) 
    return(d)
}


#' Implement the severity cascade.
#' Convert true incidence to reports, hospitalizations, critical, deaths.
severity_cascade <- function(df.inc, prms, rho, 
                             prop.crit.death     = NULL, 
                             prop.crit.death.adj = NULL,
                             i = NULL) # USED???
{
    df.pos       <- convert_true_pos(df.inc, 
                                     rho = rho, 
                                     cv = 0.2) #TODO: read from file
    df.hosp      <- convert_pos_hosp(df.pos, 
                                     h = prms$prop.pos.hosp)
    df.critical  <- convert_hosp_critical(df.hosp, 
                                          kappa = prms$prop.hosp.critical)

    # TODO: re-write when sure...    
    d = prms$prop.critical.death
    death.adj = 0
    if(!is.null(prop.crit.death)){
        d = prop.crit.death #priors$prop.crit.death[i]
        death.adj = prop.crit.death.adj #priors$prop.crit.death.adj[i]
    }
    
    df.death <- convert_critical_death(df = df.critical, d, 
                                       death.adj = death.adj)
    
    return(list(
        df.inc      = df.inc,
        df.pos      = df.pos,
        df.hosp     = df.hosp,
        df.critical = df.critical,
        df.death    = df.death
    ))
}

#' Apply lag times to all outcomes (hosp, deaths, etc.)
shift_times_all <- function(dflist, prms) {
    
    L.inf.pos    = prms$lag.inf.pos
    L.pos.hosp   = prms$lag.pos.hosp
    L.hosp.crit  = prms$lag.hosp.critical
    L.crit.death = prms$lag.critical.death
    
    df          = shift_times(df = dflist$df.inc, shift = -L.inf.pos)
    df.hosp     = shift_times(dflist$df.hosp, L.pos.hosp)
    df.critical = shift_times(dflist$df.critical, L.pos.hosp + L.hosp.crit)
    df.death    = shift_times(dflist$df.death, L.pos.hosp + L.hosp.crit + L.crit.death)
    
    return(list(
        inc      = df,
        testpos  = dflist$df.pos,
        hosp     = df.hosp,
        critical = df.critical,
        death    = df.death
    ))
}



#' Unit function for parallel computing.
#' @param i Integer. The index of the ABC iterations.
#' @param prms List of simulation and population parameters.
#' @param priors List of vectors for the sampled priors. Name is the name of the parameter.
#' @param gi List of parameter for the generation interval distribution
#' @param inc.init Matrix of initial incidence. Row is iteration, column is time.
#' @param prms.fit List of parameters used for fitting
#' @param B.fct Function. Temporal behaviour change.
#' @param B.prm List of parameters used by B.fct. 
#' 
abc_unit <- function(i, 
                     prms, 
                     priors, 
                     gi, 
                     inc.init, 
                     prms.fit, 
                     B.fct, 
                     B.prm) 
{
    pop.eff  = prms$pop.full * prms$prop.susceptible
    pop_size = pop.eff - sum(inc.init[i,])
    
    B.prm[['m1']] <- priors$lambda[i]
    
    rng.inc.init = 1:(prms.fit$t.fit.start-1)
    
    sim = resude_simulate_b(pop_size = pop_size, 
                            I.init  = inc.init[i, rng.inc.init], 
                            R0      = priors$R0[i], 
                            alpha   = priors$alpha[i], 
                            GI_span = round(gi$mean + 2*gi$variance), 
                            GI_mean = gi$mean, 
                            GI_var  = gi$variance, 
                            horizon = prms$horiz, 
                            distrib.proc = 'rgamma', 
                            phi     = gi$gam.disp ,
                            B.fct   = B.fct, 
                            B.prm   = B.prm,
                            seed    = i)
    
    tt = 1:length(sim$I)
    sim.inc = data.frame(value = sim$I, 
                         iter  = i,
                         t     = tt,
                         date  = lubridate::as_date(prms$anchor.date)+tt)
    
    # DELETE WHEN SURE:
    # From true latent incidence to positive tests:
    # Note: data frame weird, but this is to use the existing function... 
    # sim.pos = convert_true_pos(df  = sim.inc, 
    #                            rho = priors$rho, 
    #                            cv  = 0.2) #TODO: read from file...
    # 
    
    # TODO : make %>% exportable in snowfall with sFLibrary(dplyr)
    tmp = severity_cascade(df.inc = sim.inc, 
                           prms, 
                           rho = priors$rho, 
                           prop.crit.death = priors$prop.crit.death[i],
                           prop.crit.death.adj = priors$prop.crit.death.adj[i],
                           i=i)
    z = shift_times_all(dflist = tmp, prms)
    
    sim.testpos  = z$testpos
    sim.hosp     = z$hosp
    sim.critical = z$critical
    sim.death    = z$death
    
    # This should work!
    #
    # z = sim.inc %>%
    #     severity_cascade(prms, 
    #                      rho = priors$rho, 
    #                      priors = priors, 
    #                      i=i) %>%
    #     shift_times_all(prms)
    
    
    d.pos  = calc_distance(sim = z$testpos$value, # sim.pos$value, 
                           obs = prms$positives, 
                           t.fit.start = prms.fit$t.fit.start, 
                           t.fit.end   = prms.fit$t.fit.end, 
                           recency     = prms.fit$recency)
    
    d.death = calc_distance(sim = z$death$value , 
                            obs = prms$deaths, 
                            t.fit.start = prms.fit$t.fit.start, 
                            t.fit.end   = length(prms$deaths), 
                            recency     = prms.fit$recency)
    
    d = d.pos/sum(prms$positives) + d.death/sum(prms$deaths) * prms.fit$weight.death
    
    return(list(d = d, 
                sim.inc = sim.inc,
                sim.testpos  = sim.testpos,
                sim.hosp     = sim.hosp,
                sim.critical = sim.critical,
                sim.death = sim.death))
}

#' Returns ABC posterior
abc_posteriors <- function(d.abc, priors, p.abc) {
    thresh = quantile(d.abc, probs = p.abc, na.rm = TRUE)
    idx = which(d.abc <=thresh)
    
    return(list(
        idx    = idx,
        R0     = priors$R0[idx],
        alpha  = priors$alpha[idx],
        rho    = priors$rho[idx],
        lambda = priors$lambda[idx],
        prop.crit.death     = priors$prop.crit.death[idx],
        prop.crit.death.adj = priors$prop.crit.death.adj[idx]
    ))
}

#' Perform the full ABC fit.
fit_abc <- function(prms.fit, prms.prior, prms, B.fct, B.prm) {
    
    gi       = get_GI_prms()
    priors   = draw_priors(prms.prior)
    # Correct data:
    positives = prms$positives
    positives[positives<0] <- 0
    
    inc.init = positives_to_true_incidence(rho = priors$rho, 
                                           positives = positives)
    
    sfInit(parallel = prms.fit$n.cpus>1, 
           cpus = prms.fit$n.cpus)
    sfExportAll()
    x = sfLapply(1:prms.prior$n.abc, 
                 fun       = abc_unit,
                 prms      = prms, 
                 priors    = priors,
                 gi        = gi, 
                 inc.init  = inc.init, 
                 prms.fit  = prms.fit,
                 B.fct     = B.fct, 
                 B.prm     = B.prm)
    sfStop()
    
    d.abc         = sapply(x,'[[',1)
    sim.inc.abc   = lapply(x,'[[','sim.inc') %>% lapply('[[',1)
    sim.testpos.abc  = lapply(x,'[[','sim.testpos') %>% lapply('[[',1)
    sim.hosp.abc     = lapply(x,'[[','sim.hosp') %>% lapply('[[',1)
    sim.critical.abc = lapply(x,'[[','sim.critical') %>% lapply('[[',1)
    sim.death.abc    = lapply(x,'[[','sim.death') %>% lapply('[[',1)
    
    posteriors  = abc_posteriors(d.abc, priors, p.abc = prms.fit$p.abc)
    
    post.inc      = sim.inc.abc[posteriors$idx]
    post.testpos  = sim.testpos.abc[posteriors$idx]
    post.hosp     = sim.hosp.abc[posteriors$idx]
    post.critical = sim.critical.abc[posteriors$idx]
    post.death    = sim.death.abc[posteriors$idx]
    
    return(list(posteriors = posteriors, 
                post.inc = post.inc,
                post.testpos = post.testpos,
                post.hosp = post.hosp,
                post.critical = post.critical,
                post.death = post.death,
                priors = priors) )
}


plot_post_abc <- function(abcfit) {
    
    priors     = abcfit$priors
    posteriors = abcfit$posteriors
    
    df1 = data.frame(priors)
    df1$type = 'prior'
    df2 = data.frame(posteriors) %>%
        select(-idx)
    df2$type = 'posterior'
    df = rbind(df1, df2)
    
    q = df %>% 
        pivot_longer(cols = which(names(df)!='type'), 
                     names_to = 'param', values_to = 'value')
    g = q %>%
        ggplot(aes(x=value)) + 
        geom_histogram(aes(fill = type, y = ..density..), 
                       alpha=0.5, bins=30, position = 'identity') + 
        scale_fill_manual(values=c('black','grey80'))+
        facet_wrap(~param, scales = 'free')+
        theme(panel.grid = element_blank(), 
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())+
        xlab('') + ylab('')
    plot(g)
}

plot_post_traj_pos <- function(fcst, prms.fit, B.prm, positives,...) {
    
    d = fcst$digest$ss$testpos
    
    df.pos = data.frame(x = first.date+1:length(positives),
                        y = positives)
    
    fit.rng = first.date + c(prms.fit$t.fit.start, prms.fit$t.fit.end)
    
    g = d %>%
        ggplot()+
        # Fit infos:
        geom_vline(xintercept = fit.rng,
                   linetype = 'dashed')+
        geom_line(aes(x=date, y=m))+
        geom_ribbon(aes(x=date, ymin=qvlo, ymax=qvhi), alpha=0.2)+
        geom_ribbon(aes(x=date, ymin=qlo, ymax=qhi), alpha=0.2)+
        geom_step(data =df.pos, aes(x=x,y=y),colour='red', size=0.3)+
        geom_point(data =df.pos, aes(x=x,y=y), pch=21, 
                   stroke=1, colour='red2', fill='white')+
        # Behaviour markers:
        geom_vline(xintercept = first.date + B.prm$t1, linetype='dotted',
                   colour='blue')+
        scale_x_date(date_breaks = '1 month', date_labels = '%d-%b-%y')+
        scale_y_log10(breaks = c(1,10,100,1000,10000)) +
        coord_cartesian(ylim=c(1,max(d$qvhi)*1.5))+
        theme(text = element_text(size=16),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(angle=30,hjust=1))+
        ggtitle(paste('ABC fit of Positive tests incidence:',prov,'/',scenario))+
        xlab('')+ylab('')
    plot(g)
}

plot_post_traj_death <- function(fcst, prms.fit, B.prm, deaths,...) {
    
    d = fcst$digest$ss$death
    
    df.death = data.frame(x = first.date+1:length(deaths),
                          y = deaths) %>%
        mutate(cumdeath = cumsum(y))
    
    # incidence plot
    
    g = d %>%
        ggplot()+
        geom_line(aes(x=date, y=m))+
        geom_ribbon(aes(x=date, ymin=qvlo, ymax=qvhi), alpha=0.2)+
        geom_ribbon(aes(x=date, ymin=qlo, ymax=qhi), alpha=0.2)+
        geom_point(data =df.death, aes(x=x,y=y), pch=23, 
                   stroke=1, colour='red2', fill='white')+
        # geom_step(data =df.death, aes(x=x,y=y),colour='red', size=0.2, alpha=0.5)+
        # Behaviour markers:
        geom_vline(xintercept = first.date + B.prm$t1, linetype='dotted',
                   colour='blue')+
        scale_x_date(date_breaks = '1 month', date_labels = '%d-%b-%y')+
        scale_y_log10(breaks = c(1,10,seq(20,100,20),seq(200,1000,200))) +
        coord_cartesian(ylim=c(1,1000))+
        theme(text = element_text(size=12),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(angle=30,hjust=1),
              axis.text.y = element_text(size=10))+
        ggtitle(paste('ABC fit of Deaths incidence:',prov,'/',scenario))+
        xlab('')+ylab('')
    
    # cumulative plot 
    
    dc = d %>%
        mutate(cumdeath.m = cumsum(m),
               cumdeath.vlo = cumsum(qvlo),
               cumdeath.vhi = cumsum(qvhi)
        ) %>%
        select(date, starts_with('cum'))
    
    a = 0.2
    gc = dc %>% ggplot() +
        geom_line(aes(x=date, y=cumdeath.m)) + 
        geom_ribbon(aes(x=date, ymin=cumdeath.vlo, ymax=cumdeath.vhi), alpha=a)+
        geom_line(data=df.death, aes(x=x, y=cumdeath), colour='red', size=2, alpha=0.25)+
        geom_point(data=df.death, aes(x=x, y=cumdeath), colour='red', size=1, alpha=0.75)+
        # Behaviour markers:
        geom_vline(xintercept = first.date + B.prm$t1, linetype='dotted',
                   colour='blue')+
        scale_x_date(date_breaks = '1 month', date_labels = '%d-%b-%y')+
        scale_y_log10() +
        coord_cartesian(ylim=c(1,max(dc$cumdeath.vhi)*1.5))+
        theme(text = element_text(size=12),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(angle=30,hjust=1),
              axis.text.y = element_text(size=10))+
        ggtitle(paste('ABC fit of Deaths incidence:',prov,'/',scenario))+
        xlab('')+ylab('')
    
    gridExtra::grid.arrange(g,gc, ncol=2)
}



glue_list_df <- function(x, prms) {
    tmp = list()
    
    for(i in 1:length(x)){
        tt = 1:length(x[[i]])
        tmp[[i]] <- data.frame(iter = i, 
                               value = x[[i]], 
                               t = tt,
                               date = lubridate::as_date(prms$anchor.date) + tt,
                               Province = prms$prov,
                               scenario = prms$scenario)
    }
    df = do.call('rbind', tmp)
    df$Province <- as.character(df$Province)
    df$scenario <- as.character(df$scenario)
    return(df)
}


#' Reformat ABC output into dataframes
reformat_post <- function(abcfit, prms) {
    
    np = c( "post.inc", "post.testpos",
            "post.hosp", "post.critical", "post.death")
    
    a = lapply(X = abcfit[np], FUN = glue_list_df, prms = prms)
    names(a) <- stringr::str_replace(np,'post.','')
    return(a)
}


shift_times <- function(df, shift) {
    # shift = -5
    df$t = df$t + shift
    df$date = df$date + shift
    return(df)
}

#' Wrap function to calculate posteriors
#' incidence trajectories using ABC.
#' TODO: re-factor the arguments
#' 
forecast_unit_abc <- function(prms, 
                              prms.fit, 
                              prms.prior, 
                              B.fct, 
                              B.prm) {
    # Run the ABC fit:
    abcfit = fit_abc(prms.fit, 
                     prms.prior, 
                     prms, 
                     B.fct = B.fct,
                     B.prm = B.prm)
    
    # Reformat the posterior trajectories:
    sims = abcfit %>%
        reformat_post(prms = prms)
    
    ds = digest_sims(sims, 
                     anchor.date = prms$anchor.date, 
                     CI = prms$CI)
    
    return(list(sims   = sims, 
                digest = ds, 
                abcfit = abcfit))
}



save_lambda_post <- function(fcst, prov, scenario) {
    
    x = fcst$abcfit$posteriors$lambda
    write.table(round(100 * (1-mean(x))), 
                file = paste0('out/lambda-mean-',prov,'-',scenario ,'.txt'),
                quote = F, row.names = F, col.names = F, eol = '')
    write.table(round(100 * (1-quantile(x,probs = 0.025))), 
                file = paste0('out/lambda-hi-',prov,'-',scenario ,'.txt'),
                quote = F, row.names = F, col.names = F, eol = '')
    write.table(round(100 * (1-quantile(x,probs = 0.975))), 
                file = paste0('out/lambda-lo-',prov,'-',scenario ,'.txt'),
                quote = F, row.names = F, col.names = F, eol = '')
}
