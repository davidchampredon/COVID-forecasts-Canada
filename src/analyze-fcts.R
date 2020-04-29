
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2)) 
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(gridExtra))
library(scales)
theme_set(theme_bw())

options(width = 400)
source('utils-simple.R')
source('scenario.R')

dir.plot = 'plots/'
dir.out  = 'out/'


sub_df <- function(prov, df, varname) {
    df$v <- df[[varname]] # common variable name across digested dataframes
    subdf <- filter(df, Province == prov) 
    return(subdf)
}

#' Calculate the relative changes between BASELINE and other scenarios.
calc_chg_baseline <- function(subdf, relative = TRUE) {    # varname = 's'
    
    df <- subdf %>%
        group_by(type, scenario) %>%
        summarise(m = mean(v)) %>%
        pivot_wider(names_from = scenario, values_from = m)%>%
        ungroup() %>%
        as.data.frame()
    
    nm = names(df)
    col.base  <- which(nm =='BASELINE')
    col.intrv <- which(! (nm %in% c('BASELINE', 'type')) ) 
    
    chgname <- ifelse(relative,'relchg','chg')
    
    for(i in seq_along(col.intrv)){  # i=1
        
        if(relative)
            tmp <- df[,col.intrv[i]] / df[,col.base] - 1
        else
            tmp <- df[,col.intrv[i]] - df[,col.base] 
        
        df[[paste0(chgname,i)]] <- tmp
    }
    return(df)
}

calc_proba_above <- function( typ, subdf, thresh) {
    
    # Reminder: `v` is the generic variable name
    
    res <- subdf %>%
        filter( type == typ) %>%
        group_by(type, scenario) %>%
        summarise(p = mean(v > thresh[[typ]]))  %>%
        mutate(threshold = thresh[[typ]])
    return(res)
}


#' Compare intervention outcomes.
compare_intrv <- function(prov, 
                          digest.df,
                          varname,
                          outcometype,
                          thresh,
                          col.scen,
                          ci = c(0.8, 0.95)) {
    
    if(0){
        digest.df = digest.peaks
        varname = 'pkv'
        outcometype = 'Peak Values'
        prov = 'Canada'
    }
    
    subdf <- sub_df(prov, df = digest.df, varname)
    
    df.chg <- calc_chg_baseline(subdf)
    print(df.chg)
    
    scens <- unique(digest.df$scenario)
    
    # NOT SURE TO KEEP... 
    # for(s in scens){
    #     write_output_bignum(val = df[[s]][df$type=='death'],
    #                         name = paste0('cumDeath',s))
    # }
    # 
    
    
    # --- PLOT PROBABILITIES 
    
    title = paste("Proba",outcometype,"above threshold", prov)
    print(title)
    print('threshold: ')
    print(thresh)
    
    
    tmp = lapply(names(thresh), calc_proba_above, subdf=subdf, thresh=thresh)
    df.proba = do.call('rbind', tmp) %>%
        mutate(sdisp = scenario_display(scenario)) 
    
    g.proba <- df.proba %>%
        mutate(typefull = paste('Probability',outcometype,'of', 
                                toupper(type),'>',
                                format(threshold, big.mark=','))) %>%
        # mutate(z = factor(typefull, 
        #                   levels=c(paste('Probability',outcometype,'of', 
        #                                    toupper('hosp'),'>',
        #                                    format(threshold, big.mark=',')),
        #                            
        #                            paste('Probability',outcometype,'of', 
        #                                  toupper('critical'),'>',
        #                                  format(threshold, big.mark=',')),
        #                            
        #                            paste('Probability',outcometype,'of', 
        #                                  toupper('death'),'>',
    #                                  format(threshold, big.mark=','))
    #                            )
    #                   )) %>%
    ggplot() + 
        geom_bar(aes(x = sdisp, y=p, fill = sdisp), 
                 stat = 'identity')+
        facet_wrap(~typefull)+
        xlab('') + ylab('')+
        scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                           limits = c(0,1))+
        scale_fill_manual(values = col.scen)+
        scale_colour_manual(values = col.scen)+
        theme(text = element_text(size = 14),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_text(face='bold'),
              strip.text = element_text(face = 'bold'),
              strip.background = element_rect(fill='grey95'))+
        ggtitle('Probabilistic Outcomes')+
        guides(colour = F, fill = F)
    g.proba
    
    
    # --- PLOT VALUES
    
    title = paste(outcometype, "by Scenario --", prov)
    
    qt <- ci_to_quantiles(ci)
    
    ss <- subdf %>%
        group_by(type, scenario) %>%
        summarize(m = mean(v),
                  q.vlo = quantile(v, probs=qt[1], na.rm = TRUE),
                  q.lo  = quantile(v, probs=qt[2], na.rm = TRUE),
                  q.hi  = quantile(v, probs=qt[3], na.rm = TRUE),
                  q.vhi = quantile(v, probs=qt[4], na.rm = TRUE)
        ) %>%
        mutate(sdisplay = scenario_display(scenario))
    
    a = 0.3
    thick.seg = 7
    g.val <- ss %>%
        filter(type %in% c('hosp', 'critical','death')) %>%
        mutate(typeplot = factor(type, levels=c('hosp', 'critical','death'))) %>%  # OK here but not in proba plot
        ggplot(aes(x = sdisplay)) +
        geom_segment(aes(xend=sdisplay, 
                         y = q.vlo, 
                         yend = q.vhi,
                         colour = sdisplay), 
                     size = thick.seg, alpha = a)+
        geom_segment(aes(xend=sdisplay, 
                         y = q.lo, 
                         yend = q.hi,
                         colour = sdisplay), 
                     size = thick.seg, alpha = a)+
        geom_point(aes(y=m, colour = sdisplay), size= thick.seg-1)+
        scale_y_continuous(label=comma)+
        facet_wrap(~type, scales = 'free')+
        ggtitle(title, paste('Mean, ',ci[1]*100,'and',ci[2]*100,'%CI'))+
        xlab('') + ylab('')+
        scale_fill_manual(values = col.scen)+
        scale_colour_manual(values = col.scen)+
        theme(text = element_text(size = 14),
              axis.ticks.x = element_blank(),
              axis.text.x = element_text(face ='bold'), 
              panel.grid.minor.y = element_blank(),
              panel.grid.major.x = element_blank(),
              strip.text = element_text(face='bold'),
              strip.background = element_rect(fill = 'grey95'))+
        guides(colour = FALSE, fill = FALSE)
    
    return(list(val = g.val, proba = g.proba))        
}


compare_pkt_intrv <- function(prov, 
                              digest.peaks, 
                              targetdate,
                              col.scen, 
                              ci = c(0.50, 0.75, 0.99)) {
    
    if(0){ # DEBUG
        targetdate  <- "2020-05-15"
    }
    
    subdf <- sub_df(prov, digest.peaks, varname = 'pkt')
    
    # ---- CHANGE TABLE
    
    df.chg <- calc_chg_baseline(subdf, relative = FALSE)
    print(paste("Mean Peak Timing", prov))
    print(df.chg)
    
    # ---- PROBABILITIES
    
    title = paste("Proba peak date beyond target date", prov)
    print(title)
    print('target date: ')
    print(targetdate)
    
    proba.death <- subdf %>%
        filter( type == 'hosp') %>%
        group_by(scenario) %>%
        summarise(p = mean(pk.date > targetdate))
    print(proba.death)
    
    # --- Summary Stats
    
    anchor.date <- subdf$pk.date[1] - subdf$pkt[1]
    q <- ci_to_quantiles(ci)
    
    ss <- subdf %>% 
        filter(type == 'hosp') %>% # WARNING: DOING ONLY HOSPIALIZATION FOR NOW...
        group_by(scenario) %>%
        summarize(m = mean(v) + anchor.date,
                  q.lo1 = quantile(v, probs = q[1])+ anchor.date,
                  q.lo2 = quantile(v, probs = q[2])+ anchor.date,
                  q.lo3 = quantile(v, probs = q[3])+ anchor.date,
                  q.hi3 = quantile(v, probs = q[4])+ anchor.date,
                  q.hi2 = quantile(v, probs = q[5])+ anchor.date,
                  q.hi1 = quantile(v, probs = q[6])+ anchor.date
        ) %>%
        mutate(sdisp = scenario_display(scenario))
    
    # --- Changes between mean peaks:
    
    idx <- which(ss$scenario=='BASELINE')
    m <- ss$m
    names(m) <- ss$scenario
    m
    dp <- vector()
    rng <- 1:length(m)
    rng <- rng[rng!=idx]
    rng
    for(i in seq_along(rng)){
        dp[i] = m[rng[i]] - m[idx]
        names(dp)[i] = names(m)[rng[i]]
    }
    dp
    # TODO
    # Export for reports? 
    # MAke arrows between means with label of day diff
    
    # ---- PLOT PROBA
    
    title.pr = paste('Probability Hospitalization Peak Date is Beyond:',
                     format(ymd(targetdate), "%b %d"))
    g.proba <- proba.death %>%
        mutate(sdisp = scenario_display(scenario)) %>%
        ggplot(aes(x = sdisp, y = p, fill = sdisp))+
        geom_bar(stat = 'identity')+
        scale_y_continuous(limits = c(0,1), labels = scales::percent)+
        scale_fill_manual(values = col.scen)+
        theme(text = element_text(size = 16),
              axis.ticks.x = element_blank(),
              axis.text.x = element_text(face = 'bold'),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_blank())+
        xlab('') + ylab('')+
        guides(fill=F)+
        ggtitle(title.pr)
    
    # --- PLOT VALUES
    
    alpha = 0.3
    thick.seg = 16
    
    g <- ss %>%
        ggplot()+
        geom_segment(aes(x = q.lo1, xend = q.hi1, 
                         y = sdisp, yend = sdisp,
                         colour = sdisp),
                     size = thick.seg, alpha=alpha) +
        geom_segment(aes(x = q.lo2, xend = q.hi2, 
                         y = sdisp, yend = sdisp,
                         colour = sdisp),
                     size = thick.seg, alpha=alpha) +
        geom_segment(aes(x = q.lo3, xend = q.hi3, 
                         y = sdisp, yend = sdisp,
                         colour = sdisp),
                     size = thick.seg, alpha=alpha) +
        geom_point(aes(x = m, y = sdisp, colour=sdisp), 
                   size = thick.seg/2)+
        scale_x_date(date_breaks = '1 month', date_labels = '%d-%b')+
        scale_fill_manual(values = col.scen) + 
        scale_colour_manual(values = col.scen) + 
        theme(text = element_text(size = 14),
              axis.text.y = element_text(face='bold'),
              panel.grid.major.y = element_blank(),
              axis.ticks.y = element_blank())+
        ggtitle(paste("Hospitalizations Peak Date --", prov), 
                paste('Mean and',paste(ci*100,collapse = ','),'%CI'))+
        guides(colour = F, fill = F)+
        xlab('Peak Date') + ylab('')
    
    return(list(val = g, 
                proba = g.proba))
}


plot_compare_traj_intrv <- function(prov, digest.ss) {
    
    alpha= 0.15
    
    n = 1:6
    b1 <- 10^n
    b2 <- b1/2
    br <- sort(c(b1,b2))
        
    g <- sub_df(prov, df = digest.ss, varname='m') %>% 
        filter(md > 0.5, qvlo>0.5, qlo>0.5) %>%
        ggplot(aes(x=date, 
                   fill = scenario))+
        scale_y_log10(labels=scales::comma_format(accuracy=1), breaks = br)+
        scale_x_date(date_breaks = '1 month', date_labels = '%b-%d')+
        geom_line(aes(y = md, colour = scenario))+
        geom_ribbon(aes(ymin = qvlo, ymax=qvhi), 
                    alpha=alpha) + 
        geom_ribbon(aes(ymin = qlo, ymax=qhi), 
                    alpha=alpha) + 
        # facet_grid(scenario~type, scales='free_y')+
        facet_wrap(~type, nrow = 1, scales='free_y')+
        theme(axis.text.x = element_text(angle=45,hjust=1),
              panel.grid.minor.y = element_blank())+
        xlab('') + ylab('Count')+
        ggtitle(paste('Trajectory Comparison -', prov))+
        guides(colour=FALSE, fill = FALSE)
    
    plot(g)
}

#' Probability the cumul number of case for the next X days
#' will be larger than the past X days. 
proba_larger_next <- function(prov, 
                              sims, 
                              dat.all, 
                              n.days) {
    
    if(prov != 'Canada'){
        dat = dat.all %>%
            filter(Province == prov,
                   name == 'positive')    
    }
    
    if(prov == 'Canada'){
        dat = dat.all %>%
            filter(name == 'positive') %>%
            group_by(Date, name) %>%
            summarise(v = sum(value)) %>%
            rename(value = v)
    }
    
    
    last.date.obs = max(dat$Date)
    
    sum.past = sum(dat$value[dat$Date > last.date.obs - n.days])
    
    df = sims %>%
        filter(Province == prov,
               type == 'testpos',
               date > last.date.obs, 
               date <= last.date.obs + n.days) 
    z = df %>%
        group_by(scenario, iter) %>%
        summarize(greater = sum(value) > sum.past,
                  dbl  = sum(value) > 2 * sum.past,
                  half = sum(value) < 0.5 * sum.past) %>%
        group_by(scenario) %>%
        summarize(p.greater = mean(greater),
                  p.dbl  = mean(dbl),
                  p.half = mean(half))
    return(z)
}

plot_proba_larger <- function(prov, 
                              sims, 
                              dat.all, 
                              n.days,
                              col.scen) {
    p = proba_larger_next(prov, sims, dat.all, n.days) 
    g <- p %>%
        filter(scenario != 'ISO2') %>% # TODO: remove this when today is > April 30th.
        mutate(s = scenario_display(scenario)) %>%
        ggplot()+ 
        geom_bar(aes(x=s, y = p.greater, fill = s),
                 stat = 'identity')+
        scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
        scale_fill_manual(values = col.scen)+
        theme(text = element_text(size = 22),
              title = element_text(size = 18),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.y  = element_blank(),
              panel.border = element_blank(),
              axis.ticks.x  = element_blank(),
              axis.ticks.y  = element_blank(),
              axis.text.x = element_text(face='bold'))+
        guides(fill=F)+
        xlab('') + ylab('')+ 
        ggtitle(paste('Proba. cumul. confirmed cases next week will be higher than last week -', prov))
    plot(g)
}


proba_cumdeath_above <- function(prov, 
                                 sims, 
                                 thresh.value, 
                                 thresh.date) {
    
    # thresh.date = as_date('2020-04-30')
    # thresh.value = 1600
    
    d = sims %>%
        filter(type == 'death', Province == prov) %>%
        filter(date <= thresh.date) %>%
        group_by(iter, scenario) %>%
        mutate(cumdeath = cumsum(value)) %>%
        summarize(v = max(cumdeath))
    
    #TODO: Remove or do something smarter... 
    if(thresh.date <= as_date('2020-05-15')){
        d <- filter(d, scenario != 'ISO2') %>%
            filter(scenario != 'BASELINE')
    }
    
    g = ggplot(d)+
        geom_density(aes(x=v, fill=scenario),
                     color=NA,
                     alpha=0.8)+
        geom_vline(xintercept = thresh.value, linetype='dashed')+
        scale_fill_manual(values = col.scen[2]) + 
        scale_x_continuous(breaks = c(thresh.value)* c(0.5,1,1.5))+
        theme(panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              strip.background = element_rect(fill='grey95'),
              strip.text = element_text(face='bold'),
              panel.spacing.y = unit(1,"lines"),
              axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)),
              text = element_text(size=14))+
        facet_wrap(~scenario, ncol=1, scales='free_y')+
        guides(fill=FALSE)+
        ggtitle(paste('Cumulative Deaths Distribution by',
                   format(thresh.date, '%B %d'),'--',prov)) +
        ylab('') + xlab('Cumulative Number of Deaths')
    
    plot(g)
    
    # Probabilistic outputs:
    
    p.above = d %>% group_by(scenario) %>%
        summarise(p = mean(v > thresh.value))
    
    avg = d %>% group_by(scenario) %>%
        summarise(m = mean(v))
    
    for(i in 1:nrow(p.above)){
        write(round(100*p.above$p[i]), 
              file = paste0(dir.out,'cumdeath-shortterm-proba-above-',
                            prov,'-',p.above$scenario[i], '.txt'))
    }
    
    for(i in 1:nrow(avg)){
        write(round(avg$m[i]), 
              file = paste0(dir.out,'cumdeath-shortterm-mean-',
                            prov,'-',avg$scenario[i], '.txt'))
    }
    
}


plot_abc_post_prm <- function(abc.post, 
                              prm.name, 
                              scenario, 
                              ci = c(0.5, 0.95),
                              title = '',
                              one.minus = FALSE) { 
    # prm.name = 'lambda'
    # scenario = 'ISO1'
    # one.minus = T
    df = abc.post %>%
        filter(name == prm.name, scenario == !!scenario)
    
    df$y = df$value
    if(one.minus) df$y = 1 - df$value
    
    q = ci_to_quantiles(ci)
    
    dfs = df %>% 
        group_by(Province) %>%
        summarize(m = mean(y),
                  qvlo = quantile(y,probs = q[1], na.rm = TRUE),
                  qlo  = quantile(y,probs = q[2], na.rm = TRUE),
                  qhi  = quantile(y,probs = q[3], na.rm = TRUE),
                  qvhi = quantile(y,probs = q[4], na.rm = TRUE)
                  )
    
    sz = 9
    alpha = 0.20
    col = 'green4'
    g = dfs %>% 
        ggplot(aes(x=Province)) + 
        # Mean and CI
        geom_point(aes(y = m), size = sz*0.8, colour=col)+
        geom_segment(aes(xend=Province, y = qvlo, yend = qvhi), 
                     size = sz,
                     alpha=alpha, colour=col) + 
        geom_segment(aes(xend=Province, y = qlo, yend = qhi), 
                     size = sz,
                     alpha=1.5*alpha, colour=col) + 
        # Text Values
        geom_text(aes(x=Province, y=m, 
                      label=paste0(round(100*m),'%')),
                  colour = col,
                  fontface = 'bold',
                  hjust = -0.7)+
        # Cosmetics
        scale_y_continuous(labels = percent_format(accuracy = 1))+
        coord_cartesian(ylim=c(0,1))+
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              text = element_text(size=14), 
              axis.text.x = element_text(face = 'bold', size = 18),
              axis.ticks = element_blank(), 
              panel.border = element_blank())+
        ggtitle(title, "since March 20th")+
        xlab('') + ylab('Transmission rate reduction')
    plot(g)
}




