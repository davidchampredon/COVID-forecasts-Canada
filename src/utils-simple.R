suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
theme_set(theme_bw())
suppressPackageStartupMessages(library(lubridate))

library(ggrepel)

source('resude.R')


#' Converts Confidence intervals to quantiles
#' @param ci Numerical vector of CIs.
ci_to_quantiles <- function(ci) {
    qt <- c ( rev(0.5 - ci/2), 0.5 + ci/2)
    return(qt)
}


#' Converts mean and CV into shape and rate parameter
#' for the gamma distribution:
gamma_prm_conv <- function(m, cv) {
    # variance:
    v     <- (cv * m)^2
    return( c(shape = m^2/v, rate = m/v) )
}


#' Get a parameter value from its name in a CSV file
#' @param prm.name String. Name of the parameter
#' @param csv.name String. Name of the CSV file. 
get_prm <- function(prm.name, csv.name) { 
    # csv.name='prmdata.csv' ; prm.name='pop_size'
    x <- read.csv(csv.name, header = F, strip.white = T)
    val <- x[x[,1]==prm.name, 2]
    return(val)
}


#' Retrieve the first time number of tested individuals is >1.
get_idx_first_tested <- function(epi) {
    i.start <- which(epi$tested > 1)[1]
    return(i.start)
}

#' Retrieve the first time the _reported_ number of tested individuals is >0.
get_idx_first_tested_rep <- function(epi) {
    i.start <- which(epi$tested.rep > 0)[1]
    return(i.start)
}

get_idx_last_latent_incidence <- function(epi){
    return( max(which(epi$inc>0)) + 1 )
}


extract_report_times <- function(epi) {
    return(epi$time[!is.na(epi$tested.rep)])
}


#' y=x for small values, then saturates as y approaches L.
saturation_fct <- function(x, L) {
    a = -L
    b = 2*L
    c = 2/L
    return(a + b/(1+exp(-c*x)))
}
if(0){ #debug
    L <- 5000
    x <- seq(0,5*L,length.out = 1e3)
    y=saturation_fct(x, L)
    plot(x, y, typ = 'l', lwd=6)
    grid()
    abline(h=L, lty=2)
    abline(a=0, b=1, lty=2)
}


#' Crop the time range of the 
#' synthetic epidemic to a realistic range.
#' The time starts at the first reported test
#' and the length is specified by the number 
#' of reporting times.
#' 
crop_synthetic_epi <- function(n.reports) { # n.reports = 3
    
    rep.t.all <- extract_report_times(epi)
    
    t.start <- get_idx_first_tested_rep(epi)
    idx1 <- which(rep.t.all == t.start)
    t.end <- rep.t.all[idx1 + n.reports -1 ]
    
    epi.crop <- filter(epi, t.start <= time & time <=t.end)
    
    epi.crop.rep <- filter(epi.crop, !is.na(tested.rep))
    
    return(list(
        epi.crop = epi.crop,
        epi.crop.rep = epi.crop.rep
    ))
}

list_to_matrix <- function(x) {
    nr <- length(x)
    M <- matrix( unlist(x), byrow = T, nrow = nr)
    return(M)
}



#' Testing process
#' @param inc Integer vector. Latent true daily incidence.
#' @param prms Named List of testing process parameters.
#' @return Time series of tested cases (from latent incidence).
#'  
testing_process <- function(inc, prms) {
    # tau: Testing effort as a proportion of the daily latent incidence.
    # lag.test:: Testing lag. 
    tau        <- prms[['tau']]
    lag.test   <- prms[['lag.test']]
    test.d.max <- prms[['test_daily_max']]
    
    inc.lagged <- dplyr::lag(x = inc, 
                             n = lag.test, 
                             default = 0) 
    # Note: non-existant lagged values set to 0 
    # (before epidemic start incidence=0)
    n <- length(inc.lagged)
    m <- saturation_fct(tau * inc.lagged, test.d.max)
    # TODO: use `distrib.prms`
    tested <- rpois(n, lambda = m) 
    return(tested)
}


#' Testing process
#' @param tested Integer vector. Time series of tested cases
#' @param prms Named List of testing process parameters.
#' @return Time series of positive cases.
#'  
positivity_process <- function(tested, prms) { 
    # DEBUG prms = prms.positive
    # rho: Positivity rate of tested patients.
    rho <- prms[['rho']]
    n   <- length(tested)
    tested <- rpois(n, lambda = rho * tested) 
    # TODO: use `distrib.prms`
    return(tested)
}


#' Generate a synthetic epidemic.
#' 
incidence_process<- function(file.prm, distrib.proc) {
    # file.prm = 'prmdata.csv' 
    
    prmdata <- read.csv(file.prm, header = FALSE)
    message('Generating epidemic with the following parameters:')
    print(prmdata)
    pname <- trimws(as.character(prmdata[,1]))
    for(i in 1:length(pname)){ 
        assign(pname[i], prmdata[i,2], pos = 1)
    }
    
    # Generate full epidemic
    dat <- resude_simulate(pop_size = pop_size, 
                           I.init = I0,
                           R0 = R0,
                           alpha = alpha,
                           kappa = kappa, 
                           GI_span = GI_span, 
                           GI_mean = GI_mean, 
                           GI_var = GI_var,
                           horizon = horizon,
                           distrib.proc = distrib.proc,
                           seed = 1234)
    return(dat)
}


reporting_process <- function(epi, report_times) {
    
    rt <- report_times
    nr <- length(rt)
    
    # Reported tests:
    U <- numeric(length = nr)
    # Reported positive tests:
    Y <- numeric(length = nr)
    
    U[1] = sum(epi$tested[1:rt[1]])
    Y[1] = sum(epi$positive[1:rt[1]])
    
    for(k in 2:nr){
        U[k] <- sum(epi$tested[(rt[k-1]+1):rt[k]])
        Y[k] <- sum(epi$positive[(rt[k-1]+1):rt[k]])
    }    
    
    # Checks everythiing went well:    
    stopifnot(sum(U) == sum(epi$tested))
    stopifnot(sum(Y) == sum(epi$positive))
    
    epi$tested.rep       <- NA
    epi$positive.rep     <- NA
    epi$tested.rep[rt]   <- U
    epi$positive.rep[rt] <- Y
    return(epi)
}




#' Generate the synthetic epidemic.
gen_epi <- function(file.prm, 
                    prms.tested, 
                    prms.positive) {
    
    # file.prm = 'prmdata.csv' 
    
    # Latent incidence and susceptible:
    dat <- incidence_process(file.prm,
                             distrib.proc = 'rpois')
    nt <- length(dat$I)
    
    tested   <- testing_process(dat$I, prms.tested)
    positive <- positivity_process(tested, prms.positive)
    
    tmp <- data.frame(time     = 1:nt,
                      inc      = dat$I, 
                      S        = dat$S,
                      tested   = tested, 
                      positive = positive)
    
    # Reporting/aggregation process:
    report_times <- sort(sample(1:nt, size = round(nt/3), replace = F))
    epi <- reporting_process(tmp, report_times)
    
    return(epi)
}


plot_synthetic_epi <- function(epi, obs.limits = NULL) {  # obs.limits = c(15, 33)
    
    dp <- epi %>%
        pivot_longer(cols = -time , names_to = 'type')
    
    # Only reported data
    d.rep <- dp %>%
        filter(grepl(pattern = "\\.rep", type)) %>%
        filter(!is.na(value)) 
    
    tmax <- get_idx_last_latent_incidence(epi)
    
    g <- dp %>%
        filter(type != 'S') %>%
        filter(!grepl(pattern = "\\.rep", type)) %>%
        ggplot(aes(x=time, y = value, color=type)) + 
        geom_step(size=1)+
        geom_step(data = d.rep, aes(x=time, y=value, colour=type))+
        geom_point(data = d.rep, aes(x=time, y=value, colour=type))+
        scale_y_log10()+
        scale_x_continuous(limits = c(0,tmax))+
        scale_color_manual(values = c('black',
                                      'tomato1', 'tomato3',
                                      'steelblue1','steelblue3')) + 
        ggtitle('Synthetic Epidemic')
    
    if(!is.null(obs.limits)){
        g <- g + geom_vline(xintercept = obs.limits, linetype = 'dashed')
    }
    plot(g)
    plot(g + facet_wrap(~type) + theme(strip.text = element_text(size=12)))
}



#' Get all the csv files of tested individuals.
#' (get_data.R must be run beforehand)
get_all_data_csv <- function() {
    dats <- system('ls ../data/datax-*csv', intern = TRUE)
    
    tmp <- list()
    for(i in seq_along(dats)){
        tmp[[i]] <- read.csv(dats[i])
    }
    df <- do.call('rbind', tmp)
    
    df <- df %>%
        mutate(Date = lubridate::as_date(Date)) %>%
        filter(!is.na(value))
    return(df)  
}

#' Plot the dataframe returned by get_all_data().
#' 
plot_all_data <- function(df) {
    g <- df %>%
        ggplot(aes(x=Date, y=value, colour = name)) +
        geom_step(alpha=0.3, size=1) + 
        geom_point() + 
        geom_text(aes(label=value), 
                  alpha = 1,
                  size = 2.5, 
                  angle = 0, 
                  # hjust=-0.6, 
                  vjust = -1,
                  check_overlap = T,
                  fontface='bold')+
        facet_grid(name~Province)+
        scale_y_log10(limits = c(1, 1.3*max(df$value,na.rm = T))) + 
        scale_x_date(date_breaks = '7 days',
                     date_minor_breaks = '1 day', 
                     date_labels = '%b-%d')+
        theme(axis.text.x = element_text(angle=45, hjust  = 1),
              strip.text = element_text(size=12), 
              panel.grid.minor.y = element_blank())+
        ylab('Counts') + xlab('')+
        guides(colour=FALSE)+
        ggtitle('Reported Number of SARS-CoV-2 Tests Performed by Province')
    
    return(g)  
}



#' Plot the dataframe returned by get_all_data().
#' 
plot_data_report_canada <- function(df) { # df = dat.all
    
    df2 = df %>%
        mutate(Type = name)
    
    df.valtxt = df2 %>% 
        # filter(name == "positive") %>%
        group_by(Province) %>%
        filter(Date == max(Date)) 
    
    maxdate = max(df.valtxt$Date)
    mindate = max(df.valtxt$Date)
    
    col = c('tomato', 'orange2', 'steelblue')
    
    g <- df2 %>%
        ggplot(aes(x=Date, y=value, 
                   colour = Type, shape=Type)) +
        geom_step(alpha = 0.4, size=1) + 
        geom_point() + 
        geom_text(
            data = df.valtxt,
            aes(label=value, x=Date+3), 
            alpha = 1,
            size  = 3, 
            angle = 0, 
            hjust = 0.5,
            vjust = -1,
            check_overlap = T,
            fontface = 'bold',
            show.legend = FALSE) +
        # Facetting:
        facet_wrap(~Province, scales = 'free_x', nrow=2) +
        #
        # Cosmetics:
        scale_y_log10(limits = c(1, 1.3*max(df$value,na.rm = T)),
                      label = scales::comma_format(accuracy = 1)) + 
        scale_x_date(date_breaks = '7 days',
                     # date_minor_breaks = '1 day', 
                     date_labels = '%b-%d', ) +
        scale_color_manual(values = col) +
        theme(axis.text.x = element_text(angle=45, hjust  = 1),
              strip.text = element_text(size=12), 
              panel.grid = element_blank()) +
        ylab('Counts') + xlab('') +
        theme(text=element_text(size=16), 
              strip.background = element_rect(fill='grey95'),
              strip.text = element_text(face='bold'),
              legend.title = element_blank(),
              axis.text.x = element_text(size = 10)) +
        ggtitle('Reported Number of SARS-CoV-2 Tests')
    plot(g)  
}




#' Get generation interval parameters:
get_GI_prms <- function() {
    
    gi.m = get_prm('GI_mean', 'prm-model.csv')
    gi.v = get_prm('GI_var', 'prm-model.csv')  #gi.shape / gi.rate^2
    gi.phi  = gi.v / gi.m - 1
    
    return(list(
        mean       = gi.m,
        variance   = gi.v,
        gam.disp   = gi.phi,
        gam.shape  = gi.m / (1 + gi.phi),
        gam.rate   = 1 / (1 + gi.phi)
    ))
}


plot_growth_rate_cone <- function(prov=NULL) {
    
    df0 <- get_all_data_csv()
    
    df = df0
    if(!is.null(prov)) df <- filter(df0, Province==prov)
    
    df <- df %>%
        filter(name=='positive',
               value > 1 )
    date.min <- min(df$Date)
    date.max <- max(df$Date)
    
    # Calculate the exponential growth lines:
    dd <- seq(date.min, date.max, 1)
    ndd = length(dd)
    dfr <- data.frame(Date = dd, 
                      r1 = exp(log(2)/1 * 1:ndd),
                      r2 = exp(log(2)/2 * 1:ndd),
                      r3 = exp(log(2)/3 * 1:ndd),
                      r4 = exp(log(2)/4 * 1:ndd),
                      r5 = exp(log(2)/5 * 1:ndd))
    
    # Cosmetics:
    brk.y <- c((1:9)*10, (1:9)*100)
    title = paste('Daily Positives Tests Growth for',
                  ifelse(is.null(prov), 'Canada',prov))
    ndd2 <- round(ndd/2)
    ndd3 <- round(ndd*2/3)
    pos.txt <- 3
    
    g <-  df %>%
        ggplot()+
        geom_line(aes(x=Date, y=value, colour=Province), 
                  size = 2, alpha = 0.9)+
        geom_point(aes(x=Date, y=value, colour=Province),
                   pch=21, size=1, fill='white')+
        scale_x_date(date_breaks = '1 week', 
                     date_labels = '%b-%d')+
        scale_y_log10(minor_breaks = brk.y)+
        
        # reference doubling lines:
        geom_line(data=dfr, aes(x=Date, y=r3), linetype = 'dashed')+
        geom_line(data=dfr, aes(x=Date, y=r5), linetype = 'dashed')+
        geom_text(data=data.frame(x=dfr$Date[ndd-5], y=dfr$r3[ndd]), 
                  aes(x=x, y=y),
                  label="double every 3 days", angle=0)+
        
        geom_text(data=data.frame(x=dfr$Date[ndd-3], y=dfr$r5[ndd3]), 
                  aes(x=x, y=y),
                  label="double every 5 days", angle=0) +
        # theme stuff:
        xlab('')+ylab('Count')+
        theme(panel.grid.minor.x = element_blank(),
              text = element_text(size=16))+
        ggtitle(title, 'log scale')+
        guides(colour=FALSE)
    
    if(is.null(prov)){
        g <- g + facet_wrap(~Province)
    }
    plot(g)
}


plot_dbl_time <- function() {
    gi.prms <- get_GI_prms()
    dfr <- run_estim_r(g.mean = gi.prms$mean , 
                       g.disp = gi.prms$gam.disp) 
    n = nrow(dfr)
    dx = 0.2
    
    g <- dfr %>%
        ggplot(aes(x = prov))+
        geom_pointrange(aes(y=dbl, ymin=dbl.lo, ymax=dbl.hi),
                        size = 1, pch = 21, fill='white', 
                        stroke = 2)+
        scale_y_continuous(breaks = c(1:9, seq(10,20,by=2)))+
        annotate("text", x = n + dx, y = max(dfr$dbl), 
                 label = "italic(Slower)", size = 6, 
                 colour = 'green4', parse=TRUE)+
        annotate("text", x = n + dx, y = 0.5, 
                 label = "italic(Faster)", size = 6, 
                 colour = 'tomato', parse=TRUE)+
        xlab('') + ylab('Days')+
        ggtitle('Doubling Time of Daily Confirmed Cases', 'Mean and 95%CI')+
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              text = element_text(size=18),
              axis.text.x = element_text(size=22,face = 'bold'),
              axis.ticks.x = element_blank())
    plot(g)
}


plot_R0 <- function() {
    gi.prms <- get_GI_prms()
    dfr <- run_estim_r(g.mean = gi.prms$mean , 
                       g.disp = gi.prms$gam.disp) 
    n = nrow(dfr)
    col = 'tomato3'
    
    g <- dfr %>%
        ggplot(aes(x = prov))+
        geom_pointrange(aes(y=R0, ymin=R0.lo, ymax=R0.hi),
                        size = 1.5, pch = 22, fill='white', 
                        stroke = 2,
                        colour = col)+
        geom_hline(yintercept = 1, linetype = 'dashed')+
        scale_y_continuous(breaks = seq(1,4,by=0.5))+
        xlab('') + ylab('')+
        ggtitle('Basic Reproduction Number (from Confirmed Cases)', 'Mean and 95%CI')+
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              text = element_text(size=18),
              axis.text.x = element_text(size=22,face = 'bold'),
              axis.ticks.x = element_blank())+
        coord_cartesian(ylim = c(0.8,4))
    plot(g)
}


write_output_num <- function(val, rounding, name, out.dir = 'out/') {
    write.table(round(val,rounding),
                file = paste0(out.dir, name,'_', prov,'.txt'), 
                quote = F, row.names = F, col.names = F, eol = '')
}

write_output_bignum <- function(val, name, out.dir = 'out/') {
    write.table(format(round(val,0), big.mark=","),
                file = paste0(out.dir, name,'_', prov,'.txt'), 
                quote = F, row.names = F, col.names = F, eol = '')
}



plot_data_prov <- function(prov, dat.all) {
    
    # prov = 'ON'
    
    df <- dat.all %>%
        filter(Province == prov)
    sz = 2
    sz.txt = 3
    
    g <- df %>%
        ggplot(aes(x=Date, y=value, colour = name)) +
        geom_step(alpha=0.3, size=sz) + 
        geom_point(size = sz) + 
        geom_text(data = filter(df,name=='positive'),
                  aes(label=value), 
                  alpha = 1,
                  size = sz.txt, 
                  angle = 0, 
                  # hjust=-0.6, 
                  vjust = -1,
                  check_overlap = T,
                  fontface='bold')+
        scale_y_log10(limits = c(1, 2 * max(df$value,na.rm = T)),
                      labels = scales::comma_format(accuracy = 1),
                      breaks = c(1,10,50,100,500,1000,5000,10000)) + 
        scale_x_date(date_breaks = '7 days',
                     # date_minor_breaks = '1 day', 
                     date_labels = '%b-%d')+
        theme(text = element_text(size = 20),
              panel.grid.minor.x = element_blank(),
              axis.text.x = element_text(angle=0),
              panel.grid.minor.y = element_blank())+
        ylab('Counts') + xlab('')+
        guides(colour = guide_legend(title="Type"))+
        ggtitle(paste('Reported Number of SARS-CoV-2 Tests in',prov))
    
    plot(g)
}


#' Shifts the values of a vector
#' source: https://clarkrichards.org/r/timeseries/2016/02/09/a-function-to-shift-vectors/
shift_lag <- function(x, lag) {
    n <- length(x)
    xnew <- rep(NA, n)
    if (lag < 0) {
        xnew[1:(n-abs(lag))] <- x[(abs(lag)+1):n]
    } else if (lag > 0) {
        xnew[(lag+1):n] <- x[1:(n-lag)]
    } else {
        xnew <- x
    }
    return(xnew)
}

calc_dbl_time <- function(x, level) {
    x <- x[x>0]
    
    if(length(x)==0) return(NA)
    
    d = data.frame(logv = log(x), time = 1:length(x))
    q  = lm(data = d, logv ~ time)
    k  = coef(q)
    r  = max(0,k['time']  )
    # Doubling times:
    estim = log(2) / r
    
    if(is.null(level)) 
        res = c(estim = estim)
    
    if(!is.null(level)){
        ci = confint(q, level = level)
        ci = pmax(0,ci[2,])
        res = c(
            estim = estim, 
            hi    = log(2) / ci[1], 
            lo    = log(2) / ci[2]
        )
    }
    return(res)
}



plot_data_prov_dbl <- function(prov, w) {
    
    df0 = get_all_data_csv()
    df  = df0
    if(!is.null(prov)) df <- filter(df0, Province==prov)
    
    # df <- df %>%
    #     filter(value > 1 )
    # 
    date.min <- min(df$Date)
    date.max <- max(df$Date)
    
    d1 = date.max - w
    d2 = date.max
    
    # For doubling times:
    dd = df %>%
        filter(d1 <= Date & Date <= d2) %>%
        group_by(name) %>%
        summarize(d = calc_dbl_time(value, level = NULL),
                  y1 = value[Date==d1], 
                  y2 = value[Date==d2],
                  ymax = max(value)) %>%
        mutate(d1 = d1, d2 = d2 )
    
    # --- Plot
    
    y.maj = c(2.5, 5, 7.5) #seq(2,8,2)
    yy    = 10^c(0:4)
    brk.y = sort(c(2, yy, c(y.maj*10, y.maj*100, y.maj*1000) ))
    
    title = paste('Daily Count of Tests, Positive Tests and Deaths',
                  ifelse(is.null(prov), 'Canada',prov))
    
    # cosmetics:
    sz.ln = 1.5
    sz.pt = 2.0
    sz.seg = 5
    alpha.seg = 0.2
    alpha.ln  = 0.3
    
    col.ln = c('red3', 'orange2', 'steelblue3')
    
    g = df %>% 
        ggplot() + 
        # --- Trend
        geom_segment(data = dd, 
                     aes(x=d1, xend=d2, y=y1, yend=y2),
                     size = sz.seg, alpha = alpha.seg)+
        # --- Tests & deaths:
        geom_step(aes(x=Date, y=value, 
                      colour = name),
                  size = sz.ln, alpha = alpha.ln)+
        geom_point(aes(x=Date, y=value, 
                       colour=name),
                   shape = 15,
                   size = sz.pt)+
        # --- Label last count:
        geom_text(data = dd, aes(x = d2 + 5, 
                                 y = y2, 
                                 colour=name, 
                                 label=y2),
                  fontface = 'bold', show.legend = F)+
        # --- Label doubling time:
        geom_label(data = dd, aes(x=d2+3, y=ymax*1.4, 
                                  colour=name, 
                                  label=paste("Doubling time:",round(d),"days")),
                   fontface = 'bold', 
                   size = 5,
                   hjust=1, show.legend = F)+
        # --- Cosmetics:
        scale_y_log10(label = scales::comma_format(accuracy = 1),
                      breaks = brk.y)+
        scale_x_date(date_breaks = '1 week', 
                     date_labels = '%b %d', 
                     limits = c(date.min, date.max + 7))+
        scale_color_manual(values=col.ln)+
        theme(text = element_text(size=16),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.x = element_blank(),
              axis.text.x = element_text(angle = 45, hjust=1),
              axis.ticks.y = element_blank())+
        guides(colour = guide_legend(title="Type"), label=F)+
        xlab('') + ylab('Count')+
        ggtitle(title)
    
    return(g)
}



plot_dbl_evol <- function(prov, w, level) 
{
    df0 = get_all_data_csv()
    df  = df0
    if(!is.null(prov)) df <- filter(df0, Province==prov)
    
    df <- df %>%
        filter(value > 1 )
    
    date.min <- min(df$Date)
    date.max <- max(df$Date)
    
    axis.max.date = date.max + 7
    label.date  = date.max + 5
    
    d1 = date.min 
    d2 = date.min + w
    tmp = list() ; i=1
    while (d2 <= date.max) {
        tmp[[i]] = df %>%
            filter(d1 <= Date & Date <= d2) %>%
            group_by(name) %>%
            summarize(d = calc_dbl_time(value, level = NULL),
                      d.lo = calc_dbl_time(value, level)[['lo']],
                      d.hi = calc_dbl_time(value, level)[['hi']],
                      asof = d2)
        d1 = d1 + 1
        d2 = d2 + 1
        i  = i + 1
    }
    z = do.call('rbind',tmp)
    
    # z$d[z$d > 1e3]  <- NA
    # z$d.hi[z$d.hi > 1e3]  <- NA
    
    z.last = filter(z, asof == date.max)
    
    col.ln = c( 'steelblue2', 'orange2','red3')
    
    z <- within(z, name <- factor(name, 
                                  levels = c('tested', 'positive', 'death')))
    
    g = z %>%
        ggplot(aes(x = asof, y = d, 
                   colour = name)) +
        geom_ribbon(aes(ymin = d.lo, ymax=d.hi,
                        fill = name), 
                    colour = NA,
                    alpha=0.2)+
        geom_line(size = 2)+
        geom_point(aes(colour = name),
                   fill = 'white', 
                   shape = 21,
                   size = 1) +
        # --- Last value
        geom_label(data = z.last, 
                  aes(x = label.date, #asof + 3, 
                      y = pmax(5,pmin(d,50)), 
                      colour=name, 
                      fontface = 'bold',
                      label= paste(round(d),'days')))+
        # --- Fast/Slow
        annotate(geom = 'text', 
                 x = label.date, 
                 y = -2, 
                 colour = 'tomato',
                 label = 'italic(Faster)', parse=T)+
        annotate(geom = 'text', 
                 x = label.date, 
                 y = 57, 
                 colour = 'green3',
                 label = 'italic(Slower)', parse=T) +
        # --- Cosmetics:
        scale_colour_manual(values = col.ln) +
        scale_fill_manual(values = col.ln) +
        facet_wrap(~name, ncol = 1, scales='free_y') + 
        guides(colour = F, fill=F)+ 
        coord_cartesian(ylim=c(-5,60) , 
                        xlim = c(date.min+w, axis.max.date)) +
        theme(text = element_text(size = 16),
              panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(angle = 45, hjust=1),
              strip.background = element_rect(fill='grey95'),
              strip.text = element_text(face='bold'),
              plot.margin = unit(c(5.5,5.5,5.5,19), "pt"))+
        ggtitle('Evolution of Doubling Time', 
                paste0('Calculated on a ',w,'-day sliding window'))+
        xlab('Calculation Date') + 
        ylab('Doubling Time in Days')
    return(g)
}

