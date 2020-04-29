
library(EpiEstim)
source('utils-simple.R')


estim_r_one <- function(x, level) {
  
  min.n.data = 4
  
  if(nrow(x) >= min.n.data){
    x <- x[x$value>0,]
    x$logv <- log(x$value)
    q  = lm(data = x, logv ~ time)
    k  = coef(q)
    ci = confint(q, level = level)
    ci = pmax(0,ci[2,])
    r  = max(0,k['time']  )
    y0 = k['(Intercept)']  
  }
  if(nrow(x)< min.n.data){
    r = NA
    ci = c(NA,NA)
    y0 = NA
  }
  
  return(list(r = r, 
              ci = ci,
              intercept = y0))
}

estim_R0 <- function(r, g.bar, k) {
  return( (1 + k*r*g.bar)^(1/k) )
}


estim_r <- function(dat, level = 0.95) {
  df <- dat %>% 
    filter(name == 'positive')
  
  provs <- unique(df$Province)
  tmp <- list()
  
  for(i in seq_along(provs)){
    a <- df %>%
      filter(Province == provs[i]) %>%
      estim_r_one(level)
    tmp[[i]] <- data.frame(prov = provs[i],
                           r = a$r,
                           r.lo = a$ci[1],
                           r.hi = a$ci[2])
  }
  
  dfr <- do.call('rbind', tmp) %>%
    mutate(dbl.lo = log(2) / r.hi,
           dbl = log(2) / r,
           dbl.hi = log(2) / r.lo)
  
  return(dfr)
}


run_estim_r <- function(g.mean, g.disp) {  # g.mean=4.5 ; g.disp=0.56
  dat <- get_all_data_csv()
  dfr <- estim_r(dat) %>%
    mutate(R0   = estim_R0(r, g.bar = g.mean, k = g.disp),
           R0.lo = estim_R0(r.lo, g.bar = g.mean, k = g.disp),
           R0.hi = estim_R0(r.hi, g.bar = g.mean, k = g.disp))
  return(dfr)
}

estim_R0_simple <- function(prov) {
  gi.prms <- get_GI_prms()
  dfr <- run_estim_r(g.mean = gi.prms$mean , 
                     g.disp = gi.prms$gam.disp) %>% 
    rename(Province = prov) %>%
    filter(Province == prov)
  
  # TODO: change, this is rough!
  return(list(
    mean = dfr$R0,
    sd   = (dfr$R0.hi - dfr$R0.lo)/4
  ))
}


plot_r_regression <- function(prov) {
  
  gi.prms <- get_GI_prms()
  dat <- get_all_data_csv()
  
  df <- dat %>% 
    filter(name == 'positive', Province == prov)
  
  z <- estim_r_one(df,level=0.95)
  
  dfr <- run_estim_r(g.mean = gi.prms$mean , 
                     g.disp = gi.prms$gam.disp) %>% 
    rename(Province = prov) %>%
    filter(Province == prov)
  
  title = paste("Growth rate positive cases",prov,
                "dbl time:",round(dfr$dbl, 2),
                "\n R0 =",
                paste(round(dfr$R0.lo, 2),
                      round(dfr$R0,2),
                      round(dfr$R0.hi,2), 
                      sep = ' ; '))
  
  df.plot <- df[df$value>0,]
  plot(df.plot$time, log(df.plot$value), typ='b',
       main = title,
       ylab='log count', xlab = 'time', las=1)
  grid()
  abline(a=z$intercept, b=z$r, lwd=3)
  
  write_output_num(dfr$R0, 2, 'R0')
  write_output_num(dfr$R0.lo, 2, 'R0lo')
  write_output_num(dfr$R0.hi, 2, 'R0hi')
  
  write_output_num(dfr$dbl, 1, 'dbl')
  write_output_num(dfr$dbl.lo, 1, 'dbllo')
  write_output_num(dfr$dbl.hi, 1, 'dblhi')
  
  
}

#' Estimate r, doubling time  and R0  sliding 1 day backward.
run_estim_r_sliding <- function(g.mean, 
                                g.disp) {  # g.mean=4.5 ; g.disp=0.56
  dat <- get_all_data_csv() %>%
    filter(name == 'positive')
  
  dv <- sort(unique(dat$Date), decreasing = TRUE)
  
  N = 12
  
  tmp <- list()
  for(i in 1:N){
    tmp[[i]] <- dat %>%
      filter(Date <= dv[i]) %>%
      estim_r() %>%
      mutate(R0    = estim_R0(r, g.bar = g.mean, k = g.disp),
             R0.lo = estim_R0(r.lo, g.bar = g.mean, k = g.disp),
             R0.hi = estim_R0(r.hi, g.bar = g.mean, k = g.disp))
    tmp[[i]]$asof = dv[1] - i + 1
  }
  
  res = do.call('rbind', tmp)
  return(res)
}

plot_dbl_sliding <- function() {
  gi <- get_GI_prms()
  df <- run_estim_r_sliding(g.mean = gi$mean, 
                            g.disp = gi$gam.disp)
  
  df <- df[complete.cases(df),]
  
  col <- 'blue'
  
  g <- df %>%
    ggplot(aes(x=asof, y = dbl)) + 
    geom_point(colour = col) + 
    geom_line(size=1, colour = col)+
    geom_ribbon(aes(ymin=dbl.lo, ymax=dbl.hi), 
                alpha=0.2, fill = col)+
    facet_wrap(~prov, scales='free_y', nrow=2)+
    scale_y_continuous(breaks = c(1:7, 10, 14))+
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())+
    coord_cartesian(ylim=c(0,16))+
    ggtitle("Evolution of Doubling Time (confirmed cases)","mean, 95%CI")+
    xlab("Calculation Date") + ylab("Days")+
    theme(axis.text.x = element_text(angle=30, hjust=1),
          text = element_text(size=16),
          strip.text = element_text(size=18,face = 'bold'),
          strip.background = element_rect(fill='grey95'))
  plot(g)
}



plot_dbl_sliding_prov <- function(prov) {
  gi <- get_GI_prms()
  df <- run_estim_r_sliding(g.mean = gi$mean, 
                            g.disp = gi$gam.disp)
  
  df <- df[complete.cases(df),] %>%
    filter(prov == !!prov)
  
  col <- 'blue'
  title = paste('Evolution of Doubling Time for',prov)
  
  ymax <- 1.2 * max(df$dbl)
  
  g <- df %>%
    ggplot(aes(x=asof, y = dbl)) + 
    geom_line(size=1, colour = col)+
    geom_ribbon(aes(ymin=dbl.lo, ymax=dbl.hi), 
                alpha=0.2, fill = col)+
    geom_point(colour = col, pch = 21, fill='white', 
               size= 4, stroke=2) + 
    scale_y_continuous(breaks = c(1:7, 10, 14))+
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())+
    coord_cartesian(ylim=c(0, ymax))+
    ggtitle(title,
            "mean and 95%CI")+
    xlab("Calculation Date") + ylab("Days")+
    theme(axis.text.x = element_text(angle=0, hjust=0.5),
          axis.text.y = element_text(size = 18),
          text = element_text(size=14))
  
  g <- g  + 
    annotate("text", x = df$asof[2], y = ymax, 
             label = "italic(Slower)", size = 6, 
             colour = 'green4', parse=TRUE)+
    annotate("text", x = df$asof[2], y = 0.5, 
             label = "italic(Faster)", size = 6, 
             colour = 'tomato', parse=TRUE)
  
  plot(g)
}




plot_R0_sliding <- function() {
  gi <- get_GI_prms()
  df <- run_estim_r_sliding(g.mean = gi$mean, 
                            g.disp = gi$gam.disp)
  
  df <- df[complete.cases(df),]
  
  col <- 'red3'
  
  g <- df %>%
    ggplot(aes(x=asof, y = R0)) + 
    geom_point(colour = col) + 
    geom_line(size=1, colour = col)+
    geom_ribbon(aes(ymin=R0.lo, ymax=R0.hi), 
                alpha=0.2, fill = col)+
    facet_wrap(~prov, scales='free_y', nrow=2)+
    scale_y_continuous()+
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())+
    coord_cartesian(ylim=c(0.8,4))+
    ggtitle("Evolution of R0 (confirmed cases)","mean, 95%CI")+
    xlab("Calculation Date") + ylab('')+
    theme(axis.text.x = element_text(angle=30, hjust=1),
          text = element_text(size=16),
          strip.text = element_text(size=18,face = 'bold'),
          strip.background = element_rect(fill='grey95'))
  plot(g)
}

#' Estimate effective reproduction number
estim_Rt <- function() {
  
  gi <- get_GI_prms()
  
  dat <- get_all_data_csv() %>%
    filter(name=='positive')
  
  firstobs = dat %>%
    group_by(Province) %>%
    summarize(d = min(Date))
  firstobs
  
  # How far back in time to calculate:
  N = 18 # TODO: do not hard code... 
  last.d <- max(dat$Date)
  first.d = lubridate::as_date("2020-04-01")  #last.d - N
  
  df <- filter(dat, Date >= first.d)
  
  provs <- unique(dat$Province)
  tmp <- list()
  
  # Serial interval definition for EpiEstim
  # (assume same as generation interval):
  si <- make_config(list(mean_si = gi$mean, 
                         std_si  = sqrt(gi$variance)))
  
  usi <- make_config(list(
    # mean SI:
    mean_si     = gi$mean, 
    std_mean_si = 1,
    min_mean_si = gi$mean-2, 
    max_mean_si = gi$mean+2,
    # std. dev SI:
    std_si      = sqrt(gi$variance), 
    std_std_si  = 0.5,
    min_std_si  = sqrt(gi$variance)-0.5, 
    max_std_si  = sqrt(gi$variance)+0.5,
    # replicates:
    n1 = 100, n2 = 100))
  
  for(i in seq_along(provs))   # i=1
  {
    dfi = filter(dat, Province == provs[i])
    
    inc = dfi$value
    
    inc[inc<0] <- 0
    
    rt = EpiEstim::estimate_R(incid = inc,
                              method = "uncertain_si",
                              config = usi)
    
    z          = rt$R
    z$prov     = provs[i]
    z$calcdate = seq(max(dfi$Date)-nrow(z)+1, max(dfi$Date), by=1)
    
    tmp[[i]] <- z
  }
  
  res <- do.call('rbind', tmp)
  # res$calcdate = last.d - N + res$t_end -1  # TODO: double check
  return(res)
}



plot_Rt <- function() {
  df <- estim_Rt()
  
  # Cosmetics:
  a   = 0.3
  col = 'steelblue3'
  col.r1 = 'green3'
  d.shift = 3
  
  d0=min(df$calcdate) 
  d1=max(df$calcdate)
  dfpoly = data.frame(x = c(d0-2, d0-2,
                            d1+2*d.shift, d1+2*d.shift),
                      y = c(0,1,1,0))
  
  df.txt = df %>%
    group_by(prov) %>%
    filter(calcdate == max(calcdate)) %>%
    select(prov, calcdate, `Mean(R)`) %>%
    rename(y = `Mean(R)`)
  
  # Re-order West to East:
  w2e = c('BC', 'AB', 'SK', 'MB',
          'ON', 'QC', 'NS', 'NL')
  df <- within(df, 
               prov <- factor(prov, levels = w2e))
  
  
  g <- df %>% 
    ggplot(aes(x=calcdate, y=`Mean(R)`)) +
    geom_polygon(data=dfpoly, aes(x=x,y=y), 
                 fill=col.r1, alpha=0.3)+
    geom_line(colour = col, size=1.7)+
    geom_ribbon(aes(ymin=`Quantile.0.025(R)`, 
                    ymax=`Quantile.0.975(R)`),
                alpha = a, fill = col)+
    geom_point(colour = col, 
               fill = 'white',
               pch = 21, size = 1) +
    geom_label(data = df.txt, 
               aes(x = calcdate -2, 
                   y = 2.7, 
                   label = round(y,2)),
               hjust = 0.5,
               col = col, fontface = 'bold') +
    #
    # Cosmetics
    #
    scale_y_continuous(breaks = c(0:4,seq(6,10,by=2)))+
    scale_x_date(date_breaks = '1 week', date_labels = '%b %d')+
    coord_cartesian(ylim = c(0,4), xlim = c(d0,d1+d.shift))+
    facet_wrap(~prov, scales = 'fixed', nrow=2)+
    ggtitle("Evolution of the Effective Reproduction Number (from confirmed cases)",
            'Mean and 95%CI')+
    xlab('Calculation Date')  + ylab('')+
    theme(panel.grid.minor.y = element_blank(),
          strip.background = element_rect(fill='grey95'),
          strip.text = element_text(face='bold', size =12),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust=1, size=9),
          axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)))
  plot(g)
}

plot_Rt_prov <- function(prov) {
  
  df <- estim_Rt() %>%
    filter(prov == !!prov)
  
  # Cosmetics:
  a   = 0.2
  col = 'steelblue2'
  
  date.min = min(df$calcdate)
  date.max = max(df$calcdate)
  title = paste("Evolution of the Effective Reproduction Number for",prov)
  g <- df %>% 
    ggplot() +
    geom_rect(data = data.frame(xmin=date.min, xmax=date.max, ymin=0, ymax=1),
              aes(xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax), 
              fill='darkolivegreen2', alpha =0.6)+
    geom_line(aes(x=calcdate, y=`Mean(R)`), 
              colour = col, size = 2)+
    geom_ribbon(aes(x=calcdate, 
                    y=`Mean(R)`,
                    ymin=`Quantile.0.025(R)`, 
                    ymax=`Quantile.0.975(R)`),
                alpha = a, fill = col)+
    geom_point(aes(x=calcdate, y=`Mean(R)`), 
               colour = col, size =2, 
               stroke = 2, 
               pch = 21, fill='white')+ 
    scale_y_continuous(limits = c(0,max(df$`Quantile.0.975(R)`)),
                       breaks = 0:8)+
    
    ggtitle(title,'Mean and 95%CI')+
    xlab('Calculation Date')  + ylab('')+
    theme(text=element_text(size = 16),
          panel.grid.minor.y = element_blank())
  plot(g)
}