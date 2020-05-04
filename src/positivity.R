library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(ggrepel)
source('read-data.R')


pop = read.csv('../data/populations.csv', header = F)
names(pop) = c('Province', 'popsize')
pop$Province <- as.character(pop$Province)

fname <- 'https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/COVID19_Canada.csv'

x <- get_data(prov=NULL, fname) %>%
    select(Province, Date, total_testing, confirmed_positive) %>%
    group_by(Province) %>%
    mutate(tested = c(NA,diff(total_testing))) %>%
    mutate(positive = c(NA,diff(confirmed_positive)))

sort(unique(x$Province))
df  = x %>%
    filter(Province %in% c('AB',  
                           'MB',
                           'ON', 
                           'QC', 
                           'SK'
                           ))%>%   # Problems with NB, NL, etc.
    # pivot_wider(names_from = name, values_from = value) %>%
    left_join(pop, by='Province') %>%
    filter(positive>0 | tested >0) %>%
    mutate(posrate = positive / tested) %>%
    mutate(testprop = tested / popsize) %>%
    filter(posrate >0) %>%
    filter(!is.na(posrate)) %>% 
    filter(!is.infinite(posrate))

sort(unique(df$Province))

df = df %>%
    group_by(Province) %>%
    mutate(time = as.numeric(Date - min(Date)))


df2 = df %>% 
    group_by(Province) %>%
    filter(Date == max(Date))
    
pr = c(0.01, 0.02, 0.05, 0.1, 0.2)

dfpr = data.frame(x = rep(1, length(pr)),
                  y = pr,
                  pr = paste0(100*pr,'%'))

col.pr = 'steelblue3'

g = df %>%
    # group_by(Province) %>%
    # arrange(Date) %>%
    ggplot(aes(x=tested, y=positive)) + 
    geom_abline(slope = 1, 
                intercept = log10(pr), 
                colour=col.pr, size = 1, alpha=0.6) +
    geom_label(data = dfpr, aes(x=x,y=y,label = pr ), 
               colour = col.pr, fontface = 'bold')+
    geom_point(aes(fill=time), pch=21, colour='grey80') +
    geom_path(alpha=0.1)+
    geom_point(data=df2, 
               colour='red3', 
               fill='pink', 
               pch=21, size=3, stroke=2)+
    scale_fill_gradient(low='yellow', 
                         high = 'orange3',)+
    facet_wrap(~Province, scales = 'fixed')+
    scale_y_log10(label = scales::comma_format(accuracy = 1))+
    scale_x_log10()+
    guides(fill = FALSE)+
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill='grey95'),
          strip.text = element_text(face='bold'))+
    ggtitle('Daily Positivity Rate Evolution')
g



g1 = df %>%
    group_by(Province) %>%
    filter(Date >= max(Date) - 6) %>%
    # group_by(Province) %>%
    # arrange(Date) %>%
    ggplot(aes(x=tested, y=positive)) + 
    geom_abline(slope = 1, 
                intercept = log10(pr), 
                colour=col.pr, size = 1, alpha=0.8) +
    geom_label(data = dfpr, aes(x=x,y=y,label = pr ), 
               colour = col.pr, fontface = 'bold')+
    geom_point(aes(colour=Province, shape = Province), alpha=0.5, size=2) +
    geom_path(aes(colour=Province), alpha=0.1, size=1)+
    # geom_point(data=df2, aes(colour=Province, shape = Province), 
    #            size=3, stroke=2) +
    geom_label(data=df2, 
                     aes(x = tested, y=positive, 
                         label = Province,
                         colour=Province),
               fontface = 'bold',
               size = 6)+
    scale_y_log10(label = scales::comma_format(accuracy = 1))+
    scale_x_log10()+
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill='grey95'),
          strip.text = element_text(face='bold'))+
    guides(colour=FALSE, shape=F)+
    ggtitle('Daily Positivity Rate Today')
g1





g2 = df %>%
    mutate(Days.from.first.test = time) %>%
    ggplot(aes(x = testprop, y = posrate))+
    geom_path(alpha=0.1, col='orange') +
    geom_point(aes(colour=Days.from.first.test))+
    geom_point(data=df2, 
               colour='black', 
               fill='pink', 
               pch=21, size=3, stroke=2)+
    scale_color_gradient(low='yellow', high='red')+
    facet_wrap(~Province, scales = 'fixed')+
    scale_x_log10(limits = c(1e-8, 1)) +
    scale_y_log10()+
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill='grey95'),
          strip.text = element_text(face='bold'))+
    # guides(color=guide_legend(title="Days from first test"))+
    ggtitle('Positive vs. Tested')+
    xlab('Proportion tested (T/N)') + 
    ylab('Proportion positive (P/T)')
g2

pdf('tp.pdf')
plot(g2)
dev.off()
