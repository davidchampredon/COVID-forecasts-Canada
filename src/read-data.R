suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2)) 
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(gridExtra))
theme_set(theme_bw())

#' Get data 
#' @param prov String. Province name. If NULL, all provinces.
get_data <- function(prov, fname) {
    
    dat <- suppressMessages(readr::read_csv(fname)) %>%
        mutate(Date = lubridate::as_date(Date))
    
    if(!is.null(prov)) dat <- filter(dat, grepl(prov, Province))
    return(dat)
}

#' Reformat data for convenience.
#' @param test.pos Boolean. Only number tested and positives?
#' 
digest_data <- function(dat, test.pos) {  
    
    # DEBUG
    # dat = get_data('ON',fname)
    
    # If total_testing missing but 
    # positive / & negative exist, then calculate sum:
    idx <- is.na(dat$total_testing)
    if(length(idx) > 0){
        # dat$total_testing[idx] <- dat$negative[idx] + dat$confirmed_positive[idx]
        dat$total_testing[idx] <- dat$negative[idx] + dat$confirmed_positive[idx]
    }
    
    d <- dat %>%
        select(-source, -Note) %>%
        mutate(tested   = c(NA, diff(bestTotal)),
               # positive = c(NA, diff(confirmed_positive)),
               positive = newConfirmations, #c(NA, diff(newConfirmations)),
               death    = c(NA, diff(deceased)),
               time     = cumsum(c(1, diff(Date))))
    
    nm   <- names(d)
    idxc <- which (! nm %in% c('Province','Date', 'time'))
    
    res <- d %>% 
        pivot_longer(cols=idxc)  
    
    if(test.pos){
        res <- res %>%
            filter(name %in% c('tested', 'positive', 'death') )
    }
    return(res)
}


#' Plot data.
plot_data <- function(d, title = '') {
    g <- d %>%
        ggplot(aes(x=Date, y=value)) +
        geom_step(colour='darkgrey')+
        geom_point(size=2)+
        geom_text(aes(label=value), hjust=1.5, angle=90)+
        facet_wrap(~name, scales = 'free_y', ncol=1) + 
        scale_y_log10()+
        scale_x_date(date_breaks = '2 days', date_labels = '%d-%b')+
        theme(axis.text.x = element_text(angle=45, hjust = 1), 
              strip.text.x = element_text(size = 12))+
        ggtitle(title) + 
        xlab('') + ylab('')
    return(g)
}


#' Extract data of a province to a CSV file that
#' will be used as input for statistical model.
#' @param test.pos Boolean. Only number tested and positives?
extract_data_csv <- function(prov, test.pos, fname) {
    # DEBUG: prov='ON'
    
    csvname <- paste0('../data/datax-',prov,'.csv')
    
    dat <- get_data(prov, fname)
    
    d2 <- digest_data(dat = dat, 
                      test.pos = T) 
    
    d2%>%
        write.csv(file = csvname, 
                  quote = F, row.names = F)
    return(csvname)
}


#' Refomat data for using it as input for fit functions.
prep_data_for_fit <- function(dat) {
    d <- dat %>%
        filter(!is.na(value)) %>%
        pivot_wider(names_from = name, values_from = value) %>%
        filter(tested > 2)
    return(d)    
}

# ---- DEBUG ----

do.debug = FALSE  # switch to TRUE for debugging.

if(do.debug){
    
    #read_csv("https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/COVID-19_test.csv”)
    
    fname <- 'https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/COVID19_Canada.csv'
    prov = 'ON'
    
    dat <- get_data(prov, fname) %>%
        digest_data(test.pos = T)
    
    
    dat.all <- get_data(prov=NULL, fname) %>%
        digest_data(test.pos = TRUE)
    
    provs = c('ON', 'AB')
    n = length(provs)
    
    sapply(provs, extract_data_csv, test.pos = TRUE, fname=fname)
    
    g <- list()
    for(i in seq_along(provs)){
    g[[i]] <- get_data(provs[i], fname) %>%
        digest_data(test.pos = TRUE) %>%
        plot_data(title=provs[i])
    }
    h = 7
    pdf('plot-data-prov.pdf', height=h, width = 1.0*h*n)
    gridExtra::grid.arrange(grobs = g, ncol = n)
    dev.off()
}



