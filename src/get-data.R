source('read-data.R')

args <- commandArgs(trailingOnly = TRUE)
prov <- args[1]
# DEBUG: 
# prov = 'AB'   
# prov = 'ON'


fname = 'https://wzmli.github.io/COVID19-Canada/git_push/clean.Rout.csv'
# fname = 'https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/COVID19_Canada.csv'

dat <- prov %>% 
    extract_data_csv(test.pos = TRUE, fname = fname) %>%
    read.csv()

