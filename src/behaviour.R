source('utils-simple.R')


# --- PARAMS ----
# Retrieve parameters from CSV file

anchor.date = as_date("2020-03-06") #TODO: in a file? 
file.behav  = 'prm-behaviour.csv'
file.lags   = 'prm-lags.csv'

lag.inf.report = get_prm('lag_infection_to_report', file.lags)

d1 = as_date(get_prm('date_first_distancing', file.behav) )
t1 = as.numeric(d1 - anchor.date) + lag.inf.report

d_start_relax = as_date(get_prm('date_start_relax', file.behav) )
t.relax = as.numeric(d_start_relax - anchor.date) + lag.inf.report

d_backtonormal = as_date(get_prm('date_backtonormal', file.behav) )
t.bck2normal = as.numeric(d_backtonormal - anchor.date) + lag.inf.report



# --- SCENARIOS ----
# Scenario definitions

behave_NOTHING <- function(t, prm) {
    res = 1
    return(res)
}

B.NOTHING = list()


behave_BASELINE <- function(t, prm) {
    res = 1.0
    return(res)
}

B.BASELINE = list()



behave_ISO1 <- function(t, prm) {
    res = 1
    if(t >= prm$t1) res = prm$m1
    return(res)
}

B.ISO1 = list(t1 = t1, 
              m1 = 0.12345)  




behave_ISO2 <- function(t, prm) {
    
    t1 = prm$t1
    t2 = prm$t2
    t3 = prm$t3
    
    res = 1
    
    if(t1 <= t & t < t2 ) res = prm$m1
    if(t2 <= t & t < t3 ) res = (t-t3)/(t2-t3) * prm$m1 + (t-t2)/(t3-t2) * 1
    if(t3 <= t ) res = 1
    
    return(res)
}

B.ISO2 = list(t1 = t1, 
              t2 = t.relax,  
              t3 = t.bck2normal,
              m1 = 0.1234)  


# ---- DEBUG ----
if(0){
    x = 0:200
    y.ISO2 = sapply(x, behave_ISO2, B.ISO2)
    plot(x,y.ISO2, typ='l', lwd=6, las=1,ylim=c(0,1)) ; grid()
}
