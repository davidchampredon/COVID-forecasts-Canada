### 
###  ANALYZE FORECAST OUTPUTS
###

source('analyze-fcts.R')

message("Starting analysis ...")

load('digest.RData')

# ---- Thresholds ----

thresh.cum  <- list(
    Canada = c(hosp     = 150,
               critical = 50,
               death    = 30) * 1e3,
    
    ON = c(hosp     = 50,
           critical = 20,
           death    = 15) * 1e3,
    
    QC = c(hosp     = 30,
           critical = 20,
           death    = 4) * 1e3,
    
    AB = c(hosp     = 10,
           critical = 6,
           death    = 1) * 1e3,
    
    BC = c(hosp     = 0.5,
           critical = 0.25,
           death    = 0.2) * 1e3,
    
    MB = c(hosp     = 2,
           critical = 1,
           death    = 0.25) * 1e3,
    
    SK = c(hosp     = 0.2,
           critical = 0.1,
           death    = 0.1) * 1e3,
    
    NL = c(hosp     = 0.1,
           critical = 0.1,
           death    = 0.1) * 1e3,
    
    NS = c(hosp     = 0.2,
           critical = 0.1,
           death    = 0.1) * 1e3
    
    
)

thresh.pkv  <- list(
    Canada = c(hosp     = 4000,
               critical = 1500,
               death    = 200),
    
    ON = c(hosp     = 2000,
           critical = 500,
           death    = 250),
    
    QC = c(hosp     = 2000,
           critical = 1000,
           death    = 100),
    
    AB = c(hosp     = 200,
           critical = 100,
           death    = 10),
    
    BC = c(hosp     = 200,
           critical = 100,
           death    = 10),
    
    MB = c(hosp     = 50,
           critical = 25,
           death    = 10),
    
    SK = c(hosp     = 50,
           critical = 25,
           death    = 10),
    
    NL = c(hosp     = 25,
           critical = 10,
           death    = 10),
    
    NS = c(hosp     = 25,
           critical = 10,
           death    = 10)
)


cumdeath.thresh.val = c(
    Canada = 2500,
    AB = 50,
    BC = 50,
    MB = 10,
    NL = 10,
    NS = 10,
    ON = 1000,
    QC = 1600,
    SK = 10
)



targetdate = c(
    Canada = as_date('2020-05-15'),
    AB = as_date('2020-05-15'),
    BC = as_date('2020-05-15'),
    MB = as_date('2020-05-15'),
    NL = as_date('2020-05-15'),
    NS = as_date('2020-05-15'),
    ON = as_date('2020-05-15'),
    QC = as_date('2020-05-15'),
    SK = as_date('2020-05-15')
)


ci <- c(0.8, 0.95)
pdf.w = 12
pdf.h = 8
col.scen <- c('#505F90', 
              '#D4CD6A',
              '#D48E6A')

# ---- Loop Prov ----

prov.vec = c(unique(as.character(sims$Province)), 'Canada')

# Removing problems:
prov.vec = prov.vec[prov.vec != 'NL']
prov.vec = prov.vec[prov.vec != 'NS']

for( prov in prov.vec ){  
    
    # prov = 'ON'
    # prov = 'Canada'
    
    msg <- paste("=====",prov,"====")
    print("")
    print(msg)
    print("")
    
    pdf(paste0(dir.plot,'plot-analyze-pkv-',prov,'.pdf'), width = pdf.w, height=pdf.h)
    if(prov=='Canada') tmp = digest.peaks.canada
    if(prov!='Canada') tmp = digest.peaks
    g <- compare_intrv(prov = prov, 
                       digest.df = tmp,
                       varname = 'pkv',
                       outcometype = 'Daily Peak',
                       thresh = thresh.pkv[[prov]],
                       col.scen = col.scen,
                       ci = ci)
    grid.arrange(g$val, g$proba, nrow=2)
    dev.off()
    
    
    pdf(paste0(dir.plot,'plot-analyze-cum-',prov,'.pdf'), width = pdf.w, height=pdf.h)
    if(prov=='Canada') tmp = digest.cum.canada
    if(prov!='Canada') tmp = digest.cum
    g <- compare_intrv(prov = prov, 
                       digest.df = tmp,
                       varname = 's',
                       outcometype = 'Cumul. Num.',
                       thresh = thresh.cum[[prov]],
                       col.scen = col.scen,
                       ci = ci)
    grid.arrange(g$val, g$proba, nrow=2)
    dev.off()
    
    
    pdf(paste0(dir.plot,'plot-analyze-pkt-',prov,'.pdf'), width = pdf.w, height=pdf.h)
    if(prov=='Canada') tmp = digest.peaks.canada
    if(prov!='Canada') tmp = digest.peaks
    g <- compare_pkt_intrv(prov = prov, 
                           digest.peaks = tmp, 
                           targetdate = targetdate[[prov]], 
                           col.scen = col.scen)
    grid.arrange(g$val, g$proba, nrow=2)
    dev.off()
    
    pdf(paste0(dir.plot,'plot-analyze-traj-',prov,'.pdf'),width = pdf.w, height=pdf.h)
    if(prov=='Canada') tmp = digest.ss.canada
    if(prov!='Canada') tmp = digest.ss
    plot_compare_traj_intrv(prov, tmp)
    dev.off()
    
    pdf(paste0(dir.plot,'plot-analyze-probaHigher-',prov,'.pdf'),
        width = pdf.w, height=pdf.h/2)
    if(prov=='Canada') tmp = sims.canada
    if(prov!='Canada') tmp = sims
    plot_proba_larger(prov, sims = tmp, dat.all, n.days = 7, col.scen)
    dev.off()
    
    pdf(paste0(dir.plot,'plot-cumdeath-shortterm-',prov,'.pdf'),
        width = pdf.w/1.5, height=pdf.h/1.5)
    if(prov=='Canada') tmp = sims.canada
    if(prov!='Canada') tmp = sims
    proba_cumdeath_above(prov, sims = tmp, 
                         thresh.value = cumdeath.thresh.val[[prov]], 
                         thresh.date = as_date('2020-05-15'))
    dev.off()
}

pdf(paste0(dir.plot,'plot-effect-lambda.pdf'),
    width = pdf.w/1.5, height=pdf.h/2)
plot_abc_post_prm(abc.post, 
                 prm.name = 'lambda', 
                 scenario = 'ISO1', 
                 ci = c(0.5, 0.95),
                 title = 'Impact of First Social Distancing',
                 one.minus = TRUE)
dev.off()



# ------ CANADA

prov = 'Canada'

pdf(paste0(dir.plot,'plot-analyze-traj-Canada.pdf'),width = pdf.w, height=pdf.h)
plot_compare_traj_intrv(prov = 'Canada', 
                        digest.ss = digest.ss.canada)
dev.off()

pdf(paste0(dir.plot,'plot-analyze-cum-Canada.pdf'), width = pdf.w, height=pdf.h)
g <- compare_intrv(prov = prov, 
                   digest.df = digest.cum.canada,
                   varname = 's',
                   outcometype = 'Cumul. Num.',
                   thresh = thresh.cum[[prov]],
                   col.scen = col.scen,
                   ci = ci)
grid.arrange(g$val, g$proba, nrow=2)
dev.off()

pdf(paste0(dir.plot,'plot-analyze-pkt-',prov,'.pdf'), width = pdf.w, height=pdf.h)
g <- compare_pkt_intrv(prov = prov, 
                       digest.peaks = digest.peaks.canada, 
                       targetdate = targetdate[[prov]], 
                       col.scen = col.scen)
grid.arrange(g$val, g$proba, nrow=2)
dev.off()

pdf(paste0(dir.plot,'plot-analyze-pkv-',prov,'.pdf'), width = pdf.w, height=pdf.h)
g <- compare_intrv(prov = prov, 
                   digest.df = digest.peaks.canada,
                   varname = 'pkv',
                   outcometype = 'Daily Peak',
                   thresh = thresh.pkv[[prov]],
                   col.scen = col.scen,
                   ci = ci)
grid.arrange(g$val, g$proba, nrow=2)
dev.off()


