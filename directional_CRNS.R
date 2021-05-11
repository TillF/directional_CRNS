# Illustrate theoretical discernability and accuracy of CRNS signal of a directional CRNS probe
# scripts for generating the figures of the manuscript
# "Assessing the feasibility of a directional CRNS-sensor for estimating soil moisture"
# by Till Francke, Maik Heistermann, Markus Köhli, Christian Budach, Martin Schrön, Sascha E. Oswald

# set the scenario of interest here:
  scenario="optimistic" #"optimistic" or "pessimistic"
  #scenario="pessimistic" #"optimistic" or "pessimistic"
# parameters ####
  library(latex2exp)
  
  R_1  = 2100 #base count rate of area_1 [counts per hour]
  # contrast = 0.5  #relative difference in countrate in area_2 [-] #for example figure CI
  contrast = 0.2  #relative difference in countrate in area_2 [-] #for example figure CI
  
  R_2 = R_1 * (1+contrast)
  

if (scenario == "optimistic")
{
 # snr = 2.63 #optimistic detector #signal-noise-ratio of targeted to untargeted area (from Markus' simulations)
  eps = 0.1 #optimistic ## fraction of non-epithermal neutrons in counts (=noise)
  gamma_dir = 0.2 #fraction of non-albedo counts in the directional signal
  eta = 0.72 #collimation efficiency
  beta = 0.4 #fraction of total counts after installing the collimator
  
} else
if (scenario == "pessimistic")
  {
  #snr = 1.57 #actual detector #signal-noise-ratio of targeted to untargeted area (from Markus' simulations)
  eps = 0.3 #pessimistic # fraction of non-epithermal neutrons in counts (=noise)
  gamma_dir = 0.31 #fraction of non-albedo counts in the directional signal
  eta = 0.61 #collimation efficiency
  beta = 0.3 #fraction of total counts after installing the collimator
} else
{ #idealistic scenario (testing only)
  eps = 0.0  #fraction of non-epithermal neutrons in counts (=noise)
  gamma_dir = 0. #fraction of non-albedo counts in the directional signal
  eta = 0.7 #collimation efficiency
  beta = 1 #fraction of total counts after installing the collimator
} 
  gamma_no = 0.17 #fraction of non-albedo counts in the non-directional signal (0.13 - 0.26)

  #A = rbind(c(snr,1), c(1, snr)) / (1+snr) #for converting N1 and N2 to Nf1 and Nf2
  A = rbind(c(eta,1-eta), c(1-eta, eta))  #for converting N1 and N2 to Nf1 and Nf2
  
  B =  A*(beta-beta*gamma_dir) + 1/2*beta*gamma_dir
  library(matlib)
  #Ai = inv(A) #for converting Nf1 and Nf2 to N1 and N2
  Bi = inv(B) #for converting N_f to N_total
  #using analytical formulation:
  #k_1 = 2 * beta * (2 * gamma_dir*eta - gamma_dir  - 2*eta + 1) 
  #k_2=(2 * gamma_dir * eta - gamma_dir- 2 * eta)
  #Bi = 1/k_1 * rbind(c(k_2, k_2+2), c(k_2+2, k_2))  #for converting N_total to N_f
  
  #c(2000,3000) %*% A
  #c(2000,3000) %*% B

delta_t = 1 # interval for aggregating the counts [h]
p_thresh = 0.05 #threshold p-value
precision = 0.05 #threshold for width of confidence interval of a rate, relative to true rate
#ratio_precision = 0.1 #threshold for width of confidence interval of estimated ratio between two rates, relative to true ratio of rates

#selection of ranges in plots

#count rates to analyse
R_coarse = c(450, 2100, 8000, 39000, 144000)
R_fine   = exp(seq(from=log(min(R_coarse)), to = log(max(R_coarse)), length.out = 100))

#aggregation intervals rates to analyse
delta_ts_coarse = c(1, 6, 12, 24, 48)
delta_ts_fine = exp(seq(from=log(min(delta_ts_coarse)), to=log(max(delta_ts_coarse+2)), length.out=200))  #aggregation intervals to test [h]

contrasts_coarse = c(.2, .5, .8, 1.4 )
contrasts_fine = exp(seq(from=log(min(contrasts_coarse)), to=log(max(contrasts_coarse)), length.out=100))  #contrasts to test [-]
                     
source("R/compute_ci_p.R")

# 0 compute the max/min requirements for the specified singular settings
  #searching to fulfill CI requirements
  #a) required rate for 24 h
  optimize(interval = range(R_coarse),        difference_to_threshold_ci, delta_t=24, contrast=contrast,          B=B, Bi=Bi, target_sensor="R_1r_t", precision=precision)
  #b) required min aggregation interval
  optimize(interval = range(delta_ts_coarse), difference_to_threshold_ci,             contrast=contrast, R_1=R_1, B=B, Bi=Bi, target_sensor="R_1r_t", precision=precision)
  #c) max contrast interval for 24 h
  optimize(interval = range(contrasts_coarse), difference_to_threshold_ci, delta_t=24,                    R_1=R_1, B=B, Bi=Bi, target_sensor="R_1r_t", precision=precision)
  
  #searching to fulfill p-value requirements
  #d) required rate for 24 h
  optimize(interval = range(R_coarse),        difference_to_threshold_p, delta_t=24, contrast=contrast,          B=B, Bi=Bi, target_comparison="N_1r_t.N_2r_t", p_thresh=p_thresh)
  #e) required min aggregation interval
  optimize(interval = range(delta_ts_coarse), difference_to_threshold_p,             contrast=contrast, R_1=R_1, B=B, Bi=Bi, target_comparison="N_1r_t.N_2r_t", p_thresh=p_thresh)
  #f) min contrast interval resolvable for 1 h resolution
  optimize(interval = range(contrasts_coarse), difference_to_threshold_p, delta_t=1 ,                    R_1=R_1, B=B, Bi=Bi, target_comparison="N_1r_t.N_2r_t", p_thresh=p_thresh)
  #g) min contrast interval resolvable for 24 h resolution
  optimize(interval = range(contrasts_coarse), difference_to_threshold_p, delta_t=24,                    R_1=R_1, B=B, Bi=Bi, target_comparison="N_1r_t.N_2r_t", p_thresh=p_thresh)
  


# 1 simulate theoretical distributions, considering sample uncertainty, function, for a single setting ####

  res2 = t(sapply(X=delta_ts_fine, FUN=compute_ci_p2, R_1=R_1, R_2=R_2, B=B, Bi=Bi, eps=eps))
  
  #res3 = t(sapply(X=3.5, FUN=compute_ci_p2, R_1=R_1, R_2=R_2, B=B, Bi=Bi, eps=eps))
  #compute_ci_p2(delta_t = 3.5, R_1=R_1, R_2=R_2, B=B, Bi=Bi, eps=eps)
  #sqrt(R_1 * (1-eps) * 3.5)
  
  #convert to match conventions of plotting routine
  tt = c(R_1, R_2, c(R_1, R_2) %*% B) #all four count rates
  #R = t(array(data=tt, c(4, length(delta_ts_fine)), dimnames = list(c("R_1", "R_2", "R_f1", "R_f2"), NULL))) #count rates [counts per hour]
  R = t(array(data=tt, dimnames = list(c("R_1", "R_2", "R_f1", "R_f2")))) #count rates [counts per hour]
  R = c( R, R * (1-eps)) #these are the total ("_t") and epithermal ("_epi") rates
  names(R) = paste0(c("R_1", "R_2", "R_f1", "R_f2"), c(rep("_t",4), rep("_epi",4)))
  
  ci_cols = grep(dimnames(res2)[[2]], pattern = "R_", value = TRUE) #get the names of the rates contained
  ci_lohi =   array(NA, c(length(delta_ts_fine), length(ci_cols), 2), dimnames = list(NULL, ci_cols, c("lo", "hi"))) #confidence intervals of measurements
  
  ci_lohi[,,"lo"] = res2[, ci_cols]
  ci_lohi[,,"hi"] = t(2 *R[sub(ci_cols, pattern = "(\\d)r_", repl="\\1_")] - t(res2[, ci_cols])  ) #just the symetric other side
  p_value_comp    = res2[, !(dimnames(res2)[[2]] %in% ci_cols)]
  #str(res2)
  

# 1a plot counts and CIs (fig. 7 / ci_example) ####
  x11()
  par(cex=1.2)
  library(RColorBrewer)
  n = 8
  #construct palette
  palette(brewer.pal(n=n, "Dark2"))
  
  #CI
  #cols2plot  = 1:NCOL(ci_lohi)
  #ncols2plot = NCOL(ci_lohi)
  
  #plot CI as lines
  #cols2plot =  c("R_1_epi", "R_2_epi", "R_f1_epi", "R_f2_epi",  "R_1r_epi",  "R_2r_epi") #columns containing the "_epi" values
  #cols2plot =  grep(dimnames(ci_lohi)[[2]], pattern = "_t$", value = TRUE) #columns containing the "_t" (total) values
  cols2plot =  c("R_1_t", "R_2_t", "R_f1_t", "R_f2_t",  "R_1r_t",  "R_2r_t") #columns containing the "_t" values
  
  ncols2plot = length(cols2plot)
  
  ylims = pmax(0, range(ci_lohi[-(1:5), cols2plot,])) #discard the very wide interval in the first three rows
  ylims[2] = max(ylims[2], max (R)) #accommodate largest rate in plot
  #ylims = c(1800, 2500)
  
  plot(1, type="n",
       ylim=ylims, xlim=range(delta_ts_fine), ylab = "count rate R [counts/h]", xlab=expression(paste("aggregation interval ", Delta, "t [h] ")), log="x")

  #plot counts
  abline(h= R[  1:4 ], col=1:4, lty=1, lwd=2) #the "_t" rates
  #abline(h= R[-(1:4)], col=1:4, lty="dashed", lwd=2) #the "_epi" rates
  
  
  lwd = 1
  lty = ifelse(grepl(cols2plot, pattern = "\\dr_"), "dotted", "dotdash")
  
  matplot(delta_ts_fine, y= ci_lohi[, cols2plot,"hi"], type="l", add = TRUE,
          col=1:(ncols2plot-2),lty=lty, lwd=lwd)
  matplot(delta_ts_fine, y= ci_lohi[, cols2plot,"lo"], type="l", add = TRUE,
          col=1:(ncols2plot-2), lty=lty, lwd=lwd)
  
  
  #plot CIs as envelopes
  org_cols = palette() #store original palette
  transp_cols <- adjustcolor(org_cols, alpha.f = 0.3) #set transparency (for plotting envelopes)
  palette(transp_cols)
  
  for (i in 1:ncols2plot) #plot envelopes for increased visibility
    polygon(c(delta_ts_fine, rev(delta_ts_fine)), 
            c(ci_lohi[,cols2plot[i],"hi"],rev(ci_lohi[,cols2plot[i],"lo"])), border = NA, col = i)
  palette(org_cols)
  
  ci_range = apply(ci_lohi[, cols2plot,], MAR=c(1,2), diff)   
  # head(ci_range / R[ref_cols])
  # (R["R_f1_epi"] - ci_lohi[1, "R_f1_epi","lo"]) / R["R_f1_epi"] #ok
  # (R["R_f2_epi"] - ci_lohi[1, "R_f2_epi","lo"]) / R["R_f2_epi"]
  # 
  # (ci_lohi[1, "R_f1_epi","hi"] - ci_lohi[1, "R_f1_epi","lo"]) / R["R_f1_epi"] #
  # (ci_lohi[1, "R_f2_epi","hi"] - ci_lohi[1, "R_f2_epi","lo"]) / R["R_f2_epi"] # ok
  # head(t(t(ci_range) / R[ref_cols]))
  ref_cols = sub(cols2plot, pattern = "(\\d)r_", repl="\\1_") #the columns containing the respective count rates
  
  within_precision = t(t(ci_range) / R[ref_cols]) <= precision #check, where the CI is smaller than the chosen precision
  min_delta_t = apply(within_precision, MAR=2, FUN=function(x){min(which(x))}) #get indices to times when error falls below threshold
  for (i in 1:ncols2plot) #plot lines denoting minimum aggregation interval
  {  
    xplot = delta_ts_fine[min_delta_t[i]] #find position of x-axis where to plot the line
    lines(rep(xplot,2),  ci_lohi[min_delta_t[i], cols2plot[i],], col=i) #plot line illustrating the width of the interval
    lines(rep(xplot,2),  c(ci_lohi[min_delta_t[i], cols2plot[i],"lo"], 0), col=i, lty="dotted")  #plot line pointing down to the x-axis
    #points(xplot,  R[i], col=i, pch=25, crt=180 ) 
    points(xplot,  ylims[1], col=i, pch=25, bg=i) #plot marker at x-axis
    }          
  
  legend_str = sapply(FUN = TeX, X =  c("R_1", "R_2", "R_{f1}", "R_{f2}", "CI(95%) R","CI(95%), R_r"))
  legend("topright", legend=legend_str,
         lty=c(rep("solid",4), "dashed", "dotdash", "dotted"),
         col=c(1:4,rep("black",2)), 
         lwd=c(rep(2,4),rep(1,2)), ncol=1, bg="white")
  
  savePlot(filename = paste0("ci_examp_", scenario,".png"), type = "png")
  
  if (Sys.info()["sysname"]=="Windows") #pdf-export only works in Windows
    savePlot(filename = paste0("ci_examp_", scenario,".pdf"), type = "pdf")
  
  
# 1b p-values of comparison (fig. 8 / p_example)####
  x11()
  par(cex=1.2)
  #cols2plot =  dimnames(p_value_comp)[[2]]
  #cols2plot = c("N_1_epi.N_2_epi", "N_1r_epi.N_2r_epi", "N_f1_epi.N_f2_epi") #, "N_f1_t.N_f2_t"
  cols2plot = c("N_1_t.N_2_t", "N_1r_t.N_2r_t", "N_f1_t.N_f2_t") #, "N_f1_t.N_f2_t"
  
  ncols2plot = length(cols2plot)
  #ylims = range(p_value_comp[-(1:5), cols2plot]) 
    #ylims = c(1800, 2500)
  
  library(RColorBrewer)
  palette(brewer.pal(n=ncols2plot, "Dark2")) #construct palette
  
  ix_min_below_thresh = apply(p_value_comp, MAR=2, FUN=function(x){which(x <= p_thresh)[1]}) #find position, beyond which the p-value is lower than the specified threshold
  
  delta_t_subs = 1:min(max(ix_min_below_thresh, length(delta_ts_fine)/2), length(delta_ts_fine))
  
  matplot(delta_ts_fine, y= p_value_comp[, cols2plot], type="l",
          ylim=c(0, 0.1), xlim=range(delta_ts_fine[delta_t_subs]), ylab = "p-value of comparison [-]", xlab=expression(paste("aggregation time ", Delta, "t [h]")),
          col=1:ncols2plot, lty=1, lwd=2, log="x")
  abline(h=p_thresh, lty="dashed")
  
  #find indices, after which the p-value is below the threshold
  #abline(v=delta_ts_fine[ix_min_below_thresh[cols2plot]], col=1:ncols2plot, lty="dotted")
  for (i in 1:ncols2plot)
  {  
    xplot = delta_ts_fine[ix_min_below_thresh[cols2plot[i]]]
    lines(rep(xplot,2), c(-0.002, p_thresh), col=i, lty="dotted")
    points(xplot,  -0.002, col=i, pch=25, bg=i) 
  }
  
  cols2plot_l = gsub(x=cols2plot  , pattern="_epi", repl=",epi")
  cols2plot_l = gsub(x=cols2plot_l, pattern="_t", repl=",total")
  cols2plot_l = gsub(x=cols2plot_l, pattern="_([^\\.]*)", repl="_{\\1}")
  cols2plot_l = gsub(x=cols2plot_l, pattern="\\.", repl=" vs. ")
  cols2plot_l = gsub(x=cols2plot_l, pattern="N", repl="R")
  
  legend_str = sapply(FUN = TeX, X =  c(cols2plot_l, "p_{thresh} = 0.05",  paste0("$\\Delta t$ | p < p_{thresh}")))
  legend("topright", legend=legend_str, 
         lty=c(rep("solid", ncols2plot), "dashed", "dotted"),
         lwd = c(rep(2, ncols2plot), 1, 1) ,
        col=c(palette()[1:ncols2plot],"black","black"),
         bg="white" )
  savePlot(filename =  paste0("p_value_examp_", scenario,".png"), type = "png")

  if (Sys.info()["sysname"]=="Windows") #pdf-export only works in Windows
    savePlot(filename = paste0("p_value_examp_", scenario,".pdf"), type = "pdf")
  

# A1 simulate theoretical distributions for a range of settings: R vs. t_min (CI-thresh) ####
  # requires functions p_value_two_sided() and compute_ci_p2()
  
  #library(scanningCRNS)
  
delta_ts_fine = exp(seq(from=log(0.5), to=log(30), length.out=100))  #aggregation intervals to test [h]


# find the minimum required aggregation time for the CI to fall below specified precision
find_min_aggr_time_ci = function(R_1, contrast,  B,  Bi, target_sensor, precision) 
{  
  res = optimize(interval = c(0.0, 1000), difference_to_threshold_ci, R_1=R_1, contrast=contrast, B=B, Bi=Bi, target_sensor=target_sensor, precision=precision)
  return(res$minimum)
}

#
#difference_to_threshold_ci(delta_t, R_1, contrast=contrast,  B,  Bi, target_sensor="R_1r_epi", precision) 
#compute_ci_p2(delta_t = delta_t, R_1=R_1, R_2=R_1*(1+con), B=B, Bi=Bi, eps=eps)  
  
collected_min_aggr_time = array(NA, dim=c(length(contrasts_coarse), length(R_fine)))
for (i in 1:length(contrasts_coarse))
{
  contrast = contrasts_coarse[i]
  collected_min_aggr_time[i,] = sapply(X=R_fine, FUN=find_min_aggr_time_ci, contrast=contrast, B=B, Bi=Bi, target_sensor="R_1r_t", precision=precision, simplify = TRUE)
}



# A1 plot results ####
  library(RColorBrewer)
  palette(brewer.pal(n=nrow(collected_min_aggr_time), "Dark2")) #construct palette
  x11()
  par(mar=par("mar") + c(0, 0.5,0,0), cex=1.2)
  #ylims=c(0, max(collected_min_aggr_time))
  ylims=c(1,100)
  
  options(scipen=5) #prevent scientific notation in log-scale of plot
  matplot(R_fine, y= t(collected_min_aggr_time), type="l", log="x",
          ylim=ylims,
           xlab=TeX("lower countrate R_{1,total} [counts/h]"), 
           ylab=TeX("minimum aggregation time  $\\Delta t_{min}^{determ}$ [h]"),
          #main=scenario,
          col=1:4, lty=1, lwd=2)

  text(x = exp(mean(log(range(R_fine)))), y=0.9 * max(ylims), labels = paste0("\"", scenario, "\""), cex = 1.2)
              
  legend_str = sapply(FUN = TeX, X =  c("$\\Delta R$=", paste0(contrasts_coarse, "")))
  legend("topright", legend=legend_str, col=c(NA, 1:length(contrasts_coarse)), lty=1, lwd=2)
  savePlot(filename =  paste0("A1_CI_R_vs_aggrtime_", scenario,".png"), type = "png")
  if (Sys.info()["sysname"]=="Windows") #pdf-export only works in Windows
    savePlot(filename = paste0("A1_CI_R_vs_aggrtime_", scenario,".pdf"), type = "pdf")

# A2 simulate theoretical distributions for a range of settings: R vs. t_min (p-thresh) ####
 # requires functions p_value_two_sided() and compute_ci_p2()
library(scanningCRNS)
  
#delta_ts_fine = exp(seq(from=log(0.5), to=log(30), length.out=100))  #aggregation intervals to test [h]




# find the minimum required aggregation time for the p-value to fall below specified p_thresh
find_min_aggr_time_p = function(R_1, contrast,  B,  Bi, target_comparison, p_thresh) 
{  
  minv = optimize(interval = c(0.01, 48), difference_to_threshold_p, R_1=R_1, contrast=contrast, B=B, Bi=Bi, target_comparison="N_1r_t.N_2r_t", p_thresh=p_thresh)$minimum
  return(minv)
}

#minimum aggregation time required for example setting 
#(delta_t = find_min_aggr_time_p(R_1 = R_1, contrast=contrast, B=B, Bi=Bi, target_comparison="N_1r_t.N_2r_t", p_thresh=p_thresh))
#difference_to_threshold_p(delta_t = 1, R_1=R_1, contrast=contrast, B=B, Bi=Bi, target_comparison="N_1r_t.N_2r_t", p_thresh=p_thresh)
#compute_ci_p2(delta_t = 7.97, R_1=R_1, R_2=R_2, B=B, Bi=Bi)  

collected_min_aggr_time = array(NA, dim=c(length(contrasts_coarse), length(R_fine)))

for (i in 1:length(contrasts_coarse))
{
  contrast = contrasts_coarse[i]
  collected_min_aggr_time[i,] = sapply(X=R_fine, FUN=find_min_aggr_time_p, contrast=contrast, B=B, Bi=Bi, target_comparison="N_1r_t.N_2r_t", p_thresh=p_thresh, simplify = TRUE)
}


# A2 plot results ####
library(RColorBrewer)
palette(brewer.pal(n=nrow(collected_min_aggr_time), "Dark2")) #construct palette

x11()
par(mar=par("mar") + c(0, 0.5,0,0), cex=1.2)
matplot(R_fine, y= t(collected_min_aggr_time), type="l", log="x",
        xlab=TeX("lower countrate R_{1,total} [counts/h]"), 
        ylab=TeX("minimum aggregation time $ $ $ \\Delta t_{min}^{disting}$ [h]"),
#        main=scenario,
                col=1:4, lty=1, lwd=2)


text(x = exp(mean(log(range(R_fine)))), y=0.9 * max(collected_min_aggr_time), labels = paste0("\"", scenario, "\""), cex = 1.2)

legend_str = sapply(FUN = TeX, X =  c("$\\Delta R$=", paste0(contrasts_coarse, "")))
legend("topright", legend=legend_str, col=c(NA, 1:length(contrasts_coarse)), lty=1, lwd=2)
savePlot(filename = paste0("A2_p_R_vs_aggrtime_", scenario,".png"), type = "png")

if (Sys.info()["sysname"]=="Windows") #pdf-export only works in Windows
  savePlot(filename = paste0("A2_p_R_vs_aggrtime_", scenario,".pdf"), type = "pdf")

# B1 simulate theoretical distributions for a range of settings: R vs. contrasts (CI-thresh) ####
# requires functions p_value_two_sided() and compute_ci_p2()

# find the minimum required contrast for the CI to fall below specified precision
find_max_contrast = function(R_1, delta_t,  B,  Bi, target_sensor, precision) 
{  
  res = optimize(interval = c(0.0, 10), difference_to_threshold_ci, R_1=R_1, delta_t=delta_t, B=B, Bi=Bi, target_sensor, precision)
  
  #optimize(interval = range(contrasts_coarse), difference_to_threshold_ci, R_1=2000, R_2=2000+100, B=B, Bi=Bi, target_sensor, precision)
  if (res$objective > precision/100) res$minimum = NA #if the requested precision cannot be reached, return NA
  
  return(res$minimum)
}
#max contrast allowed for example setting 
#find_max_contrast(delta_t=0.1, R_1 = R_1, B=B, Bi=Bi, target_sensor="R_1r_t", precision=precision) #test

  
collected_max_contrast = array(NA, dim=c(length(delta_ts_coarse), length(R_fine)))

for (i in 1:length(delta_ts_coarse))
{
  delta_t = delta_ts_coarse[i]
  collected_max_contrast[i,] = sapply(X=R_fine, FUN=find_max_contrast, delta_t=delta_t, B=B, Bi=Bi, target_sensor="R_1r_t", precision=precision, simplify = TRUE)
}
# R_1=5000
# con=0.1
# #delta_t=9
# delta_t=6.350596
# 
# (con = find_max_contrast(delta_t=delta_t, R_1 = R_1, B=B, Bi=Bi, target_sensor="R_1r_epi", precision=precision)) #test
# difference_to_threshold_ci(delta_t, R_1, contrast=con,  B,  Bi, target_sensor="R_1r_epi", precision) 
# compute_ci_p2(delta_t = delta_t, R_1=R_1, R_2=R_1*(1+con), B=B, Bi=Bi)  

# B1 plot results ####
library(RColorBrewer)
palette(brewer.pal(n=nrow(collected_max_contrast), "Paired")) #construct palette

x11()
par(mar=par("mar") + c(0, 0.5,0,0), cex=1.2)
options(scipen=5) #prevent scientific notation in log-scale of plot
ylims=c(0, 2*max(contrasts_coarse))

matplot(R_fine, y= t(collected_max_contrast), type="l",
          xlab=TeX("lower countrate R_{1,epi} [counts/h]"), 
          ylab=TeX("max contrast $ $ $\\Delta R$ [-]"),
          col=1:length(delta_ts_coarse), lty=1, lwd=2, log="x",
          ylim=ylims
        #, main=scenario
        )

text(x = exp(mean(log(range(R_fine)))), y=0.9 * max(ylims), labels = paste0("\"", scenario, "\""), cex = 1.2)

legend_str = sapply(FUN = TeX, X =  c("$\\Delta t$=", paste0(delta_ts_coarse, " h")))
legend("topright", legend=legend_str, 
         col=c(NA, 1:length(delta_ts_coarse)), lty=1, lwd=2, bg = "white")
savePlot(filename = paste0("B1_ci_R_vs_contrast_", scenario,".png"), type = "png")
if (Sys.info()["sysname"]=="Windows") #pdf-export only works in Windows
  savePlot(filename = paste0("B1_ci_R_vs_contrast_", scenario,".pdf"), type = "pdf")

# B2 simulate theoretical distributions for a range of settings: R vs. contrasts (p-thresh) ####
# requires functions p_value_two_sided() and compute_ci_p()

# find the minimum required aggregation time for the p-value to fall below specified p_thresh
find_min_contrast = function(R_1, delta_t,  B,  Bi, target_comparison, p_thresh) 
{  
  res = optimize(interval = c(0.001, 100), difference_to_threshold_p, R_1=R_1, delta_t=delta_t, B=B, Bi=Bi, target_comparison, p_thresh)
  if (res$objective > p_thresh/100) res$minimum = NA #if the requested p_value cannot be reached, return NA
  
  return(res$minimum)
}
# R_1=5000
# con=0.1
# #delta_t=9
# delta_t=6.350596
# 
# (con = find_min_contrast (R_1, delta_t,  B,  Bi, target_comparison="N_1r_epi.N_2r_epi", p_thresh) )
# difference_to_threshold_p(delta_t, R_1, con,  B,  Bi,  target_comparison="N_1r_epi.N_2r_epi", p_thresh) 
# compute_ci_p2(delta_t = delta_t, R_1=R_1, R_2=R_1*(1+con), B=B, Bi=Bi, eps=eps)  

  
collected_min_contrast = array(NA, dim=c(length(delta_ts_coarse), length(R_fine)))

for (i in 1:length(delta_ts_coarse))
{
  delta_t = delta_ts_coarse[i]
  collected_min_contrast[i,] = sapply(X=R_fine, FUN=find_min_contrast, delta_t=delta_t, B=B, Bi=Bi, target_comparison="N_1r_t.N_2r_t", p_thresh=p_thresh, simplify = TRUE)
}

#find_min_contrast (R_1, delta_t,  B,  Bi, target_comparison="N_1r_t.N_2r_t", p_thresh) 


# B2 plot results ####
library(RColorBrewer)
palette(brewer.pal(n=nrow(collected_min_contrast), "Paired")) #construct palette

x11()
par(mar=par("mar") + c(0, 0.5,0,0), cex=1.2)
ylims = c(0, 0.5*max(contrasts_coarse))
matplot(R_fine, y= t(collected_min_contrast), type="l",
        xlab=TeX("lower countrate R_{1} [counts/h]"), 
        ylab=TeX("min contrast $ $ $\\Delta R$ [-]"),
        col=1:length(delta_ts_coarse), lty=1, lwd=2 , log="x",
        ylim=ylims
        #, main=scenario
        )

text(x = exp(mean(log(range(R_fine)))), y=0.9 * max(ylims), labels = paste0("\"", scenario, "\""), cex = 1.2)


legend_str = sapply(FUN = TeX, X =  c("$\\Delta t$=", paste0(delta_ts_coarse, " h")))
legend("topright", legend=legend_str, col=c(NA, 1:length(delta_ts_coarse)), lty=1, lwd=2)
savePlot(filename = paste0("B2_p_R_vs_contrast_", scenario,".png"), type = "png")
if (Sys.info()["sysname"]=="Windows") #pdf-export only works in Windows
  savePlot(filename = paste0("B2_p_R_vs_contrast_", scenario,".pdf"), type = "pdf")



# C1 simulate theoretical distributions for a range of settings: aggr_time vs. contrasts (CI-thresh) ####
# requires functions p_value_two_sided() and compute_ci_p2()

# find the minimum required aggregation time for the CI to fall below specified precision
delta_ts_fine2 = exp(seq(from=log(min(delta_ts_coarse)), to=log(2*max(delta_ts_coarse)), length.out=200))  #aggregation intervals to test [h]
#here, we need some finer resolution to produce continuous lines


find_max_contrast = function(R_1, delta_t,  B,  Bi, target_sensor, precision) 
{  
  res = optimize(interval = c(0, 15), difference_to_threshold_ci, delta_t=delta_t, R_1=R_1, B=B, Bi=Bi, target_sensor=target_sensor, precision=precision)
  
  if (res$objective > precision/100) res$minimum = NA #if the requested precision cannot be reached, return NA
  
  return(res$minimum)
}
#find_max_contrast(delta_t = 1, R_1 = R_coarse[1], B=B, Bi=Bi, target_sensor="R_1r_t", precision=precision) #test


#required rate for 24 h and low contrast
optimize(interval = range(R_coarse), difference_to_threshold_ci, delta_t=24, contrast=min(contrasts_coarse), B=B, Bi=Bi, target_sensor="R_1r_t", precision=precision)
#required rate for 24 h and high contrast
optimize(interval = range(R_coarse), difference_to_threshold_ci, delta_t=24, contrast=max(contrasts_coarse), B=B, Bi=Bi, target_sensor="R_1r_t", precision=precision)




collected_max_contrast =array(NA, dim=c(length(R_coarse), length(delta_ts_fine2)))

for (i in 1:length(R_coarse))
{
  R_1 = R_coarse[i]
  collected_max_contrast[i,] = sapply(X=delta_ts_fine2, FUN=find_max_contrast, R_1=R_1, B=B, Bi=Bi, target_sensor="R_1r_t", precision=precision, simplify = TRUE)
}


# C1 plot results ####
library(RColorBrewer)
palette(brewer.pal(n=nrow(collected_max_contrast), "Paired")) #construct palette

x11()
par(mar=par("mar") + c(0, 0.5,0,0), cex=1.2)
options(scipen=5) #prevent scientific notation in log-scale of plot
ylims = c(0.5*min(contrasts_coarse), 2*max(contrasts_coarse))
matplot(delta_ts_fine2, y= t(collected_max_contrast), type="l",
        xlab=TeX("aggregation time  $\\Delta t$ [h]"),
        ylab=TeX("max contrast $ $ $\\Delta R$ [-]"),
        
        col=1:length(R_coarse), lty=1, lwd=2, log="xy",
        ylim=ylims
               #, main=scenario
        )

text(x = exp(mean(log(range(delta_ts_fine2)))), y=0.9 * max(ylims), labels = paste0("\"", scenario, "\""), cex = 1.2)

legend_str = sapply(FUN = TeX, X =  c("R_1 [counts/h]=", paste0(R_coarse, "")))
legend("bottomright", legend=legend_str, 
       col=c(NA, 1:length(R_coarse)), lty=1, lwd=2, bg = "white")

savePlot(filename = paste0("C1_ci_aggrtime_vs_contrast_", scenario,".png"), type = "png")
if (Sys.info()["sysname"]=="Windows") #pdf-export only works in Windows
  savePlot(filename = paste0("C1_ci_aggrtime_vs_contrast_", scenario,".pdf"), type = "pdf")


# C2 simulate theoretical distributions for a range of settings: R vs. contrasts (p-thresh) ####
# requires functions p_value_two_sided() and compute_ci_p()

# find the minimum required contrast for the p-value to fall below specified p_thresh
find_min_contrast = function(R_1, delta_t,  B,  Bi, target_comparison, p_thresh) 
{  
  res = optimize(interval = c(0.001, 10), difference_to_threshold_p, R_1=R_1, delta_t=delta_t, B=B, Bi=Bi, target_comparison, p_thresh)
  return(res$minimum)
}


collected_min_contrast = array(NA, dim=c(length(R_coarse), length(delta_ts_fine)))

for (i in 1:length(R_coarse))
{
  R_1 = R_coarse[i]
  collected_min_contrast[i,] = sapply(X=delta_ts_fine, FUN=find_min_contrast, R_1=R_1, B=B, Bi=Bi, target_comparison="N_1r_t.N_2r_t", p_thresh=p_thresh, simplify = TRUE)
}


#find_min_contrast(delta_ts_fine[i], R_1=R_1, B=B, Bi=Bi, target_comparison="N_1r_t.N_2r_t", p_thresh=p_thresh)


#min required rate for 1 h and low contrast
optimize(interval = c(1,100000), difference_to_threshold_p, delta_t=1, contrast=min(contrasts_coarse), B=B, Bi=Bi, target_comparison="N_1r_t.N_2r_t", p_thresh=p_thresh)
#required rate for 1 h and high contrast
optimize(interval = c(1,100000), difference_to_threshold_p, delta_t=1, contrast=max(contrasts_coarse), B=B, Bi=Bi, target_comparison="N_1r_t.N_2r_t", p_thresh=p_thresh)
#required rate for 1 h and contrast 0.4
optimize(interval = c(1,100000), difference_to_threshold_p, delta_t=1, contrast=0.4, B=B, Bi=Bi, target_comparison="N_1r_t.N_2r_t", p_thresh=p_thresh)



# C2 plot results ####
library(RColorBrewer)
palette(brewer.pal(n=nrow(collected_min_contrast), "Paired")) #construct palette

x11()
par(mar=par("mar") + c(0, 0.5,0,0), cex=1.2)
options(scipen=5) #prevent scientific notation in log-scale of plot
ylims = c(1e-3, 0.5*max(contrasts_coarse))
matplot(delta_ts_fine, y= t(collected_min_contrast), type="l",
        xlab=TeX("aggregation time  $\\Delta t$ [h]"),
        ylab=TeX("min contrast $ $ $\\Delta R$ [-]"),
        col=1:length(R_coarse), lty=1, lwd=2, log="xy",
        ylim=ylims
        #, main=scenario
)

text(x = exp(mean(log(range(delta_ts_fine)))), y=0.9 * max(ylims), labels = paste0("\"", scenario, "\""), cex = 1.2)

legend_str = sapply(FUN = TeX, X =  c("R_1 [counts/h]=", paste0(R_coarse, "")))
legend("topright", legend=legend_str, 
       col=c(NA, 1:length(R_coarse)), lty=1, lwd=2, bg = "white")
savePlot(filename = paste0("C2_p_aggrtime_vs_contrast_", scenario,".png"), type = "png")
if (Sys.info()["sysname"]=="Windows") #pdf-export only works in Windows
  savePlot(filename = paste0("C2_p_aggrtime_vs_contrast_", scenario,".pdf"), type = "pdf")  