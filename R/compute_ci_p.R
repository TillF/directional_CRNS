 #https://statisticsbyjim.com/hypothesis-testing/confidence-intervals-compare-means/
  #http://www.stat.yale.edu/Courses/1997-98/101/meancomp.htm
  p_value_two_sided = function(mu1, mu2, se1, se2)  #compute p-value for two-sided comparison of two means with the given standard errors
  {
    se_diff = sqrt(se1^2 + se2^2) #get standard error of difference by addition of variance of standard errors
    pv = 2*pnorm(q = 0, mean = abs(mu1-mu2), sd = se_diff) # corresponds to two-sided comparison
    return(pv)
  }
  
  
  #rev. 4
  
  #compute error in reconstructed counts in the half-planes from directional measurements
  compute_s_r_total = function(N_total, B, Bi)
  {
  #N_total: vector of 2; total counts that would result within half-planes (including non-epithermal neutrons (N_non_epi))
  #B: matrix for converting counts in the half-planes (N_1, N_2) to counts registered by the directional scanner (N_f1, N_f2)
  #Bi: inverse of B
  #eps: fraction of noise in N_total (N_non_epi = N_total * eps)
    s_r_total = sqrt(Bi^2 %*% B  %*% N_total) 
    return(s_r_total)
  }  
  
  #compute confidence intervals of count rates and p-values for pairwise comparison of these
  compute_ci_p2 = function(delta_t, R_1_t, R_2_t, B, Bi, eps = 0) 
  {
    #R_1_t, R_2_t: effective neutron count rates in the half-planes  (including non-epithermal neutrons (N_non_epi))
    #B: matrix for converting counts in the half-planes (N_1, N_2) to counts registered by the directional scanner (N_f1, N_f2)
    #Bi: inverse of B
    #eps: fraction of noise in N_total (N_non_epi = N_total * eps)
    
    ci_lohi      = array() #, 6, dimnames = list(c("R_1", "R_2", "R_f1", "R_f2", "R_1r", "R_2r"))) #confidence intervals of measurements
    p_value_comp = array() #NA, 6, dimnames = list(c("N_1.N_2", "N_1.N_f1", "N_1.N_f2", "N_f1.N_f2", "N_f1.N_2", "N_2.N_f2"))) #p-values of two-sided comparison
    
    
    N_1_t = R_1_t * delta_t #theoretical number of counts expected
    N_2_t = R_2_t * delta_t

    tt = round(c(N_1_t, N_2_t) %*% B) #convert to mixed signals (i.e. counts expected for the directional operation)
    N_f1_t = tt[1]
    N_f2_t = tt[2]
    
    N_1_epi =  (1-eps) * N_1_t
    N_2_epi =  (1-eps) * N_2_t
    
    N_f1_epi = (1-eps) * N_f1_t
    N_f2_epi = (1-eps) * N_f2_t
    
    N_1_t = round(N_1_t)  # as counts are integer, we need to round here (we didn't round in the first place to avoid carrying the effect of rounding into N_f)
    N_2_t = round(N_2_t)
    
    #compute standard deviations (errors) of count rates, assuming Poisson characteristics
    s_1_t   = sqrt(N_1_t) 
    s_2_t   = sqrt(N_2_t)
    s_1_epi   = sqrt(N_1_epi) 
    s_2_epi   = sqrt(N_2_epi)
    s_f1_t  = sqrt(N_f1_t)
    s_f2_t  = sqrt(N_f2_t)
    s_f1_epi  = sqrt(N_f1_epi)
    s_f2_epi  = sqrt(N_f2_epi)
    # browser()
    # 
    # s_f1_epi/N_f1_epi
    # s_f2_epi/N_f2_epi

    
    #assume addition of Gaussian errors for computing total errors (we approximate Poisson quantiles with Gaussian, so we can add them this way)
  
    # errors of reconstructed values for N_1 and N_2
    tt = compute_s_r_total(c(N_1_t, N_2_t), B, Bi)
    s_1r_t = tt[1]
    s_2r_t = tt[2]
    
    #for the (virtual) epithermal rates, this error is just scaled proportionally (probably not used anyway) 
    s_1r_epi = sqrt((1-gamma_no) * s_1r_t^2)
    s_2r_epi = sqrt((1-gamma_no) * s_2r_t^2)
    
    #confidence intervals of measurement from errors computed above
    #approximate poisson quantiles with gaussian, so we can add them as above
    #qpois(p = c(p_thresh/2, 1-p_thresh/2), lambda =  N_1) #is the same for large N
    #as gaussian implies symmetry, we only compute lower side
    
    ci_lohi["R_1_t"]     = qnorm(p = p_thresh/2, mean =  N_1_t,  sd=s_1_t  ) 
    ci_lohi["R_2_t"]     = qnorm(p = p_thresh/2, mean =  N_2_t,  sd=s_2_t  )
    ci_lohi["R_1_epi"]     = qnorm(p = p_thresh/2, mean =  N_1_epi,  sd=s_1_epi  ) 
    ci_lohi["R_2_epi"]     = qnorm(p = p_thresh/2, mean =  N_2_epi,  sd=s_2_epi  )
    
    ci_lohi["R_1r_t"]    = qnorm(p = p_thresh/2, mean =  N_1_t,  sd=s_1r_t  ) 
    ci_lohi["R_2r_t"]    = qnorm(p = p_thresh/2, mean =  N_2_t,  sd=s_2r_t  )
    ci_lohi["R_1r_epi"]    = qnorm(p = p_thresh/2, mean =  N_1_epi,  sd=s_1r_epi  ) 
    ci_lohi["R_2r_epi"]    = qnorm(p = p_thresh/2, mean =  N_2_epi,  sd=s_2r_epi  )
    
    ci_lohi["R_f1_t"]    = qnorm(p = p_thresh/2, mean =  N_f1_t, sd=s_f1_t )
    ci_lohi["R_f2_t"]    = qnorm(p = p_thresh/2, mean =  N_f2_t, sd=s_f2_t )
    ci_lohi["R_f1_epi"]    = qnorm(p = p_thresh/2, mean =  N_f1_epi, sd=s_f1_epi )
    ci_lohi["R_f2_epi"]    = qnorm(p = p_thresh/2, mean =  N_f2_epi, sd=s_f2_epi )
    
    (N_f1_epi - ci_lohi["R_f1_epi"]) /  N_f1_epi
    (N_f2_epi - ci_lohi["R_f2_epi"]) /  N_f2_epi
    
    ci_lohi = ci_lohi /delta_t #convert counts to count rates
    
    #do statistical tests
    p_value_comp["N_1_t.N_2_t"]    = p_value_two_sided(mu1=N_1_t,  mu2=N_2_t,  se1=s_1_t,  se2=s_2_t) 
    p_value_comp["N_1_epi.N_2_epi"]    = p_value_two_sided(mu1=N_1_epi,  mu2=N_2_epi,  se1=s_1_epi,  se2=s_2_epi) 
    
    p_value_comp["N_1r_t.N_2r_t"]    = p_value_two_sided(mu1=N_1_t,  mu2=N_2_t,  se1=s_1r_t,  se2=s_2r_t) 
    p_value_comp["N_1r_epi.N_2r_epi"]    = p_value_two_sided(mu1=N_1_epi,  mu2=N_2_epi,  se1=s_1r_epi,  se2=s_2r_epi) 
    
    p_value_comp["N_f1_t.N_f2_t"]  = p_value_two_sided(mu1=N_f1_t,  mu2=N_f2_t,  se1=s_f1_t,  se2=s_f2_t)  
    p_value_comp["N_f1_epi.N_f2_epi"]  = p_value_two_sided(mu1=N_f1_epi,  mu2=N_f2_epi,  se1=s_f1_epi,  se2=s_f2_epi)  
    
    #discard "NA" value at first position; relic of initialisation
    ci_lohi = ci_lohi[-1]
    p_value_comp = p_value_comp[-1]
    
    return(c(ci_lohi, p_value_comp))
  }
  
  #auxiliary function used in the optimization to find threshold values
  difference_to_threshold_ci = function(delta_t, R_1, contrast, B, Bi, target_sensor, precision) 
  {
    R_2=(1+contrast)*R_1
    #if (delta_t<=0) return(9999) #reject out-of-range-values
    res = compute_ci_p2(delta_t = delta_t, R_1=R_1, R_2=R_2, B=B, Bi=Bi, eps = eps)  
    true_values = array(NA, dim = 6, dimnames=list(c("R_1", "R_2","R_f1", "R_f2","R_1r", "R_2r")))
    true_values[c("R_1",   "R_2")] = c(R_1, R_2)
    true_values[c("R_1r", "R_2r")] = c(R_1, R_2)
    true_values[c("R_f1", "R_f2")] = true_values[1:2] %*% B #compute true values for R_f1 and R_f2
    true_values = c(true_values, true_values * (1-eps) )
    names(true_values) = paste0(names(true_values), c(rep("_t",6), rep("_epi",6)))
    
    current_prec = 2 * (true_values[target_sensor] - res[target_sensor]) #assume symmetry of confidence interval
    current_prec_relative = abs(current_prec) / true_values[target_sensor]  #convert to relative precision
    obj_fun = abs(current_prec_relative - precision)  #return deviation of target precision
    return(obj_fun)
  }
  
  #auxiliary function used in the optimization to find threshold values
  difference_to_threshold_p = function(delta_t, R_1, contrast, B, Bi, target_comparison, p_thresh) 
  {
    R_2=(1+contrast)*R_1
    #if (delta_t<=0) return(9999) #reject out-of-range-values
    res = compute_ci_p2(delta_t = delta_t, R_1=R_1, R_2=R_2, B=B, Bi=Bi)  
    
    current_p =  res[target_comparison]
    if (current_p < 1e-5) 
      return(contrast*10 + delta_t + R_1/100) #help to avoid the flat area where no difference can be detected for numerical reasons)
                        #drive optimization towards lower aggregation times, contrasts and count rates
    obj_fun = abs(current_p - p_thresh)  #return deviation of target p-value
    return(obj_fun)
  }
  
  