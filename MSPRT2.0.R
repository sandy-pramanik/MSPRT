
library(ggplot2)
library(ggpubr)
library(doParallel)

#### Type II error function of fixed-design tests ####

## arguments
# theta         : parameter value where the Type II error probability is evaluated
# test.type     : test type; oneProp, oneZ, oneT, twoZ, twoT
# side          : direction of H1; 'right', 'left'
# theta0        : hypothesized value under H0
# N             : available sample size in fixed design for one-sample tests
# N1            : available sample size from Group 1 in fixed design for two-sample tests
# N2            : available sample size from Group 2 in fixed design for two-sample tests
# Type1         : Type I error probability
# sigma         : known sd in one-sample z test
# sigma1        : known sd of Group 1 in two-sample z test
# sigma2        : known sd of Group 2 in two-sample z test

## returns Type II error probability at theta

Type2.fixed_design = function(theta, test.type, side = "right", theta0, N, N1, N2, 
                              Type1 = 0.005, sigma = 1, sigma1 = 1, sigma2 = 1){
  
  if((test.type!="oneProp") & (test.type!="oneZ") & (test.type!="oneT") &
     (test.type!="twoZ") & (test.type!="twoT")){
    return(print("Unknown 'test type'. Has to be one of 'oneProp', 'oneZ', 'oneT', 'twoZ' or 'twoT'."))
  }
  
  if(test.type=="oneProp"){
    
    ####################### one-sample proportion test #######################
    
    if(missing(theta0)) theta0 = 0.5
    
    if(side=="right"){
      
      c.alpha = qbinom(p = Type1, size = N, prob = theta0, lower.tail = F)
      return(pbinom(q = c.alpha, size = N, prob = theta))
      
    }else if(side=="left"){
      
      c.alpha = qbinom(p = Type1, size = N, prob = theta0)
      return(pbinom(q = c.alpha-1, size = N, prob = theta, lower.tail = F))
    }
    
  }else if(test.type=="oneZ"){
    
    ####################### one-sample z test #######################
    
    if(missing(theta0)) theta0 = 0
    
    if(side=="right"){
      
      z.alpha = qnorm(p = Type1, lower.tail = F)
      return(pnorm(q = theta0 + (z.alpha*sigma)/sqrt(N), mean = theta, sd = sigma/sqrt(N)))
      
    }else if(side=="left"){
      
      z.alpha = qnorm(p = Type1, lower.tail = F)
      return(pnorm(q = theta0 - (z.alpha*sigma)/sqrt(N), mean = theta, sd = sigma/sqrt(N),
                   lower.tail = F))
    }
    
  }else if(test.type=="oneT"){
    
    ####################### one-sample t test #######################
    
    if(missing(theta0)) theta0 = 0
    
    if(side=="right"){
      
      t.alpha = qt(p = Type1, df = N-1, lower.tail = F)
      return(pt(q = t.alpha, df = N-1, ncp = sqrt(N)*(theta - theta0)))
      
    }else if(side=="left"){
      
      t.alpha = qt(p = Type1, df = N-1, lower.tail = F)
      return(pt(q = -t.alpha, df = N-1, ncp = sqrt(N)*(theta - theta0),
                lower.tail = F))
    }
    
  }else if(test.type=="twoZ"){
    
    ####################### two-sample z test #######################
    
    if(missing(theta0)) theta0 = 0
    
    if(side=="right"){
      
      z.alpha = qnorm(p = Type1, lower.tail = F)
      sigmaD = sqrt((sigma1^2)/N1 + (sigma2^2)/N2)
      return(pnorm(q = theta0 + z.alpha*sigmaD, mean = theta, sd = sigmaD))
      
    }else if(side=="left"){
      
      z.alpha = qnorm(p = Type1, lower.tail = F)
      sigmaD = sqrt((sigma1^2)/N1 + (sigma2^2)/N2)
      return(pnorm(q = theta0 - z.alpha*sigmaD, mean = theta, sd = sigmaD,
                   lower.tail = F))
    }
    
  }else if(test.type=="twoT"){
    
    ####################### two-sample t test #######################
    
    if(missing(theta0)) theta0 = 0
    
    if(side=="right"){
      
      t.alpha = qt(p = Type1, df = N1 + N2 - 2, lower.tail = F)
      return(pt(q = t.alpha, df = N1 + N2 - 2, 
                ncp = (theta - theta0)/sqrt(1/N1 + 1/N2)))
      
    }else if(side=="left"){
      
      t.alpha = qt(p = Type1, df = N1 + N2 - 2, lower.tail = F)
      return(pt(q = -t.alpha, df = N1 + N2 - 2, 
                ncp = (theta - theta0)/sqrt(1/N1 + 1/N2), lower.tail = F))
    }
  }
}


#### Fixed-design alternative ####

## arguments
# test.type     : test type; oneProp, oneZ, oneT, twoZ, twoT
# side          : direction of H1; 'right', 'left'
# theta0        : hypothesized value under H0
# N             : available sample size in fixed design for one-sample tests
# N1            : available sample size from Group 1 in fixed design for two-sample tests
# N2            : available sample size from Group 2 in fixed design for two-sample tests
# Type1         : Type I error probability
# Type2         : Type II error probability
# sigma         : known sd in one-sample z test
# sigma1        : known sd of Group 1 in two-sample z test
# sigma2        : known sd of Group 2 in two-sample z test

## returns the fixed-design alternative

fixed_design.alt = function(test.type, side = "right", theta0, N, N1, N2, 
                            Type1 = 0.005, Type2 = .2, sigma = 1, sigma1 = 1, sigma2 = 1){
  
  if((test.type!="oneProp") & (test.type!="oneZ") & (test.type!="oneT") &
     (test.type!="twoZ") & (test.type!="twoT")){
    return(print("Unknown 'test type'. Has to be one of 'oneProp', 'oneZ', 'oneT', 'twoZ' or 'twoT'."))
  }
  
  if(test.type=="oneProp"){
    
    ####################### one-sample proportion test #######################
    
    if(missing(theta0)) theta0 = 0.5
    
    if(side=="right"){
      
      c.alpha = qbinom(p = Type1, size = N, prob = theta0, lower.tail = F)
      solve.out = uniroot(interval = c(theta0, 1),
                          f = function(x){
                            
                            pbinom(q = c.alpha, size = N, prob = x) - Type2
                          })
      
      return(solve.out$root)
      
    }else if(side=="left"){
      
      c.alpha = qbinom(p = Type1, size = N, prob = theta0)
      solve.out = uniroot(interval = c(0, theta0),
                          f = function(x){
                            
                            pbinom(q = c.alpha-1, size = N, prob = x, 
                                   lower.tail = F) - Type2
                          })
      
      return(solve.out$root)
    }
    
  }else if(test.type=="oneZ"){
    
    ####################### one-sample z test #######################
    
    if(missing(theta0)==T) theta0 = 0
    
    z.alpha = qnorm(p = Type1, lower.tail = F)
    if(side=="right"){
      
      return(theta0 - ((qnorm(p = Type2) - z.alpha)*sigma)/sqrt(N))
      
    }else if(side=="left"){
      
      return(theta0 - ((qnorm(p = 1-Type2) + z.alpha)*sigma)/sqrt(N))
    }
    
  }else if(test.type=="oneT"){
    
    ####################### one-sample t test #######################
    
    if(missing(theta0)==T) theta0 = 0
    
    t.alpha = qt(p = Type1, df = N-1, lower.tail = F)
    if(side=="right"){
      
      solve.out = uniroot(interval = c(theta0, .Machine$integer.max),
                          f = function(x){
                            
                            pt(q = t.alpha, df = N-1, ncp = sqrt(N)*(x - theta0)) -
                              Type2
                          })
      
      return(solve.out$root)
      
    }else if(side=="left"){
      
      solve.out = uniroot(interval = c(-.Machine$integer.max, theta0),
                          f = function(x){
                            
                            pt(q = -t.alpha, df = N-1, ncp = sqrt(N)*(x - theta0),
                               lower.tail = F) - Type2
                          })
      
      return(solve.out$root)
    }
    
  }else if(test.type=="twoZ"){
    
    ####################### two-sample z test #######################
    
    if(missing(theta0)==T) theta0 = 0
    
    z.alpha = qnorm(p = Type1, lower.tail = F)
    sigmaD = sqrt((sigma1^2)/N1 + (sigma2^2)/N2)
    if(side=="right"){
      
      return(theta0 - (qnorm(p = Type2) - z.alpha)*sigmaD)
      
    }else if(side=="left"){
      
      return(theta0 - (qnorm(p = 1-Type2) + z.alpha)*sigmaD)
    }
    
  }else if(test.type=="twoT"){
    
    ####################### two-sample t test #######################
    
    if(missing(theta0)==T) theta0 = 0
    
    t.alpha = qt(p = Type1, df = N1 + N2 - 2, lower.tail = F)
    if(side=="right"){
      
      solve.out = uniroot(interval = c(theta0, .Machine$integer.max),
                          f = function(x){
                            
                            pt(q = t.alpha, df = N1 + N2 - 2, 
                               ncp = (x - theta0)/sqrt(1/N1 + 1/N2)) - Type2
                          })
      
      return(solve.out$root)
      
    }else if(side=="left"){
      
      solve.out = uniroot(interval = c(-.Machine$integer.max, theta0),
                          f = function(x){
                            
                            pt(q = -t.alpha, df = N1 + N2 - 2, 
                               ncp = (x - theta0)/sqrt(1/N1 + 1/N2), 
                               lower.tail = F) - Type2
                          })
      
      return(solve.out$root)
    }
  }
}


#### UMPBT alternative  ####

## arguments
# test.type     : test type; oneProp, oneZ, oneT, twoZ, twoT
# side          : direction of H1; 'right', 'left'
# theta0        : hypothesized value under H0
# N             : available sample size in fixed design for one-sample tests
# N1            : available sample size from Group 1 in fixed design for two-sample tests
# N2            : available sample size from Group 2 in fixed design for two-sample tests
# Type1         : Type I error probability
# sigma         : known sd in one-sample z test
# sigma1        : known sd of Group 1 in two-sample z test
# sigma2        : known sd of Group 2 in two-sample z test
# obs           : observations in one-sample t-test
# sd.obs        : sd (divisor N-1) of observations in one-sample t-test
# obs1          : observations from Group 1 in two-sample t-test
# obs2          : observations from Group 2 in two-sample t-test
# pooled.sd     : pooled sd (divisor N1 + N2 - 2) observations in two-sample t-test

## returns the UMPBT alternative

UMPBT.alt = function(test.type, side = "right", theta0, N, N1, N2, 
                     Type1 = 0.005, sigma = 1, sigma1 = 1, sigma2 = 1,
                     obs, sd.obs, obs1, obs2, pooled.sd){
  
  if((test.type!="oneProp") & (test.type!="oneZ") & (test.type!="oneT") &
     (test.type!="twoZ") & (test.type!="twoT")){
    return(print("Unknown 'test type'. Has to be one of 'oneProp', 'oneZ', 'oneT', 'twoZ' or 'twoT'."))
  }
  
  if(test.type=="oneProp"){
    
    ####################### one-sample proportion test #######################
    
    if(missing(theta0)) theta0 = 0.5
    
    if(side=="right"){
      
      # fixed design cutoff; rejection region (c.alpha, N]
      c.alpha = qbinom(p = Type1, size = N, prob = theta0, lower.tail = F)
      
      #### finding the outer UMPBT alternative ####
      # solving for BF threshold in UMPBT
      solve.delta.outer = 
        nleqslv::nleqslv(x = 3,
                         fn = function(delta){
                           
                           out.optimize.UMPBTobjective = 
                             optimize(interval = c(theta0, 1),
                                      f = function(p){
                                        
                                        (log(delta) - N*(log(1 - p) - log(1 - theta0)))/
                                          (log(p/(1 - p)) - log(theta0/(1 - theta0)))
                                      })
                           
                           out.optimize.UMPBTobjective$objective - c.alpha
                         })
      
      # the outer UMPBT alternative
      out.optimize.UMPBTobjective.outer = 
        optimize(interval = c(theta0, 1),
                 f = function(p){
                   
                   (log(solve.delta.outer$x) - N*(log(1 - p) - log(1 - theta0)))/
                     (log(p/(1 - p)) - log(theta0/(1 - theta0)))
                 })
      
      #### finding the inner UMPBT alternative ####
      # solving for BF threshold in UMPBT
      solve.delta.inner = 
        nleqslv::nleqslv(x = 3,
                         fn = function(delta){
                           
                           out.optimize.UMPBTobjective = 
                             optimize(interval = c(theta0, 1),
                                      f = function(p){
                                        
                                        (log(delta) - N*(log(1 - p) - log(1 - theta0)))/
                                          (log(p/(1 - p)) - log(theta0/(1 - theta0)))
                                      })
                           
                           out.optimize.UMPBTobjective$objective - (c.alpha - 1)
                         })
      
      # the inner UMPBT alternative
      out.optimize.UMPBTobjective.inner = 
        optimize(interval = c(theta0, 1),
                 f = function(p){
                   
                   (log(solve.delta.inner$x) - N*(log(1 - p) - log(1 - theta0)))/
                     (log(p/(1 - p)) - log(theta0/(1 - theta0)))
                 })
      
      #### probability corresponding to the outer UMPBT alternative ####
      mix.prob = (Type1 - pbinom(q = c.alpha, size = N, prob = theta0, lower.tail = F))/
        dbinom(x = c.alpha, size = N, prob = theta0)
      
      return(list("theta" = c(out.optimize.UMPBTobjective.inner$minimum,
                              out.optimize.UMPBTobjective.outer$minimum),
                  "mix.prob" = c(mix.prob, 1-mix.prob)))
      
    }else if(side=="left"){
      
      # fixed design cutoff; rejection region [0, c.alpha)
      c.alpha = qbinom(p = Type1, size = N, prob = theta0)
      
      #### finding the outer UMPBT alternative ####
      # solving for BF threshold in UMPBT
      solve.delta.outer = 
        nleqslv::nleqslv(x = 3,
                         fn = function(delta){
                           
                           out.optimize.UMPBTobjective = 
                             optimize(interval = c(0, theta0), maximum = T,
                                      f = function(p){
                                        
                                        (log(delta) - N*(log(1 - p) - log(1 - theta0)))/
                                          (log(p/(1 - p)) - log(theta0/(1 - theta0)))
                                      })
                           
                           out.optimize.UMPBTobjective$objective - c.alpha
                         })
      
      # the outer UMPBT alternative
      out.optimize.UMPBTobjective.outer = 
        optimize(interval = c(0, theta0), maximum = T,
                 f = function(p){
                   
                   (log(solve.delta.outer$x) - N*(log(1 - p) - log(1 - theta0)))/
                     (log(p/(1 - p)) - log(theta0/(1 - theta0)))
                 })
      
      #### finding the inner UMPBT alternative ####
      # solving for BF threshold in UMPBT
      solve.delta.inner = 
        nleqslv::nleqslv(x = 3,
                         fn = function(delta){
                           
                           out.optimize.UMPBTobjective = 
                             optimize(interval = c(0, theta0), maximum = T,
                                      f = function(p){
                                        
                                        (log(delta) - N*(log(1 - p) - log(1 - theta0)))/
                                          (log(p/(1 - p)) - log(theta0/(1 - theta0)))
                                      })
                           
                           out.optimize.UMPBTobjective$objective - (c.alpha + 1)
                         })
      
      # the inner UMPBT alternative
      out.optimize.UMPBTobjective.inner = 
        optimize(interval = c(0, theta0), maximum = T,
                 f = function(p){
                   
                   (log(solve.delta.inner$x) - N*(log(1 - p) - log(1 - theta0)))/
                     (log(p/(1 - p)) - log(theta0/(1 - theta0)))
                 })
      
      #### probability corresponding to the outer UMPBT alternative ####
      mix.prob = (Type1 - pbinom(q = c.alpha-1, size = N, prob = theta0))/
        dbinom(x = c.alpha, size = N, prob = theta0)
      
      return(list("theta" = c(out.optimize.UMPBTobjective.inner$maximum,
                              out.optimize.UMPBTobjective.outer$maximum),
                  "mix.prob" = c(mix.prob, 1-mix.prob)))
    }
    
  }else if(test.type=="oneZ"){
    
    ####################### one-sample z test #######################
    
    if(missing(theta0)==T) theta0 = 0
    
    z.alpha = qnorm(p = Type1, lower.tail = F)
    if(side=="right"){
      
      return(theta0 + (z.alpha*sigma)/sqrt(N))
      
    }else if(side=="left"){
      
      return(theta0 - (z.alpha*sigma)/sqrt(N))
    }
    
  }else if(test.type=="oneT"){
    
    ####################### one-sample t test #######################
    
    if(missing(theta0)) theta0 = 0
    
    if(missing(sd.obs)){
      
      if(missing(obs)){
        
        return("Need to provide either 'sd.obs' or 'obs'.")
        
      }else{
        
        sd.obs = sd(obs)
      }
      
    }else{
      
      if((!missing(obs))&&(round(sd(obs), 5)!=sd.obs)){
        
        return("The 'sd.obs' provided by the user doesn't match with the sd (divisor (n-1)) calculated from 'obs'.")
        
      }else{
        
        print("Ignoring 'obs'. Specifying 'sd.obs' is enough.")
      }
    }
    
    t.alpha = qt(p = Type1, df = N-1, lower.tail = F)
    if(side=="right"){
      
      return(theta0 + (t.alpha*sd.obs)/sqrt(N))
      
    }else if(side=="left"){
      
      return(theta0 - (t.alpha*sd.obs)/sqrt(N))
    }
    
  }else if(test.type=="twoZ"){
    
    ####################### two-sample z test #######################
    
    if(missing(theta0)) theta0 = 0
    
    z.alpha = qnorm(p = Type1, lower.tail = F)
    if(side=="right"){
      
      return(theta0 + z.alpha*sqrt((sigma1^2)/N1 + (sigma2^2)/N2))
      
    }else if(side=="left"){
      
      return(theta0 - z.alpha*sqrt((sigma1^2)/N1 + (sigma2^2)/N2))
    }
    
  }else if(test.type=="twoT"){
    
    ####################### two-sample t test #######################
    
    if(missing(theta0)) theta0 = 0
    
    if(missing(pooled.sd)){
      
      if(any(c(!missing(obs1),!missing(obs1)))){
        
        return("Need to provide either 'pooled.sd' or both 'obs1' and 'obs2.")
        
      }else{
        
        pooled.sd = sqrt(((length(obs1)-1)*var(obs1) + (length(obs2)-1)*var(obs2))/
                           (length(obs1) + length(obs2) - 2))
      }
      
    }else{
      
      if((!missing(obs1))&&(!missing(obs2))&&
         (round(sqrt(((length(obs1)-1)*var(obs1) + (length(obs2)-1)*var(obs2))/
                     (length(obs1) + length(obs2) - 2)), 5)!=pooled.sd)){
        
        return("The 'pooled.sd' provided by the user doesn't match with the pooled sd (divisor (N1 + N2 - 2)) calculated from 'obs1' and 'obs2.")
        
      }else{
        
        print("Ignoring 'obs1' and/or 'obs2. Specifying 'pooled.sd' is enough.")
      }
    }
    
    t.alpha = qt(p = Type1, df = N1 + N2 - 2, lower.tail = F)
    if(side=="right"){
      
      return(theta0 + t.alpha*pooled.sd*sqrt(1/N1 + 1/N2))
      
    }else if(side=="left"){
      
      return(theta0 - t.alpha*pooled.sd*sqrt(1/N1 + 1/N2))
    }
  }
}


################################# designing the MSPRT #################################

#### one-sample proportion test ####
design.MSPRT_oneProp = function(side = 'right', theta0 = 0.5, theta1 = T,
                                Type1.target =.005, Type2.target = .2,
                                N.max, batch.size,
                                nReplicate = 1e+6, verbose = T, seed = 1){
  
  if(side!='both'){
    
    ########################### one-sample proportion (right/left sided) ###########################
    
    ## batch sizes and N.max
    if(missing(batch.size)){
      
      if(missing(N.max)){
        
        return("Either 'batch.size' or 'N.max' needs to be specified")
        
      }else{batch.size = rep(1, N.max)}
      
    }else{
      
      if(missing(N.max)){
        
        N.max = sum(batch.size)
        
      }else{
        
        if(sum(batch.size)!=N.max) return("Sum of batch sizes should add up to N.max")
      }
    }
    
    nAnalyses = length(batch.size)
    
    ## msg
    if(verbose){
      
      if(any(batch.size>1)){
        
        cat('\n')
        print("=========================================================================")
        print("Designing the group sequential MSPRT for a one-sample proportion test:")
        print("=========================================================================")
        
      }else{
        
        cat('\n')
        print("=========================================================================")
        print("Designing the sequential MSPRT for a one-sample proportion test:")
        print("=========================================================================")
      }
      
      print(paste("Maximum available sample size: ", N.max, sep = ""))
      print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
      print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
      print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
      print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
      print(paste("Hypothesized value under H0: ", theta0, sep = ""))
      print(paste("Direction of the H1: ", side, sep = ""))
    }
    
    batch.size = c(0, cumsum(batch.size))
    
    
    if(is.logical(theta1)&&(theta1==F)){
      
      ################ no alternative comparison ################
      
      ################ UMPBT alternative ################
      UMPBT = UMPBT.alt(test.type = 'oneProp', side = side, theta0 = theta0,
                        N = N.max, Type1 = Type1.target)
      
      # msg
      if(verbose==T){
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(UMPBT$theta[1], 3), " & ",
                    round(UMPBT$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$mix.prob[1], 3), " & ", 1 - round(UMPBT$mix.prob[1], 3), sep = ''))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # required storages
      cumsum0_n = LR0_n = numeric(nReplicate)
      type1.error.AR = rep(F, nReplicate)
      N0.AR = rep(N.max, nReplicate)
      not.reached.decisionH0.AR = 1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # sum of observations at step n
          sum0_n = rbinom(length(not.reached.decisionH0.AR),
                          batch.size[n+1]-batch.size[n], theta0)
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + sum0_n
          
          # likelihood ratio of observations until step n
          LR0_n[not.reached.decisionH0.AR] = 
            UMPBT$mix.prob[1]*(((1 - UMPBT$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$theta[1])))^cumsum0_n[not.reached.decisionH0.AR] +
            (1 - UMPBT$mix.prob[2])*(((1 - UMPBT$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$theta[2])))^cumsum0_n[not.reached.decisionH0.AR]
          
          # comparing with the thresholds
          AcceptedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]<=Accept.threshold)
          RejectedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]>=Reject.threshold)
          reached.decisionH0_n.AR = union(AcceptedH0.underH0_n.AR, RejectedH0.underH0_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH0_n.AR)>0){
            
            N0.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch.size[n+1]
            type1.error.AR[not.reached.decisionH0.AR[RejectedH0.underH0_n.AR]] = T
            not.reached.decisionH0.AR = not.reached.decisionH0.AR[-reached.decisionH0_n.AR]
          }
        }
        
        setTxtProgressBar(pb, n)
      }
      
      # determining termination threshold
      # H0 is rejected if LR or (BF) is >= termination threshold
      nNot.reached.decisionH0.AR = length(not.reached.decisionH0.AR)
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 Reject.threshold -
                                                   max(LR0_n[not.reached.decisionH0.AR]))))
          termination.threshold.AR = (floor(max(LR0_n[not.reached.decisionH0.AR])*
                                              (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 min(LR0_n[not.reached.decisionH0.AR]) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(LR0_n[not.reached.decisionH0.AR]))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(LR0_n[not.reached.decisionH0.AR]))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(min(cumRejFreq_not.reached.decisionH0.AR)>nNewRejects.AR){
            
            nDecimal.accuracy = 
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR + 
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      
      # Expected sample sizes
      EN0 = mean(N0.AR)
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", termination.threshold.AR))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  "N0" = list('H0' = N0.AR), "EN" = EN0,
                  "UMPBT" = UMPBT, "Type2.fixed.design" = Type2.target, 
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneProp', 'side' = side, 'theta0' = theta0, 
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
      
    }else if(is.logical(theta1)&&(theta1==T)){
      
      ################ comparison at the fixed-design alternative (default) ################
      theta1 = fixed_design.alt(test.type = 'oneProp', side = side, theta0 = theta0,
                                N = N.max, Type1 = Type1.target, Type2 = Type2.target)
      
      ################ UMPBT alternative ################
      UMPBT = UMPBT.alt(test.type = 'oneProp', side = side, theta0 = theta0,
                        N = N.max, Type1 = Type1.target)
      
      # msg
      if(verbose==T){
        
        print(paste("Alternative under comparison: ", round(theta1, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(UMPBT$theta[1], 3), " & ",
                    round(UMPBT$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$mix.prob[1], 3), " & ", 1 - round(UMPBT$mix.prob[1], 3), sep = ''))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # required storages
      cumsum0_n = cumsum1_n = LR0_n = LR1_n = numeric(nReplicate)
      type1.error.AR = type2.error.AR = rep(F, nReplicate)
      N0.AR = N1.AR = rep(N.max, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH1.AR = 1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # sum of observations at step n
          sum0_n = rbinom(length(not.reached.decisionH0.AR),
                          batch.size[n+1]-batch.size[n], theta0)
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + sum0_n
          
          # likelihood ratio of observations until step n
          LR0_n[not.reached.decisionH0.AR] = 
            UMPBT$mix.prob[1]*(((1 - UMPBT$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$theta[1])))^cumsum0_n[not.reached.decisionH0.AR] +
            (1 - UMPBT$mix.prob[2])*(((1 - UMPBT$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$theta[2])))^cumsum0_n[not.reached.decisionH0.AR]
          
          # comparing with the thresholds
          AcceptedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]<=Accept.threshold)
          RejectedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]>=Reject.threshold)
          reached.decisionH0_n.AR = union(AcceptedH0.underH0_n.AR, RejectedH0.underH0_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH0_n.AR)>0){
            
            N0.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch.size[n+1]
            type1.error.AR[not.reached.decisionH0.AR[RejectedH0.underH0_n.AR]] = T
            not.reached.decisionH0.AR = not.reached.decisionH0.AR[-reached.decisionH0_n.AR]
          }
        }
        
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          # sum of observations at step n
          sum1_n = rbinom(length(not.reached.decisionH1.AR),
                          batch.size[n+1]-batch.size[n], theta1)
          
          # sum of observations until step n
          cumsum1_n[not.reached.decisionH1.AR] = 
            cumsum1_n[not.reached.decisionH1.AR] + sum1_n
          
          # likelihood ratio of observations until step n
          LR1_n[not.reached.decisionH1.AR] = 
            UMPBT$mix.prob[1]*(((1 - UMPBT$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$theta[1])))^cumsum1_n[not.reached.decisionH1.AR] +
            (1 - UMPBT$mix.prob[2])*(((1 - UMPBT$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$theta[2])))^cumsum1_n[not.reached.decisionH1.AR]
          
          # comparing with the thresholds
          AcceptedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]<=Accept.threshold)
          RejectedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]>=Reject.threshold)
          reached.decisionH1_n.AR = union(AcceptedH0.underH1_n.AR, RejectedH0.underH1_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH1_n.AR)>0){
            
            N1.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch.size[n+1]
            type2.error.AR[not.reached.decisionH1.AR[AcceptedH0.underH1_n.AR]] = T
            not.reached.decisionH1.AR = not.reached.decisionH1.AR[-reached.decisionH1_n.AR]
          }
        }
        
        setTxtProgressBar(pb, n)
      }
      
      # determining termination threshold
      # H0 is rejected if LR or (BF) is >= termination threshold
      nNot.reached.decisionH0.AR = length(not.reached.decisionH0.AR)
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 Reject.threshold -
                                                   max(LR0_n[not.reached.decisionH0.AR]))))
          termination.threshold.AR = (floor(max(LR0_n[not.reached.decisionH0.AR])*
                                              (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 min(LR0_n[not.reached.decisionH0.AR]) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(LR0_n[not.reached.decisionH0.AR]))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(LR0_n[not.reached.decisionH0.AR]))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(min(cumRejFreq_not.reached.decisionH0.AR)>nNewRejects.AR){
            
            nDecimal.accuracy = 
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR + 
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      # attained Type II error probability
      actual.type2.error.AR = mean(type2.error.AR) +
        sum(LR1_n[not.reached.decisionH1.AR]<termination.threshold.AR)/nReplicate
      
      
      
      # Expected sample sizes
      EN0 = mean(N0.AR)
      EN1 = mean(N1.AR)
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", termination.threshold.AR))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Attained Type II error probability: ", round(actual.type2.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print(paste("Expected sample size at the alternative: ", round(EN1, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR, 
                  "Type2.attained" = actual.type2.error.AR, 
                  "N0" = list('H0' = N0.AR, 'H1' = N1.AR), "EN" = c(EN0, EN1),
                  "UMPBT" = UMPBT,
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target, 
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneProp', 'side' = side, 'theta0' = theta0, 
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
      
    }else{
      
      ################ comparison at user provided point alternative ################
      
      ################ UMPBT alternative ################
      UMPBT = UMPBT.alt(test.type = 'oneProp', side = side, theta0 = theta0,
                        N = N.max, Type1 = Type1.target)
      
      # msg
      if(verbose==T){
        
        print(paste("Alternative under comparison: ", round(theta1, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(UMPBT$theta[1], 3), " & ",
                    round(UMPBT$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$mix.prob[1], 3), " & ", 1 - round(UMPBT$mix.prob[1], 3), sep = ''))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # required storages
      cumsum0_n = cumsum1_n = LR0_n = LR1_n = numeric(nReplicate)
      type1.error.AR = type2.error.AR = rep(F, nReplicate)
      N0.AR = N1.AR = rep(N.max, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH1.AR = 1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # sum of observations at step n
          sum0_n = rbinom(length(not.reached.decisionH0.AR),
                          batch.size[n+1]-batch.size[n], theta0)
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + sum0_n
          
          # likelihood ratio of observations until step n
          LR0_n[not.reached.decisionH0.AR] = 
            UMPBT$mix.prob[1]*(((1 - UMPBT$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$theta[1])))^cumsum0_n[not.reached.decisionH0.AR] +
            (1 - UMPBT$mix.prob[2])*(((1 - UMPBT$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$theta[2])))^cumsum0_n[not.reached.decisionH0.AR]
          
          # comparing with the thresholds
          AcceptedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]<=Accept.threshold)
          RejectedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]>=Reject.threshold)
          reached.decisionH0_n.AR = union(AcceptedH0.underH0_n.AR, RejectedH0.underH0_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH0_n.AR)>0){
            
            N0.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch.size[n+1]
            type1.error.AR[not.reached.decisionH0.AR[RejectedH0.underH0_n.AR]] = T
            not.reached.decisionH0.AR = not.reached.decisionH0.AR[-reached.decisionH0_n.AR]
          }
        }
        
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          # sum of observations at step n
          sum1_n = rbinom(length(not.reached.decisionH1.AR),
                          batch.size[n+1]-batch.size[n], theta1)
          
          # sum of observations until step n
          cumsum1_n[not.reached.decisionH1.AR] = 
            cumsum1_n[not.reached.decisionH1.AR] + sum1_n
          
          # likelihood ratio of observations until step n
          LR1_n[not.reached.decisionH1.AR] = 
            UMPBT$mix.prob[1]*(((1 - UMPBT$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$theta[1])))^cumsum1_n[not.reached.decisionH1.AR] +
            (1 - UMPBT$mix.prob[2])*(((1 - UMPBT$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$theta[2])))^cumsum1_n[not.reached.decisionH1.AR]
          
          # comparing with the thresholds
          AcceptedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]<=Accept.threshold)
          RejectedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]>=Reject.threshold)
          reached.decisionH1_n.AR = union(AcceptedH0.underH1_n.AR, RejectedH0.underH1_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH1_n.AR)>0){
            
            N1.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch.size[n+1]
            type2.error.AR[not.reached.decisionH1.AR[AcceptedH0.underH1_n.AR]] = T
            not.reached.decisionH1.AR = not.reached.decisionH1.AR[-reached.decisionH1_n.AR]
          }
        }
        
        setTxtProgressBar(pb, n)
      }
      
      # determining termination threshold
      # H0 is rejected if LR or (BF) is >= termination threshold
      nNot.reached.decisionH0.AR = length(not.reached.decisionH0.AR)
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 Reject.threshold -
                                                   max(LR0_n[not.reached.decisionH0.AR]))))
          termination.threshold.AR = (floor(max(LR0_n[not.reached.decisionH0.AR])*
                                              (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 min(LR0_n[not.reached.decisionH0.AR]) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(LR0_n[not.reached.decisionH0.AR]))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(LR0_n[not.reached.decisionH0.AR]))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(min(cumRejFreq_not.reached.decisionH0.AR)>nNewRejects.AR){
            
            nDecimal.accuracy = 
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR + 
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      # attained Type II error probability
      actual.type2.error.AR = mean(type2.error.AR) +
        sum(LR1_n[not.reached.decisionH1.AR]<termination.threshold.AR)/nReplicate
      
      
      
      # Expected sample sizes
      EN0 = mean(N0.AR)
      EN1 = mean(N1.AR)
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", termination.threshold.AR))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Attained Type II error probability: ", round(actual.type2.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print(paste("Expected sample size at the alternative: ", round(EN1, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR, 
                  "Type2.attained" = actual.type2.error.AR, 
                  "N0" = list('H0' = N0.AR, 'H1' = N1.AR), "EN" = c(EN0, EN1),
                  "UMPBT" = UMPBT,
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target, 
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneProp', 'side' = side, 'theta0' = theta0, 
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
    }
    
  }else if(side=='both'){
    
    ########################### one-sample proportion (both sided) ###########################
    
    ## batch sizes and N.max
    if(missing(batch.size)){
      
      if(missing(N.max)){
        
        return("Either 'batch.size' or 'N.max' needs to be specified")
        
      }else{batch.size = rep(1, N.max)}
      
    }else{
      
      if(missing(N.max)){
        
        N.max = sum(batch.size)
        
      }else{
        
        if(sum(batch.size)!=N.max) return("Sum of batch sizes should add up to N.max")
      }
    }
    
    nAnalyses = length(batch.size)
    
    ## msg
    if(verbose){
      
      if(any(batch.size>1)){
        
        cat('\n')
        print("=========================================================================")
        print("Designing the group sequential MSPRT for a one-sample proportion test:")
        print("=========================================================================")
        
      }else{
        
        cat('\n')
        print("=========================================================================")
        print("Designing the sequential MSPRT for a one-sample proportion test:")
        print("=========================================================================")
      }
      
      print(paste("Maximum available sample size: ", N.max, sep = ""))
      print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
      print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
      print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
      print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
      print(paste("Hypothesized value under H0: ", theta0, sep = ""))
      print(paste("Direction of the H1: ", side, sep = ""))
    }
    
    batch.size = c(0, cumsum(batch.size))
    
    
    if(is.logical(theta1)&&(theta1==F)){
      
      ################ no fixed-design alternative ################
      
      ################ UMPBT alternative ################
      UMPBT = list('right' = UMPBT.alt(test.type = 'oneProp', side = 'right', 
                                       theta0 = theta0, N = N.max, Type1 = Type1.target/2),
                   'left' = UMPBT.alt(test.type = 'oneProp', side = 'left',
                                      theta0 = theta0, N = N.max, Type1 = Type1.target/2))
      
      # msg
      if(verbose==T){
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(UMPBT$right$theta[1], 3), " & ",
                    round(UMPBT$right$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$right$mix.prob[1], 3), " & ", 1 - round(UMPBT$right$mix.prob[1], 3),
                    sep = ""))
        print(paste(' On the left: ', round(UMPBT$left$theta[1], 3), " & ",
                    round(UMPBT$left$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$left$mix.prob[1], 3), " & ", 1 - round(UMPBT$left$mix.prob[1], 3),
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # required storages
      cumsum0_n = LR0_n.r = LR0_n.l = numeric(nReplicate)
      type1.error.AR = rep(F, nReplicate)
      N0.AR = N0.AR.r = N0.AR.l = rep(N.max, nReplicate)
      decision.underH0.AR.r = decision.underH0.AR.l = rep(NA, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.l =
        1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # sum of observations at step n
          sum0_n = rbinom(length(not.reached.decisionH0.AR),
                          batch.size[n+1]-batch.size[n], theta0)
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + sum0_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR0_n.r[not.reached.decisionH0.AR.r] = 
            UMPBT$right$mix.prob[1]*(((1 - UMPBT$right$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[1])))^cumsum0_n[not.reached.decisionH0.AR.r] +
            (1 - UMPBT$right$mix.prob[2])*(((1 - UMPBT$right$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[2])))^cumsum0_n[not.reached.decisionH0.AR.r]
          
          # for left sided check
          LR0_n.l[not.reached.decisionH0.AR.l] = 
            UMPBT$left$mix.prob[1]*(((1 - UMPBT$left$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[1])))^cumsum0_n[not.reached.decisionH0.AR.l] +
            (1 - UMPBT$left$mix.prob[2])*(((1 - UMPBT$left$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[2])))^cumsum0_n[not.reached.decisionH0.AR.l]
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]<=Accept.threshold
          RejectedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]>=Reject.threshold
          reached.decisionH0_n.AR.r = AcceptedH0.underH0_n.AR.r|RejectedH0.underH0_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.r)){
            
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[AcceptedH0.underH0_n.AR.r]] = 'A'
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[RejectedH0.underH0_n.AR.r]] = 'R'
            N0.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.r[!reached.decisionH0_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]<=Accept.threshold
          RejectedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]>=Reject.threshold
          reached.decisionH0_n.AR.l = AcceptedH0.underH0_n.AR.l|RejectedH0.underH0_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.l)){
            
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[AcceptedH0.underH0_n.AR.l]] = 'A'
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[RejectedH0.underH0_n.AR.l]] = 'R'
            N0.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH0.AR.l = not.reached.decisionH0.AR.l[!reached.decisionH0_n.AR.l]
          }
          
          not.reached.decisionH0.AR = union(not.reached.decisionH0.AR.r,
                                            not.reached.decisionH0.AR.l)
        }
        
        setTxtProgressBar(pb, n)
      }
      
      
      ### both-sided checking
      ## under H0
      # accepted or rejected ones
      accepted.by.both0 = intersect(which(decision.underH0.AR.r=='A'),
                                    which(decision.underH0.AR.l=='A'))
      onlyrejected.by.right0 = intersect(which(decision.underH0.AR.r=='R'),
                                         which(decision.underH0.AR.l!='R'))
      onlyrejected.by.left0 = intersect(which(decision.underH0.AR.r!='R'),
                                        which(decision.underH0.AR.l=='R'))
      rejected.by.both0 = intersect(which(decision.underH0.AR.r=='R'),
                                    which(decision.underH0.AR.l=='R'))
      
      # sample sizes required
      N0.AR[accepted.by.both0] = pmax(N0.AR.r[accepted.by.both0],
                                      N0.AR.l[accepted.by.both0])
      N0.AR[onlyrejected.by.right0] = N0.AR.r[onlyrejected.by.right0]
      N0.AR[onlyrejected.by.left0] = N0.AR.l[onlyrejected.by.left0]
      N0.AR[rejected.by.both0] = pmin(N0.AR.r[rejected.by.both0],
                                      N0.AR.l[rejected.by.both0])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right0 = intersect(which(decision.underH0.AR.r=='A'),
                                         which(is.na(decision.underH0.AR.l)))
      onlyaccepted.by.left0 = intersect(which(is.na(decision.underH0.AR.r)),
                                        which(decision.underH0.AR.l=='A'))
      both.inconclusive0 = intersect(which(is.na(decision.underH0.AR.r)),
                                     which(is.na(decision.underH0.AR.l)))
      all.inconclusive0 = c(onlyaccepted.by.right0, onlyaccepted.by.left0,
                            both.inconclusive0)
      nNot.reached.decisionH0.AR = length(all.inconclusive0)
      
      # Type I error probability
      type1.error.AR[c(onlyrejected.by.right0, onlyrejected.by.left0,
                       rejected.by.both0)] = T
      
      
      ## determining termination threshold
      ## H0 is rejected if LR or (BF) is >= termination threshold
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        term.thresh.possible.choices =
          c(LR0_n.r[onlyaccepted.by.left0],
            LR0_n.l[onlyaccepted.by.right0],
            pmin(LR0_n.r[both.inconclusive0], LR0_n.l[both.inconclusive0]))
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          max.LR0_n = max(term.thresh.possible.choices)
          nDecimal.accuracy = ceiling(-log10(min(0.01, Reject.threshold - max.LR0_n)))
          termination.threshold.AR = (floor(max.LR0_n*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01, min(term.thresh.possible.choices) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(term.thresh.possible.choices))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(term.thresh.possible.choices))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(cumRejFreq_not.reached.decisionH0.AR[1]>nNewRejects.AR){
            
            nDecimal.accuracy =
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR =
              (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                       (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR +
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      ## Expected sample sizes
      EN0 = mean(N0.AR)     # under H0
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", round(termination.threshold.AR, 3)))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  'N' = list('H0' = N0.AR), 'EN' = EN0, "UMPBT" = UMPBT,
                  "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneProp', 'side' = side, 'theta0' = theta0, 'sigma' = sigma,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
      
    }else if(is.logical(theta1)&&(theta1==T)){
      
      ################ comparison at the fixed-design alternative (default) ################
      
      theta1 = list('right' = fixed_design.alt(test.type = 'oneProp', side = 'right', 
                                               theta0 = theta0, N = N.max, 
                                               Type1 = Type1.target/2, Type2 = Type2.target),
                    'left' = fixed_design.alt(test.type = 'oneProp', side = 'left', 
                                              theta0 = theta0, N = N.max, 
                                              Type1 = Type1.target/2, Type2 = Type2.target))
      
      ################ UMPBT alternative ################
      UMPBT = list('right' = UMPBT.alt(test.type = 'oneProp', side = 'right', 
                                       theta0 = theta0, N = N.max, Type1 = Type1.target/2),
                   'left' = UMPBT.alt(test.type = 'oneProp', side = 'left',
                                      theta0 = theta0, N = N.max, Type1 = Type1.target/2))
      
      # msg
      if(verbose==T){
        
        print("Alternative under comparison:")
        print(paste(' On the right:, ', round(theta1$right, 3), sep = ""))
        print(paste(' On the left:, ', round(theta1$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(UMPBT$right$theta[1], 3), " & ",
                    round(UMPBT$right$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$right$mix.prob[1], 3), " & ", 1 - round(UMPBT$right$mix.prob[1], 3),
                    sep = ""))
        print(paste(' On the left: ', round(UMPBT$left$theta[1], 3), " & ",
                    round(UMPBT$left$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$left$mix.prob[1], 3), " & ", 1 - round(UMPBT$left$mix.prob[1], 3),
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # required storages
      cumsum0_n = cumsum1r_n = cumsum1l_n =
        LR0_n.r = LR0_n.l = LR1r_n.r = LR1r_n.l = LR1l_n.r = LR1l_n.l = numeric(nReplicate)
      type1.error.AR = PowerH1r.AR = PowerH1l.AR = rep(F, nReplicate)
      N0.AR = N0.AR.r = N0.AR.l = N1r.AR = N1r.AR.r = N1r.AR.l = 
        N1l.AR = N1l.AR.r = N1l.AR.l = rep(N.max, nReplicate)
      decision.underH0.AR.r = decision.underH0.AR.l = 
        decision.underH1r.AR.r = decision.underH1r.AR.l = 
        decision.underH1l.AR.r = decision.underH1l.AR.l = rep(NA, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.l =
        not.reached.decisionH1r.AR = not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.l =
        not.reached.decisionH1l.AR = not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.l =
        1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # sum of observations at step n
          sum0_n = rbinom(length(not.reached.decisionH0.AR),
                          batch.size[n+1]-batch.size[n], theta0)
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + sum0_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR0_n.r[not.reached.decisionH0.AR.r] = 
            UMPBT$right$mix.prob[1]*(((1 - UMPBT$right$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[1])))^cumsum0_n[not.reached.decisionH0.AR.r] +
            (1 - UMPBT$right$mix.prob[2])*(((1 - UMPBT$right$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[2])))^cumsum0_n[not.reached.decisionH0.AR.r]
          
          # for left sided check
          LR0_n.l[not.reached.decisionH0.AR.l] = 
            UMPBT$left$mix.prob[1]*(((1 - UMPBT$left$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[1])))^cumsum0_n[not.reached.decisionH0.AR.l] +
            (1 - UMPBT$left$mix.prob[2])*(((1 - UMPBT$left$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[2])))^cumsum0_n[not.reached.decisionH0.AR.l]
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]<=Accept.threshold
          RejectedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]>=Reject.threshold
          reached.decisionH0_n.AR.r = AcceptedH0.underH0_n.AR.r|RejectedH0.underH0_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.r)){
            
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[AcceptedH0.underH0_n.AR.r]] = 'A'
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[RejectedH0.underH0_n.AR.r]] = 'R'
            N0.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.r[!reached.decisionH0_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]<=Accept.threshold
          RejectedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]>=Reject.threshold
          reached.decisionH0_n.AR.l = AcceptedH0.underH0_n.AR.l|RejectedH0.underH0_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.l)){
            
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[AcceptedH0.underH0_n.AR.l]] = 'A'
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[RejectedH0.underH0_n.AR.l]] = 'R'
            N0.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH0.AR.l = not.reached.decisionH0.AR.l[!reached.decisionH0_n.AR.l]
          }
          
          not.reached.decisionH0.AR = union(not.reached.decisionH0.AR.r,
                                            not.reached.decisionH0.AR.l)
        }
        
        
        ## under right-sided H1
        if(length(not.reached.decisionH1r.AR)>0){
          
          # sum of observations at step n
          sum1r_n = rbinom(length(not.reached.decisionH1r.AR),
                           batch.size[n+1]-batch.size[n], theta1$right)
          
          # sum of observations until step n
          cumsum1r_n[not.reached.decisionH1r.AR] =
            cumsum1r_n[not.reached.decisionH1r.AR] + sum1r_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1r_n.r[not.reached.decisionH1r.AR.r] = 
            UMPBT$right$mix.prob[1]*(((1 - UMPBT$right$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[1])))^cumsum1r_n[not.reached.decisionH1r.AR.r] +
            (1 - UMPBT$right$mix.prob[2])*(((1 - UMPBT$right$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[2])))^cumsum1r_n[not.reached.decisionH1r.AR.r]
          
          # for left sided check
          LR1r_n.l[not.reached.decisionH1r.AR.l] = 
            UMPBT$left$mix.prob[1]*(((1 - UMPBT$left$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[1])))^cumsum1r_n[not.reached.decisionH1r.AR.l] +
            (1 - UMPBT$left$mix.prob[2])*(((1 - UMPBT$left$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[2])))^cumsum1r_n[not.reached.decisionH1r.AR.l]
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]<=Accept.threshold
          RejectedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]>=Reject.threshold
          reached.decisionH1r_n.AR.r = AcceptedH0.underH1r_n.AR.r|RejectedH0.underH1r_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.r)){
            
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[AcceptedH0.underH1r_n.AR.r]] = 'A'
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[RejectedH0.underH1r_n.AR.r]] = 'R'
            N1r.AR.r[not.reached.decisionH1r.AR.r[reached.decisionH1r_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.r[!reached.decisionH1r_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]<=Accept.threshold
          RejectedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]>=Reject.threshold
          reached.decisionH1r_n.AR.l = AcceptedH0.underH1r_n.AR.l|RejectedH0.underH1r_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.l)){
            
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[AcceptedH0.underH1r_n.AR.l]] = 'A'
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[RejectedH0.underH1r_n.AR.l]] = 'R'
            N1r.AR.l[not.reached.decisionH1r.AR.l[reached.decisionH1r_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH1r.AR.l = not.reached.decisionH1r.AR.l[!reached.decisionH1r_n.AR.l]
          }
          
          not.reached.decisionH1r.AR = union(not.reached.decisionH1r.AR.r,
                                             not.reached.decisionH1r.AR.l)
        }
        
        
        ## under left-sided H1
        if(length(not.reached.decisionH1l.AR)>0){
          
          # sum of observations at step n
          sum1l_n = rbinom(length(not.reached.decisionH1l.AR),
                           batch.size[n+1]-batch.size[n], theta1$left)
          
          # sum of observations until step n
          cumsum1l_n[not.reached.decisionH1l.AR] =
            cumsum1l_n[not.reached.decisionH1l.AR] + sum1l_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1l_n.r[not.reached.decisionH1l.AR.r] = 
            UMPBT$right$mix.prob[1]*(((1 - UMPBT$right$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[1])))^cumsum1l_n[not.reached.decisionH1l.AR.r] +
            (1 - UMPBT$right$mix.prob[2])*(((1 - UMPBT$right$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[2])))^cumsum1l_n[not.reached.decisionH1l.AR.r]
          
          # for left sided check
          LR1l_n.l[not.reached.decisionH1l.AR.l] = 
            UMPBT$left$mix.prob[1]*(((1 - UMPBT$left$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[1])))^cumsum1l_n[not.reached.decisionH1l.AR.l] +
            (1 - UMPBT$left$mix.prob[2])*(((1 - UMPBT$left$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[2])))^cumsum1l_n[not.reached.decisionH1l.AR.l]
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]<=Accept.threshold
          RejectedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]>=Reject.threshold
          reached.decisionH1l_n.AR.r = AcceptedH0.underH1l_n.AR.r|RejectedH0.underH1l_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.r)){
            
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[AcceptedH0.underH1l_n.AR.r]] = 'A'
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[RejectedH0.underH1l_n.AR.r]] = 'R'
            N1l.AR.r[not.reached.decisionH1l.AR.r[reached.decisionH1l_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.r[!reached.decisionH1l_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]<=Accept.threshold
          RejectedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]>=Reject.threshold
          reached.decisionH1l_n.AR.l = AcceptedH0.underH1l_n.AR.l|RejectedH0.underH1l_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.l)){
            
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[AcceptedH0.underH1l_n.AR.l]] = 'A'
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[RejectedH0.underH1l_n.AR.l]] = 'R'
            N1l.AR.l[not.reached.decisionH1l.AR.l[reached.decisionH1l_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH1l.AR.l = not.reached.decisionH1l.AR.l[!reached.decisionH1l_n.AR.l]
          }
          
          not.reached.decisionH1l.AR = union(not.reached.decisionH1l.AR.r,
                                             not.reached.decisionH1l.AR.l)
        }
        
        setTxtProgressBar(pb, n)
      }
      
      
      ### both-sided checking
      ## under H0
      # accepted or rejected ones
      accepted.by.both0 = intersect(which(decision.underH0.AR.r=='A'),
                                    which(decision.underH0.AR.l=='A'))
      onlyrejected.by.right0 = intersect(which(decision.underH0.AR.r=='R'),
                                         which(decision.underH0.AR.l!='R'))
      onlyrejected.by.left0 = intersect(which(decision.underH0.AR.r!='R'),
                                        which(decision.underH0.AR.l=='R'))
      rejected.by.both0 = intersect(which(decision.underH0.AR.r=='R'),
                                    which(decision.underH0.AR.l=='R'))
      
      # sample sizes required
      N0.AR[accepted.by.both0] = pmax(N0.AR.r[accepted.by.both0],
                                      N0.AR.l[accepted.by.both0])
      N0.AR[onlyrejected.by.right0] = N0.AR.r[onlyrejected.by.right0]
      N0.AR[onlyrejected.by.left0] = N0.AR.l[onlyrejected.by.left0]
      N0.AR[rejected.by.both0] = pmin(N0.AR.r[rejected.by.both0],
                                      N0.AR.l[rejected.by.both0])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right0 = intersect(which(decision.underH0.AR.r=='A'),
                                         which(is.na(decision.underH0.AR.l)))
      onlyaccepted.by.left0 = intersect(which(is.na(decision.underH0.AR.r)),
                                        which(decision.underH0.AR.l=='A'))
      both.inconclusive0 = intersect(which(is.na(decision.underH0.AR.r)),
                                     which(is.na(decision.underH0.AR.l)))
      all.inconclusive0 = c(onlyaccepted.by.right0, onlyaccepted.by.left0,
                            both.inconclusive0)
      nNot.reached.decisionH0.AR = length(all.inconclusive0)
      
      # Type I error probability
      type1.error.AR[c(onlyrejected.by.right0, onlyrejected.by.left0,
                       rejected.by.both0)] = T
      
      
      ## under right-sided H1
      # accepted or rejected ones
      accepted.by.both1r = intersect(which(decision.underH1r.AR.r=='A'),
                                     which(decision.underH1r.AR.l=='A'))
      onlyrejected.by.right1r = intersect(which(decision.underH1r.AR.r=='R'),
                                          which(decision.underH1r.AR.l!='R'))
      onlyrejected.by.left1r = intersect(which(decision.underH1r.AR.r!='R'),
                                         which(decision.underH1r.AR.l=='R'))
      rejected.by.both1r = intersect(which(decision.underH1r.AR.r=='R'),
                                     which(decision.underH1r.AR.l=='R'))
      
      # sample sizes required
      N1r.AR[accepted.by.both1r] = pmax(N1r.AR.r[accepted.by.both1r],
                                        N1r.AR.l[accepted.by.both1r])
      N1r.AR[onlyrejected.by.right1r] = N1r.AR.r[onlyrejected.by.right1r]
      N1r.AR[onlyrejected.by.left1r] = N1r.AR.l[onlyrejected.by.left1r]
      N1r.AR[rejected.by.both1r] = pmin(N1r.AR.r[rejected.by.both1r],
                                        N1r.AR.l[rejected.by.both1r])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1r = intersect(which(decision.underH1r.AR.r=='A'),
                                          which(is.na(decision.underH1r.AR.l)))
      onlyaccepted.by.left1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                         which(decision.underH1r.AR.l=='A'))
      both.inconclusive1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                      which(is.na(decision.underH1r.AR.l)))
      all.inconclusive1r = c(onlyaccepted.by.right1r, onlyaccepted.by.left1r,
                             both.inconclusive1r)
      nNot.reached.decisionH1r.AR = length(all.inconclusive1r)
      
      # Type I error probability
      PowerH1r.AR[c(onlyrejected.by.right1r, onlyrejected.by.left1r,
                    rejected.by.both1r)] = T
      
      
      ## under left-sided H1
      # accepted or rejected ones
      accepted.by.both1l = intersect(which(decision.underH1l.AR.r=='A'),
                                     which(decision.underH1l.AR.l=='A'))
      onlyrejected.by.right1l = intersect(which(decision.underH1l.AR.r=='R'),
                                          which(decision.underH1l.AR.l!='R'))
      onlyrejected.by.left1l = intersect(which(decision.underH1l.AR.r!='R'),
                                         which(decision.underH1l.AR.l=='R'))
      rejected.by.both1l = intersect(which(decision.underH1l.AR.r=='R'),
                                     which(decision.underH1l.AR.l=='R'))
      
      # sample sizes required
      N1l.AR[accepted.by.both1l] = pmax(N1l.AR.r[accepted.by.both1l],
                                        N1l.AR.l[accepted.by.both1l])
      N1l.AR[onlyrejected.by.right1l] = N1l.AR.r[onlyrejected.by.right1l]
      N1l.AR[onlyrejected.by.left1l] = N1l.AR.l[onlyrejected.by.left1l]
      N1l.AR[rejected.by.both1l] = pmin(N1l.AR.r[rejected.by.both1l],
                                        N1l.AR.l[rejected.by.both1l])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1l = intersect(which(decision.underH1l.AR.r=='A'),
                                          which(is.na(decision.underH1l.AR.l)))
      onlyaccepted.by.left1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                         which(decision.underH1l.AR.l=='A'))
      both.inconclusive1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                      which(is.na(decision.underH1l.AR.l)))
      all.inconclusive1l = c(onlyaccepted.by.right1l, onlyaccepted.by.left1l,
                             both.inconclusive1l)
      nNot.reached.decisionH1l.AR = length(all.inconclusive1l)
      
      # Type I error probability
      PowerH1l.AR[c(onlyrejected.by.right1l, onlyrejected.by.left1l,
                    rejected.by.both1l)] = T
      
      
      ## determining termination threshold
      ## H0 is rejected if LR or (BF) is >= termination threshold
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        term.thresh.possible.choices =
          c(LR0_n.r[onlyaccepted.by.left0],
            LR0_n.l[onlyaccepted.by.right0],
            pmin(LR0_n.r[both.inconclusive0], LR0_n.l[both.inconclusive0]))
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          max.LR0_n = max(term.thresh.possible.choices)
          nDecimal.accuracy = ceiling(-log10(min(0.01, Reject.threshold - max.LR0_n)))
          termination.threshold.AR = (floor(max.LR0_n*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01, min(term.thresh.possible.choices) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(term.thresh.possible.choices))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(term.thresh.possible.choices))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(cumRejFreq_not.reached.decisionH0.AR[1]>nNewRejects.AR){
            
            nDecimal.accuracy =
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR =
              (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                       (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR +
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      ## attained Type II error probability
      # right-sided H1
      actual.PowerH1r.AR.r = mean(PowerH1r.AR) +
        sum(c(LR1r_n.r[onlyaccepted.by.left1r],
              LR1r_n.l[onlyaccepted.by.right1r],
              pmax(LR1r_n.r[both.inconclusive1r], LR1r_n.l[both.inconclusive1r]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1r.AR = 1 - actual.PowerH1r.AR.r
      
      # left-sided H1
      actual.PowerH1l.AR.r = mean(PowerH1l.AR) +
        sum(c(LR1l_n.r[onlyaccepted.by.left1l],
              LR1l_n.l[onlyaccepted.by.right1l],
              pmax(LR1l_n.r[both.inconclusive1l], LR1l_n.l[both.inconclusive1l]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1l.AR = 1 - actual.PowerH1l.AR.r
      
      ## Expected sample sizes
      EN0 = mean(N0.AR)     # under H0
      EN1r = mean(N1r.AR)   # under right-sided H1
      EN1l = mean(N1l.AR)   # under left-sided H1
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", round(termination.threshold.AR, 3)))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print("Attained Type II error probability:")
        print(paste(" On the right: ", round(actual.type2.errorH1r.AR, 4)))
        print(paste(" On the left: ", round(actual.type2.errorH1l.AR, 4)))
        print("Expected sample size at the alternatives:")
        print(paste(" On the right: ", round(EN1r, 2)))
        print(paste(" On the left: ", round(EN1l, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  "Type2.attained" = c(actual.type2.errorH1r.AR, actual.type2.errorH1l.AR),
                  'N' = list('H0' = N0.AR, 'right' = N1r.AR, 'left' = N1l.AR),
                  'EN' = c(EN0, EN1r, EN1l), "UMPBT" = UMPBT,
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneProp', 'side' = side, 'theta0' = theta0, 'sigma' = sigma,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
      
    }else{
      
      ################ comparison at user provided point alternative ################
      
      ################ UMPBT alternative ################
      UMPBT = list('right' = UMPBT.alt(test.type = 'oneProp', side = 'right', 
                                       theta0 = theta0, N = N.max, Type1 = Type1.target/2),
                   'left' = UMPBT.alt(test.type = 'oneProp', side = 'left',
                                      theta0 = theta0, N = N.max, Type1 = Type1.target/2))
      
      # msg
      if(verbose==T){
        
        print("Alternative under comparison:")
        print(paste(' On the right: ', round(theta1$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta1$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(UMPBT$right$theta[1], 3), " & ",
                    round(UMPBT$right$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$right$mix.prob[1], 3), " & ", 1 - round(UMPBT$right$mix.prob[1], 3),
                    sep = ""))
        print(paste(' On the left: ', round(UMPBT$left$theta[1], 3), " & ",
                    round(UMPBT$left$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$left$mix.prob[1], 3), " & ", 1 - round(UMPBT$left$mix.prob[1], 3),
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # required storages
      cumsum0_n = cumsum1r_n = cumsum1l_n =
        LR0_n.r = LR0_n.l = LR1r_n.r = LR1r_n.l = LR1l_n.r = LR1l_n.l = numeric(nReplicate)
      type1.error.AR = PowerH1r.AR = PowerH1l.AR = rep(F, nReplicate)
      N0.AR = N0.AR.r = N0.AR.l = N1r.AR = N1r.AR.r = N1r.AR.l = 
        N1l.AR = N1l.AR.r = N1l.AR.l = rep(N.max, nReplicate)
      decision.underH0.AR.r = decision.underH0.AR.l = 
        decision.underH1r.AR.r = decision.underH1r.AR.l = 
        decision.underH1l.AR.r = decision.underH1l.AR.l = rep(NA, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.l =
        not.reached.decisionH1r.AR = not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.l =
        not.reached.decisionH1l.AR = not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.l =
        1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # sum of observations at step n
          sum0_n = rbinom(length(not.reached.decisionH0.AR),
                          batch.size[n+1]-batch.size[n], theta0)
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + sum0_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR0_n.r[not.reached.decisionH0.AR.r] = 
            UMPBT$right$mix.prob[1]*(((1 - UMPBT$right$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[1])))^cumsum0_n[not.reached.decisionH0.AR.r] +
            (1 - UMPBT$right$mix.prob[2])*(((1 - UMPBT$right$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[2])))^cumsum0_n[not.reached.decisionH0.AR.r]
          
          # for left sided check
          LR0_n.l[not.reached.decisionH0.AR.l] = 
            UMPBT$left$mix.prob[1]*(((1 - UMPBT$left$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[1])))^cumsum0_n[not.reached.decisionH0.AR.l] +
            (1 - UMPBT$left$mix.prob[2])*(((1 - UMPBT$left$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[2])))^cumsum0_n[not.reached.decisionH0.AR.l]
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]<=Accept.threshold
          RejectedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]>=Reject.threshold
          reached.decisionH0_n.AR.r = AcceptedH0.underH0_n.AR.r|RejectedH0.underH0_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.r)){
            
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[AcceptedH0.underH0_n.AR.r]] = 'A'
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[RejectedH0.underH0_n.AR.r]] = 'R'
            N0.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.r[!reached.decisionH0_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]<=Accept.threshold
          RejectedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]>=Reject.threshold
          reached.decisionH0_n.AR.l = AcceptedH0.underH0_n.AR.l|RejectedH0.underH0_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.l)){
            
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[AcceptedH0.underH0_n.AR.l]] = 'A'
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[RejectedH0.underH0_n.AR.l]] = 'R'
            N0.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH0.AR.l = not.reached.decisionH0.AR.l[!reached.decisionH0_n.AR.l]
          }
          
          not.reached.decisionH0.AR = union(not.reached.decisionH0.AR.r,
                                            not.reached.decisionH0.AR.l)
        }
        
        
        ## under right-sided H1
        if(length(not.reached.decisionH1r.AR)>0){
          
          # sum of observations at step n
          sum1r_n = rbinom(length(not.reached.decisionH1r.AR),
                           batch.size[n+1]-batch.size[n], theta1$right)
          
          # sum of observations until step n
          cumsum1r_n[not.reached.decisionH1r.AR] =
            cumsum1r_n[not.reached.decisionH1r.AR] + sum1r_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1r_n.r[not.reached.decisionH1r.AR.r] = 
            UMPBT$right$mix.prob[1]*(((1 - UMPBT$right$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[1])))^cumsum1r_n[not.reached.decisionH1r.AR.r] +
            (1 - UMPBT$right$mix.prob[2])*(((1 - UMPBT$right$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[2])))^cumsum1r_n[not.reached.decisionH1r.AR.r]
          
          # for left sided check
          LR1r_n.l[not.reached.decisionH1r.AR.l] = 
            UMPBT$left$mix.prob[1]*(((1 - UMPBT$left$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[1])))^cumsum1r_n[not.reached.decisionH1r.AR.l] +
            (1 - UMPBT$left$mix.prob[2])*(((1 - UMPBT$left$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[2])))^cumsum1r_n[not.reached.decisionH1r.AR.l]
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]<=Accept.threshold
          RejectedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]>=Reject.threshold
          reached.decisionH1r_n.AR.r = AcceptedH0.underH1r_n.AR.r|RejectedH0.underH1r_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.r)){
            
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[AcceptedH0.underH1r_n.AR.r]] = 'A'
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[RejectedH0.underH1r_n.AR.r]] = 'R'
            N1r.AR.r[not.reached.decisionH1r.AR.r[reached.decisionH1r_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.r[!reached.decisionH1r_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]<=Accept.threshold
          RejectedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]>=Reject.threshold
          reached.decisionH1r_n.AR.l = AcceptedH0.underH1r_n.AR.l|RejectedH0.underH1r_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.l)){
            
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[AcceptedH0.underH1r_n.AR.l]] = 'A'
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[RejectedH0.underH1r_n.AR.l]] = 'R'
            N1r.AR.l[not.reached.decisionH1r.AR.l[reached.decisionH1r_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH1r.AR.l = not.reached.decisionH1r.AR.l[!reached.decisionH1r_n.AR.l]
          }
          
          not.reached.decisionH1r.AR = union(not.reached.decisionH1r.AR.r,
                                             not.reached.decisionH1r.AR.l)
        }
        
        
        ## under left-sided H1
        if(length(not.reached.decisionH1l.AR)>0){
          
          # sum of observations at step n
          sum1l_n = rbinom(length(not.reached.decisionH1l.AR),
                           batch.size[n+1]-batch.size[n], theta1$left)
          
          # sum of observations until step n
          cumsum1l_n[not.reached.decisionH1l.AR] =
            cumsum1l_n[not.reached.decisionH1l.AR] + sum1l_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1l_n.r[not.reached.decisionH1l.AR.r] = 
            UMPBT$right$mix.prob[1]*(((1 - UMPBT$right$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[1])))^cumsum1l_n[not.reached.decisionH1l.AR.r] +
            (1 - UMPBT$right$mix.prob[2])*(((1 - UMPBT$right$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[2])))^cumsum1l_n[not.reached.decisionH1l.AR.r]
          
          # for left sided check
          LR1l_n.l[not.reached.decisionH1l.AR.l] = 
            UMPBT$left$mix.prob[1]*(((1 - UMPBT$left$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[1])))^cumsum1l_n[not.reached.decisionH1l.AR.l] +
            (1 - UMPBT$left$mix.prob[2])*(((1 - UMPBT$left$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[2])))^cumsum1l_n[not.reached.decisionH1l.AR.l]
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]<=Accept.threshold
          RejectedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]>=Reject.threshold
          reached.decisionH1l_n.AR.r = AcceptedH0.underH1l_n.AR.r|RejectedH0.underH1l_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.r)){
            
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[AcceptedH0.underH1l_n.AR.r]] = 'A'
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[RejectedH0.underH1l_n.AR.r]] = 'R'
            N1l.AR.r[not.reached.decisionH1l.AR.r[reached.decisionH1l_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.r[!reached.decisionH1l_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]<=Accept.threshold
          RejectedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]>=Reject.threshold
          reached.decisionH1l_n.AR.l = AcceptedH0.underH1l_n.AR.l|RejectedH0.underH1l_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.l)){
            
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[AcceptedH0.underH1l_n.AR.l]] = 'A'
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[RejectedH0.underH1l_n.AR.l]] = 'R'
            N1l.AR.l[not.reached.decisionH1l.AR.l[reached.decisionH1l_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH1l.AR.l = not.reached.decisionH1l.AR.l[!reached.decisionH1l_n.AR.l]
          }
          
          not.reached.decisionH1l.AR = union(not.reached.decisionH1l.AR.r,
                                             not.reached.decisionH1l.AR.l)
        }
        
        setTxtProgressBar(pb, n)
      }
      
      
      ### both-sided checking
      ## under H0
      # accepted or rejected ones
      accepted.by.both0 = intersect(which(decision.underH0.AR.r=='A'),
                                    which(decision.underH0.AR.l=='A'))
      onlyrejected.by.right0 = intersect(which(decision.underH0.AR.r=='R'),
                                         which(decision.underH0.AR.l!='R'))
      onlyrejected.by.left0 = intersect(which(decision.underH0.AR.r!='R'),
                                        which(decision.underH0.AR.l=='R'))
      rejected.by.both0 = intersect(which(decision.underH0.AR.r=='R'),
                                    which(decision.underH0.AR.l=='R'))
      
      # sample sizes required
      N0.AR[accepted.by.both0] = pmax(N0.AR.r[accepted.by.both0],
                                      N0.AR.l[accepted.by.both0])
      N0.AR[onlyrejected.by.right0] = N0.AR.r[onlyrejected.by.right0]
      N0.AR[onlyrejected.by.left0] = N0.AR.l[onlyrejected.by.left0]
      N0.AR[rejected.by.both0] = pmin(N0.AR.r[rejected.by.both0],
                                      N0.AR.l[rejected.by.both0])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right0 = intersect(which(decision.underH0.AR.r=='A'),
                                         which(is.na(decision.underH0.AR.l)))
      onlyaccepted.by.left0 = intersect(which(is.na(decision.underH0.AR.r)),
                                        which(decision.underH0.AR.l=='A'))
      both.inconclusive0 = intersect(which(is.na(decision.underH0.AR.r)),
                                     which(is.na(decision.underH0.AR.l)))
      all.inconclusive0 = c(onlyaccepted.by.right0, onlyaccepted.by.left0,
                            both.inconclusive0)
      nNot.reached.decisionH0.AR = length(all.inconclusive0)
      
      # Type I error probability
      type1.error.AR[c(onlyrejected.by.right0, onlyrejected.by.left0,
                       rejected.by.both0)] = T
      
      
      ## under right-sided H1
      # accepted or rejected ones
      accepted.by.both1r = intersect(which(decision.underH1r.AR.r=='A'),
                                     which(decision.underH1r.AR.l=='A'))
      onlyrejected.by.right1r = intersect(which(decision.underH1r.AR.r=='R'),
                                          which(decision.underH1r.AR.l!='R'))
      onlyrejected.by.left1r = intersect(which(decision.underH1r.AR.r!='R'),
                                         which(decision.underH1r.AR.l=='R'))
      rejected.by.both1r = intersect(which(decision.underH1r.AR.r=='R'),
                                     which(decision.underH1r.AR.l=='R'))
      
      # sample sizes required
      N1r.AR[accepted.by.both1r] = pmax(N1r.AR.r[accepted.by.both1r],
                                        N1r.AR.l[accepted.by.both1r])
      N1r.AR[onlyrejected.by.right1r] = N1r.AR.r[onlyrejected.by.right1r]
      N1r.AR[onlyrejected.by.left1r] = N1r.AR.l[onlyrejected.by.left1r]
      N1r.AR[rejected.by.both1r] = pmin(N1r.AR.r[rejected.by.both1r],
                                        N1r.AR.l[rejected.by.both1r])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1r = intersect(which(decision.underH1r.AR.r=='A'),
                                          which(is.na(decision.underH1r.AR.l)))
      onlyaccepted.by.left1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                         which(decision.underH1r.AR.l=='A'))
      both.inconclusive1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                      which(is.na(decision.underH1r.AR.l)))
      all.inconclusive1r = c(onlyaccepted.by.right1r, onlyaccepted.by.left1r,
                             both.inconclusive1r)
      nNot.reached.decisionH1r.AR = length(all.inconclusive1r)
      
      # Type I error probability
      PowerH1r.AR[c(onlyrejected.by.right1r, onlyrejected.by.left1r,
                    rejected.by.both1r)] = T
      
      
      ## under left-sided H1
      # accepted or rejected ones
      accepted.by.both1l = intersect(which(decision.underH1l.AR.r=='A'),
                                     which(decision.underH1l.AR.l=='A'))
      onlyrejected.by.right1l = intersect(which(decision.underH1l.AR.r=='R'),
                                          which(decision.underH1l.AR.l!='R'))
      onlyrejected.by.left1l = intersect(which(decision.underH1l.AR.r!='R'),
                                         which(decision.underH1l.AR.l=='R'))
      rejected.by.both1l = intersect(which(decision.underH1l.AR.r=='R'),
                                     which(decision.underH1l.AR.l=='R'))
      
      # sample sizes required
      N1l.AR[accepted.by.both1l] = pmax(N1l.AR.r[accepted.by.both1l],
                                        N1l.AR.l[accepted.by.both1l])
      N1l.AR[onlyrejected.by.right1l] = N1l.AR.r[onlyrejected.by.right1l]
      N1l.AR[onlyrejected.by.left1l] = N1l.AR.l[onlyrejected.by.left1l]
      N1l.AR[rejected.by.both1l] = pmin(N1l.AR.r[rejected.by.both1l],
                                        N1l.AR.l[rejected.by.both1l])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1l = intersect(which(decision.underH1l.AR.r=='A'),
                                          which(is.na(decision.underH1l.AR.l)))
      onlyaccepted.by.left1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                         which(decision.underH1l.AR.l=='A'))
      both.inconclusive1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                      which(is.na(decision.underH1l.AR.l)))
      all.inconclusive1l = c(onlyaccepted.by.right1l, onlyaccepted.by.left1l,
                             both.inconclusive1l)
      nNot.reached.decisionH1l.AR = length(all.inconclusive1l)
      
      # Type I error probability
      PowerH1l.AR[c(onlyrejected.by.right1l, onlyrejected.by.left1l,
                    rejected.by.both1l)] = T
      
      
      ## determining termination threshold
      ## H0 is rejected if LR or (BF) is >= termination threshold
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        term.thresh.possible.choices =
          c(LR0_n.r[onlyaccepted.by.left0],
            LR0_n.l[onlyaccepted.by.right0],
            pmin(LR0_n.r[both.inconclusive0], LR0_n.l[both.inconclusive0]))
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          max.LR0_n = max(term.thresh.possible.choices)
          nDecimal.accuracy = ceiling(-log10(min(0.01, Reject.threshold - max.LR0_n)))
          termination.threshold.AR = (floor(max.LR0_n*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01, min(term.thresh.possible.choices) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(term.thresh.possible.choices))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(term.thresh.possible.choices))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(cumRejFreq_not.reached.decisionH0.AR[1]>nNewRejects.AR){
            
            nDecimal.accuracy =
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR =
              (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                       (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR +
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      ## attained Type II error probability
      # right-sided H1
      actual.PowerH1r.AR.r = mean(PowerH1r.AR) +
        sum(c(LR1r_n.r[onlyaccepted.by.left1r],
              LR1r_n.l[onlyaccepted.by.right1r],
              pmax(LR1r_n.r[both.inconclusive1r], LR1r_n.l[both.inconclusive1r]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1r.AR = 1 - actual.PowerH1r.AR.r
      
      # left-sided H1
      actual.PowerH1l.AR.r = mean(PowerH1l.AR) +
        sum(c(LR1l_n.r[onlyaccepted.by.left1l],
              LR1l_n.l[onlyaccepted.by.right1l],
              pmax(LR1l_n.r[both.inconclusive1l], LR1l_n.l[both.inconclusive1l]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1l.AR = 1 - actual.PowerH1l.AR.r
      
      ## Expected sample sizes
      EN0 = mean(N0.AR)     # under H0
      EN1r = mean(N1r.AR)   # under right-sided H1
      EN1l = mean(N1l.AR)   # under left-sided H1
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", round(termination.threshold.AR, 3)))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print("Attained Type II error probability:")
        print(paste(" On the right: ", round(actual.type2.errorH1r.AR, 4)))
        print(paste(" On the left: ", round(actual.type2.errorH1l.AR, 4)))
        print("Expected sample size at the alternatives:")
        print(paste(" On the right: ", round(EN1r, 2)))
        print(paste(" On the left: ", round(EN1l, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  "Type2.attained" = c(actual.type2.errorH1r.AR, actual.type2.errorH1l.AR),
                  'N' = list('H0' = N0.AR, 'right' = N1r.AR, 'left' = N1l.AR),
                  'EN' = c(EN0, EN1r, EN1l), "UMPBT" = UMPBT,
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneProp', 'side' = side, 'theta0' = theta0, 'sigma' = sigma,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
    }
  }
}

#### one-sample z test ####
design.MSPRT_oneZ = function(side = 'right', theta0 = 0, theta1 = T,
                             Type1.target =.005, Type2.target = .2,
                             N.max, batch.size, sigma = 1,
                             nReplicate = 1e+6, verbose = T, seed = 1){
  
  
  if(side!='both'){
    
    ################################# one-sample z (right/left sided) #################################
    
    ## batch sizes and N.max
    if(missing(batch.size)){
      
      if(missing(N.max)){
        
        return("Either 'batch.size' or 'N.max' needs to be specified")
        
      }else{batch.size = rep(1, N.max)}
      
    }else{
      
      if(missing(N.max)){
        
        N.max = sum(batch.size)
        
      }else{
        
        if(sum(batch.size)!=N.max) return("Sum of batch sizes should add up to N.max")
      }
    }
    
    nAnalyses = length(batch.size)
    
    ## msg
    if(verbose){
      
      if(any(batch.size>1)){
        
        cat('\n')
        print("=========================================================================")
        print("Designing the group sequential MSPRT for a one-sample z test:")
        print("=========================================================================")
        
      }else{
        
        cat('\n')
        print("=========================================================================")
        print("Designing the sequential MSPRT for a one-sample z test:")
        print("=========================================================================")
      }
      
      print(paste("Maximum available sample size: ", N.max, sep = ""))
      print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
      print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
      print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
      print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
      print(paste("Hypothesized value under H0: ", theta0, sep = ""))
      print(paste("Direction of the H1: ", side, sep = ""))
      print(paste("Known standard deviation: ", sigma, sep = ""))
    }
    
    batch.size = c(0, cumsum(batch.size))
    
    
    if(is.logical(theta1)&&(theta1==F)){
      
      ################ no fixed-design alternative ################
      
      ################ UMPBT alternative ################
      theta.UMPBT = UMPBT.alt(test.type = 'oneZ', side = side, theta0 = theta0,
                              N = N.max, Type1 = Type1.target, sigma = sigma)
      
      # msg
      if(verbose==T){
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(theta.UMPBT, 3)))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # required storages
      cumsum0_n = LR0_n = numeric(nReplicate)
      type1.error.AR = rep(F, nReplicate)
      N0.AR = rep(N.max, nReplicate)
      not.reached.decisionH0.AR = 1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # sum of observations at step n
          sum0_n = rnorm(length(not.reached.decisionH0.AR),
                         (batch.size[n+1]-batch.size[n])*theta0,
                         sqrt(batch.size[n+1]-batch.size[n])*sigma)
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + sum0_n
          
          # likelihood ratio of observations until step n
          LR0_n[not.reached.decisionH0.AR] = 
            exp((cumsum0_n[not.reached.decisionH0.AR]*(theta.UMPBT - theta0) - 
                   ((batch.size[n+1]*((theta.UMPBT^2) - (theta0^2)))/2))/(sigma^2))
          
          # comparing with the thresholds
          AcceptedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]<=Accept.threshold)
          RejectedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]>=Reject.threshold)
          reached.decisionH0_n.AR = union(AcceptedH0.underH0_n.AR, RejectedH0.underH0_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH0_n.AR)>0){
            
            N0.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch.size[n+1]
            type1.error.AR[not.reached.decisionH0.AR[RejectedH0.underH0_n.AR]] = T
            not.reached.decisionH0.AR = not.reached.decisionH0.AR[-reached.decisionH0_n.AR]
          }
        }
        
        setTxtProgressBar(pb, n)
      }
      
      # determining termination threshold
      # H0 is rejected if LR or (BF) is >= termination threshold
      nNot.reached.decisionH0.AR = length(not.reached.decisionH0.AR)
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 Reject.threshold -
                                                   max(LR0_n[not.reached.decisionH0.AR]))))
          termination.threshold.AR = (floor(max(LR0_n[not.reached.decisionH0.AR])*
                                              (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 min(LR0_n[not.reached.decisionH0.AR]) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(LR0_n[not.reached.decisionH0.AR]))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(LR0_n[not.reached.decisionH0.AR]))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(min(cumRejFreq_not.reached.decisionH0.AR)>nNewRejects.AR){
            
            nDecimal.accuracy = 
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR + 
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      # Expected sample sizes
      EN0 = mean(N0.AR)
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", Reject.threshold))
        print(paste("Termination threshold: ", termination.threshold.AR))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("TypeI.attained" = actual.type1.error.AR,
                  "N" = list('H0' = N0.AR), "EN" = EN0,
                  "theta.UMPBT" = theta.UMPBT, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneZ', 'side' = side, 'theta0' = theta0, 'sigma' = sigma,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
      
    }else if(is.logical(theta1)&&(theta1==T)){
      
      ################ comparison at the fixed-design alternative (default) ################
      theta1 = fixed_design.alt(test.type = 'oneZ', side = side, theta0 = theta0,
                                N = N.max, Type1 = Type1.target, Type2 = Type2.target,
                                sigma = sigma)
      
      ################ UMPBT alternative ################
      theta.UMPBT = UMPBT.alt(test.type = 'oneZ', side = side, theta0 = theta0,
                              N = N.max, Type1 = Type1.target, sigma = sigma)
      
      # msg
      if(verbose==T){
        
        print(paste("Alternative under comparison: ", round(theta1, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(theta.UMPBT, 3)))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # required storages
      cumsum0_n = cumsum1_n = LR0_n = LR1_n = numeric(nReplicate)
      type1.error.AR = type2.error.AR = rep(F, nReplicate)
      N0.AR = N1.AR = rep(N.max, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH1.AR = 1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # sum of observations at step n
          sum0_n = rnorm(length(not.reached.decisionH0.AR),
                         (batch.size[n+1]-batch.size[n])*theta0,
                         sqrt(batch.size[n+1]-batch.size[n])*sigma)
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + sum0_n
          
          # likelihood ratio of observations until step n
          LR0_n[not.reached.decisionH0.AR] = 
            exp((cumsum0_n[not.reached.decisionH0.AR]*(theta.UMPBT - theta0) - 
                   ((batch.size[n+1]*((theta.UMPBT^2) - (theta0^2)))/2))/(sigma^2))
          
          # comparing with the thresholds
          AcceptedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]<=Accept.threshold)
          RejectedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]>=Reject.threshold)
          reached.decisionH0_n.AR = union(AcceptedH0.underH0_n.AR, RejectedH0.underH0_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH0_n.AR)>0){
            
            N0.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch.size[n+1]
            type1.error.AR[not.reached.decisionH0.AR[RejectedH0.underH0_n.AR]] = T
            not.reached.decisionH0.AR = not.reached.decisionH0.AR[-reached.decisionH0_n.AR]
          }
        }
        
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          # sum of observations at step n
          sum1_n = rnorm(length(not.reached.decisionH1.AR),
                         (batch.size[n+1]-batch.size[n])*theta1,
                         sqrt(batch.size[n+1]-batch.size[n])*sigma)
          
          # sum of observations until step n
          cumsum1_n[not.reached.decisionH1.AR] = 
            cumsum1_n[not.reached.decisionH1.AR] + sum1_n
          
          # likelihood ratio of observations until step n
          LR1_n[not.reached.decisionH1.AR] = 
            exp((cumsum1_n[not.reached.decisionH1.AR]*(theta.UMPBT - theta0) - 
                   ((batch.size[n+1]*((theta.UMPBT^2) - (theta0^2)))/2))/(sigma^2))
          
          # comparing with the thresholds
          AcceptedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]<=Accept.threshold)
          RejectedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]>=Reject.threshold)
          reached.decisionH1_n.AR = union(AcceptedH0.underH1_n.AR, RejectedH0.underH1_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH1_n.AR)>0){
            
            N1.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch.size[n+1]
            type2.error.AR[not.reached.decisionH1.AR[AcceptedH0.underH1_n.AR]] = T
            not.reached.decisionH1.AR = not.reached.decisionH1.AR[-reached.decisionH1_n.AR]
          }
        }
        
        setTxtProgressBar(pb, n)
      }
      
      # determining termination threshold
      # H0 is rejected if LR or (BF) is >= termination threshold
      nNot.reached.decisionH0.AR = length(not.reached.decisionH0.AR)
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 Reject.threshold -
                                                   max(LR0_n[not.reached.decisionH0.AR]))))
          termination.threshold.AR = (floor(max(LR0_n[not.reached.decisionH0.AR])*
                                              (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 min(LR0_n[not.reached.decisionH0.AR]) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(LR0_n[not.reached.decisionH0.AR]))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(LR0_n[not.reached.decisionH0.AR]))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(min(cumRejFreq_not.reached.decisionH0.AR)>nNewRejects.AR){
            
            nDecimal.accuracy = 
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR + 
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      # attained Type II error probability
      actual.type2.error.AR = mean(type2.error.AR) +
        sum(LR1_n[not.reached.decisionH1.AR]<termination.threshold.AR)/nReplicate
      
      # Expected sample sizes
      EN0 = mean(N0.AR)
      EN1 = mean(N1.AR)
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", Reject.threshold))
        print(paste("Termination threshold: ", termination.threshold.AR))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Attained Type II error probability: ", round(actual.type2.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print(paste("Expected sample size at the alternative: ", round(EN1, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("TypeI.attained" = actual.type1.error.AR,
                  "TypeII.attained" = actual.type2.error.AR,
                  "N" = list('H0' = N0.AR, 'H1' = N1.AR), "EN" = c(EN0, EN1),
                  "theta.UMPBT" = theta.UMPBT, "theta1" = theta1,
                  "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneZ', 'side' = side, 'theta0' = theta0, 'sigma' = sigma,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
      
    }else{
      
      ################ comparison at user provided point alternative ################
      
      ################ UMPBT alternative ################
      theta.UMPBT = UMPBT.alt(test.type = 'oneZ', side = side, theta0 = theta0,
                              N = N.max, Type1 = Type1.target, sigma = sigma)
      
      # msg
      if(verbose==T){
        
        print(paste("Alternative under comparison: ", round(theta1, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(theta.UMPBT, 3)))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # required storages
      cumsum0_n = cumsum1_n = LR0_n = LR1_n = numeric(nReplicate)
      type1.error.AR = type2.error.AR = rep(F, nReplicate)
      N0.AR = N1.AR = rep(N.max, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH1.AR = 1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # sum of observations at step n
          sum0_n = rnorm(length(not.reached.decisionH0.AR),
                         (batch.size[n+1]-batch.size[n])*theta0,
                         sqrt(batch.size[n+1]-batch.size[n])*sigma)
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + sum0_n
          
          # likelihood ratio of observations until step n
          LR0_n[not.reached.decisionH0.AR] = 
            exp((cumsum0_n[not.reached.decisionH0.AR]*(theta.UMPBT - theta0) - 
                   ((batch.size[n+1]*((theta.UMPBT^2) - (theta0^2)))/2))/(sigma^2))
          
          # comparing with the thresholds
          AcceptedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]<=Accept.threshold)
          RejectedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]>=Reject.threshold)
          reached.decisionH0_n.AR = union(AcceptedH0.underH0_n.AR, RejectedH0.underH0_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH0_n.AR)>0){
            
            N0.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch.size[n+1]
            type1.error.AR[not.reached.decisionH0.AR[RejectedH0.underH0_n.AR]] = T
            not.reached.decisionH0.AR = not.reached.decisionH0.AR[-reached.decisionH0_n.AR]
          }
        }
        
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          # sum of observations at step n
          sum1_n = rnorm(length(not.reached.decisionH1.AR),
                         (batch.size[n+1]-batch.size[n])*theta1,
                         sqrt(batch.size[n+1]-batch.size[n])*sigma)
          
          # sum of observations until step n
          cumsum1_n[not.reached.decisionH1.AR] = 
            cumsum1_n[not.reached.decisionH1.AR] + sum1_n
          
          # likelihood ratio of observations until step n
          LR1_n[not.reached.decisionH1.AR] = 
            exp((cumsum1_n[not.reached.decisionH1.AR]*(theta.UMPBT - theta0) - 
                   ((batch.size[n+1]*((theta.UMPBT^2) - (theta0^2)))/2))/(sigma^2))
          
          # comparing with the thresholds
          AcceptedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]<=Accept.threshold)
          RejectedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]>=Reject.threshold)
          reached.decisionH1_n.AR = union(AcceptedH0.underH1_n.AR, RejectedH0.underH1_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH1_n.AR)>0){
            
            N1.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch.size[n+1]
            type2.error.AR[not.reached.decisionH1.AR[AcceptedH0.underH1_n.AR]] = T
            not.reached.decisionH1.AR = not.reached.decisionH1.AR[-reached.decisionH1_n.AR]
          }
        }
        
        setTxtProgressBar(pb, n)
      }
      
      # determining termination threshold
      # H0 is rejected if LR or (BF) is >= termination threshold
      nNot.reached.decisionH0.AR = length(not.reached.decisionH0.AR)
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 Reject.threshold -
                                                   max(LR0_n[not.reached.decisionH0.AR]))))
          termination.threshold.AR = (floor(max(LR0_n[not.reached.decisionH0.AR])*
                                              (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 min(LR0_n[not.reached.decisionH0.AR]) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(LR0_n[not.reached.decisionH0.AR]))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(LR0_n[not.reached.decisionH0.AR]))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(min(cumRejFreq_not.reached.decisionH0.AR)>nNewRejects.AR){
            
            nDecimal.accuracy = 
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR + 
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      # attained Type II error probability
      actual.type2.error.AR = mean(type2.error.AR) +
        sum(LR1_n[not.reached.decisionH1.AR]<termination.threshold.AR)/nReplicate
      
      # Expected sample sizes
      EN0 = mean(N0.AR)
      EN1 = mean(N1.AR)
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", Reject.threshold))
        print(paste("Termination threshold: ", termination.threshold.AR))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Attained Type II error probability: ", round(actual.type2.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print(paste("Expected sample size at the alternative: ", round(EN1, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("TypeI.attained" = actual.type1.error.AR,
                  "TypeII.attained" = actual.type2.error.AR,
                  "N" = list('H0' = N0.AR, 'H1' = N1.AR), "EN" = c(EN0, EN1),
                  "theta.UMPBT" = theta.UMPBT, "theta1" = theta1,
                  "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneZ', 'side' = side, 'theta0' = theta0, 'sigma' = sigma,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
    }
    
  }else{
    
    ################################# one-sample z (both sided) #################################
    
    ## batch sizes and N.max
    if(missing(batch.size)){
      
      if(missing(N.max)){
        
        return("Either 'batch.size' or 'N.max' needs to be specified")
        
      }else{batch.size = rep(1, N.max)}
      
    }else{
      
      if(missing(N.max)){
        
        N.max = sum(batch.size)
        
      }else{
        
        if(sum(batch.size)!=N.max) return("Sum of batch sizes should add up to N.max")
      }
    }
    
    nAnalyses = length(batch.size)
    
    ## msg
    if(verbose){
      
      if(any(batch.size>1)){
        
        cat('\n')
        print("=========================================================================")
        print("Designing the group sequential MSPRT for a one-sample z test:")
        print("=========================================================================")
        
      }else{
        
        cat('\n')
        print("=========================================================================")
        print("Designing the sequential MSPRT for a one-sample z test:")
        print("=========================================================================")
      }
      
      print(paste("Maximum available sample size: ", N.max, sep = ""))
      print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
      print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
      print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
      print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
      print(paste("Hypothesized value under H0: ", theta0, sep = ""))
      print(paste("Direction of the H1: ", side, sep = ""))
      print(paste("Known standard deviation: ", sigma, sep = ""))
    }
    
    batch.size = c(0, cumsum(batch.size))
    
    
    if(is.logical(theta1)&&(theta1==F)){
      
      ################ no fixed-design alternative ################
      
      ################ UMPBT alternative ################
      theta.UMPBT = list('right' = UMPBT.alt(test.type = 'oneZ', side = 'right', 
                                             theta0 = theta0, N = N.max, 
                                             Type1 = Type1.target/2, sigma = sigma),
                         'left' = UMPBT.alt(test.type = 'oneZ', side = 'left', 
                                            theta0 = theta0, N = N.max,
                                            Type1 = Type1.target/2, sigma = sigma))
      
      # msg
      if(verbose==T){
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(theta.UMPBT$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta.UMPBT$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # required storages
      cumsum0_n = LR0_n.r = LR0_n.l = numeric(nReplicate)
      type1.error.AR = rep(F, nReplicate)
      N0.AR = N0.AR.r = N0.AR.l = rep(N.max, nReplicate)
      decision.underH0.AR.r = decision.underH0.AR.l = rep(NA, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.l =
        1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # sum of observations at step n
          sum0_n = rnorm(length(not.reached.decisionH0.AR),
                         (batch.size[n+1]-batch.size[n])*theta0,
                         sqrt(batch.size[n+1]-batch.size[n])*sigma)
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + sum0_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR0_n.r[not.reached.decisionH0.AR.r] = 
            exp((cumsum0_n[not.reached.decisionH0.AR.r]*(theta.UMPBT$right - theta0) - 
                   ((batch.size[n+1]*((theta.UMPBT$right^2) - (theta0^2)))/2))/(sigma^2))
          
          # for left sided check
          LR0_n.l[not.reached.decisionH0.AR.l] = 
            exp((cumsum0_n[not.reached.decisionH0.AR.l]*(theta.UMPBT$left - theta0) - 
                   ((batch.size[n+1]*((theta.UMPBT$left^2) - (theta0^2)))/2))/(sigma^2))
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]<=Accept.threshold
          RejectedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]>=Reject.threshold
          reached.decisionH0_n.AR.r = AcceptedH0.underH0_n.AR.r|RejectedH0.underH0_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.r)){
            
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[AcceptedH0.underH0_n.AR.r]] = 'A'
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[RejectedH0.underH0_n.AR.r]] = 'R'
            N0.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.r[!reached.decisionH0_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]<=Accept.threshold
          RejectedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]>=Reject.threshold
          reached.decisionH0_n.AR.l = AcceptedH0.underH0_n.AR.l|RejectedH0.underH0_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.l)){
            
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[AcceptedH0.underH0_n.AR.l]] = 'A'
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[RejectedH0.underH0_n.AR.l]] = 'R'
            N0.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH0.AR.l = not.reached.decisionH0.AR.l[!reached.decisionH0_n.AR.l]
          }
          
          not.reached.decisionH0.AR = union(not.reached.decisionH0.AR.r,
                                            not.reached.decisionH0.AR.l)
        }
        
        setTxtProgressBar(pb, n)
      }
      
      
      ### both-sided checking
      ## under H0
      # accepted or rejected ones
      accepted.by.both0 = intersect(which(decision.underH0.AR.r=='A'),
                                    which(decision.underH0.AR.l=='A'))
      onlyrejected.by.right0 = intersect(which(decision.underH0.AR.r=='R'),
                                         which(decision.underH0.AR.l!='R'))
      onlyrejected.by.left0 = intersect(which(decision.underH0.AR.r!='R'),
                                        which(decision.underH0.AR.l=='R'))
      rejected.by.both0 = intersect(which(decision.underH0.AR.r=='R'),
                                    which(decision.underH0.AR.l=='R'))
      
      # sample sizes required
      N0.AR[accepted.by.both0] = pmax(N0.AR.r[accepted.by.both0],
                                      N0.AR.l[accepted.by.both0])
      N0.AR[onlyrejected.by.right0] = N0.AR.r[onlyrejected.by.right0]
      N0.AR[onlyrejected.by.left0] = N0.AR.l[onlyrejected.by.left0]
      N0.AR[rejected.by.both0] = pmin(N0.AR.r[rejected.by.both0],
                                      N0.AR.l[rejected.by.both0])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right0 = intersect(which(decision.underH0.AR.r=='A'),
                                         which(is.na(decision.underH0.AR.l)))
      onlyaccepted.by.left0 = intersect(which(is.na(decision.underH0.AR.r)),
                                        which(decision.underH0.AR.l=='A'))
      both.inconclusive0 = intersect(which(is.na(decision.underH0.AR.r)),
                                     which(is.na(decision.underH0.AR.l)))
      all.inconclusive0 = c(onlyaccepted.by.right0, onlyaccepted.by.left0,
                            both.inconclusive0)
      nNot.reached.decisionH0.AR = length(all.inconclusive0)
      
      # Type I error probability
      type1.error.AR[c(onlyrejected.by.right0, onlyrejected.by.left0,
                       rejected.by.both0)] = T
      
      
      ## determining termination threshold
      ## H0 is rejected if LR or (BF) is >= termination threshold
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        term.thresh.possible.choices =
          c(LR0_n.r[onlyaccepted.by.left0],
            LR0_n.l[onlyaccepted.by.right0],
            pmin(LR0_n.r[both.inconclusive0], LR0_n.l[both.inconclusive0]))
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          max.LR0_n = max(term.thresh.possible.choices)
          nDecimal.accuracy = ceiling(-log10(min(0.01, Reject.threshold - max.LR0_n)))
          termination.threshold.AR = (floor(max.LR0_n*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01, min(term.thresh.possible.choices) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(term.thresh.possible.choices))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(term.thresh.possible.choices))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(cumRejFreq_not.reached.decisionH0.AR[1]>nNewRejects.AR){
            
            nDecimal.accuracy =
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR =
              (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                       (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR +
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      ## Expected sample sizes
      EN0 = mean(N0.AR)     # under H0
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", round(termination.threshold.AR, 3)))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  'N' = list('H0' = N0.AR),
                  'EN' = EN0, "theta.UMPBT" = theta.UMPBT, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneZ', 'side' = side, 'theta0' = theta0, 'sigma' = sigma,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
      
    }else if(is.logical(theta1)&&(theta1==T)){
      
      ################ comparison at the fixed-design alternative (default) ################
      
      theta1 = list('right' = fixed_design.alt(test.type = 'oneZ', side = 'right', 
                                               theta0 = theta0, N = N.max, 
                                               Type1 = Type1.target/2, Type2 = Type2.target,
                                               sigma = sigma),
                    'left' = fixed_design.alt(test.type = 'oneZ', side = 'left',
                                              theta0 = theta0, N = N.max, 
                                              Type1 = Type1.target/2, Type2 = Type2.target,
                                              sigma = sigma))
      
      ################ UMPBT alternative ################
      theta.UMPBT = list('right' = UMPBT.alt(test.type = 'oneZ', side = 'right', 
                                             theta0 = theta0, N = N.max, 
                                             Type1 = Type1.target/2, sigma = sigma),
                         'left' = UMPBT.alt(test.type = 'oneZ', side = 'left', 
                                            theta0 = theta0, N = N.max,
                                            Type1 = Type1.target/2, sigma = sigma))
      
      # msg
      if(verbose==T){
        
        print("Alternative under comparison:")
        print(paste(' On the right: ', round(theta1$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta1$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(theta.UMPBT$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta.UMPBT$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # required storages
      cumsum0_n = cumsum1r_n = cumsum1l_n =
        LR0_n.r = LR0_n.l = LR1r_n.r = LR1r_n.l = LR1l_n.r = LR1l_n.l = numeric(nReplicate)
      type1.error.AR = PowerH1r.AR = PowerH1l.AR = rep(F, nReplicate)
      N0.AR = N0.AR.r = N0.AR.l = N1r.AR = N1r.AR.r = N1r.AR.l = 
        N1l.AR = N1l.AR.r = N1l.AR.l = rep(N.max, nReplicate)
      decision.underH0.AR.r = decision.underH0.AR.l = 
        decision.underH1r.AR.r = decision.underH1r.AR.l = 
        decision.underH1l.AR.r = decision.underH1l.AR.l = rep(NA, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.l =
        not.reached.decisionH1r.AR = not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.l =
        not.reached.decisionH1l.AR = not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.l =
        1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # sum of observations at step n
          sum0_n = rnorm(length(not.reached.decisionH0.AR),
                         (batch.size[n+1]-batch.size[n])*theta0,
                         sqrt(batch.size[n+1]-batch.size[n])*sigma)
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + sum0_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR0_n.r[not.reached.decisionH0.AR.r] = 
            exp((cumsum0_n[not.reached.decisionH0.AR.r]*(theta.UMPBT$right - theta0) - 
                   ((batch.size[n+1]*((theta.UMPBT$right^2) - (theta0^2)))/2))/(sigma^2))
          
          # for left sided check
          LR0_n.l[not.reached.decisionH0.AR.l] = 
            exp((cumsum0_n[not.reached.decisionH0.AR.l]*(theta.UMPBT$left - theta0) - 
                   ((batch.size[n+1]*((theta.UMPBT$left^2) - (theta0^2)))/2))/(sigma^2))
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]<=Accept.threshold
          RejectedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]>=Reject.threshold
          reached.decisionH0_n.AR.r = AcceptedH0.underH0_n.AR.r|RejectedH0.underH0_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.r)){
            
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[AcceptedH0.underH0_n.AR.r]] = 'A'
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[RejectedH0.underH0_n.AR.r]] = 'R'
            N0.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.r[!reached.decisionH0_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]<=Accept.threshold
          RejectedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]>=Reject.threshold
          reached.decisionH0_n.AR.l = AcceptedH0.underH0_n.AR.l|RejectedH0.underH0_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.l)){
            
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[AcceptedH0.underH0_n.AR.l]] = 'A'
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[RejectedH0.underH0_n.AR.l]] = 'R'
            N0.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH0.AR.l = not.reached.decisionH0.AR.l[!reached.decisionH0_n.AR.l]
          }
          
          not.reached.decisionH0.AR = union(not.reached.decisionH0.AR.r,
                                            not.reached.decisionH0.AR.l)
        }
        
        
        ## under right-sided H1
        if(length(not.reached.decisionH1r.AR)>0){
          
          # sum of observations at step n
          sum1r_n = rnorm(length(not.reached.decisionH1r.AR),
                          (batch.size[n+1]-batch.size[n])*theta1$right,
                          sqrt(batch.size[n+1]-batch.size[n])*sigma)
          
          # sum of observations until step n
          cumsum1r_n[not.reached.decisionH1r.AR] =
            cumsum1r_n[not.reached.decisionH1r.AR] + sum1r_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1r_n.r[not.reached.decisionH1r.AR.r] =
            exp((cumsum1r_n[not.reached.decisionH1r.AR.r]*(theta.UMPBT$right - theta0) -
                   ((batch.size[n+1]*((theta.UMPBT$right^2) - (theta0^2)))/2))/(sigma^2))
          
          # for left sided check
          LR1r_n.l[not.reached.decisionH1r.AR.l] =
            exp((cumsum1r_n[not.reached.decisionH1r.AR.l]*(theta.UMPBT$left - theta0) -
                   ((batch.size[n+1]*((theta.UMPBT$left^2) - (theta0^2)))/2))/(sigma^2))
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]<=Accept.threshold
          RejectedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]>=Reject.threshold
          reached.decisionH1r_n.AR.r = AcceptedH0.underH1r_n.AR.r|RejectedH0.underH1r_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.r)){
            
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[AcceptedH0.underH1r_n.AR.r]] = 'A'
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[RejectedH0.underH1r_n.AR.r]] = 'R'
            N1r.AR.r[not.reached.decisionH1r.AR.r[reached.decisionH1r_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.r[!reached.decisionH1r_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]<=Accept.threshold
          RejectedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]>=Reject.threshold
          reached.decisionH1r_n.AR.l = AcceptedH0.underH1r_n.AR.l|RejectedH0.underH1r_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.l)){
            
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[AcceptedH0.underH1r_n.AR.l]] = 'A'
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[RejectedH0.underH1r_n.AR.l]] = 'R'
            N1r.AR.l[not.reached.decisionH1r.AR.l[reached.decisionH1r_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH1r.AR.l = not.reached.decisionH1r.AR.l[!reached.decisionH1r_n.AR.l]
          }
          
          not.reached.decisionH1r.AR = union(not.reached.decisionH1r.AR.r,
                                             not.reached.decisionH1r.AR.l)
        }
        
        
        ## under left-sided H1
        if(length(not.reached.decisionH1l.AR)>0){
          
          # sum of observations at step n
          sum1l_n = rnorm(length(not.reached.decisionH1l.AR),
                          (batch.size[n+1]-batch.size[n])*theta1$left,
                          sqrt(batch.size[n+1]-batch.size[n])*sigma)
          
          # sum of observations until step n
          cumsum1l_n[not.reached.decisionH1l.AR] =
            cumsum1l_n[not.reached.decisionH1l.AR] + sum1l_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1l_n.r[not.reached.decisionH1l.AR.r] =
            exp((cumsum1l_n[not.reached.decisionH1l.AR.r]*(theta.UMPBT$right - theta0) -
                   ((batch.size[n+1]*((theta.UMPBT$right^2) - (theta0^2)))/2))/(sigma^2))
          
          # for left sided check
          LR1l_n.l[not.reached.decisionH1l.AR.l] =
            exp((cumsum1l_n[not.reached.decisionH1l.AR.l]*(theta.UMPBT$left - theta0) -
                   ((batch.size[n+1]*((theta.UMPBT$left^2) - (theta0^2)))/2))/(sigma^2))
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]<=Accept.threshold
          RejectedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]>=Reject.threshold
          reached.decisionH1l_n.AR.r = AcceptedH0.underH1l_n.AR.r|RejectedH0.underH1l_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.r)){
            
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[AcceptedH0.underH1l_n.AR.r]] = 'A'
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[RejectedH0.underH1l_n.AR.r]] = 'R'
            N1l.AR.r[not.reached.decisionH1l.AR.r[reached.decisionH1l_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.r[!reached.decisionH1l_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]<=Accept.threshold
          RejectedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]>=Reject.threshold
          reached.decisionH1l_n.AR.l = AcceptedH0.underH1l_n.AR.l|RejectedH0.underH1l_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.l)){
            
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[AcceptedH0.underH1l_n.AR.l]] = 'A'
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[RejectedH0.underH1l_n.AR.l]] = 'R'
            N1l.AR.l[not.reached.decisionH1l.AR.l[reached.decisionH1l_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH1l.AR.l = not.reached.decisionH1l.AR.l[!reached.decisionH1l_n.AR.l]
          }
          
          not.reached.decisionH1l.AR = union(not.reached.decisionH1l.AR.r,
                                             not.reached.decisionH1l.AR.l)
        }
        
        setTxtProgressBar(pb, n)
      }
      
      
      ### both-sided checking
      ## under H0
      # accepted or rejected ones
      accepted.by.both0 = intersect(which(decision.underH0.AR.r=='A'),
                                    which(decision.underH0.AR.l=='A'))
      onlyrejected.by.right0 = intersect(which(decision.underH0.AR.r=='R'),
                                         which(decision.underH0.AR.l!='R'))
      onlyrejected.by.left0 = intersect(which(decision.underH0.AR.r!='R'),
                                        which(decision.underH0.AR.l=='R'))
      rejected.by.both0 = intersect(which(decision.underH0.AR.r=='R'),
                                    which(decision.underH0.AR.l=='R'))
      
      # sample sizes required
      N0.AR[accepted.by.both0] = pmax(N0.AR.r[accepted.by.both0],
                                      N0.AR.l[accepted.by.both0])
      N0.AR[onlyrejected.by.right0] = N0.AR.r[onlyrejected.by.right0]
      N0.AR[onlyrejected.by.left0] = N0.AR.l[onlyrejected.by.left0]
      N0.AR[rejected.by.both0] = pmin(N0.AR.r[rejected.by.both0],
                                      N0.AR.l[rejected.by.both0])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right0 = intersect(which(decision.underH0.AR.r=='A'),
                                         which(is.na(decision.underH0.AR.l)))
      onlyaccepted.by.left0 = intersect(which(is.na(decision.underH0.AR.r)),
                                        which(decision.underH0.AR.l=='A'))
      both.inconclusive0 = intersect(which(is.na(decision.underH0.AR.r)),
                                     which(is.na(decision.underH0.AR.l)))
      all.inconclusive0 = c(onlyaccepted.by.right0, onlyaccepted.by.left0,
                            both.inconclusive0)
      nNot.reached.decisionH0.AR = length(all.inconclusive0)
      
      # Type I error probability
      type1.error.AR[c(onlyrejected.by.right0, onlyrejected.by.left0,
                       rejected.by.both0)] = T
      
      
      ## under right-sided H1
      # accepted or rejected ones
      accepted.by.both1r = intersect(which(decision.underH1r.AR.r=='A'),
                                     which(decision.underH1r.AR.l=='A'))
      onlyrejected.by.right1r = intersect(which(decision.underH1r.AR.r=='R'),
                                          which(decision.underH1r.AR.l!='R'))
      onlyrejected.by.left1r = intersect(which(decision.underH1r.AR.r!='R'),
                                         which(decision.underH1r.AR.l=='R'))
      rejected.by.both1r = intersect(which(decision.underH1r.AR.r=='R'),
                                     which(decision.underH1r.AR.l=='R'))
      
      # sample sizes required
      N1r.AR[accepted.by.both1r] = pmax(N1r.AR.r[accepted.by.both1r],
                                        N1r.AR.l[accepted.by.both1r])
      N1r.AR[onlyrejected.by.right1r] = N1r.AR.r[onlyrejected.by.right1r]
      N1r.AR[onlyrejected.by.left1r] = N1r.AR.l[onlyrejected.by.left1r]
      N1r.AR[rejected.by.both1r] = pmin(N1r.AR.r[rejected.by.both1r],
                                        N1r.AR.l[rejected.by.both1r])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1r = intersect(which(decision.underH1r.AR.r=='A'),
                                          which(is.na(decision.underH1r.AR.l)))
      onlyaccepted.by.left1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                         which(decision.underH1r.AR.l=='A'))
      both.inconclusive1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                      which(is.na(decision.underH1r.AR.l)))
      all.inconclusive1r = c(onlyaccepted.by.right1r, onlyaccepted.by.left1r,
                             both.inconclusive1r)
      nNot.reached.decisionH1r.AR = length(all.inconclusive1r)
      
      # Type I error probability
      PowerH1r.AR[c(onlyrejected.by.right1r, onlyrejected.by.left1r,
                    rejected.by.both1r)] = T
      
      
      ## under left-sided H1
      # accepted or rejected ones
      accepted.by.both1l = intersect(which(decision.underH1l.AR.r=='A'),
                                     which(decision.underH1l.AR.l=='A'))
      onlyrejected.by.right1l = intersect(which(decision.underH1l.AR.r=='R'),
                                          which(decision.underH1l.AR.l!='R'))
      onlyrejected.by.left1l = intersect(which(decision.underH1l.AR.r!='R'),
                                         which(decision.underH1l.AR.l=='R'))
      rejected.by.both1l = intersect(which(decision.underH1l.AR.r=='R'),
                                     which(decision.underH1l.AR.l=='R'))
      
      # sample sizes required
      N1l.AR[accepted.by.both1l] = pmax(N1l.AR.r[accepted.by.both1l],
                                        N1l.AR.l[accepted.by.both1l])
      N1l.AR[onlyrejected.by.right1l] = N1l.AR.r[onlyrejected.by.right1l]
      N1l.AR[onlyrejected.by.left1l] = N1l.AR.l[onlyrejected.by.left1l]
      N1l.AR[rejected.by.both1l] = pmin(N1l.AR.r[rejected.by.both1l],
                                        N1l.AR.l[rejected.by.both1l])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1l = intersect(which(decision.underH1l.AR.r=='A'),
                                          which(is.na(decision.underH1l.AR.l)))
      onlyaccepted.by.left1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                         which(decision.underH1l.AR.l=='A'))
      both.inconclusive1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                      which(is.na(decision.underH1l.AR.l)))
      all.inconclusive1l = c(onlyaccepted.by.right1l, onlyaccepted.by.left1l,
                             both.inconclusive1l)
      nNot.reached.decisionH1l.AR = length(all.inconclusive1l)
      
      # Type I error probability
      PowerH1l.AR[c(onlyrejected.by.right1l, onlyrejected.by.left1l,
                    rejected.by.both1l)] = T
      
      
      ## determining termination threshold
      ## H0 is rejected if LR or (BF) is >= termination threshold
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        term.thresh.possible.choices =
          c(LR0_n.r[onlyaccepted.by.left0],
            LR0_n.l[onlyaccepted.by.right0],
            pmin(LR0_n.r[both.inconclusive0], LR0_n.l[both.inconclusive0]))
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          max.LR0_n = max(term.thresh.possible.choices)
          nDecimal.accuracy = ceiling(-log10(min(0.01, Reject.threshold - max.LR0_n)))
          termination.threshold.AR = (floor(max.LR0_n*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01, min(term.thresh.possible.choices) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(term.thresh.possible.choices))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(term.thresh.possible.choices))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(cumRejFreq_not.reached.decisionH0.AR[1]>nNewRejects.AR){
            
            nDecimal.accuracy =
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR =
              (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                       (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR +
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      ## attained Type II error probability
      # right-sided H1
      actual.PowerH1r.AR.r = mean(PowerH1r.AR) +
        sum(c(LR1r_n.r[onlyaccepted.by.left1r],
              LR1r_n.l[onlyaccepted.by.right1r],
              pmax(LR1r_n.r[both.inconclusive1r], LR1r_n.l[both.inconclusive1r]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1r.AR = 1 - actual.PowerH1r.AR.r
      
      # left-sided H1
      actual.PowerH1l.AR.r = mean(PowerH1l.AR) +
        sum(c(LR1l_n.r[onlyaccepted.by.left1l],
              LR1l_n.l[onlyaccepted.by.right1l],
              pmax(LR1l_n.r[both.inconclusive1l], LR1l_n.l[both.inconclusive1l]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1l.AR = 1 - actual.PowerH1l.AR.r
      
      ## Expected sample sizes
      EN0 = mean(N0.AR)     # under H0
      EN1r = mean(N1r.AR)   # under right-sided H1
      EN1l = mean(N1l.AR)   # under left-sided H1
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", round(termination.threshold.AR, 3)))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print("Attained Type II error probability:")
        print(paste(" On the right: ", round(actual.type2.errorH1r.AR, 4)))
        print(paste(" On the left: ", round(actual.type2.errorH1l.AR, 4)))
        print("Expected sample size at the alternatives:")
        print(paste(" On the right: ", round(EN1r, 2)))
        print(paste(" On the left: ", round(EN1l, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  "Type2.attained" = c(actual.type2.errorH1r.AR, actual.type2.errorH1l.AR),
                  'N' = list('H0' = N0.AR, 'right' = N1r.AR, 'left' = N1l.AR),
                  'EN' = c(EN0, EN1r, EN1l), "theta.UMPBT" = theta.UMPBT,
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneZ', 'side' = side, 'theta0' = theta0, 'sigma' = sigma,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
      
    }else{
      
      ################ comparison at user provided point alternative ################
      
      ################ UMPBT alternative ################
      theta.UMPBT = list('right' = UMPBT.alt(test.type = 'oneZ', side = 'right', 
                                             theta0 = theta0, N = N.max, 
                                             Type1 = Type1.target/2, sigma = sigma),
                         'left' = UMPBT.alt(test.type = 'oneZ', side = 'left', 
                                            theta0 = theta0, N = N.max,
                                            Type1 = Type1.target/2, sigma = sigma))
      
      # msg
      if(verbose==T){
        
        print("Alternative under comparison:")
        print(paste(' On the right: ', round(theta1$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta1$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(theta.UMPBT$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta.UMPBT$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # required storages
      cumsum0_n = cumsum1r_n = cumsum1l_n =
        LR0_n.r = LR0_n.l = LR1r_n.r = LR1r_n.l = LR1l_n.r = LR1l_n.l = numeric(nReplicate)
      type1.error.AR = PowerH1r.AR = PowerH1l.AR = rep(F, nReplicate)
      N0.AR = N0.AR.r = N0.AR.l = N1r.AR = N1r.AR.r = N1r.AR.l = 
        N1l.AR = N1l.AR.r = N1l.AR.l = rep(N.max, nReplicate)
      decision.underH0.AR.r = decision.underH0.AR.l = 
        decision.underH1r.AR.r = decision.underH1r.AR.l = 
        decision.underH1l.AR.r = decision.underH1l.AR.l = rep(NA, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.l =
        not.reached.decisionH1r.AR = not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.l =
        not.reached.decisionH1l.AR = not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.l =
        1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # sum of observations at step n
          sum0_n = rnorm(length(not.reached.decisionH0.AR),
                         (batch.size[n+1]-batch.size[n])*theta0,
                         sqrt(batch.size[n+1]-batch.size[n])*sigma)
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + sum0_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR0_n.r[not.reached.decisionH0.AR.r] = 
            exp((cumsum0_n[not.reached.decisionH0.AR.r]*(theta.UMPBT$right - theta0) - 
                   ((batch.size[n+1]*((theta.UMPBT$right^2) - (theta0^2)))/2))/(sigma^2))
          
          # for left sided check
          LR0_n.l[not.reached.decisionH0.AR.l] = 
            exp((cumsum0_n[not.reached.decisionH0.AR.l]*(theta.UMPBT$left - theta0) - 
                   ((batch.size[n+1]*((theta.UMPBT$left^2) - (theta0^2)))/2))/(sigma^2))
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]<=Accept.threshold
          RejectedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]>=Reject.threshold
          reached.decisionH0_n.AR.r = AcceptedH0.underH0_n.AR.r|RejectedH0.underH0_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.r)){
            
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[AcceptedH0.underH0_n.AR.r]] = 'A'
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[RejectedH0.underH0_n.AR.r]] = 'R'
            N0.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.r[!reached.decisionH0_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]<=Accept.threshold
          RejectedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]>=Reject.threshold
          reached.decisionH0_n.AR.l = AcceptedH0.underH0_n.AR.l|RejectedH0.underH0_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.l)){
            
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[AcceptedH0.underH0_n.AR.l]] = 'A'
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[RejectedH0.underH0_n.AR.l]] = 'R'
            N0.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH0.AR.l = not.reached.decisionH0.AR.l[!reached.decisionH0_n.AR.l]
          }
          
          not.reached.decisionH0.AR = union(not.reached.decisionH0.AR.r,
                                            not.reached.decisionH0.AR.l)
        }
        
        
        ## under right-sided H1
        if(length(not.reached.decisionH1r.AR)>0){
          
          # sum of observations at step n
          sum1r_n = rnorm(length(not.reached.decisionH1r.AR),
                          (batch.size[n+1]-batch.size[n])*theta1$right,
                          sqrt(batch.size[n+1]-batch.size[n])*sigma)
          
          # sum of observations until step n
          cumsum1r_n[not.reached.decisionH1r.AR] =
            cumsum1r_n[not.reached.decisionH1r.AR] + sum1r_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1r_n.r[not.reached.decisionH1r.AR.r] =
            exp((cumsum1r_n[not.reached.decisionH1r.AR.r]*(theta.UMPBT$right - theta0) -
                   ((batch.size[n+1]*((theta.UMPBT$right^2) - (theta0^2)))/2))/(sigma^2))
          
          # for left sided check
          LR1r_n.l[not.reached.decisionH1r.AR.l] =
            exp((cumsum1r_n[not.reached.decisionH1r.AR.l]*(theta.UMPBT$left - theta0) -
                   ((batch.size[n+1]*((theta.UMPBT$left^2) - (theta0^2)))/2))/(sigma^2))
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]<=Accept.threshold
          RejectedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]>=Reject.threshold
          reached.decisionH1r_n.AR.r = AcceptedH0.underH1r_n.AR.r|RejectedH0.underH1r_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.r)){
            
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[AcceptedH0.underH1r_n.AR.r]] = 'A'
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[RejectedH0.underH1r_n.AR.r]] = 'R'
            N1r.AR.r[not.reached.decisionH1r.AR.r[reached.decisionH1r_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.r[!reached.decisionH1r_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]<=Accept.threshold
          RejectedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]>=Reject.threshold
          reached.decisionH1r_n.AR.l = AcceptedH0.underH1r_n.AR.l|RejectedH0.underH1r_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.l)){
            
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[AcceptedH0.underH1r_n.AR.l]] = 'A'
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[RejectedH0.underH1r_n.AR.l]] = 'R'
            N1r.AR.l[not.reached.decisionH1r.AR.l[reached.decisionH1r_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH1r.AR.l = not.reached.decisionH1r.AR.l[!reached.decisionH1r_n.AR.l]
          }
          
          not.reached.decisionH1r.AR = union(not.reached.decisionH1r.AR.r,
                                             not.reached.decisionH1r.AR.l)
        }
        
        
        ## under left-sided H1
        if(length(not.reached.decisionH1l.AR)>0){
          
          # sum of observations at step n
          sum1l_n = rnorm(length(not.reached.decisionH1l.AR),
                          (batch.size[n+1]-batch.size[n])*theta1$left,
                          sqrt(batch.size[n+1]-batch.size[n])*sigma)
          
          # sum of observations until step n
          cumsum1l_n[not.reached.decisionH1l.AR] =
            cumsum1l_n[not.reached.decisionH1l.AR] + sum1l_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1l_n.r[not.reached.decisionH1l.AR.r] =
            exp((cumsum1l_n[not.reached.decisionH1l.AR.r]*(theta.UMPBT$right - theta0) -
                   ((batch.size[n+1]*((theta.UMPBT$right^2) - (theta0^2)))/2))/(sigma^2))
          
          # for left sided check
          LR1l_n.l[not.reached.decisionH1l.AR.l] =
            exp((cumsum1l_n[not.reached.decisionH1l.AR.l]*(theta.UMPBT$left - theta0) -
                   ((batch.size[n+1]*((theta.UMPBT$left^2) - (theta0^2)))/2))/(sigma^2))
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]<=Accept.threshold
          RejectedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]>=Reject.threshold
          reached.decisionH1l_n.AR.r = AcceptedH0.underH1l_n.AR.r|RejectedH0.underH1l_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.r)){
            
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[AcceptedH0.underH1l_n.AR.r]] = 'A'
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[RejectedH0.underH1l_n.AR.r]] = 'R'
            N1l.AR.r[not.reached.decisionH1l.AR.r[reached.decisionH1l_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.r[!reached.decisionH1l_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]<=Accept.threshold
          RejectedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]>=Reject.threshold
          reached.decisionH1l_n.AR.l = AcceptedH0.underH1l_n.AR.l|RejectedH0.underH1l_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.l)){
            
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[AcceptedH0.underH1l_n.AR.l]] = 'A'
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[RejectedH0.underH1l_n.AR.l]] = 'R'
            N1l.AR.l[not.reached.decisionH1l.AR.l[reached.decisionH1l_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH1l.AR.l = not.reached.decisionH1l.AR.l[!reached.decisionH1l_n.AR.l]
          }
          
          not.reached.decisionH1l.AR = union(not.reached.decisionH1l.AR.r,
                                             not.reached.decisionH1l.AR.l)
        }
        
        setTxtProgressBar(pb, n)
      }
      
      
      ### both-sided checking
      ## under H0
      # accepted or rejected ones
      accepted.by.both0 = intersect(which(decision.underH0.AR.r=='A'),
                                    which(decision.underH0.AR.l=='A'))
      onlyrejected.by.right0 = intersect(which(decision.underH0.AR.r=='R'),
                                         which(decision.underH0.AR.l!='R'))
      onlyrejected.by.left0 = intersect(which(decision.underH0.AR.r!='R'),
                                        which(decision.underH0.AR.l=='R'))
      rejected.by.both0 = intersect(which(decision.underH0.AR.r=='R'),
                                    which(decision.underH0.AR.l=='R'))
      
      # sample sizes required
      N0.AR[accepted.by.both0] = pmax(N0.AR.r[accepted.by.both0],
                                      N0.AR.l[accepted.by.both0])
      N0.AR[onlyrejected.by.right0] = N0.AR.r[onlyrejected.by.right0]
      N0.AR[onlyrejected.by.left0] = N0.AR.l[onlyrejected.by.left0]
      N0.AR[rejected.by.both0] = pmin(N0.AR.r[rejected.by.both0],
                                      N0.AR.l[rejected.by.both0])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right0 = intersect(which(decision.underH0.AR.r=='A'),
                                         which(is.na(decision.underH0.AR.l)))
      onlyaccepted.by.left0 = intersect(which(is.na(decision.underH0.AR.r)),
                                        which(decision.underH0.AR.l=='A'))
      both.inconclusive0 = intersect(which(is.na(decision.underH0.AR.r)),
                                     which(is.na(decision.underH0.AR.l)))
      all.inconclusive0 = c(onlyaccepted.by.right0, onlyaccepted.by.left0,
                            both.inconclusive0)
      nNot.reached.decisionH0.AR = length(all.inconclusive0)
      
      # Type I error probability
      type1.error.AR[c(onlyrejected.by.right0, onlyrejected.by.left0,
                       rejected.by.both0)] = T
      
      
      ## under right-sided H1
      # accepted or rejected ones
      accepted.by.both1r = intersect(which(decision.underH1r.AR.r=='A'),
                                     which(decision.underH1r.AR.l=='A'))
      onlyrejected.by.right1r = intersect(which(decision.underH1r.AR.r=='R'),
                                          which(decision.underH1r.AR.l!='R'))
      onlyrejected.by.left1r = intersect(which(decision.underH1r.AR.r!='R'),
                                         which(decision.underH1r.AR.l=='R'))
      rejected.by.both1r = intersect(which(decision.underH1r.AR.r=='R'),
                                     which(decision.underH1r.AR.l=='R'))
      
      # sample sizes required
      N1r.AR[accepted.by.both1r] = pmax(N1r.AR.r[accepted.by.both1r],
                                        N1r.AR.l[accepted.by.both1r])
      N1r.AR[onlyrejected.by.right1r] = N1r.AR.r[onlyrejected.by.right1r]
      N1r.AR[onlyrejected.by.left1r] = N1r.AR.l[onlyrejected.by.left1r]
      N1r.AR[rejected.by.both1r] = pmin(N1r.AR.r[rejected.by.both1r],
                                        N1r.AR.l[rejected.by.both1r])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1r = intersect(which(decision.underH1r.AR.r=='A'),
                                          which(is.na(decision.underH1r.AR.l)))
      onlyaccepted.by.left1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                         which(decision.underH1r.AR.l=='A'))
      both.inconclusive1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                      which(is.na(decision.underH1r.AR.l)))
      all.inconclusive1r = c(onlyaccepted.by.right1r, onlyaccepted.by.left1r,
                             both.inconclusive1r)
      nNot.reached.decisionH1r.AR = length(all.inconclusive1r)
      
      # Type I error probability
      PowerH1r.AR[c(onlyrejected.by.right1r, onlyrejected.by.left1r,
                    rejected.by.both1r)] = T
      
      
      ## under left-sided H1
      # accepted or rejected ones
      accepted.by.both1l = intersect(which(decision.underH1l.AR.r=='A'),
                                     which(decision.underH1l.AR.l=='A'))
      onlyrejected.by.right1l = intersect(which(decision.underH1l.AR.r=='R'),
                                          which(decision.underH1l.AR.l!='R'))
      onlyrejected.by.left1l = intersect(which(decision.underH1l.AR.r!='R'),
                                         which(decision.underH1l.AR.l=='R'))
      rejected.by.both1l = intersect(which(decision.underH1l.AR.r=='R'),
                                     which(decision.underH1l.AR.l=='R'))
      
      # sample sizes required
      N1l.AR[accepted.by.both1l] = pmax(N1l.AR.r[accepted.by.both1l],
                                        N1l.AR.l[accepted.by.both1l])
      N1l.AR[onlyrejected.by.right1l] = N1l.AR.r[onlyrejected.by.right1l]
      N1l.AR[onlyrejected.by.left1l] = N1l.AR.l[onlyrejected.by.left1l]
      N1l.AR[rejected.by.both1l] = pmin(N1l.AR.r[rejected.by.both1l],
                                        N1l.AR.l[rejected.by.both1l])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1l = intersect(which(decision.underH1l.AR.r=='A'),
                                          which(is.na(decision.underH1l.AR.l)))
      onlyaccepted.by.left1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                         which(decision.underH1l.AR.l=='A'))
      both.inconclusive1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                      which(is.na(decision.underH1l.AR.l)))
      all.inconclusive1l = c(onlyaccepted.by.right1l, onlyaccepted.by.left1l,
                             both.inconclusive1l)
      nNot.reached.decisionH1l.AR = length(all.inconclusive1l)
      
      # Type I error probability
      PowerH1l.AR[c(onlyrejected.by.right1l, onlyrejected.by.left1l,
                    rejected.by.both1l)] = T
      
      
      ## determining termination threshold
      ## H0 is rejected if LR or (BF) is >= termination threshold
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        term.thresh.possible.choices =
          c(LR0_n.r[onlyaccepted.by.left0],
            LR0_n.l[onlyaccepted.by.right0],
            pmin(LR0_n.r[both.inconclusive0], LR0_n.l[both.inconclusive0]))
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          max.LR0_n = max(term.thresh.possible.choices)
          nDecimal.accuracy = ceiling(-log10(min(0.01, Reject.threshold - max.LR0_n)))
          termination.threshold.AR = (floor(max.LR0_n*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01, min(term.thresh.possible.choices) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(term.thresh.possible.choices))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(term.thresh.possible.choices))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(cumRejFreq_not.reached.decisionH0.AR[1]>nNewRejects.AR){
            
            nDecimal.accuracy =
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR =
              (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                       (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR +
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      ## attained Type II error probability
      # right-sided H1
      actual.PowerH1r.AR.r = mean(PowerH1r.AR) +
        sum(c(LR1r_n.r[onlyaccepted.by.left1r],
              LR1r_n.l[onlyaccepted.by.right1r],
              pmax(LR1r_n.r[both.inconclusive1r], LR1r_n.l[both.inconclusive1r]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1r.AR = 1 - actual.PowerH1r.AR.r
      
      # left-sided H1
      actual.PowerH1l.AR.r = mean(PowerH1l.AR) +
        sum(c(LR1l_n.r[onlyaccepted.by.left1l],
              LR1l_n.l[onlyaccepted.by.right1l],
              pmax(LR1l_n.r[both.inconclusive1l], LR1l_n.l[both.inconclusive1l]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1l.AR = 1 - actual.PowerH1l.AR.r
      
      ## Expected sample sizes
      EN0 = mean(N0.AR)     # under H0
      EN1r = mean(N1r.AR)   # under right-sided H1
      EN1l = mean(N1l.AR)   # under left-sided H1
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", round(termination.threshold.AR, 3)))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print("Attained Type II error probability:")
        print(paste(" On the right: ", round(actual.type2.errorH1r.AR, 4)))
        print(paste(" On the left: ", round(actual.type2.errorH1l.AR, 4)))
        print("Expected sample size at the alternatives:")
        print(paste(" On the right: ", round(EN1r, 2)))
        print(paste(" On the left: ", round(EN1l, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  "Type2.attained" = c(actual.type2.errorH1r.AR, actual.type2.errorH1l.AR),
                  'N' = list('H0' = N0.AR, 'right' = N1r.AR, 'left' = N1l.AR),
                  'EN' = c(EN0, EN1r, EN1l), "theta.UMPBT" = theta.UMPBT,
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneZ', 'side' = side, 'theta0' = theta0, 'sigma' = sigma,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
    }
  }
}

#### one-sample t test ####
design.MSPRT_oneT = function(side = 'right', theta0 = 0, theta1 = T,
                             Type1.target =.005, Type2.target = .2,
                             N.max, batch.size,
                             nReplicate = 1e+6, verbose = T, seed = 1){
  
  if(side!='both'){
    
    ################################# one-sample t (right/left sided) #################################
    
    ## batch sizes and N.max
    if(missing(batch.size)){
      
      if(missing(N.max)){
        
        return("Either 'batch.size' or 'N.max' needs to be specified")
        
      }else{batch.size = c(2, rep(1, N.max-2))}
      
    }else{
      
      if(batch.size[1]<2){
        
        return("First batch size should be at least 2")
        
      }else{
        
        if(missing(N.max)){
          
          N.max = sum(batch.size)
          
        }else{
          
          if(sum(batch.size)!=N.max) return("Sum of batch.size should add up to N.max")
        }
      }
    }
    
    nAnalyses = length(batch.size)
    
    # msg
    if(verbose){
      
      if((batch.size[1]>2)||any(batch.size[-1]>1)){
        
        cat('\n')
        print("=========================================================================")
        print("Designing the group sequential MSPRT for a one-sample t test:")
        print("=========================================================================")
        
      }else{
        
        cat('\n')
        print("=========================================================================")
        print("Designing the sequential MSPRT for a one-sample t test:")
        print("=========================================================================")
      }
      
      print(paste("Maximum available sample size: ", N.max, sep = ""))
      print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
      print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
      print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
      print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
      print(paste("Hypothesized value under H0: ", theta0, sep = ""))
      print(paste("Direction of the H1: ", side, sep = ""))
    }
    
    batch.size = c(0, cumsum(batch.size))
    
    
    if(is.logical(theta1)&&(theta1==F)){
      
      ################ no fixed-design alternative ################
      
      # msg
      if(verbose==T){
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # cut-off (with sign) in fixed design one-sample t test
      signed_t.alpha = (2*(side=='right')-1)*qt(Type1.target, df = N.max -1, lower.tail = F)
      
      # required storages
      cumSS0_n = cumsum0_n = LR0_n = numeric(nReplicate)
      type1.error.AR = rep(F, nReplicate)
      N0.AR = rep(N.max, nReplicate)
      not.reached.decisionH0.AR = 1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # observations at step n
          if(length(not.reached.decisionH0.AR)>1){
            
            obs0_n = mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                            FUN = function(X){
                              
                              rnorm(length(not.reached.decisionH0.AR), theta0, 1)
                            })
            
          }else{
            
            obs0_n = matrix(mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                                   FUN = function(X){
                                     
                                     rnorm(length(not.reached.decisionH0.AR), theta0, 1)
                                     
                                   }), nrow = 1, ncol = batch.size[n+1]-batch.size[n],
                            byrow = T)
          }
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + rowSums(obs0_n)
          
          # sum of squares of observations until step n
          cumSS0_n[not.reached.decisionH0.AR] = 
            cumSS0_n[not.reached.decisionH0.AR] + rowSums(obs0_n^2)
          
          # xbar and (n-1)*(s^2) until step n
          xbar0_n = cumsum0_n[not.reached.decisionH0.AR]/batch.size[n+1]
          divisor.s0_n.sq = 
            cumSS0_n[not.reached.decisionH0.AR] - ((cumsum0_n[not.reached.decisionH0.AR])^2)/batch.size[n+1]
          
          # likelihood ratio of observations until step n
          LR0_n[not.reached.decisionH0.AR] = 
            ((1 + (batch.size[n+1]*((xbar0_n - theta0)^2))/divisor.s0_n.sq)/
               (1 + (batch.size[n+1]*((xbar0_n - (theta0 + signed_t.alpha*
                                                    sqrt(divisor.s0_n.sq/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s0_n.sq))^(batch.size[n+1]/2)
          
          # comparing with the thresholds
          AcceptedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]<=Accept.threshold)
          RejectedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]>=Reject.threshold)
          reached.decisionH0_n.AR = union(AcceptedH0.underH0_n.AR, RejectedH0.underH0_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH0_n.AR)>0){
            
            N0.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch.size[n+1]
            type1.error.AR[not.reached.decisionH0.AR[RejectedH0.underH0_n.AR]] = T
            not.reached.decisionH0.AR = not.reached.decisionH0.AR[-reached.decisionH0_n.AR]
          }
        }
        
        setTxtProgressBar(pb, n)
      }
      
      # determining termination threshold
      # H0 is rejected if LR or (BF) is >= termination threshold
      nNot.reached.decisionH0.AR = length(not.reached.decisionH0.AR)
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 Reject.threshold -
                                                   max(LR0_n[not.reached.decisionH0.AR]))))
          termination.threshold.AR = (floor(max(LR0_n[not.reached.decisionH0.AR])*
                                              (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 min(LR0_n[not.reached.decisionH0.AR]) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(LR0_n[not.reached.decisionH0.AR]))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(LR0_n[not.reached.decisionH0.AR]))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(min(cumRejFreq_not.reached.decisionH0.AR)>nNewRejects.AR){
            
            nDecimal.accuracy = 
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR + 
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      # Expected sample sizes
      EN0 = mean(N0.AR)
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", Reject.threshold))
        print(paste("Termination threshold: ", termination.threshold.AR))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("TypeI.attained" = actual.type1.error.AR, 
                  "N" = list('H0' = N0.AR), "EN" = EN0, "Type2.fixed.design" = Type2.target, 
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneT', 'side' = side, 'theta0' = theta0,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
      
    }else if(is.logical(theta1)&&(theta1==T)){
      
      ################ comparison at the fixed-design alternative (default) ################
      theta1 = fixed_design.alt(test.type = 'oneT', side = side, theta0 = theta0,
                                N = N.max, Type1 = Type1.target, Type2 = Type2.target)
      
      # msg
      if(verbose==T){
        print(paste("Alternative under comparison: ", round(theta1, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # cut-off (with sign) in fixed design one-sample t test
      signed_t.alpha = (2*(side=='right')-1)*qt(Type1.target, df = N.max -1, lower.tail = F)
      
      # required storages
      cumSS0_n = cumSS1_n = cumsum0_n = cumsum1_n = LR0_n = LR1_n = numeric(nReplicate)
      type1.error.AR = type2.error.AR = rep(F, nReplicate)
      N0.AR = N1.AR = rep(N.max, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH1.AR = 1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # observations at step n
          if(length(not.reached.decisionH0.AR)>1){
            
            obs0_n = mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                            FUN = function(X){
                              
                              rnorm(length(not.reached.decisionH0.AR), theta0, 1)
                            })
            
          }else{
            
            obs0_n = matrix(mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                                   FUN = function(X){
                                     
                                     rnorm(length(not.reached.decisionH0.AR), theta0, 1)
                                     
                                   }), nrow = 1, ncol = batch.size[n+1]-batch.size[n],
                            byrow = T)
          }
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + rowSums(obs0_n)
          
          # sum of squares of observations until step n
          cumSS0_n[not.reached.decisionH0.AR] = 
            cumSS0_n[not.reached.decisionH0.AR] + rowSums(obs0_n^2)
          
          # xbar and (n-1)*(s^2) until step n
          xbar0_n = cumsum0_n[not.reached.decisionH0.AR]/batch.size[n+1]
          divisor.s0_n.sq = 
            cumSS0_n[not.reached.decisionH0.AR] - ((cumsum0_n[not.reached.decisionH0.AR])^2)/batch.size[n+1]
          
          # likelihood ratio of observations until step n
          LR0_n[not.reached.decisionH0.AR] = 
            ((1 + (batch.size[n+1]*((xbar0_n - theta0)^2))/divisor.s0_n.sq)/
               (1 + (batch.size[n+1]*((xbar0_n - (theta0 + signed_t.alpha*
                                                    sqrt(divisor.s0_n.sq/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s0_n.sq))^(batch.size[n+1]/2)
          
          # comparing with the thresholds
          AcceptedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]<=Accept.threshold)
          RejectedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]>=Reject.threshold)
          reached.decisionH0_n.AR = union(AcceptedH0.underH0_n.AR, RejectedH0.underH0_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH0_n.AR)>0){
            
            N0.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch.size[n+1]
            type1.error.AR[not.reached.decisionH0.AR[RejectedH0.underH0_n.AR]] = T
            not.reached.decisionH0.AR = not.reached.decisionH0.AR[-reached.decisionH0_n.AR]
          }
        }
        
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          # observations at step n
          if(length(not.reached.decisionH1.AR)>1){
            
            obs1_n = mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                            FUN = function(X){
                              
                              rnorm(length(not.reached.decisionH1.AR), theta1, 1)
                            })
            
          }else{
            
            obs1_n = matrix(mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                                   FUN = function(X){
                                     
                                     rnorm(length(not.reached.decisionH1.AR), theta1, 1)
                                     
                                   }), nrow = 1, ncol = batch.size[n+1]-batch.size[n],
                            byrow = T)
          }
          
          # sum of observations until step n
          cumsum1_n[not.reached.decisionH1.AR] =
            cumsum1_n[not.reached.decisionH1.AR] + rowSums(obs1_n)
          
          # sum of squares of observations until step n
          cumSS1_n[not.reached.decisionH1.AR] = 
            cumSS1_n[not.reached.decisionH1.AR] + rowSums(obs1_n^2)
          
          # xbar and (n-1)*(s^2) until step n
          xbar1_n = cumsum1_n[not.reached.decisionH1.AR]/batch.size[n+1]
          divisor.s1_n.sq = 
            cumSS1_n[not.reached.decisionH1.AR] - ((cumsum1_n[not.reached.decisionH1.AR])^2)/batch.size[n+1]
          
          # likelihood ratio of observations until step n
          LR1_n[not.reached.decisionH1.AR] = 
            ((1 + (batch.size[n+1]*((xbar1_n - theta0)^2))/divisor.s1_n.sq)/
               (1 + (batch.size[n+1]*((xbar1_n - (theta0 + signed_t.alpha*
                                                    sqrt(divisor.s1_n.sq/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s1_n.sq))^(batch.size[n+1]/2)
          
          # comparing with the thresholds
          AcceptedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]<=Accept.threshold)
          RejectedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]>=Reject.threshold)
          reached.decisionH1_n.AR = union(AcceptedH0.underH1_n.AR, RejectedH0.underH1_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH1_n.AR)>0){
            
            N1.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch.size[n+1]
            type2.error.AR[not.reached.decisionH1.AR[AcceptedH0.underH1_n.AR]] = T
            not.reached.decisionH1.AR = not.reached.decisionH1.AR[-reached.decisionH1_n.AR]
          }
        }
        
        setTxtProgressBar(pb, n)
      }
      
      # determining termination threshold
      # H0 is rejected if LR or (BF) is >= termination threshold
      nNot.reached.decisionH0.AR = length(not.reached.decisionH0.AR)
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 Reject.threshold -
                                                   max(LR0_n[not.reached.decisionH0.AR]))))
          termination.threshold.AR = (floor(max(LR0_n[not.reached.decisionH0.AR])*
                                              (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 min(LR0_n[not.reached.decisionH0.AR]) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(LR0_n[not.reached.decisionH0.AR]))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(LR0_n[not.reached.decisionH0.AR]))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(min(cumRejFreq_not.reached.decisionH0.AR)>nNewRejects.AR){
            
            nDecimal.accuracy = 
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR + 
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      # attained Type II error probability
      actual.type2.error.AR = mean(type2.error.AR) +
        sum(LR1_n[not.reached.decisionH1.AR]<termination.threshold.AR)/nReplicate
      
      # Expected sample sizes
      EN0 = mean(N0.AR)
      EN1 = mean(N1.AR)
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", Reject.threshold))
        print(paste("Termination threshold: ", termination.threshold.AR))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Attained Type II error probability: ", round(actual.type2.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print(paste("Expected sample size at the alternative: ", round(EN1, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("TypeI.attained" = actual.type1.error.AR, 
                  "TypeII.attained" = actual.type2.error.AR, 
                  "N" = list('H0' = N0.AR, 'H1' = N1.AR), "EN" = c(EN0, EN1),
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target, 
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneT', 'side' = side, 'theta0' = theta0,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
      
    }else{
      
      ################ comparison at user provided point alternative ################
      
      # msg
      if(verbose==T){
        print(paste("Alternative under comparison: ", round(theta1, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # cut-off (with sign) in fixed design one-sample t test
      signed_t.alpha = (2*(side=='right')-1)*qt(Type1.target, df = N.max -1, lower.tail = F)
      
      # required storages
      cumSS0_n = cumSS1_n = cumsum0_n = cumsum1_n = LR0_n = LR1_n = numeric(nReplicate)
      type1.error.AR = type2.error.AR = rep(F, nReplicate)
      N0.AR = N1.AR = rep(N.max, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH1.AR = 1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # observations at step n
          if(length(not.reached.decisionH0.AR)>1){
            
            obs0_n = mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                            FUN = function(X){
                              
                              rnorm(length(not.reached.decisionH0.AR), theta0, 1)
                            })
            
          }else{
            
            obs0_n = matrix(mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                                   FUN = function(X){
                                     
                                     rnorm(length(not.reached.decisionH0.AR), theta0, 1)
                                     
                                   }), nrow = 1, ncol = batch.size[n+1]-batch.size[n],
                            byrow = T)
          }
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + rowSums(obs0_n)
          
          # sum of squares of observations until step n
          cumSS0_n[not.reached.decisionH0.AR] = 
            cumSS0_n[not.reached.decisionH0.AR] + rowSums(obs0_n^2)
          
          # xbar and (n-1)*(s^2) until step n
          xbar0_n = cumsum0_n[not.reached.decisionH0.AR]/batch.size[n+1]
          divisor.s0_n.sq = 
            cumSS0_n[not.reached.decisionH0.AR] - ((cumsum0_n[not.reached.decisionH0.AR])^2)/batch.size[n+1]
          
          # likelihood ratio of observations until step n
          LR0_n[not.reached.decisionH0.AR] = 
            ((1 + (batch.size[n+1]*((xbar0_n - theta0)^2))/divisor.s0_n.sq)/
               (1 + (batch.size[n+1]*((xbar0_n - (theta0 + signed_t.alpha*
                                                    sqrt(divisor.s0_n.sq/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s0_n.sq))^(batch.size[n+1]/2)
          
          # comparing with the thresholds
          AcceptedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]<=Accept.threshold)
          RejectedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]>=Reject.threshold)
          reached.decisionH0_n.AR = union(AcceptedH0.underH0_n.AR, RejectedH0.underH0_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH0_n.AR)>0){
            
            N0.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch.size[n+1]
            type1.error.AR[not.reached.decisionH0.AR[RejectedH0.underH0_n.AR]] = T
            not.reached.decisionH0.AR = not.reached.decisionH0.AR[-reached.decisionH0_n.AR]
          }
        }
        
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          # observations at step n
          if(length(not.reached.decisionH1.AR)>1){
            
            obs1_n = mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                            FUN = function(X){
                              
                              rnorm(length(not.reached.decisionH1.AR), theta1, 1)
                            })
            
          }else{
            
            obs1_n = matrix(mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                                   FUN = function(X){
                                     
                                     rnorm(length(not.reached.decisionH1.AR), theta1, 1)
                                     
                                   }), nrow = 1, ncol = batch.size[n+1]-batch.size[n],
                            byrow = T)
          }
          
          # sum of observations until step n
          cumsum1_n[not.reached.decisionH1.AR] =
            cumsum1_n[not.reached.decisionH1.AR] + rowSums(obs1_n)
          
          # sum of squares of observations until step n
          cumSS1_n[not.reached.decisionH1.AR] = 
            cumSS1_n[not.reached.decisionH1.AR] + rowSums(obs1_n^2)
          
          # xbar and (n-1)*(s^2) until step n
          xbar1_n = cumsum1_n[not.reached.decisionH1.AR]/batch.size[n+1]
          divisor.s1_n.sq = 
            cumSS1_n[not.reached.decisionH1.AR] - ((cumsum1_n[not.reached.decisionH1.AR])^2)/batch.size[n+1]
          
          # likelihood ratio of observations until step n
          LR1_n[not.reached.decisionH1.AR] = 
            ((1 + (batch.size[n+1]*((xbar1_n - theta0)^2))/divisor.s1_n.sq)/
               (1 + (batch.size[n+1]*((xbar1_n - (theta0 + signed_t.alpha*
                                                    sqrt(divisor.s1_n.sq/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s1_n.sq))^(batch.size[n+1]/2)
          
          # comparing with the thresholds
          AcceptedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]<=Accept.threshold)
          RejectedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]>=Reject.threshold)
          reached.decisionH1_n.AR = union(AcceptedH0.underH1_n.AR, RejectedH0.underH1_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH1_n.AR)>0){
            
            N1.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch.size[n+1]
            type2.error.AR[not.reached.decisionH1.AR[AcceptedH0.underH1_n.AR]] = T
            not.reached.decisionH1.AR = not.reached.decisionH1.AR[-reached.decisionH1_n.AR]
          }
        }
        
        setTxtProgressBar(pb, n)
      }
      
      # determining termination threshold
      # H0 is rejected if LR or (BF) is >= termination threshold
      nNot.reached.decisionH0.AR = length(not.reached.decisionH0.AR)
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 Reject.threshold -
                                                   max(LR0_n[not.reached.decisionH0.AR]))))
          termination.threshold.AR = (floor(max(LR0_n[not.reached.decisionH0.AR])*
                                              (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 min(LR0_n[not.reached.decisionH0.AR]) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(LR0_n[not.reached.decisionH0.AR]))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(LR0_n[not.reached.decisionH0.AR]))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(min(cumRejFreq_not.reached.decisionH0.AR)>nNewRejects.AR){
            
            nDecimal.accuracy = 
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR + 
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      # attained Type II error probability
      actual.type2.error.AR = mean(type2.error.AR) +
        sum(LR1_n[not.reached.decisionH1.AR]<termination.threshold.AR)/nReplicate
      
      # Expected sample sizes
      EN0 = mean(N0.AR)
      EN1 = mean(N1.AR)
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", Reject.threshold))
        print(paste("Termination threshold: ", termination.threshold.AR))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Attained Type II error probability: ", round(actual.type2.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print(paste("Expected sample size at the alternative: ", round(EN1, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("TypeI.attained" = actual.type1.error.AR, 
                  "TypeII.attained" = actual.type2.error.AR, 
                  "N" = list('H0' = N0.AR, 'H1' = N1.AR), "EN" = c(EN0, EN1),
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target, 
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneT', 'side' = side, 'theta0' = theta0,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
    }
    
  }else{
    
    ################################# one-sample t (both sided) #################################
    
    ## batch sizes and N.max
    if(missing(batch.size)){
      
      if(missing(N.max)){
        
        return("Either 'batch.size' or 'N.max' needs to be specified")
        
      }else{batch.size = c(2, rep(1, N.max-2))}
      
    }else{
      
      if(batch.size[1]<2){
        
        return("First batch size should be at least 2")
        
      }else{
        
        if(missing(N.max)){
          
          N.max = sum(batch.size)
          
        }else{
          
          if(sum(batch.size)!=N.max) return("Sum of batch.size should add up to N.max")
        }
      }
    }
    
    nAnalyses = length(batch.size)
    
    # msg
    if(verbose){
      
      if((batch.size[1]>2)||any(batch.size[-1]>1)){
        
        cat('\n')
        print("=========================================================================")
        print("Designing the group sequential MSPRT for a one-sample t test:")
        print("=========================================================================")
        
      }else{
        
        cat('\n')
        print("=========================================================================")
        print("Designing the sequential MSPRT for a one-sample t test:")
        print("=========================================================================")
      }
      
      print(paste("Maximum available sample size: ", N.max, sep = ""))
      print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
      print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
      print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
      print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
      print(paste("Hypothesized value under H0: ", theta0, sep = ""))
      print(paste("Direction of the H1: ", side, sep = ""))
    }
    
    batch.size = c(0, cumsum(batch.size))
    
    
    if(is.logical(theta1)&&(theta1==F)){
      
      ################ no fixed-design alternative ################
      
      # msg
      if(verbose==T){
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # cut-off (with sign) in fixed design one-sample t test
      t.alpha = qt(Type1.target/2, df = N.max -1, lower.tail = F)
      
      # required storages
      cumSS0_n = cumsum0_n = LR0_n.r = LR0_n.l = numeric(nReplicate)
      type1.error.AR = rep(F, nReplicate)
      N0.AR = N0.AR.r = N0.AR.l = rep(N.max, nReplicate)
      decision.underH0.AR.r = decision.underH0.AR.l = rep(NA, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.l =
        1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # observations at step n
          if(length(not.reached.decisionH0.AR)>1){
            
            obs0_n = mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                            FUN = function(X){
                              
                              rnorm(length(not.reached.decisionH0.AR), theta0, 1)
                            })
            
          }else{
            
            obs0_n = matrix(mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                                   FUN = function(X){
                                     
                                     rnorm(length(not.reached.decisionH0.AR), theta0, 1)
                                     
                                   }), nrow = 1, ncol = batch.size[n+1]-batch.size[n],
                            byrow = T)
          }
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + rowSums(obs0_n)
          
          # sum of squares of observations until step n
          cumSS0_n[not.reached.decisionH0.AR] = 
            cumSS0_n[not.reached.decisionH0.AR] + rowSums(obs0_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          # for right sided check
          xbar0_n.r = cumsum0_n[not.reached.decisionH0.AR.r]/batch.size[n+1]
          divisor.s0_n.sq.r = 
            cumSS0_n[not.reached.decisionH0.AR.r] - 
            ((cumsum0_n[not.reached.decisionH0.AR.r])^2)/batch.size[n+1]
          
          # for left sided check
          xbar0_n.l = cumsum0_n[not.reached.decisionH0.AR.l]/batch.size[n+1]
          divisor.s0_n.sq.l = 
            cumSS0_n[not.reached.decisionH0.AR.l] - 
            ((cumsum0_n[not.reached.decisionH0.AR.l])^2)/batch.size[n+1]
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR0_n.r[not.reached.decisionH0.AR.r] = 
            ((1 + (batch.size[n+1]*((xbar0_n.r - theta0)^2))/divisor.s0_n.sq.r)/
               (1 + (batch.size[n+1]*((xbar0_n.r - 
                                         (theta0 + t.alpha*
                                            sqrt(divisor.s0_n.sq.r/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s0_n.sq.r))^(batch.size[n+1]/2)
          
          # for left sided check
          LR0_n.l[not.reached.decisionH0.AR.l] = 
            ((1 + (batch.size[n+1]*((xbar0_n.l - theta0)^2))/divisor.s0_n.sq.l)/
               (1 + (batch.size[n+1]*((xbar0_n.l - 
                                         (theta0 - t.alpha*
                                            sqrt(divisor.s0_n.sq.l/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s0_n.sq.l))^(batch.size[n+1]/2)
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]<=Accept.threshold
          RejectedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]>=Reject.threshold
          reached.decisionH0_n.AR.r = AcceptedH0.underH0_n.AR.r|RejectedH0.underH0_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.r)){
            
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[AcceptedH0.underH0_n.AR.r]] = 'A'
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[RejectedH0.underH0_n.AR.r]] = 'R'
            N0.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.r[!reached.decisionH0_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]<=Accept.threshold
          RejectedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]>=Reject.threshold
          reached.decisionH0_n.AR.l = AcceptedH0.underH0_n.AR.l|RejectedH0.underH0_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.l)){
            
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[AcceptedH0.underH0_n.AR.l]] = 'A'
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[RejectedH0.underH0_n.AR.l]] = 'R'
            N0.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH0.AR.l = not.reached.decisionH0.AR.l[!reached.decisionH0_n.AR.l]
          }
          
          not.reached.decisionH0.AR = union(not.reached.decisionH0.AR.r,
                                            not.reached.decisionH0.AR.l)
        }
        
        setTxtProgressBar(pb, n)
      }
      
      
      ### both-sided checking
      ## under H0
      # accepted or rejected ones
      accepted.by.both0 = intersect(which(decision.underH0.AR.r=='A'),
                                    which(decision.underH0.AR.l=='A'))
      onlyrejected.by.right0 = intersect(which(decision.underH0.AR.r=='R'),
                                         which(decision.underH0.AR.l!='R'))
      onlyrejected.by.left0 = intersect(which(decision.underH0.AR.r!='R'),
                                        which(decision.underH0.AR.l=='R'))
      rejected.by.both0 = intersect(which(decision.underH0.AR.r=='R'),
                                    which(decision.underH0.AR.l=='R'))
      
      # sample sizes required
      N0.AR[accepted.by.both0] = pmax(N0.AR.r[accepted.by.both0],
                                      N0.AR.l[accepted.by.both0])
      N0.AR[onlyrejected.by.right0] = N0.AR.r[onlyrejected.by.right0]
      N0.AR[onlyrejected.by.left0] = N0.AR.l[onlyrejected.by.left0]
      N0.AR[rejected.by.both0] = pmin(N0.AR.r[rejected.by.both0],
                                      N0.AR.l[rejected.by.both0])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right0 = intersect(which(decision.underH0.AR.r=='A'),
                                         which(is.na(decision.underH0.AR.l)))
      onlyaccepted.by.left0 = intersect(which(is.na(decision.underH0.AR.r)),
                                        which(decision.underH0.AR.l=='A'))
      both.inconclusive0 = intersect(which(is.na(decision.underH0.AR.r)),
                                     which(is.na(decision.underH0.AR.l)))
      all.inconclusive0 = c(onlyaccepted.by.right0, onlyaccepted.by.left0,
                            both.inconclusive0)
      nNot.reached.decisionH0.AR = length(all.inconclusive0)
      
      # Type I error probability
      type1.error.AR[c(onlyrejected.by.right0, onlyrejected.by.left0,
                       rejected.by.both0)] = T
      
      
      ## determining termination threshold
      ## H0 is rejected if LR or (BF) is >= termination threshold
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        term.thresh.possible.choices =
          c(LR0_n.r[onlyaccepted.by.left0],
            LR0_n.l[onlyaccepted.by.right0],
            pmin(LR0_n.r[both.inconclusive0], LR0_n.l[both.inconclusive0]))
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          max.LR0_n = max(term.thresh.possible.choices)
          nDecimal.accuracy = ceiling(-log10(min(0.01, Reject.threshold - max.LR0_n)))
          termination.threshold.AR = (floor(max.LR0_n*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01, min(term.thresh.possible.choices) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(term.thresh.possible.choices))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(term.thresh.possible.choices))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(cumRejFreq_not.reached.decisionH0.AR[1]>nNewRejects.AR){
            
            nDecimal.accuracy =
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR =
              (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                       (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR +
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      ## Expected sample sizes
      EN0 = mean(N0.AR)     # under H0
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", round(termination.threshold.AR, 3)))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  'N' = list('H0' = N0.AR), 'EN' = EN0, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneT', 'side' = side, 'theta0' = theta0, 'sigma' = sigma,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
      
    }else if(is.logical(theta1)&&(theta1==T)){
      
      ################ fixed-design alternative ################
      theta1 = list('right' = fixed_design.alt(test.type = 'oneT', side = 'right', 
                                               theta0 = theta0, N = N.max, 
                                               Type1 = Type1.target/2, Type2 = Type2.target),
                    'left' = fixed_design.alt(test.type = 'oneT', side = 'left', 
                                              theta0 = theta0, N = N.max, 
                                              Type1 = Type1.target/2, Type2 = Type2.target))
      
      # msg
      if(verbose==T){
        print("Alternative under comparison:")
        print(paste(' On the right: ', round(theta1$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta1$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # cut-off (with sign) in fixed design one-sample t test
      t.alpha = qt(Type1.target/2, df = N.max -1, lower.tail = F)
      
      # required storages
      cumSS0_n = cumSS1r_n = cumSS1l_n = 
        cumsum0_n = cumsum1r_n = cumsum1l_n = 
        LR0_n.r = LR0_n.l = LR1r_n.r = LR1r_n.l = LR1l_n.r = LR1l_n.l = numeric(nReplicate)
      type1.error.AR = PowerH1r.AR = PowerH1l.AR = rep(F, nReplicate)
      N0.AR = N0.AR.r = N0.AR.l = N1r.AR = N1r.AR.r = N1r.AR.l = 
        N1l.AR = N1l.AR.r = N1l.AR.l = rep(N.max, nReplicate)
      decision.underH0.AR.r = decision.underH0.AR.l = 
        decision.underH1r.AR.r = decision.underH1r.AR.l = 
        decision.underH1l.AR.r = decision.underH1l.AR.l = rep(NA, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.l =
        not.reached.decisionH1r.AR = not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.l =
        not.reached.decisionH1l.AR = not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.l =
        1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # observations at step n
          if(length(not.reached.decisionH0.AR)>1){
            
            obs0_n = mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                            FUN = function(X){
                              
                              rnorm(length(not.reached.decisionH0.AR), theta0, 1)
                            })
            
          }else{
            
            obs0_n = matrix(mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                                   FUN = function(X){
                                     
                                     rnorm(length(not.reached.decisionH0.AR), theta0, 1)
                                     
                                   }), nrow = 1, ncol = batch.size[n+1]-batch.size[n],
                            byrow = T)
          }
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + rowSums(obs0_n)
          
          # sum of squares of observations until step n
          cumSS0_n[not.reached.decisionH0.AR] = 
            cumSS0_n[not.reached.decisionH0.AR] + rowSums(obs0_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          # for right sided check
          xbar0_n.r = cumsum0_n[not.reached.decisionH0.AR.r]/batch.size[n+1]
          divisor.s0_n.sq.r = 
            cumSS0_n[not.reached.decisionH0.AR.r] - 
            ((cumsum0_n[not.reached.decisionH0.AR.r])^2)/batch.size[n+1]
          
          # for left sided check
          xbar0_n.l = cumsum0_n[not.reached.decisionH0.AR.l]/batch.size[n+1]
          divisor.s0_n.sq.l = 
            cumSS0_n[not.reached.decisionH0.AR.l] - 
            ((cumsum0_n[not.reached.decisionH0.AR.l])^2)/batch.size[n+1]
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR0_n.r[not.reached.decisionH0.AR.r] = 
            ((1 + (batch.size[n+1]*((xbar0_n.r - theta0)^2))/divisor.s0_n.sq.r)/
               (1 + (batch.size[n+1]*((xbar0_n.r - 
                                         (theta0 + t.alpha*
                                            sqrt(divisor.s0_n.sq.r/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s0_n.sq.r))^(batch.size[n+1]/2)
          
          # for left sided check
          LR0_n.l[not.reached.decisionH0.AR.l] = 
            ((1 + (batch.size[n+1]*((xbar0_n.l - theta0)^2))/divisor.s0_n.sq.l)/
               (1 + (batch.size[n+1]*((xbar0_n.l - 
                                         (theta0 - t.alpha*
                                            sqrt(divisor.s0_n.sq.l/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s0_n.sq.l))^(batch.size[n+1]/2)
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]<=Accept.threshold
          RejectedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]>=Reject.threshold
          reached.decisionH0_n.AR.r = AcceptedH0.underH0_n.AR.r|RejectedH0.underH0_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.r)){
            
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[AcceptedH0.underH0_n.AR.r]] = 'A'
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[RejectedH0.underH0_n.AR.r]] = 'R'
            N0.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.r[!reached.decisionH0_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]<=Accept.threshold
          RejectedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]>=Reject.threshold
          reached.decisionH0_n.AR.l = AcceptedH0.underH0_n.AR.l|RejectedH0.underH0_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.l)){
            
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[AcceptedH0.underH0_n.AR.l]] = 'A'
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[RejectedH0.underH0_n.AR.l]] = 'R'
            N0.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH0.AR.l = not.reached.decisionH0.AR.l[!reached.decisionH0_n.AR.l]
          }
          
          not.reached.decisionH0.AR = union(not.reached.decisionH0.AR.r,
                                            not.reached.decisionH0.AR.l)
        }
        
        
        ## under right-sided H1
        if(length(not.reached.decisionH1r.AR)>0){
          
          # observations at step n
          if(length(not.reached.decisionH1r.AR)>1){
            
            obs1r_n = mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                             FUN = function(X){
                               
                               rnorm(length(not.reached.decisionH1r.AR), theta1$right, 1)
                             })
            
          }else{
            
            obs1r_n = matrix(mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                                    FUN = function(X){
                                      
                                      rnorm(length(not.reached.decisionH1r.AR), theta1$right, 1)
                                      
                                    }), nrow = 1, ncol = batch.size[n+1]-batch.size[n],
                             byrow = T)
          }
          
          # sum of observations until step n
          cumsum1r_n[not.reached.decisionH1r.AR] = 
            cumsum1r_n[not.reached.decisionH1r.AR] + rowSums(obs1r_n)
          
          # sum of squares of observations until step n
          cumSS1r_n[not.reached.decisionH1r.AR] = 
            cumSS1r_n[not.reached.decisionH1r.AR] + rowSums(obs1r_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          # for right sided check
          xbar1r_n.r = cumsum1r_n[not.reached.decisionH1r.AR.r]/batch.size[n+1]
          divisor.s1r_n.sq.r = 
            cumSS1r_n[not.reached.decisionH1r.AR.r] - 
            ((cumsum1r_n[not.reached.decisionH1r.AR.r])^2)/batch.size[n+1]
          
          # for left sided check
          xbar1r_n.l = cumsum1r_n[not.reached.decisionH1r.AR.l]/batch.size[n+1]
          divisor.s1r_n.sq.l = 
            cumSS1r_n[not.reached.decisionH1r.AR.l] - 
            ((cumsum1r_n[not.reached.decisionH1r.AR.l])^2)/batch.size[n+1]
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1r_n.r[not.reached.decisionH1r.AR.r] = 
            ((1 + (batch.size[n+1]*((xbar1r_n.r - theta0)^2))/divisor.s1r_n.sq.r)/
               (1 + (batch.size[n+1]*((xbar1r_n.r - 
                                         (theta0 + t.alpha*
                                            sqrt(divisor.s1r_n.sq.r/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s1r_n.sq.r))^(batch.size[n+1]/2)
          
          # for left sided check
          LR1r_n.l[not.reached.decisionH1r.AR.l] = 
            ((1 + (batch.size[n+1]*((xbar1r_n.l - theta0)^2))/divisor.s1r_n.sq.l)/
               (1 + (batch.size[n+1]*((xbar1r_n.l - 
                                         (theta0 - t.alpha*
                                            sqrt(divisor.s1r_n.sq.l/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s1r_n.sq.l))^(batch.size[n+1]/2)
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]<=Accept.threshold
          RejectedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]>=Reject.threshold
          reached.decisionH1r_n.AR.r = AcceptedH0.underH1r_n.AR.r|RejectedH0.underH1r_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.r)){
            
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[AcceptedH0.underH1r_n.AR.r]] = 'A'
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[RejectedH0.underH1r_n.AR.r]] = 'R'
            N1r.AR.r[not.reached.decisionH1r.AR.r[reached.decisionH1r_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.r[!reached.decisionH1r_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]<=Accept.threshold
          RejectedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]>=Reject.threshold
          reached.decisionH1r_n.AR.l = AcceptedH0.underH1r_n.AR.l|RejectedH0.underH1r_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.l)){
            
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[AcceptedH0.underH1r_n.AR.l]] = 'A'
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[RejectedH0.underH1r_n.AR.l]] = 'R'
            N1r.AR.l[not.reached.decisionH1r.AR.l[reached.decisionH1r_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH1r.AR.l = not.reached.decisionH1r.AR.l[!reached.decisionH1r_n.AR.l]
          }
          
          not.reached.decisionH1r.AR = union(not.reached.decisionH1r.AR.r,
                                             not.reached.decisionH1r.AR.l)
        }
        
        
        ## under left-sided H1
        if(length(not.reached.decisionH1l.AR)>0){
          
          # observations at step n
          if(length(not.reached.decisionH1l.AR)>1){
            
            obs1l_n = mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                             FUN = function(X){
                               
                               rnorm(length(not.reached.decisionH1l.AR), theta1$left, 1)
                             })
            
          }else{
            
            obs1l_n = matrix(mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                                    FUN = function(X){
                                      
                                      rnorm(length(not.reached.decisionH1l.AR), theta1$left, 1)
                                      
                                    }), nrow = 1, ncol = batch.size[n+1]-batch.size[n],
                             byrow = T)
          }
          
          # sum of observations until step n
          cumsum1l_n[not.reached.decisionH1l.AR] = 
            cumsum1l_n[not.reached.decisionH1l.AR] + rowSums(obs1l_n)
          
          # sum of squares of observations until step n
          cumSS1l_n[not.reached.decisionH1l.AR] = 
            cumSS1l_n[not.reached.decisionH1l.AR] + rowSums(obs1l_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          # for right sided check
          xbar1l_n.r = cumsum1l_n[not.reached.decisionH1l.AR.r]/batch.size[n+1]
          divisor.s1l_n.sq.r = 
            cumSS1l_n[not.reached.decisionH1l.AR.r] - 
            ((cumsum1l_n[not.reached.decisionH1l.AR.r])^2)/batch.size[n+1]
          
          # for left sided check
          xbar1l_n.l = cumsum1l_n[not.reached.decisionH1l.AR.l]/batch.size[n+1]
          divisor.s1l_n.sq.l = 
            cumSS1l_n[not.reached.decisionH1l.AR.l] - 
            ((cumsum1l_n[not.reached.decisionH1l.AR.l])^2)/batch.size[n+1]
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1l_n.r[not.reached.decisionH1l.AR.r] = 
            ((1 + (batch.size[n+1]*((xbar1l_n.r - theta0)^2))/divisor.s1l_n.sq.r)/
               (1 + (batch.size[n+1]*((xbar1l_n.r - 
                                         (theta0 + t.alpha*
                                            sqrt(divisor.s1l_n.sq.r/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s1l_n.sq.r))^(batch.size[n+1]/2)
          
          # for left sided check
          LR1l_n.l[not.reached.decisionH1l.AR.l] = 
            ((1 + (batch.size[n+1]*((xbar1l_n.l - theta0)^2))/divisor.s1l_n.sq.l)/
               (1 + (batch.size[n+1]*((xbar1l_n.l - 
                                         (theta0 - t.alpha*
                                            sqrt(divisor.s1l_n.sq.l/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s1l_n.sq.l))^(batch.size[n+1]/2)
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]<=Accept.threshold
          RejectedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]>=Reject.threshold
          reached.decisionH1l_n.AR.r = AcceptedH0.underH1l_n.AR.r|RejectedH0.underH1l_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.r)){
            
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[AcceptedH0.underH1l_n.AR.r]] = 'A'
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[RejectedH0.underH1l_n.AR.r]] = 'R'
            N1l.AR.r[not.reached.decisionH1l.AR.r[reached.decisionH1l_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.r[!reached.decisionH1l_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]<=Accept.threshold
          RejectedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]>=Reject.threshold
          reached.decisionH1l_n.AR.l = AcceptedH0.underH1l_n.AR.l|RejectedH0.underH1l_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.l)){
            
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[AcceptedH0.underH1l_n.AR.l]] = 'A'
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[RejectedH0.underH1l_n.AR.l]] = 'R'
            N1l.AR.l[not.reached.decisionH1l.AR.l[reached.decisionH1l_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH1l.AR.l = not.reached.decisionH1l.AR.l[!reached.decisionH1l_n.AR.l]
          }
          
          not.reached.decisionH1l.AR = union(not.reached.decisionH1l.AR.r,
                                             not.reached.decisionH1l.AR.l)
        }
        
        setTxtProgressBar(pb, n)
      }
      
      
      ### both-sided checking
      ## under H0
      # accepted or rejected ones
      accepted.by.both0 = intersect(which(decision.underH0.AR.r=='A'),
                                    which(decision.underH0.AR.l=='A'))
      onlyrejected.by.right0 = intersect(which(decision.underH0.AR.r=='R'),
                                         which(decision.underH0.AR.l!='R'))
      onlyrejected.by.left0 = intersect(which(decision.underH0.AR.r!='R'),
                                        which(decision.underH0.AR.l=='R'))
      rejected.by.both0 = intersect(which(decision.underH0.AR.r=='R'),
                                    which(decision.underH0.AR.l=='R'))
      
      # sample sizes required
      N0.AR[accepted.by.both0] = pmax(N0.AR.r[accepted.by.both0],
                                      N0.AR.l[accepted.by.both0])
      N0.AR[onlyrejected.by.right0] = N0.AR.r[onlyrejected.by.right0]
      N0.AR[onlyrejected.by.left0] = N0.AR.l[onlyrejected.by.left0]
      N0.AR[rejected.by.both0] = pmin(N0.AR.r[rejected.by.both0],
                                      N0.AR.l[rejected.by.both0])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right0 = intersect(which(decision.underH0.AR.r=='A'),
                                         which(is.na(decision.underH0.AR.l)))
      onlyaccepted.by.left0 = intersect(which(is.na(decision.underH0.AR.r)),
                                        which(decision.underH0.AR.l=='A'))
      both.inconclusive0 = intersect(which(is.na(decision.underH0.AR.r)),
                                     which(is.na(decision.underH0.AR.l)))
      all.inconclusive0 = c(onlyaccepted.by.right0, onlyaccepted.by.left0,
                            both.inconclusive0)
      nNot.reached.decisionH0.AR = length(all.inconclusive0)
      
      # Type I error probability
      type1.error.AR[c(onlyrejected.by.right0, onlyrejected.by.left0,
                       rejected.by.both0)] = T
      
      
      ## under right-sided H1
      # accepted or rejected ones
      accepted.by.both1r = intersect(which(decision.underH1r.AR.r=='A'),
                                     which(decision.underH1r.AR.l=='A'))
      onlyrejected.by.right1r = intersect(which(decision.underH1r.AR.r=='R'),
                                          which(decision.underH1r.AR.l!='R'))
      onlyrejected.by.left1r = intersect(which(decision.underH1r.AR.r!='R'),
                                         which(decision.underH1r.AR.l=='R'))
      rejected.by.both1r = intersect(which(decision.underH1r.AR.r=='R'),
                                     which(decision.underH1r.AR.l=='R'))
      
      # sample sizes required
      N1r.AR[accepted.by.both1r] = pmax(N1r.AR.r[accepted.by.both1r],
                                        N1r.AR.l[accepted.by.both1r])
      N1r.AR[onlyrejected.by.right1r] = N1r.AR.r[onlyrejected.by.right1r]
      N1r.AR[onlyrejected.by.left1r] = N1r.AR.l[onlyrejected.by.left1r]
      N1r.AR[rejected.by.both1r] = pmin(N1r.AR.r[rejected.by.both1r],
                                        N1r.AR.l[rejected.by.both1r])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1r = intersect(which(decision.underH1r.AR.r=='A'),
                                          which(is.na(decision.underH1r.AR.l)))
      onlyaccepted.by.left1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                         which(decision.underH1r.AR.l=='A'))
      both.inconclusive1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                      which(is.na(decision.underH1r.AR.l)))
      all.inconclusive1r = c(onlyaccepted.by.right1r, onlyaccepted.by.left1r,
                             both.inconclusive1r)
      nNot.reached.decisionH1r.AR = length(all.inconclusive1r)
      
      # Type I error probability
      PowerH1r.AR[c(onlyrejected.by.right1r, onlyrejected.by.left1r,
                    rejected.by.both1r)] = T
      
      
      ## under left-sided H1
      # accepted or rejected ones
      accepted.by.both1l = intersect(which(decision.underH1l.AR.r=='A'),
                                     which(decision.underH1l.AR.l=='A'))
      onlyrejected.by.right1l = intersect(which(decision.underH1l.AR.r=='R'),
                                          which(decision.underH1l.AR.l!='R'))
      onlyrejected.by.left1l = intersect(which(decision.underH1l.AR.r!='R'),
                                         which(decision.underH1l.AR.l=='R'))
      rejected.by.both1l = intersect(which(decision.underH1l.AR.r=='R'),
                                     which(decision.underH1l.AR.l=='R'))
      
      # sample sizes required
      N1l.AR[accepted.by.both1l] = pmax(N1l.AR.r[accepted.by.both1l],
                                        N1l.AR.l[accepted.by.both1l])
      N1l.AR[onlyrejected.by.right1l] = N1l.AR.r[onlyrejected.by.right1l]
      N1l.AR[onlyrejected.by.left1l] = N1l.AR.l[onlyrejected.by.left1l]
      N1l.AR[rejected.by.both1l] = pmin(N1l.AR.r[rejected.by.both1l],
                                        N1l.AR.l[rejected.by.both1l])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1l = intersect(which(decision.underH1l.AR.r=='A'),
                                          which(is.na(decision.underH1l.AR.l)))
      onlyaccepted.by.left1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                         which(decision.underH1l.AR.l=='A'))
      both.inconclusive1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                      which(is.na(decision.underH1l.AR.l)))
      all.inconclusive1l = c(onlyaccepted.by.right1l, onlyaccepted.by.left1l,
                             both.inconclusive1l)
      nNot.reached.decisionH1l.AR = length(all.inconclusive1l)
      
      # Type I error probability
      PowerH1l.AR[c(onlyrejected.by.right1l, onlyrejected.by.left1l,
                    rejected.by.both1l)] = T
      
      
      ## determining termination threshold
      ## H0 is rejected if LR or (BF) is >= termination threshold
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        term.thresh.possible.choices =
          c(LR0_n.r[onlyaccepted.by.left0],
            LR0_n.l[onlyaccepted.by.right0],
            pmin(LR0_n.r[both.inconclusive0], LR0_n.l[both.inconclusive0]))
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          max.LR0_n = max(term.thresh.possible.choices)
          nDecimal.accuracy = ceiling(-log10(min(0.01, Reject.threshold - max.LR0_n)))
          termination.threshold.AR = (floor(max.LR0_n*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01, min(term.thresh.possible.choices) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(term.thresh.possible.choices))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(term.thresh.possible.choices))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(cumRejFreq_not.reached.decisionH0.AR[1]>nNewRejects.AR){
            
            nDecimal.accuracy =
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR =
              (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                       (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR +
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      ## attained Type II error probability
      # right-sided H1
      actual.PowerH1r.AR.r = mean(PowerH1r.AR) +
        sum(c(LR1r_n.r[onlyaccepted.by.left1r],
              LR1r_n.l[onlyaccepted.by.right1r],
              pmax(LR1r_n.r[both.inconclusive1r], LR1r_n.l[both.inconclusive1r]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1r.AR = 1 - actual.PowerH1r.AR.r
      
      # left-sided H1
      actual.PowerH1l.AR.r = mean(PowerH1l.AR) +
        sum(c(LR1l_n.r[onlyaccepted.by.left1l],
              LR1l_n.l[onlyaccepted.by.right1l],
              pmax(LR1l_n.r[both.inconclusive1l], LR1l_n.l[both.inconclusive1l]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1l.AR = 1 - actual.PowerH1l.AR.r
      
      ## Expected sample sizes
      EN0 = mean(N0.AR)     # under H0
      EN1r = mean(N1r.AR)   # under right-sided H1
      EN1l = mean(N1l.AR)   # under left-sided H1
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", round(termination.threshold.AR, 3)))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print("Attained Type II error probability:")
        print(paste(" On the right: ", round(actual.type2.errorH1r.AR, 4)))
        print(paste(" On the left: ", round(actual.type2.errorH1l.AR, 4)))
        print("Expected sample size at the alternatives:")
        print(paste(" On the right: ", round(EN1r, 2)))
        print(paste(" On the left: ", round(EN1l, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  "Type2.attained" = c(actual.type2.errorH1r.AR, actual.type2.errorH1l.AR),
                  'N' = list('H0' = N0.AR, 'right' = N1r.AR, 'left' = N1l.AR),
                  'EN' = c(EN0, EN1r, EN1l),
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneT', 'side' = side, 'theta0' = theta0, 'sigma' = sigma,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
      
    }else{
      
      ################ comparison at user provided point alternative ################
      
      # msg
      if(verbose==T){
        print("Alternative under comparison:")
        print(paste(' On the right: ', round(theta1$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta1$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # cut-off (with sign) in fixed design one-sample t test
      t.alpha = qt(Type1.target/2, df = N.max -1, lower.tail = F)
      
      # required storages
      cumSS0_n = cumSS1r_n = cumSS1l_n = 
        cumsum0_n = cumsum1r_n = cumsum1l_n = 
        LR0_n.r = LR0_n.l = LR1r_n.r = LR1r_n.l = LR1l_n.r = LR1l_n.l = numeric(nReplicate)
      type1.error.AR = PowerH1r.AR = PowerH1l.AR = rep(F, nReplicate)
      N0.AR = N0.AR.r = N0.AR.l = N1r.AR = N1r.AR.r = N1r.AR.l = 
        N1l.AR = N1l.AR.r = N1l.AR.l = rep(N.max, nReplicate)
      decision.underH0.AR.r = decision.underH0.AR.l = 
        decision.underH1r.AR.r = decision.underH1r.AR.l = 
        decision.underH1l.AR.r = decision.underH1l.AR.l = rep(NA, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.l =
        not.reached.decisionH1r.AR = not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.l =
        not.reached.decisionH1l.AR = not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.l =
        1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          # observations at step n
          if(length(not.reached.decisionH0.AR)>1){
            
            obs0_n = mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                            FUN = function(X){
                              
                              rnorm(length(not.reached.decisionH0.AR), theta0, 1)
                            })
            
          }else{
            
            obs0_n = matrix(mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                                   FUN = function(X){
                                     
                                     rnorm(length(not.reached.decisionH0.AR), theta0, 1)
                                     
                                   }), nrow = 1, ncol = batch.size[n+1]-batch.size[n],
                            byrow = T)
          }
          
          # sum of observations until step n
          cumsum0_n[not.reached.decisionH0.AR] = 
            cumsum0_n[not.reached.decisionH0.AR] + rowSums(obs0_n)
          
          # sum of squares of observations until step n
          cumSS0_n[not.reached.decisionH0.AR] = 
            cumSS0_n[not.reached.decisionH0.AR] + rowSums(obs0_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          # for right sided check
          xbar0_n.r = cumsum0_n[not.reached.decisionH0.AR.r]/batch.size[n+1]
          divisor.s0_n.sq.r = 
            cumSS0_n[not.reached.decisionH0.AR.r] - 
            ((cumsum0_n[not.reached.decisionH0.AR.r])^2)/batch.size[n+1]
          
          # for left sided check
          xbar0_n.l = cumsum0_n[not.reached.decisionH0.AR.l]/batch.size[n+1]
          divisor.s0_n.sq.l = 
            cumSS0_n[not.reached.decisionH0.AR.l] - 
            ((cumsum0_n[not.reached.decisionH0.AR.l])^2)/batch.size[n+1]
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR0_n.r[not.reached.decisionH0.AR.r] = 
            ((1 + (batch.size[n+1]*((xbar0_n.r - theta0)^2))/divisor.s0_n.sq.r)/
               (1 + (batch.size[n+1]*((xbar0_n.r - 
                                         (theta0 + t.alpha*
                                            sqrt(divisor.s0_n.sq.r/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s0_n.sq.r))^(batch.size[n+1]/2)
          
          # for left sided check
          LR0_n.l[not.reached.decisionH0.AR.l] = 
            ((1 + (batch.size[n+1]*((xbar0_n.l - theta0)^2))/divisor.s0_n.sq.l)/
               (1 + (batch.size[n+1]*((xbar0_n.l - 
                                         (theta0 - t.alpha*
                                            sqrt(divisor.s0_n.sq.l/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s0_n.sq.l))^(batch.size[n+1]/2)
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]<=Accept.threshold
          RejectedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]>=Reject.threshold
          reached.decisionH0_n.AR.r = AcceptedH0.underH0_n.AR.r|RejectedH0.underH0_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.r)){
            
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[AcceptedH0.underH0_n.AR.r]] = 'A'
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[RejectedH0.underH0_n.AR.r]] = 'R'
            N0.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.r[!reached.decisionH0_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]<=Accept.threshold
          RejectedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]>=Reject.threshold
          reached.decisionH0_n.AR.l = AcceptedH0.underH0_n.AR.l|RejectedH0.underH0_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.l)){
            
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[AcceptedH0.underH0_n.AR.l]] = 'A'
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[RejectedH0.underH0_n.AR.l]] = 'R'
            N0.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH0.AR.l = not.reached.decisionH0.AR.l[!reached.decisionH0_n.AR.l]
          }
          
          not.reached.decisionH0.AR = union(not.reached.decisionH0.AR.r,
                                            not.reached.decisionH0.AR.l)
        }
        
        
        ## under right-sided H1
        if(length(not.reached.decisionH1r.AR)>0){
          
          # observations at step n
          if(length(not.reached.decisionH1r.AR)>1){
            
            obs1r_n = mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                             FUN = function(X){
                               
                               rnorm(length(not.reached.decisionH1r.AR), theta1$right, 1)
                             })
            
          }else{
            
            obs1r_n = matrix(mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                                    FUN = function(X){
                                      
                                      rnorm(length(not.reached.decisionH1r.AR), theta1$right, 1)
                                      
                                    }), nrow = 1, ncol = batch.size[n+1]-batch.size[n],
                             byrow = T)
          }
          
          # sum of observations until step n
          cumsum1r_n[not.reached.decisionH1r.AR] = 
            cumsum1r_n[not.reached.decisionH1r.AR] + rowSums(obs1r_n)
          
          # sum of squares of observations until step n
          cumSS1r_n[not.reached.decisionH1r.AR] = 
            cumSS1r_n[not.reached.decisionH1r.AR] + rowSums(obs1r_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          # for right sided check
          xbar1r_n.r = cumsum1r_n[not.reached.decisionH1r.AR.r]/batch.size[n+1]
          divisor.s1r_n.sq.r = 
            cumSS1r_n[not.reached.decisionH1r.AR.r] - 
            ((cumsum1r_n[not.reached.decisionH1r.AR.r])^2)/batch.size[n+1]
          
          # for left sided check
          xbar1r_n.l = cumsum1r_n[not.reached.decisionH1r.AR.l]/batch.size[n+1]
          divisor.s1r_n.sq.l = 
            cumSS1r_n[not.reached.decisionH1r.AR.l] - 
            ((cumsum1r_n[not.reached.decisionH1r.AR.l])^2)/batch.size[n+1]
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1r_n.r[not.reached.decisionH1r.AR.r] = 
            ((1 + (batch.size[n+1]*((xbar1r_n.r - theta0)^2))/divisor.s1r_n.sq.r)/
               (1 + (batch.size[n+1]*((xbar1r_n.r - 
                                         (theta0 + t.alpha*
                                            sqrt(divisor.s1r_n.sq.r/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s1r_n.sq.r))^(batch.size[n+1]/2)
          
          # for left sided check
          LR1r_n.l[not.reached.decisionH1r.AR.l] = 
            ((1 + (batch.size[n+1]*((xbar1r_n.l - theta0)^2))/divisor.s1r_n.sq.l)/
               (1 + (batch.size[n+1]*((xbar1r_n.l - 
                                         (theta0 - t.alpha*
                                            sqrt(divisor.s1r_n.sq.l/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s1r_n.sq.l))^(batch.size[n+1]/2)
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]<=Accept.threshold
          RejectedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]>=Reject.threshold
          reached.decisionH1r_n.AR.r = AcceptedH0.underH1r_n.AR.r|RejectedH0.underH1r_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.r)){
            
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[AcceptedH0.underH1r_n.AR.r]] = 'A'
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[RejectedH0.underH1r_n.AR.r]] = 'R'
            N1r.AR.r[not.reached.decisionH1r.AR.r[reached.decisionH1r_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.r[!reached.decisionH1r_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]<=Accept.threshold
          RejectedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]>=Reject.threshold
          reached.decisionH1r_n.AR.l = AcceptedH0.underH1r_n.AR.l|RejectedH0.underH1r_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.l)){
            
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[AcceptedH0.underH1r_n.AR.l]] = 'A'
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[RejectedH0.underH1r_n.AR.l]] = 'R'
            N1r.AR.l[not.reached.decisionH1r.AR.l[reached.decisionH1r_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH1r.AR.l = not.reached.decisionH1r.AR.l[!reached.decisionH1r_n.AR.l]
          }
          
          not.reached.decisionH1r.AR = union(not.reached.decisionH1r.AR.r,
                                             not.reached.decisionH1r.AR.l)
        }
        
        
        ## under left-sided H1
        if(length(not.reached.decisionH1l.AR)>0){
          
          # observations at step n
          if(length(not.reached.decisionH1l.AR)>1){
            
            obs1l_n = mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                             FUN = function(X){
                               
                               rnorm(length(not.reached.decisionH1l.AR), theta1$left, 1)
                             })
            
          }else{
            
            obs1l_n = matrix(mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                                    FUN = function(X){
                                      
                                      rnorm(length(not.reached.decisionH1l.AR), theta1$left, 1)
                                      
                                    }), nrow = 1, ncol = batch.size[n+1]-batch.size[n],
                             byrow = T)
          }
          
          # sum of observations until step n
          cumsum1l_n[not.reached.decisionH1l.AR] = 
            cumsum1l_n[not.reached.decisionH1l.AR] + rowSums(obs1l_n)
          
          # sum of squares of observations until step n
          cumSS1l_n[not.reached.decisionH1l.AR] = 
            cumSS1l_n[not.reached.decisionH1l.AR] + rowSums(obs1l_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          # for right sided check
          xbar1l_n.r = cumsum1l_n[not.reached.decisionH1l.AR.r]/batch.size[n+1]
          divisor.s1l_n.sq.r = 
            cumSS1l_n[not.reached.decisionH1l.AR.r] - 
            ((cumsum1l_n[not.reached.decisionH1l.AR.r])^2)/batch.size[n+1]
          
          # for left sided check
          xbar1l_n.l = cumsum1l_n[not.reached.decisionH1l.AR.l]/batch.size[n+1]
          divisor.s1l_n.sq.l = 
            cumSS1l_n[not.reached.decisionH1l.AR.l] - 
            ((cumsum1l_n[not.reached.decisionH1l.AR.l])^2)/batch.size[n+1]
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1l_n.r[not.reached.decisionH1l.AR.r] = 
            ((1 + (batch.size[n+1]*((xbar1l_n.r - theta0)^2))/divisor.s1l_n.sq.r)/
               (1 + (batch.size[n+1]*((xbar1l_n.r - 
                                         (theta0 + t.alpha*
                                            sqrt(divisor.s1l_n.sq.r/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s1l_n.sq.r))^(batch.size[n+1]/2)
          
          # for left sided check
          LR1l_n.l[not.reached.decisionH1l.AR.l] = 
            ((1 + (batch.size[n+1]*((xbar1l_n.l - theta0)^2))/divisor.s1l_n.sq.l)/
               (1 + (batch.size[n+1]*((xbar1l_n.l - 
                                         (theta0 - t.alpha*
                                            sqrt(divisor.s1l_n.sq.l/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s1l_n.sq.l))^(batch.size[n+1]/2)
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]<=Accept.threshold
          RejectedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]>=Reject.threshold
          reached.decisionH1l_n.AR.r = AcceptedH0.underH1l_n.AR.r|RejectedH0.underH1l_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.r)){
            
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[AcceptedH0.underH1l_n.AR.r]] = 'A'
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[RejectedH0.underH1l_n.AR.r]] = 'R'
            N1l.AR.r[not.reached.decisionH1l.AR.r[reached.decisionH1l_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.r[!reached.decisionH1l_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]<=Accept.threshold
          RejectedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]>=Reject.threshold
          reached.decisionH1l_n.AR.l = AcceptedH0.underH1l_n.AR.l|RejectedH0.underH1l_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.l)){
            
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[AcceptedH0.underH1l_n.AR.l]] = 'A'
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[RejectedH0.underH1l_n.AR.l]] = 'R'
            N1l.AR.l[not.reached.decisionH1l.AR.l[reached.decisionH1l_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH1l.AR.l = not.reached.decisionH1l.AR.l[!reached.decisionH1l_n.AR.l]
          }
          
          not.reached.decisionH1l.AR = union(not.reached.decisionH1l.AR.r,
                                             not.reached.decisionH1l.AR.l)
        }
        
        setTxtProgressBar(pb, n)
      }
      
      
      ### both-sided checking
      ## under H0
      # accepted or rejected ones
      accepted.by.both0 = intersect(which(decision.underH0.AR.r=='A'),
                                    which(decision.underH0.AR.l=='A'))
      onlyrejected.by.right0 = intersect(which(decision.underH0.AR.r=='R'),
                                         which(decision.underH0.AR.l!='R'))
      onlyrejected.by.left0 = intersect(which(decision.underH0.AR.r!='R'),
                                        which(decision.underH0.AR.l=='R'))
      rejected.by.both0 = intersect(which(decision.underH0.AR.r=='R'),
                                    which(decision.underH0.AR.l=='R'))
      
      # sample sizes required
      N0.AR[accepted.by.both0] = pmax(N0.AR.r[accepted.by.both0],
                                      N0.AR.l[accepted.by.both0])
      N0.AR[onlyrejected.by.right0] = N0.AR.r[onlyrejected.by.right0]
      N0.AR[onlyrejected.by.left0] = N0.AR.l[onlyrejected.by.left0]
      N0.AR[rejected.by.both0] = pmin(N0.AR.r[rejected.by.both0],
                                      N0.AR.l[rejected.by.both0])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right0 = intersect(which(decision.underH0.AR.r=='A'),
                                         which(is.na(decision.underH0.AR.l)))
      onlyaccepted.by.left0 = intersect(which(is.na(decision.underH0.AR.r)),
                                        which(decision.underH0.AR.l=='A'))
      both.inconclusive0 = intersect(which(is.na(decision.underH0.AR.r)),
                                     which(is.na(decision.underH0.AR.l)))
      all.inconclusive0 = c(onlyaccepted.by.right0, onlyaccepted.by.left0,
                            both.inconclusive0)
      nNot.reached.decisionH0.AR = length(all.inconclusive0)
      
      # Type I error probability
      type1.error.AR[c(onlyrejected.by.right0, onlyrejected.by.left0,
                       rejected.by.both0)] = T
      
      
      ## under right-sided H1
      # accepted or rejected ones
      accepted.by.both1r = intersect(which(decision.underH1r.AR.r=='A'),
                                     which(decision.underH1r.AR.l=='A'))
      onlyrejected.by.right1r = intersect(which(decision.underH1r.AR.r=='R'),
                                          which(decision.underH1r.AR.l!='R'))
      onlyrejected.by.left1r = intersect(which(decision.underH1r.AR.r!='R'),
                                         which(decision.underH1r.AR.l=='R'))
      rejected.by.both1r = intersect(which(decision.underH1r.AR.r=='R'),
                                     which(decision.underH1r.AR.l=='R'))
      
      # sample sizes required
      N1r.AR[accepted.by.both1r] = pmax(N1r.AR.r[accepted.by.both1r],
                                        N1r.AR.l[accepted.by.both1r])
      N1r.AR[onlyrejected.by.right1r] = N1r.AR.r[onlyrejected.by.right1r]
      N1r.AR[onlyrejected.by.left1r] = N1r.AR.l[onlyrejected.by.left1r]
      N1r.AR[rejected.by.both1r] = pmin(N1r.AR.r[rejected.by.both1r],
                                        N1r.AR.l[rejected.by.both1r])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1r = intersect(which(decision.underH1r.AR.r=='A'),
                                          which(is.na(decision.underH1r.AR.l)))
      onlyaccepted.by.left1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                         which(decision.underH1r.AR.l=='A'))
      both.inconclusive1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                      which(is.na(decision.underH1r.AR.l)))
      all.inconclusive1r = c(onlyaccepted.by.right1r, onlyaccepted.by.left1r,
                             both.inconclusive1r)
      nNot.reached.decisionH1r.AR = length(all.inconclusive1r)
      
      # Type I error probability
      PowerH1r.AR[c(onlyrejected.by.right1r, onlyrejected.by.left1r,
                    rejected.by.both1r)] = T
      
      
      ## under left-sided H1
      # accepted or rejected ones
      accepted.by.both1l = intersect(which(decision.underH1l.AR.r=='A'),
                                     which(decision.underH1l.AR.l=='A'))
      onlyrejected.by.right1l = intersect(which(decision.underH1l.AR.r=='R'),
                                          which(decision.underH1l.AR.l!='R'))
      onlyrejected.by.left1l = intersect(which(decision.underH1l.AR.r!='R'),
                                         which(decision.underH1l.AR.l=='R'))
      rejected.by.both1l = intersect(which(decision.underH1l.AR.r=='R'),
                                     which(decision.underH1l.AR.l=='R'))
      
      # sample sizes required
      N1l.AR[accepted.by.both1l] = pmax(N1l.AR.r[accepted.by.both1l],
                                        N1l.AR.l[accepted.by.both1l])
      N1l.AR[onlyrejected.by.right1l] = N1l.AR.r[onlyrejected.by.right1l]
      N1l.AR[onlyrejected.by.left1l] = N1l.AR.l[onlyrejected.by.left1l]
      N1l.AR[rejected.by.both1l] = pmin(N1l.AR.r[rejected.by.both1l],
                                        N1l.AR.l[rejected.by.both1l])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1l = intersect(which(decision.underH1l.AR.r=='A'),
                                          which(is.na(decision.underH1l.AR.l)))
      onlyaccepted.by.left1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                         which(decision.underH1l.AR.l=='A'))
      both.inconclusive1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                      which(is.na(decision.underH1l.AR.l)))
      all.inconclusive1l = c(onlyaccepted.by.right1l, onlyaccepted.by.left1l,
                             both.inconclusive1l)
      nNot.reached.decisionH1l.AR = length(all.inconclusive1l)
      
      # Type I error probability
      PowerH1l.AR[c(onlyrejected.by.right1l, onlyrejected.by.left1l,
                    rejected.by.both1l)] = T
      
      
      ## determining termination threshold
      ## H0 is rejected if LR or (BF) is >= termination threshold
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        term.thresh.possible.choices =
          c(LR0_n.r[onlyaccepted.by.left0],
            LR0_n.l[onlyaccepted.by.right0],
            pmin(LR0_n.r[both.inconclusive0], LR0_n.l[both.inconclusive0]))
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          max.LR0_n = max(term.thresh.possible.choices)
          nDecimal.accuracy = ceiling(-log10(min(0.01, Reject.threshold - max.LR0_n)))
          termination.threshold.AR = (floor(max.LR0_n*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01, min(term.thresh.possible.choices) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(term.thresh.possible.choices))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(term.thresh.possible.choices))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(cumRejFreq_not.reached.decisionH0.AR[1]>nNewRejects.AR){
            
            nDecimal.accuracy =
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR =
              (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                       (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR +
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      ## attained Type II error probability
      # right-sided H1
      actual.PowerH1r.AR.r = mean(PowerH1r.AR) +
        sum(c(LR1r_n.r[onlyaccepted.by.left1r],
              LR1r_n.l[onlyaccepted.by.right1r],
              pmax(LR1r_n.r[both.inconclusive1r], LR1r_n.l[both.inconclusive1r]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1r.AR = 1 - actual.PowerH1r.AR.r
      
      # left-sided H1
      actual.PowerH1l.AR.r = mean(PowerH1l.AR) +
        sum(c(LR1l_n.r[onlyaccepted.by.left1l],
              LR1l_n.l[onlyaccepted.by.right1l],
              pmax(LR1l_n.r[both.inconclusive1l], LR1l_n.l[both.inconclusive1l]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1l.AR = 1 - actual.PowerH1l.AR.r
      
      ## Expected sample sizes
      EN0 = mean(N0.AR)     # under H0
      EN1r = mean(N1r.AR)   # under right-sided H1
      EN1l = mean(N1l.AR)   # under left-sided H1
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", round(termination.threshold.AR, 3)))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: ", round(EN0, 2)))
        print("Attained Type II error probability:")
        print(paste(" On the right: ", round(actual.type2.errorH1r.AR, 4)))
        print(paste(" On the left: ", round(actual.type2.errorH1l.AR, 4)))
        print("Expected sample size at the alternatives:")
        print(paste(" On the right: ", round(EN1r, 2)))
        print(paste(" On the left: ", round(EN1l, 2)))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  "Type2.attained" = c(actual.type2.errorH1r.AR, actual.type2.errorH1l.AR),
                  'N' = list('H0' = N0.AR, 'right' = N1r.AR, 'left' = N1l.AR),
                  'EN' = c(EN0, EN1r, EN1l),
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'oneT', 'side' = side, 'theta0' = theta0, 'sigma' = sigma,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'N.max' = N.max, 'batch.size' = diff(batch.size), 'nAnalyses' = nAnalyses,
                  'nReplicate' = nReplicate, 'seed' = seed))
    }
  }
}

#### two-sample z test ####
design.MSPRT_twoZ = function(side = 'right', theta0 = 0, theta1 = T,
                             Type1.target =.005, Type2.target = .2,
                             N1.max, N2.max, sigma1 = 1, sigma2 = 1,
                             batch1.size, batch2.size,
                             nReplicate = 1e+6, verbose = T, seed = 1){
  
  
  if(side!='both'){
    
    ########################### two-sample z (right/left sided) ###########################
    
    ## checking if length(batch1.size) and length(batch2.size) are equal
    if((!missing(batch1.size)) && (!missing(batch2.size)) &&
       (length(batch1.size)!=length(batch2.size))) return("Lenghts of batch1.size and batch2.size should be same")
    
    ## batch sizes and N for group 1
    if(missing(batch1.size)){
      
      if(missing(N1.max)){
        
        return(print("Either 'batch1.size' or 'N1.max' needs to be specified"))
        
      }else{batch1.size = rep(1, N1.max)}
      
    }else{
      
      if(missing(N1.max)){
        
        N1.max = sum(batch1.size)
        
      }else{
        
        if(sum(batch1.size)!=N1.max) return(print("Sum of batch1.size should add up to N1.max"))
      }
    }
    
    ## batch sizes and N for group 2
    if(missing(batch2.size)){
      
      if(missing(N2.max)){
        
        return(print("Either 'batch2.size' or 'N2.max' needs to be specified"))
        
      }else{batch2.size = rep(1, N2.max)}
      
    }else{
      
      if(missing(N2.max)){
        
        N2.max = sum(batch2.size)
        
      }else{
        
        if(sum(batch2.size)!=N1.max) return(print("Sum of batch2.size should add up to N2.max"))
      }
    }
    
    nAnalyses = length(batch1.size)
    
    ## msg
    if(verbose){
      
      if(any(batch1.size>1)||any(batch2.size>1)){
        
        cat('\n')
        print("=========================================================================")
        print("Designing the group sequential MSPRT for a two-sample z test:")
        print("=========================================================================")
        
      }else{
        
        cat('\n')
        print("=========================================================================")
        print("Designing the sequential MSPRT for a two-sample z test:")
        print("=========================================================================")
      }
      
      print("Group 1:")
      print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
      print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
      print(paste(" Known standard deviation: ", sigma1, sep = ""))
      print("Group 2:")
      print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
      print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
      print(paste(" Known standard deviation: ", sigma2, sep = ""))
      print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
      print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
      print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
      print(paste("Hypothesized value under H0: ", theta0, sep = ""))
      print(paste("Direction of the H1: ", side, sep = ""))
    }
    
    batch1.size = c(0, cumsum(batch1.size))
    batch2.size = c(0, cumsum(batch2.size))
    
    
    if(is.logical(theta1)&&(theta1==F)){
      
      ################ no alternative  comparison ################
      
      ################ UMPBT alternative ################
      theta.UMPBT = UMPBT.alt(test.type = 'twoZ', side = side, theta0 = theta0,
                              N1 = N1.max, N2 = N2.max, Type1 = Type1.target,
                              sigma1 = sigma1, sigma2 = sigma2)
      
      # msg
      if(verbose==T){
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(theta.UMPBT, 3)))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # required storages
      cumsum10_n = cumsum20_n = LR0_n = numeric(nReplicate)
      type1.error.AR = rep(F, nReplicate)
      N10.AR = rep(N1.max, nReplicate)
      N20.AR = rep(N2.max, nReplicate)
      not.reached.decisionH0.AR = 1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          ## sum of observations at step n
          # Group 1
          sum10_n = rnorm(length(not.reached.decisionH0.AR),
                          (batch1.size[n+1]-batch1.size[n])*(theta0/2),
                          sqrt(batch1.size[n+1]-batch1.size[n])*sigma1)
          
          # Group 2
          sum20_n = rnorm(length(not.reached.decisionH0.AR),
                          -(batch2.size[n+1]-batch2.size[n])*(theta0/2),
                          sqrt(batch2.size[n+1]-batch2.size[n])*sigma2)
          
          ## sum of observations until step n
          # Group 1
          cumsum10_n[not.reached.decisionH0.AR] = 
            cumsum10_n[not.reached.decisionH0.AR] + sum10_n
          
          # Group 2
          cumsum20_n[not.reached.decisionH0.AR] = 
            cumsum20_n[not.reached.decisionH0.AR] + sum20_n
          
          # likelihood ratio of observations until step n
          LR0_n[not.reached.decisionH0.AR] = 
            exp(-(((theta.UMPBT^2) - (theta0^2)) - 
                    2*(theta.UMPBT - theta0)*
                    (cumsum10_n[not.reached.decisionH0.AR]/batch1.size[n+1] - 
                       cumsum20_n[not.reached.decisionH0.AR]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          # comparing with the thresholds
          AcceptedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]<=Accept.threshold)
          RejectedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]>=Reject.threshold)
          reached.decisionH0_n.AR = union(AcceptedH0.underH0_n.AR, RejectedH0.underH0_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH0_n.AR)>0){
            
            N10.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch1.size[n+1]
            N20.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch2.size[n+1]
            type1.error.AR[not.reached.decisionH0.AR[RejectedH0.underH0_n.AR]] = T
            not.reached.decisionH0.AR = not.reached.decisionH0.AR[-reached.decisionH0_n.AR]
          }
        }
        
        setTxtProgressBar(pb, n)
      }
      
      # determining termination threshold
      # H0 is rejected if LR or (BF) is >= termination threshold
      nNot.reached.decisionH0.AR = length(not.reached.decisionH0.AR)
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 Reject.threshold -
                                                   max(LR0_n[not.reached.decisionH0.AR]))))
          termination.threshold.AR = (floor(max(LR0_n[not.reached.decisionH0.AR])*
                                              (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 min(LR0_n[not.reached.decisionH0.AR]) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(LR0_n[not.reached.decisionH0.AR]))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(LR0_n[not.reached.decisionH0.AR]))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(min(cumRejFreq_not.reached.decisionH0.AR)>nNewRejects.AR){
            
            nDecimal.accuracy = 
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR + 
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      # Expected sample sizes
      EN10 = mean(N10.AR)
      EN20 = mean(N20.AR)
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", Reject.threshold))
        print(paste("Termination threshold: ", termination.threshold.AR))
        print(paste("Attained Type I error probability: ", actual.type1.error.AR))
        print(paste(" Expected sample size under H0: Group 1 - ", round(EN10, 2),
                    ", Group 2 - ", round(EN20, 2), sep = ''))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  'N' = list('H0' = list('Group1' = N10.AR, 'Group2' = N20.AR)),
                  'EN' = list('H0' = list('Group1' = EN10, 'Group2' = EN20)),
                  "theta.UMPBT" = theta.UMPBT, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'twoZ', 'side' = side, 'theta0' = theta0, 
                  'sigma1' = sigma1, 'sigma2' = sigma2, 'N1.max' = N1.max, 'N2.max' = N2.max,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'batch1.size' = diff(batch1.size), 'batch2.size' = diff(batch2.size),
                  'nAnalyses' = nAnalyses, 'nReplicate' = nReplicate, 'seed' = seed))
      
    }else if(is.logical(theta1)&&(theta1==T)){
      
      ################ comparison at the fixed-design alternative ################
      theta1 = fixed_design.alt(test.type = 'twoZ', side = side, theta0 = theta0,
                                N1 = N1.max, N2 = N2.max, Type1 = Type1.target,
                                Type2 = Type2.target, sigma1 = 1, sigma2 = 1)
      
      ################ UMPBT alternative ################
      theta.UMPBT = UMPBT.alt(test.type = 'twoZ', side = side, theta0 = theta0,
                              N1 = N1.max, N2 = N2.max, Type1 = Type1.target,
                              sigma1 = sigma1, sigma2 = sigma2)
      
      # msg
      if(verbose==T){
        
        print(paste("Alternative under comparison: ", round(theta1, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(theta.UMPBT, 3)))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # required storages
      cumsum10_n = cumsum20_n = cumsum11_n = cumsum21_n = LR0_n = LR1_n = numeric(nReplicate)
      type1.error.AR = type2.error.AR = rep(F, nReplicate)
      N10.AR = N11.AR = rep(N1.max, nReplicate)
      N20.AR = N21.AR = rep(N2.max, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH1.AR = 1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          ## sum of observations at step n
          # Group 1
          sum10_n = rnorm(length(not.reached.decisionH0.AR),
                          (batch1.size[n+1]-batch1.size[n])*(theta0/2),
                          sqrt(batch1.size[n+1]-batch1.size[n])*sigma1)
          
          # Group 2
          sum20_n = rnorm(length(not.reached.decisionH0.AR),
                          -(batch2.size[n+1]-batch2.size[n])*(theta0/2),
                          sqrt(batch2.size[n+1]-batch2.size[n])*sigma2)
          
          ## sum of observations until step n
          # Group 1
          cumsum10_n[not.reached.decisionH0.AR] = 
            cumsum10_n[not.reached.decisionH0.AR] + sum10_n
          
          # Group 2
          cumsum20_n[not.reached.decisionH0.AR] = 
            cumsum20_n[not.reached.decisionH0.AR] + sum20_n
          
          # likelihood ratio of observations until step n
          LR0_n[not.reached.decisionH0.AR] = 
            exp(-(((theta.UMPBT^2) - (theta0^2)) - 
                    2*(theta.UMPBT - theta0)*
                    (cumsum10_n[not.reached.decisionH0.AR]/batch1.size[n+1] - 
                       cumsum20_n[not.reached.decisionH0.AR]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          # comparing with the thresholds
          AcceptedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]<=Accept.threshold)
          RejectedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]>=Reject.threshold)
          reached.decisionH0_n.AR = union(AcceptedH0.underH0_n.AR, RejectedH0.underH0_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH0_n.AR)>0){
            
            N10.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch1.size[n+1]
            N20.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch2.size[n+1]
            type1.error.AR[not.reached.decisionH0.AR[RejectedH0.underH0_n.AR]] = T
            not.reached.decisionH0.AR = not.reached.decisionH0.AR[-reached.decisionH0_n.AR]
          }
        }
        
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          ## sum of observations at step n
          # Group 1
          sum11_n = rnorm(length(not.reached.decisionH1.AR),
                          (batch1.size[n+1]-batch1.size[n])*(theta1/2),
                          sqrt(batch1.size[n+1]-batch1.size[n])*sigma1)
          
          # Group 2
          sum21_n = rnorm(length(not.reached.decisionH1.AR),
                          -(batch2.size[n+1]-batch2.size[n])*(theta1/2),
                          sqrt(batch2.size[n+1]-batch2.size[n])*sigma2)
          
          ## sum of observations until step n
          # Group 1
          cumsum11_n[not.reached.decisionH1.AR] = 
            cumsum11_n[not.reached.decisionH1.AR] + sum11_n
          
          # Group 2
          cumsum21_n[not.reached.decisionH1.AR] = 
            cumsum21_n[not.reached.decisionH1.AR] + sum21_n
          
          # likelihood ratio of observations until step n
          LR1_n[not.reached.decisionH1.AR] = 
            exp(-(((theta.UMPBT^2) - (theta0^2)) - 
                    2*(theta.UMPBT - theta0)*
                    (cumsum11_n[not.reached.decisionH1.AR]/batch1.size[n+1] - 
                       cumsum21_n[not.reached.decisionH1.AR]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          # comparing with the thresholds
          AcceptedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]<=Accept.threshold)
          RejectedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]>=Reject.threshold)
          reached.decisionH1_n.AR = union(AcceptedH0.underH1_n.AR, RejectedH0.underH1_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH1_n.AR)>0){
            
            N11.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch1.size[n+1]
            N21.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch2.size[n+1]
            type2.error.AR[not.reached.decisionH1.AR[AcceptedH0.underH1_n.AR]] = T
            not.reached.decisionH1.AR = not.reached.decisionH1.AR[-reached.decisionH1_n.AR]
          }
        }
        
        setTxtProgressBar(pb, n)
      }
      
      # determining termination threshold
      # H0 is rejected if LR or (BF) is >= termination threshold
      nNot.reached.decisionH0.AR = length(not.reached.decisionH0.AR)
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 Reject.threshold -
                                                   max(LR0_n[not.reached.decisionH0.AR]))))
          termination.threshold.AR = (floor(max(LR0_n[not.reached.decisionH0.AR])*
                                              (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 min(LR0_n[not.reached.decisionH0.AR]) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(LR0_n[not.reached.decisionH0.AR]))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(LR0_n[not.reached.decisionH0.AR]))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(min(cumRejFreq_not.reached.decisionH0.AR)>nNewRejects.AR){
            
            nDecimal.accuracy = 
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR + 
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      # attained Type II error probability
      actual.type2.error.AR = mean(type2.error.AR) +
        sum(LR1_n[not.reached.decisionH1.AR]<termination.threshold.AR)/nReplicate
      
      # Expected sample sizes
      EN10 = mean(N10.AR)
      EN20 = mean(N20.AR)
      EN11 = mean(N11.AR)
      EN21 = mean(N21.AR)
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", Reject.threshold))
        print(paste("Termination threshold: ", termination.threshold.AR))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Attained Type II error probability: ", round(actual.type2.error.AR, 4)))
        print(paste(" Expected sample size under H0: Group 1 - ", round(EN10, 2),
                    ", Group 2 - ", round(EN20, 2), sep = ''))
        print(paste(" Expected sample size at the alternative: Group 1 - ", round(EN11, 2),
                    ", Group 2 - ", round(EN21, 2), sep = ''))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  "Type2.attained" = actual.type2.error.AR,
                  'N' = list('H0' = list('Group1' = N10.AR, 'Group2' = N20.AR),
                             'H1' = list('Group1' = N11.AR, 'Group2' = N21.AR)),
                  'EN' = list('H0' = list('Group1' = EN10, 'Group2' = EN20),
                              'H1' = list('Group1' = EN11, 'Group2' = EN21)),
                  "theta.UMPBT" = theta.UMPBT,
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'twoZ', 'side' = side, 'theta0' = theta0, 
                  'sigma1' = sigma1, 'sigma2' = sigma2, 'N1.max' = N1.max, 'N2.max' = N2.max,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'batch1.size' = diff(batch1.size), 'batch2.size' = diff(batch2.size),
                  'nAnalyses' = nAnalyses, 'nReplicate' = nReplicate, 'seed' = seed))
      
    }else{
      
      ################ comparison at user provided point alternative ################
      
      ################ UMPBT alternative ################
      theta.UMPBT = UMPBT.alt(test.type = 'twoZ', side = side, theta0 = theta0,
                              N1 = N1.max, N2 = N2.max, Type1 = Type1.target,
                              sigma1 = sigma1, sigma2 = sigma2)
      
      # msg
      if(verbose==T){
        
        print(paste("Alternative under comparison: ", round(theta1, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(theta.UMPBT, 3)))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # required storages
      cumsum10_n = cumsum20_n = cumsum11_n = cumsum21_n = LR0_n = LR1_n = numeric(nReplicate)
      type1.error.AR = type2.error.AR = rep(F, nReplicate)
      N10.AR = N11.AR = rep(N1.max, nReplicate)
      N20.AR = N21.AR = rep(N2.max, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH1.AR = 1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          ## sum of observations at step n
          # Group 1
          sum10_n = rnorm(length(not.reached.decisionH0.AR),
                          (batch1.size[n+1]-batch1.size[n])*(theta0/2),
                          sqrt(batch1.size[n+1]-batch1.size[n])*sigma1)
          
          # Group 2
          sum20_n = rnorm(length(not.reached.decisionH0.AR),
                          -(batch2.size[n+1]-batch2.size[n])*(theta0/2),
                          sqrt(batch2.size[n+1]-batch2.size[n])*sigma2)
          
          ## sum of observations until step n
          # Group 1
          cumsum10_n[not.reached.decisionH0.AR] = 
            cumsum10_n[not.reached.decisionH0.AR] + sum10_n
          
          # Group 2
          cumsum20_n[not.reached.decisionH0.AR] = 
            cumsum20_n[not.reached.decisionH0.AR] + sum20_n
          
          # likelihood ratio of observations until step n
          LR0_n[not.reached.decisionH0.AR] = 
            exp(-(((theta.UMPBT^2) - (theta0^2)) - 
                    2*(theta.UMPBT - theta0)*
                    (cumsum10_n[not.reached.decisionH0.AR]/batch1.size[n+1] - 
                       cumsum20_n[not.reached.decisionH0.AR]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          # comparing with the thresholds
          AcceptedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]<=Accept.threshold)
          RejectedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]>=Reject.threshold)
          reached.decisionH0_n.AR = union(AcceptedH0.underH0_n.AR, RejectedH0.underH0_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH0_n.AR)>0){
            
            N10.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch1.size[n+1]
            N20.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch2.size[n+1]
            type1.error.AR[not.reached.decisionH0.AR[RejectedH0.underH0_n.AR]] = T
            not.reached.decisionH0.AR = not.reached.decisionH0.AR[-reached.decisionH0_n.AR]
          }
        }
        
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          ## sum of observations at step n
          # Group 1
          sum11_n = rnorm(length(not.reached.decisionH1.AR),
                          (batch1.size[n+1]-batch1.size[n])*(theta1/2),
                          sqrt(batch1.size[n+1]-batch1.size[n])*sigma1)
          
          # Group 2
          sum21_n = rnorm(length(not.reached.decisionH1.AR),
                          -(batch2.size[n+1]-batch2.size[n])*(theta1/2),
                          sqrt(batch2.size[n+1]-batch2.size[n])*sigma2)
          
          ## sum of observations until step n
          # Group 1
          cumsum11_n[not.reached.decisionH1.AR] = 
            cumsum11_n[not.reached.decisionH1.AR] + sum11_n
          
          # Group 2
          cumsum21_n[not.reached.decisionH1.AR] = 
            cumsum21_n[not.reached.decisionH1.AR] + sum21_n
          
          # likelihood ratio of observations until step n
          LR1_n[not.reached.decisionH1.AR] = 
            exp(-(((theta.UMPBT^2) - (theta0^2)) - 
                    2*(theta.UMPBT - theta0)*
                    (cumsum11_n[not.reached.decisionH1.AR]/batch1.size[n+1] - 
                       cumsum21_n[not.reached.decisionH1.AR]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          # comparing with the thresholds
          AcceptedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]<=Accept.threshold)
          RejectedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]>=Reject.threshold)
          reached.decisionH1_n.AR = union(AcceptedH0.underH1_n.AR, RejectedH0.underH1_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH1_n.AR)>0){
            
            N11.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch1.size[n+1]
            N21.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch2.size[n+1]
            type2.error.AR[not.reached.decisionH1.AR[AcceptedH0.underH1_n.AR]] = T
            not.reached.decisionH1.AR = not.reached.decisionH1.AR[-reached.decisionH1_n.AR]
          }
        }
        
        setTxtProgressBar(pb, n)
      }
      
      # determining termination threshold
      # H0 is rejected if LR or (BF) is >= termination threshold
      nNot.reached.decisionH0.AR = length(not.reached.decisionH0.AR)
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 Reject.threshold -
                                                   max(LR0_n[not.reached.decisionH0.AR]))))
          termination.threshold.AR = (floor(max(LR0_n[not.reached.decisionH0.AR])*
                                              (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 min(LR0_n[not.reached.decisionH0.AR]) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(LR0_n[not.reached.decisionH0.AR]))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(LR0_n[not.reached.decisionH0.AR]))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(min(cumRejFreq_not.reached.decisionH0.AR)>nNewRejects.AR){
            
            nDecimal.accuracy = 
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR + 
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      # attained Type II error probability
      actual.type2.error.AR = mean(type2.error.AR) +
        sum(LR1_n[not.reached.decisionH1.AR]<termination.threshold.AR)/nReplicate
      
      # Expected sample sizes
      EN10 = mean(N10.AR)
      EN20 = mean(N20.AR)
      EN11 = mean(N11.AR)
      EN21 = mean(N21.AR)
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", Reject.threshold))
        print(paste("Termination threshold: ", termination.threshold.AR))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Attained Type II error probability: ", round(actual.type2.error.AR, 4)))
        print(paste(" Expected sample size under H0: Group 1 - ", round(EN10, 2),
                    ", Group 2 - ", round(EN20, 2), sep = ''))
        print(paste(" Expected sample size at the alternative: Group 1 - ", round(EN11, 2),
                    ", Group 2 - ", round(EN21, 2), sep = ''))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  "Type2.attained" = actual.type2.error.AR,
                  'N' = list('H0' = list('Group1' = N10.AR, 'Group2' = N20.AR),
                             'H1' = list('Group1' = N11.AR, 'Group2' = N21.AR)),
                  'EN' = list('H0' = list('Group1' = EN10, 'Group2' = EN20),
                              'H1' = list('Group1' = EN11, 'Group2' = EN21)),
                  "theta.UMPBT" = theta.UMPBT,
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'twoZ', 'side' = side, 'theta0' = theta0, 
                  'sigma1' = sigma1, 'sigma2' = sigma2, 'N1.max' = N1.max, 'N2.max' = N2.max,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'batch1.size' = diff(batch1.size), 'batch2.size' = diff(batch2.size),
                  'nAnalyses' = nAnalyses, 'nReplicate' = nReplicate, 'seed' = seed))
    }
    
  }else{
    
    ################################# two-sample z (both sided) #################################
    
    ## checking if length(batch1.size) and length(batch2.size) are equal
    if((!missing(batch1.size)) && (!missing(batch2.size)) &&
       (length(batch1.size)!=length(batch2.size))) return("Lenghts of batch1.size and batch2.size should be same")
    
    ## batch sizes and N for group 1
    if(missing(batch1.size)){
      
      if(missing(N1.max)){
        
        return(print("Either 'batch1.size' or 'N1.max' needs to be specified"))
        
      }else{batch1.size = rep(1, N1.max)}
      
    }else{
      
      if(missing(N1.max)){
        
        N1.max = sum(batch1.size)
        
      }else{
        
        if(sum(batch1.size)!=N1.max) return(print("Sum of batch1.size should add up to N1.max"))
      }
    }
    
    ## batch sizes and N for group 2
    if(missing(batch2.size)){
      
      if(missing(N2.max)){
        
        return(print("Either 'batch2.size' or 'N2.max' needs to be specified"))
        
      }else{batch2.size = rep(1, N2.max)}
      
    }else{
      
      if(missing(N2.max)){
        
        N2.max = sum(batch2.size)
        
      }else{
        
        if(sum(batch2.size)!=N1.max) return(print("Sum of batch2.size should add up to N2.max"))
      }
    }
    
    nAnalyses = length(batch1.size)
    
    ## msg
    if(verbose){
      
      if(any(batch1.size>1)||any(batch2.size>1)){
        
        cat('\n')
        print("=========================================================================")
        print("Designing the group sequential MSPRT for a two-sample z test:")
        print("=========================================================================")
        
      }else{
        
        cat('\n')
        print("=========================================================================")
        print("Designing the sequential MSPRT for a two-sample z test:")
        print("=========================================================================")
      }
      
      print("Group 1:")
      print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
      print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
      print(paste(" Known standard deviation: ", sigma1, sep = ""))
      print("Group 2:")
      print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
      print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
      print(paste(" Known standard deviation: ", sigma2, sep = ""))
      print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
      print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
      print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
      print(paste("Hypothesized value under H0: ", theta0, sep = ""))
      print(paste("Direction of the H1: ", side, sep = ""))
    }
    
    batch1.size = c(0, cumsum(batch1.size))
    batch2.size = c(0, cumsum(batch2.size))
    
    
    if(is.logical(theta1)&&(theta1==F)){
      
      ################ no alternative comparison ################
      
      ################ UMPBT alternative ################
      theta.UMPBT = list('right' = UMPBT.alt(test.type = 'twoZ', side = 'right',
                                             theta0 = theta0, N1 = N1.max, N2 = N2.max, 
                                             Type1 = Type1.target/2,
                                             sigma1 = sigma1, sigma2 = sigma2),
                         'left' = UMPBT.alt(test.type = 'twoZ', side = 'left',
                                            theta0 = theta0, N1 = N1.max, N2 = N2.max, 
                                            Type1 = Type1.target/2,
                                            sigma1 = sigma1, sigma2 = sigma2))
      
      # msg
      if(verbose==T){
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(theta.UMPBT$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta.UMPBT$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # required storages
      cumsum10_n = cumsum20_n = LR0_n.r = LR0_n.l = numeric(nReplicate)
      type1.error.AR = rep(F, nReplicate)
      N10.AR = N10.AR.r = N10.AR.l = rep(N1.max, nReplicate)
      N20.AR = N20.AR.r = N20.AR.l = rep(N2.max, nReplicate)
      decision.underH0.AR.r = decision.underH0.AR.l = rep(NA, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.l =
        1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          ## sum of observations at step n
          # Group 1
          sum10_n = rnorm(length(not.reached.decisionH0.AR),
                          (batch1.size[n+1]-batch1.size[n])*(theta0/2),
                          sqrt(batch1.size[n+1]-batch1.size[n])*sigma1)
          
          # Group 2
          sum20_n = rnorm(length(not.reached.decisionH0.AR),
                          -(batch2.size[n+1]-batch2.size[n])*(theta0/2),
                          sqrt(batch2.size[n+1]-batch2.size[n])*sigma2)
          
          ## sum of observations until step n
          # Group 1
          cumsum10_n[not.reached.decisionH0.AR] = 
            cumsum10_n[not.reached.decisionH0.AR] + sum10_n
          
          # Group 2
          cumsum20_n[not.reached.decisionH0.AR] = 
            cumsum20_n[not.reached.decisionH0.AR] + sum20_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR0_n.r[not.reached.decisionH0.AR.r] = 
            exp(-(((theta.UMPBT$right^2) - (theta0^2)) - 
                    2*(theta.UMPBT$right - theta0)*
                    (cumsum10_n[not.reached.decisionH0.AR.r]/batch1.size[n+1] - 
                       cumsum20_n[not.reached.decisionH0.AR.r]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          # for left sided check
          LR0_n.l[not.reached.decisionH0.AR.l] = 
            exp(-(((theta.UMPBT$left^2) - (theta0^2)) - 
                    2*(theta.UMPBT$left - theta0)*
                    (cumsum10_n[not.reached.decisionH0.AR.l]/batch1.size[n+1] - 
                       cumsum20_n[not.reached.decisionH0.AR.l]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]<=Accept.threshold
          RejectedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]>=Reject.threshold
          reached.decisionH0_n.AR.r = AcceptedH0.underH0_n.AR.r|RejectedH0.underH0_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.r)){
            
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[AcceptedH0.underH0_n.AR.r]] = 'A'
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[RejectedH0.underH0_n.AR.r]] = 'R'
            N10.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch1.size[n+1]
            N20.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch2.size[n+1]
            not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.r[!reached.decisionH0_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]<=Accept.threshold
          RejectedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]>=Reject.threshold
          reached.decisionH0_n.AR.l = AcceptedH0.underH0_n.AR.l|RejectedH0.underH0_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.l)){
            
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[AcceptedH0.underH0_n.AR.l]] = 'A'
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[RejectedH0.underH0_n.AR.l]] = 'R'
            N10.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch1.size[n+1]
            N20.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch2.size[n+1]
            not.reached.decisionH0.AR.l = not.reached.decisionH0.AR.l[!reached.decisionH0_n.AR.l]
          }
          
          not.reached.decisionH0.AR = union(not.reached.decisionH0.AR.r,
                                            not.reached.decisionH0.AR.l)
        }
        
        setTxtProgressBar(pb, n)
      }
      
      
      ### both-sided checking
      ## under H0
      # accepted or rejected ones
      accepted.by.both0 = intersect(which(decision.underH0.AR.r=='A'),
                                    which(decision.underH0.AR.l=='A'))
      onlyrejected.by.right0 = intersect(which(decision.underH0.AR.r=='R'),
                                         which(decision.underH0.AR.l!='R'))
      onlyrejected.by.left0 = intersect(which(decision.underH0.AR.r!='R'),
                                        which(decision.underH0.AR.l=='R'))
      rejected.by.both0 = intersect(which(decision.underH0.AR.r=='R'),
                                    which(decision.underH0.AR.l=='R'))
      
      ## sample sizes required
      # Group 1
      N10.AR[accepted.by.both0] = pmax(N10.AR.r[accepted.by.both0],
                                       N10.AR.l[accepted.by.both0])
      N10.AR[onlyrejected.by.right0] = N10.AR.r[onlyrejected.by.right0]
      N10.AR[onlyrejected.by.left0] = N10.AR.l[onlyrejected.by.left0]
      N10.AR[rejected.by.both0] = pmin(N10.AR.r[rejected.by.both0],
                                       N10.AR.l[rejected.by.both0])
      
      # Group 2
      N20.AR[accepted.by.both0] = pmax(N20.AR.r[accepted.by.both0],
                                       N20.AR.l[accepted.by.both0])
      N20.AR[onlyrejected.by.right0] = N20.AR.r[onlyrejected.by.right0]
      N20.AR[onlyrejected.by.left0] = N20.AR.l[onlyrejected.by.left0]
      N20.AR[rejected.by.both0] = pmin(N20.AR.r[rejected.by.both0],
                                       N20.AR.l[rejected.by.both0])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right0 = intersect(which(decision.underH0.AR.r=='A'),
                                         which(is.na(decision.underH0.AR.l)))
      onlyaccepted.by.left0 = intersect(which(is.na(decision.underH0.AR.r)),
                                        which(decision.underH0.AR.l=='A'))
      both.inconclusive0 = intersect(which(is.na(decision.underH0.AR.r)),
                                     which(is.na(decision.underH0.AR.l)))
      all.inconclusive0 = c(onlyaccepted.by.right0, onlyaccepted.by.left0,
                            both.inconclusive0)
      nNot.reached.decisionH0.AR = length(all.inconclusive0)
      
      # Type I error probability
      type1.error.AR[c(onlyrejected.by.right0, onlyrejected.by.left0,
                       rejected.by.both0)] = T
      
      
      ## determining termination threshold
      ## H0 is rejected if LR or (BF) is >= termination threshold
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        term.thresh.possible.choices =
          c(LR0_n.r[onlyaccepted.by.left0],
            LR0_n.l[onlyaccepted.by.right0],
            pmin(LR0_n.r[both.inconclusive0], LR0_n.l[both.inconclusive0]))
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          max.LR0_n = max(term.thresh.possible.choices)
          nDecimal.accuracy = ceiling(-log10(min(0.01, Reject.threshold - max.LR0_n)))
          termination.threshold.AR = (floor(max.LR0_n*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01, min(term.thresh.possible.choices) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(term.thresh.possible.choices))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(term.thresh.possible.choices))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(cumRejFreq_not.reached.decisionH0.AR[1]>nNewRejects.AR){
            
            nDecimal.accuracy =
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR =
              (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                       (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR +
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      ## Expected sample sizes
      # Group 1
      EN10 = mean(N10.AR)     # under H0
      
      # Group 2
      EN20 = mean(N20.AR)     # under H0
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", round(termination.threshold.AR, 3)))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: Group 1 - ", round(EN10, 2), 
                    ', Group 2 - ', round(EN20, 2), sep = ''))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  'N' = list('H0' = list('Group1' = N10.AR, 'Group2' = N20.AR)),
                  'EN' = list('H0' = list('Group1' = EN10, 'Group2' = EN20)),
                  "theta.UMPBT" = theta.UMPBT, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'twoZ', 'side' = side, 'theta0' = theta0, 
                  'sigma1' = sigma1, 'sigma2' = sigma2, 'N1.max' = N1.max, 'N2.max' = N2.max,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'batch1.size' = diff(batch1.size), 'batch2.size' = diff(batch2.size),
                  'nAnalyses' = nAnalyses, 'nReplicate' = nReplicate, 'seed' = seed))
      
    }else if(is.logical(theta1)&&(theta1==T)){
      
      ################ comparison at the fixed-design alternative (default) ################
      
      theta1 = list('right' = fixed_design.alt(test.type = 'twoZ', side = 'right',
                                               theta0 = theta0, N1 = N1.max, N2 = N2.max,
                                               Type1 = Type1.target/2,
                                               Type2 = Type2.target, sigma1 = sigma1, sigma2 = sigma2),
                    'left' = fixed_design.alt(test.type = 'twoZ', side = 'left',
                                              theta0 = theta0, N1 = N1.max, N2 = N2.max,
                                              Type1 = Type1.target/2,
                                              Type2 = Type2.target, sigma1 = sigma1, sigma2 = sigma2))
      
      ################ UMPBT alternative ################
      theta.UMPBT = list('right' = UMPBT.alt(test.type = 'twoZ', side = 'right',
                                             theta0 = theta0, N1 = N1.max, N2 = N2.max, 
                                             Type1 = Type1.target/2,
                                             sigma1 = sigma1, sigma2 = sigma2),
                         'left' = UMPBT.alt(test.type = 'twoZ', side = 'left',
                                            theta0 = theta0, N1 = N1.max, N2 = N2.max, 
                                            Type1 = Type1.target/2,
                                            sigma1 = sigma1, sigma2 = sigma2))
      
      # msg
      if(verbose==T){
        
        print("Alternative under comparison:")
        print(paste(' On the right: ', round(theta1$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta1$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(theta.UMPBT$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta.UMPBT$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # required storages
      cumsum10_n = cumsum20_n = cumsum11r_n = cumsum21r_n = cumsum11l_n = cumsum21l_n = 
        LR0_n.r = LR0_n.l = LR1r_n.r = LR1r_n.l = LR1l_n.r = LR1l_n.l = numeric(nReplicate)
      type1.error.AR = PowerH1r.AR = PowerH1l.AR = rep(F, nReplicate)
      N10.AR = N10.AR.r = N10.AR.l = 
        N11r.AR = N11r.AR.r = N11r.AR.l = 
        N11l.AR = N11l.AR.r = N11l.AR.l = rep(N1.max, nReplicate)
      N20.AR = N20.AR.r = N20.AR.l = 
        N21r.AR = N21r.AR.r = N21r.AR.l = 
        N21l.AR = N21l.AR.r = N21l.AR.l = rep(N2.max, nReplicate)
      decision.underH0.AR.r = decision.underH0.AR.l = 
        decision.underH1r.AR.r = decision.underH1r.AR.l = 
        decision.underH1l.AR.r = decision.underH1l.AR.l = rep(NA, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.l =
        not.reached.decisionH1r.AR = not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.l =
        not.reached.decisionH1l.AR = not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.l =
        1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          ## sum of observations at step n
          # Group 1
          sum10_n = rnorm(length(not.reached.decisionH0.AR),
                          (batch1.size[n+1]-batch1.size[n])*(theta0/2),
                          sqrt(batch1.size[n+1]-batch1.size[n])*sigma1)
          
          # Group 2
          sum20_n = rnorm(length(not.reached.decisionH0.AR),
                          -(batch2.size[n+1]-batch2.size[n])*(theta0/2),
                          sqrt(batch2.size[n+1]-batch2.size[n])*sigma2)
          
          ## sum of observations until step n
          # Group 1
          cumsum10_n[not.reached.decisionH0.AR] = 
            cumsum10_n[not.reached.decisionH0.AR] + sum10_n
          
          # Group 2
          cumsum20_n[not.reached.decisionH0.AR] = 
            cumsum20_n[not.reached.decisionH0.AR] + sum20_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR0_n.r[not.reached.decisionH0.AR.r] = 
            exp(-(((theta.UMPBT$right^2) - (theta0^2)) - 
                    2*(theta.UMPBT$right - theta0)*
                    (cumsum10_n[not.reached.decisionH0.AR.r]/batch1.size[n+1] - 
                       cumsum20_n[not.reached.decisionH0.AR.r]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          # for left sided check
          LR0_n.l[not.reached.decisionH0.AR.l] = 
            exp(-(((theta.UMPBT$left^2) - (theta0^2)) - 
                    2*(theta.UMPBT$left - theta0)*
                    (cumsum10_n[not.reached.decisionH0.AR.l]/batch1.size[n+1] - 
                       cumsum20_n[not.reached.decisionH0.AR.l]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]<=Accept.threshold
          RejectedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]>=Reject.threshold
          reached.decisionH0_n.AR.r = AcceptedH0.underH0_n.AR.r|RejectedH0.underH0_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.r)){
            
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[AcceptedH0.underH0_n.AR.r]] = 'A'
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[RejectedH0.underH0_n.AR.r]] = 'R'
            N10.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch1.size[n+1]
            N20.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch2.size[n+1]
            not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.r[!reached.decisionH0_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]<=Accept.threshold
          RejectedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]>=Reject.threshold
          reached.decisionH0_n.AR.l = AcceptedH0.underH0_n.AR.l|RejectedH0.underH0_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.l)){
            
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[AcceptedH0.underH0_n.AR.l]] = 'A'
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[RejectedH0.underH0_n.AR.l]] = 'R'
            N10.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch1.size[n+1]
            N20.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch2.size[n+1]
            not.reached.decisionH0.AR.l = not.reached.decisionH0.AR.l[!reached.decisionH0_n.AR.l]
          }
          
          not.reached.decisionH0.AR = union(not.reached.decisionH0.AR.r,
                                            not.reached.decisionH0.AR.l)
        }
        
        
        ## under right-sided H1
        if(length(not.reached.decisionH1r.AR)>0){
          
          ## sum of observations at step n
          # Group 1
          sum11r_n = rnorm(length(not.reached.decisionH1r.AR),
                           (batch1.size[n+1]-batch1.size[n])*(theta1$right/2),
                           sqrt(batch1.size[n+1]-batch1.size[n])*sigma1)
          
          # Group 2
          sum21r_n = rnorm(length(not.reached.decisionH1r.AR),
                           -(batch2.size[n+1]-batch2.size[n])*(theta1$right/2),
                           sqrt(batch2.size[n+1]-batch2.size[n])*sigma2)
          
          ## sum of observations until step n
          # Group 1
          cumsum11r_n[not.reached.decisionH1r.AR] = 
            cumsum11r_n[not.reached.decisionH1r.AR] + sum11r_n
          
          # Group 2
          cumsum21r_n[not.reached.decisionH1r.AR] = 
            cumsum21r_n[not.reached.decisionH1r.AR] + sum21r_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1r_n.r[not.reached.decisionH1r.AR.r] = 
            exp(-(((theta.UMPBT$right^2) - (theta0^2)) - 
                    2*(theta.UMPBT$right - theta0)*
                    (cumsum11r_n[not.reached.decisionH1r.AR.r]/batch1.size[n+1] - 
                       cumsum21r_n[not.reached.decisionH1r.AR.r]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          # for left sided check
          LR1r_n.l[not.reached.decisionH1r.AR.l] = 
            exp(-(((theta.UMPBT$left^2) - (theta0^2)) - 
                    2*(theta.UMPBT$left - theta0)*
                    (cumsum11r_n[not.reached.decisionH1r.AR.l]/batch1.size[n+1] - 
                       cumsum21r_n[not.reached.decisionH1r.AR.l]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]<=Accept.threshold
          RejectedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]>=Reject.threshold
          reached.decisionH1r_n.AR.r = AcceptedH0.underH1r_n.AR.r|RejectedH0.underH1r_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.r)){
            
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[AcceptedH0.underH1r_n.AR.r]] = 'A'
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[RejectedH0.underH1r_n.AR.r]] = 'R'
            N11r.AR.r[not.reached.decisionH1r.AR.r[reached.decisionH1r_n.AR.r]] = batch1.size[n+1]
            N21r.AR.r[not.reached.decisionH1r.AR.r[reached.decisionH1r_n.AR.r]] = batch2.size[n+1]
            not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.r[!reached.decisionH1r_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]<=Accept.threshold
          RejectedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]>=Reject.threshold
          reached.decisionH1r_n.AR.l = AcceptedH0.underH1r_n.AR.l|RejectedH0.underH1r_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.l)){
            
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[AcceptedH0.underH1r_n.AR.l]] = 'A'
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[RejectedH0.underH1r_n.AR.l]] = 'R'
            N11r.AR.l[not.reached.decisionH1r.AR.l[reached.decisionH1r_n.AR.l]] = batch1.size[n+1]
            N21r.AR.l[not.reached.decisionH1r.AR.l[reached.decisionH1r_n.AR.l]] = batch2.size[n+1]
            not.reached.decisionH1r.AR.l = not.reached.decisionH1r.AR.l[!reached.decisionH1r_n.AR.l]
          }
          
          not.reached.decisionH1r.AR = union(not.reached.decisionH1r.AR.r,
                                             not.reached.decisionH1r.AR.l)
        }
        
        
        ## under left-sided H1
        if(length(not.reached.decisionH1l.AR)>0){
          
          ## sum of observations at step n
          # Group 1
          sum11l_n = rnorm(length(not.reached.decisionH1l.AR),
                           (batch1.size[n+1]-batch1.size[n])*(theta1$left/2),
                           sqrt(batch1.size[n+1]-batch1.size[n])*sigma1)
          
          # Group 2
          sum21l_n = rnorm(length(not.reached.decisionH1l.AR),
                           -(batch2.size[n+1]-batch2.size[n])*(theta1$left/2),
                           sqrt(batch2.size[n+1]-batch2.size[n])*sigma2)
          
          ## sum of observations until step n
          # Group 1
          cumsum11l_n[not.reached.decisionH1l.AR] = 
            cumsum11l_n[not.reached.decisionH1l.AR] + sum11l_n
          
          # Group 2
          cumsum21l_n[not.reached.decisionH1l.AR] = 
            cumsum21l_n[not.reached.decisionH1l.AR] + sum21l_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1l_n.r[not.reached.decisionH1l.AR.r] = 
            exp(-(((theta.UMPBT$right^2) - (theta0^2)) - 
                    2*(theta.UMPBT$right - theta0)*
                    (cumsum11l_n[not.reached.decisionH1l.AR.r]/batch1.size[n+1] - 
                       cumsum21l_n[not.reached.decisionH1l.AR.r]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          # for left sided check
          LR1l_n.l[not.reached.decisionH1l.AR.l] = 
            exp(-(((theta.UMPBT$left^2) - (theta0^2)) - 
                    2*(theta.UMPBT$left - theta0)*
                    (cumsum11l_n[not.reached.decisionH1l.AR.l]/batch1.size[n+1] - 
                       cumsum21l_n[not.reached.decisionH1l.AR.l]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]<=Accept.threshold
          RejectedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]>=Reject.threshold
          reached.decisionH1l_n.AR.r = AcceptedH0.underH1l_n.AR.r|RejectedH0.underH1l_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.r)){
            
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[AcceptedH0.underH1l_n.AR.r]] = 'A'
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[RejectedH0.underH1l_n.AR.r]] = 'R'
            N11l.AR.r[not.reached.decisionH1l.AR.r[reached.decisionH1l_n.AR.r]] = batch1.size[n+1]
            N21l.AR.r[not.reached.decisionH1l.AR.r[reached.decisionH1l_n.AR.r]] = batch2.size[n+1]
            not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.r[!reached.decisionH1l_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]<=Accept.threshold
          RejectedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]>=Reject.threshold
          reached.decisionH1l_n.AR.l = AcceptedH0.underH1l_n.AR.l|RejectedH0.underH1l_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.l)){
            
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[AcceptedH0.underH1l_n.AR.l]] = 'A'
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[RejectedH0.underH1l_n.AR.l]] = 'R'
            N11l.AR.l[not.reached.decisionH1l.AR.l[reached.decisionH1l_n.AR.l]] = batch1.size[n+1]
            N21l.AR.l[not.reached.decisionH1l.AR.l[reached.decisionH1l_n.AR.l]] = batch2.size[n+1]
            not.reached.decisionH1l.AR.l = not.reached.decisionH1l.AR.l[!reached.decisionH1l_n.AR.l]
          }
          
          not.reached.decisionH1l.AR = union(not.reached.decisionH1l.AR.r,
                                             not.reached.decisionH1l.AR.l)
        }
        
        setTxtProgressBar(pb, n)
      }
      
      
      ### both-sided checking
      ## under H0
      # accepted or rejected ones
      accepted.by.both0 = intersect(which(decision.underH0.AR.r=='A'),
                                    which(decision.underH0.AR.l=='A'))
      onlyrejected.by.right0 = intersect(which(decision.underH0.AR.r=='R'),
                                         which(decision.underH0.AR.l!='R'))
      onlyrejected.by.left0 = intersect(which(decision.underH0.AR.r!='R'),
                                        which(decision.underH0.AR.l=='R'))
      rejected.by.both0 = intersect(which(decision.underH0.AR.r=='R'),
                                    which(decision.underH0.AR.l=='R'))
      
      ## sample sizes required
      # Group 1
      N10.AR[accepted.by.both0] = pmax(N10.AR.r[accepted.by.both0],
                                       N10.AR.l[accepted.by.both0])
      N10.AR[onlyrejected.by.right0] = N10.AR.r[onlyrejected.by.right0]
      N10.AR[onlyrejected.by.left0] = N10.AR.l[onlyrejected.by.left0]
      N10.AR[rejected.by.both0] = pmin(N10.AR.r[rejected.by.both0],
                                       N10.AR.l[rejected.by.both0])
      
      # Group 2
      N20.AR[accepted.by.both0] = pmax(N20.AR.r[accepted.by.both0],
                                       N20.AR.l[accepted.by.both0])
      N20.AR[onlyrejected.by.right0] = N20.AR.r[onlyrejected.by.right0]
      N20.AR[onlyrejected.by.left0] = N20.AR.l[onlyrejected.by.left0]
      N20.AR[rejected.by.both0] = pmin(N20.AR.r[rejected.by.both0],
                                       N20.AR.l[rejected.by.both0])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right0 = intersect(which(decision.underH0.AR.r=='A'),
                                         which(is.na(decision.underH0.AR.l)))
      onlyaccepted.by.left0 = intersect(which(is.na(decision.underH0.AR.r)),
                                        which(decision.underH0.AR.l=='A'))
      both.inconclusive0 = intersect(which(is.na(decision.underH0.AR.r)),
                                     which(is.na(decision.underH0.AR.l)))
      all.inconclusive0 = c(onlyaccepted.by.right0, onlyaccepted.by.left0,
                            both.inconclusive0)
      nNot.reached.decisionH0.AR = length(all.inconclusive0)
      
      # Type I error probability
      type1.error.AR[c(onlyrejected.by.right0, onlyrejected.by.left0,
                       rejected.by.both0)] = T
      
      
      ## under right-sided H1
      # accepted or rejected ones
      accepted.by.both1r = intersect(which(decision.underH1r.AR.r=='A'),
                                     which(decision.underH1r.AR.l=='A'))
      onlyrejected.by.right1r = intersect(which(decision.underH1r.AR.r=='R'),
                                          which(decision.underH1r.AR.l!='R'))
      onlyrejected.by.left1r = intersect(which(decision.underH1r.AR.r!='R'),
                                         which(decision.underH1r.AR.l=='R'))
      rejected.by.both1r = intersect(which(decision.underH1r.AR.r=='R'),
                                     which(decision.underH1r.AR.l=='R'))
      
      ## sample sizes required
      # Group 1
      N11r.AR[accepted.by.both1r] = pmax(N11r.AR.r[accepted.by.both1r],
                                         N11r.AR.l[accepted.by.both1r])
      N11r.AR[onlyrejected.by.right1r] = N11r.AR.r[onlyrejected.by.right1r]
      N11r.AR[onlyrejected.by.left1r] = N11r.AR.l[onlyrejected.by.left1r]
      N11r.AR[rejected.by.both1r] = pmin(N11r.AR.r[rejected.by.both1r],
                                         N11r.AR.l[rejected.by.both1r])
      
      # Group 2
      N21r.AR[accepted.by.both1r] = pmax(N21r.AR.r[accepted.by.both1r],
                                         N21r.AR.l[accepted.by.both1r])
      N21r.AR[onlyrejected.by.right1r] = N21r.AR.r[onlyrejected.by.right1r]
      N21r.AR[onlyrejected.by.left1r] = N21r.AR.l[onlyrejected.by.left1r]
      N21r.AR[rejected.by.both1r] = pmin(N21r.AR.r[rejected.by.both1r],
                                         N21r.AR.l[rejected.by.both1r])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1r = intersect(which(decision.underH1r.AR.r=='A'),
                                          which(is.na(decision.underH1r.AR.l)))
      onlyaccepted.by.left1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                         which(decision.underH1r.AR.l=='A'))
      both.inconclusive1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                      which(is.na(decision.underH1r.AR.l)))
      all.inconclusive1r = c(onlyaccepted.by.right1r, onlyaccepted.by.left1r,
                             both.inconclusive1r)
      nNot.reached.decisionH1r.AR = length(all.inconclusive1r)
      
      # Type I error probability
      PowerH1r.AR[c(onlyrejected.by.right1r, onlyrejected.by.left1r,
                    rejected.by.both1r)] = T
      
      
      ## under left-sided H1
      # accepted or rejected ones
      accepted.by.both1l = intersect(which(decision.underH1l.AR.r=='A'),
                                     which(decision.underH1l.AR.l=='A'))
      onlyrejected.by.right1l = intersect(which(decision.underH1l.AR.r=='R'),
                                          which(decision.underH1l.AR.l!='R'))
      onlyrejected.by.left1l = intersect(which(decision.underH1l.AR.r!='R'),
                                         which(decision.underH1l.AR.l=='R'))
      rejected.by.both1l = intersect(which(decision.underH1l.AR.r=='R'),
                                     which(decision.underH1l.AR.l=='R'))
      
      ## sample sizes required
      # Group 1
      N11l.AR[accepted.by.both1l] = pmax(N11l.AR.r[accepted.by.both1l],
                                         N11l.AR.l[accepted.by.both1l])
      N11l.AR[onlyrejected.by.right1l] = N11l.AR.r[onlyrejected.by.right1l]
      N11l.AR[onlyrejected.by.left1l] = N11l.AR.l[onlyrejected.by.left1l]
      N11l.AR[rejected.by.both1l] = pmin(N11l.AR.r[rejected.by.both1l],
                                         N11l.AR.l[rejected.by.both1l])
      
      # Group 2
      N21l.AR[accepted.by.both1l] = pmax(N21l.AR.r[accepted.by.both1l],
                                         N21l.AR.l[accepted.by.both1l])
      N21l.AR[onlyrejected.by.right1l] = N21l.AR.r[onlyrejected.by.right1l]
      N21l.AR[onlyrejected.by.left1l] = N21l.AR.l[onlyrejected.by.left1l]
      N21l.AR[rejected.by.both1l] = pmin(N21l.AR.r[rejected.by.both1l],
                                         N21l.AR.l[rejected.by.both1l])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1l = intersect(which(decision.underH1l.AR.r=='A'),
                                          which(is.na(decision.underH1l.AR.l)))
      onlyaccepted.by.left1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                         which(decision.underH1l.AR.l=='A'))
      both.inconclusive1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                      which(is.na(decision.underH1l.AR.l)))
      all.inconclusive1l = c(onlyaccepted.by.right1l, onlyaccepted.by.left1l,
                             both.inconclusive1l)
      nNot.reached.decisionH1l.AR = length(all.inconclusive1l)
      
      # Type I error probability
      PowerH1l.AR[c(onlyrejected.by.right1l, onlyrejected.by.left1l,
                    rejected.by.both1l)] = T
      
      
      ## determining termination threshold
      ## H0 is rejected if LR or (BF) is >= termination threshold
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        term.thresh.possible.choices =
          c(LR0_n.r[onlyaccepted.by.left0],
            LR0_n.l[onlyaccepted.by.right0],
            pmin(LR0_n.r[both.inconclusive0], LR0_n.l[both.inconclusive0]))
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          max.LR0_n = max(term.thresh.possible.choices)
          nDecimal.accuracy = ceiling(-log10(min(0.01, Reject.threshold - max.LR0_n)))
          termination.threshold.AR = (floor(max.LR0_n*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01, min(term.thresh.possible.choices) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(term.thresh.possible.choices))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(term.thresh.possible.choices))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(cumRejFreq_not.reached.decisionH0.AR[1]>nNewRejects.AR){
            
            nDecimal.accuracy =
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR =
              (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                       (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR +
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      ## attained Type II error probability
      # right-sided H1
      actual.PowerH1r.AR.r = mean(PowerH1r.AR) +
        sum(c(LR1r_n.r[onlyaccepted.by.left1r],
              LR1r_n.l[onlyaccepted.by.right1r],
              pmax(LR1r_n.r[both.inconclusive1r], LR1r_n.l[both.inconclusive1r]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1r.AR = 1 - actual.PowerH1r.AR.r
      
      # left-sided H1
      actual.PowerH1l.AR.r = mean(PowerH1l.AR) +
        sum(c(LR1l_n.r[onlyaccepted.by.left1l],
              LR1l_n.l[onlyaccepted.by.right1l],
              pmax(LR1l_n.r[both.inconclusive1l], LR1l_n.l[both.inconclusive1l]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1l.AR = 1 - actual.PowerH1l.AR.r
      
      ## Expected sample sizes
      # Group 1
      EN10 = mean(N10.AR)     # under H0
      EN11r = mean(N11r.AR)   # under right-sided H1
      EN11l = mean(N11l.AR)   # under left-sided H1
      
      # Group 2
      EN20 = mean(N20.AR)     # under H0
      EN21r = mean(N21r.AR)   # under right-sided H1
      EN21l = mean(N21l.AR)   # under left-sided H1
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", round(termination.threshold.AR, 3)))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: Group 1 - ", round(EN10, 2), 
                    ', Group 2 - ', round(EN20, 2), sep = ''))
        print("Attained Type II error probability:")
        print(paste(" On the right: ", round(actual.type2.errorH1r.AR, 4)))
        print(paste(" On the left: ", round(actual.type2.errorH1l.AR, 4)))
        print("Expected sample size at the alternatives:")
        print(paste(" On the right: Group 1 - ", round(EN11r, 2), 
                    ', Group 2 - ', round(EN21r, 2), sep = ''))
        print(paste(" On the left: Group 1 - ", round(EN11l, 2), 
                    ', Group 2 - ', round(EN21l, 2), sep = ''))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  "Type2.attained" = c(actual.type2.errorH1r.AR, actual.type2.errorH1l.AR),
                  'N' = list('H0' = list('Group1' = N10.AR, 'Group2' = N20.AR),
                             'right' = list('Group1' = N11r.AR, 'Group2' = N21r.AR),
                             'left' = list('Group1' = N11l.AR, 'Group2' = N21l.AR)),
                  'EN' = list('H0' = list('Group1' = EN10, 'Group2' = EN20),
                              'right' = list('Group1' = EN11r, 'Group2' = EN21r),
                              'left' = list('Group1' = EN11l, 'Group2' = EN21l)),
                  "theta.UMPBT" = theta.UMPBT,
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'twoZ', 'side' = side, 'theta0' = theta0, 
                  'sigma1' = sigma1, 'sigma2' = sigma2, 'N1.max' = N1.max, 'N2.max' = N2.max,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'batch1.size' = diff(batch1.size), 'batch2.size' = diff(batch2.size),
                  'nAnalyses' = nAnalyses, 'nReplicate' = nReplicate, 'seed' = seed))
      
    }else{
      
      ################ comparison at user specified point alternative ################
      
      ################ UMPBT alternative ################
      theta.UMPBT = list('right' = UMPBT.alt(test.type = 'twoZ', side = 'right',
                                             theta0 = theta0, N1 = N1.max, N2 = N2.max, 
                                             Type1 = Type1.target/2,
                                             sigma1 = sigma1, sigma2 = sigma2),
                         'left' = UMPBT.alt(test.type = 'twoZ', side = 'left',
                                            theta0 = theta0, N1 = N1.max, N2 = N2.max, 
                                            Type1 = Type1.target/2,
                                            sigma1 = sigma1, sigma2 = sigma2))
      
      # msg
      if(verbose==T){
        
        print("Alternative under comparison:")
        print(paste(' On the right: ', round(theta1$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta1$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(theta.UMPBT$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta.UMPBT$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # required storages
      cumsum10_n = cumsum20_n = cumsum11r_n = cumsum21r_n = cumsum11l_n = cumsum21l_n = 
        LR0_n.r = LR0_n.l = LR1r_n.r = LR1r_n.l = LR1l_n.r = LR1l_n.l = numeric(nReplicate)
      type1.error.AR = PowerH1r.AR = PowerH1l.AR = rep(F, nReplicate)
      N10.AR = N10.AR.r = N10.AR.l = 
        N11r.AR = N11r.AR.r = N11r.AR.l = 
        N11l.AR = N11l.AR.r = N11l.AR.l = rep(N1.max, nReplicate)
      N20.AR = N20.AR.r = N20.AR.l = 
        N21r.AR = N21r.AR.r = N21r.AR.l = 
        N21l.AR = N21l.AR.r = N21l.AR.l = rep(N2.max, nReplicate)
      decision.underH0.AR.r = decision.underH0.AR.l = 
        decision.underH1r.AR.r = decision.underH1r.AR.l = 
        decision.underH1l.AR.r = decision.underH1l.AR.l = rep(NA, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.l =
        not.reached.decisionH1r.AR = not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.l =
        not.reached.decisionH1l.AR = not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.l =
        1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          ## sum of observations at step n
          # Group 1
          sum10_n = rnorm(length(not.reached.decisionH0.AR),
                          (batch1.size[n+1]-batch1.size[n])*(theta0/2),
                          sqrt(batch1.size[n+1]-batch1.size[n])*sigma1)
          
          # Group 2
          sum20_n = rnorm(length(not.reached.decisionH0.AR),
                          -(batch2.size[n+1]-batch2.size[n])*(theta0/2),
                          sqrt(batch2.size[n+1]-batch2.size[n])*sigma2)
          
          ## sum of observations until step n
          # Group 1
          cumsum10_n[not.reached.decisionH0.AR] = 
            cumsum10_n[not.reached.decisionH0.AR] + sum10_n
          
          # Group 2
          cumsum20_n[not.reached.decisionH0.AR] = 
            cumsum20_n[not.reached.decisionH0.AR] + sum20_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR0_n.r[not.reached.decisionH0.AR.r] = 
            exp(-(((theta.UMPBT$right^2) - (theta0^2)) - 
                    2*(theta.UMPBT$right - theta0)*
                    (cumsum10_n[not.reached.decisionH0.AR.r]/batch1.size[n+1] - 
                       cumsum20_n[not.reached.decisionH0.AR.r]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          # for left sided check
          LR0_n.l[not.reached.decisionH0.AR.l] = 
            exp(-(((theta.UMPBT$left^2) - (theta0^2)) - 
                    2*(theta.UMPBT$left - theta0)*
                    (cumsum10_n[not.reached.decisionH0.AR.l]/batch1.size[n+1] - 
                       cumsum20_n[not.reached.decisionH0.AR.l]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]<=Accept.threshold
          RejectedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]>=Reject.threshold
          reached.decisionH0_n.AR.r = AcceptedH0.underH0_n.AR.r|RejectedH0.underH0_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.r)){
            
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[AcceptedH0.underH0_n.AR.r]] = 'A'
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[RejectedH0.underH0_n.AR.r]] = 'R'
            N10.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch1.size[n+1]
            N20.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch2.size[n+1]
            not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.r[!reached.decisionH0_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]<=Accept.threshold
          RejectedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]>=Reject.threshold
          reached.decisionH0_n.AR.l = AcceptedH0.underH0_n.AR.l|RejectedH0.underH0_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.l)){
            
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[AcceptedH0.underH0_n.AR.l]] = 'A'
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[RejectedH0.underH0_n.AR.l]] = 'R'
            N10.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch1.size[n+1]
            N20.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch2.size[n+1]
            not.reached.decisionH0.AR.l = not.reached.decisionH0.AR.l[!reached.decisionH0_n.AR.l]
          }
          
          not.reached.decisionH0.AR = union(not.reached.decisionH0.AR.r,
                                            not.reached.decisionH0.AR.l)
        }
        
        
        ## under right-sided H1
        if(length(not.reached.decisionH1r.AR)>0){
          
          ## sum of observations at step n
          # Group 1
          sum11r_n = rnorm(length(not.reached.decisionH1r.AR),
                           (batch1.size[n+1]-batch1.size[n])*(theta1$right/2),
                           sqrt(batch1.size[n+1]-batch1.size[n])*sigma1)
          
          # Group 2
          sum21r_n = rnorm(length(not.reached.decisionH1r.AR),
                           -(batch2.size[n+1]-batch2.size[n])*(theta1$right/2),
                           sqrt(batch2.size[n+1]-batch2.size[n])*sigma2)
          
          ## sum of observations until step n
          # Group 1
          cumsum11r_n[not.reached.decisionH1r.AR] = 
            cumsum11r_n[not.reached.decisionH1r.AR] + sum11r_n
          
          # Group 2
          cumsum21r_n[not.reached.decisionH1r.AR] = 
            cumsum21r_n[not.reached.decisionH1r.AR] + sum21r_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1r_n.r[not.reached.decisionH1r.AR.r] = 
            exp(-(((theta.UMPBT$right^2) - (theta0^2)) - 
                    2*(theta.UMPBT$right - theta0)*
                    (cumsum11r_n[not.reached.decisionH1r.AR.r]/batch1.size[n+1] - 
                       cumsum21r_n[not.reached.decisionH1r.AR.r]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          # for left sided check
          LR1r_n.l[not.reached.decisionH1r.AR.l] = 
            exp(-(((theta.UMPBT$left^2) - (theta0^2)) - 
                    2*(theta.UMPBT$left - theta0)*
                    (cumsum11r_n[not.reached.decisionH1r.AR.l]/batch1.size[n+1] - 
                       cumsum21r_n[not.reached.decisionH1r.AR.l]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]<=Accept.threshold
          RejectedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]>=Reject.threshold
          reached.decisionH1r_n.AR.r = AcceptedH0.underH1r_n.AR.r|RejectedH0.underH1r_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.r)){
            
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[AcceptedH0.underH1r_n.AR.r]] = 'A'
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[RejectedH0.underH1r_n.AR.r]] = 'R'
            N11r.AR.r[not.reached.decisionH1r.AR.r[reached.decisionH1r_n.AR.r]] = batch1.size[n+1]
            N21r.AR.r[not.reached.decisionH1r.AR.r[reached.decisionH1r_n.AR.r]] = batch2.size[n+1]
            not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.r[!reached.decisionH1r_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]<=Accept.threshold
          RejectedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]>=Reject.threshold
          reached.decisionH1r_n.AR.l = AcceptedH0.underH1r_n.AR.l|RejectedH0.underH1r_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.l)){
            
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[AcceptedH0.underH1r_n.AR.l]] = 'A'
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[RejectedH0.underH1r_n.AR.l]] = 'R'
            N11r.AR.l[not.reached.decisionH1r.AR.l[reached.decisionH1r_n.AR.l]] = batch1.size[n+1]
            N21r.AR.l[not.reached.decisionH1r.AR.l[reached.decisionH1r_n.AR.l]] = batch2.size[n+1]
            not.reached.decisionH1r.AR.l = not.reached.decisionH1r.AR.l[!reached.decisionH1r_n.AR.l]
          }
          
          not.reached.decisionH1r.AR = union(not.reached.decisionH1r.AR.r,
                                             not.reached.decisionH1r.AR.l)
        }
        
        
        ## under left-sided H1
        if(length(not.reached.decisionH1l.AR)>0){
          
          ## sum of observations at step n
          # Group 1
          sum11l_n = rnorm(length(not.reached.decisionH1l.AR),
                           (batch1.size[n+1]-batch1.size[n])*(theta1$left/2),
                           sqrt(batch1.size[n+1]-batch1.size[n])*sigma1)
          
          # Group 2
          sum21l_n = rnorm(length(not.reached.decisionH1l.AR),
                           -(batch2.size[n+1]-batch2.size[n])*(theta1$left/2),
                           sqrt(batch2.size[n+1]-batch2.size[n])*sigma2)
          
          ## sum of observations until step n
          # Group 1
          cumsum11l_n[not.reached.decisionH1l.AR] = 
            cumsum11l_n[not.reached.decisionH1l.AR] + sum11l_n
          
          # Group 2
          cumsum21l_n[not.reached.decisionH1l.AR] = 
            cumsum21l_n[not.reached.decisionH1l.AR] + sum21l_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1l_n.r[not.reached.decisionH1l.AR.r] = 
            exp(-(((theta.UMPBT$right^2) - (theta0^2)) - 
                    2*(theta.UMPBT$right - theta0)*
                    (cumsum11l_n[not.reached.decisionH1l.AR.r]/batch1.size[n+1] - 
                       cumsum21l_n[not.reached.decisionH1l.AR.r]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          # for left sided check
          LR1l_n.l[not.reached.decisionH1l.AR.l] = 
            exp(-(((theta.UMPBT$left^2) - (theta0^2)) - 
                    2*(theta.UMPBT$left - theta0)*
                    (cumsum11l_n[not.reached.decisionH1l.AR.l]/batch1.size[n+1] - 
                       cumsum21l_n[not.reached.decisionH1l.AR.l]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]<=Accept.threshold
          RejectedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]>=Reject.threshold
          reached.decisionH1l_n.AR.r = AcceptedH0.underH1l_n.AR.r|RejectedH0.underH1l_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.r)){
            
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[AcceptedH0.underH1l_n.AR.r]] = 'A'
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[RejectedH0.underH1l_n.AR.r]] = 'R'
            N11l.AR.r[not.reached.decisionH1l.AR.r[reached.decisionH1l_n.AR.r]] = batch1.size[n+1]
            N21l.AR.r[not.reached.decisionH1l.AR.r[reached.decisionH1l_n.AR.r]] = batch2.size[n+1]
            not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.r[!reached.decisionH1l_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]<=Accept.threshold
          RejectedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]>=Reject.threshold
          reached.decisionH1l_n.AR.l = AcceptedH0.underH1l_n.AR.l|RejectedH0.underH1l_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.l)){
            
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[AcceptedH0.underH1l_n.AR.l]] = 'A'
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[RejectedH0.underH1l_n.AR.l]] = 'R'
            N11l.AR.l[not.reached.decisionH1l.AR.l[reached.decisionH1l_n.AR.l]] = batch1.size[n+1]
            N21l.AR.l[not.reached.decisionH1l.AR.l[reached.decisionH1l_n.AR.l]] = batch2.size[n+1]
            not.reached.decisionH1l.AR.l = not.reached.decisionH1l.AR.l[!reached.decisionH1l_n.AR.l]
          }
          
          not.reached.decisionH1l.AR = union(not.reached.decisionH1l.AR.r,
                                             not.reached.decisionH1l.AR.l)
        }
        
        setTxtProgressBar(pb, n)
      }
      
      
      ### both-sided checking
      ## under H0
      # accepted or rejected ones
      accepted.by.both0 = intersect(which(decision.underH0.AR.r=='A'),
                                    which(decision.underH0.AR.l=='A'))
      onlyrejected.by.right0 = intersect(which(decision.underH0.AR.r=='R'),
                                         which(decision.underH0.AR.l!='R'))
      onlyrejected.by.left0 = intersect(which(decision.underH0.AR.r!='R'),
                                        which(decision.underH0.AR.l=='R'))
      rejected.by.both0 = intersect(which(decision.underH0.AR.r=='R'),
                                    which(decision.underH0.AR.l=='R'))
      
      ## sample sizes required
      # Group 1
      N10.AR[accepted.by.both0] = pmax(N10.AR.r[accepted.by.both0],
                                       N10.AR.l[accepted.by.both0])
      N10.AR[onlyrejected.by.right0] = N10.AR.r[onlyrejected.by.right0]
      N10.AR[onlyrejected.by.left0] = N10.AR.l[onlyrejected.by.left0]
      N10.AR[rejected.by.both0] = pmin(N10.AR.r[rejected.by.both0],
                                       N10.AR.l[rejected.by.both0])
      
      # Group 2
      N20.AR[accepted.by.both0] = pmax(N20.AR.r[accepted.by.both0],
                                       N20.AR.l[accepted.by.both0])
      N20.AR[onlyrejected.by.right0] = N20.AR.r[onlyrejected.by.right0]
      N20.AR[onlyrejected.by.left0] = N20.AR.l[onlyrejected.by.left0]
      N20.AR[rejected.by.both0] = pmin(N20.AR.r[rejected.by.both0],
                                       N20.AR.l[rejected.by.both0])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right0 = intersect(which(decision.underH0.AR.r=='A'),
                                         which(is.na(decision.underH0.AR.l)))
      onlyaccepted.by.left0 = intersect(which(is.na(decision.underH0.AR.r)),
                                        which(decision.underH0.AR.l=='A'))
      both.inconclusive0 = intersect(which(is.na(decision.underH0.AR.r)),
                                     which(is.na(decision.underH0.AR.l)))
      all.inconclusive0 = c(onlyaccepted.by.right0, onlyaccepted.by.left0,
                            both.inconclusive0)
      nNot.reached.decisionH0.AR = length(all.inconclusive0)
      
      # Type I error probability
      type1.error.AR[c(onlyrejected.by.right0, onlyrejected.by.left0,
                       rejected.by.both0)] = T
      
      
      ## under right-sided H1
      # accepted or rejected ones
      accepted.by.both1r = intersect(which(decision.underH1r.AR.r=='A'),
                                     which(decision.underH1r.AR.l=='A'))
      onlyrejected.by.right1r = intersect(which(decision.underH1r.AR.r=='R'),
                                          which(decision.underH1r.AR.l!='R'))
      onlyrejected.by.left1r = intersect(which(decision.underH1r.AR.r!='R'),
                                         which(decision.underH1r.AR.l=='R'))
      rejected.by.both1r = intersect(which(decision.underH1r.AR.r=='R'),
                                     which(decision.underH1r.AR.l=='R'))
      
      ## sample sizes required
      # Group 1
      N11r.AR[accepted.by.both1r] = pmax(N11r.AR.r[accepted.by.both1r],
                                         N11r.AR.l[accepted.by.both1r])
      N11r.AR[onlyrejected.by.right1r] = N11r.AR.r[onlyrejected.by.right1r]
      N11r.AR[onlyrejected.by.left1r] = N11r.AR.l[onlyrejected.by.left1r]
      N11r.AR[rejected.by.both1r] = pmin(N11r.AR.r[rejected.by.both1r],
                                         N11r.AR.l[rejected.by.both1r])
      
      # Group 2
      N21r.AR[accepted.by.both1r] = pmax(N21r.AR.r[accepted.by.both1r],
                                         N21r.AR.l[accepted.by.both1r])
      N21r.AR[onlyrejected.by.right1r] = N21r.AR.r[onlyrejected.by.right1r]
      N21r.AR[onlyrejected.by.left1r] = N21r.AR.l[onlyrejected.by.left1r]
      N21r.AR[rejected.by.both1r] = pmin(N21r.AR.r[rejected.by.both1r],
                                         N21r.AR.l[rejected.by.both1r])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1r = intersect(which(decision.underH1r.AR.r=='A'),
                                          which(is.na(decision.underH1r.AR.l)))
      onlyaccepted.by.left1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                         which(decision.underH1r.AR.l=='A'))
      both.inconclusive1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                      which(is.na(decision.underH1r.AR.l)))
      all.inconclusive1r = c(onlyaccepted.by.right1r, onlyaccepted.by.left1r,
                             both.inconclusive1r)
      nNot.reached.decisionH1r.AR = length(all.inconclusive1r)
      
      # Type I error probability
      PowerH1r.AR[c(onlyrejected.by.right1r, onlyrejected.by.left1r,
                    rejected.by.both1r)] = T
      
      
      ## under left-sided H1
      # accepted or rejected ones
      accepted.by.both1l = intersect(which(decision.underH1l.AR.r=='A'),
                                     which(decision.underH1l.AR.l=='A'))
      onlyrejected.by.right1l = intersect(which(decision.underH1l.AR.r=='R'),
                                          which(decision.underH1l.AR.l!='R'))
      onlyrejected.by.left1l = intersect(which(decision.underH1l.AR.r!='R'),
                                         which(decision.underH1l.AR.l=='R'))
      rejected.by.both1l = intersect(which(decision.underH1l.AR.r=='R'),
                                     which(decision.underH1l.AR.l=='R'))
      
      ## sample sizes required
      # Group 1
      N11l.AR[accepted.by.both1l] = pmax(N11l.AR.r[accepted.by.both1l],
                                         N11l.AR.l[accepted.by.both1l])
      N11l.AR[onlyrejected.by.right1l] = N11l.AR.r[onlyrejected.by.right1l]
      N11l.AR[onlyrejected.by.left1l] = N11l.AR.l[onlyrejected.by.left1l]
      N11l.AR[rejected.by.both1l] = pmin(N11l.AR.r[rejected.by.both1l],
                                         N11l.AR.l[rejected.by.both1l])
      
      # Group 2
      N21l.AR[accepted.by.both1l] = pmax(N21l.AR.r[accepted.by.both1l],
                                         N21l.AR.l[accepted.by.both1l])
      N21l.AR[onlyrejected.by.right1l] = N21l.AR.r[onlyrejected.by.right1l]
      N21l.AR[onlyrejected.by.left1l] = N21l.AR.l[onlyrejected.by.left1l]
      N21l.AR[rejected.by.both1l] = pmin(N21l.AR.r[rejected.by.both1l],
                                         N21l.AR.l[rejected.by.both1l])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1l = intersect(which(decision.underH1l.AR.r=='A'),
                                          which(is.na(decision.underH1l.AR.l)))
      onlyaccepted.by.left1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                         which(decision.underH1l.AR.l=='A'))
      both.inconclusive1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                      which(is.na(decision.underH1l.AR.l)))
      all.inconclusive1l = c(onlyaccepted.by.right1l, onlyaccepted.by.left1l,
                             both.inconclusive1l)
      nNot.reached.decisionH1l.AR = length(all.inconclusive1l)
      
      # Type I error probability
      PowerH1l.AR[c(onlyrejected.by.right1l, onlyrejected.by.left1l,
                    rejected.by.both1l)] = T
      
      
      ## determining termination threshold
      ## H0 is rejected if LR or (BF) is >= termination threshold
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        term.thresh.possible.choices =
          c(LR0_n.r[onlyaccepted.by.left0],
            LR0_n.l[onlyaccepted.by.right0],
            pmin(LR0_n.r[both.inconclusive0], LR0_n.l[both.inconclusive0]))
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          max.LR0_n = max(term.thresh.possible.choices)
          nDecimal.accuracy = ceiling(-log10(min(0.01, Reject.threshold - max.LR0_n)))
          termination.threshold.AR = (floor(max.LR0_n*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01, min(term.thresh.possible.choices) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(term.thresh.possible.choices))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(term.thresh.possible.choices))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(cumRejFreq_not.reached.decisionH0.AR[1]>nNewRejects.AR){
            
            nDecimal.accuracy =
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR =
              (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                       (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR +
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      ## attained Type II error probability
      # right-sided H1
      actual.PowerH1r.AR.r = mean(PowerH1r.AR) +
        sum(c(LR1r_n.r[onlyaccepted.by.left1r],
              LR1r_n.l[onlyaccepted.by.right1r],
              pmax(LR1r_n.r[both.inconclusive1r], LR1r_n.l[both.inconclusive1r]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1r.AR = 1 - actual.PowerH1r.AR.r
      
      # left-sided H1
      actual.PowerH1l.AR.r = mean(PowerH1l.AR) +
        sum(c(LR1l_n.r[onlyaccepted.by.left1l],
              LR1l_n.l[onlyaccepted.by.right1l],
              pmax(LR1l_n.r[both.inconclusive1l], LR1l_n.l[both.inconclusive1l]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1l.AR = 1 - actual.PowerH1l.AR.r
      
      ## Expected sample sizes
      # Group 1
      EN10 = mean(N10.AR)     # under H0
      EN11r = mean(N11r.AR)   # under right-sided H1
      EN11l = mean(N11l.AR)   # under left-sided H1
      
      # Group 2
      EN20 = mean(N20.AR)     # under H0
      EN21r = mean(N21r.AR)   # under right-sided H1
      EN21l = mean(N21l.AR)   # under left-sided H1
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", round(termination.threshold.AR, 3)))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: Group 1 - ", round(EN10, 2), 
                    ', Group 2 - ', round(EN20, 2), sep = ''))
        print("Attained Type II error probability:")
        print(paste(" On the right: ", round(actual.type2.errorH1r.AR, 4)))
        print(paste(" On the left: ", round(actual.type2.errorH1l.AR, 4)))
        print("Expected sample size at the alternatives:")
        print(paste(" On the right: Group 1 - ", round(EN11r, 2), 
                    ', Group 2 - ', round(EN21r, 2), sep = ''))
        print(paste(" On the left: Group 1 - ", round(EN11l, 2), 
                    ', Group 2 - ', round(EN21l, 2), sep = ''))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  "Type2.attained" = c(actual.type2.errorH1r.AR, actual.type2.errorH1l.AR),
                  'N' = list('H0' = list('Group1' = N10.AR, 'Group2' = N20.AR),
                             'right' = list('Group1' = N11r.AR, 'Group2' = N21r.AR),
                             'left' = list('Group1' = N11l.AR, 'Group2' = N21l.AR)),
                  'EN' = list('H0' = list('Group1' = EN10, 'Group2' = EN20),
                              'right' = list('Group1' = EN11r, 'Group2' = EN21r),
                              'left' = list('Group1' = EN11l, 'Group2' = EN21l)),
                  "theta.UMPBT" = theta.UMPBT,
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'twoZ', 'side' = side, 'theta0' = theta0, 
                  'sigma1' = sigma1, 'sigma2' = sigma2, 'N1.max' = N1.max, 'N2.max' = N2.max,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'batch1.size' = diff(batch1.size), 'batch2.size' = diff(batch2.size),
                  'nAnalyses' = nAnalyses, 'nReplicate' = nReplicate, 'seed' = seed))
    }
  }
}

#### two-sample t test ####
design.MSPRT_twoT = function(side = 'right', theta0 = 0, theta1 = T,
                             Type1.target =.005, Type2.target = .2,
                             N1.max, N2.max, batch1.size, batch2.size,
                             nReplicate = 1e+6, verbose = T, seed = 1){
  
  
  if(side!='both'){
    
    ################################# two-sample t (right/left sided) #################################
    
    ## checking if length(batch1.size) and length(batch2.size) are equal
    if((!missing(batch1.size)) && (!missing(batch2.size)) &&
       (length(batch1.size)!=length(batch2.size))) return("Lenghts of batch1.size and batch2.size should be same")
    
    ## batch sizes and N for group 1
    if(missing(batch1.size)){
      
      if(missing(N1.max)){
        
        return("Either 'batch1.size' or 'N1.max' needs to be specified")
        
      }else{batch1.size = c(2, rep(1, N1.max-2))}
      
    }else{
      
      if(batch1.size[1]<2){
        
        return("First batch size in Group 1 should be at least 2")
        
      }else{
        
        if(missing(N1.max)){
          
          N1.max = sum(batch1.size)
          
        }else{
          
          if(sum(batch1.size)!=N1.max) return("Sum of batch1.size should add up to N1.max")
        }
      }
    }
    
    ## batch sizes and N for group 2
    if(missing(batch2.size)){
      
      if(missing(N2.max)){
        
        return("Either 'batch2.size' or 'N2.max' needs to be specified")
        
      }else{batch2.size = c(2, rep(1, N2.max-2))}
      
    }else{
      
      if(batch2.size[1]<2){
        
        return("First batch size in Group 2 should be at least 2")
        
      }else{
        
        if(missing(N2.max)){
          
          N2.max = sum(batch2.size)
          
        }else{
          
          if(sum(batch2.size)!=N2.max) return("Sum of batch2.size should add up to N2.max")
        }
      }
    }
    
    nAnalyses = length(batch1.size)
    
    ## msg
    if(verbose){
      
      if((batch1.size[1]>2)||any(batch1.size[-1]>1)||
         (batch2.size[1]>2)||any(batch2.size[-1]>1)){
        
        cat('\n')
        print("=========================================================================")
        print("Designing the group sequential MSPRT for a two-sample t test:")
        print("=========================================================================")
        
      }else{
        
        cat('\n')
        print("=========================================================================")
        print("Designing the sequential MSPRT for a two-sample t test:")
        print("=========================================================================")
      }
      
      print("Group 1:")
      print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
      print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
      print("Group 2:")
      print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
      print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
      print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
      print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
      print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
      print(paste("Hypothesized value under H0: ", theta0, sep = ""))
      print(paste("Direction of the H1: ", side, sep = ""))
    }
    
    batch1.size = c(0, cumsum(batch1.size))
    batch2.size = c(0, cumsum(batch2.size))
    
    
    if(is.logical(theta1)&&(theta1==F)){
      
      ################ no alternative comparison ################
      
      # msg
      if(verbose==T){
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # cut-off (with sign) in fixed design one-sample t test
      signed_t.alpha = (2*(side=='right')-1)*
        qt(Type1.target, df = N1.max + N2.max -2, lower.tail = F)
      
      # required storages
      cumSS10_n = cumSS20_n = cumsum10_n = cumsum20_n = 
        LR0_n = LR1_n = numeric(nReplicate)
      type1.error.AR = rep(F, nReplicate)
      N10.AR = rep(N1.max, nReplicate)
      N20.AR = rep(N2.max, nReplicate)
      not.reached.decisionH0.AR = 1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          ## observations at step n
          # Group 1
          obs10_n = mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                           FUN = function(X){
                             
                             rnorm(length(not.reached.decisionH0.AR), theta0/2, 1)
                           })
          # Group 2
          obs20_n = mapply(X = 1:(batch2.size[n+1]-batch2.size[n]), 
                           FUN = function(X){
                             
                             rnorm(length(not.reached.decisionH0.AR), -theta0/2, 1)
                           })
          
          ## sum of observations until step n
          # Group 1
          cumsum10_n[not.reached.decisionH0.AR] = 
            cumsum10_n[not.reached.decisionH0.AR] + rowSums(obs10_n)
          # Group 2
          cumsum20_n[not.reached.decisionH0.AR] = 
            cumsum20_n[not.reached.decisionH0.AR] + rowSums(obs20_n)
          
          ## sum of squares of observations until step n
          # Group 1
          cumSS10_n[not.reached.decisionH0.AR] = 
            cumSS10_n[not.reached.decisionH0.AR] + rowSums(obs10_n^2)
          # Group 2
          cumSS20_n[not.reached.decisionH0.AR] = 
            cumSS20_n[not.reached.decisionH0.AR] + rowSums(obs20_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          xbar.diff0_n = cumsum10_n[not.reached.decisionH0.AR]/batch1.size[n+1] -
            cumsum20_n[not.reached.decisionH0.AR]/batch2.size[n+1]
          divisor.pooled.sd0_n.sq = 
            cumSS10_n[not.reached.decisionH0.AR] - ((cumsum10_n[not.reached.decisionH0.AR])^2)/batch1.size[n+1] +
            cumSS20_n[not.reached.decisionH0.AR] - ((cumsum20_n[not.reached.decisionH0.AR])^2)/batch2.size[n+1]
          
          # likelihood ratio of observations until step n
          LR0_n[not.reached.decisionH0.AR] = 
            ((1 + ((xbar.diff0_n - theta0)^2)/(divisor.pooled.sd0_n.sq*
                                                 (1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff0_n - 
                        (theta0 + signed_t.alpha*
                           sqrt((divisor.pooled.sd0_n.sq/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd0_n.sq*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          # comparing with the thresholds
          AcceptedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]<=Accept.threshold)
          RejectedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]>=Reject.threshold)
          reached.decisionH0_n.AR = union(AcceptedH0.underH0_n.AR, RejectedH0.underH0_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH0_n.AR)>0){
            
            N10.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch1.size[n+1]
            N20.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch2.size[n+1]
            type1.error.AR[not.reached.decisionH0.AR[RejectedH0.underH0_n.AR]] = T
            not.reached.decisionH0.AR = not.reached.decisionH0.AR[-reached.decisionH0_n.AR]
          }
        }
        
        setTxtProgressBar(pb, n)
      }
      
      # determining termination threshold
      # H0 is rejected if LR or (BF) is >= termination threshold
      nNot.reached.decisionH0.AR = length(not.reached.decisionH0.AR)
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 Reject.threshold -
                                                   max(LR0_n[not.reached.decisionH0.AR]))))
          termination.threshold.AR = (floor(max(LR0_n[not.reached.decisionH0.AR])*
                                              (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 min(LR0_n[not.reached.decisionH0.AR]) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(LR0_n[not.reached.decisionH0.AR]))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(LR0_n[not.reached.decisionH0.AR]))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(min(cumRejFreq_not.reached.decisionH0.AR)>nNewRejects.AR){
            
            nDecimal.accuracy = 
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR + 
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      # Expected sample sizes
      EN10 = mean(N10.AR)
      EN20 = mean(N20.AR)
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", Reject.threshold))
        print(paste("Termination threshold: ", termination.threshold.AR))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste(" Expected sample size under H0: Group 1 - ", round(EN10, 2),
                    ", Group 2 - ", round(EN20, 2), sep = ''))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  'N' = list('H0' = list('Group1' = N10.AR, 'Group2' = N20.AR)),
                  'EN' = list('H0' = list('Group1' = EN10, 'Group2' = EN20)),
                  "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'twoT', 'side' = side, 'theta0' = theta0, 
                  'N1.max' = N1.max, 'N2.max' = N2.max,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'batch1.size' = diff(batch1.size), 'batch2.size' = diff(batch2.size),
                  'nAnalyses' = nAnalyses, 'nReplicate' = nReplicate, 'seed' = seed))
      
    }else if(is.logical(theta1)&&(theta1==T)){
      
      ################ comparison at the fixed-design alternative ################
      
      theta1 = fixed_design.alt(test.type = 'twoT', side = side, theta0 = theta0,
                                N1 = N1.max, N2 = N2.max, Type1 = Type1.target,
                                Type2 = Type2.target)
      
      # msg
      if(verbose==T){
        print(paste("Alternative under comparison: ", round(theta1, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # cut-off (with sign) in fixed design one-sample t test
      signed_t.alpha = (2*(side=='right')-1)*
        qt(Type1.target, df = N1.max + N2.max -2, lower.tail = F)
      
      # required storages
      cumSS10_n = cumSS20_n = cumSS11_n = cumSS21_n = 
        cumsum10_n = cumsum20_n = cumsum11_n = cumsum21_n = 
        LR0_n = LR1_n = numeric(nReplicate)
      type1.error.AR = type2.error.AR = rep(F, nReplicate)
      N10.AR = N11.AR = rep(N1.max, nReplicate)
      N20.AR = N21.AR = rep(N2.max, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH1.AR = 1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          ## observations at step n
          # Group 1
          obs10_n = mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                           FUN = function(X){
                             
                             rnorm(length(not.reached.decisionH0.AR), theta0/2, 1)
                           })
          # Group 2
          obs20_n = mapply(X = 1:(batch2.size[n+1]-batch2.size[n]), 
                           FUN = function(X){
                             
                             rnorm(length(not.reached.decisionH0.AR), -theta0/2, 1)
                           })
          
          ## sum of observations until step n
          # Group 1
          cumsum10_n[not.reached.decisionH0.AR] = 
            cumsum10_n[not.reached.decisionH0.AR] + rowSums(obs10_n)
          # Group 2
          cumsum20_n[not.reached.decisionH0.AR] = 
            cumsum20_n[not.reached.decisionH0.AR] + rowSums(obs20_n)
          
          ## sum of squares of observations until step n
          # Group 1
          cumSS10_n[not.reached.decisionH0.AR] = 
            cumSS10_n[not.reached.decisionH0.AR] + rowSums(obs10_n^2)
          # Group 2
          cumSS20_n[not.reached.decisionH0.AR] = 
            cumSS20_n[not.reached.decisionH0.AR] + rowSums(obs20_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          xbar.diff0_n = cumsum10_n[not.reached.decisionH0.AR]/batch1.size[n+1] -
            cumsum20_n[not.reached.decisionH0.AR]/batch2.size[n+1]
          divisor.pooled.sd0_n.sq = 
            cumSS10_n[not.reached.decisionH0.AR] - ((cumsum10_n[not.reached.decisionH0.AR])^2)/batch1.size[n+1] +
            cumSS20_n[not.reached.decisionH0.AR] - ((cumsum20_n[not.reached.decisionH0.AR])^2)/batch2.size[n+1]
          
          # likelihood ratio of observations until step n
          LR0_n[not.reached.decisionH0.AR] = 
            ((1 + ((xbar.diff0_n - theta0)^2)/(divisor.pooled.sd0_n.sq*
                                                 (1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff0_n - 
                        (theta0 + signed_t.alpha*
                           sqrt((divisor.pooled.sd0_n.sq/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd0_n.sq*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          # comparing with the thresholds
          AcceptedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]<=Accept.threshold)
          RejectedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]>=Reject.threshold)
          reached.decisionH0_n.AR = union(AcceptedH0.underH0_n.AR, RejectedH0.underH0_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH0_n.AR)>0){
            
            N10.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch1.size[n+1]
            N20.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch2.size[n+1]
            type1.error.AR[not.reached.decisionH0.AR[RejectedH0.underH0_n.AR]] = T
            not.reached.decisionH0.AR = not.reached.decisionH0.AR[-reached.decisionH0_n.AR]
          }
        }
        
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          ## observations at step n
          # Group 1
          obs11_n = mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                           FUN = function(X){
                             
                             rnorm(length(not.reached.decisionH1.AR), theta1/2, 1)
                           })
          # Group 2
          obs21_n = mapply(X = 1:(batch2.size[n+1]-batch2.size[n]), 
                           FUN = function(X){
                             
                             rnorm(length(not.reached.decisionH1.AR), -theta1/2, 1)
                           })
          
          ## sum of observations until step n
          # Group 1
          cumsum11_n[not.reached.decisionH1.AR] = 
            cumsum11_n[not.reached.decisionH1.AR] + rowSums(obs11_n)
          # Group 2
          cumsum21_n[not.reached.decisionH1.AR] = 
            cumsum21_n[not.reached.decisionH1.AR] + rowSums(obs21_n)
          
          ## sum of squares of observations until step n
          # Group 1
          cumSS11_n[not.reached.decisionH1.AR] = 
            cumSS11_n[not.reached.decisionH1.AR] + rowSums(obs11_n^2)
          # Group 2
          cumSS21_n[not.reached.decisionH1.AR] = 
            cumSS21_n[not.reached.decisionH1.AR] + rowSums(obs21_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          xbar.diff1_n = cumsum11_n[not.reached.decisionH1.AR]/batch1.size[n+1] -
            cumsum21_n[not.reached.decisionH1.AR]/batch2.size[n+1]
          divisor.pooled.sd1_n.sq = 
            cumSS11_n[not.reached.decisionH1.AR] - ((cumsum11_n[not.reached.decisionH1.AR])^2)/batch1.size[n+1] +
            cumSS21_n[not.reached.decisionH1.AR] - ((cumsum21_n[not.reached.decisionH1.AR])^2)/batch2.size[n+1]
          
          # likelihood ratio of observations until step n
          LR1_n[not.reached.decisionH1.AR] = 
            ((1 + ((xbar.diff1_n - theta0)^2)/(divisor.pooled.sd1_n.sq*
                                                 (1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff1_n - 
                        (theta0 + signed_t.alpha*
                           sqrt((divisor.pooled.sd1_n.sq/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd1_n.sq*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          # comparing with the thresholds
          AcceptedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]<=Accept.threshold)
          RejectedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]>=Reject.threshold)
          reached.decisionH1_n.AR = union(AcceptedH0.underH1_n.AR, RejectedH0.underH1_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH1_n.AR)>0){
            
            N11.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch1.size[n+1]
            N21.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch2.size[n+1]
            type2.error.AR[not.reached.decisionH1.AR[AcceptedH0.underH1_n.AR]] = T
            not.reached.decisionH1.AR = not.reached.decisionH1.AR[-reached.decisionH1_n.AR]
          }
        }
        
        setTxtProgressBar(pb, n)
      }
      
      # determining termination threshold
      # H0 is rejected if LR or (BF) is >= termination threshold
      nNot.reached.decisionH0.AR = length(not.reached.decisionH0.AR)
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 Reject.threshold -
                                                   max(LR0_n[not.reached.decisionH0.AR]))))
          termination.threshold.AR = (floor(max(LR0_n[not.reached.decisionH0.AR])*
                                              (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 min(LR0_n[not.reached.decisionH0.AR]) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(LR0_n[not.reached.decisionH0.AR]))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(LR0_n[not.reached.decisionH0.AR]))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(min(cumRejFreq_not.reached.decisionH0.AR)>nNewRejects.AR){
            
            nDecimal.accuracy = 
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR + 
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      # attained Type II error probability
      actual.type2.error.AR = mean(type2.error.AR) +
        sum(LR1_n[not.reached.decisionH1.AR]<termination.threshold.AR)/nReplicate
      
      # Expected sample sizes
      EN10 = mean(N10.AR)
      EN20 = mean(N20.AR)
      EN11 = mean(N11.AR)
      EN21 = mean(N21.AR)
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", Reject.threshold))
        print(paste("Termination threshold: ", termination.threshold.AR))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Attained Type II error probability: ", round(actual.type2.error.AR, 4)))
        print(paste(" Expected sample size under H0: Group 1 - ", round(EN10, 2),
                    ", Group 2 - ", round(EN20, 2), sep = ''))
        print(paste(" Expected sample size at the alternative: Group 1 - ", round(EN11, 2),
                    ", Group 2 - ", round(EN21, 2), sep = ''))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  "Type2.attained" = actual.type2.error.AR,
                  'N' = list('H0' = list('Group1' = N10.AR, 'Group2' = N20.AR),
                             'H1' = list('Group1' = N11.AR, 'Group2' = N21.AR)),
                  'EN' = list('H0' = list('Group1' = EN10, 'Group2' = EN20),
                              'H1' = list('Group1' = EN11, 'Group2' = EN21)),
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'twoT', 'side' = side, 'theta0' = theta0, 
                  'N1.max' = N1.max, 'N2.max' = N2.max,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'batch1.size' = diff(batch1.size), 'batch2.size' = diff(batch2.size),
                  'nAnalyses' = nAnalyses, 'nReplicate' = nReplicate, 'seed' = seed))
      
    }else{
      
      ################ comparison at the use specified point alternative ################
      
      # msg
      if(verbose==T){
        print(paste("Alternative under comparison: ", round(theta1, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # cut-off (with sign) in fixed design one-sample t test
      signed_t.alpha = (2*(side=='right')-1)*
        qt(Type1.target, df = N1.max + N2.max -2, lower.tail = F)
      
      # required storages
      cumSS10_n = cumSS20_n = cumSS11_n = cumSS21_n = 
        cumsum10_n = cumsum20_n = cumsum11_n = cumsum21_n = 
        LR0_n = LR1_n = numeric(nReplicate)
      type1.error.AR = type2.error.AR = rep(F, nReplicate)
      N10.AR = N11.AR = rep(N1.max, nReplicate)
      N20.AR = N21.AR = rep(N2.max, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH1.AR = 1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          ## observations at step n
          # Group 1
          obs10_n = mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                           FUN = function(X){
                             
                             rnorm(length(not.reached.decisionH0.AR), theta0/2, 1)
                           })
          # Group 2
          obs20_n = mapply(X = 1:(batch2.size[n+1]-batch2.size[n]), 
                           FUN = function(X){
                             
                             rnorm(length(not.reached.decisionH0.AR), -theta0/2, 1)
                           })
          
          ## sum of observations until step n
          # Group 1
          cumsum10_n[not.reached.decisionH0.AR] = 
            cumsum10_n[not.reached.decisionH0.AR] + rowSums(obs10_n)
          # Group 2
          cumsum20_n[not.reached.decisionH0.AR] = 
            cumsum20_n[not.reached.decisionH0.AR] + rowSums(obs20_n)
          
          ## sum of squares of observations until step n
          # Group 1
          cumSS10_n[not.reached.decisionH0.AR] = 
            cumSS10_n[not.reached.decisionH0.AR] + rowSums(obs10_n^2)
          # Group 2
          cumSS20_n[not.reached.decisionH0.AR] = 
            cumSS20_n[not.reached.decisionH0.AR] + rowSums(obs20_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          xbar.diff0_n = cumsum10_n[not.reached.decisionH0.AR]/batch1.size[n+1] -
            cumsum20_n[not.reached.decisionH0.AR]/batch2.size[n+1]
          divisor.pooled.sd0_n.sq = 
            cumSS10_n[not.reached.decisionH0.AR] - ((cumsum10_n[not.reached.decisionH0.AR])^2)/batch1.size[n+1] +
            cumSS20_n[not.reached.decisionH0.AR] - ((cumsum20_n[not.reached.decisionH0.AR])^2)/batch2.size[n+1]
          
          # likelihood ratio of observations until step n
          LR0_n[not.reached.decisionH0.AR] = 
            ((1 + ((xbar.diff0_n - theta0)^2)/(divisor.pooled.sd0_n.sq*
                                                 (1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff0_n - 
                        (theta0 + signed_t.alpha*
                           sqrt((divisor.pooled.sd0_n.sq/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd0_n.sq*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          # comparing with the thresholds
          AcceptedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]<=Accept.threshold)
          RejectedH0.underH0_n.AR = which(LR0_n[not.reached.decisionH0.AR]>=Reject.threshold)
          reached.decisionH0_n.AR = union(AcceptedH0.underH0_n.AR, RejectedH0.underH0_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH0_n.AR)>0){
            
            N10.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch1.size[n+1]
            N20.AR[not.reached.decisionH0.AR[reached.decisionH0_n.AR]] = batch2.size[n+1]
            type1.error.AR[not.reached.decisionH0.AR[RejectedH0.underH0_n.AR]] = T
            not.reached.decisionH0.AR = not.reached.decisionH0.AR[-reached.decisionH0_n.AR]
          }
        }
        
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          ## observations at step n
          # Group 1
          obs11_n = mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                           FUN = function(X){
                             
                             rnorm(length(not.reached.decisionH1.AR), theta1/2, 1)
                           })
          # Group 2
          obs21_n = mapply(X = 1:(batch2.size[n+1]-batch2.size[n]), 
                           FUN = function(X){
                             
                             rnorm(length(not.reached.decisionH1.AR), -theta1/2, 1)
                           })
          
          ## sum of observations until step n
          # Group 1
          cumsum11_n[not.reached.decisionH1.AR] = 
            cumsum11_n[not.reached.decisionH1.AR] + rowSums(obs11_n)
          # Group 2
          cumsum21_n[not.reached.decisionH1.AR] = 
            cumsum21_n[not.reached.decisionH1.AR] + rowSums(obs21_n)
          
          ## sum of squares of observations until step n
          # Group 1
          cumSS11_n[not.reached.decisionH1.AR] = 
            cumSS11_n[not.reached.decisionH1.AR] + rowSums(obs11_n^2)
          # Group 2
          cumSS21_n[not.reached.decisionH1.AR] = 
            cumSS21_n[not.reached.decisionH1.AR] + rowSums(obs21_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          xbar.diff1_n = cumsum11_n[not.reached.decisionH1.AR]/batch1.size[n+1] -
            cumsum21_n[not.reached.decisionH1.AR]/batch2.size[n+1]
          divisor.pooled.sd1_n.sq = 
            cumSS11_n[not.reached.decisionH1.AR] - ((cumsum11_n[not.reached.decisionH1.AR])^2)/batch1.size[n+1] +
            cumSS21_n[not.reached.decisionH1.AR] - ((cumsum21_n[not.reached.decisionH1.AR])^2)/batch2.size[n+1]
          
          # likelihood ratio of observations until step n
          LR1_n[not.reached.decisionH1.AR] = 
            ((1 + ((xbar.diff1_n - theta0)^2)/(divisor.pooled.sd1_n.sq*
                                                 (1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff1_n - 
                        (theta0 + signed_t.alpha*
                           sqrt((divisor.pooled.sd1_n.sq/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd1_n.sq*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          # comparing with the thresholds
          AcceptedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]<=Accept.threshold)
          RejectedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]>=Reject.threshold)
          reached.decisionH1_n.AR = union(AcceptedH0.underH1_n.AR, RejectedH0.underH1_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH1_n.AR)>0){
            
            N11.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch1.size[n+1]
            N21.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch2.size[n+1]
            type2.error.AR[not.reached.decisionH1.AR[AcceptedH0.underH1_n.AR]] = T
            not.reached.decisionH1.AR = not.reached.decisionH1.AR[-reached.decisionH1_n.AR]
          }
        }
        
        setTxtProgressBar(pb, n)
      }
      
      # determining termination threshold
      # H0 is rejected if LR or (BF) is >= termination threshold
      nNot.reached.decisionH0.AR = length(not.reached.decisionH0.AR)
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 Reject.threshold -
                                                   max(LR0_n[not.reached.decisionH0.AR]))))
          termination.threshold.AR = (floor(max(LR0_n[not.reached.decisionH0.AR])*
                                              (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                 min(LR0_n[not.reached.decisionH0.AR]) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(LR0_n[not.reached.decisionH0.AR]))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(LR0_n[not.reached.decisionH0.AR]))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(min(cumRejFreq_not.reached.decisionH0.AR)>nNewRejects.AR){
            
            nDecimal.accuracy = 
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR + 
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      # attained Type II error probability
      actual.type2.error.AR = mean(type2.error.AR) +
        sum(LR1_n[not.reached.decisionH1.AR]<termination.threshold.AR)/nReplicate
      
      # Expected sample sizes
      EN10 = mean(N10.AR)
      EN20 = mean(N20.AR)
      EN11 = mean(N11.AR)
      EN21 = mean(N21.AR)
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", Reject.threshold))
        print(paste("Termination threshold: ", termination.threshold.AR))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Attained Type II error probability: ", round(actual.type2.error.AR, 4)))
        print(paste(" Expected sample size under H0: Group 1 - ", round(EN10, 2),
                    ", Group 2 - ", round(EN20, 2), sep = ''))
        print(paste(" Expected sample size at the alternative: Group 1 - ", round(EN11, 2),
                    ", Group 2 - ", round(EN21, 2), sep = ''))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  "Type2.attained" = actual.type2.errorH1r.AR,
                  'N' = list('H0' = list('Group1' = N10.AR, 'Group2' = N20.AR),
                             'H1' = list('Group1' = N11.AR, 'Group2' = N21.AR)),
                  'EN' = list('H0' = list('Group1' = EN10, 'Group2' = EN20),
                              'H1' = list('Group1' = EN11, 'Group2' = EN21)),
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'twoT', 'side' = side, 'theta0' = theta0, 
                  'N1.max' = N1.max, 'N2.max' = N2.max,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'batch1.size' = diff(batch1.size), 'batch2.size' = diff(batch2.size),
                  'nAnalyses' = nAnalyses, 'nReplicate' = nReplicate, 'seed' = seed))
    }
    
    
  }else{
    
    ################################# two-sample t (both sided) #################################
    
    ## checking if length(batch1.size) and length(batch2.size) are equal
    if((!missing(batch1.size)) && (!missing(batch2.size)) &&
       (length(batch1.size)!=length(batch2.size))) return("Lenghts of batch1.size and batch2.size should be same")
    
    ## batch sizes and N for group 1
    if(missing(batch1.size)){
      
      if(missing(N1.max)){
        
        return("Either 'batch1.size' or 'N1.max' needs to be specified")
        
      }else{batch1.size = c(2, rep(1, N1.max-2))}
      
    }else{
      
      if(batch1.size[1]<2){
        
        return("First batch size in Group 1 should be at least 2")
        
      }else{
        
        if(missing(N1.max)){
          
          N1.max = sum(batch1.size)
          
        }else{
          
          if(sum(batch1.size)!=N1.max) return("Sum of batch1.size should add up to N1.max")
        }
      }
    }
    
    ## batch sizes and N for group 2
    if(missing(batch2.size)){
      
      if(missing(N2.max)){
        
        return("Either 'batch2.size' or 'N2.max' needs to be specified")
        
      }else{batch2.size = c(2, rep(1, N2.max-2))}
      
    }else{
      
      if(batch2.size[1]<2){
        
        return("First batch size in Group 2 should be at least 2")
        
      }else{
        
        if(missing(N2.max)){
          
          N2.max = sum(batch2.size)
          
        }else{
          
          if(sum(batch2.size)!=N2.max) return("Sum of batch2.size should add up to N2.max")
        }
      }
    }
    
    nAnalyses = length(batch1.size)
    
    ## msg
    if(verbose){
      
      if((batch1.size[1]>2)||any(batch1.size[-1]>1)||
         (batch2.size[1]>2)||any(batch2.size[-1]>1)){
        
        cat('\n')
        print("=========================================================================")
        print("Designing the group sequential MSPRT for a two-sample t test:")
        print("=========================================================================")
        
      }else{
        
        cat('\n')
        print("=========================================================================")
        print("Designing the sequential MSPRT for a two-sample t test:")
        print("=========================================================================")
      }
      
      print("Group 1:")
      print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
      print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
      print("Group 2:")
      print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
      print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
      print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
      print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
      print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
      print(paste("Hypothesized value under H0: ", theta0, sep = ""))
      print(paste("Direction of the H1: ", side, sep = ""))
    }
    
    batch1.size = c(0, cumsum(batch1.size))
    batch2.size = c(0, cumsum(batch2.size))
    
    
    if(is.logical(theta1)&&(theta1==F)){
      
      ################ no alternative comparison ################
      
      # msg
      if(verbose==T){
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # cut-off (with sign) in fixed design one-sample t test
      t.alpha = qt(Type1.target/2, df = N1.max + N2.max -2, lower.tail = F)
      
      # required storages
      cumSS10_n = cumSS20_n = cumsum10_n = cumsum20_n = LR0_n.r = LR0_n.l = numeric(nReplicate)
      type1.error.AR = rep(F, nReplicate)
      N10.AR = N10.AR.r = N10.AR.l = rep(N1.max, nReplicate)
      N20.AR = N20.AR.r = N20.AR.l = rep(N2.max, nReplicate)
      decision.underH0.AR.r = decision.underH0.AR.l = rep(NA, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.l =
        1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          ## observations at step n
          # Group 1
          if(length(not.reached.decisionH0.AR)>1){
            
            obs10_n = mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                             FUN = function(X){
                               
                               rnorm(length(not.reached.decisionH0.AR), theta0/2, 1)
                             })
            
          }else{
            
            obs10_n = matrix(mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                                    FUN = function(X){
                                      
                                      rnorm(length(not.reached.decisionH0.AR), theta0/2, 1)
                                      
                                    }), nrow = 1, ncol = batch1.size[n+1]-batch1.size[n], 
                             byrow = T)
          }
          
          # Group 2
          if(length(not.reached.decisionH0.AR)>1){
            
            obs20_n = mapply(X = 1:(batch2.size[n+1]-batch2.size[n]), 
                             FUN = function(X){
                               
                               rnorm(length(not.reached.decisionH0.AR), -theta0/2, 1)
                             })
            
          }else{
            
            obs20_n = matrix(mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                                    FUN = function(X){
                                      
                                      rnorm(length(not.reached.decisionH0.AR), -theta0/2, 1)
                                      
                                    }), nrow = 1, ncol = batch1.size[n+1]-batch1.size[n], 
                             byrow = T)
          }
          
          ## sum of observations until step n
          # Group 1
          cumsum10_n[not.reached.decisionH0.AR] = 
            cumsum10_n[not.reached.decisionH0.AR] + rowSums(obs10_n)
          # Group 2
          cumsum20_n[not.reached.decisionH0.AR] = 
            cumsum20_n[not.reached.decisionH0.AR] + rowSums(obs20_n)
          
          ## sum of squares of observations until step n
          # Group 1
          cumSS10_n[not.reached.decisionH0.AR] = 
            cumSS10_n[not.reached.decisionH0.AR] + rowSums(obs10_n^2)
          # Group 2
          cumSS20_n[not.reached.decisionH0.AR] = 
            cumSS20_n[not.reached.decisionH0.AR] + rowSums(obs20_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          # for right sided check
          xbar.diff0_n.r = cumsum10_n[not.reached.decisionH0.AR.r]/batch1.size[n+1] -
            cumsum20_n[not.reached.decisionH0.AR.r]/batch2.size[n+1]
          divisor.pooled.sd0_n.sq.r = 
            cumSS10_n[not.reached.decisionH0.AR.r] - 
            ((cumsum10_n[not.reached.decisionH0.AR.r])^2)/batch1.size[n+1] +
            cumSS20_n[not.reached.decisionH0.AR.r] - 
            ((cumsum20_n[not.reached.decisionH0.AR.r])^2)/batch2.size[n+1]
          
          # for left sided check
          xbar.diff0_n.l = cumsum10_n[not.reached.decisionH0.AR.l]/batch1.size[n+1] -
            cumsum20_n[not.reached.decisionH0.AR.l]/batch2.size[n+1]
          divisor.pooled.sd0_n.sq.l = 
            cumSS10_n[not.reached.decisionH0.AR.l] -
            ((cumsum10_n[not.reached.decisionH0.AR.l])^2)/batch1.size[n+1] +
            cumSS20_n[not.reached.decisionH0.AR.l] - 
            ((cumsum20_n[not.reached.decisionH0.AR.l])^2)/batch2.size[n+1]
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR0_n.r[not.reached.decisionH0.AR.r] = 
            ((1 + ((xbar.diff0_n.r - theta0)^2)/
                (divisor.pooled.sd0_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff0_n.r - 
                        (theta0 + t.alpha*
                           sqrt((divisor.pooled.sd0_n.sq.r/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd0_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          # for left sided check
          LR0_n.l[not.reached.decisionH0.AR.l] = 
            ((1 + ((xbar.diff0_n.l - theta0)^2)/
                (divisor.pooled.sd0_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff0_n.l - 
                        (theta0 - t.alpha*
                           sqrt((divisor.pooled.sd0_n.sq.l/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd0_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]<=Accept.threshold
          RejectedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]>=Reject.threshold
          reached.decisionH0_n.AR.r = AcceptedH0.underH0_n.AR.r|RejectedH0.underH0_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.r)){
            
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[AcceptedH0.underH0_n.AR.r]] = 'A'
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[RejectedH0.underH0_n.AR.r]] = 'R'
            N10.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch1.size[n+1]
            N20.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch2.size[n+1]
            not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.r[!reached.decisionH0_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]<=Accept.threshold
          RejectedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]>=Reject.threshold
          reached.decisionH0_n.AR.l = AcceptedH0.underH0_n.AR.l|RejectedH0.underH0_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.l)){
            
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[AcceptedH0.underH0_n.AR.l]] = 'A'
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[RejectedH0.underH0_n.AR.l]] = 'R'
            N10.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch1.size[n+1]
            N20.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch2.size[n+1]
            not.reached.decisionH0.AR.l = not.reached.decisionH0.AR.l[!reached.decisionH0_n.AR.l]
          }
          
          not.reached.decisionH0.AR = union(not.reached.decisionH0.AR.r,
                                            not.reached.decisionH0.AR.l)
        }
        
        setTxtProgressBar(pb, n)
      }
      
      
      ### both-sided checking
      ## under H0
      # accepted or rejected ones
      accepted.by.both0 = intersect(which(decision.underH0.AR.r=='A'),
                                    which(decision.underH0.AR.l=='A'))
      onlyrejected.by.right0 = intersect(which(decision.underH0.AR.r=='R'),
                                         which(decision.underH0.AR.l!='R'))
      onlyrejected.by.left0 = intersect(which(decision.underH0.AR.r!='R'),
                                        which(decision.underH0.AR.l=='R'))
      rejected.by.both0 = intersect(which(decision.underH0.AR.r=='R'),
                                    which(decision.underH0.AR.l=='R'))
      
      ## sample sizes required
      # Group 1
      N10.AR[accepted.by.both0] = pmax(N10.AR.r[accepted.by.both0],
                                       N10.AR.l[accepted.by.both0])
      N10.AR[onlyrejected.by.right0] = N10.AR.r[onlyrejected.by.right0]
      N10.AR[onlyrejected.by.left0] = N10.AR.l[onlyrejected.by.left0]
      N10.AR[rejected.by.both0] = pmin(N10.AR.r[rejected.by.both0],
                                       N10.AR.l[rejected.by.both0])
      
      # Group 2
      N20.AR[accepted.by.both0] = pmax(N20.AR.r[accepted.by.both0],
                                       N20.AR.l[accepted.by.both0])
      N20.AR[onlyrejected.by.right0] = N20.AR.r[onlyrejected.by.right0]
      N20.AR[onlyrejected.by.left0] = N20.AR.l[onlyrejected.by.left0]
      N20.AR[rejected.by.both0] = pmin(N20.AR.r[rejected.by.both0],
                                       N20.AR.l[rejected.by.both0])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right0 = intersect(which(decision.underH0.AR.r=='A'),
                                         which(is.na(decision.underH0.AR.l)))
      onlyaccepted.by.left0 = intersect(which(is.na(decision.underH0.AR.r)),
                                        which(decision.underH0.AR.l=='A'))
      both.inconclusive0 = intersect(which(is.na(decision.underH0.AR.r)),
                                     which(is.na(decision.underH0.AR.l)))
      all.inconclusive0 = c(onlyaccepted.by.right0, onlyaccepted.by.left0,
                            both.inconclusive0)
      nNot.reached.decisionH0.AR = length(all.inconclusive0)
      
      # Type I error probability
      type1.error.AR[c(onlyrejected.by.right0, onlyrejected.by.left0,
                       rejected.by.both0)] = T
      
      
      ## determining termination threshold
      ## H0 is rejected if LR or (BF) is >= termination threshold
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        term.thresh.possible.choices =
          c(LR0_n.r[onlyaccepted.by.left0],
            LR0_n.l[onlyaccepted.by.right0],
            pmin(LR0_n.r[both.inconclusive0], LR0_n.l[both.inconclusive0]))
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          max.LR0_n = max(term.thresh.possible.choices)
          nDecimal.accuracy = ceiling(-log10(min(0.01, Reject.threshold - max.LR0_n)))
          termination.threshold.AR = (floor(max.LR0_n*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01, min(term.thresh.possible.choices) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(term.thresh.possible.choices))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(term.thresh.possible.choices))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(cumRejFreq_not.reached.decisionH0.AR[1]>nNewRejects.AR){
            
            nDecimal.accuracy =
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR =
              (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                       (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR +
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      ## Expected sample sizes
      # Group 1
      EN10 = mean(N10.AR)     # under H0
      
      # Group 2
      EN20 = mean(N20.AR)     # under H0
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", round(termination.threshold.AR, 3)))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: Group 1 - ", round(EN10, 2), 
                    ', Group 2 - ', round(EN20, 2), sep = ''))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  'N' = list('H0' = list('Group1' = N10.AR, 'Group2' = N20.AR)),
                  'EN' = list('H0' = list('Group1' = EN10, 'Group2' = EN20)),
                  "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'twoT', 'side' = side, 'theta0' = theta0, 
                  'N1.max' = N1.max, 'N2.max' = N2.max,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'batch1.size' = diff(batch1.size), 'batch2.size' = diff(batch2.size),
                  'nAnalyses' = nAnalyses, 'nReplicate' = nReplicate, 'seed' = seed))
      
    }else if(is.logical(theta1)&&(theta1==T)){
      
      ################ comparison at the fixed-design alternative ################
      
      theta1 = list('right' = fixed_design.alt(test.type = 'twoT', side = 'right',
                                               theta0 = theta0, N1 = N1.max, N2 = N2.max,
                                               Type1 = Type1.target/2, Type2 = Type2.target),
                    'left' = fixed_design.alt(test.type = 'twoT', side = 'left',
                                              theta0 = theta0, N1 = N1.max, N2 = N2.max,
                                              Type1 = Type1.target/2, Type2 = Type2.target))
      
      # msg
      if(verbose==T){
        print("Alternative under comparison:")
        print(paste(' On the right: ', round(theta1$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta1$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # cut-off (with sign) in fixed design one-sample t test
      t.alpha = qt(Type1.target/2, df = N1.max + N2.max -2, lower.tail = F)
      
      # required storages
      cumSS10_n = cumSS20_n = cumSS11r_n = cumSS21r_n = cumSS11l_n = cumSS21l_n = 
        cumsum10_n = cumsum20_n = cumsum11r_n = cumsum21r_n = cumsum11l_n = cumsum21l_n = 
        LR0_n.r = LR0_n.l = LR1r_n.r = LR1r_n.l = LR1l_n.r = LR1l_n.l = numeric(nReplicate)
      type1.error.AR = PowerH1r.AR = PowerH1l.AR = rep(F, nReplicate)
      N10.AR = N10.AR.r = N10.AR.l = 
        N11r.AR = N11r.AR.r = N11r.AR.l = 
        N11l.AR = N11l.AR.r = N11l.AR.l = rep(N1.max, nReplicate)
      N20.AR = N20.AR.r = N20.AR.l = 
        N21r.AR = N21r.AR.r = N21r.AR.l = 
        N21l.AR = N21l.AR.r = N21l.AR.l = rep(N2.max, nReplicate)
      decision.underH0.AR.r = decision.underH0.AR.l = 
        decision.underH1r.AR.r = decision.underH1r.AR.l = 
        decision.underH1l.AR.r = decision.underH1l.AR.l = rep(NA, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.l =
        not.reached.decisionH1r.AR = not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.l =
        not.reached.decisionH1l.AR = not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.l =
        1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          ## observations at step n
          # Group 1
          if(length(not.reached.decisionH0.AR)>1){
            
            obs10_n = mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                             FUN = function(X){
                               
                               rnorm(length(not.reached.decisionH0.AR), theta0/2, 1)
                             })
            
          }else{
            
            obs10_n = matrix(mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                                    FUN = function(X){
                                      
                                      rnorm(length(not.reached.decisionH0.AR), theta0/2, 1)
                                      
                                    }), nrow = 1, ncol = batch1.size[n+1]-batch1.size[n], 
                             byrow = T)
          }
          
          # Group 2
          if(length(not.reached.decisionH0.AR)>1){
            
            obs20_n = mapply(X = 1:(batch2.size[n+1]-batch2.size[n]), 
                             FUN = function(X){
                               
                               rnorm(length(not.reached.decisionH0.AR), -theta0/2, 1)
                             })
            
          }else{
            
            obs20_n = matrix(mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                                    FUN = function(X){
                                      
                                      rnorm(length(not.reached.decisionH0.AR), -theta0/2, 1)
                                      
                                    }), nrow = 1, ncol = batch1.size[n+1]-batch1.size[n], 
                             byrow = T)
          }
          
          ## sum of observations until step n
          # Group 1
          cumsum10_n[not.reached.decisionH0.AR] = 
            cumsum10_n[not.reached.decisionH0.AR] + rowSums(obs10_n)
          # Group 2
          cumsum20_n[not.reached.decisionH0.AR] = 
            cumsum20_n[not.reached.decisionH0.AR] + rowSums(obs20_n)
          
          ## sum of squares of observations until step n
          # Group 1
          cumSS10_n[not.reached.decisionH0.AR] = 
            cumSS10_n[not.reached.decisionH0.AR] + rowSums(obs10_n^2)
          # Group 2
          cumSS20_n[not.reached.decisionH0.AR] = 
            cumSS20_n[not.reached.decisionH0.AR] + rowSums(obs20_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          # for right sided check
          xbar.diff0_n.r = cumsum10_n[not.reached.decisionH0.AR.r]/batch1.size[n+1] -
            cumsum20_n[not.reached.decisionH0.AR.r]/batch2.size[n+1]
          divisor.pooled.sd0_n.sq.r = 
            cumSS10_n[not.reached.decisionH0.AR.r] - 
            ((cumsum10_n[not.reached.decisionH0.AR.r])^2)/batch1.size[n+1] +
            cumSS20_n[not.reached.decisionH0.AR.r] - 
            ((cumsum20_n[not.reached.decisionH0.AR.r])^2)/batch2.size[n+1]
          
          # for left sided check
          xbar.diff0_n.l = cumsum10_n[not.reached.decisionH0.AR.l]/batch1.size[n+1] -
            cumsum20_n[not.reached.decisionH0.AR.l]/batch2.size[n+1]
          divisor.pooled.sd0_n.sq.l = 
            cumSS10_n[not.reached.decisionH0.AR.l] -
            ((cumsum10_n[not.reached.decisionH0.AR.l])^2)/batch1.size[n+1] +
            cumSS20_n[not.reached.decisionH0.AR.l] - 
            ((cumsum20_n[not.reached.decisionH0.AR.l])^2)/batch2.size[n+1]
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR0_n.r[not.reached.decisionH0.AR.r] = 
            ((1 + ((xbar.diff0_n.r - theta0)^2)/
                (divisor.pooled.sd0_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff0_n.r - 
                        (theta0 + t.alpha*
                           sqrt((divisor.pooled.sd0_n.sq.r/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd0_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          # for left sided check
          LR0_n.l[not.reached.decisionH0.AR.l] = 
            ((1 + ((xbar.diff0_n.l - theta0)^2)/
                (divisor.pooled.sd0_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff0_n.l - 
                        (theta0 - t.alpha*
                           sqrt((divisor.pooled.sd0_n.sq.l/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd0_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]<=Accept.threshold
          RejectedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]>=Reject.threshold
          reached.decisionH0_n.AR.r = AcceptedH0.underH0_n.AR.r|RejectedH0.underH0_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.r)){
            
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[AcceptedH0.underH0_n.AR.r]] = 'A'
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[RejectedH0.underH0_n.AR.r]] = 'R'
            N10.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch1.size[n+1]
            N20.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch2.size[n+1]
            not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.r[!reached.decisionH0_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]<=Accept.threshold
          RejectedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]>=Reject.threshold
          reached.decisionH0_n.AR.l = AcceptedH0.underH0_n.AR.l|RejectedH0.underH0_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.l)){
            
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[AcceptedH0.underH0_n.AR.l]] = 'A'
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[RejectedH0.underH0_n.AR.l]] = 'R'
            N10.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch1.size[n+1]
            N20.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch2.size[n+1]
            not.reached.decisionH0.AR.l = not.reached.decisionH0.AR.l[!reached.decisionH0_n.AR.l]
          }
          
          not.reached.decisionH0.AR = union(not.reached.decisionH0.AR.r,
                                            not.reached.decisionH0.AR.l)
        }
        
        
        ## under right-sided H1
        if(length(not.reached.decisionH1r.AR)>0){
          
          ## observations at step n
          # Group 1
          if(length(not.reached.decisionH1r.AR)>1){
            
            obs11r_n = mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                              FUN = function(X){
                                
                                rnorm(length(not.reached.decisionH1r.AR), theta1$right/2, 1)
                              })
            
          }else{
            
            obs11r_n = matrix(mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                                     FUN = function(X){
                                       
                                       rnorm(length(not.reached.decisionH1r.AR), theta1$right/2, 1)
                                       
                                     }), nrow = 1, ncol = batch1.size[n+1]-batch1.size[n], 
                              byrow = T)
          }
          
          # Group 2
          if(length(not.reached.decisionH1r.AR)>1){
            
            obs21r_n = mapply(X = 1:(batch2.size[n+1]-batch2.size[n]), 
                              FUN = function(X){
                                
                                rnorm(length(not.reached.decisionH1r.AR), -theta1$right/2, 1)
                              })
            
          }else{
            
            obs21r_n = matrix(mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                                     FUN = function(X){
                                       
                                       rnorm(length(not.reached.decisionH1r.AR), -theta1$right/2, 1)
                                       
                                     }), nrow = 1, ncol = batch1.size[n+1]-batch1.size[n], 
                              byrow = T)
          }
          
          ## sum of observations until step n
          # Group 1
          cumsum11r_n[not.reached.decisionH1r.AR] = 
            cumsum11r_n[not.reached.decisionH1r.AR] + rowSums(obs11r_n)
          # Group 2
          cumsum21r_n[not.reached.decisionH1r.AR] = 
            cumsum21r_n[not.reached.decisionH1r.AR] + rowSums(obs21r_n)
          
          ## sum of squares of observations until step n
          # Group 1
          cumSS11r_n[not.reached.decisionH1r.AR] = 
            cumSS11r_n[not.reached.decisionH1r.AR] + rowSums(obs11r_n^2)
          # Group 2
          cumSS21r_n[not.reached.decisionH1r.AR] = 
            cumSS21r_n[not.reached.decisionH1r.AR] + rowSums(obs21r_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          # for right sided check
          xbar.diff1r_n.r = cumsum11r_n[not.reached.decisionH1r.AR.r]/batch1.size[n+1] -
            cumsum21r_n[not.reached.decisionH1r.AR.r]/batch2.size[n+1]
          divisor.pooled.sd1r_n.sq.r = 
            cumSS11r_n[not.reached.decisionH1r.AR.r] - 
            ((cumsum11r_n[not.reached.decisionH1r.AR.r])^2)/batch1.size[n+1] +
            cumSS21r_n[not.reached.decisionH1r.AR.r] - 
            ((cumsum21r_n[not.reached.decisionH1r.AR.r])^2)/batch2.size[n+1]
          
          # for left sided check
          xbar.diff1r_n.l = cumsum11r_n[not.reached.decisionH1r.AR.l]/batch1.size[n+1] -
            cumsum21r_n[not.reached.decisionH1r.AR.l]/batch2.size[n+1]
          divisor.pooled.sd1r_n.sq.l = 
            cumSS11r_n[not.reached.decisionH1r.AR.l] -
            ((cumsum11r_n[not.reached.decisionH1r.AR.l])^2)/batch1.size[n+1] +
            cumSS21r_n[not.reached.decisionH1r.AR.l] - 
            ((cumsum21r_n[not.reached.decisionH1r.AR.l])^2)/batch2.size[n+1]
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1r_n.r[not.reached.decisionH1r.AR.r] = 
            ((1 + ((xbar.diff1r_n.r - theta0)^2)/
                (divisor.pooled.sd1r_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff1r_n.r - 
                        (theta0 + t.alpha*
                           sqrt((divisor.pooled.sd1r_n.sq.r/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd1r_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          # for left sided check
          LR1r_n.l[not.reached.decisionH1r.AR.l] = 
            ((1 + ((xbar.diff1r_n.l - theta0)^2)/
                (divisor.pooled.sd1r_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff1r_n.l - 
                        (theta0 - t.alpha*
                           sqrt((divisor.pooled.sd1r_n.sq.l/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd1r_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]<=Accept.threshold
          RejectedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]>=Reject.threshold
          reached.decisionH1r_n.AR.r = AcceptedH0.underH1r_n.AR.r|RejectedH0.underH1r_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.r)){
            
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[AcceptedH0.underH1r_n.AR.r]] = 'A'
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[RejectedH0.underH1r_n.AR.r]] = 'R'
            N11r.AR.r[not.reached.decisionH1r.AR.r[reached.decisionH1r_n.AR.r]] = batch1.size[n+1]
            N21r.AR.r[not.reached.decisionH1r.AR.r[reached.decisionH1r_n.AR.r]] = batch2.size[n+1]
            not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.r[!reached.decisionH1r_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]<=Accept.threshold
          RejectedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]>=Reject.threshold
          reached.decisionH1r_n.AR.l = AcceptedH0.underH1r_n.AR.l|RejectedH0.underH1r_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.l)){
            
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[AcceptedH0.underH1r_n.AR.l]] = 'A'
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[RejectedH0.underH1r_n.AR.l]] = 'R'
            N11r.AR.l[not.reached.decisionH1r.AR.l[reached.decisionH1r_n.AR.l]] = batch1.size[n+1]
            N21r.AR.l[not.reached.decisionH1r.AR.l[reached.decisionH1r_n.AR.l]] = batch2.size[n+1]
            not.reached.decisionH1r.AR.l = not.reached.decisionH1r.AR.l[!reached.decisionH1r_n.AR.l]
          }
          
          not.reached.decisionH1r.AR = union(not.reached.decisionH1r.AR.r,
                                             not.reached.decisionH1r.AR.l)
        }
        
        
        ## under left-sided H1
        if(length(not.reached.decisionH1l.AR)>0){
          
          ## observations at step n
          # Group 1
          if(length(not.reached.decisionH1l.AR)>1){
            
            obs11l_n = mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                              FUN = function(X){
                                
                                rnorm(length(not.reached.decisionH1l.AR), theta1$left/2, 1)
                              })
            
          }else{
            
            obs11l_n = matrix(mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                                     FUN = function(X){
                                       
                                       rnorm(length(not.reached.decisionH1l.AR), theta1$left/2, 1)
                                       
                                     }), nrow = 1, ncol = batch1.size[n+1]-batch1.size[n], 
                              byrow = T)
          }
          
          # Group 2
          if(length(not.reached.decisionH1l.AR)>1){
            
            obs21l_n = mapply(X = 1:(batch2.size[n+1]-batch2.size[n]), 
                              FUN = function(X){
                                
                                rnorm(length(not.reached.decisionH1l.AR), -theta1$left/2, 1)
                              })
            
          }else{
            
            obs21l_n = matrix(mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                                     FUN = function(X){
                                       
                                       rnorm(length(not.reached.decisionH1l.AR), -theta1$left/2, 1)
                                       
                                     }), nrow = 1, ncol = batch1.size[n+1]-batch1.size[n], 
                              byrow = T)
          }
          
          ## sum of observations until step n
          # Group 1
          cumsum11l_n[not.reached.decisionH1l.AR] = 
            cumsum11l_n[not.reached.decisionH1l.AR] + rowSums(obs11l_n)
          # Group 2
          cumsum21l_n[not.reached.decisionH1l.AR] = 
            cumsum21l_n[not.reached.decisionH1l.AR] + rowSums(obs21l_n)
          
          ## sum of squares of observations until step n
          # Group 1
          cumSS11l_n[not.reached.decisionH1l.AR] = 
            cumSS11l_n[not.reached.decisionH1l.AR] + rowSums(obs11l_n^2)
          # Group 2
          cumSS21l_n[not.reached.decisionH1l.AR] = 
            cumSS21l_n[not.reached.decisionH1l.AR] + rowSums(obs21l_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          # for right sided check
          xbar.diff1l_n.r = cumsum11l_n[not.reached.decisionH1l.AR.r]/batch1.size[n+1] -
            cumsum21l_n[not.reached.decisionH1l.AR.r]/batch2.size[n+1]
          divisor.pooled.sd1l_n.sq.r = 
            cumSS11l_n[not.reached.decisionH1l.AR.r] - 
            ((cumsum11l_n[not.reached.decisionH1l.AR.r])^2)/batch1.size[n+1] +
            cumSS21l_n[not.reached.decisionH1l.AR.r] - 
            ((cumsum21l_n[not.reached.decisionH1l.AR.r])^2)/batch2.size[n+1]
          
          # for left sided check
          xbar.diff1l_n.l = cumsum11l_n[not.reached.decisionH1l.AR.l]/batch1.size[n+1] -
            cumsum21l_n[not.reached.decisionH1l.AR.l]/batch2.size[n+1]
          divisor.pooled.sd1l_n.sq.l = 
            cumSS11l_n[not.reached.decisionH1l.AR.l] -
            ((cumsum11l_n[not.reached.decisionH1l.AR.l])^2)/batch1.size[n+1] +
            cumSS21l_n[not.reached.decisionH1l.AR.l] - 
            ((cumsum21l_n[not.reached.decisionH1l.AR.l])^2)/batch2.size[n+1]
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1l_n.r[not.reached.decisionH1l.AR.r] = 
            ((1 + ((xbar.diff1l_n.r - theta0)^2)/
                (divisor.pooled.sd1l_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff1l_n.r - 
                        (theta0 + t.alpha*
                           sqrt((divisor.pooled.sd1l_n.sq.r/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd1l_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          # for left sided check
          LR1l_n.l[not.reached.decisionH1l.AR.l] = 
            ((1 + ((xbar.diff1l_n.l - theta0)^2)/
                (divisor.pooled.sd1l_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff1l_n.l - 
                        (theta0 - t.alpha*
                           sqrt((divisor.pooled.sd1l_n.sq.l/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd1l_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]<=Accept.threshold
          RejectedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]>=Reject.threshold
          reached.decisionH1l_n.AR.r = AcceptedH0.underH1l_n.AR.r|RejectedH0.underH1l_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.r)){
            
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[AcceptedH0.underH1l_n.AR.r]] = 'A'
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[RejectedH0.underH1l_n.AR.r]] = 'R'
            N11l.AR.r[not.reached.decisionH1l.AR.r[reached.decisionH1l_n.AR.r]] = batch1.size[n+1]
            N21l.AR.r[not.reached.decisionH1l.AR.r[reached.decisionH1l_n.AR.r]] = batch2.size[n+1]
            not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.r[!reached.decisionH1l_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]<=Accept.threshold
          RejectedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]>=Reject.threshold
          reached.decisionH1l_n.AR.l = AcceptedH0.underH1l_n.AR.l|RejectedH0.underH1l_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.l)){
            
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[AcceptedH0.underH1l_n.AR.l]] = 'A'
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[RejectedH0.underH1l_n.AR.l]] = 'R'
            N11l.AR.l[not.reached.decisionH1l.AR.l[reached.decisionH1l_n.AR.l]] = batch1.size[n+1]
            N21l.AR.l[not.reached.decisionH1l.AR.l[reached.decisionH1l_n.AR.l]] = batch2.size[n+1]
            not.reached.decisionH1l.AR.l = not.reached.decisionH1l.AR.l[!reached.decisionH1l_n.AR.l]
          }
          
          not.reached.decisionH1l.AR = union(not.reached.decisionH1l.AR.r,
                                             not.reached.decisionH1l.AR.l)
        }
        
        setTxtProgressBar(pb, n)
      }
      
      
      ### both-sided checking
      ## under H0
      # accepted or rejected ones
      accepted.by.both0 = intersect(which(decision.underH0.AR.r=='A'),
                                    which(decision.underH0.AR.l=='A'))
      onlyrejected.by.right0 = intersect(which(decision.underH0.AR.r=='R'),
                                         which(decision.underH0.AR.l!='R'))
      onlyrejected.by.left0 = intersect(which(decision.underH0.AR.r!='R'),
                                        which(decision.underH0.AR.l=='R'))
      rejected.by.both0 = intersect(which(decision.underH0.AR.r=='R'),
                                    which(decision.underH0.AR.l=='R'))
      
      ## sample sizes required
      # Group 1
      N10.AR[accepted.by.both0] = pmax(N10.AR.r[accepted.by.both0],
                                       N10.AR.l[accepted.by.both0])
      N10.AR[onlyrejected.by.right0] = N10.AR.r[onlyrejected.by.right0]
      N10.AR[onlyrejected.by.left0] = N10.AR.l[onlyrejected.by.left0]
      N10.AR[rejected.by.both0] = pmin(N10.AR.r[rejected.by.both0],
                                       N10.AR.l[rejected.by.both0])
      
      # Group 2
      N20.AR[accepted.by.both0] = pmax(N20.AR.r[accepted.by.both0],
                                       N20.AR.l[accepted.by.both0])
      N20.AR[onlyrejected.by.right0] = N20.AR.r[onlyrejected.by.right0]
      N20.AR[onlyrejected.by.left0] = N20.AR.l[onlyrejected.by.left0]
      N20.AR[rejected.by.both0] = pmin(N20.AR.r[rejected.by.both0],
                                       N20.AR.l[rejected.by.both0])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right0 = intersect(which(decision.underH0.AR.r=='A'),
                                         which(is.na(decision.underH0.AR.l)))
      onlyaccepted.by.left0 = intersect(which(is.na(decision.underH0.AR.r)),
                                        which(decision.underH0.AR.l=='A'))
      both.inconclusive0 = intersect(which(is.na(decision.underH0.AR.r)),
                                     which(is.na(decision.underH0.AR.l)))
      all.inconclusive0 = c(onlyaccepted.by.right0, onlyaccepted.by.left0,
                            both.inconclusive0)
      nNot.reached.decisionH0.AR = length(all.inconclusive0)
      
      # Type I error probability
      type1.error.AR[c(onlyrejected.by.right0, onlyrejected.by.left0,
                       rejected.by.both0)] = T
      
      
      ## under right-sided H1
      # accepted or rejected ones
      accepted.by.both1r = intersect(which(decision.underH1r.AR.r=='A'),
                                     which(decision.underH1r.AR.l=='A'))
      onlyrejected.by.right1r = intersect(which(decision.underH1r.AR.r=='R'),
                                          which(decision.underH1r.AR.l!='R'))
      onlyrejected.by.left1r = intersect(which(decision.underH1r.AR.r!='R'),
                                         which(decision.underH1r.AR.l=='R'))
      rejected.by.both1r = intersect(which(decision.underH1r.AR.r=='R'),
                                     which(decision.underH1r.AR.l=='R'))
      
      ## sample sizes required
      # Group 1
      N11r.AR[accepted.by.both1r] = pmax(N11r.AR.r[accepted.by.both1r],
                                         N11r.AR.l[accepted.by.both1r])
      N11r.AR[onlyrejected.by.right1r] = N11r.AR.r[onlyrejected.by.right1r]
      N11r.AR[onlyrejected.by.left1r] = N11r.AR.l[onlyrejected.by.left1r]
      N11r.AR[rejected.by.both1r] = pmin(N11r.AR.r[rejected.by.both1r],
                                         N11r.AR.l[rejected.by.both1r])
      
      # Group 2
      N21r.AR[accepted.by.both1r] = pmax(N21r.AR.r[accepted.by.both1r],
                                         N21r.AR.l[accepted.by.both1r])
      N21r.AR[onlyrejected.by.right1r] = N21r.AR.r[onlyrejected.by.right1r]
      N21r.AR[onlyrejected.by.left1r] = N21r.AR.l[onlyrejected.by.left1r]
      N21r.AR[rejected.by.both1r] = pmin(N21r.AR.r[rejected.by.both1r],
                                         N21r.AR.l[rejected.by.both1r])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1r = intersect(which(decision.underH1r.AR.r=='A'),
                                          which(is.na(decision.underH1r.AR.l)))
      onlyaccepted.by.left1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                         which(decision.underH1r.AR.l=='A'))
      both.inconclusive1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                      which(is.na(decision.underH1r.AR.l)))
      all.inconclusive1r = c(onlyaccepted.by.right1r, onlyaccepted.by.left1r,
                             both.inconclusive1r)
      nNot.reached.decisionH1r.AR = length(all.inconclusive1r)
      
      # Type I error probability
      PowerH1r.AR[c(onlyrejected.by.right1r, onlyrejected.by.left1r,
                    rejected.by.both1r)] = T
      
      
      ## under left-sided H1
      # accepted or rejected ones
      accepted.by.both1l = intersect(which(decision.underH1l.AR.r=='A'),
                                     which(decision.underH1l.AR.l=='A'))
      onlyrejected.by.right1l = intersect(which(decision.underH1l.AR.r=='R'),
                                          which(decision.underH1l.AR.l!='R'))
      onlyrejected.by.left1l = intersect(which(decision.underH1l.AR.r!='R'),
                                         which(decision.underH1l.AR.l=='R'))
      rejected.by.both1l = intersect(which(decision.underH1l.AR.r=='R'),
                                     which(decision.underH1l.AR.l=='R'))
      
      ## sample sizes required
      # Group 1
      N11l.AR[accepted.by.both1l] = pmax(N11l.AR.r[accepted.by.both1l],
                                         N11l.AR.l[accepted.by.both1l])
      N11l.AR[onlyrejected.by.right1l] = N11l.AR.r[onlyrejected.by.right1l]
      N11l.AR[onlyrejected.by.left1l] = N11l.AR.l[onlyrejected.by.left1l]
      N11l.AR[rejected.by.both1l] = pmin(N11l.AR.r[rejected.by.both1l],
                                         N11l.AR.l[rejected.by.both1l])
      
      # Group 2
      N21l.AR[accepted.by.both1l] = pmax(N21l.AR.r[accepted.by.both1l],
                                         N21l.AR.l[accepted.by.both1l])
      N21l.AR[onlyrejected.by.right1l] = N21l.AR.r[onlyrejected.by.right1l]
      N21l.AR[onlyrejected.by.left1l] = N21l.AR.l[onlyrejected.by.left1l]
      N21l.AR[rejected.by.both1l] = pmin(N21l.AR.r[rejected.by.both1l],
                                         N21l.AR.l[rejected.by.both1l])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1l = intersect(which(decision.underH1l.AR.r=='A'),
                                          which(is.na(decision.underH1l.AR.l)))
      onlyaccepted.by.left1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                         which(decision.underH1l.AR.l=='A'))
      both.inconclusive1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                      which(is.na(decision.underH1l.AR.l)))
      all.inconclusive1l = c(onlyaccepted.by.right1l, onlyaccepted.by.left1l,
                             both.inconclusive1l)
      nNot.reached.decisionH1l.AR = length(all.inconclusive1l)
      
      # Type I error probability
      PowerH1l.AR[c(onlyrejected.by.right1l, onlyrejected.by.left1l,
                    rejected.by.both1l)] = T
      
      
      ## determining termination threshold
      ## H0 is rejected if LR or (BF) is >= termination threshold
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        term.thresh.possible.choices =
          c(LR0_n.r[onlyaccepted.by.left0],
            LR0_n.l[onlyaccepted.by.right0],
            pmin(LR0_n.r[both.inconclusive0], LR0_n.l[both.inconclusive0]))
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          max.LR0_n = max(term.thresh.possible.choices)
          nDecimal.accuracy = ceiling(-log10(min(0.01, Reject.threshold - max.LR0_n)))
          termination.threshold.AR = (floor(max.LR0_n*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01, min(term.thresh.possible.choices) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(term.thresh.possible.choices))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(term.thresh.possible.choices))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(cumRejFreq_not.reached.decisionH0.AR[1]>nNewRejects.AR){
            
            nDecimal.accuracy =
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR =
              (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                       (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR +
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      ## attained Type II error probability
      # right-sided H1
      actual.PowerH1r.AR.r = mean(PowerH1r.AR) +
        sum(c(LR1r_n.r[onlyaccepted.by.left1r],
              LR1r_n.l[onlyaccepted.by.right1r],
              pmax(LR1r_n.r[both.inconclusive1r], LR1r_n.l[both.inconclusive1r]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1r.AR = 1 - actual.PowerH1r.AR.r
      
      # left-sided H1
      actual.PowerH1l.AR.r = mean(PowerH1l.AR) +
        sum(c(LR1l_n.r[onlyaccepted.by.left1l],
              LR1l_n.l[onlyaccepted.by.right1l],
              pmax(LR1l_n.r[both.inconclusive1l], LR1l_n.l[both.inconclusive1l]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1l.AR = 1 - actual.PowerH1l.AR.r
      
      ## Expected sample sizes
      # Group 1
      EN10 = mean(N10.AR)     # under H0
      EN11r = mean(N11r.AR)   # under right-sided H1
      EN11l = mean(N11l.AR)   # under left-sided H1
      
      # Group 2
      EN20 = mean(N20.AR)     # under H0
      EN21r = mean(N21r.AR)   # under right-sided H1
      EN21l = mean(N21l.AR)   # under left-sided H1
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", round(termination.threshold.AR, 3)))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: Group 1 - ", round(EN10, 2), 
                    ', Group 2 - ', round(EN20, 2), sep = ''))
        print("Attained Type II error probability:")
        print(paste(" On the right: ", round(actual.type2.errorH1r.AR, 4)))
        print(paste(" On the left: ", round(actual.type2.errorH1l.AR, 4)))
        print("Expected sample size at the alternatives:")
        print(paste(" On the right: Group 1 - ", round(EN11r, 2), 
                    ', Group 2 - ', round(EN21r, 2), sep = ''))
        print(paste(" On the left: Group 1 - ", round(EN11l, 2), 
                    ', Group 2 - ', round(EN21l, 2), sep = ''))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  "Type2.attained" = c(actual.type2.errorH1r.AR, actual.type2.errorH1l.AR),
                  'N' = list('H0' = list('Group1' = N10.AR, 'Group2' = N20.AR),
                             'right' = list('Group1' = N11r.AR, 'Group2' = N21r.AR),
                             'left' = list('Group1' = N11l.AR, 'Group2' = N21l.AR)),
                  'EN' = list('H0' = list('Group1' = EN10, 'Group2' = EN20),
                              'right' = list('Group1' = EN11r, 'Group2' = EN21r),
                              'left' = list('Group1' = EN11l, 'Group2' = EN21l)),
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'twoT', 'side' = side, 'theta0' = theta0, 
                  'N1.max' = N1.max, 'N2.max' = N2.max,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'batch1.size' = diff(batch1.size), 'batch2.size' = diff(batch2.size),
                  'nAnalyses' = nAnalyses, 'nReplicate' = nReplicate, 'seed' = seed))
      
    }else{
      
      ################ comparison at the user-specified point alternative ################
      
      # msg
      if(verbose==T){
        print("Alternative under comparison: ")
        print("-------------------------------------------------------------------------")
        print(paste(' On the right: ', round(theta1$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta1$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the Termination threshold ...")
      }
      
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # cut-off (with sign) in fixed design one-sample t test
      t.alpha = qt(Type1.target/2, df = N1.max + N2.max -2, lower.tail = F)
      
      # required storages
      cumSS10_n = cumSS20_n = cumSS11r_n = cumSS21r_n = cumSS11l_n = cumSS21l_n = 
        cumsum10_n = cumsum20_n = cumsum11r_n = cumsum21r_n = cumsum11l_n = cumsum21l_n = 
        LR0_n.r = LR0_n.l = LR1r_n.r = LR1r_n.l = LR1l_n.r = LR1l_n.l = numeric(nReplicate)
      type1.error.AR = PowerH1r.AR = PowerH1l.AR = rep(F, nReplicate)
      N10.AR = N10.AR.r = N10.AR.l = 
        N11r.AR = N11r.AR.r = N11r.AR.l = 
        N11l.AR = N11l.AR.r = N11l.AR.l = rep(N1.max, nReplicate)
      N20.AR = N20.AR.r = N20.AR.l = 
        N21r.AR = N21r.AR.r = N21r.AR.l = 
        N21l.AR = N21l.AR.r = N21l.AR.l = rep(N2.max, nReplicate)
      decision.underH0.AR.r = decision.underH0.AR.l = 
        decision.underH1r.AR.r = decision.underH1r.AR.l = 
        decision.underH1l.AR.r = decision.underH1l.AR.l = rep(NA, nReplicate)
      not.reached.decisionH0.AR = not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.l =
        not.reached.decisionH1r.AR = not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.l =
        not.reached.decisionH1l.AR = not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.l =
        1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under H0
        if(length(not.reached.decisionH0.AR)>0){
          
          ## observations at step n
          # Group 1
          if(length(not.reached.decisionH0.AR)>1){
            
            obs10_n = mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                             FUN = function(X){
                               
                               rnorm(length(not.reached.decisionH0.AR), theta0/2, 1)
                             })
            
          }else{
            
            obs10_n = matrix(mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                                    FUN = function(X){
                                      
                                      rnorm(length(not.reached.decisionH0.AR), theta0/2, 1)
                                      
                                    }), nrow = 1, ncol = batch1.size[n+1]-batch1.size[n], 
                             byrow = T)
          }
          
          # Group 2
          if(length(not.reached.decisionH0.AR)>1){
            
            obs20_n = mapply(X = 1:(batch2.size[n+1]-batch2.size[n]), 
                             FUN = function(X){
                               
                               rnorm(length(not.reached.decisionH0.AR), -theta0/2, 1)
                             })
            
          }else{
            
            obs20_n = matrix(mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                                    FUN = function(X){
                                      
                                      rnorm(length(not.reached.decisionH0.AR), -theta0/2, 1)
                                      
                                    }), nrow = 1, ncol = batch1.size[n+1]-batch1.size[n], 
                             byrow = T)
          }
          
          ## sum of observations until step n
          # Group 1
          cumsum10_n[not.reached.decisionH0.AR] = 
            cumsum10_n[not.reached.decisionH0.AR] + rowSums(obs10_n)
          # Group 2
          cumsum20_n[not.reached.decisionH0.AR] = 
            cumsum20_n[not.reached.decisionH0.AR] + rowSums(obs20_n)
          
          ## sum of squares of observations until step n
          # Group 1
          cumSS10_n[not.reached.decisionH0.AR] = 
            cumSS10_n[not.reached.decisionH0.AR] + rowSums(obs10_n^2)
          # Group 2
          cumSS20_n[not.reached.decisionH0.AR] = 
            cumSS20_n[not.reached.decisionH0.AR] + rowSums(obs20_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          # for right sided check
          xbar.diff0_n.r = cumsum10_n[not.reached.decisionH0.AR.r]/batch1.size[n+1] -
            cumsum20_n[not.reached.decisionH0.AR.r]/batch2.size[n+1]
          divisor.pooled.sd0_n.sq.r = 
            cumSS10_n[not.reached.decisionH0.AR.r] - 
            ((cumsum10_n[not.reached.decisionH0.AR.r])^2)/batch1.size[n+1] +
            cumSS20_n[not.reached.decisionH0.AR.r] - 
            ((cumsum20_n[not.reached.decisionH0.AR.r])^2)/batch2.size[n+1]
          
          # for left sided check
          xbar.diff0_n.l = cumsum10_n[not.reached.decisionH0.AR.l]/batch1.size[n+1] -
            cumsum20_n[not.reached.decisionH0.AR.l]/batch2.size[n+1]
          divisor.pooled.sd0_n.sq.l = 
            cumSS10_n[not.reached.decisionH0.AR.l] -
            ((cumsum10_n[not.reached.decisionH0.AR.l])^2)/batch1.size[n+1] +
            cumSS20_n[not.reached.decisionH0.AR.l] - 
            ((cumsum20_n[not.reached.decisionH0.AR.l])^2)/batch2.size[n+1]
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR0_n.r[not.reached.decisionH0.AR.r] = 
            ((1 + ((xbar.diff0_n.r - theta0)^2)/
                (divisor.pooled.sd0_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff0_n.r - 
                        (theta0 + t.alpha*
                           sqrt((divisor.pooled.sd0_n.sq.r/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd0_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          # for left sided check
          LR0_n.l[not.reached.decisionH0.AR.l] = 
            ((1 + ((xbar.diff0_n.l - theta0)^2)/
                (divisor.pooled.sd0_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff0_n.l - 
                        (theta0 - t.alpha*
                           sqrt((divisor.pooled.sd0_n.sq.l/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd0_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]<=Accept.threshold
          RejectedH0.underH0_n.AR.r = LR0_n.r[not.reached.decisionH0.AR.r]>=Reject.threshold
          reached.decisionH0_n.AR.r = AcceptedH0.underH0_n.AR.r|RejectedH0.underH0_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.r)){
            
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[AcceptedH0.underH0_n.AR.r]] = 'A'
            decision.underH0.AR.r[not.reached.decisionH0.AR.r[RejectedH0.underH0_n.AR.r]] = 'R'
            N10.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch1.size[n+1]
            N20.AR.r[not.reached.decisionH0.AR.r[reached.decisionH0_n.AR.r]] = batch2.size[n+1]
            not.reached.decisionH0.AR.r = not.reached.decisionH0.AR.r[!reached.decisionH0_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]<=Accept.threshold
          RejectedH0.underH0_n.AR.l = LR0_n.l[not.reached.decisionH0.AR.l]>=Reject.threshold
          reached.decisionH0_n.AR.l = AcceptedH0.underH0_n.AR.l|RejectedH0.underH0_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH0_n.AR.l)){
            
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[AcceptedH0.underH0_n.AR.l]] = 'A'
            decision.underH0.AR.l[not.reached.decisionH0.AR.l[RejectedH0.underH0_n.AR.l]] = 'R'
            N10.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch1.size[n+1]
            N20.AR.l[not.reached.decisionH0.AR.l[reached.decisionH0_n.AR.l]] = batch2.size[n+1]
            not.reached.decisionH0.AR.l = not.reached.decisionH0.AR.l[!reached.decisionH0_n.AR.l]
          }
          
          not.reached.decisionH0.AR = union(not.reached.decisionH0.AR.r,
                                            not.reached.decisionH0.AR.l)
        }
        
        
        ## under right-sided H1
        if(length(not.reached.decisionH1r.AR)>0){
          
          ## observations at step n
          # Group 1
          if(length(not.reached.decisionH1r.AR)>1){
            
            obs11r_n = mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                              FUN = function(X){
                                
                                rnorm(length(not.reached.decisionH1r.AR), theta1$right/2, 1)
                              })
            
          }else{
            
            obs11r_n = matrix(mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                                     FUN = function(X){
                                       
                                       rnorm(length(not.reached.decisionH1r.AR), theta1$right/2, 1)
                                       
                                     }), nrow = 1, ncol = batch1.size[n+1]-batch1.size[n], 
                              byrow = T)
          }
          
          # Group 2
          if(length(not.reached.decisionH1r.AR)>1){
            
            obs21r_n = mapply(X = 1:(batch2.size[n+1]-batch2.size[n]), 
                              FUN = function(X){
                                
                                rnorm(length(not.reached.decisionH1r.AR), -theta1$right/2, 1)
                              })
            
          }else{
            
            obs21r_n = matrix(mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                                     FUN = function(X){
                                       
                                       rnorm(length(not.reached.decisionH1r.AR), -theta1$right/2, 1)
                                       
                                     }), nrow = 1, ncol = batch1.size[n+1]-batch1.size[n], 
                              byrow = T)
          }
          
          ## sum of observations until step n
          # Group 1
          cumsum11r_n[not.reached.decisionH1r.AR] = 
            cumsum11r_n[not.reached.decisionH1r.AR] + rowSums(obs11r_n)
          # Group 2
          cumsum21r_n[not.reached.decisionH1r.AR] = 
            cumsum21r_n[not.reached.decisionH1r.AR] + rowSums(obs21r_n)
          
          ## sum of squares of observations until step n
          # Group 1
          cumSS11r_n[not.reached.decisionH1r.AR] = 
            cumSS11r_n[not.reached.decisionH1r.AR] + rowSums(obs11r_n^2)
          # Group 2
          cumSS21r_n[not.reached.decisionH1r.AR] = 
            cumSS21r_n[not.reached.decisionH1r.AR] + rowSums(obs21r_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          # for right sided check
          xbar.diff1r_n.r = cumsum11r_n[not.reached.decisionH1r.AR.r]/batch1.size[n+1] -
            cumsum21r_n[not.reached.decisionH1r.AR.r]/batch2.size[n+1]
          divisor.pooled.sd1r_n.sq.r = 
            cumSS11r_n[not.reached.decisionH1r.AR.r] - 
            ((cumsum11r_n[not.reached.decisionH1r.AR.r])^2)/batch1.size[n+1] +
            cumSS21r_n[not.reached.decisionH1r.AR.r] - 
            ((cumsum21r_n[not.reached.decisionH1r.AR.r])^2)/batch2.size[n+1]
          
          # for left sided check
          xbar.diff1r_n.l = cumsum11r_n[not.reached.decisionH1r.AR.l]/batch1.size[n+1] -
            cumsum21r_n[not.reached.decisionH1r.AR.l]/batch2.size[n+1]
          divisor.pooled.sd1r_n.sq.l = 
            cumSS11r_n[not.reached.decisionH1r.AR.l] -
            ((cumsum11r_n[not.reached.decisionH1r.AR.l])^2)/batch1.size[n+1] +
            cumSS21r_n[not.reached.decisionH1r.AR.l] - 
            ((cumsum21r_n[not.reached.decisionH1r.AR.l])^2)/batch2.size[n+1]
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1r_n.r[not.reached.decisionH1r.AR.r] = 
            ((1 + ((xbar.diff1r_n.r - theta0)^2)/
                (divisor.pooled.sd1r_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff1r_n.r - 
                        (theta0 + t.alpha*
                           sqrt((divisor.pooled.sd1r_n.sq.r/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd1r_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          # for left sided check
          LR1r_n.l[not.reached.decisionH1r.AR.l] = 
            ((1 + ((xbar.diff1r_n.l - theta0)^2)/
                (divisor.pooled.sd1r_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff1r_n.l - 
                        (theta0 - t.alpha*
                           sqrt((divisor.pooled.sd1r_n.sq.l/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd1r_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]<=Accept.threshold
          RejectedH0.underH1r_n.AR.r = LR1r_n.r[not.reached.decisionH1r.AR.r]>=Reject.threshold
          reached.decisionH1r_n.AR.r = AcceptedH0.underH1r_n.AR.r|RejectedH0.underH1r_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.r)){
            
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[AcceptedH0.underH1r_n.AR.r]] = 'A'
            decision.underH1r.AR.r[not.reached.decisionH1r.AR.r[RejectedH0.underH1r_n.AR.r]] = 'R'
            N11r.AR.r[not.reached.decisionH1r.AR.r[reached.decisionH1r_n.AR.r]] = batch1.size[n+1]
            N21r.AR.r[not.reached.decisionH1r.AR.r[reached.decisionH1r_n.AR.r]] = batch2.size[n+1]
            not.reached.decisionH1r.AR.r = not.reached.decisionH1r.AR.r[!reached.decisionH1r_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]<=Accept.threshold
          RejectedH0.underH1r_n.AR.l = LR1r_n.l[not.reached.decisionH1r.AR.l]>=Reject.threshold
          reached.decisionH1r_n.AR.l = AcceptedH0.underH1r_n.AR.l|RejectedH0.underH1r_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1r_n.AR.l)){
            
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[AcceptedH0.underH1r_n.AR.l]] = 'A'
            decision.underH1r.AR.l[not.reached.decisionH1r.AR.l[RejectedH0.underH1r_n.AR.l]] = 'R'
            N11r.AR.l[not.reached.decisionH1r.AR.l[reached.decisionH1r_n.AR.l]] = batch1.size[n+1]
            N21r.AR.l[not.reached.decisionH1r.AR.l[reached.decisionH1r_n.AR.l]] = batch2.size[n+1]
            not.reached.decisionH1r.AR.l = not.reached.decisionH1r.AR.l[!reached.decisionH1r_n.AR.l]
          }
          
          not.reached.decisionH1r.AR = union(not.reached.decisionH1r.AR.r,
                                             not.reached.decisionH1r.AR.l)
        }
        
        
        ## under left-sided H1
        if(length(not.reached.decisionH1l.AR)>0){
          
          ## observations at step n
          # Group 1
          if(length(not.reached.decisionH1l.AR)>1){
            
            obs11l_n = mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                              FUN = function(X){
                                
                                rnorm(length(not.reached.decisionH1l.AR), theta1$left/2, 1)
                              })
            
          }else{
            
            obs11l_n = matrix(mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                                     FUN = function(X){
                                       
                                       rnorm(length(not.reached.decisionH1l.AR), theta1$left/2, 1)
                                       
                                     }), nrow = 1, ncol = batch1.size[n+1]-batch1.size[n], 
                              byrow = T)
          }
          
          # Group 2
          if(length(not.reached.decisionH1l.AR)>1){
            
            obs21l_n = mapply(X = 1:(batch2.size[n+1]-batch2.size[n]), 
                              FUN = function(X){
                                
                                rnorm(length(not.reached.decisionH1l.AR), -theta1$left/2, 1)
                              })
            
          }else{
            
            obs21l_n = matrix(mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                                     FUN = function(X){
                                       
                                       rnorm(length(not.reached.decisionH1l.AR), -theta1$left/2, 1)
                                       
                                     }), nrow = 1, ncol = batch1.size[n+1]-batch1.size[n], 
                              byrow = T)
          }
          
          ## sum of observations until step n
          # Group 1
          cumsum11l_n[not.reached.decisionH1l.AR] = 
            cumsum11l_n[not.reached.decisionH1l.AR] + rowSums(obs11l_n)
          # Group 2
          cumsum21l_n[not.reached.decisionH1l.AR] = 
            cumsum21l_n[not.reached.decisionH1l.AR] + rowSums(obs21l_n)
          
          ## sum of squares of observations until step n
          # Group 1
          cumSS11l_n[not.reached.decisionH1l.AR] = 
            cumSS11l_n[not.reached.decisionH1l.AR] + rowSums(obs11l_n^2)
          # Group 2
          cumSS21l_n[not.reached.decisionH1l.AR] = 
            cumSS21l_n[not.reached.decisionH1l.AR] + rowSums(obs21l_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          # for right sided check
          xbar.diff1l_n.r = cumsum11l_n[not.reached.decisionH1l.AR.r]/batch1.size[n+1] -
            cumsum21l_n[not.reached.decisionH1l.AR.r]/batch2.size[n+1]
          divisor.pooled.sd1l_n.sq.r = 
            cumSS11l_n[not.reached.decisionH1l.AR.r] - 
            ((cumsum11l_n[not.reached.decisionH1l.AR.r])^2)/batch1.size[n+1] +
            cumSS21l_n[not.reached.decisionH1l.AR.r] - 
            ((cumsum21l_n[not.reached.decisionH1l.AR.r])^2)/batch2.size[n+1]
          
          # for left sided check
          xbar.diff1l_n.l = cumsum11l_n[not.reached.decisionH1l.AR.l]/batch1.size[n+1] -
            cumsum21l_n[not.reached.decisionH1l.AR.l]/batch2.size[n+1]
          divisor.pooled.sd1l_n.sq.l = 
            cumSS11l_n[not.reached.decisionH1l.AR.l] -
            ((cumsum11l_n[not.reached.decisionH1l.AR.l])^2)/batch1.size[n+1] +
            cumSS21l_n[not.reached.decisionH1l.AR.l] - 
            ((cumsum21l_n[not.reached.decisionH1l.AR.l])^2)/batch2.size[n+1]
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1l_n.r[not.reached.decisionH1l.AR.r] = 
            ((1 + ((xbar.diff1l_n.r - theta0)^2)/
                (divisor.pooled.sd1l_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff1l_n.r - 
                        (theta0 + t.alpha*
                           sqrt((divisor.pooled.sd1l_n.sq.r/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd1l_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          # for left sided check
          LR1l_n.l[not.reached.decisionH1l.AR.l] = 
            ((1 + ((xbar.diff1l_n.l - theta0)^2)/
                (divisor.pooled.sd1l_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff1l_n.l - 
                        (theta0 - t.alpha*
                           sqrt((divisor.pooled.sd1l_n.sq.l/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd1l_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]<=Accept.threshold
          RejectedH0.underH1l_n.AR.r = LR1l_n.r[not.reached.decisionH1l.AR.r]>=Reject.threshold
          reached.decisionH1l_n.AR.r = AcceptedH0.underH1l_n.AR.r|RejectedH0.underH1l_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.r)){
            
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[AcceptedH0.underH1l_n.AR.r]] = 'A'
            decision.underH1l.AR.r[not.reached.decisionH1l.AR.r[RejectedH0.underH1l_n.AR.r]] = 'R'
            N11l.AR.r[not.reached.decisionH1l.AR.r[reached.decisionH1l_n.AR.r]] = batch1.size[n+1]
            N21l.AR.r[not.reached.decisionH1l.AR.r[reached.decisionH1l_n.AR.r]] = batch2.size[n+1]
            not.reached.decisionH1l.AR.r = not.reached.decisionH1l.AR.r[!reached.decisionH1l_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]<=Accept.threshold
          RejectedH0.underH1l_n.AR.l = LR1l_n.l[not.reached.decisionH1l.AR.l]>=Reject.threshold
          reached.decisionH1l_n.AR.l = AcceptedH0.underH1l_n.AR.l|RejectedH0.underH1l_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1l_n.AR.l)){
            
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[AcceptedH0.underH1l_n.AR.l]] = 'A'
            decision.underH1l.AR.l[not.reached.decisionH1l.AR.l[RejectedH0.underH1l_n.AR.l]] = 'R'
            N11l.AR.l[not.reached.decisionH1l.AR.l[reached.decisionH1l_n.AR.l]] = batch1.size[n+1]
            N21l.AR.l[not.reached.decisionH1l.AR.l[reached.decisionH1l_n.AR.l]] = batch2.size[n+1]
            not.reached.decisionH1l.AR.l = not.reached.decisionH1l.AR.l[!reached.decisionH1l_n.AR.l]
          }
          
          not.reached.decisionH1l.AR = union(not.reached.decisionH1l.AR.r,
                                             not.reached.decisionH1l.AR.l)
        }
        
        setTxtProgressBar(pb, n)
      }
      
      
      ### both-sided checking
      ## under H0
      # accepted or rejected ones
      accepted.by.both0 = intersect(which(decision.underH0.AR.r=='A'),
                                    which(decision.underH0.AR.l=='A'))
      onlyrejected.by.right0 = intersect(which(decision.underH0.AR.r=='R'),
                                         which(decision.underH0.AR.l!='R'))
      onlyrejected.by.left0 = intersect(which(decision.underH0.AR.r!='R'),
                                        which(decision.underH0.AR.l=='R'))
      rejected.by.both0 = intersect(which(decision.underH0.AR.r=='R'),
                                    which(decision.underH0.AR.l=='R'))
      
      ## sample sizes required
      # Group 1
      N10.AR[accepted.by.both0] = pmax(N10.AR.r[accepted.by.both0],
                                       N10.AR.l[accepted.by.both0])
      N10.AR[onlyrejected.by.right0] = N10.AR.r[onlyrejected.by.right0]
      N10.AR[onlyrejected.by.left0] = N10.AR.l[onlyrejected.by.left0]
      N10.AR[rejected.by.both0] = pmin(N10.AR.r[rejected.by.both0],
                                       N10.AR.l[rejected.by.both0])
      
      # Group 2
      N20.AR[accepted.by.both0] = pmax(N20.AR.r[accepted.by.both0],
                                       N20.AR.l[accepted.by.both0])
      N20.AR[onlyrejected.by.right0] = N20.AR.r[onlyrejected.by.right0]
      N20.AR[onlyrejected.by.left0] = N20.AR.l[onlyrejected.by.left0]
      N20.AR[rejected.by.both0] = pmin(N20.AR.r[rejected.by.both0],
                                       N20.AR.l[rejected.by.both0])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right0 = intersect(which(decision.underH0.AR.r=='A'),
                                         which(is.na(decision.underH0.AR.l)))
      onlyaccepted.by.left0 = intersect(which(is.na(decision.underH0.AR.r)),
                                        which(decision.underH0.AR.l=='A'))
      both.inconclusive0 = intersect(which(is.na(decision.underH0.AR.r)),
                                     which(is.na(decision.underH0.AR.l)))
      all.inconclusive0 = c(onlyaccepted.by.right0, onlyaccepted.by.left0,
                            both.inconclusive0)
      nNot.reached.decisionH0.AR = length(all.inconclusive0)
      
      # Type I error probability
      type1.error.AR[c(onlyrejected.by.right0, onlyrejected.by.left0,
                       rejected.by.both0)] = T
      
      
      ## under right-sided H1
      # accepted or rejected ones
      accepted.by.both1r = intersect(which(decision.underH1r.AR.r=='A'),
                                     which(decision.underH1r.AR.l=='A'))
      onlyrejected.by.right1r = intersect(which(decision.underH1r.AR.r=='R'),
                                          which(decision.underH1r.AR.l!='R'))
      onlyrejected.by.left1r = intersect(which(decision.underH1r.AR.r!='R'),
                                         which(decision.underH1r.AR.l=='R'))
      rejected.by.both1r = intersect(which(decision.underH1r.AR.r=='R'),
                                     which(decision.underH1r.AR.l=='R'))
      
      ## sample sizes required
      # Group 1
      N11r.AR[accepted.by.both1r] = pmax(N11r.AR.r[accepted.by.both1r],
                                         N11r.AR.l[accepted.by.both1r])
      N11r.AR[onlyrejected.by.right1r] = N11r.AR.r[onlyrejected.by.right1r]
      N11r.AR[onlyrejected.by.left1r] = N11r.AR.l[onlyrejected.by.left1r]
      N11r.AR[rejected.by.both1r] = pmin(N11r.AR.r[rejected.by.both1r],
                                         N11r.AR.l[rejected.by.both1r])
      
      # Group 2
      N21r.AR[accepted.by.both1r] = pmax(N21r.AR.r[accepted.by.both1r],
                                         N21r.AR.l[accepted.by.both1r])
      N21r.AR[onlyrejected.by.right1r] = N21r.AR.r[onlyrejected.by.right1r]
      N21r.AR[onlyrejected.by.left1r] = N21r.AR.l[onlyrejected.by.left1r]
      N21r.AR[rejected.by.both1r] = pmin(N21r.AR.r[rejected.by.both1r],
                                         N21r.AR.l[rejected.by.both1r])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1r = intersect(which(decision.underH1r.AR.r=='A'),
                                          which(is.na(decision.underH1r.AR.l)))
      onlyaccepted.by.left1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                         which(decision.underH1r.AR.l=='A'))
      both.inconclusive1r = intersect(which(is.na(decision.underH1r.AR.r)),
                                      which(is.na(decision.underH1r.AR.l)))
      all.inconclusive1r = c(onlyaccepted.by.right1r, onlyaccepted.by.left1r,
                             both.inconclusive1r)
      nNot.reached.decisionH1r.AR = length(all.inconclusive1r)
      
      # Type I error probability
      PowerH1r.AR[c(onlyrejected.by.right1r, onlyrejected.by.left1r,
                    rejected.by.both1r)] = T
      
      
      ## under left-sided H1
      # accepted or rejected ones
      accepted.by.both1l = intersect(which(decision.underH1l.AR.r=='A'),
                                     which(decision.underH1l.AR.l=='A'))
      onlyrejected.by.right1l = intersect(which(decision.underH1l.AR.r=='R'),
                                          which(decision.underH1l.AR.l!='R'))
      onlyrejected.by.left1l = intersect(which(decision.underH1l.AR.r!='R'),
                                         which(decision.underH1l.AR.l=='R'))
      rejected.by.both1l = intersect(which(decision.underH1l.AR.r=='R'),
                                     which(decision.underH1l.AR.l=='R'))
      
      ## sample sizes required
      # Group 1
      N11l.AR[accepted.by.both1l] = pmax(N11l.AR.r[accepted.by.both1l],
                                         N11l.AR.l[accepted.by.both1l])
      N11l.AR[onlyrejected.by.right1l] = N11l.AR.r[onlyrejected.by.right1l]
      N11l.AR[onlyrejected.by.left1l] = N11l.AR.l[onlyrejected.by.left1l]
      N11l.AR[rejected.by.both1l] = pmin(N11l.AR.r[rejected.by.both1l],
                                         N11l.AR.l[rejected.by.both1l])
      
      # Group 2
      N21l.AR[accepted.by.both1l] = pmax(N21l.AR.r[accepted.by.both1l],
                                         N21l.AR.l[accepted.by.both1l])
      N21l.AR[onlyrejected.by.right1l] = N21l.AR.r[onlyrejected.by.right1l]
      N21l.AR[onlyrejected.by.left1l] = N21l.AR.l[onlyrejected.by.left1l]
      N21l.AR[rejected.by.both1l] = pmin(N21l.AR.r[rejected.by.both1l],
                                         N21l.AR.l[rejected.by.both1l])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1l = intersect(which(decision.underH1l.AR.r=='A'),
                                          which(is.na(decision.underH1l.AR.l)))
      onlyaccepted.by.left1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                         which(decision.underH1l.AR.l=='A'))
      both.inconclusive1l = intersect(which(is.na(decision.underH1l.AR.r)),
                                      which(is.na(decision.underH1l.AR.l)))
      all.inconclusive1l = c(onlyaccepted.by.right1l, onlyaccepted.by.left1l,
                             both.inconclusive1l)
      nNot.reached.decisionH1l.AR = length(all.inconclusive1l)
      
      # Type I error probability
      PowerH1l.AR[c(onlyrejected.by.right1l, onlyrejected.by.left1l,
                    rejected.by.both1l)] = T
      
      
      ## determining termination threshold
      ## H0 is rejected if LR or (BF) is >= termination threshold
      type1.error.spent.AR = mean(type1.error.AR) # type 1 error already spent
      if(nNot.reached.decisionH0.AR==0){
        
        nDecimal.accuracy = 2
        termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
          (10^nDecimal.accuracy)
        
        actual.type1.error.AR = type1.error.spent.AR
        
      }else{
        
        term.thresh.possible.choices =
          c(LR0_n.r[onlyaccepted.by.left0],
            LR0_n.l[onlyaccepted.by.right0],
            pmin(LR0_n.r[both.inconclusive0], LR0_n.l[both.inconclusive0]))
        
        type1.error.max.AR = type1.error.spent.AR + nNot.reached.decisionH0.AR/nReplicate
        if(type1.error.spent.AR>Type1.target){
          
          max.LR0_n = max(term.thresh.possible.choices)
          nDecimal.accuracy = ceiling(-log10(min(0.01, Reject.threshold - max.LR0_n)))
          termination.threshold.AR = (floor(max.LR0_n*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.spent.AR
          
        }else if(type1.error.max.AR<=Type1.target){
          
          nDecimal.accuracy = ceiling(-log10(min(0.01, min(term.thresh.possible.choices) -
                                                   Accept.threshold)))
          termination.threshold.AR = (floor(Accept.threshold*(10^nDecimal.accuracy)) + 1)/
            (10^nDecimal.accuracy)
          
          actual.type1.error.AR = type1.error.max.AR
          
        }else{
          
          uniqLR0.not.reached.decisionH0.inc.AR = sort(unique(term.thresh.possible.choices))
          cumRejFreq_not.reached.decisionH0.AR = cumsum(rev(as.numeric(table(term.thresh.possible.choices))))
          nNewRejects.AR = floor(nReplicate*(Type1.target - type1.error.spent.AR)) # max new rejects
          
          if(cumRejFreq_not.reached.decisionH0.AR[1]>nNewRejects.AR){
            
            nDecimal.accuracy =
              ceiling(-log10(min(0.01, Reject.threshold -
                                   uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)])))
            termination.threshold.AR =
              (floor(uniqLR0.not.reached.decisionH0.inc.AR[length(uniqLR0.not.reached.decisionH0.inc.AR)]*
                       (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR
            
          }else{
            
            opt.indx.AR = max(which(cumRejFreq_not.reached.decisionH0.AR<=nNewRejects.AR))
            min.rej.indx.AR = length(uniqLR0.not.reached.decisionH0.inc.AR) - (opt.indx.AR - 1)
            
            nDecimal.accuracy = ceiling(-log10(min(0.01,
                                                   uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR] -
                                                     uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1])))
            termination.threshold.AR = (floor(uniqLR0.not.reached.decisionH0.inc.AR[min.rej.indx.AR-1]*
                                                (10^nDecimal.accuracy)) + 1)/(10^nDecimal.accuracy)
            
            actual.type1.error.AR = type1.error.spent.AR +
              cumRejFreq_not.reached.decisionH0.AR[opt.indx.AR]/nReplicate
          }
        }
      }
      
      ## attained Type II error probability
      # right-sided H1
      actual.PowerH1r.AR.r = mean(PowerH1r.AR) +
        sum(c(LR1r_n.r[onlyaccepted.by.left1r],
              LR1r_n.l[onlyaccepted.by.right1r],
              pmax(LR1r_n.r[both.inconclusive1r], LR1r_n.l[both.inconclusive1r]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1r.AR = 1 - actual.PowerH1r.AR.r
      
      # left-sided H1
      actual.PowerH1l.AR.r = mean(PowerH1l.AR) +
        sum(c(LR1l_n.r[onlyaccepted.by.left1l],
              LR1l_n.l[onlyaccepted.by.right1l],
              pmax(LR1l_n.r[both.inconclusive1l], LR1l_n.l[both.inconclusive1l]))>=
              termination.threshold.AR)/nReplicate
      actual.type2.errorH1l.AR = 1 - actual.PowerH1l.AR.r
      
      ## Expected sample sizes
      # Group 1
      EN10 = mean(N10.AR)     # under H0
      EN11r = mean(N11r.AR)   # under right-sided H1
      EN11l = mean(N11l.AR)   # under left-sided H1
      
      # Group 2
      EN20 = mean(N20.AR)     # under H0
      EN21r = mean(N21r.AR)   # under right-sided H1
      EN21l = mean(N21l.AR)   # under left-sided H1
      
      # msg
      if(verbose==T){
        cat('\n')
        print('Done.')
        print("-------------------------------------------------------------------------")
        cat('\n\n')
        print("=========================================================================")
        print("Performance summary:")
        print("=========================================================================")
        print(paste("H0 acceptance threshold: ", round(Accept.threshold, 3)))
        print(paste("H0 rejection threshold: ", round(Reject.threshold, 3)))
        print(paste("Termination threshold: ", round(termination.threshold.AR, 3)))
        print(paste("Attained Type I error probability: ", round(actual.type1.error.AR, 4)))
        print(paste("Expected sample size under H0: Group 1 - ", round(EN10, 2), 
                    ', Group 2 - ', round(EN20, 2), sep = ''))
        print("Attained Type II error probability:")
        print(paste(" On the right: ", round(actual.type2.errorH1r.AR, 4)))
        print(paste(" On the left: ", round(actual.type2.errorH1l.AR, 4)))
        print("Expected sample size at the alternatives:")
        print(paste(" On the right: Group 1 - ", round(EN11r, 2), 
                    ', Group 2 - ', round(EN21r, 2), sep = ''))
        print(paste(" On the left: Group 1 - ", round(EN11l, 2), 
                    ', Group 2 - ', round(EN21l, 2), sep = ''))
        print("=========================================================================")
        cat('\n')
      }
      
      return(list("Type1.attained" = actual.type1.error.AR,
                  "Type2.attained" = c(actual.type2.errorH1r.AR, actual.type2.errorH1l.AR),
                  'N' = list('H0' = list('Group1' = N10.AR, 'Group2' = N20.AR),
                             'right' = list('Group1' = N11r.AR, 'Group2' = N21r.AR),
                             'left' = list('Group1' = N11l.AR, 'Group2' = N21l.AR)),
                  'EN' = list('H0' = list('Group1' = EN10, 'Group2' = EN20),
                              'right' = list('Group1' = EN11r, 'Group2' = EN21r),
                              'left' = list('Group1' = EN11l, 'Group2' = EN21l)),
                  "theta1" = theta1, "Type2.fixed.design" = Type2.target,
                  "Reject.threshold" = Reject.threshold, "Accept.threshold" = Accept.threshold,
                  "termination.threshold" = termination.threshold.AR,
                  'test.type' = 'twoT', 'side' = side, 'theta0' = theta0, 
                  'N1.max' = N1.max, 'N2.max' = N2.max,
                  'Type1.target' = Type1.target, 'Type2.target' = Type2.target,
                  'batch1.size' = diff(batch1.size), 'batch2.size' = diff(batch2.size),
                  'nAnalyses' = nAnalyses, 'nReplicate' = nReplicate, 'seed' = seed))
    }
  }
}

#### designing the MSPRT combined for all ####
design.MSPRT = function(test.type, side = 'right', theta0, theta1 = T,
                        Type1.target = .005, Type2.target = .2,
                        N.max, N1.max, N2.max,
                        sigma = 1, sigma1 = 1, sigma2 = 1,
                        batch.size, batch1.size, batch2.size,
                        nReplicate = 1e+6, verbose = T, seed = 1){
  
  if(test.type=='oneProp'){
    
    ## ignoring batch1.seq & batch2.seq
    if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
    if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
    
    ## ignoring N1.max & N2.max
    if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
    if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
    
    ## batch sizes and N.max
    if(missing(batch.size)){
      
      if(missing(N.max)){
        
        return("Either 'batch.size' or 'N.max' needs to be specified")
        
      }else{batch.size = rep(1, N.max)}
      
    }else{
      
      if(missing(N.max)){
        
        N.max = sum(batch.size)
        
      }else{
        
        if(sum(batch.size)!=N.max) return("Sum of batch sizes should add up to N.max")
      }
    }
    
    if(missing(theta0)) theta0 = 0.5
    
    return(design.MSPRT_oneProp(side = side, theta0 = theta0, theta1 = theta1,
                                Type1.target = Type1.target,
                                Type2.target = Type2.target,
                                N.max = N.max, batch.size = batch.size,
                                nReplicate = nReplicate, 
                                verbose = verbose, seed = seed))
    
  }else if(test.type=='oneZ'){
    
    ## ignoring batch1.seq & batch2.seq
    if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
    if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
    
    ## ignoring N1.max & N2.max
    if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
    if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
    
    ## batch sizes and N.max
    if(missing(batch.size)){
      
      if(missing(N.max)){
        
        return("Either 'batch.size' or 'N.max' needs to be specified")
        
      }else{batch.size = rep(1, N.max)}
      
    }else{
      
      if(missing(N.max)){
        
        N.max = sum(batch.size)
        
      }else{
        
        if(sum(batch.size)!=N.max) return("Sum of batch sizes should add up to N.max")
      }
    }
    
    if(missing(theta0)) theta0 = 0
    
    return(design.MSPRT_oneZ(side = side, theta0 = theta0, theta1 = theta1,
                             Type1.target = Type1.target,
                             Type2.target = Type2.target,
                             N.max = N.max, sigma = sigma,
                             batch.size = batch.size,
                             nReplicate = nReplicate, 
                             verbose = verbose, seed = seed))
    
  }else if(test.type=='oneT'){
    
    ## ignoring batch1.seq & batch2.seq
    if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
    if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
    
    ## ignoring N1.max & N2.max
    if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
    if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
    
    ## batch sizes and N.max
    if(missing(batch.size)){
      
      if(missing(N.max)){
        
        return("Either 'batch.size' or 'N.max' needs to be specified")
        
      }else{batch.size = c(2, rep(1, N.max-2))}
      
    }else{
      
      if(batch.size[1]<2){
        
        return("First batch size should be at least 2")
        
      }else{
        
        if(missing(N.max)){
          
          N.max = sum(batch.size)
          
        }else{
          
          if(sum(batch.size)!=N.max) return("Sum of batch.size should add up to N.max")
        }
      }
    }
    
    if(missing(theta0)) theta0 = 0
    
    return(design.MSPRT_oneT(side = side, theta0 = theta0, theta1 = theta1,
                             Type1.target = Type1.target,
                             Type2.target = Type2.target,
                             N.max = N.max, batch.size = batch.size,
                             nReplicate = nReplicate, 
                             verbose = verbose, seed = seed))
    
  }else if(test.type=='twoZ'){
    
    ## ignoring batch.size
    if(!missing(batch.size)) print("'batch.size' is ignored. Not required in two-sample tests.")
    
    ## ignoring N.max
    if(!missing(N.max)) print("'N.max' is ignored. Not required in two-sample tests.")
    
    ## checking if length(batch1.size) and length(batch2.size) are equal
    if((!missing(batch1.size)) && (!missing(batch2.size)) &&
       (length(batch1.size)!=length(batch2.size))) return("Lenghts of batch1.size and batch2.size should be same")
    
    ## batch sizes and N for group 1
    if(missing(batch1.size)){
      
      if(missing(N1.max)){
        
        return(print("Either 'batch1.size' or 'N1.max' needs to be specified"))
        
      }else{batch1.size = rep(1, N1.max)}
      
    }else{
      
      if(missing(N1.max)){
        
        N1.max = sum(batch1.size)
        
      }else{
        
        if(sum(batch1.size)!=N1.max) return(print("Sum of batch1.size should add up to N1.max"))
      }
    }
    
    ## batch sizes and N for group 2
    if(missing(batch2.size)){
      
      if(missing(N2.max)){
        
        return(print("Either 'batch2.size' or 'N2.max' needs to be specified"))
        
      }else{batch2.size = rep(1, N2.max)}
      
    }else{
      
      if(missing(N2.max)){
        
        N2.max = sum(batch2.size)
        
      }else{
        
        if(sum(batch2.size)!=N1.max) return(print("Sum of batch2.size should add up to N2.max"))
      }
    }
    
    if(missing(theta0)) theta0 = 0
    
    return(design.MSPRT_twoZ(side = side, theta0 = theta0, theta1 = theta1,
                             Type1.target = Type1.target, Type2.target = Type2.target,
                             N1.max = N1.max, N2.max = N2.max,
                             sigma1 = sigma1, sigma2 = sigma2,
                             batch1.size = batch1.size, batch2.size = batch2.size,
                             nReplicate = nReplicate, verbose = verbose, seed = seed))
    
  }else if(test.type=='twoT'){
    
    ## ignoring batch.size
    if(!missing(batch.size)) print("'batch.size' is ignored. Not required in two-sample tests.")
    
    ## ignoring N.max
    if(!missing(N.max)) print("'N.max' is ignored. Not required in two-sample tests.")
    
    ## checking if length(batch1.size) and length(batch2.size) are equal
    if((!missing(batch1.size)) && (!missing(batch2.size)) &&
       (length(batch1.size)!=length(batch2.size))) return("Lenghts of batch1.size and batch2.size should be same")
    
    ## batch sizes and N for group 1
    if(missing(batch1.size)){
      
      if(missing(N1.max)){
        
        return("Either 'batch1.size' or 'N1.max' needs to be specified")
        
      }else{batch1.size = c(2, rep(1, N1.max-2))}
      
    }else{
      
      if(batch1.size[1]<2){
        
        return("First batch size in Group 1 should be at least 2")
        
      }else{
        
        if(missing(N1.max)){
          
          N1.max = sum(batch1.size)
          
        }else{
          
          if(sum(batch1.size)!=N1.max) return("Sum of batch1.size should add up to N1.max")
        }
      }
    }
    
    ## batch sizes and N for group 2
    if(missing(batch2.size)){
      
      if(missing(N2.max)){
        
        return("Either 'batch2.size' or 'N2.max' needs to be specified")
        
      }else{batch2.size = c(2, rep(1, N2.max-2))}
      
    }else{
      
      if(batch2.size[1]<2){
        
        return("First batch size in Group 2 should be at least 2")
        
      }else{
        
        if(missing(N2.max)){
          
          N2.max = sum(batch2.size)
          
        }else{
          
          if(sum(batch2.size)!=N2.max) return("Sum of batch2.size should add up to N2.max")
        }
      }
    }
    
    if(missing(theta0)) theta0 = 0
    
    return(design.MSPRT_twoT(side = side, theta0 = theta0, theta1 = theta1,
                             Type1.target = Type1.target, Type2.target = Type2.target,
                             N1.max = N1.max, N2.max = N2.max,
                             batch1.size = batch1.size, batch2.size = batch2.size,
                             nReplicate = nReplicate, verbose = verbose, seed = seed))
    
  }
}


################################### OC and ASN of the MSPRT ###################################

#### one-sample proportion test ####
OCandASN.MSPRT_oneProp = function(theta, design.MSPRT.object, 
                                  termination.threshold,
                                  side = 'right', theta0 = 0.5, 
                                  Type1.target =.005, Type2.target = .2,
                                  N.max, batch.size,
                                  nReplicate = 1e+6, nCore = detectCores() - 1,
                                  verbose = T, seed = 1){
  
  # side
  if(!missing(design.MSPRT.object)) side = design.MSPRT.object$side
  
  if(side!='both'){
    
    #################### one-sample proportion (right/left sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch.size = design.MSPRT.object$batch.size
      N.max = design.MSPRT.object$N.max
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      termination.threshold = design.MSPRT.object$termination.threshold
      UMPBT = design.MSPRT.object$UMPBT
      nAnalyses = design.MSPRT.object$nAnalyses
      nReplicate = design.MSPRT.object$nReplicate
      
      # msg
      if(verbose){
        
        if(any(batch.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the group sequential MSPRT for a one-sample proportion test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the sequential MSPRT for a one-sample proportion test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", design.MSPRT.object$termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(UMPBT$theta[1], 3), " & ",
                    round(UMPBT$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$mix.prob[1], 3), " & ", 1 - round(UMPBT$mix.prob[1], 3), sep = ''))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch.size = c(0, cumsum(batch.size))
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
      if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
      
      ## ignoring N1.max & N2.max
      if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
      if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
      
      ## batch sizes and N.max
      if(missing(batch.size)){
        
        if(missing(N.max)){
          
          return("Either 'batch.size' or 'N.max' needs to be specified")
          
        }else{batch.size = rep(1, N.max)}
        
      }else{
        
        if(missing(N.max)){
          
          N.max = sum(batch.size)
          
        }else{
          
          if(sum(batch.size)!=N.max) return("Sum of batch sizes should add up to N.max")
        }
      }
      
      nAnalyses = length(batch.size)
      
      ######################## UMPBT alternative ########################
      UMPBT = UMPBT.alt(test.type = 'oneProp', side = side, theta0 = theta0,
                        N = N.max, Type1 = Type1.target)
      
      
      # msg
      if(verbose){
        
        if(any(batch.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the group sequential MSPRT for a one-sample proportion test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the sequential MSPRT for a one-sample proportion test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", design.MSPRT.object$termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(UMPBT$theta[1], 3), " & ",
                    round(UMPBT$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$mix.prob[1], 3), " & ", 1 - round(UMPBT$mix.prob[1], 3), sep = ''))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch.size = c(0, cumsum(batch.size))
    }
    
    
    registerDoParallel(cores = nCore)
    out.OCandASN = foreach(theta1 = theta, .combine = 'rbind') %dopar% {
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # required storages
      cumsum1_n = LR1_n = numeric(nReplicate)
      type2.error.AR = rep(F, nReplicate)
      N1.AR = rep(N.max, nReplicate)
      not.reached.decisionH1.AR = 1:nReplicate
      
      set.seed(seed)
      for(n in 1:nAnalyses){
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          # sum of observations at step n
          sum1_n = rbinom(length(not.reached.decisionH1.AR),
                          batch.size[n+1]-batch.size[n], theta1)
          
          # sum of observations until step n
          cumsum1_n[not.reached.decisionH1.AR] = 
            cumsum1_n[not.reached.decisionH1.AR] + sum1_n
          
          # likelihood ratio of observations until step n
          LR1_n[not.reached.decisionH1.AR] = 
            UMPBT$mix.prob[1]*(((1 - UMPBT$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$theta[1])))^cumsum1_n[not.reached.decisionH1.AR] +
            (1 - UMPBT$mix.prob[2])*(((1 - UMPBT$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$theta[2])))^cumsum1_n[not.reached.decisionH1.AR]
          
          # comparing with the thresholds
          AcceptedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]<=Accept.threshold)
          RejectedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]>=Reject.threshold)
          reached.decisionH1_n.AR = union(AcceptedH0.underH1_n.AR, RejectedH0.underH1_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH1_n.AR)>0){
            
            N1.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch.size[n+1]
            type2.error.AR[not.reached.decisionH1.AR[AcceptedH0.underH1_n.AR]] = T
            not.reached.decisionH1.AR = not.reached.decisionH1.AR[-reached.decisionH1_n.AR]
          }
        }
      }
      
      # attained Type II error probability
      actual.type2.error.AR = mean(type2.error.AR) +
        sum(LR1_n[not.reached.decisionH1.AR]<termination.threshold)/nReplicate
      
      # Expected sample sizes
      EN1 = mean(N1.AR)
      
      c(theta1, actual.type2.error.AR, EN1)
    }
    
    if(length(theta)==1) out.OCandASN = matrix(data = out.OCandASN, nrow = 1,
                                               ncol = 3, byrow = T)
    
    out.OCandASN = as.data.frame(out.OCandASN)
    colnames(out.OCandASN) = c('theta', 'acceptH0.prob', 'EN')
    
    # msg
    if(verbose==T){
      cat('\n')
      print('Done.')
      print("-------------------------------------------------------------------------")
      cat('\n\n')
      print("=========================================================================")
      print("Performance summary:")
      print("=========================================================================")
      print(paste('Parameter value(s): ', paste(round(theta, 3), collapse = ', '), sep = ''))
      print(paste('Probability of accepting H0: ',
                  paste(round(out.OCandASN$acceptH0.prob, 3), collapse = ', '), sep = ''))
      print(paste('Expected sample size: ',
                  paste(round(out.OCandASN$EN, 2), collapse = ', '), sep = ''))
      print("=========================================================================")
      cat('\n')
    }
    
    return(out.OCandASN)
    
    # end one-sided oneProp
    
  }else{
    
    #################### one-sample proportion (both sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch.size = design.MSPRT.object$batch.size
      N.max = design.MSPRT.object$N.max
      nAnalyses = design.MSPRT.object$nAnalyses
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      termination.threshold = design.MSPRT.object$termination.threshold
      nReplicate = design.MSPRT.object$nReplicate
      UMPBT = design.MSPRT.object$UMPBT
      
      # msg
      if(verbose){
        
        if(any(batch.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the group sequential MSPRT for a one-sample proportion test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the sequential MSPRT for a one-sample proportion test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", design.MSPRT.object$termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right:, ', round(UMPBT$right$theta[1], 3), " & ",
                    round(UMPBT$right$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$right$mix.prob[1], 3), " & ", 1 - round(UMPBT$right$mix.prob[1], 3),
                    sep = ""))
        print(paste(' On the left:, ', round(UMPBT$left$theta[1], 3), " & ",
                    round(UMPBT$left$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$left$mix.prob[1], 3), " & ", 1 - round(UMPBT$left$mix.prob[1], 3),
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch.size = c(0, cumsum(batch.size))
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
      if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
      
      ## ignoring N1.max & N2.max
      if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
      if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
      
      ## batch sizes and N.max
      if(missing(batch.size)){
        
        if(missing(N.max)){
          
          return("Either 'batch.size' or 'N.max' needs to be specified")
          
        }else{batch.size = rep(1, N.max)}
        
      }else{
        
        if(missing(N.max)){
          
          N.max = sum(batch.size)
          
        }else{
          
          if(sum(batch.size)!=N.max) return("Sum of batch sizes should add up to N.max")
        }
      }
      
      nAnalyses = length(batch.size)
      
      ######################## UMPBT alternative ########################
      UMPBT = list('right' = UMPBT.alt(test.type = 'oneProp', side = 'right', 
                                       theta0 = theta0, N = N.max, Type1 = Type1.target/2),
                   'left' = UMPBT.alt(test.type = 'oneProp', side = 'left',
                                      theta0 = theta0, N = N.max, Type1 = Type1.target/2))
      
      # msg
      if(verbose){
        
        if(any(batch.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the group sequential MSPRT for a one-sample proportion test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the sequential MSPRT for a one-sample proportion test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", design.MSPRT.object$termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(UMPBT$right$theta[1], 3), " & ",
                    round(UMPBT$right$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$right$mix.prob[1], 3), " & ", 1 - round(UMPBT$right$mix.prob[1], 3),
                    sep = ""))
        print(paste(' On the left: ', round(UMPBT$left$theta[1], 3), " & ",
                    round(UMPBT$left$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$left$mix.prob[1], 3), " & ", 1 - round(UMPBT$left$mix.prob[1], 3),
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch.size = c(0, cumsum(batch.size))
    }
    
    
    registerDoParallel(cores = nCore)
    out.OCandASN = foreach(theta1 = theta, .combine = 'rbind') %dopar% {
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # required storages
      cumsum1_n = LR1_n.r = LR1_n.l = numeric(nReplicate)
      PowerH1.AR = rep(F, nReplicate)
      N1.AR = N1.AR.r = N1.AR.l = rep(N.max, nReplicate)
      decision.underH1.AR.r = decision.underH1.AR.l = rep(NA, nReplicate)
      not.reached.decisionH1.AR = not.reached.decisionH1.AR.r = not.reached.decisionH1.AR.l =
        1:nReplicate
      
      set.seed(seed)
      pb = txtProgressBar(min = 1, max = nAnalyses, style = 3)
      for(n in 1:nAnalyses){
        
        ## under right-sided H1
        if(length(not.reached.decisionH1.AR)>0){
          
          # sum of observations at step n
          sum1_n = rbinom(length(not.reached.decisionH1.AR),
                          batch.size[n+1]-batch.size[n], theta1)
          
          # sum of observations until step n
          cumsum1_n[not.reached.decisionH1.AR] =
            cumsum1_n[not.reached.decisionH1.AR] + sum1_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1_n.r[not.reached.decisionH1.AR.r] = 
            UMPBT$right$mix.prob[1]*(((1 - UMPBT$right$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[1])))^cumsum1_n[not.reached.decisionH1.AR.r] +
            (1 - UMPBT$right$mix.prob[2])*(((1 - UMPBT$right$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[2])))^cumsum1_n[not.reached.decisionH1.AR.r]
          
          # for left sided check
          LR1_n.l[not.reached.decisionH1.AR.l] = 
            UMPBT$left$mix.prob[1]*(((1 - UMPBT$left$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[1])))^cumsum1_n[not.reached.decisionH1.AR.l] +
            (1 - UMPBT$left$mix.prob[2])*(((1 - UMPBT$left$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[2])))^cumsum1_n[not.reached.decisionH1.AR.l]
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1_n.AR.r = LR1_n.r[not.reached.decisionH1.AR.r]<=Accept.threshold
          RejectedH0.underH1_n.AR.r = LR1_n.r[not.reached.decisionH1.AR.r]>=Reject.threshold
          reached.decisionH1_n.AR.r = AcceptedH0.underH1_n.AR.r|RejectedH0.underH1_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1_n.AR.r)){
            
            decision.underH1.AR.r[not.reached.decisionH1.AR.r[AcceptedH0.underH1_n.AR.r]] = 'A'
            decision.underH1.AR.r[not.reached.decisionH1.AR.r[RejectedH0.underH1_n.AR.r]] = 'R'
            N1.AR.r[not.reached.decisionH1.AR.r[reached.decisionH1_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH1.AR.r = not.reached.decisionH1.AR.r[!reached.decisionH1_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1_n.AR.l = LR1_n.l[not.reached.decisionH1.AR.l]<=Accept.threshold
          RejectedH0.underH1_n.AR.l = LR1_n.l[not.reached.decisionH1.AR.l]>=Reject.threshold
          reached.decisionH1_n.AR.l = AcceptedH0.underH1_n.AR.l|RejectedH0.underH1_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1_n.AR.l)){
            
            decision.underH1.AR.l[not.reached.decisionH1.AR.l[AcceptedH0.underH1_n.AR.l]] = 'A'
            decision.underH1.AR.l[not.reached.decisionH1.AR.l[RejectedH0.underH1_n.AR.l]] = 'R'
            N1.AR.l[not.reached.decisionH1.AR.l[reached.decisionH1_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH1.AR.l = not.reached.decisionH1.AR.l[!reached.decisionH1_n.AR.l]
          }
          
          not.reached.decisionH1.AR = union(not.reached.decisionH1.AR.r,
                                            not.reached.decisionH1.AR.l)
        }
      }
      
      
      ### both-sided checking
      ## under H1
      # accepted or rejected ones
      accepted.by.both1 = intersect(which(decision.underH1.AR.r=='A'),
                                    which(decision.underH1.AR.l=='A'))
      onlyrejected.by.right1 = intersect(which(decision.underH1.AR.r=='R'),
                                         which(decision.underH1.AR.l!='R'))
      onlyrejected.by.left1 = intersect(which(decision.underH1.AR.r!='R'),
                                        which(decision.underH1.AR.l=='R'))
      rejected.by.both1 = intersect(which(decision.underH1.AR.r=='R'),
                                    which(decision.underH1.AR.l=='R'))
      
      # sample sizes required
      N1.AR[accepted.by.both1] = pmax(N1.AR.r[accepted.by.both1],
                                      N1.AR.l[accepted.by.both1])
      N1.AR[onlyrejected.by.right1] = N1.AR.r[onlyrejected.by.right1]
      N1.AR[onlyrejected.by.left1] = N1.AR.l[onlyrejected.by.left1]
      N1.AR[rejected.by.both1] = pmin(N1.AR.r[rejected.by.both1],
                                      N1.AR.l[rejected.by.both1])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1 = intersect(which(decision.underH1.AR.r=='A'),
                                         which(is.na(decision.underH1.AR.l)))
      onlyaccepted.by.left1 = intersect(which(is.na(decision.underH1.AR.r)),
                                        which(decision.underH1.AR.l=='A'))
      both.inconclusive1 = intersect(which(is.na(decision.underH1.AR.r)),
                                     which(is.na(decision.underH1.AR.l)))
      all.inconclusive1 = c(onlyaccepted.by.right1, onlyaccepted.by.left1,
                            both.inconclusive1)
      nNot.reached.decisionH1.AR = length(all.inconclusive1)
      
      # Type I error probability
      PowerH1.AR[c(onlyrejected.by.right1, onlyrejected.by.left1,
                   rejected.by.both1)] = T
      
      ## attained Type II error probability
      actual.PowerH1.AR.r = mean(PowerH1.AR) +
        sum(c(LR1_n.r[onlyaccepted.by.left1],
              LR1_n.l[onlyaccepted.by.right1],
              pmax(LR1_n.r[both.inconclusive1], LR1_n.l[both.inconclusive1]))>=
              termination.threshold)/nReplicate
      actual.type2.errorH1.AR = 1 - actual.PowerH1.AR.r
      
      ## Expected sample sizes
      EN1 = mean(N1.AR)   # under right-sided H1
      
      c(theta1, actual.type2.errorH1.AR, EN1)
    }
    
    if(length(theta)==1) out.OCandASN = matrix(data = out.OCandASN, nrow = 1,
                                               ncol = 3, byrow = T)
    
    out.OCandASN = as.data.frame(out.OCandASN)
    colnames(out.OCandASN) = c('theta', 'acceptH0.prob', 'EN')
    
    # msg
    if(verbose==T){
      cat('\n')
      print('Done.')
      print("-------------------------------------------------------------------------")
      cat('\n\n')
      print("=========================================================================")
      print("Performance summary:")
      print("=========================================================================")
      print(paste('Parameter value(s): ', paste(round(theta, 3), collapse = ', '), sep = ''))
      print(paste('Probability of accepting H0: ',
                  paste(round(out.OCandASN$acceptH0.prob, 3), collapse = ', '), sep = ''))
      print(paste('Expected sample size: ',
                  paste(round(out.OCandASN$EN, 2), collapse = ', '), sep = ''))
      print("=========================================================================")
      cat('\n')
    }
    
    return(out.OCandASN)
    
    # end both-sided oneProp
  }
}

#### one-sample z test ####
OCandASN.MSPRT_oneZ = function(theta, design.MSPRT.object, 
                               termination.threshold,
                               side = 'right', theta0 = 0, 
                               Type1.target =.005, Type2.target = .2,
                               N.max, sigma = 1, batch.size,
                               nReplicate = 1e+6, nCore = detectCores() - 1,
                               verbose = T, seed = 1){
  
  # side
  if(!missing(design.MSPRT.object)) side = design.MSPRT.object$side
  
  if(side!='both'){
    
    #################### one-sample z (right/left sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch.size = design.MSPRT.object$batch.size
      N.max = design.MSPRT.object$N.max
      nAnalyses = design.MSPRT.object$nAnalyses
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      sigma = design.MSPRT.object$sigma
      termination.threshold = design.MSPRT.object$termination.threshold
      nReplicate = design.MSPRT.object$nReplicate
      theta.UMPBT = design.MSPRT.object$theta.UMPBT
      
      # msg
      if(verbose){
        
        if(any(batch.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the group sequential MSPRT for a one-sample z test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the sequential MSPRT for a one-sample z test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Known standard deviation: ", sigma, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(theta.UMPBT, 3)))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch.size = c(0, cumsum(batch.size))
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
      if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
      
      ## ignoring N1.max & N2.max
      if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
      if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
      
      ## batch sizes and N.max
      if(missing(batch.size)){
        
        if(missing(N.max)){
          
          return("Either 'batch.size' or 'N.max' needs to be specified")
          
        }else{batch.size = rep(1, N.max)}
        
      }else{
        
        if(missing(N.max)){
          
          N.max = sum(batch.size)
          
        }else{
          
          if(sum(batch.size)!=N.max) return("Sum of batch sizes should add up to N.max")
        }
      }
      
      nAnalyses = length(batch.size)
      
      ######################## UMPBT alternative ########################
      theta.UMPBT = UMPBT.alt(test.type = 'oneZ', side = side, theta0 = theta0,
                              N = N.max, Type1 = Type1.target, sigma = sigma)
      
      
      # msg
      if(verbose){
        
        if(any(batch.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the group sequential MSPRT for a one-sample z test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the sequential MSPRT for a one-sample z test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Known standard deviation: ", sigma, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(theta.UMPBT, 3)))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch.size = c(0, cumsum(batch.size))
    }
    
    
    registerDoParallel(cores = nCore)
    out.OCandASN = foreach(theta1 = theta, .combine = 'rbind') %dopar% {
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # required storages
      cumsum1_n = LR1_n = numeric(nReplicate)
      type2.error.AR = rep(F, nReplicate)
      N1.AR = rep(N.max, nReplicate)
      not.reached.decisionH1.AR = 1:nReplicate
      
      set.seed(seed)
      for(n in 1:nAnalyses){
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          # sum of observations at step n
          sum1_n = rnorm(length(not.reached.decisionH1.AR),
                         (batch.size[n+1]-batch.size[n])*theta1,
                         sqrt(batch.size[n+1]-batch.size[n])*sigma)
          
          # sum of observations until step n
          cumsum1_n[not.reached.decisionH1.AR] = 
            cumsum1_n[not.reached.decisionH1.AR] + sum1_n
          
          # likelihood ratio of observations until step n
          LR1_n[not.reached.decisionH1.AR] = 
            exp((cumsum1_n[not.reached.decisionH1.AR]*(theta.UMPBT - theta0) - 
                   ((batch.size[n+1]*((theta.UMPBT^2) - (theta0^2)))/2))/(sigma^2))
          
          # comparing with the thresholds
          AcceptedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]<=Accept.threshold)
          RejectedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]>=Reject.threshold)
          reached.decisionH1_n.AR = union(AcceptedH0.underH1_n.AR, RejectedH0.underH1_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH1_n.AR)>0){
            
            N1.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch.size[n+1]
            type2.error.AR[not.reached.decisionH1.AR[AcceptedH0.underH1_n.AR]] = T
            not.reached.decisionH1.AR = not.reached.decisionH1.AR[-reached.decisionH1_n.AR]
          }
        }
      }
      
      # attained Type II error probability
      actual.type2.error.AR = mean(type2.error.AR) +
        sum(LR1_n[not.reached.decisionH1.AR]<termination.threshold)/nReplicate
      
      # Expected sample sizes
      EN1 = mean(N1.AR)
      
      c(theta1, actual.type2.error.AR, EN1)
    }
    
    if(length(theta)==1) out.OCandASN = matrix(data = out.OCandASN, nrow = 1,
                                               ncol = 3, byrow = T)
    
    out.OCandASN = as.data.frame(out.OCandASN)
    colnames(out.OCandASN) = c('theta', 'acceptH0.prob', 'EN')
    
    # msg
    if(verbose==T){
      cat('\n')
      print('Done.')
      print("-------------------------------------------------------------------------")
      cat('\n\n')
      print("=========================================================================")
      print("Performance summary:")
      print("=========================================================================")
      print(paste('Parameter value(s): ', paste(round(theta, 3), collapse = ', '), sep = ''))
      print(paste('Probability of accepting H0: ',
                  paste(round(out.OCandASN$acceptH0.prob, 3), collapse = ', '), sep = ''))
      print(paste('Expected sample size: ',
                  paste(round(out.OCandASN$EN, 2), collapse = ', '), sep = ''))
      print("=========================================================================")
      cat('\n')
    }
    
    return(out.OCandASN)
    
    # end one-sided oneZ
  }else{
    
    #################### one-sample z (both sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch.size = design.MSPRT.object$batch.size
      N.max = design.MSPRT.object$N.max
      nAnalyses = design.MSPRT.object$nAnalyses
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      sigma = design.MSPRT.object$sigma
      termination.threshold = design.MSPRT.object$termination.threshold
      nReplicate = design.MSPRT.object$nReplicate
      theta.UMPBT = design.MSPRT.object$theta.UMPBT
      
      # msg
      if(verbose){
        
        if(any(batch.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the group sequential MSPRT for a one-sample z test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the sequential MSPRT for a one-sample z test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Known standard deviation: ", sigma, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(theta.UMPBT$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta.UMPBT$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch.size = c(0, cumsum(batch.size))
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
      if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
      
      ## ignoring N1.max & N2.max
      if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
      if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
      
      ## batch sizes and N.max
      if(missing(batch.size)){
        
        if(missing(N.max)){
          
          return("Either 'batch.size' or 'N.max' needs to be specified")
          
        }else{batch.size = rep(1, N.max)}
        
      }else{
        
        if(missing(N.max)){
          
          N.max = sum(batch.size)
          
        }else{
          
          if(sum(batch.size)!=N.max) return("Sum of batch sizes should add up to N.max")
        }
      }
      
      nAnalyses = length(batch.size)
      
      ######################## UMPBT alternative ########################
      theta.UMPBT = list('right' = UMPBT.alt(test.type = 'oneZ', side = 'right', 
                                             theta0 = theta0, N = N.max, 
                                             Type1 = Type1.target/2, sigma = sigma),
                         'left' = UMPBT.alt(test.type = 'oneZ', side = 'left', 
                                            theta0 = theta0, N = N.max,
                                            Type1 = Type1.target/2, sigma = sigma))
      
      # msg
      if(verbose){
        
        if(any(batch.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the group sequential MSPRT for a one-sample z test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the sequential MSPRT for a one-sample z test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Known standard deviation: ", sigma, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(theta.UMPBT$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta.UMPBT$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch.size = c(0, cumsum(batch.size))
    }
    
    
    registerDoParallel(cores = nCore)
    out.OCandASN = foreach(theta1 = theta, .combine = 'rbind') %dopar% {
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # required storages
      cumsum1_n = LR1_n.r = LR1_n.l = numeric(nReplicate)
      PowerH1.AR = rep(F, nReplicate)
      N1.AR = N1.AR.r = N1.AR.l = rep(N.max, nReplicate)
      decision.underH1.AR.r = decision.underH1.AR.l = rep(NA, nReplicate)
      not.reached.decisionH1.AR = not.reached.decisionH1.AR.r = not.reached.decisionH1.AR.l =
        1:nReplicate
      
      set.seed(seed)
      for(n in 1:nAnalyses){
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          # sum of observations at step n
          sum1_n = rnorm(length(not.reached.decisionH1.AR),
                         (batch.size[n+1]-batch.size[n])*theta1,
                         sqrt(batch.size[n+1]-batch.size[n])*sigma)
          
          # sum of observations until step n
          cumsum1_n[not.reached.decisionH1.AR] =
            cumsum1_n[not.reached.decisionH1.AR] + sum1_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1_n.r[not.reached.decisionH1.AR.r] =
            exp((cumsum1_n[not.reached.decisionH1.AR.r]*(theta.UMPBT$right - theta0) -
                   ((batch.size[n+1]*((theta.UMPBT$right^2) - (theta0^2)))/2))/(sigma^2))
          
          # for left sided check
          LR1_n.l[not.reached.decisionH1.AR.l] =
            exp((cumsum1_n[not.reached.decisionH1.AR.l]*(theta.UMPBT$left - theta0) -
                   ((batch.size[n+1]*((theta.UMPBT$left^2) - (theta0^2)))/2))/(sigma^2))
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1_n.AR.r = LR1_n.r[not.reached.decisionH1.AR.r]<=Accept.threshold
          RejectedH0.underH1_n.AR.r = LR1_n.r[not.reached.decisionH1.AR.r]>=Reject.threshold
          reached.decisionH1_n.AR.r = AcceptedH0.underH1_n.AR.r|RejectedH0.underH1_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1_n.AR.r)){
            
            decision.underH1.AR.r[not.reached.decisionH1.AR.r[AcceptedH0.underH1_n.AR.r]] = 'A'
            decision.underH1.AR.r[not.reached.decisionH1.AR.r[RejectedH0.underH1_n.AR.r]] = 'R'
            N1.AR.r[not.reached.decisionH1.AR.r[reached.decisionH1_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH1.AR.r = not.reached.decisionH1.AR.r[!reached.decisionH1_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1_n.AR.l = LR1_n.l[not.reached.decisionH1.AR.l]<=Accept.threshold
          RejectedH0.underH1_n.AR.l = LR1_n.l[not.reached.decisionH1.AR.l]>=Reject.threshold
          reached.decisionH1_n.AR.l = AcceptedH0.underH1_n.AR.l|RejectedH0.underH1_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1_n.AR.l)){
            
            decision.underH1.AR.l[not.reached.decisionH1.AR.l[AcceptedH0.underH1_n.AR.l]] = 'A'
            decision.underH1.AR.l[not.reached.decisionH1.AR.l[RejectedH0.underH1_n.AR.l]] = 'R'
            N1.AR.l[not.reached.decisionH1.AR.l[reached.decisionH1_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH1.AR.l = not.reached.decisionH1.AR.l[!reached.decisionH1_n.AR.l]
          }
          
          not.reached.decisionH1.AR = union(not.reached.decisionH1.AR.r,
                                            not.reached.decisionH1.AR.l)
        }
      }
      
      
      ### both-sided checking
      ## under H1
      # accepted or rejected ones
      accepted.by.both1 = intersect(which(decision.underH1.AR.r=='A'),
                                    which(decision.underH1.AR.l=='A'))
      onlyrejected.by.right1 = intersect(which(decision.underH1.AR.r=='R'),
                                         which(decision.underH1.AR.l!='R'))
      onlyrejected.by.left1 = intersect(which(decision.underH1.AR.r!='R'),
                                        which(decision.underH1.AR.l=='R'))
      rejected.by.both1 = intersect(which(decision.underH1.AR.r=='R'),
                                    which(decision.underH1.AR.l=='R'))
      
      # sample sizes required
      N1.AR[accepted.by.both1] = pmax(N1.AR.r[accepted.by.both1],
                                      N1.AR.l[accepted.by.both1])
      N1.AR[onlyrejected.by.right1] = N1.AR.r[onlyrejected.by.right1]
      N1.AR[onlyrejected.by.left1] = N1.AR.l[onlyrejected.by.left1]
      N1.AR[rejected.by.both1] = pmin(N1.AR.r[rejected.by.both1],
                                      N1.AR.l[rejected.by.both1])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1 = intersect(which(decision.underH1.AR.r=='A'),
                                         which(is.na(decision.underH1.AR.l)))
      onlyaccepted.by.left1 = intersect(which(is.na(decision.underH1.AR.r)),
                                        which(decision.underH1.AR.l=='A'))
      both.inconclusive1 = intersect(which(is.na(decision.underH1.AR.r)),
                                     which(is.na(decision.underH1.AR.l)))
      all.inconclusive1 = c(onlyaccepted.by.right1, onlyaccepted.by.left1,
                            both.inconclusive1)
      nNot.reached.decisionH1.AR = length(all.inconclusive1)
      
      # power
      PowerH1.AR[c(onlyrejected.by.right1, onlyrejected.by.left1,
                   rejected.by.both1)] = T
      
      ## attained Type II error probability
      # right-sided H1
      actual.PowerH1.AR = mean(PowerH1.AR) +
        sum(c(LR1_n.r[onlyaccepted.by.left1],
              LR1_n.l[onlyaccepted.by.right1],
              pmax(LR1_n.r[both.inconclusive1], LR1_n.l[both.inconclusive1]))>=
              termination.threshold)/nReplicate
      actual.type2.errorH1.AR = 1 - actual.PowerH1.AR
      
      ## Expected sample sizes
      EN1 = mean(N1.AR)
      
      c(theta1, actual.type2.errorH1.AR, EN1)
    }
    
    if(length(theta)==1) out.OCandASN = matrix(data = out.OCandASN, nrow = 1,
                                               ncol = 3, byrow = T)
    
    out.OCandASN = as.data.frame(out.OCandASN)
    colnames(out.OCandASN) = c('theta', 'acceptH0.prob', 'EN')
    
    # msg
    if(verbose==T){
      cat('\n')
      print('Done.')
      print("-------------------------------------------------------------------------")
      cat('\n\n')
      print("=========================================================================")
      print("Performance summary:")
      print("=========================================================================")
      print(paste('Parameter value(s): ', paste(round(theta, 3), collapse = ', '), sep = ''))
      print(paste('Probability of accepting H0: ',
                  paste(round(out.OCandASN$acceptH0.prob, 3), collapse = ', '), sep = ''))
      print(paste('Expected sample size: ',
                  paste(round(out.OCandASN$EN, 2), collapse = ', '), sep = ''))
      print("=========================================================================")
      cat('\n')
    }
    
    return(out.OCandASN)
    
  } # end both-sided oneZ
}

#### one-sample t test ####
OCandASN.MSPRT_oneT = function(theta, design.MSPRT.object, 
                               termination.threshold,
                               side = 'right', theta0 = 0, 
                               Type1.target =.005, Type2.target = .2,
                               N.max, batch.size,
                               nReplicate = 1e+6, nCore = detectCores() - 1,
                               verbose = T, seed = 1){
  
  # side
  if(!missing(design.MSPRT.object)) side = design.MSPRT.object$side
  
  if(side!='both'){
    
    #################### one-sample t (right/left sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch.size = design.MSPRT.object$batch.size
      N.max = design.MSPRT.object$N.max
      nAnalyses = design.MSPRT.object$nAnalyses
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      termination.threshold = design.MSPRT.object$termination.threshold
      nReplicate = design.MSPRT.object$nReplicate
      
      # msg
      if(verbose){
        
        if((batch.size[1]>2)||any(batch.size[-1]>1)){
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the group sequential MSPRT for a one-sample t test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the sequential MSPRT for a one-sample t test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch.size = c(0, cumsum(batch.size))
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
      if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
      
      ## ignoring N1.max & N2.max
      if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
      if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
      
      ## batch sizes and N.max
      if(missing(batch.size)){
        
        if(missing(N.max)){
          
          return("Either 'batch.size' or 'N.max' needs to be specified")
          
        }else{batch.size = c(2, rep(1, N.max-2))}
        
      }else{
        
        if(batch.size[1]<2){
          
          return("First batch size should be at least 2")
          
        }else{
          
          if(missing(N.max)){
            
            N.max = sum(batch.size)
            
          }else{
            
            if(sum(batch.size)!=N.max) return("Sum of batch.size should add up to N.max")
          }
        }
      }
      
      nAnalyses = length(batch.size)
      
      # msg
      if(verbose){
        
        if((batch.size[1]>2)||any(batch.size[-1]>1)){
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the group sequential MSPRT for a one-sample t test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the sequential MSPRT for a one-sample t test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch.size = c(0, cumsum(batch.size))
    }
    
    
    registerDoParallel(cores = nCore)
    out.OCandASN = foreach(theta1 = theta, .combine = 'rbind') %dopar% {
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # cut-off (with sign) in fixed design one-sample t test
      signed_t.alpha = (2*(side=='right')-1)*qt(Type1.target, df = N.max -1, lower.tail = F)
      
      # required storages
      cumSS1_n = cumsum1_n = LR1_n = numeric(nReplicate)
      type2.error.AR = rep(F, nReplicate)
      N1.AR = rep(N.max, nReplicate)
      not.reached.decisionH1.AR = 1:nReplicate
      
      set.seed(seed)
      for(n in 1:nAnalyses){
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          # observations at step n
          if(length(not.reached.decisionH1.AR)>1){
            
            obs1_n = mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                            FUN = function(X){
                              
                              rnorm(length(not.reached.decisionH1.AR), theta1, 1)
                            })
            
          }else{
            
            obs1_n = matrix(mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                                   FUN = function(X){
                                     
                                     rnorm(length(not.reached.decisionH1.AR), theta1, 1)
                                     
                                   }), nrow = 1, ncol = batch.size[n+1]-batch.size[n],
                            byrow = T)
          }
          
          # sum of observations until step n
          cumsum1_n[not.reached.decisionH1.AR] =
            cumsum1_n[not.reached.decisionH1.AR] + rowSums(obs1_n)
          
          # sum of squares of observations until step n
          cumSS1_n[not.reached.decisionH1.AR] = 
            cumSS1_n[not.reached.decisionH1.AR] + rowSums(obs1_n^2)
          
          # xbar and (n-1)*(s^2) until step n
          xbar1_n = cumsum1_n[not.reached.decisionH1.AR]/batch.size[n+1]
          divisor.s1_n.sq = 
            cumSS1_n[not.reached.decisionH1.AR] - ((cumsum1_n[not.reached.decisionH1.AR])^2)/batch.size[n+1]
          
          # likelihood ratio of observations until step n
          LR1_n[not.reached.decisionH1.AR] = 
            ((1 + (batch.size[n+1]*((xbar1_n - theta0)^2))/divisor.s1_n.sq)/
               (1 + (batch.size[n+1]*((xbar1_n - (theta0 + signed_t.alpha*
                                                    sqrt(divisor.s1_n.sq/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s1_n.sq))^(batch.size[n+1]/2)
          
          # comparing with the thresholds
          AcceptedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]<=Accept.threshold)
          RejectedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]>=Reject.threshold)
          reached.decisionH1_n.AR = union(AcceptedH0.underH1_n.AR, RejectedH0.underH1_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH1_n.AR)>0){
            
            N1.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch.size[n+1]
            type2.error.AR[not.reached.decisionH1.AR[AcceptedH0.underH1_n.AR]] = T
            not.reached.decisionH1.AR = not.reached.decisionH1.AR[-reached.decisionH1_n.AR]
          }
        }
      }
      
      # attained Type II error probability
      actual.type2.error.AR = mean(type2.error.AR) +
        sum(LR1_n[not.reached.decisionH1.AR]<termination.threshold)/nReplicate
      
      # Expected sample sizes
      EN1 = mean(N1.AR)
      
      c(theta1, actual.type2.error.AR, EN1)
    }
    
    if(length(theta)==1) out.OCandASN = matrix(data = out.OCandASN, nrow = 1,
                                               ncol = 3, byrow = T)
    
    out.OCandASN = as.data.frame(out.OCandASN)
    colnames(out.OCandASN) = c('theta', 'acceptH0.prob', 'EN')
    
    # msg
    if(verbose==T){
      cat('\n')
      print('Done.')
      print("-------------------------------------------------------------------------")
      cat('\n\n')
      print("=========================================================================")
      print("Performance summary:")
      print("=========================================================================")
      print(paste('Parameter value(s): ', paste(round(theta, 3), collapse = ', '), sep = ''))
      print(paste('Probability of accepting H0: ',
                  paste(round(out.OCandASN$acceptH0.prob, 3), collapse = ', '), sep = ''))
      print(paste('Expected sample size: ',
                  paste(round(out.OCandASN$EN, 2), collapse = ', '), sep = ''))
      print("=========================================================================")
      cat('\n')
    }
    
    return(out.OCandASN)
    
    # end one-sided oneT
  }else{
    
    #################### one-sample t (both sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch.size = design.MSPRT.object$batch.size
      N.max = design.MSPRT.object$N.max
      nAnalyses = design.MSPRT.object$nAnalyses
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      termination.threshold = design.MSPRT.object$termination.threshold
      nReplicate = design.MSPRT.object$nReplicate
      
      # msg
      if(verbose){
        
        if((batch.size[1]>2)||any(batch.size[-1]>1)){
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the group sequential MSPRT for a one-sample t test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the sequential MSPRT for a one-sample t test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch.size = c(0, cumsum(batch.size))
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
      if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
      
      ## ignoring N1.max & N2.max
      if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
      if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
      
      ## batch sizes and N.max
      if(missing(batch.size)){
        
        if(missing(N.max)){
          
          return("Either 'batch.size' or 'N.max' needs to be specified")
          
        }else{batch.size = c(2, rep(1, N.max-2))}
        
      }else{
        
        if(batch.size[1]<2){
          
          return("First batch size should be at least 2")
          
        }else{
          
          if(missing(N.max)){
            
            N.max = sum(batch.size)
            
          }else{
            
            if(sum(batch.size)!=N.max) return("Sum of batch.size should add up to N.max")
          }
        }
      }
      
      nAnalyses = length(batch.size)
      
      # msg
      if(verbose){
        
        if((batch.size[1]>2)||any(batch.size[-1]>1)){
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the group sequential MSPRT for a one-sample t test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("OC and ASN of the sequential MSPRT for a one-sample t test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch.size = c(0, cumsum(batch.size))
    }
    
    
    registerDoParallel(cores = nCore)
    out.OCandASN = foreach(theta1 = theta, .combine = 'rbind') %dopar% {
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # cut-off (with sign) in fixed design one-sample t test
      t.alpha = qt(Type1.target/2, df = N.max -1, lower.tail = F)
      
      # required storages
      cumSS1_n = cumsum1_n = LR1_n.r = LR1_n.l = numeric(nReplicate)
      PowerH1.AR = rep(F, nReplicate)
      N1.AR = N1.AR.r = N1.AR.l = rep(N.max, nReplicate)
      decision.underH1.AR.r = decision.underH1.AR.l = rep(NA, nReplicate)
      not.reached.decisionH1.AR = not.reached.decisionH1.AR.r = not.reached.decisionH1.AR.l =
        1:nReplicate
      
      set.seed(seed)
      for(n in 1:nAnalyses){
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          # observations at step n
          if(length(not.reached.decisionH1.AR)>1){
            
            obs1_n = mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                            FUN = function(X){
                              
                              rnorm(length(not.reached.decisionH1.AR), theta1, 1)
                            })
            
          }else{
            
            obs1_n = matrix(mapply(X = 1:(batch.size[n+1]-batch.size[n]), 
                                   FUN = function(X){
                                     
                                     rnorm(length(not.reached.decisionH1.AR), theta1, 1)
                                     
                                   }), nrow = 1, ncol = batch.size[n+1]-batch.size[n],
                            byrow = T)
          }
          
          # sum of observations until step n
          cumsum1_n[not.reached.decisionH1.AR] = 
            cumsum1_n[not.reached.decisionH1.AR] + rowSums(obs1_n)
          
          # sum of squares of observations until step n
          cumSS1_n[not.reached.decisionH1.AR] = 
            cumSS1_n[not.reached.decisionH1.AR] + rowSums(obs1_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          # for right sided check
          xbar1_n.r = cumsum1_n[not.reached.decisionH1.AR.r]/batch.size[n+1]
          divisor.s1_n.sq.r = 
            cumSS1_n[not.reached.decisionH1.AR.r] - 
            ((cumsum1_n[not.reached.decisionH1.AR.r])^2)/batch.size[n+1]
          
          # for left sided check
          xbar1_n.l = cumsum1_n[not.reached.decisionH1.AR.l]/batch.size[n+1]
          divisor.s1_n.sq.l = 
            cumSS1_n[not.reached.decisionH1.AR.l] - 
            ((cumsum1_n[not.reached.decisionH1.AR.l])^2)/batch.size[n+1]
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1_n.r[not.reached.decisionH1.AR.r] = 
            ((1 + (batch.size[n+1]*((xbar1_n.r - theta0)^2))/divisor.s1_n.sq.r)/
               (1 + (batch.size[n+1]*((xbar1_n.r - 
                                         (theta0 + t.alpha*
                                            sqrt(divisor.s1_n.sq.r/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s1_n.sq.r))^(batch.size[n+1]/2)
          
          # for left sided check
          LR1_n.l[not.reached.decisionH1.AR.l] = 
            ((1 + (batch.size[n+1]*((xbar1_n.l - theta0)^2))/divisor.s1_n.sq.l)/
               (1 + (batch.size[n+1]*((xbar1_n.l - 
                                         (theta0 - t.alpha*
                                            sqrt(divisor.s1_n.sq.l/(N.max*(batch.size[n+1]-1)))))^2))/
                  divisor.s1_n.sq.l))^(batch.size[n+1]/2)
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1_n.AR.r = LR1_n.r[not.reached.decisionH1.AR.r]<=Accept.threshold
          RejectedH0.underH1_n.AR.r = LR1_n.r[not.reached.decisionH1.AR.r]>=Reject.threshold
          reached.decisionH1_n.AR.r = AcceptedH0.underH1_n.AR.r|RejectedH0.underH1_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1_n.AR.r)){
            
            decision.underH1.AR.r[not.reached.decisionH1.AR.r[AcceptedH0.underH1_n.AR.r]] = 'A'
            decision.underH1.AR.r[not.reached.decisionH1.AR.r[RejectedH0.underH1_n.AR.r]] = 'R'
            N1.AR.r[not.reached.decisionH1.AR.r[reached.decisionH1_n.AR.r]] = batch.size[n+1]
            not.reached.decisionH1.AR.r = not.reached.decisionH1.AR.r[!reached.decisionH1_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1_n.AR.l = LR1_n.l[not.reached.decisionH1.AR.l]<=Accept.threshold
          RejectedH0.underH1_n.AR.l = LR1_n.l[not.reached.decisionH1.AR.l]>=Reject.threshold
          reached.decisionH1_n.AR.l = AcceptedH0.underH1_n.AR.l|RejectedH0.underH1_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1_n.AR.l)){
            
            decision.underH1.AR.l[not.reached.decisionH1.AR.l[AcceptedH0.underH1_n.AR.l]] = 'A'
            decision.underH1.AR.l[not.reached.decisionH1.AR.l[RejectedH0.underH1_n.AR.l]] = 'R'
            N1.AR.l[not.reached.decisionH1.AR.l[reached.decisionH1_n.AR.l]] = batch.size[n+1]
            not.reached.decisionH1.AR.l = not.reached.decisionH1.AR.l[!reached.decisionH1_n.AR.l]
          }
          
          not.reached.decisionH1.AR = union(not.reached.decisionH1.AR.r,
                                            not.reached.decisionH1.AR.l)
        }
      }
      
      
      ### both-sided checking
      ## under H1
      # accepted or rejected ones
      accepted.by.both1 = intersect(which(decision.underH1.AR.r=='A'),
                                    which(decision.underH1.AR.l=='A'))
      onlyrejected.by.right1 = intersect(which(decision.underH1.AR.r=='R'),
                                         which(decision.underH1.AR.l!='R'))
      onlyrejected.by.left1 = intersect(which(decision.underH1.AR.r!='R'),
                                        which(decision.underH1.AR.l=='R'))
      rejected.by.both1 = intersect(which(decision.underH1.AR.r=='R'),
                                    which(decision.underH1.AR.l=='R'))
      
      # sample sizes required
      N1.AR[accepted.by.both1] = pmax(N1.AR.r[accepted.by.both1],
                                      N1.AR.l[accepted.by.both1])
      N1.AR[onlyrejected.by.right1] = N1.AR.r[onlyrejected.by.right1]
      N1.AR[onlyrejected.by.left1] = N1.AR.l[onlyrejected.by.left1]
      N1.AR[rejected.by.both1] = pmin(N1.AR.r[rejected.by.both1],
                                      N1.AR.l[rejected.by.both1])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1 = intersect(which(decision.underH1.AR.r=='A'),
                                         which(is.na(decision.underH1.AR.l)))
      onlyaccepted.by.left1 = intersect(which(is.na(decision.underH1.AR.r)),
                                        which(decision.underH1.AR.l=='A'))
      both.inconclusive1 = intersect(which(is.na(decision.underH1.AR.r)),
                                     which(is.na(decision.underH1.AR.l)))
      all.inconclusive1 = c(onlyaccepted.by.right1, onlyaccepted.by.left1,
                            both.inconclusive1)
      nNot.reached.decisionH1.AR = length(all.inconclusive1)
      
      # Type I error probability
      PowerH1.AR[c(onlyrejected.by.right1, onlyrejected.by.left1,
                   rejected.by.both1)] = T
      
      ## attained Type II error probability
      actual.PowerH1.AR.r = mean(PowerH1.AR) +
        sum(c(LR1_n.r[onlyaccepted.by.left1],
              LR1_n.l[onlyaccepted.by.right1],
              pmax(LR1_n.r[both.inconclusive1], LR1_n.l[both.inconclusive1]))>=
              termination.threshold)/nReplicate
      actual.type2.errorH1.AR = 1 - actual.PowerH1.AR.r
      
      ## Expected sample sizes
      EN1 = mean(N1.AR)
      
      c(theta1, actual.type2.errorH1.AR, EN1)
    }
    
    if(length(theta)==1) out.OCandASN = matrix(data = out.OCandASN, nrow = 1,
                                               ncol = 3, byrow = T)
    
    out.OCandASN = as.data.frame(out.OCandASN)
    colnames(out.OCandASN) = c('theta', 'acceptH0.prob', 'EN')
    
    # msg
    if(verbose==T){
      cat('\n')
      print('Done.')
      print("-------------------------------------------------------------------------")
      cat('\n\n')
      print("=========================================================================")
      print("Performance summary:")
      print("=========================================================================")
      print(paste('Parameter value(s): ', paste(round(theta, 3), collapse = ', '), sep = ''))
      print(paste('Probability of accepting H0: ',
                  paste(round(out.OCandASN$acceptH0.prob, 3), collapse = ', '), sep = ''))
      print(paste('Expected sample size: ',
                  paste(round(out.OCandASN$EN, 2), collapse = ', '), sep = ''))
      print("=========================================================================")
      cat('\n')
    }
    
    return(out.OCandASN)
    
  } # end both-sided oneT
  
  # end oneT
}

#### two-sample z test ####
OCandASN.MSPRT_twoZ = function(theta, design.MSPRT.object, 
                               termination.threshold,
                               side = 'right', theta0 = 0, 
                               Type1.target =.005, Type2.target = .2,
                               N1.max, N2.max, sigma1 = 1, sigma2 = 1,
                               batch1.size, batch2.size,
                               nReplicate = 1e+6, nCore = detectCores() - 1,
                               verbose = T, seed = 1){
  
  # side
  if(!missing(design.MSPRT.object)) side = design.MSPRT.object$side
  
  if(side!='both'){
    
    #################### two-sample z (right/left sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch1.size = design.MSPRT.object$batch1.size
      batch2.size = design.MSPRT.object$batch2.size
      N1.max = design.MSPRT.object$N1.max
      N2.max = design.MSPRT.object$N2.max
      nAnalyses = design.MSPRT.object$nAnalyses
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      sigma1 = design.MSPRT.object$sigma1
      sigma2 = design.MSPRT.object$sigma2
      termination.threshold = design.MSPRT.object$termination.threshold
      theta.UMPBT = design.MSPRT.object$theta.UMPBT
      nReplicate = design.MSPRT.object$nReplicate
      
      # msg
      if(verbose){
        
        if(any(batch1.size>1)||any(batch2.size>1)){
          
          cat('\n')
          print("=========================================================================")
          print("OC and ASN of the group sequential MSPRT for a two-sample z test:")
          print("=========================================================================")
          
        }else{
          
          cat('\n')
          print("=========================================================================")
          print("OC and ASN of the sequential MSPRT for a two-sample z test:")
          print("=========================================================================")
        }
        
        print("Group 1:")
        print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
        print(paste(" Known standard deviation: ", sigma1, sep = ""))
        print("Group 2:")
        print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
        print(paste(" Known standard deviation: ", sigma2, sep = ""))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(theta.UMPBT, 3)))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch1.size = c(0, cumsum(batch1.size))
      batch2.size = c(0, cumsum(batch2.size))
      
    }else{
      
      ## ignoring batch.seq
      if(!missing(batch.size)) print("'batch.size' is ignored. Not required in two-sample tests.")
      
      ## ignoring N.max
      if(!missing(N.max)) print("'N.max' is ignored. Not required in two-sample tests.")
      
      ## checking if length(batch1.size) and length(batch2.size) are equal
      if((!missing(batch1.size)) && (!missing(batch2.size)) &&
         (length(batch1.size)!=length(batch2.size))) return("Lenghts of batch1.size and batch2.size should be same")
      
      ## batch sizes and N for group 1
      if(missing(batch1.size)){
        
        if(missing(N1.max)){
          
          return(print("Either 'batch1.size' or 'N1.max' needs to be specified"))
          
        }else{batch1.size = rep(1, N1.max)}
        
      }else{
        
        if(missing(N1.max)){
          
          N1.max = sum(batch1.size)
          
        }else{
          
          if(sum(batch1.size)!=N1.max) return(print("Sum of batch1.size should add up to N1.max"))
        }
      }
      
      ## batch sizes and N for group 2
      if(missing(batch2.size)){
        
        if(missing(N2.max)){
          
          return(print("Either 'batch2.size' or 'N2.max' needs to be specified"))
          
        }else{batch2.size = rep(1, N2.max)}
        
      }else{
        
        if(missing(N2.max)){
          
          N2.max = sum(batch2.size)
          
        }else{
          
          if(sum(batch2.size)!=N1.max) return(print("Sum of batch2.size should add up to N2.max"))
        }
      }
      
      nAnalyses = length(batch1.size)
      
      ################ UMPBT alternative ################
      theta.UMPBT = UMPBT.alt(test.type = 'twoZ', side = side, theta0 = theta0,
                              N1 = N1.max, N2 = N2.max, Type1 = Type1.target,
                              sigma1 = sigma1, sigma2 = sigma2)
      
      # msg
      if(verbose){
        
        if(any(batch1.size>1)||any(batch2.size>1)){
          
          cat('\n')
          print("=========================================================================")
          print("Designing the group sequential MSPRT for a two-sample z test:")
          print("=========================================================================")
          
        }else{
          
          cat('\n')
          print("=========================================================================")
          print("Designing the sequential MSPRT for a two-sample z test:")
          print("=========================================================================")
        }
        
        print("Group 1:")
        print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
        print(paste(" Known standard deviation: ", sigma1, sep = ""))
        print("Group 2:")
        print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
        print(paste(" Known standard deviation: ", sigma2, sep = ""))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(theta.UMPBT, 3)))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch1.size = c(0, cumsum(batch1.size))
      batch2.size = c(0, cumsum(batch2.size))
    }
    
    
    registerDoParallel(cores = nCore)
    out.OCandASN = foreach(theta1 = theta, .combine = 'rbind') %dopar% {
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # required storages
      cumsum11_n = cumsum21_n = LR1_n = numeric(nReplicate)
      type2.error.AR = rep(F, nReplicate)
      N11.AR = rep(N1.max, nReplicate)
      N21.AR = rep(N2.max, nReplicate)
      not.reached.decisionH1.AR = 1:nReplicate
      
      set.seed(seed)
      for(n in 1:nAnalyses){
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          ## sum of observations at step n
          # Group 1
          sum11_n = rnorm(length(not.reached.decisionH1.AR),
                          (batch1.size[n+1]-batch1.size[n])*(theta1/2),
                          sqrt(batch1.size[n+1]-batch1.size[n])*sigma1)
          
          # Group 2
          sum21_n = rnorm(length(not.reached.decisionH1.AR),
                          -(batch2.size[n+1]-batch2.size[n])*(theta1/2),
                          sqrt(batch2.size[n+1]-batch2.size[n])*sigma2)
          
          ## sum of observations until step n
          # Group 1
          cumsum11_n[not.reached.decisionH1.AR] = 
            cumsum11_n[not.reached.decisionH1.AR] + sum11_n
          
          # Group 2
          cumsum21_n[not.reached.decisionH1.AR] = 
            cumsum21_n[not.reached.decisionH1.AR] + sum21_n
          
          # likelihood ratio of observations until step n
          LR1_n[not.reached.decisionH1.AR] = 
            exp(-(((theta.UMPBT^2) - (theta0^2)) - 
                    2*(theta.UMPBT - theta0)*
                    (cumsum11_n[not.reached.decisionH1.AR]/batch1.size[n+1] - 
                       cumsum21_n[not.reached.decisionH1.AR]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          # comparing with the thresholds
          AcceptedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]<=Accept.threshold)
          RejectedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]>=Reject.threshold)
          reached.decisionH1_n.AR = union(AcceptedH0.underH1_n.AR, RejectedH0.underH1_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH1_n.AR)>0){
            
            N11.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch1.size[n+1]
            N21.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch2.size[n+1]
            type2.error.AR[not.reached.decisionH1.AR[AcceptedH0.underH1_n.AR]] = T
            not.reached.decisionH1.AR = not.reached.decisionH1.AR[-reached.decisionH1_n.AR]
          }
        }
      }
      
      # attained Type II error probability
      actual.type2.error.AR = mean(type2.error.AR) +
        sum(LR1_n[not.reached.decisionH1.AR]<termination.threshold)/nReplicate
      
      # Expected sample sizes
      EN11 = mean(N11.AR)
      EN21 = mean(N21.AR)
      
      c(theta1, actual.type2.error.AR, EN11, EN21)
    }
    
    if(length(theta)==1) out.OCandASN = matrix(data = out.OCandASN, nrow = 1,
                                               ncol = 4, byrow = T)
    
    out.OCandASN = as.data.frame(out.OCandASN)
    colnames(out.OCandASN) = c('theta', 'acceptH0.prob', 'EN1', 'EN2')
    
    # msg
    if(verbose==T){
      cat('\n')
      print('Done.')
      print("-------------------------------------------------------------------------")
      cat('\n\n')
      print("=========================================================================")
      print("Performance summary:")
      print("=========================================================================")
      print(paste('Parameter value(s): ', paste(round(theta, 3), collapse = ', '), sep = ''))
      print(paste('Probability of accepting H0: ',
                  paste(round(out.OCandASN$acceptH0.prob, 3), collapse = ', '), sep = ''))
      print("Expected sample size:")
      print(paste(' Group 1 - ', paste(round(out.OCandASN$EN1, 2), collapse = ', '), sep = ''))
      print(paste(' Group 2 - ', paste(round(out.OCandASN$EN2, 2), collapse = ', '), sep = ''))
      print("=========================================================================")
      cat('\n')
    }
    
    return(out.OCandASN)
    
    # end one-sided twoZ
  }else{
    
    #################### two-sample z (both sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch1.size = design.MSPRT.object$batch1.size
      batch2.size = design.MSPRT.object$batch2.size
      N1.max = design.MSPRT.object$N1.max
      N2.max = design.MSPRT.object$N2.max
      nAnalyses = design.MSPRT.object$nAnalyses
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      sigma1 = design.MSPRT.object$sigma1
      sigma2 = design.MSPRT.object$sigma2
      termination.threshold = design.MSPRT.object$termination.threshold
      theta.UMPBT = design.MSPRT.object$theta.UMPBT
      nReplicate = design.MSPRT.object$nReplicate
      
      # msg
      if(verbose){
        
        if(any(batch1.size>1)||any(batch2.size>1)){
          
          cat('\n')
          print("=========================================================================")
          print("OC and ASN of the group sequential MSPRT for a two-sample z test:")
          print("=========================================================================")
          
        }else{
          
          cat('\n')
          print("=========================================================================")
          print("OC and ASN of the sequential MSPRT for a two-sample z test:")
          print("=========================================================================")
        }
        
        print("Group 1:")
        print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
        print(paste(" Known standard deviation: ", sigma1, sep = ""))
        print("Group 2:")
        print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
        print(paste(" Known standard deviation: ", sigma2, sep = ""))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(theta.UMPBT$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta.UMPBT$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch1.size = c(0, cumsum(batch1.size))
      batch2.size = c(0, cumsum(batch2.size))
      
    }else{
      
      ## ignoring batch.size
      if(!missing(batch.size)) print("'batch.size' is ignored. Not required in two-sample tests.")
      
      ## ignoring N.max
      if(!missing(N.max)) print("'N.max' is ignored. Not required in two-sample tests.")
      
      ## checking if length(batch1.size) and length(batch2.size) are equal
      if((!missing(batch1.size)) && (!missing(batch2.size)) &&
         (length(batch1.size)!=length(batch2.size))) return("Lenghts of batch1.size and batch2.size should be same")
      
      ## batch sizes and N for group 1
      if(missing(batch1.size)){
        
        if(missing(N1.max)){
          
          return(print("Either 'batch1.size' or 'N1.max' needs to be specified"))
          
        }else{batch1.size = rep(1, N1.max)}
        
      }else{
        
        if(missing(N1.max)){
          
          N1.max = sum(batch1.size)
          
        }else{
          
          if(sum(batch1.size)!=N1.max) return(print("Sum of batch1.size should add up to N1.max"))
        }
      }
      
      ## batch sizes and N for group 2
      if(missing(batch2.size)){
        
        if(missing(N2.max)){
          
          return(print("Either 'batch2.size' or 'N2.max' needs to be specified"))
          
        }else{batch2.size = rep(1, N2.max)}
        
      }else{
        
        if(missing(N2.max)){
          
          N2.max = sum(batch2.size)
          
        }else{
          
          if(sum(batch2.size)!=N1.max) return(print("Sum of batch2.size should add up to N2.max"))
        }
      }
      
      nAnalyses = length(batch1.size)
      
      ################ UMPBT alternative ################
      theta.UMPBT = list('right' = UMPBT.alt(test.type = 'twoZ', side = 'right',
                                             theta0 = theta0, N1 = N1.max, N2 = N2.max, 
                                             Type1 = Type1.target/2,
                                             sigma1 = sigma1, sigma2 = sigma2),
                         'left' = UMPBT.alt(test.type = 'twoZ', side = 'left',
                                            theta0 = theta0, N1 = N1.max, N2 = N2.max, 
                                            Type1 = Type1.target/2,
                                            sigma1 = sigma1, sigma2 = sigma2))
      
      # msg
      if(verbose){
        
        if(any(batch1.size>1)||any(batch2.size>1)){
          
          cat('\n')
          print("=========================================================================")
          print("OC and ASN of the group sequential MSPRT for a two-sample z test:")
          print("=========================================================================")
          
        }else{
          
          cat('\n')
          print("=========================================================================")
          print("OC and ASN of the sequential MSPRT for a two-sample z test:")
          print("=========================================================================")
        }
        
        print("Group 1:")
        print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
        print(paste(" Known standard deviation: ", sigma1, sep = ""))
        print("Group 2:")
        print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
        print(paste(" Known standard deviation: ", sigma2, sep = ""))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(theta.UMPBT$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta.UMPBT$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch1.size = c(0, cumsum(batch1.size))
      batch2.size = c(0, cumsum(batch2.size))
    }
    
    
    registerDoParallel(cores = nCore)
    out.OCandASN = foreach(theta1 = theta, .combine = 'rbind') %dopar% {
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # required storages
      cumsum11_n = cumsum21_n = LR1_n.r = LR1_n.l = numeric(nReplicate)
      PowerH1.AR = rep(F, nReplicate)
      N11.AR = N11.AR.r = N11.AR.l = rep(N1.max, nReplicate)
      N21.AR = N21.AR.r = N21.AR.l = rep(N2.max, nReplicate)
      decision.underH1.AR.r = decision.underH1.AR.l = rep(NA, nReplicate)
      not.reached.decisionH1.AR = not.reached.decisionH1.AR.r = not.reached.decisionH1.AR.l =
        1:nReplicate
      
      set.seed(seed)
      for(n in 1:nAnalyses){
        
        ## under right-sided H1
        if(length(not.reached.decisionH1.AR)>0){
          
          ## sum of observations at step n
          # Group 1
          sum11_n = rnorm(length(not.reached.decisionH1.AR),
                          (batch1.size[n+1]-batch1.size[n])*(theta1/2),
                          sqrt(batch1.size[n+1]-batch1.size[n])*sigma1)
          
          # Group 2
          sum21_n = rnorm(length(not.reached.decisionH1.AR),
                          -(batch2.size[n+1]-batch2.size[n])*(theta1/2),
                          sqrt(batch2.size[n+1]-batch2.size[n])*sigma2)
          
          ## sum of observations until step n
          # Group 1
          cumsum11_n[not.reached.decisionH1.AR] = 
            cumsum11_n[not.reached.decisionH1.AR] + sum11_n
          
          # Group 2
          cumsum21_n[not.reached.decisionH1.AR] = 
            cumsum21_n[not.reached.decisionH1.AR] + sum21_n
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1_n.r[not.reached.decisionH1.AR.r] = 
            exp(-(((theta.UMPBT$right^2) - (theta0^2)) - 
                    2*(theta.UMPBT$right - theta0)*
                    (cumsum11_n[not.reached.decisionH1.AR.r]/batch1.size[n+1] - 
                       cumsum21_n[not.reached.decisionH1.AR.r]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          # for left sided check
          LR1_n.l[not.reached.decisionH1.AR.l] = 
            exp(-(((theta.UMPBT$left^2) - (theta0^2)) - 
                    2*(theta.UMPBT$left - theta0)*
                    (cumsum11_n[not.reached.decisionH1.AR.l]/batch1.size[n+1] - 
                       cumsum21_n[not.reached.decisionH1.AR.l]/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1_n.AR.r = LR1_n.r[not.reached.decisionH1.AR.r]<=Accept.threshold
          RejectedH0.underH1_n.AR.r = LR1_n.r[not.reached.decisionH1.AR.r]>=Reject.threshold
          reached.decisionH1_n.AR.r = AcceptedH0.underH1_n.AR.r|RejectedH0.underH1_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1_n.AR.r)){
            
            decision.underH1.AR.r[not.reached.decisionH1.AR.r[AcceptedH0.underH1_n.AR.r]] = 'A'
            decision.underH1.AR.r[not.reached.decisionH1.AR.r[RejectedH0.underH1_n.AR.r]] = 'R'
            N11.AR.r[not.reached.decisionH1.AR.r[reached.decisionH1_n.AR.r]] = batch1.size[n+1]
            N21.AR.r[not.reached.decisionH1.AR.r[reached.decisionH1_n.AR.r]] = batch2.size[n+1]
            not.reached.decisionH1.AR.r = not.reached.decisionH1.AR.r[!reached.decisionH1_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1_n.AR.l = LR1_n.l[not.reached.decisionH1.AR.l]<=Accept.threshold
          RejectedH0.underH1_n.AR.l = LR1_n.l[not.reached.decisionH1.AR.l]>=Reject.threshold
          reached.decisionH1_n.AR.l = AcceptedH0.underH1_n.AR.l|RejectedH0.underH1_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1_n.AR.l)){
            
            decision.underH1.AR.l[not.reached.decisionH1.AR.l[AcceptedH0.underH1_n.AR.l]] = 'A'
            decision.underH1.AR.l[not.reached.decisionH1.AR.l[RejectedH0.underH1_n.AR.l]] = 'R'
            N11.AR.l[not.reached.decisionH1.AR.l[reached.decisionH1_n.AR.l]] = batch1.size[n+1]
            N21.AR.l[not.reached.decisionH1.AR.l[reached.decisionH1_n.AR.l]] = batch2.size[n+1]
            not.reached.decisionH1.AR.l = not.reached.decisionH1.AR.l[!reached.decisionH1_n.AR.l]
          }
          
          not.reached.decisionH1.AR = union(not.reached.decisionH1.AR.r,
                                            not.reached.decisionH1.AR.l)
        }
      }
      
      
      ### both-sided checking
      ## under H1
      # accepted or rejected ones
      accepted.by.both1 = intersect(which(decision.underH1.AR.r=='A'),
                                    which(decision.underH1.AR.l=='A'))
      onlyrejected.by.right1 = intersect(which(decision.underH1.AR.r=='R'),
                                         which(decision.underH1.AR.l!='R'))
      onlyrejected.by.left1 = intersect(which(decision.underH1.AR.r!='R'),
                                        which(decision.underH1.AR.l=='R'))
      rejected.by.both1 = intersect(which(decision.underH1.AR.r=='R'),
                                    which(decision.underH1.AR.l=='R'))
      
      ## sample sizes required
      # Group 1
      N11.AR[accepted.by.both1] = pmax(N11.AR.r[accepted.by.both1],
                                       N11.AR.l[accepted.by.both1])
      N11.AR[onlyrejected.by.right1] = N11.AR.r[onlyrejected.by.right1]
      N11.AR[onlyrejected.by.left1] = N11.AR.l[onlyrejected.by.left1]
      N11.AR[rejected.by.both1] = pmin(N11.AR.r[rejected.by.both1],
                                       N11.AR.l[rejected.by.both1])
      
      # Group 2
      N21.AR[accepted.by.both1] = pmax(N21.AR.r[accepted.by.both1],
                                       N21.AR.l[accepted.by.both1])
      N21.AR[onlyrejected.by.right1] = N21.AR.r[onlyrejected.by.right1]
      N21.AR[onlyrejected.by.left1] = N21.AR.l[onlyrejected.by.left1]
      N21.AR[rejected.by.both1] = pmin(N21.AR.r[rejected.by.both1],
                                       N21.AR.l[rejected.by.both1])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1 = intersect(which(decision.underH1.AR.r=='A'),
                                         which(is.na(decision.underH1.AR.l)))
      onlyaccepted.by.left1 = intersect(which(is.na(decision.underH1.AR.r)),
                                        which(decision.underH1.AR.l=='A'))
      both.inconclusive1 = intersect(which(is.na(decision.underH1.AR.r)),
                                     which(is.na(decision.underH1.AR.l)))
      all.inconclusive1 = c(onlyaccepted.by.right1, onlyaccepted.by.left1,
                            both.inconclusive1)
      nNot.reached.decisionH1.AR = length(all.inconclusive1)
      
      # Type I error probability
      PowerH1.AR[c(onlyrejected.by.right1, onlyrejected.by.left1,
                   rejected.by.both1)] = T
      
      ## attained Type II error probability
      actual.PowerH1.AR = mean(PowerH1.AR) +
        sum(c(LR1_n.r[onlyaccepted.by.left1],
              LR1_n.l[onlyaccepted.by.right1],
              pmax(LR1_n.r[both.inconclusive1], LR1_n.l[both.inconclusive1]))>=
              termination.threshold)/nReplicate
      actual.type2.errorH1.AR = 1 - actual.PowerH1.AR
      
      ## Expected sample sizes
      # Group 1
      EN11 = mean(N11.AR)
      
      # Group 2
      EN21 = mean(N21.AR)
      
      c(theta1, actual.type2.errorH1.AR, EN11, EN21)
    }
    
    if(length(theta)==1) out.OCandASN = matrix(data = out.OCandASN, nrow = 1,
                                               ncol = 4, byrow = T)
    
    out.OCandASN = as.data.frame(out.OCandASN)
    colnames(out.OCandASN) = c('theta', 'acceptH0.prob', 'EN1', 'EN2')
    
    # msg
    if(verbose==T){
      cat('\n')
      print('Done.')
      print("-------------------------------------------------------------------------")
      cat('\n\n')
      print("=========================================================================")
      print("Performance summary:")
      print("=========================================================================")
      print(paste('Parameter value(s): ', paste(round(theta, 3), collapse = ', '), sep = ''))
      print(paste('Probability of accepting H0: ',
                  paste(round(out.OCandASN$acceptH0.prob, 3), collapse = ', '), sep = ''))
      print("Expected sample size:")
      print(paste(' Group 1 - ', paste(round(out.OCandASN$EN1, 2), collapse = ', '), sep = ''))
      print(paste(' Group 2 - ', paste(round(out.OCandASN$EN2, 2), collapse = ', '), sep = ''))
      print("=========================================================================")
      cat('\n')
    }
    
    return(out.OCandASN)
    
  } # end both-sided twoZ
  
  # end twoZ
}

#### two-sample t test ####
OCandASN.MSPRT_twoT = function(theta, design.MSPRT.object, 
                               termination.threshold,
                               side = 'right', theta0 = 0, 
                               Type1.target =.005, Type2.target = .2,
                               N1.max, N2.max, batch1.size, batch2.size,
                               nReplicate = 1e+6, nCore = detectCores() - 1,
                               verbose = T, seed = 1){
  
  # side
  if(!missing(design.MSPRT.object)) side = design.MSPRT.object$side
  
  if(side!='both'){
    
    #################### two-sample t (right/left sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch1.size = design.MSPRT.object$batch1.size
      batch2.size = design.MSPRT.object$batch2.size
      N1.max = design.MSPRT.object$N1.max
      N2.max = design.MSPRT.object$N2.max
      nAnalyses = design.MSPRT.object$nAnalyses
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      termination.threshold = design.MSPRT.object$termination.threshold
      nReplicate = design.MSPRT.object$nReplicate
      
      # msg
      if(verbose){
        
        if(any(batch1.size>1)||any(batch2.size>1)){
          
          cat('\n')
          print("=========================================================================")
          print("Designing the group sequential MSPRT for a two-sample t test:")
          print("=========================================================================")
          
        }else{
          
          cat('\n')
          print("=========================================================================")
          print("Designing the sequential MSPRT for a two-sample t test:")
          print("=========================================================================")
        }
        
        print("Group 1:")
        print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
        print("Group 2:")
        print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch1.size = c(0, cumsum(batch1.size))
      batch2.size = c(0, cumsum(batch2.size))
      
    }else{
      
      ## ignoring batch.size
      if(!missing(batch.size)) print("'batch.size' is ignored. Not required in two-sample tests.")
      
      ## ignoring N.max
      if(!missing(N.max)) print("'N.max' is ignored. Not required in two-sample tests.")
      
      ## checking if length(batch1.size) and length(batch2.size) are equal
      if((!missing(batch1.size)) && (!missing(batch2.size)) &&
         (length(batch1.size)!=length(batch2.size))) return("Lenghts of batch1.size and batch2.size should be same")
      
      ## batch sizes and N for group 1
      if(missing(batch1.size)){
        
        if(missing(N1.max)){
          
          return("Either 'batch1.size' or 'N1.max' needs to be specified")
          
        }else{batch1.size = c(2, rep(1, N1.max-2))}
        
      }else{
        
        if(batch1.size[1]<2){
          
          return("First batch size in Group 1 should be at least 2")
          
        }else{
          
          if(missing(N1.max)){
            
            N1.max = sum(batch1.size)
            
          }else{
            
            if(sum(batch1.size)!=N1.max) return("Sum of batch1.size should add up to N1.max")
          }
        }
      }
      
      ## batch sizes and N for group 2
      if(missing(batch2.size)){
        
        if(missing(N2.max)){
          
          return("Either 'batch2.size' or 'N2.max' needs to be specified")
          
        }else{batch2.size = c(2, rep(1, N2.max-2))}
        
      }else{
        
        if(batch2.size[1]<2){
          
          return("First batch size in Group 2 should be at least 2")
          
        }else{
          
          if(missing(N2.max)){
            
            N2.max = sum(batch2.size)
            
          }else{
            
            if(sum(batch2.size)!=N2.max) return("Sum of batch2.size should add up to N2.max")
          }
        }
      }
      
      nAnalyses = length(batch1.size)
      
      # msg
      if(verbose){
        
        if((batch1.size[1]>2)||any(batch1.size[-1]>1)||
           (batch2.size[1]>2)||any(batch2.size[-1]>1)){
          
          cat('\n')
          print("=========================================================================")
          print("OC and ASN of the group sequential MSPRT for a two-sample t test:")
          print("=========================================================================")
          
        }else{
          
          cat('\n')
          print("=========================================================================")
          print("OC and ASN of the sequential MSPRT for a two-sample t test:")
          print("=========================================================================")
        }
        
        print("Group 1:")
        print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
        print("Group 2:")
        print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch1.size = c(0, cumsum(batch1.size))
      batch2.size = c(0, cumsum(batch2.size))
    }
    
    
    registerDoParallel(cores = nCore)
    out.OCandASN = foreach(theta1 = theta, .combine = 'rbind') %dopar% {
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target)
      Reject.threshold = (1 - Type2.target)/Type1.target
      
      # cut-off (with sign) in fixed design one-sample t test
      signed_t.alpha = (2*(side=='right')-1)*
        qt(Type1.target, df = N1.max + N2.max -2, lower.tail = F)
      
      # required storages
      cumSS11_n = cumSS21_n = cumsum11_n = cumsum21_n = LR1_n = numeric(nReplicate)
      type2.error.AR = rep(F, nReplicate)
      N11.AR = rep(N1.max, nReplicate)
      N21.AR = rep(N2.max, nReplicate)
      not.reached.decisionH1.AR = 1:nReplicate
      
      set.seed(seed)
      for(n in 1:nAnalyses){
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          ## observations at step n
          # Group 1
          obs11_n = mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                           FUN = function(X){
                             
                             rnorm(length(not.reached.decisionH1.AR), theta1/2, 1)
                           })
          # Group 2
          obs21_n = mapply(X = 1:(batch2.size[n+1]-batch2.size[n]), 
                           FUN = function(X){
                             
                             rnorm(length(not.reached.decisionH1.AR), -theta1/2, 1)
                           })
          
          ## sum of observations until step n
          # Group 1
          cumsum11_n[not.reached.decisionH1.AR] = 
            cumsum11_n[not.reached.decisionH1.AR] + rowSums(obs11_n)
          # Group 2
          cumsum21_n[not.reached.decisionH1.AR] = 
            cumsum21_n[not.reached.decisionH1.AR] + rowSums(obs21_n)
          
          ## sum of squares of observations until step n
          # Group 1
          cumSS11_n[not.reached.decisionH1.AR] = 
            cumSS11_n[not.reached.decisionH1.AR] + rowSums(obs11_n^2)
          # Group 2
          cumSS21_n[not.reached.decisionH1.AR] = 
            cumSS21_n[not.reached.decisionH1.AR] + rowSums(obs21_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          xbar.diff1_n = cumsum11_n[not.reached.decisionH1.AR]/batch1.size[n+1] -
            cumsum21_n[not.reached.decisionH1.AR]/batch2.size[n+1]
          divisor.pooled.sd1_n.sq = 
            cumSS11_n[not.reached.decisionH1.AR] - ((cumsum11_n[not.reached.decisionH1.AR])^2)/batch1.size[n+1] +
            cumSS21_n[not.reached.decisionH1.AR] - ((cumsum21_n[not.reached.decisionH1.AR])^2)/batch2.size[n+1]
          
          # likelihood ratio of observations until step n
          LR1_n[not.reached.decisionH1.AR] = 
            ((1 + ((xbar.diff1_n - theta0)^2)/(divisor.pooled.sd1_n.sq*
                                                 (1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff1_n - 
                        (theta0 + signed_t.alpha*
                           sqrt((divisor.pooled.sd1_n.sq/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd1_n.sq*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          # comparing with the thresholds
          AcceptedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]<=Accept.threshold)
          RejectedH0.underH1_n.AR = which(LR1_n[not.reached.decisionH1.AR]>=Reject.threshold)
          reached.decisionH1_n.AR = union(AcceptedH0.underH1_n.AR, RejectedH0.underH1_n.AR)
          
          # tracking those reaching/not reaching a decision at step n
          if(length(reached.decisionH1_n.AR)>0){
            
            N11.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch1.size[n+1]
            N21.AR[not.reached.decisionH1.AR[reached.decisionH1_n.AR]] = batch2.size[n+1]
            type2.error.AR[not.reached.decisionH1.AR[AcceptedH0.underH1_n.AR]] = T
            not.reached.decisionH1.AR = not.reached.decisionH1.AR[-reached.decisionH1_n.AR]
          }
        }
      }
      
      # attained Type II error probability
      actual.type2.error.AR = mean(type2.error.AR) +
        sum(LR1_n[not.reached.decisionH1.AR]<termination.threshold)/nReplicate
      
      # Expected sample sizes
      EN11 = mean(N11.AR)
      EN21 = mean(N21.AR)
      
      c(theta1, actual.type2.error.AR, EN11, EN21)
    }
    
    if(length(theta)==1) out.OCandASN = matrix(data = out.OCandASN, nrow = 1,
                                               ncol = 4, byrow = T)
    
    out.OCandASN = as.data.frame(out.OCandASN)
    colnames(out.OCandASN) = c('theta', 'acceptH0.prob', 'EN1', 'EN2')
    
    # msg
    if(verbose==T){
      cat('\n')
      print('Done.')
      print("-------------------------------------------------------------------------")
      cat('\n\n')
      print("=========================================================================")
      print("Performance summary:")
      print("=========================================================================")
      print(paste('Parameter value(s): ', paste(round(theta, 3), collapse = ', '), sep = ''))
      print(paste('Probability of accepting H0: ',
                  paste(round(out.OCandASN$acceptH0.prob, 3), collapse = ', '), sep = ''))
      print("Expected sample size:")
      print(paste(' Group 1 - ', paste(round(out.OCandASN$EN1, 2), collapse = ', '), sep = ''))
      print(paste(' Group 2 - ', paste(round(out.OCandASN$EN2, 2), collapse = ', '), sep = ''))
      print("=========================================================================")
      cat('\n')
    }
    
    return(out.OCandASN)
    
    # end one-sided twoT
  }else{
    
    #################### two-sample t (both sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch1.size = design.MSPRT.object$batch1.size
      batch2.size = design.MSPRT.object$batch2.size
      N1.max = design.MSPRT.object$N1.max
      N2.max = design.MSPRT.object$N2.max
      nAnalyses = design.MSPRT.object$nAnalyses
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      termination.threshold = design.MSPRT.object$termination.threshold
      nReplicate = design.MSPRT.object$nReplicate
      
      # msg
      if(verbose){
        
        if((batch1.size[1]>2)||any(batch1.size[-1]>1)||
           (batch2.size[1]>2)||any(batch2.size[-1]>1)){
          
          cat('\n')
          print("=========================================================================")
          print("OC and ASN of the group sequential MSPRT for a two-sample t test:")
          print("=========================================================================")
          
        }else{
          
          cat('\n')
          print("=========================================================================")
          print("OC and ASN of the sequential MSPRT for a two-sample t test:")
          print("=========================================================================")
        }
        
        print("Group 1:")
        print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
        print("Group 2:")
        print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch1.size = c(0, cumsum(batch1.size))
      batch2.size = c(0, cumsum(batch2.size))
      
    }else{
      
      ## ignoring batch.size
      if(!missing(batch.size)) print("'batch.size' is ignored. Not required in two-sample tests.")
      
      ## ignoring N.max
      if(!missing(N.max)) print("'N.max' is ignored. Not required in two-sample tests.")
      
      ## checking if length(batch1.size) and length(batch2.size) are equal
      if((!missing(batch1.size)) && (!missing(batch2.size)) &&
         (length(batch1.size)!=length(batch2.size))) return("Lenghts of batch1.size and batch2.size should be same")
      
      ## batch sizes and N for group 1
      if(missing(batch1.size)){
        
        if(missing(N1.max)){
          
          return("Either 'batch1.size' or 'N1.max' needs to be specified")
          
        }else{batch1.size = c(2, rep(1, N1.max-2))}
        
      }else{
        
        if(batch1.size[1]<2){
          
          return("First batch size in Group 1 should be at least 2")
          
        }else{
          
          if(missing(N1.max)){
            
            N1.max = sum(batch1.size)
            
          }else{
            
            if(sum(batch1.size)!=N1.max) return("Sum of batch1.size should add up to N1.max")
          }
        }
      }
      
      ## batch sizes and N for group 2
      if(missing(batch2.size)){
        
        if(missing(N2.max)){
          
          return("Either 'batch2.size' or 'N2.max' needs to be specified")
          
        }else{batch2.size = c(2, rep(1, N2.max-2))}
        
      }else{
        
        if(batch2.size[1]<2){
          
          return("First batch size in Group 2 should be at least 2")
          
        }else{
          
          if(missing(N2.max)){
            
            N2.max = sum(batch2.size)
            
          }else{
            
            if(sum(batch2.size)!=N2.max) return("Sum of batch2.size should add up to N2.max")
          }
        }
      }
      
      nAnalyses = length(batch1.size)
      
      # msg
      if(verbose){
        
        if((batch1.size[1]>2)||any(batch1.size[-1]>1)||
           (batch2.size[1]>2)||any(batch2.size[-1]>1)){
          
          cat('\n')
          print("=========================================================================")
          print("OC and ASN of the group sequential MSPRT for a two-sample t test:")
          print("=========================================================================")
          
        }else{
          
          cat('\n')
          print("=========================================================================")
          print("OC and ASN of the sequential MSPRT for a two-sample t test:")
          print("=========================================================================")
        }
        
        print("Group 1:")
        print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
        print("Group 2:")
        print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
        print(paste("Total number of sequential analyses: ", nAnalyses, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste('Parameter value(s) where OC and ASN is desired: ',
                    paste(round(theta, 3), collapse = ', '), sep = ''))
        print(paste("Termination threshold: ", termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
        print("Calculating the OC and ASN ...")
      }
      
      batch1.size = c(0, cumsum(batch1.size))
      batch2.size = c(0, cumsum(batch2.size))
    }
    
    
    registerDoParallel(cores = nCore)
    out.OCandASN = foreach(theta1 = theta, .combine = 'rbind') %dopar% {
      
      #### simulating data, calculating likelihood ratio, and finding the termination threshold ####
      
      # Wald's thresholds
      Accept.threshold = Type2.target/(1 - Type1.target/2)
      Reject.threshold = (1 - Type2.target)/(Type1.target/2)
      
      # cut-off (with sign) in fixed design one-sample t test
      t.alpha = qt(Type1.target/2, df = N1.max + N2.max -2, lower.tail = F)
      
      # required storages
      cumSS11_n = cumSS21_n = cumsum11_n = cumsum21_n = 
        LR1_n.r = LR1_n.l = numeric(nReplicate)
      PowerH1.AR = rep(F, nReplicate)
      N11.AR = N11.AR.r = N11.AR.l = rep(N1.max, nReplicate)
      N21.AR = N21.AR.r = N21.AR.l = rep(N2.max, nReplicate)
      decision.underH1.AR.r = decision.underH1.AR.l = rep(NA, nReplicate)
      not.reached.decisionH1.AR = not.reached.decisionH1.AR.r = not.reached.decisionH1.AR.l =
        1:nReplicate
      
      set.seed(seed)
      for(n in 1:nAnalyses){
        
        ## under H1
        if(length(not.reached.decisionH1.AR)>0){
          
          ## observations at step n
          # Group 1
          if(length(not.reached.decisionH1.AR)>1){
            
            obs11_n = mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                             FUN = function(X){
                               
                               rnorm(length(not.reached.decisionH1.AR), theta1/2, 1)
                             })
            
          }else{
            
            obs11_n = matrix(mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                                    FUN = function(X){
                                      
                                      rnorm(length(not.reached.decisionH1.AR), theta1/2, 1)
                                      
                                    }), nrow = 1, ncol = batch1.size[n+1]-batch1.size[n], 
                             byrow = T)
          }
          
          # Group 2
          if(length(not.reached.decisionH1.AR)>1){
            
            obs21_n = mapply(X = 1:(batch2.size[n+1]-batch2.size[n]), 
                             FUN = function(X){
                               
                               rnorm(length(not.reached.decisionH1.AR), -theta1/2, 1)
                             })
            
          }else{
            
            obs21_n = matrix(mapply(X = 1:(batch1.size[n+1]-batch1.size[n]), 
                                    FUN = function(X){
                                      
                                      rnorm(length(not.reached.decisionH1.AR), -theta1/2, 1)
                                      
                                    }), nrow = 1, ncol = batch1.size[n+1]-batch1.size[n], 
                             byrow = T)
          }
          
          ## sum of observations until step n
          # Group 1
          cumsum11_n[not.reached.decisionH1.AR] = 
            cumsum11_n[not.reached.decisionH1.AR] + rowSums(obs11_n)
          # Group 2
          cumsum21_n[not.reached.decisionH1.AR] = 
            cumsum21_n[not.reached.decisionH1.AR] + rowSums(obs21_n)
          
          ## sum of squares of observations until step n
          # Group 1
          cumSS11_n[not.reached.decisionH1.AR] = 
            cumSS11_n[not.reached.decisionH1.AR] + rowSums(obs11_n^2)
          # Group 2
          cumSS21_n[not.reached.decisionH1.AR] = 
            cumSS21_n[not.reached.decisionH1.AR] + rowSums(obs21_n^2)
          
          ## xbar and (n-1)*(s^2) until step n
          # for right sided check
          xbar.diff1_n.r = cumsum11_n[not.reached.decisionH1.AR.r]/batch1.size[n+1] -
            cumsum21_n[not.reached.decisionH1.AR.r]/batch2.size[n+1]
          divisor.pooled.sd1_n.sq.r = 
            cumSS11_n[not.reached.decisionH1.AR.r] - 
            ((cumsum11_n[not.reached.decisionH1.AR.r])^2)/batch1.size[n+1] +
            cumSS21_n[not.reached.decisionH1.AR.r] - 
            ((cumsum21_n[not.reached.decisionH1.AR.r])^2)/batch2.size[n+1]
          
          # for left sided check
          xbar.diff1_n.l = cumsum11_n[not.reached.decisionH1.AR.l]/batch1.size[n+1] -
            cumsum21_n[not.reached.decisionH1.AR.l]/batch2.size[n+1]
          divisor.pooled.sd1_n.sq.l = 
            cumSS11_n[not.reached.decisionH1.AR.l] -
            ((cumsum11_n[not.reached.decisionH1.AR.l])^2)/batch1.size[n+1] +
            cumSS21_n[not.reached.decisionH1.AR.l] - 
            ((cumsum21_n[not.reached.decisionH1.AR.l])^2)/batch2.size[n+1]
          
          ## likelihood ratio of observations until step n
          # for right sided check
          LR1_n.r[not.reached.decisionH1.AR.r] = 
            ((1 + ((xbar.diff1_n.r - theta0)^2)/
                (divisor.pooled.sd1_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff1_n.r - 
                        (theta0 + t.alpha*
                           sqrt((divisor.pooled.sd1_n.sq.r/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd1_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          # for left sided check
          LR1_n.l[not.reached.decisionH1.AR.l] = 
            ((1 + ((xbar.diff1_n.l - theta0)^2)/
                (divisor.pooled.sd1_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff1_n.l - 
                        (theta0 - t.alpha*
                           sqrt((divisor.pooled.sd1_n.sq.l/(batch1.size[n+1] + batch2.size[n+1] -2))*
                                  (1/N1.max + 1/N2.max))))^2)/
                  (divisor.pooled.sd1_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          ### comparing with the thresholds
          ## for right sided check
          AcceptedH0.underH1_n.AR.r = LR1_n.r[not.reached.decisionH1.AR.r]<=Accept.threshold
          RejectedH0.underH1_n.AR.r = LR1_n.r[not.reached.decisionH1.AR.r]>=Reject.threshold
          reached.decisionH1_n.AR.r = AcceptedH0.underH1_n.AR.r|RejectedH0.underH1_n.AR.r
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1_n.AR.r)){
            
            decision.underH1.AR.r[not.reached.decisionH1.AR.r[AcceptedH0.underH1_n.AR.r]] = 'A'
            decision.underH1.AR.r[not.reached.decisionH1.AR.r[RejectedH0.underH1_n.AR.r]] = 'R'
            N11.AR.r[not.reached.decisionH1.AR.r[reached.decisionH1_n.AR.r]] = batch1.size[n+1]
            N21.AR.r[not.reached.decisionH1.AR.r[reached.decisionH1_n.AR.r]] = batch2.size[n+1]
            not.reached.decisionH1.AR.r = not.reached.decisionH1.AR.r[!reached.decisionH1_n.AR.r]
          }
          
          ## for left sided check
          AcceptedH0.underH1_n.AR.l = LR1_n.l[not.reached.decisionH1.AR.l]<=Accept.threshold
          RejectedH0.underH1_n.AR.l = LR1_n.l[not.reached.decisionH1.AR.l]>=Reject.threshold
          reached.decisionH1_n.AR.l = AcceptedH0.underH1_n.AR.l|RejectedH0.underH1_n.AR.l
          
          # tracking those reaching/not reaching a decision at step n
          if(any(reached.decisionH1_n.AR.l)){
            
            decision.underH1.AR.l[not.reached.decisionH1.AR.l[AcceptedH0.underH1_n.AR.l]] = 'A'
            decision.underH1.AR.l[not.reached.decisionH1.AR.l[RejectedH0.underH1_n.AR.l]] = 'R'
            N11.AR.l[not.reached.decisionH1.AR.l[reached.decisionH1_n.AR.l]] = batch1.size[n+1]
            N21.AR.l[not.reached.decisionH1.AR.l[reached.decisionH1_n.AR.l]] = batch2.size[n+1]
            not.reached.decisionH1.AR.l = not.reached.decisionH1.AR.l[!reached.decisionH1_n.AR.l]
          }
          
          not.reached.decisionH1.AR = union(not.reached.decisionH1.AR.r,
                                            not.reached.decisionH1.AR.l)
        }
      }
      
      
      ### both-sided checking
      # accepted or rejected ones
      accepted.by.both1 = intersect(which(decision.underH1.AR.r=='A'),
                                    which(decision.underH1.AR.l=='A'))
      onlyrejected.by.right1 = intersect(which(decision.underH1.AR.r=='R'),
                                         which(decision.underH1.AR.l!='R'))
      onlyrejected.by.left1 = intersect(which(decision.underH1.AR.r!='R'),
                                        which(decision.underH1.AR.l=='R'))
      rejected.by.both1 = intersect(which(decision.underH1.AR.r=='R'),
                                    which(decision.underH1.AR.l=='R'))
      
      ## sample sizes required
      # Group 1
      N11.AR[accepted.by.both1] = pmax(N11.AR.r[accepted.by.both1],
                                       N11.AR.l[accepted.by.both1])
      N11.AR[onlyrejected.by.right1] = N11.AR.r[onlyrejected.by.right1]
      N11.AR[onlyrejected.by.left1] = N11.AR.l[onlyrejected.by.left1]
      N11.AR[rejected.by.both1] = pmin(N11.AR.r[rejected.by.both1],
                                       N11.AR.l[rejected.by.both1])
      
      # Group 2
      N21.AR[accepted.by.both1] = pmax(N21.AR.r[accepted.by.both1],
                                       N21.AR.l[accepted.by.both1])
      N21.AR[onlyrejected.by.right1] = N21.AR.r[onlyrejected.by.right1]
      N21.AR[onlyrejected.by.left1] = N21.AR.l[onlyrejected.by.left1]
      N21.AR[rejected.by.both1] = pmin(N21.AR.r[rejected.by.both1],
                                       N21.AR.l[rejected.by.both1])
      
      # inconclusive cases after both sided checking
      onlyaccepted.by.right1 = intersect(which(decision.underH1.AR.r=='A'),
                                         which(is.na(decision.underH1.AR.l)))
      onlyaccepted.by.left1 = intersect(which(is.na(decision.underH1.AR.r)),
                                        which(decision.underH1.AR.l=='A'))
      both.inconclusive1 = intersect(which(is.na(decision.underH1.AR.r)),
                                     which(is.na(decision.underH1.AR.l)))
      all.inconclusive1 = c(onlyaccepted.by.right1, onlyaccepted.by.left1,
                            both.inconclusive1)
      nNot.reached.decisionH1.AR = length(all.inconclusive1)
      
      # Type I error probability
      PowerH1.AR[c(onlyrejected.by.right1, onlyrejected.by.left1,
                   rejected.by.both1)] = T
      
      ## attained Type II error probability
      actual.PowerH1.AR = mean(PowerH1.AR) +
        sum(c(LR1_n.r[onlyaccepted.by.left1],
              LR1_n.l[onlyaccepted.by.right1],
              pmax(LR1_n.r[both.inconclusive1], LR1_n.l[both.inconclusive1]))>=
              termination.threshold)/nReplicate
      actual.type2.errorH1.AR = 1 - actual.PowerH1.AR
      
      ## Expected sample sizes
      # Group 1
      EN11 = mean(N11.AR)
      
      # Group 2
      EN21 = mean(N21.AR)
      
      c(theta1, actual.type2.errorH1.AR, EN11, EN21)
    }
    
    if(length(theta)==1) out.OCandASN = matrix(data = out.OCandASN, nrow = 1,
                                               ncol = 4, byrow = T)
    out.OCandASN = as.data.frame(out.OCandASN)
    colnames(out.OCandASN) = c('theta', 'acceptH0.prob', 'EN1', 'EN2')
    
    # msg
    if(verbose==T){
      cat('\n')
      print('Done.')
      print("-------------------------------------------------------------------------")
      cat('\n\n')
      print("=========================================================================")
      print("Performance summary:")
      print("=========================================================================")
      print(paste('Parameter value(s): ', paste(round(theta, 3), collapse = ', '), sep = ''))
      print(paste('Probability of accepting H0: ',
                  paste(round(out.OCandASN$acceptH0.prob, 3), collapse = ', '), sep = ''))
      print("Expected sample size:")
      print(paste(' Group 1 - ', paste(round(out.OCandASN$EN1, 2), collapse = ', '), sep = ''))
      print(paste(' Group 2 - ', paste(round(out.OCandASN$EN2, 2), collapse = ', '), sep = ''))
      print("=========================================================================")
      cat('\n')
    }
    
    return(out.OCandASN)
    
  } # end both-sided twoT
  
  # end twoT
}


#### OC and ASN of the MSPRT combined for all ####
OCandASN.MSPRT = function(theta, design.MSPRT.object, 
                          termination.threshold, test.type, 
                          side = 'right', theta0, 
                          Type1.target =.005, Type2.target = .2,
                          N.max, N1.max, N2.max,
                          sigma = 1, sigma1 = 1, sigma2 = 1,
                          batch.size, batch1.size, batch2.size,
                          nReplicate = 1e+6, nCore = detectCores() - 1,
                          verbose = T, seed = 1){
  
  if(!missing(design.MSPRT.object)) test.type = design.MSPRT.object$test.type
  
  if(test.type=='oneProp'){
    
    if(!missing(design.MSPRT.object)){
      
      return(OCandASN.MSPRT_oneProp(theta = theta,
                                    design.MSPRT.object = design.MSPRT.object,
                                    nReplicate = nReplicate, nCore = nCore,
                                    verbose = verbose, seed = seed))
      
    }else{
      
      return(OCandASN.MSPRT_oneProp(theta = theta,
                                    termination.threshold = termination.threshold,
                                    side = side, theta0 = theta0, 
                                    Type1.target = Type1.target,
                                    Type2.target = Type2.target,
                                    N.max = N.max, batch.size = batch.size,
                                    nReplicate = nReplicate, nCore = nCore,
                                    verbose = verbose, seed = seed))
    }
    
  }else if(test.type=='oneZ'){
    
    if(!missing(design.MSPRT.object)){
      
      return(OCandASN.MSPRT_oneZ(theta = theta,
                                 design.MSPRT.object = design.MSPRT.object,
                                 nReplicate = nReplicate, nCore = nCore,
                                 verbose = verbose, seed = seed))
      
    }else{
      
      return(OCandASN.MSPRT_oneZ(theta = theta,
                                 termination.threshold = termination.threshold,
                                 side = side, theta0 = theta0, sigma = sigma,
                                 Type1.target = Type1.target,
                                 Type2.target = Type2.target,
                                 N.max = N.max, batch.size = batch.size,
                                 nReplicate = nReplicate, nCore = nCore,
                                 verbose = verbose, seed = seed))
    }
    
  }else if(test.type=='oneT'){
    
    if(!missing(design.MSPRT.object)){
      
      return(OCandASN.MSPRT_oneT(theta = theta,
                                 design.MSPRT.object = design.MSPRT.object,
                                 nReplicate = nReplicate, nCore = nCore,
                                 verbose = verbose, seed = seed))
      
    }else{
      
      return(OCandASN.MSPRT_oneT(theta = theta,
                                 termination.threshold = termination.threshold,
                                 side = side, theta0 = theta0,
                                 Type1.target = Type1.target,
                                 Type2.target = Type2.target,
                                 N.max = N.max, batch.size = batch.size,
                                 nReplicate = nReplicate, nCore = nCore,
                                 verbose = verbose, seed = seed))
    }
    
  }else if(test.type=='twoZ'){
    
    if(!missing(design.MSPRT.object)){
      
      return(OCandASN.MSPRT_twoZ(theta = theta,
                                 design.MSPRT.object = design.MSPRT.object,
                                 nReplicate = nReplicate, nCore = nCore,
                                 verbose = verbose, seed = seed))
      
    }else{
      
      return(OCandASN.MSPRT_twoZ(theta = theta,
                                 termination.threshold = termination.threshold,
                                 side = side, theta0 = theta0,
                                 Type1.target = Type1.target,
                                 Type2.target = Type2.target,
                                 N1.max = N1.max, N2.max = N2.max,
                                 sigma1 = sigma1, sigma2 = sigma2,
                                 batch1.size = batch1.size,
                                 batch2.size = batch2.size,
                                 nReplicate = nReplicate, nCore = nCore,
                                 verbose = verbose, seed = seed))
    }
    
  }else if(test.type=='twoT'){
    
    if(!missing(design.MSPRT.object)){
      
      return(OCandASN.MSPRT_twoT(theta = theta,
                                 design.MSPRT.object = design.MSPRT.object,
                                 nReplicate = nReplicate, nCore = nCore,
                                 verbose = verbose, seed = seed))
      
    }else{
      
      return(OCandASN.MSPRT_twoT(theta = theta,
                                 termination.threshold = termination.threshold,
                                 side = side, theta0 = theta0,
                                 Type1.target = Type1.target,
                                 Type2.target = Type2.target,
                                 N1.max = N1.max, N2.max = N2.max,
                                 batch1.size = batch1.size,
                                 batch2.size = batch2.size,
                                 nReplicate = nReplicate, nCore = nCore,
                                 verbose = verbose, seed = seed))
    }
  }
}



################# implementing the MSPRT #################

#### one-sample proportion test ####
implement.MSPRT_oneProp = function(obs, design.MSPRT.object, 
                                   termination.threshold,
                                   side = 'right', theta0 = 0.5, 
                                   Type1.target =.005, Type2.target = .2,
                                   N.max, batch.size,
                                   verbose = T, plot.it = T){
  
  # side
  if(!missing(design.MSPRT.object)) side = design.MSPRT.object$side
  
  if(side!='both'){
    
    #################### one-sample proportion (right/left sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch.size = design.MSPRT.object$batch.size
      N.max = design.MSPRT.object$N.max
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      termination.threshold = design.MSPRT.object$termination.threshold
      UMPBT = design.MSPRT.object$UMPBT
      nAnalyses.max = design.MSPRT.object$nAnalyses
      
      # msg
      if(verbose){
        
        if(any(batch.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the group sequential MSPRT for a one-sample proportion test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the sequential MSPRT for a one-sample proportion test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(UMPBT$theta[1], 3), " & ",
                    round(UMPBT$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$mix.prob[1], 3), " & ", 1 - round(UMPBT$mix.prob[1], 3), sep = ''))
        print("-------------------------------------------------------------------------")
      }
      
      batch.size = c(0, cumsum(batch.size))
      nAnalyses = max(which(batch.size<=length(obs))) - 1
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(obs1)) print("'obs1' is ignored. Not required in one-sample tests.")
      if(!missing(obs2)) print("'obs2' is ignored. Not required in one-sample tests.")
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
      if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
      
      ## ignoring N1.max & N2.max
      if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
      if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
      
      ## batch sizes and N.max
      if(missing(batch.size)){
        
        if(missing(N.max)){
          
          return("Either 'batch.size' or 'N.max' needs to be specified")
          
        }else{batch.size = rep(1, N.max)}
        
      }else{
        
        if(missing(N.max)){
          
          N.max = sum(batch.size)
          
        }else{
          
          if(sum(batch.size)!=N.max) return("Sum of batch sizes should add up to N.max")
        }
      }
      
      nAnalyses.max = length(batch.size)
      
      ## point H0
      if(missing(theta0)) theta0 = 0.5
      
      ######################## UMPBT alternative ########################
      UMPBT = UMPBT.alt(test.type = 'oneProp', side = side, theta0 = theta0,
                        N = N.max, Type1 = Type1.target)
      
      # msg
      if(verbose){
        
        if(any(batch.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the group sequential MSPRT for a one-sample proportion test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the sequential MSPRT for a one-sample proportion test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(UMPBT$theta[1], 3), " & ",
                    round(UMPBT$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$mix.prob[1], 3), " & ", 1 - round(UMPBT$mix.prob[1], 3), sep = ''))
        print("-------------------------------------------------------------------------")
      }
      
      batch.size = c(0, cumsum(batch.size))
      nAnalyses = max(which(batch.size<=length(obs))) - 1
    }
    
    
    ###################### sequential comparison ######################
    # Wald's thresholds
    Accept.threshold = Type2.target/(1 - Type1.target)
    Reject.threshold = (1 - Type2.target)/Type1.target
    
    # required storages
    cumsum_n = 0
    reached.decision = F
    rejectH0 = NA
    LR = rep(NA, nAnalyses)
    
    for(n in 1:nAnalyses){
      
      if(!reached.decision){
        
        # sum of observations until step n
        cumsum_n = cumsum_n + sum(obs[(batch.size[n]+1):batch.size[n+1]])
        
        # likelihood ratio of observations until step n
        LR[n] = 
          UMPBT$mix.prob[1]*(((1 - UMPBT$theta[1])/(1 - theta0))^batch.size[n+1])*
          ((UMPBT$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$theta[1])))^cumsum_n +
          (1 - UMPBT$mix.prob[2])*(((1 - UMPBT$theta[2])/(1 - theta0))^batch.size[n+1])*
          ((UMPBT$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$theta[2])))^cumsum_n
        
        # comparing with the thresholds
        AcceptedH0_n = LR[n]<=Accept.threshold
        RejectedH0_n = LR[n]>=Reject.threshold
        reached.decision = AcceptedH0_n||RejectedH0_n
        if(reached.decision){
          
          n0 = batch.size[n+1]
          rejectH0 = RejectedH0_n
          if(rejectH0){
            
            decision = 'reject'
            
          }else{decision = 'accept'}
        }
      }
    }
    
    # inconclusive cases
    if(!reached.decision){
      
      if(nAnalyses==nAnalyses.max){
        
        n0 = N.max
        rejectH0 = LR[nAnalyses]>=termination.threshold
        
        if(rejectH0){
          
          decision = 'reject'
          
        }else if(!rejectH0){decision = 'accept'}
        reached.decision = T
        
      }else{
        
        n0 = batch.size[nAnalyses+1]
        decision = 'continue'
      }
    }
    
    # msg
    if(verbose==T){
      
      if(decision=='continue'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Continue sampling')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='reject'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Reject the null hypothesis')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='accept'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Accept the null hypothesis')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
      }
    }
    
    
    nAnalyses.test = max(which(!is.na(LR)))
    
    ## plotting
    if(plot.it==T){
      
      # title and decision
      if(side=='right'){
        
        testname = 'Right-sided one-sample proportion test'
        
      }else{testname = 'Left-sided one-sample proportion test'}
      
      if(decision=="accept"){
        
        plot.subtitle = 
          paste('Accept the null (n = ', n0, ')', sep = '')
        
        min.LR = min(LR, na.rm = T)
        if(min.LR<Accept.threshold){
          
          ylow = min.LR
          yup = max(LR, na.rm = T)
          
        }else if((min.LR>=Accept.threshold)&&
                 (min.LR<(termination.threshold + Accept.threshold)/2)){
          
          ylow = Accept.threshold
          yup = termination.threshold
          
        }else{
          
          ylow = Accept.threshold
          yup = Reject.threshold
        }
        
      }else if(decision=="reject"){
        
        plot.subtitle = 
          paste('Reject the null (n = ', n0, ')', sep = '')
        
        ylow = Accept.threshold
        yup = max(LR, na.rm = T)
        
      }else if(decision=="continue"){
        
        plot.subtitle = 
          paste('Continue sampling (n = ', n0, ')', sep = '')
        
        last.LR = LR[max(which(!is.na(LR)))]
        if(last.LR<(termination.threshold + Accept.threshold)/2){
          
          ylow = Accept.threshold
          yup = termination.threshold
          
        }else{
          
          ylow = Accept.threshold
          yup = Reject.threshold
        }
      }
      
      df = rbind.data.frame(data.frame('xval' = 1:nAnalyses.test,
                                       'yval' = Accept.threshold,
                                       'type' = 'A'),
                            data.frame('xval' = 1:nAnalyses.test,
                                       'yval' = Reject.threshold,
                                       'type' = 'R'),
                            data.frame('xval' = 1:nAnalyses.test,
                                       'yval' = LR[1:nAnalyses.test],
                                       'type' = 'LR'),
                            data.frame('xval' = design.MSPRT.object$nAnalyses,
                                       'yval' = termination.threshold,
                                       'type' = 'term.thresh'))
      
      df$type = factor(as.character(df$type),
                       levels = c('A', 'R', 'LR', 'term.thresh'))
      
      seqcompare = ggplot(data = df,
                          aes(x = xval, y = yval, group = type)) + 
        geom_point(aes(colour = type), size = 2) +
        geom_line(aes(colour = type), size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Accept.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="forestgreen", size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Reject.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="red2", size=1) +
        geom_point(aes(x = length(batch.size) - 1, y = termination.threshold),
                   colour = "black", size = 2) +
        xlim(c(1, nAnalyses.test)) +
        ylim(c(ylow, yup)) +
        scale_color_manual(labels = c('Acceptance threshold', 'Rejection threshold',
                                      'Wtd. likelihood ratio', 'Termination threshold'),
                           values = c('A' = 'forestgreen', 'R' = 'red2',
                                      'LR' = 'dodgerblue', 'term.thresh' = 'black'),
                           guide = guide_legend(nrow = 2, 
                                                override.aes = list(linetype=c(1,1,1,0), 
                                                                    shape=16))) +
        theme(plot.title = element_text(size=22, face="bold"),
              plot.subtitle = element_text(size=22, face="bold"),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              axis.text.x = element_text(color = "black", size = 22),
              axis.text.y = element_text(color = "black", size = 22),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1.2),
              panel.background = element_blank(), legend.title=element_blank(),
              legend.key.width = unit(1.4, "cm"),
              legend.spacing.x = unit(1, 'cm'), legend.text=element_text(size=22),
              legend.position = 'bottom') +
        labs(title = testname, subtitle = plot.subtitle,
             y = 'Wtd. likelihood ratio',
             x = 'Steps in sequential analyses')
      
      suppressWarnings(print(seqcompare))
    }
    
    return(list('n' = n0, 'decision' = decision, 
                'Accept.threshold' = Accept.threshold, 'Reject.threshold' = Reject.threshold,
                'LR' = LR[1:nAnalyses.test], 'UMPBT' = UMPBT))
    
    # end one-sided oneProp
    
  }else{
    
    #################### one-sample proportion (both sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch.size = design.MSPRT.object$batch.size
      N.max = design.MSPRT.object$N.max
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      termination.threshold = design.MSPRT.object$termination.threshold
      UMPBT = design.MSPRT.object$UMPBT
      nAnalyses.max = design.MSPRT.object$nAnalyses
      
      # msg
      if(verbose){
        
        if(any(batch.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the group sequential MSPRT for a one-sample proportion test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the sequential MSPRT for a one-sample proportion test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Maximum number of sequential analyses: ", design.MSPRT.object$nAnalyses,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right:, ', round(UMPBT$right$theta[1], 3), " & ",
                    round(UMPBT$right$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$right$mix.prob[1], 3), " & ", 1 - round(UMPBT$right$mix.prob[1], 3),
                    sep = ""))
        print(paste(' On the left:, ', round(UMPBT$left$theta[1], 3), " & ",
                    round(UMPBT$left$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$left$mix.prob[1], 3), " & ", 1 - round(UMPBT$left$mix.prob[1], 3),
                    sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
      batch.size = c(0, cumsum(batch.size))
      nAnalyses = max(which(batch.size<=length(obs))) - 1
      
    }else{
      
      ## ignoring obs1 & obs2
      if(!missing(obs1)) print("'obs1' is ignored. Not required in one-sample tests.")
      if(!missing(obs2)) print("'obs2' is ignored. Not required in one-sample tests.")
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
      if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
      
      ## ignoring N1.max & N2.max
      if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
      if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
      
      ## batch sizes and N.max
      if(missing(batch.size)){
        
        if(missing(N.max)){
          
          return("Either 'batch.size' or 'N.max' needs to be specified")
          
        }else{batch.size = rep(1, N.max)}
        
      }else{
        
        if(missing(N.max)){
          
          N.max = sum(batch.size)
          
        }else{
          
          if(sum(batch.size)!=N.max) return("Sum of batch sizes should add up to N.max")
        }
      }
      
      nAnalyses.max = length(batch.size)
      
      ## point H0
      if(missing(theta0)) theta0 = 0.5
      
      ######################## UMPBT alternative ########################
      UMPBT = list('right' = UMPBT.alt(test.type = 'oneProp', side = 'right', 
                                       theta0 = theta0, N = N.max, Type1 = Type1.target/2),
                   'left' = UMPBT.alt(test.type = 'oneProp', side = 'left',
                                      theta0 = theta0, N = N.max, Type1 = Type1.target/2))
      
      # msg
      if(verbose){
        
        if(any(batch.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the group sequential MSPRT for a one-sample proportion test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the sequential MSPRT for a one-sample proportion test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(UMPBT$right$theta[1], 3), " & ",
                    round(UMPBT$right$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$right$mix.prob[1], 3), " & ", 1 - round(UMPBT$right$mix.prob[1], 3),
                    sep = ""))
        print(paste(' On the left: ', round(UMPBT$left$theta[1], 3), " & ",
                    round(UMPBT$left$theta[2], 3), " with respective probabilities ",
                    round(UMPBT$left$mix.prob[1], 3), " & ", 1 - round(UMPBT$left$mix.prob[1], 3),
                    sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
      batch.size = c(0, cumsum(batch.size))
      nAnalyses = max(which(batch.size<=length(obs))) - 1
    }
    
    
    ###################### sequential comparison ######################
    # Wald's thresholds
    Accept.threshold = Type2.target/(1 - Type1.target/2)
    Reject.threshold = (1 - Type2.target)/(Type1.target/2)
    
    # required storages
    cumsum_n = 0
    reached.decision.r = reached.decision.l = reached.decision = F
    rejectH0 = NA
    LR.r = LR.l = rep(NA, nAnalyses)
    
    for(n in 1:nAnalyses){
      
      if(!reached.decision){
        
        # sum of observations until step n
        cumsum_n = cumsum_n + sum(obs[(batch.size[n]+1):batch.size[n+1]])
        
        ## for right sided check
        if(!reached.decision.r){
          
          # likelihood ratio of observations until step n
          LR.r[n] = 
            UMPBT$right$mix.prob[1]*(((1 - UMPBT$right$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[1])))^cumsum_n +
            (1 - UMPBT$right$mix.prob[2])*(((1 - UMPBT$right$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$right$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$right$theta[2])))^cumsum_n
          
          # comparing with the thresholds
          AcceptedH0_n.r = LR.r[n]<=Accept.threshold
          RejectedH0_n.r = LR.r[n]>=Reject.threshold
          reached.decision.r = AcceptedH0_n.r||RejectedH0_n.r
        }
        
        ## for left sided check
        if(!reached.decision.l){
          
          # likelihood ratio of observations until step n
          LR.l[n] = 
            UMPBT$left$mix.prob[1]*(((1 - UMPBT$left$theta[1])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[1]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[1])))^cumsum_n +
            (1 - UMPBT$left$mix.prob[2])*(((1 - UMPBT$left$theta[2])/(1 - theta0))^batch.size[n+1])*
            ((UMPBT$left$theta[2]*(1 - theta0))/(theta0*(1 - UMPBT$left$theta[2])))^cumsum_n
          
          # comparing with the thresholds
          AcceptedH0_n.l = LR.l[n]<=Accept.threshold
          RejectedH0_n.l = LR.l[n]>=Reject.threshold
          reached.decision.l = AcceptedH0_n.l||RejectedH0_n.l
        }
        
        ## both-sided check
        if(AcceptedH0_n.r&&AcceptedH0_n.l){
          
          rejectH0 = F
          decision = 'accept'
          n0 = batch.size[n+1]
          reached.decision = T
          
        }else if(RejectedH0_n.r||RejectedH0_n.l){
          
          rejectH0 = T
          decision = 'reject'
          n0 = batch.size[n+1]
          reached.decision = T
        }
      }
    }
    
    # inconclusive cases
    if(!reached.decision){
      
      if(nAnalyses==nAnalyses.max){
        
        n0 = N.max
        if(AcceptedH0_n.l&&(!reached.decision.r)){
          
          rejectH0 = LR.r[nAnalyses]>=termination.threshold
          
        }else if(AcceptedH0_n.r&&(!reached.decision.l)){
          
          rejectH0 = LR.l[nAnalyses]>=termination.threshold
          
        }else if((!reached.decision.r)&&(!reached.decision.l)){
          
          rejectH0 = max(LR.r[nAnalyses], LR.l[nAnalyses])>=termination.threshold
          
        }else{rejectH0 = F}
        
        if(rejectH0){
          
          decision = 'reject'
          
        }else if(!rejectH0){decision = 'accept'}
        reached.decision = T
        
      }else{
        
        n0 = batch.size[nAnalyses + 1]
        decision = 'continue'
      }
    }
    
    # msg
    if(verbose==T){
      
      if(decision=='continue'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Continue sampling')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='reject'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Reject the null hypothesis')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='accept'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Accept the null hypothesis')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
      }
    }
    
    
    nAnalyses.r = max(which(!is.na(LR.r)))
    nAnalyses.l = max(which(!is.na(LR.l)))
    
    ## plotting
    if(plot.it==T){
      
      # title and decision
      testname = 'Two-sided one-sample proportion test'
      if(decision=="accept"){
        
        plot.subtitle = 
          paste('Accept the null (n = ', n0, ')', sep = '')
        
      }else if(decision=="reject"){
        
        plot.subtitle = 
          paste('Reject the null (n = ', n0, ')', sep = '')
        
      }else if(decision=="continue"){
        
        plot.subtitle = 
          paste('Continue sampling (n = ', n0, ')', sep = '')
      }
      
      # right sided test
      if(AcceptedH0_n.r){
        
        min.LR.r = min(LR.r, na.rm = T)
        if(min.LR.r<Accept.threshold){
          
          ylow.r = min.LR.r
          yup.r = max(LR.r, na.rm = T)
          
        }else if((min.LR.r>=Accept.threshold)&&
                 (min.LR.r<(termination.threshold + Accept.threshold)/2)){
          
          ylow.r = Accept.threshold
          yup.r = termination.threshold
          
        }else{
          
          ylow.r = Accept.threshold
          yup.r = Reject.threshold
        }
        
      }else if(RejectedH0_n.r){
        
        ylow.r = Accept.threshold
        yup.r = max(LR.r, na.rm = T)
        
      }else if(!reached.decision.r){
        
        last.LR.r = LR.r[max(which(!is.na(LR.r)))]
        if(last.LR.r<(termination.threshold + Accept.threshold)/2){
          
          ylow.r = Accept.threshold
          yup.r = termination.threshold
          
        }else{
          
          ylow.r = Accept.threshold
          yup.r = Reject.threshold
        }
      }
      
      df.r = rbind.data.frame(data.frame('xval' = 1:nAnalyses.r,
                                         'yval' = Accept.threshold,
                                         'type' = 'A'),
                              data.frame('xval' = 1:nAnalyses.r,
                                         'yval' = Reject.threshold,
                                         'type' = 'R'),
                              data.frame('xval' = 1:nAnalyses.r,
                                         'yval' = LR.r[1:nAnalyses.r],
                                         'type' = 'LR'),
                              data.frame('xval' = design.MSPRT.object$nAnalyses,
                                         'yval' = termination.threshold,
                                         'type' = 'term.thresh'))
      
      df.r$type = factor(as.character(df.r$type),
                         levels = c('A', 'R', 'LR', 'term.thresh'))
      
      seqcompare.r = ggplot(data = df.r,
                            aes(x = xval, y = yval, group = type)) + 
        geom_point(aes(colour = type), size = 2) +
        geom_line(aes(colour = type), size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Accept.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="forestgreen", size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Reject.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="red2", size=1) +
        geom_point(aes(x = length(batch.size) - 1, y = termination.threshold),
                   colour = "black", size = 2) +
        xlim(c(1, nAnalyses.r)) +
        ylim(c(ylow.r, yup.r)) +
        scale_color_manual(labels = c('Acceptance threshold', 'Rejection threshold',
                                      'Wtd. likelihood ratio', 'Termination threshold'),
                           values = c('A' = 'forestgreen', 'R' = 'red2',
                                      'LR' = 'dodgerblue', 'term.thresh' = 'black'),
                           guide = guide_legend(nrow = 2,
                                                override.aes = list(linetype=c(1,1,1,0), 
                                                                    shape=16))) +
        theme(plot.title = element_text(size=22, face="bold"),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              axis.text.x = element_text(color = "black", size = 22),
              axis.text.y = element_text(color = "black", size = 22),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1.2),
              panel.background = element_blank(), legend.title=element_blank(),
              legend.key.width = unit(1.4, "cm"),
              legend.spacing.x = unit(1, 'cm'), legend.text=element_text(size=22),
              legend.position = 'bottom') +
        labs(title = paste('Right-sided test at level ', Type1.target/2, sep = ''),
             y = 'Wtd. likelihood ratio',
             x = 'Steps in sequential analyses')
      
      # left sided test
      if(AcceptedH0_n.l){
        
        min.LR.l = min(LR.l, na.rm = T)
        if(min.LR.l<Accept.threshold){
          
          ylow.l = min.LR.l
          yup.l = max(LR.l, na.rm = T)
          
        }else if((min.LR.l>=Accept.threshold)&&
                 (min.LR.l<(termination.threshold + Accept.threshold)/2)){
          
          ylow.l = Accept.threshold
          yup.l = termination.threshold
          
        }else{
          
          ylow.l = Accept.threshold
          yup.l = Reject.threshold
        }
        
      }else if(RejectedH0_n.l){
        
        ylow.l = Accept.threshold
        yup.l = max(LR.l, na.rm = T)
        
      }else if(!reached.decision.l){
        
        last.LR.l = LR.l[max(which(!is.na(LR.l)))]
        if(last.LR.l<(termination.threshold + Accept.threshold)/2){
          
          ylow.l = Accept.threshold
          yup.l = termination.threshold
          
        }else{
          
          ylow.l = Accept.threshold
          yup.l = Reject.threshold
        }
      }
      
      df.l = rbind.data.frame(data.frame('xval' = 1:nAnalyses.l,
                                         'yval' = Accept.threshold,
                                         'type' = 'A'),
                              data.frame('xval' = 1:nAnalyses.l,
                                         'yval' = Reject.threshold,
                                         'type' = 'R'),
                              data.frame('xval' = 1:nAnalyses.l,
                                         'yval' = LR.l[1:nAnalyses.l],
                                         'type' = 'LR'),
                              data.frame('xval' = design.MSPRT.object$nAnalyses,
                                         'yval' = termination.threshold,
                                         'type' = 'term.thresh'))
      
      df.l$type = factor(as.character(df.l$type),
                         levels = c('A', 'R', 'LR', 'term.thresh'))
      
      seqcompare.l = ggplot(data = df.l,
                            aes(x = xval, y = yval, group = type)) + 
        geom_point(aes(colour = type), size = 2) +
        geom_line(aes(colour = type), size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Accept.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="forestgreen", size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Reject.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="red2", size=1) +
        geom_point(aes(x = length(batch.size) - 1, y = termination.threshold),
                   colour = "black", size = 2) +
        xlim(c(1, nAnalyses.l)) +
        ylim(c(ylow.l, yup.l)) +
        scale_color_manual(labels = c('Acceptance threshold', 'Rejection threshold',
                                      'Wtd. likelihood ratio', 'Termination threshold'),
                           values = c('A' = 'forestgreen', 'R' = 'red2',
                                      'LR' = 'dodgerblue', 'term.thresh' = 'black'),
                           guide = guide_legend(nrow = 2)) +
        theme(plot.title = element_text(size=22, face="bold"),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              axis.text.x = element_text(color = "black", size = 22),
              axis.text.y = element_text(color = "black", size = 22),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1.2),
              panel.background = element_blank(), legend.title=element_blank(),
              legend.key.width = unit(1.4, "cm"),
              legend.spacing.x = unit(1, 'cm'), legend.text=element_text(size=22),
              legend.position = 'bottom') +
        labs(title = paste('Left-sided test at level ', Type1.target/2, sep = ''),
             y = 'Wtd. likelihood ratio',
             x = 'Steps in sequential analyses')
      
      suppressWarnings(print(annotate_figure(ggarrange(seqcompare.l, seqcompare.r, 
                                                       nrow = 1, ncol = 2,
                                                       legend = 'bottom', common.legend = T),
                                             top = text_grob(paste(testname, ': ', plot.subtitle, '\n'), face = "bold", 
                                                             size = 25, hjust = .5))))
    }
    
    return(list('n' = n0, 'decision' = decision, 
                'Accept.threshold' = Accept.threshold, 'Reject.threshold' = Reject.threshold,
                'LR' = list('right' = LR.r[1:nAnalyses.r], 'left' = LR.l[1:nAnalyses.l]),
                'UMPBT' = UMPBT))
    
  } # end both-sided oneProp
}

#### one-sample z test ####
implement.MSPRT_oneZ = function(obs, design.MSPRT.object, 
                                termination.threshold,
                                side = 'right', theta0 = 0, 
                                Type1.target =.005, Type2.target = .2,
                                N.max, sigma = 1, batch.size,
                                verbose = T, plot.it = T){
  
  # side
  if(!missing(design.MSPRT.object)) side = design.MSPRT.object$side
  
  if(side!='both'){
    
    #################### one-sample z (right/left sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch.size = design.MSPRT.object$batch.size
      N.max = design.MSPRT.object$N.max
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      sigma = design.MSPRT.object$sigma
      termination.threshold = design.MSPRT.object$termination.threshold
      theta.UMPBT = design.MSPRT.object$theta.UMPBT
      nAnalyses.max = design.MSPRT.object$nAnalyses
      
      # msg
      if(verbose){
        
        if(any(batch.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the group sequential MSPRT for a one-sample z test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the sequential MSPRT for a one-sample z test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Known standard deviation: ", sigma, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(theta.UMPBT, 3)))
        print("-------------------------------------------------------------------------")
      }
      
      batch.size = c(0, cumsum(batch.size))
      nAnalyses = max(which(batch.size<=length(obs))) - 1
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(obs1)) print("'obs1' is ignored. Not required in one-sample tests.")
      if(!missing(obs2)) print("'obs2' is ignored. Not required in one-sample tests.")
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
      if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
      
      ## ignoring N1.max & N2.max
      if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
      if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
      
      ## batch sizes and N.max
      if(missing(batch.size)){
        
        if(missing(N.max)){
          
          return("Either 'batch.size' or 'N.max' needs to be specified")
          
        }else{batch.size = rep(1, N.max)}
        
      }else{
        
        if(missing(N.max)){
          
          N.max = sum(batch.size)
          
        }else{
          
          if(sum(batch.size)!=N.max) return("Sum of batch sizes should add up to N.max")
        }
      }
      
      nAnalyses.max = length(batch.size)
      
      ## point H0
      if(missing(theta0)) theta0 = 0
      
      ######################## UMPBT alternative ########################
      theta.UMPBT = UMPBT.alt(test.type = 'oneZ', side = side, theta0 = theta0,
                              N = N.max, Type1 = Type1.target, sigma = sigma)
      
      # msg
      if(verbose){
        
        if(any(batch.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the group sequential MSPRT for a one-sample z test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the sequential MSPRT for a one-sample z test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Known standard deviation: ", sigma, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(theta.UMPBT, 3)))
        print("-------------------------------------------------------------------------")
      }
      
      batch.size = c(0, cumsum(batch.size))
      nAnalyses = max(which(batch.size<=length(obs))) - 1
    }
    
    
    ###################### sequential comparison ######################
    # Wald's thresholds
    Accept.threshold = Type2.target/(1 - Type1.target)
    Reject.threshold = (1 - Type2.target)/Type1.target
    
    # required storages
    cumsum_n = 0
    reached.decision = F
    rejectH0 = NA
    LR = rep(NA, nAnalyses)
    
    for(n in 1:nAnalyses){
      
      if(!reached.decision){
        
        # sum of observations until step n
        cumsum_n = cumsum_n + sum(obs[(batch.size[n]+1):batch.size[n+1]])
        
        # likelihood ratio of observations until step n
        LR[n] = 
          exp((cumsum_n*(theta.UMPBT - theta0) - 
                 ((batch.size[n+1]*((theta.UMPBT^2) - (theta0^2)))/2))/(sigma^2))
        
        # comparing with the thresholds
        AcceptedH0_n = LR[n]<=Accept.threshold
        RejectedH0_n = LR[n]>=Reject.threshold
        reached.decision = AcceptedH0_n||RejectedH0_n
        if(reached.decision){
          
          n0 = batch.size[n+1]
          rejectH0 = RejectedH0_n
          if(rejectH0){
            
            decision = 'reject'
            
          }else{decision = 'accept'}
        }
      }
    }
    
    # inconclusive cases
    if(!reached.decision){
      
      if(nAnalyses==nAnalyses.max){
        
        n0 = N.max
        rejectH0 = LR[nAnalyses]>=termination.threshold
        
        if(rejectH0){
          
          decision = 'reject'
          
        }else if(!rejectH0){decision = 'accept'}
        reached.decision = T
        
      }else{
        
        n0 = batch.size[nAnalyses+1]
        decision = 'continue'
      }
    }
    
    # msg
    if(verbose==T){
      
      if(decision=='continue'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Continue sampling')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='reject'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Reject the null hypothesis')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='accept'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Accept the null hypothesis')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
      }
    }
    
    
    nAnalyses.test = max(which(!is.na(LR)))
    
    ## plotting
    if(plot.it==T){
      
      # title and decision
      if(side=='right'){
        
        testname = 'Right-sided one-sample z test'
        
      }else{testname = 'Left-sided one-sample z test'}
      
      if(decision=="accept"){
        
        plot.subtitle = 
          paste('Accept the null (n = ', n0, ')', sep = '')
        
        min.LR = min(LR, na.rm = T)
        if(min.LR<Accept.threshold){
          
          ylow = min.LR
          yup = max(LR, na.rm = T)
          
        }else if((min.LR>=Accept.threshold)&&
                 (min.LR<(termination.threshold + Accept.threshold)/2)){
          
          ylow = Accept.threshold
          yup = termination.threshold
          
        }else{
          
          ylow = Accept.threshold
          yup = Reject.threshold
        }
        
      }else if(decision=="reject"){
        
        plot.subtitle = 
          paste('Reject the null (n = ', n0, ')', sep = '')
        
        ylow = Accept.threshold
        yup = max(LR, na.rm = T)
        
      }else if(decision=="continue"){
        
        plot.subtitle = 
          paste('Continue sampling (n = ', n0, ')', sep = '')
        
        last.LR = LR[max(which(!is.na(LR)))]
        if(last.LR<(termination.threshold + Accept.threshold)/2){
          
          ylow = Accept.threshold
          yup = termination.threshold
          
        }else{
          
          ylow = Accept.threshold
          yup = Reject.threshold
        }
      }
      
      df = rbind.data.frame(data.frame('xval' = 1:nAnalyses.test,
                                       'yval' = Accept.threshold,
                                       'type' = 'A'),
                            data.frame('xval' = 1:nAnalyses.test,
                                       'yval' = Reject.threshold,
                                       'type' = 'R'),
                            data.frame('xval' = 1:nAnalyses.test,
                                       'yval' = LR[1:nAnalyses.test],
                                       'type' = 'LR'),
                            data.frame('xval' = design.MSPRT.object$nAnalyses,
                                       'yval' = termination.threshold,
                                       'type' = 'term.thresh'))
      
      df$type = factor(as.character(df$type),
                       levels = c('A', 'R', 'LR', 'term.thresh'))
      
      seqcompare = ggplot(data = df,
                          aes(x = xval, y = yval, group = type)) + 
        geom_point(aes(colour = type), size = 2) +
        geom_line(aes(colour = type), size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Accept.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="forestgreen", size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Reject.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="red2", size=1) +
        geom_point(aes(x = length(batch.size) - 1, y = termination.threshold),
                   colour = "black", size = 2) +
        xlim(c(1, nAnalyses.test)) +
        ylim(c(ylow, yup)) +
        scale_color_manual(labels = c('Acceptance threshold', 'Rejection threshold',
                                      'Likelihood ratio', 'Termination threshold'),
                           values = c('A' = 'forestgreen', 'R' = 'red2',
                                      'LR' = 'dodgerblue', 'term.thresh' = 'black'),
                           guide = guide_legend(nrow = 2, 
                                                override.aes = list(linetype=c(1,1,1,0), 
                                                                    shape=16))) +
        theme(plot.title = element_text(size=22, face="bold"),
              plot.subtitle = element_text(size=22, face="bold"),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              axis.text.x = element_text(color = "black", size = 22),
              axis.text.y = element_text(color = "black", size = 22),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1.2),
              panel.background = element_blank(), legend.title=element_blank(),
              legend.key.width = unit(1.4, "cm"),
              legend.spacing.x = unit(1, 'cm'), legend.text=element_text(size=22),
              legend.position = 'bottom') +
        labs(title = testname, subtitle = plot.subtitle,
             y = 'Likelihood ratio',
             x = 'Steps in sequential analyses')
      
      suppressWarnings(print(seqcompare))
    }
    
    return(list('n' = n0, 'decision' = decision, 
                'Accept.threshold' = Accept.threshold, 'Reject.threshold' = Reject.threshold,
                'LR' = LR[1:nAnalyses.test], 'theta.UMPBT' = theta.UMPBT))
    
    # end one-sided oneZ
    
  }else{
    
    #################### one-sample z (both sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch.size = design.MSPRT.object$batch.size
      N.max = design.MSPRT.object$N.max
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      sigma = design.MSPRT.object$sigma
      termination.threshold = design.MSPRT.object$termination.threshold
      theta.UMPBT = design.MSPRT.object$theta.UMPBT
      nAnalyses.max = design.MSPRT.object$nAnalyses
      
      # msg
      if(verbose){
        
        if(any(batch.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the group sequential MSPRT for a one-sample z test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the sequential MSPRT for a one-sample z test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Known standard deviation: ", sigma, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(theta.UMPBT$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta.UMPBT$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
      batch.size = c(0, cumsum(batch.size))
      nAnalyses = max(which(batch.size<=length(obs))) - 1
      
    }else{
      
      ## ignoring obs1 & obs2
      if(!missing(obs1)) print("'obs1' is ignored. Not required in one-sample tests.")
      if(!missing(obs2)) print("'obs2' is ignored. Not required in one-sample tests.")
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
      if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
      
      ## ignoring N1.max & N2.max
      if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
      if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
      
      ## batch sizes and N.max
      if(missing(batch.size)){
        
        if(missing(N.max)){
          
          return("Either 'batch.size' or 'N.max' needs to be specified")
          
        }else{batch.size = rep(1, N.max)}
        
      }else{
        
        if(missing(N.max)){
          
          N.max = sum(batch.size)
          
        }else{
          
          if(sum(batch.size)!=N.max) return("Sum of batch sizes should add up to N.max")
        }
      }
      
      nAnalyses.max = length(batch.size)
      
      ## point H0
      if(missing(theta0)) theta0 = 0
      
      ######################## UMPBT alternative ########################
      theta.UMPBT = list('right' = UMPBT.alt(test.type = 'oneZ', side = 'right', 
                                             theta0 = theta0, N = N.max, 
                                             Type1 = Type1.target/2, sigma = sigma),
                         'left' = UMPBT.alt(test.type = 'oneZ', side = 'left', 
                                            theta0 = theta0, N = N.max,
                                            Type1 = Type1.target/2, sigma = sigma))
      
      # msg
      if(verbose){
        
        if(any(batch.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the group sequential MSPRT for a one-sample z test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the sequential MSPRT for a one-sample z test:")
          print("==========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Known standard deviation: ", sigma, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(theta.UMPBT$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta.UMPBT$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
      batch.size = c(0, cumsum(batch.size))
      nAnalyses = max(which(batch.size<=length(obs))) - 1
    }
    
    
    ###################### sequential comparison ######################
    # Wald's thresholds
    Accept.threshold = Type2.target/(1 - Type1.target/2)
    Reject.threshold = (1 - Type2.target)/(Type1.target/2)
    
    # required storages
    cumsum_n = 0
    reached.decision.r = reached.decision.l = reached.decision = F
    rejectH0 = NA
    LR.r = LR.l = rep(NA, nAnalyses)
    
    for(n in 1:nAnalyses){
      
      if(!reached.decision){
        
        # sum of observations until step n
        cumsum_n = cumsum_n + sum(obs[(batch.size[n]+1):batch.size[n+1]])
        
        ## for right sided check
        if(!reached.decision.r){
          
          # likelihood ratio of observations until step n
          LR.r[n] = 
            exp((cumsum_n*(theta.UMPBT$right - theta0) -
                   ((batch.size[n+1]*((theta.UMPBT$right^2) - (theta0^2)))/2))/(sigma^2))
          
          # comparing with the thresholds
          AcceptedH0_n.r = LR.r[n]<=Accept.threshold
          RejectedH0_n.r = LR.r[n]>=Reject.threshold
          reached.decision.r = AcceptedH0_n.r||RejectedH0_n.r
        }
        
        ## for left sided check
        if(!reached.decision.l){
          
          # likelihood ratio of observations until step n
          LR.l[n] = 
            exp((cumsum_n*(theta.UMPBT$left - theta0) -
                   ((batch.size[n+1]*((theta.UMPBT$left^2) - (theta0^2)))/2))/(sigma^2))
          
          # comparing with the thresholds
          AcceptedH0_n.l = LR.l[n]<=Accept.threshold
          RejectedH0_n.l = LR.l[n]>=Reject.threshold
          reached.decision.l = AcceptedH0_n.l||RejectedH0_n.l
        }
        
        ## both-sided check
        if(AcceptedH0_n.r&&AcceptedH0_n.l){
          
          rejectH0 = F
          decision = 'accept'
          n0 = batch.size[n+1]
          reached.decision = T
          
        }else if(RejectedH0_n.r||RejectedH0_n.l){
          
          rejectH0 = T
          decision = 'reject'
          n0 = batch.size[n+1]
          reached.decision = T
        }
      }
    }
    
    # inconclusive cases
    if(!reached.decision){
      
      if(nAnalyses==nAnalyses.max){
        
        n0 = N.max
        if(AcceptedH0_n.l&&(!reached.decision.r)){
          
          rejectH0 = LR.r[nAnalyses]>=termination.threshold
          
        }else if(AcceptedH0_n.r&&(!reached.decision.l)){
          
          rejectH0 = LR.l[nAnalyses]>=termination.threshold
          
        }else if((!reached.decision.r)&&(!reached.decision.l)){
          
          rejectH0 = max(LR.r[nAnalyses], LR.l[nAnalyses])>=termination.threshold
          
        }else{rejectH0 = F}
        
        if(rejectH0){
          
          decision = 'reject'
          
        }else if(!rejectH0){decision = 'accept'}
        reached.decision = T
        
      }else{
        
        n0 = batch.size[nAnalyses + 1]
        decision = 'continue'
      }
    }
    
    # msg
    if(verbose==T){
      
      if(decision=='continue'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Continue sampling')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='reject'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Reject the null hypothesis')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='accept'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Accept the null hypothesis')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
      }
    }
    
    
    nAnalyses.r = max(which(!is.na(LR.r)))
    nAnalyses.l = max(which(!is.na(LR.l)))
    
    ## plotting
    if(plot.it==T){
      
      # title and decision
      testname = 'Two-sided one-sample z test'
      if(decision=="accept"){
        
        plot.subtitle = 
          paste('Accept the null (n = ', n0, ')', sep = '')
        
      }else if(decision=="reject"){
        
        plot.subtitle = 
          paste('Reject the null (n = ', n0, ')', sep = '')
        
      }else if(decision=="continue"){
        
        plot.subtitle = 
          paste('Continue sampling (n = ', n0, ')', sep = '')
      }
      
      # right sided test
      if(AcceptedH0_n.r){
        
        min.LR.r = min(LR.r, na.rm = T)
        if(min.LR.r<Accept.threshold){
          
          ylow.r = min.LR.r
          yup.r = max(LR.r, na.rm = T)
          
        }else if((min.LR.r>=Accept.threshold)&&
                 (min.LR.r<(termination.threshold + Accept.threshold)/2)){
          
          ylow.r = Accept.threshold
          yup.r = termination.threshold
          
        }else{
          
          ylow.r = Accept.threshold
          yup.r = Reject.threshold
        }
        
      }else if(RejectedH0_n.r){
        
        ylow.r = Accept.threshold
        yup.r = max(LR.r, na.rm = T)
        
      }else if(!reached.decision.r){
        
        last.LR.r = LR.r[max(which(!is.na(LR.r)))]
        if(last.LR.r<(termination.threshold + Accept.threshold)/2){
          
          ylow.r = Accept.threshold
          yup.r = termination.threshold
          
        }else{
          
          ylow.r = Accept.threshold
          yup.r = Reject.threshold
        }
      }
      
      df.r = rbind.data.frame(data.frame('xval' = 1:nAnalyses.r,
                                         'yval' = Accept.threshold,
                                         'type' = 'A'),
                              data.frame('xval' = 1:nAnalyses.r,
                                         'yval' = Reject.threshold,
                                         'type' = 'R'),
                              data.frame('xval' = 1:nAnalyses.r,
                                         'yval' = LR.r[1:nAnalyses.r],
                                         'type' = 'LR'),
                              data.frame('xval' = design.MSPRT.object$nAnalyses,
                                         'yval' = termination.threshold,
                                         'type' = 'term.thresh'))
      
      df.r$type = factor(as.character(df.r$type),
                         levels = c('A', 'R', 'LR', 'term.thresh'))
      
      seqcompare.r = ggplot(data = df.r,
                            aes(x = xval, y = yval, group = type)) + 
        geom_point(aes(colour = type), size = 2) +
        geom_line(aes(colour = type), size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Accept.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="forestgreen", size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Reject.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="red2", size=1) +
        geom_point(aes(x = length(batch.size) - 1, y = termination.threshold),
                   colour = "black", size = 2) +
        xlim(c(1, nAnalyses.r)) +
        ylim(c(ylow.r, yup.r)) +
        scale_color_manual(labels = c('Acceptance threshold', 'Rejection threshold',
                                      'Likelihood ratio', 'Termination threshold'),
                           values = c('A' = 'forestgreen', 'R' = 'red2',
                                      'LR' = 'dodgerblue', 'term.thresh' = 'black'),
                           guide = guide_legend(nrow = 2,
                                                override.aes = list(linetype=c(1,1,1,0), 
                                                                    shape=16))) +
        theme(plot.title = element_text(size=22, face="bold"),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              axis.text.x = element_text(color = "black", size = 22),
              axis.text.y = element_text(color = "black", size = 22),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1.2),
              panel.background = element_blank(), legend.title=element_blank(),
              legend.key.width = unit(1.4, "cm"),
              legend.spacing.x = unit(1, 'cm'), legend.text=element_text(size=22),
              legend.position = 'bottom') +
        labs(title = paste('Right-sided test at level ', Type1.target/2, sep = ''),
             y = 'Likelihood ratio',
             x = 'Steps in sequential analyses')
      
      # left sided test
      if(AcceptedH0_n.l){
        
        min.LR.l = min(LR.l, na.rm = T)
        if(min.LR.l<Accept.threshold){
          
          ylow.l = min.LR.l
          yup.l = max(LR.l, na.rm = T)
          
        }else if((min.LR.l>=Accept.threshold)&&
                 (min.LR.l<(termination.threshold + Accept.threshold)/2)){
          
          ylow.l = Accept.threshold
          yup.l = termination.threshold
          
        }else{
          
          ylow.l = Accept.threshold
          yup.l = Reject.threshold
        }
        
      }else if(RejectedH0_n.l){
        
        ylow.l = Accept.threshold
        yup.l = max(LR.l, na.rm = T)
        
      }else if(!reached.decision.l){
        
        last.LR.l = LR.l[max(which(!is.na(LR.l)))]
        if(last.LR.l<(termination.threshold + Accept.threshold)/2){
          
          ylow.l = Accept.threshold
          yup.l = termination.threshold
          
        }else{
          
          ylow.l = Accept.threshold
          yup.l = Reject.threshold
        }
      }
      
      df.l = rbind.data.frame(data.frame('xval' = 1:nAnalyses.l,
                                         'yval' = Accept.threshold,
                                         'type' = 'A'),
                              data.frame('xval' = 1:nAnalyses.l,
                                         'yval' = Reject.threshold,
                                         'type' = 'R'),
                              data.frame('xval' = 1:nAnalyses.l,
                                         'yval' = LR.l[1:nAnalyses.l],
                                         'type' = 'LR'),
                              data.frame('xval' = design.MSPRT.object$nAnalyses,
                                         'yval' = termination.threshold,
                                         'type' = 'term.thresh'))
      
      df.l$type = factor(as.character(df.l$type),
                         levels = c('A', 'R', 'LR', 'term.thresh'))
      
      seqcompare.l = ggplot(data = df.l,
                            aes(x = xval, y = yval, group = type)) + 
        geom_point(aes(colour = type), size = 2) +
        geom_line(aes(colour = type), size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Accept.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="forestgreen", size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Reject.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="red2", size=1) +
        geom_point(aes(x = length(batch.size) - 1, y = termination.threshold),
                   colour = "black", size = 2) +
        xlim(c(1, nAnalyses.l)) +
        ylim(c(ylow.l, yup.l)) +
        scale_color_manual(labels = c('Acceptance threshold', 'Rejection threshold',
                                      'Likelihood ratio', 'Termination threshold'),
                           values = c('A' = 'forestgreen', 'R' = 'red2',
                                      'LR' = 'dodgerblue', 'term.thresh' = 'black'),
                           guide = guide_legend(nrow = 2)) +
        theme(plot.title = element_text(size=22, face="bold"),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              axis.text.x = element_text(color = "black", size = 22),
              axis.text.y = element_text(color = "black", size = 22),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1.2),
              panel.background = element_blank(), legend.title=element_blank(),
              legend.key.width = unit(1.4, "cm"),
              legend.spacing.x = unit(1, 'cm'), legend.text=element_text(size=22),
              legend.position = 'bottom') +
        labs(title = paste('Left-sided test at level ', Type1.target/2, sep = ''),
             y = 'Likelihood ratio',
             x = 'Steps in sequential analyses')
      
      suppressWarnings(print(annotate_figure(ggarrange(seqcompare.l, seqcompare.r, 
                                                       nrow = 1, ncol = 2,
                                                       legend = 'bottom', common.legend = T),
                                             top = text_grob(paste(testname, ': ', plot.subtitle, '\n'), face = "bold", 
                                                             size = 25, hjust = .5))))
    }
    
    return(list('n' = n0, 'decision' = decision, 
                'Accept.threshold' = Accept.threshold, 'Reject.threshold' = Reject.threshold,
                'LR' = list('right' = LR.r[1:nAnalyses.r], 'left' = LR.l[1:nAnalyses.l]),
                'theta.UMPBT' = theta.UMPBT))
    
  } # end both-sided oneZ
}

#### one-sample t test ####
implement.MSPRT_oneT = function(obs, design.MSPRT.object, 
                                termination.threshold,
                                side = 'right', theta0 = 0, 
                                Type1.target =.005, Type2.target = .2,
                                N.max, batch.size,
                                verbose = T, plot.it = T){
  
  # side
  if(!missing(design.MSPRT.object)) side = design.MSPRT.object$side
  
  if(side!='both'){
    
    #################### one-sample t (right/left sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch.size = design.MSPRT.object$batch.size
      N.max = design.MSPRT.object$N.max
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      termination.threshold = design.MSPRT.object$termination.threshold
      nAnalyses.max = design.MSPRT.object$nAnalyses
      
      # msg
      if(verbose){
        
        if((batch.size[1]>2)||any(batch.size[-1]>1)){
          
          cat('\n')
          print("=========================================================================")
          print("Implementing the group sequential MSPRT for a one-sample t test:")
          print("=========================================================================")
          
        }else{
          
          cat('\n')
          print("=========================================================================")
          print("Implementing the sequential MSPRT for a one-sample t test:")
          print("=========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Termination threshold: ", design.MSPRT.object$termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
      batch.size = c(0, cumsum(batch.size))
      nAnalyses = max(which(batch.size<=length(obs))) - 1
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(obs1)) print("'obs1' is ignored. Not required in one-sample tests.")
      if(!missing(obs2)) print("'obs2' is ignored. Not required in one-sample tests.")
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
      if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
      
      ## ignoring N1.max & N2.max
      if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
      if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
      
      ## batch sizes and N.max
      if(missing(batch.size)){
        
        if(missing(N.max)){
          
          return("Either 'batch.size' or 'N.max' needs to be specified")
          
        }else{batch.size = c(2, rep(1, N.max-2))}
        
      }else{
        
        if(batch.size[1]<2){
          
          return("First batch size should be at least 2")
          
        }else{
          
          if(missing(N.max)){
            
            N.max = sum(batch.size)
            
          }else{
            
            if(sum(batch.size)!=N.max) return("Sum of batch.size should add up to N.max")
          }
        }
      }
      
      nAnalyses.max = length(batch.size)
      
      ## point H0
      if(missing(theta0)) theta0 = 0
      
      # msg
      if(verbose){
        
        if((batch.size[1]>2)||any(batch.size[-1]>1)){
          
          cat('\n')
          print("=========================================================================")
          print("Implementing the group sequential MSPRT for a one-sample t test:")
          print("=========================================================================")
          
        }else{
          
          cat('\n')
          print("=========================================================================")
          print("Implementing the sequential MSPRT for a one-sample t test:")
          print("=========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Termination threshold: ", termination.threshold,
                    sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
      batch.size = c(0, cumsum(batch.size))
      nAnalyses = max(which(batch.size<=length(obs))) - 1
    }
    
    
    ###################### sequential comparison ######################
    # Wald's thresholds
    Accept.threshold = Type2.target/(1 - Type1.target)
    Reject.threshold = (1 - Type2.target)/Type1.target
    
    # cut-off (with sign) in fixed design one-sample t test
    signed_t.alpha = (2*(side=='right')-1)*qt(Type1.target, df = N.max -1, 
                                              lower.tail = F)
    
    # required storages
    cumSS_n = cumsum_n = 0
    reached.decision = F
    rejectH0 = NA
    theta.UMPBT = LR = rep(NA, nAnalyses)
    
    for(n in 1:nAnalyses){
      
      if(!reached.decision){
        
        # sum of observations until step n
        cumsum_n = cumsum_n + sum(obs[(batch.size[n]+1):batch.size[n+1]])
        
        # sum of squares of observations until step n
        cumSS_n = cumSS_n + sum(obs[(batch.size[n]+1):batch.size[n+1]]^2)
        
        # xbar and (n-1)*(s^2) until step n
        xbar_n = cumsum_n/batch.size[n+1]
        divisor.s_n.sq = cumSS_n - (cumsum_n^2)/batch.size[n+1]
        
        # UMPBT alternative
        theta.UMPBT[n] = theta0 + signed_t.alpha*
          sqrt(divisor.s_n.sq/(N.max*(batch.size[n+1]-1)))
        
        # likelihood ratio of observations until step n
        LR[n] = 
          ((1 + (batch.size[n+1]*((xbar_n - theta0)^2))/divisor.s_n.sq)/
             (1 + (batch.size[n+1]*((xbar_n - theta.UMPBT[n])^2))/
                divisor.s_n.sq))^(batch.size[n+1]/2)
        
        # comparing with the thresholds
        AcceptedH0_n = LR[n]<=Accept.threshold
        RejectedH0_n = LR[n]>=Reject.threshold
        reached.decision = AcceptedH0_n||RejectedH0_n
        if(reached.decision){
          
          n0 = batch.size[n+1]
          rejectH0 = RejectedH0_n
          if(rejectH0){
            
            decision = 'reject'
            
          }else{decision = 'accept'}
        }
      }
    }
    
    # inconclusive cases
    if(!reached.decision){
      
      if(nAnalyses==nAnalyses.max){
        
        n0 = N.max
        rejectH0 = LR[nAnalyses]>=termination.threshold
        
        if(rejectH0){
          
          decision = 'reject'
          
        }else if(!rejectH0){decision = 'accept'}
        reached.decision = T
        
      }else{
        
        n0 = batch.size[nAnalyses+1]
        decision = 'continue'
      }
    }
    
    # msg
    if(verbose==T){
      
      if(decision=='continue'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Continue sampling')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='reject'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Reject the null hypothesis')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='accept'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Accept the null hypothesis')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
      }
    }
    
    
    nAnalyses.test = max(which(!is.na(LR)))
    
    ## plotting
    if(plot.it==T){
      
      # title and decision
      if(side=='right'){
        
        testname = 'Right-sided one-sample t test'
        
      }else{testname = 'Left-sided one-sample t test'}
      
      if(decision=="accept"){
        
        plot.subtitle = 
          paste('Accept the null (n = ', n0, ')', sep = '')
        
        min.LR = min(LR, na.rm = T)
        if(min.LR<Accept.threshold){
          
          ylow = min.LR
          yup = max(LR, na.rm = T)
          
        }else if((min.LR>=Accept.threshold)&&
                 (min.LR<(termination.threshold + Accept.threshold)/2)){
          
          ylow = Accept.threshold
          yup = termination.threshold
          
        }else{
          
          ylow = Accept.threshold
          yup = Reject.threshold
        }
        
      }else if(decision=="reject"){
        
        plot.subtitle = 
          paste('Reject the null (n = ', n0, ')', sep = '')
        
        ylow = Accept.threshold
        yup = max(LR, na.rm = T)
        
      }else if(decision=="continue"){
        
        plot.subtitle = 
          paste('Continue sampling (n = ', n0, ')', sep = '')
        
        last.LR = LR[max(which(!is.na(LR)))]
        if(last.LR<(termination.threshold + Accept.threshold)/2){
          
          ylow = Accept.threshold
          yup = termination.threshold
          
        }else{
          
          ylow = Accept.threshold
          yup = Reject.threshold
        }
      }
      
      df = rbind.data.frame(data.frame('xval' = 1:nAnalyses.test,
                                       'yval' = Accept.threshold,
                                       'type' = 'A'),
                            data.frame('xval' = 1:nAnalyses.test,
                                       'yval' = Reject.threshold,
                                       'type' = 'R'),
                            data.frame('xval' = 1:nAnalyses.test,
                                       'yval' = LR[1:nAnalyses.test],
                                       'type' = 'LR'),
                            data.frame('xval' = design.MSPRT.object$nAnalyses,
                                       'yval' = termination.threshold,
                                       'type' = 'term.thresh'))
      
      df$type = factor(as.character(df$type),
                       levels = c('A', 'R', 'LR', 'term.thresh'))
      
      seqcompare = ggplot(data = df,
                          aes(x = xval, y = yval, group = type)) + 
        geom_point(aes(colour = type), size = 2) +
        geom_line(aes(colour = type), size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Accept.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="forestgreen", size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Reject.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="red2", size=1) +
        geom_point(aes(x = length(batch.size) - 1, y = termination.threshold),
                   colour = "black", size = 2) +
        xlim(c(1, nAnalyses.test)) +
        ylim(c(ylow, yup)) +
        scale_color_manual(labels = c('Acceptance threshold', 'Rejection threshold',
                                      'Bayes factor', 'Termination threshold'),
                           values = c('A' = 'forestgreen', 'R' = 'red2',
                                      'LR' = 'dodgerblue', 'term.thresh' = 'black'),
                           guide = guide_legend(nrow = 2, 
                                                override.aes = list(linetype=c(1,1,1,0), 
                                                                    shape=16))) +
        theme(plot.title = element_text(size=22, face="bold"),
              plot.subtitle = element_text(size=18, face="bold"),
              axis.title.x = element_text(size=15),
              axis.title.y = element_text(size=15),
              axis.text.x = element_text(color = "black", size = 15),
              axis.text.y = element_text(color = "black", size = 15),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1.2),
              panel.background = element_blank(), legend.title=element_blank(),
              legend.key.width = unit(1.4, "cm"),
              legend.spacing.x = unit(1, 'cm'), legend.text=element_text(size=20),
              legend.position = 'bottom') +
        labs(title = testname, subtitle = plot.subtitle,
             y = 'Bayes factor',
             x = 'Steps in sequential analyses')
      
      suppressWarnings(print(seqcompare))
    }
    
    return(list('n' = n0, 'decision' = decision, 
                'Accept.threshold' = Accept.threshold, 'Reject.threshold' = Reject.threshold,
                'LR' = LR[1:nAnalyses.test], 'theta.UMPBT' = theta.UMPBT[1:nAnalyses.test]))
    
    # end one-sided oneT
    
  }else{
    
    #################### one-sample t (both sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch.size = design.MSPRT.object$batch.size
      N.max = design.MSPRT.object$N.max
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      termination.threshold = design.MSPRT.object$termination.threshold
      nAnalyses.max = design.MSPRT.object$nAnalyses
      
      # msg
      if(verbose){
        
        if((batch.size[1]>2)||any(batch.size[-1]>1)){
          
          cat('\n')
          print("=========================================================================")
          print("Implementing the group sequential MSPRT for a one-sample t test:")
          print("=========================================================================")
          
        }else{
          
          cat('\n')
          print("=========================================================================")
          print("Implementing the sequential MSPRT for a one-sample t test:")
          print("=========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
      batch.size = c(0, cumsum(batch.size))
      nAnalyses = max(which(batch.size<=length(obs))) - 1
      
    }else{
      
      ## ignoring obs1 & obs2
      if(!missing(obs1)) print("'obs1' is ignored. Not required in one-sample tests.")
      if(!missing(obs2)) print("'obs2' is ignored. Not required in one-sample tests.")
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
      if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
      
      ## ignoring N1.max & N2.max
      if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
      if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
      
      ## batch sizes and N.max
      if(missing(batch.size)){
        
        if(missing(N.max)){
          
          return("Either 'batch.size' or 'N.max' needs to be specified")
          
        }else{batch.size = c(2, rep(1, N.max-2))}
        
      }else{
        
        if(batch.size[1]<2){
          
          return("First batch size should be at least 2")
          
        }else{
          
          if(missing(N.max)){
            
            N.max = sum(batch.size)
            
          }else{
            
            if(sum(batch.size)!=N.max) return("Sum of batch.size should add up to N.max")
          }
        }
      }
      
      nAnalyses.max = length(batch.size)
      
      ## point H0
      if(missing(theta0)) theta0 = 0
      
      # msg
      if(verbose){
        
        if((batch.size[1]>2)||any(batch.size[-1]>1)){
          
          cat('\n')
          print("=========================================================================")
          print("Implementing the group sequential MSPRT for a one-sample t test:")
          print("=========================================================================")
          
        }else{
          
          cat('\n')
          print("=========================================================================")
          print("Implementing the sequential MSPRT for a one-sample t test:")
          print("=========================================================================")
        }
        
        print(paste("Maximum available sample size: ", N.max, sep = ""))
        print(paste('Batch sizes: ', paste(batch.size, collapse = ', '), sep = ''))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max, sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
      batch.size = c(0, cumsum(batch.size))
      nAnalyses = max(which(batch.size<=length(obs))) - 1
    }
    
    
    ###################### sequential comparison ######################
    # Wald's thresholds
    Accept.threshold = Type2.target/(1 - Type1.target/2)
    Reject.threshold = (1 - Type2.target)/(Type1.target/2)
    
    # cut-off (with sign) in fixed design one-sample t test
    t.alpha = qt(Type1.target/2, df = N.max -1, lower.tail = F)
    
    # required storages
    cumSS_n = cumsum_n = 0
    reached.decision.r = reached.decision.l = reached.decision = F
    rejectH0 = NA
    theta.UMPBT.r = theta.UMPBT.l = LR.r = LR.l = rep(NA, nAnalyses)
    
    
    for(n in 1:nAnalyses){
      
      if(!reached.decision){
        
        # sum of observations until step n
        cumsum_n = cumsum_n + sum(obs[(batch.size[n]+1):batch.size[n+1]])
        
        # sum of squares of observations until step n
        cumSS_n = cumSS_n + sum(obs[(batch.size[n]+1):batch.size[n+1]]^2)
        
        ## for right sided check
        if(!reached.decision.r){
          
          # xbar and (n-1)*(s^2) until step n
          xbar_n.r = cumsum_n/batch.size[n+1]
          divisor.s_n.sq.r = cumSS_n - ((cumsum_n)^2)/batch.size[n+1]
          
          # UMPBT alternative
          theta.UMPBT.r[n] = theta0 + t.alpha*
            sqrt(divisor.s_n.sq.r/(N.max*(batch.size[n+1]-1)))
          
          # likelihood ratio of observations until step n
          LR.r[n] = 
            ((1 + (batch.size[n+1]*((xbar_n.r - theta0)^2))/divisor.s_n.sq.r)/
               (1 + (batch.size[n+1]*((xbar_n.r - theta.UMPBT.r[n])^2))/
                  divisor.s_n.sq.r))^(batch.size[n+1]/2)
          
          # comparing with the thresholds
          AcceptedH0_n.r = LR.r[n]<=Accept.threshold
          RejectedH0_n.r = LR.r[n]>=Reject.threshold
          reached.decision.r = AcceptedH0_n.r||RejectedH0_n.r
        }
        
        ## for left sided check
        if(!reached.decision.l){
          
          # xbar and (n-1)*(s^2) until step n
          xbar_n.l = cumsum_n/batch.size[n+1]
          divisor.s_n.sq.l = cumSS_n - ((cumsum_n)^2)/batch.size[n+1]
          
          # UMPBT alternative
          theta.UMPBT.l[n] = theta0 - t.alpha*
            sqrt(divisor.s_n.sq.l/(N.max*(batch.size[n+1]-1)))
          
          # likelihood ratio of observations until step n
          LR.l[n] = 
            ((1 + (batch.size[n+1]*((xbar_n.l - theta0)^2))/divisor.s_n.sq.l)/
               (1 + (batch.size[n+1]*((xbar_n.l - theta.UMPBT.l[n])^2))/
                  divisor.s_n.sq.l))^(batch.size[n+1]/2)
          
          # comparing with the thresholds
          AcceptedH0_n.l = LR.l[n]<=Accept.threshold
          RejectedH0_n.l = LR.l[n]>=Reject.threshold
          reached.decision.l = AcceptedH0_n.l||RejectedH0_n.l
        }
        
        ## both-sided check
        if(AcceptedH0_n.r&&AcceptedH0_n.l){
          
          rejectH0 = F
          decision = 'accept'
          n0 = batch.size[n+1]
          reached.decision = T
          
        }else if(RejectedH0_n.r||RejectedH0_n.l){
          
          rejectH0 = T
          decision = 'reject'
          n0 = batch.size[n+1]
          reached.decision = T
        }
      }
    }
    
    # inconclusive cases
    if(!reached.decision){
      
      if(nAnalyses==nAnalyses.max){
        
        n0 = N.max
        if(AcceptedH0_n.l&&(!reached.decision.r)){
          
          rejectH0 = LR.r[nAnalyses]>=termination.threshold
          
        }else if(AcceptedH0_n.r&&(!reached.decision.l)){
          
          rejectH0 = LR.l[nAnalyses]>=termination.threshold
          
        }else if((!reached.decision.r)&&(!reached.decision.l)){
          
          rejectH0 = max(LR.r[nAnalyses], LR.l[nAnalyses])>=termination.threshold
          
        }else{rejectH0 = F}
        
        if(rejectH0){
          
          decision = 'reject'
          
        }else if(!rejectH0){decision = 'accept'}
        reached.decision = T
        
      }else{
        
        n0 = batch.size[nAnalyses + 1]
        decision = 'continue'
      }
    }
    
    # msg
    if(verbose==T){
      
      if(decision=='continue'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Continue sampling')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='reject'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Reject the null hypothesis')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='accept'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Accept the null hypothesis')
        print(paste("Sample size used: ", n0, sep = ''))
        print("=========================================================================")
        cat('\n')
      }
    }
    
    
    nAnalyses.r = max(which(!is.na(LR.r)))
    nAnalyses.l = max(which(!is.na(LR.l)))
    
    ## plotting
    if(plot.it==T){
      
      # title and decision
      testname = 'Two-sided one-sample t test'
      if(decision=="accept"){
        
        plot.subtitle = 
          paste('Accept the null (n = ', n0, ')', sep = '')
        
      }else if(decision=="reject"){
        
        plot.subtitle = 
          paste('Reject the null (n = ', n0, ')', sep = '')
        
      }else if(decision=="continue"){
        
        plot.subtitle = 
          paste('Continue sampling (n = ', n0, ')', sep = '')
      }
      
      # right sided test
      if(AcceptedH0_n.r){
        
        min.LR.r = min(LR.r, na.rm = T)
        if(min.LR.r<Accept.threshold){
          
          ylow.r = min.LR.r
          yup.r = max(LR.r, na.rm = T)
          
        }else if((min.LR.r>=Accept.threshold)&&
                 (min.LR.r<(termination.threshold + Accept.threshold)/2)){
          
          ylow.r = Accept.threshold
          yup.r = termination.threshold
          
        }else{
          
          ylow.r = Accept.threshold
          yup.r = Reject.threshold
        }
        
      }else if(RejectedH0_n.r){
        
        ylow.r = Accept.threshold
        yup.r = max(LR.r, na.rm = T)
        
      }else if(!reached.decision.r){
        
        last.LR.r = LR.r[max(which(!is.na(LR.r)))]
        if(last.LR.r<(termination.threshold + Accept.threshold)/2){
          
          ylow.r = Accept.threshold
          yup.r = termination.threshold
          
        }else{
          
          ylow.r = Accept.threshold
          yup.r = Reject.threshold
        }
      }
      
      df.r = rbind.data.frame(data.frame('xval' = 1:nAnalyses.r,
                                         'yval' = Accept.threshold,
                                         'type' = 'A'),
                              data.frame('xval' = 1:nAnalyses.r,
                                         'yval' = Reject.threshold,
                                         'type' = 'R'),
                              data.frame('xval' = 1:nAnalyses.r,
                                         'yval' = LR.r[1:nAnalyses.r],
                                         'type' = 'LR'),
                              data.frame('xval' = design.MSPRT.object$nAnalyses,
                                         'yval' = termination.threshold,
                                         'type' = 'term.thresh'))
      
      df.r$type = factor(as.character(df.r$type),
                         levels = c('A', 'R', 'LR', 'term.thresh'))
      
      seqcompare.r = ggplot(data = df.r,
                            aes(x = xval, y = yval, group = type)) + 
        geom_point(aes(colour = type), size = 2) +
        geom_line(aes(colour = type), size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Accept.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="forestgreen", size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Reject.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="red2", size=1) +
        geom_point(aes(x = length(batch.size) - 1, y = termination.threshold),
                   colour = "black", size = 2) +
        xlim(c(1, nAnalyses.r)) +
        ylim(c(ylow.r, yup.r)) +
        scale_color_manual(labels = c('Acceptance threshold', 'Rejection threshold',
                                      'Bayes factor', 'Termination threshold'),
                           values = c('A' = 'forestgreen', 'R' = 'red2',
                                      'LR' = 'dodgerblue', 'term.thresh' = 'black'),
                           guide = guide_legend(nrow = 2,
                                                override.aes = list(linetype=c(1,1,1,0), 
                                                                    shape=16))) +
        theme(plot.title = element_text(size=22, face="bold"),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              axis.text.x = element_text(color = "black", size = 22),
              axis.text.y = element_text(color = "black", size = 22),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1.2),
              panel.background = element_blank(), legend.title=element_blank(),
              legend.key.width = unit(1.4, "cm"),
              legend.spacing.x = unit(1, 'cm'), legend.text=element_text(size=22),
              legend.position = 'bottom') +
        labs(title = paste('Right-sided test at level ', Type1.target/2, sep = ''),
             y = 'Bayes factor',
             x = 'Steps in sequential analyses')
      
      # left sided test
      if(AcceptedH0_n.l){
        
        min.LR.l = min(LR.l, na.rm = T)
        if(min.LR.l<Accept.threshold){
          
          ylow.l = min.LR.l
          yup.l = max(LR.l, na.rm = T)
          
        }else if((min.LR.l>=Accept.threshold)&&
                 (min.LR.l<(termination.threshold + Accept.threshold)/2)){
          
          ylow.l = Accept.threshold
          yup.l = termination.threshold
          
        }else{
          
          ylow.l = Accept.threshold
          yup.l = Reject.threshold
        }
        
      }else if(RejectedH0_n.l){
        
        ylow.l = Accept.threshold
        yup.l = max(LR.l, na.rm = T)
        
      }else if(!reached.decision.l){
        
        last.LR.l = LR.l[max(which(!is.na(LR.l)))]
        if(last.LR.l<(termination.threshold + Accept.threshold)/2){
          
          ylow.l = Accept.threshold
          yup.l = termination.threshold
          
        }else{
          
          ylow.l = Accept.threshold
          yup.l = Reject.threshold
        }
      }
      
      df.l = rbind.data.frame(data.frame('xval' = 1:nAnalyses.l,
                                         'yval' = Accept.threshold,
                                         'type' = 'A'),
                              data.frame('xval' = 1:nAnalyses.l,
                                         'yval' = Reject.threshold,
                                         'type' = 'R'),
                              data.frame('xval' = 1:nAnalyses.l,
                                         'yval' = LR.l[1:nAnalyses.l],
                                         'type' = 'LR'),
                              data.frame('xval' = design.MSPRT.object$nAnalyses,
                                         'yval' = termination.threshold,
                                         'type' = 'term.thresh'))
      
      df.l$type = factor(as.character(df.l$type),
                         levels = c('A', 'R', 'LR', 'term.thresh'))
      
      seqcompare.l = ggplot(data = df.l,
                            aes(x = xval, y = yval, group = type)) + 
        geom_point(aes(colour = type), size = 2) +
        geom_line(aes(colour = type), size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Accept.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="forestgreen", size = 1) +
        geom_segment(aes(x = length(batch.size) - 1, y = Reject.threshold,
                         xend = length(batch.size) - 1, yend = termination.threshold),
                     color="red2", size=1) +
        geom_point(aes(x = length(batch.size) - 1, y = termination.threshold),
                   colour = "black", size = 2) +
        xlim(c(1, nAnalyses.l)) +
        ylim(c(ylow.l, yup.l)) +
        scale_color_manual(labels = c('Acceptance threshold', 'Rejection threshold',
                                      'Bayes factor', 'Termination threshold'),
                           values = c('A' = 'forestgreen', 'R' = 'red2',
                                      'LR' = 'dodgerblue', 'term.thresh' = 'black'),
                           guide = guide_legend(nrow = 2)) +
        theme(plot.title = element_text(size=22, face="bold"),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              axis.text.x = element_text(color = "black", size = 22),
              axis.text.y = element_text(color = "black", size = 22),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1.2),
              panel.background = element_blank(), legend.title=element_blank(),
              legend.key.width = unit(1.4, "cm"),
              legend.spacing.x = unit(1, 'cm'), legend.text=element_text(size=22),
              legend.position = 'bottom') +
        labs(title = paste('Left-sided test at level ', Type1.target/2, sep = ''),
             y = 'Bayes factor',
             x = 'Steps in sequential analyses')
      
      suppressWarnings(print(annotate_figure(ggarrange(seqcompare.l, seqcompare.r, 
                                                       nrow = 1, ncol = 2,
                                                       legend = 'bottom', common.legend = T),
                                             top = text_grob(paste(testname, ': ', plot.subtitle, '\n'), face = "bold", 
                                                             size = 25, hjust = .5))))
    }
    
    return(list('n' = n0, 'decision' = decision, 
                'Accept.threshold' = Accept.threshold, 'Reject.threshold' = Reject.threshold,
                'LR' = list('right' = LR.r[1:nAnalyses.r], 'left' = LR.l[1:nAnalyses.l]),
                'theta.UMPBT' = list('right' = theta.UMPBT.r[1:nAnalyses.r],
                                     'left' = theta.UMPBT.l[1:nAnalyses.l])))
    
  } # end both-sided oneT
}


#### two-sample z test ####
implement.MSPRT_twoZ = function(obs1, obs2, design.MSPRT.object, 
                                termination.threshold,
                                side = 'right', theta0 = 0, 
                                Type1.target =.005, Type2.target = .2,
                                N1.max, N2.max,
                                sigma1 = 1, sigma2 = 1,
                                batch1.size, batch2.size,
                                verbose = T, plot.it = T){
  
  # side
  if(!missing(design.MSPRT.object)) side = design.MSPRT.object$side
  
  if(side!='both'){
    
    #################### two-sample z (right/left sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch1.size = design.MSPRT.object$batch1.size
      batch2.size = design.MSPRT.object$batch2.size
      N1.max = design.MSPRT.object$N1.max
      N2.max = design.MSPRT.object$N2.max
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      sigma1 = design.MSPRT.object$sigma1
      sigma2 = design.MSPRT.object$sigma2
      termination.threshold = design.MSPRT.object$termination.threshold
      theta.UMPBT = design.MSPRT.object$theta.UMPBT
      nAnalyses.max = design.MSPRT.object$nAnalyses
      
      # msg
      if(verbose){
        
        if(any(batch1.size>1)||any(batch2.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the group sequential MSPRT for a two-sample z test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the sequential MSPRT for a two-sample z test:")
          print("==========================================================================")
        }
        
        print("Group 1:")
        print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
        print(paste(" Known standard deviation: ", sigma1, sep = ""))
        print("Group 2:")
        print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
        print(paste(" Known standard deviation: ", sigma2, sep = ""))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(theta.UMPBT, 3)))
        print("-------------------------------------------------------------------------")
      }
      
      batch1.size = c(0, cumsum(batch1.size))
      batch2.size = c(0, cumsum(batch2.size))
      
      nAnalyses = min(max(which(batch1.size<=length(obs1))),
                      max(which(batch2.size<=length(obs2)))) - 1
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(obs)) print("'obs' is ignored. Not required in two-sample tests.")
      
      ## ignoring batch.seq
      if(!missing(batch.size)) print("'batch.size' is ignored. Not required in two-sample tests.")
      
      ## ignoring N.max
      if(!missing(N.max)) print("'N.max' is ignored. Not required in two-sample tests.")
      
      ## checking if length(batch1.size) and length(batch2.size) are equal
      if((!missing(batch1.size)) && (!missing(batch2.size)) &&
         (length(batch1.size)!=length(batch2.size))) return("Lenghts of batch1.size and batch2.size should be same")
      
      ## batch sizes and N for group 1
      if(missing(batch1.size)){
        
        if(missing(N1.max)){
          
          return(print("Either 'batch1.size' or 'N1.max' needs to be specified"))
          
        }else{batch1.size = rep(1, N1.max)}
        
      }else{
        
        if(missing(N1.max)){
          
          N1.max = sum(batch1.size)
          
        }else{
          
          if(sum(batch1.size)!=N1.max) return(print("Sum of batch1.size should add up to N1.max"))
        }
      }
      
      ## batch sizes and N for group 2
      if(missing(batch2.size)){
        
        if(missing(N2.max)){
          
          return(print("Either 'batch2.size' or 'N2.max' needs to be specified"))
          
        }else{batch2.size = rep(1, N2.max)}
        
      }else{
        
        if(missing(N2.max)){
          
          N2.max = sum(batch2.size)
          
        }else{
          
          if(sum(batch2.size)!=N1.max) return(print("Sum of batch2.size should add up to N2.max"))
        }
      }
      
      nAnalyses.max = length(batch1.size)
      
      ## point H0
      if(missing(theta0)) theta0 = 0
      
      ######################## UMPBT alternative ########################
      theta.UMPBT = UMPBT.alt(test.type = 'twoZ', side = side, theta0 = theta0,
                              N1 = N1.max, N2 = N2.max, Type1 = Type1.target,
                              sigma1 = sigma1, sigma2 = sigma2)
      
      # msg
      if(verbose){
        
        if(any(batch1.size>1)||any(batch2.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the group sequential MSPRT for a two-sample z test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the sequential MSPRT for a two-sample z test:")
          print("==========================================================================")
        }
        
        print("Group 1:")
        print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
        print(paste(" Known standard deviation: ", sigma1, sep = ""))
        print("Group 2:")
        print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
        print(paste(" Known standard deviation: ", sigma2, sep = ""))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
        print(paste("The UMPBT alternative is: ", round(theta.UMPBT, 3)))
        print("-------------------------------------------------------------------------")
      }
      
      batch1.size = c(0, cumsum(batch1.size))
      batch2.size = c(0, cumsum(batch2.size))
      
      nAnalyses = min(max(which(batch1.size<=length(obs1))),
                      max(which(batch2.size<=length(obs2)))) - 1
    }
    
    
    ###################### sequential comparison ######################
    # Wald's thresholds
    Accept.threshold = Type2.target/(1 - Type1.target)
    Reject.threshold = (1 - Type2.target)/Type1.target
    
    # required storages
    cumsum1_n = cumsum2_n = 0
    reached.decision = F
    rejectH0 = NA
    LR = rep(NA, nAnalyses)
    
    for(n in 1:nAnalyses){
      
      if(!reached.decision){
        
        ## sum of observations until step n
        # Group 1
        cumsum1_n = cumsum1_n + sum(obs1[(batch1.size[n]+1):batch1.size[n+1]])
        
        # Group 2
        cumsum2_n = cumsum2_n + sum(obs2[(batch2.size[n]+1):batch2.size[n+1]])
        
        # likelihood ratio of observations until step n
        LR[n] = 
          exp(-(((theta.UMPBT^2) - (theta0^2)) - 2*(theta.UMPBT - theta0)*
                  (cumsum1_n/batch1.size[n+1] - cumsum2_n/batch2.size[n+1]))/
                (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
        
        # comparing with the thresholds
        AcceptedH0_n = LR[n]<=Accept.threshold
        RejectedH0_n = LR[n]>=Reject.threshold
        reached.decision = AcceptedH0_n||RejectedH0_n
        if(reached.decision){
          
          n1 = batch1.size[n+1]
          n2 = batch2.size[n+1]
          rejectH0 = RejectedH0_n
          if(rejectH0){
            
            decision = 'reject'
            
          }else{decision = 'accept'}
        }
      }
    }
    
    # inconclusive cases
    if(!reached.decision){
      
      if(nAnalyses==nAnalyses.max){
        
        n1 = N1.max
        n2 = N2.max
        rejectH0 = LR[nAnalyses]>=termination.threshold
        
        if(rejectH0){
          
          decision = 'reject'
          
        }else if(!rejectH0){decision = 'accept'}
        reached.decision = T
        
      }else{
        
        n1 = batch1.size[nAnalyses+1]
        n2 = batch2.size[nAnalyses+1]
        decision = 'continue'
      }
    }
    
    # msg
    if(verbose==T){
      
      if(decision=='continue'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Continue sampling')
        print(paste("Sample size used: Group 1 - ", n1, 
                    ", Group 2 - ", n2, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='reject'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Reject the null hypothesis')
        print(paste("Sample size used: Group 1 - ", n1, 
                    ", Group 2 - ", n2, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='accept'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Accept the null hypothesis')
        print(paste("Sample size used: Group 1 - ", n1, 
                    ", Group 2 - ", n2, sep = ''))
        print("=========================================================================")
        cat('\n')
      }
    }
    
    
    nAnalyses.test = max(which(!is.na(LR)))
    
    ## plotting
    if(plot.it==T){
      
      # title and decision
      if(side=='right'){
        
        testname = 'Right-sided two-sample z test'
        
      }else{testname = 'Left-sided two-sample z test'}
      
      if(decision=="accept"){
        
        plot.subtitle = 
          paste('Accept the null (n1 = ', n1, 
                ', n2 = ', n2, ')', sep = '')
        
        min.LR = min(LR, na.rm = T)
        if(min.LR<Accept.threshold){
          
          ylow = min.LR
          yup = max(LR, na.rm = T)
          
        }else if((min.LR>=Accept.threshold)&&
                 (min.LR<(termination.threshold + Accept.threshold)/2)){
          
          ylow = Accept.threshold
          yup = termination.threshold
          
        }else{
          
          ylow = Accept.threshold
          yup = Reject.threshold
        }
        
      }else if(decision=="reject"){
        
        plot.subtitle = 
          paste('Reject the null (n1 = ', n1, 
                ', n2 = ', n2, ')', sep = '')
        
        ylow = Accept.threshold
        yup = max(LR, na.rm = T)
        
      }else if(decision=="continue"){
        
        plot.subtitle = 
          paste('Continue sampling (n1 = ', n1, 
                ', n2 = ', n2, ')', sep = '')
        
        last.LR = LR[max(which(!is.na(LR)))]
        if(last.LR<(termination.threshold + Accept.threshold)/2){
          
          ylow = Accept.threshold
          yup = termination.threshold
          
        }else{
          
          ylow = Accept.threshold
          yup = Reject.threshold
        }
      }
      
      df = rbind.data.frame(data.frame('xval' = 1:nAnalyses.test,
                                       'yval' = Accept.threshold,
                                       'type' = 'A'),
                            data.frame('xval' = 1:nAnalyses.test,
                                       'yval' = Reject.threshold,
                                       'type' = 'R'),
                            data.frame('xval' = 1:nAnalyses.test,
                                       'yval' = LR[1:nAnalyses.test],
                                       'type' = 'LR'),
                            data.frame('xval' = nAnalyses.max,
                                       'yval' = termination.threshold,
                                       'type' = 'term.thresh'))
      
      df$type = factor(as.character(df$type),
                       levels = c('A', 'R', 'LR', 'term.thresh'))
      
      seqcompare = ggplot(data = df,
                          aes(x = xval, y = yval, group = type)) + 
        geom_point(aes(colour = type), size = 2) +
        geom_line(aes(colour = type), size = 1) +
        geom_segment(aes(x = nAnalyses.max, y = Accept.threshold,
                         xend = nAnalyses.max, yend = termination.threshold),
                     color="forestgreen", size = 1) +
        geom_segment(aes(x = nAnalyses.max, y = Reject.threshold,
                         xend = nAnalyses.max, yend = termination.threshold),
                     color="red2", size=1) +
        geom_point(aes(x = nAnalyses.max, y = termination.threshold),
                   colour = "black", size = 2) +
        xlim(c(1, nAnalyses.test)) +
        ylim(c(ylow, yup)) +
        scale_color_manual(labels = c('Acceptance threshold', 'Rejection threshold',
                                      'Likelihood ratio', 'Termination threshold'),
                           values = c('A' = 'forestgreen', 'R' = 'red2',
                                      'LR' = 'dodgerblue', 'term.thresh' = 'black'),
                           guide = guide_legend(nrow = 2, 
                                                override.aes = list(linetype=c(1,1,1,0), 
                                                                    shape=16))) +
        theme(plot.title = element_text(size=22, face="bold"),
              plot.subtitle = element_text(size=18, face="bold"),
              axis.title.x = element_text(size=15),
              axis.title.y = element_text(size=15),
              axis.text.x = element_text(color = "black", size = 15),
              axis.text.y = element_text(color = "black", size = 15),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1.2),
              panel.background = element_blank(), legend.title=element_blank(),
              legend.key.width = unit(1.4, "cm"),
              legend.spacing.x = unit(1, 'cm'), legend.text=element_text(size=20),
              legend.position = 'bottom') +
        labs(title = testname, subtitle = plot.subtitle,
             y = 'Likelihood ratio',
             x = 'Steps in sequential analyses')
      
      suppressWarnings(print(seqcompare))
    }
    
    return(list('n1' = n1, 'n2' = n2, 'decision' = decision, 
                'Accept.threshold' = Accept.threshold, 'Reject.threshold' = Reject.threshold,
                'LR' = LR[1:nAnalyses.test], 'theta.UMPBT' = theta.UMPBT))
    
    # end one-sided twoZ
    
  }else{
    
    #################### two-sample z (both sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch1.size = design.MSPRT.object$batch1.size
      batch2.size = design.MSPRT.object$batch2.size
      N1.max = design.MSPRT.object$N1.max
      N2.max = design.MSPRT.object$N2.max
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      sigma1 = design.MSPRT.object$sigma1
      sigma2 = design.MSPRT.object$sigma2
      termination.threshold = design.MSPRT.object$termination.threshold
      theta.UMPBT = design.MSPRT.object$theta.UMPBT
      nAnalyses.max = design.MSPRT.object$nAnalyses
      
      # msg
      if(verbose){
        
        if(any(batch1.size>1)||any(batch2.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the group sequential MSPRT for a two-sample z test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the sequential MSPRT for a two-sample z test:")
          print("==========================================================================")
        }
        
        print("Group 1:")
        print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
        print(paste(" Known standard deviation: ", sigma1, sep = ""))
        print("Group 2:")
        print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
        print(paste(" Known standard deviation: ", sigma2, sep = ""))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(theta.UMPBT$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta.UMPBT$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
      batch1.size = c(0, cumsum(batch1.size))
      batch2.size = c(0, cumsum(batch2.size))
      
      nAnalyses = min(max(which(batch1.size<=length(obs1))),
                      max(which(batch2.size<=length(obs2)))) - 1
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(obs)) print("'obs' is ignored. Not required in two-sample tests.")
      
      ## ignoring batch.seq
      if(!missing(batch.size)) print("'batch.size' is ignored. Not required in two-sample tests.")
      
      ## ignoring N.max
      if(!missing(N.max)) print("'N.max' is ignored. Not required in two-sample tests.")
      
      ## checking if length(batch1.size) and length(batch2.size) are equal
      if((!missing(batch1.size)) && (!missing(batch2.size)) &&
         (length(batch1.size)!=length(batch2.size))) return("Lenghts of batch1.size and batch2.size should be same")
      
      ## batch sizes and N for group 1
      if(missing(batch1.size)){
        
        if(missing(N1.max)){
          
          return(print("Either 'batch1.size' or 'N1.max' needs to be specified"))
          
        }else{batch1.size = rep(1, N1.max)}
        
      }else{
        
        if(missing(N1.max)){
          
          N1.max = sum(batch1.size)
          
        }else{
          
          if(sum(batch1.size)!=N1.max) return(print("Sum of batch1.size should add up to N1.max"))
        }
      }
      
      ## batch sizes and N for group 2
      if(missing(batch2.size)){
        
        if(missing(N2.max)){
          
          return(print("Either 'batch2.size' or 'N2.max' needs to be specified"))
          
        }else{batch2.size = rep(1, N2.max)}
        
      }else{
        
        if(missing(N2.max)){
          
          N2.max = sum(batch2.size)
          
        }else{
          
          if(sum(batch2.size)!=N1.max) return(print("Sum of batch2.size should add up to N2.max"))
        }
      }
      
      nAnalyses.max = length(batch1.size)
      
      ## point H0
      if(missing(theta0)) theta0 = 0
      
      ######################## UMPBT alternative ########################
      theta.UMPBT = list('right' = UMPBT.alt(test.type = 'twoZ', side = 'right',
                                             theta0 = theta0, N1 = N1.max, N2 = N2.max, 
                                             Type1 = Type1.target/2,
                                             sigma1 = sigma1, sigma2 = sigma2),
                         'left' = UMPBT.alt(test.type = 'twoZ', side = 'left',
                                            theta0 = theta0, N1 = N1.max, N2 = N2.max, 
                                            Type1 = Type1.target/2,
                                            sigma1 = sigma1, sigma2 = sigma2))
      
      # msg
      if(verbose){
        
        if(any(batch1.size>1)||any(batch2.size>1)){
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the group sequential MSPRT for a two-sample z test:")
          print("==========================================================================")
          
        }else{
          
          cat('\n')
          print("==========================================================================")
          print("Implementing the sequential MSPRT for a two-sample z test:")
          print("==========================================================================")
        }
        
        print("Group 1:")
        print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
        print(paste(" Known standard deviation: ", sigma1, sep = ""))
        print("Group 2:")
        print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
        print(paste(" Known standard deviation: ", sigma2, sep = ""))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
        print("The UMPBT alternative:")
        print(paste(' On the right: ', round(theta.UMPBT$right, 3), sep = ""))
        print(paste(' On the left: ', round(theta.UMPBT$left, 3), sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
      batch1.size = c(0, cumsum(batch1.size))
      batch2.size = c(0, cumsum(batch2.size))
      
      nAnalyses = min(max(which(batch1.size<=length(obs1))),
                      max(which(batch2.size<=length(obs2)))) - 1
    }
    
    
    ###################### sequential comparison ######################
    # Wald's thresholds
    Accept.threshold = Type2.target/(1 - Type1.target/2)
    Reject.threshold = (1 - Type2.target)/(Type1.target/2)
    
    # required storages
    cumsum1_n = cumsum2_n = 0
    reached.decision.r = reached.decision.l = reached.decision = F
    rejectH0 = NA
    LR.r = LR.l = rep(NA, nAnalyses)
    
    for(n in 1:nAnalyses){
      
      if(!reached.decision){
        
        ## sum of observations until step n
        # Group 1
        cumsum1_n = cumsum1_n + sum(obs1[(batch1.size[n]+1):batch1.size[n+1]])
        
        # Group 2
        cumsum2_n = cumsum2_n + sum(obs2[(batch2.size[n]+1):batch2.size[n+1]])
        
        ## for right sided check
        if(!reached.decision.r){
          
          # likelihood ratio of observations until step n
          LR.r[n] = 
            exp(-(((theta.UMPBT$right^2) - (theta0^2)) - 2*(theta.UMPBT$right - theta0)*
                    (cumsum1_n/batch1.size[n+1] - cumsum2_n/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          # comparing with the thresholds
          AcceptedH0_n.r = LR.r[n]<=Accept.threshold
          RejectedH0_n.r = LR.r[n]>=Reject.threshold
          reached.decision.r = AcceptedH0_n.r||RejectedH0_n.r
        }
        
        ## for left sided check
        if(!reached.decision.l){
          
          # likelihood ratio of observations until step n
          LR.l[n] = 
            exp(-(((theta.UMPBT$left^2) - (theta0^2)) - 2*(theta.UMPBT$left - theta0)*
                    (cumsum1_n/batch1.size[n+1] - cumsum2_n/batch2.size[n+1]))/
                  (2*((sigma1^2)/batch1.size[n+1] + (sigma2^2)/batch2.size[n+1])))
          
          # comparing with the thresholds
          AcceptedH0_n.l = LR.l[n]<=Accept.threshold
          RejectedH0_n.l = LR.l[n]>=Reject.threshold
          reached.decision.l = AcceptedH0_n.l||RejectedH0_n.l
        }
        
        ## both-sided check
        if(AcceptedH0_n.r&&AcceptedH0_n.l){
          
          rejectH0 = F
          decision = 'accept'
          n1 = batch1.size[n+1]
          n2 = batch2.size[n+1]
          reached.decision = T
          
        }else if(RejectedH0_n.r||RejectedH0_n.l){
          
          rejectH0 = T
          decision = 'reject'
          n1 = batch1.size[n+1]
          n2 = batch2.size[n+1]
          reached.decision = T
        }
      }
    }
    
    # inconclusive cases
    if(!reached.decision){
      
      if(nAnalyses==nAnalyses.max){
        
        n1 = N1.max
        n2 = N2.max
        if(AcceptedH0_n.l&&(!reached.decision.r)){
          
          rejectH0 = LR.r[nAnalyses]>=termination.threshold
          
        }else if(AcceptedH0_n.r&&(!reached.decision.l)){
          
          rejectH0 = LR.l[nAnalyses]>=termination.threshold
          
        }else if((!reached.decision.r)&&(!reached.decision.l)){
          
          rejectH0 = max(LR.r[nAnalyses], LR.l[nAnalyses])>=termination.threshold
          
        }else{rejectH0 = F}
        
        if(rejectH0){
          
          decision = 'reject'
          
        }else if(!rejectH0){decision = 'accept'}
        reached.decision = T
        
      }else{
        
        n1 = batch1.size[nAnalyses+1]
        n2 = batch2.size[nAnalyses+1]
        decision = 'continue'
      }
    }
    
    # msg
    if(verbose==T){
      
      if(decision=='continue'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Continue sampling')
        print(paste("Sample size used: Group 1 - ", n1, 
                    ", Group 2 - ", n2, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='reject'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Reject the null hypothesis')
        print(paste("Sample size used: Group 1 - ", n1, 
                    ", Group 2 - ", n2, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='accept'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Accept the null hypothesis')
        print(paste("Sample size used: Group 1 - ", n1, 
                    ", Group 2 - ", n2, sep = ''))
        print("=========================================================================")
        cat('\n')
      }
    }
    
    
    nAnalyses.r = max(which(!is.na(LR.r)))
    nAnalyses.l = max(which(!is.na(LR.l)))
    
    ## plotting
    if(plot.it==T){
      
      # title and decision
      testname = 'Two-sided two-sample z test'
      if(decision=="accept"){
        
        plot.subtitle = 
          paste('Accept the null (n1 = ', n1, 
                ', n2 = ', n2, ')', sep = '')
        
      }else if(decision=="reject"){
        
        plot.subtitle = 
          paste('Reject the null (n1 = ', n1, 
                ', n2 = ', n2, ')', sep = '')
        
      }else if(decision=="continue"){
        
        plot.subtitle = 
          paste('Continue sampling (n1 = ', n1, 
                ', n2 = ', n2, ')', sep = '')
      }
      
      # right sided test
      if(AcceptedH0_n.r){
        
        min.LR.r = min(LR.r, na.rm = T)
        if(min.LR.r<Accept.threshold){
          
          ylow.r = min.LR.r
          yup.r = max(LR.r, na.rm = T)
          
        }else if((min.LR.r>=Accept.threshold)&&
                 (min.LR.r<(termination.threshold + Accept.threshold)/2)){
          
          ylow.r = Accept.threshold
          yup.r = termination.threshold
          
        }else{
          
          ylow.r = Accept.threshold
          yup.r = Reject.threshold
        }
        
      }else if(RejectedH0_n.r){
        
        ylow.r = Accept.threshold
        yup.r = max(LR.r, na.rm = T)
        
      }else if(!reached.decision.r){
        
        last.LR.r = LR.r[max(which(!is.na(LR.r)))]
        if(last.LR.r<(termination.threshold + Accept.threshold)/2){
          
          ylow.r = Accept.threshold
          yup.r = termination.threshold
          
        }else{
          
          ylow.r = Accept.threshold
          yup.r = Reject.threshold
        }
      }
      
      df.r = rbind.data.frame(data.frame('xval' = 1:nAnalyses.r,
                                         'yval' = Accept.threshold,
                                         'type' = 'A'),
                              data.frame('xval' = 1:nAnalyses.r,
                                         'yval' = Reject.threshold,
                                         'type' = 'R'),
                              data.frame('xval' = 1:nAnalyses.r,
                                         'yval' = LR.r[1:nAnalyses.r],
                                         'type' = 'LR'),
                              data.frame('xval' = nAnalyses.max,
                                         'yval' = termination.threshold,
                                         'type' = 'term.thresh'))
      
      df.r$type = factor(as.character(df.r$type),
                         levels = c('A', 'R', 'LR', 'term.thresh'))
      
      seqcompare.r = ggplot(data = df.r,
                            aes(x = xval, y = yval, group = type)) + 
        geom_point(aes(colour = type), size = 2) +
        geom_line(aes(colour = type), size = 1) +
        geom_segment(aes(x = nAnalyses.max, y = Accept.threshold,
                         xend = nAnalyses.max, yend = termination.threshold),
                     color="forestgreen", size = 1) +
        geom_segment(aes(x = nAnalyses.max, y = Reject.threshold,
                         xend = nAnalyses.max, yend = termination.threshold),
                     color="red2", size=1) +
        geom_point(aes(x = nAnalyses.max, y = termination.threshold),
                   colour = "black", size = 2) +
        xlim(c(1, nAnalyses.r)) +
        ylim(c(ylow.r, yup.r)) +
        scale_color_manual(labels = c('Acceptance threshold', 'Rejection threshold',
                                      'Likelihood ratio', 'Termination threshold'),
                           values = c('A' = 'forestgreen', 'R' = 'red2',
                                      'LR' = 'dodgerblue', 'term.thresh' = 'black'),
                           guide = guide_legend(nrow = 2,
                                                override.aes = list(linetype=c(1,1,1,0), 
                                                                    shape=16))) +
        theme(plot.title = element_text(size=22, face="bold"),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              axis.text.x = element_text(color = "black", size = 22),
              axis.text.y = element_text(color = "black", size = 22),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1.2),
              panel.background = element_blank(), legend.title=element_blank(),
              legend.key.width = unit(1.4, "cm"),
              legend.spacing.x = unit(1, 'cm'), legend.text=element_text(size=22),
              legend.position = 'bottom') +
        labs(title = paste('Right-sided test at level ', Type1.target/2, sep = ''),
             y = 'Likelihood ratio',
             x = 'Steps in sequential analyses')
      
      # left sided test
      if(AcceptedH0_n.l){
        
        min.LR.l = min(LR.l, na.rm = T)
        if(min.LR.l<Accept.threshold){
          
          ylow.l = min.LR.l
          yup.l = max(LR.l, na.rm = T)
          
        }else if((min.LR.l>=Accept.threshold)&&
                 (min.LR.l<(termination.threshold + Accept.threshold)/2)){
          
          ylow.l = Accept.threshold
          yup.l = termination.threshold
          
        }else{
          
          ylow.l = Accept.threshold
          yup.l = Reject.threshold
        }
        
      }else if(RejectedH0_n.l){
        
        ylow.l = Accept.threshold
        yup.l = max(LR.l, na.rm = T)
        
      }else if(!reached.decision.l){
        
        last.LR.l = LR.l[max(which(!is.na(LR.l)))]
        if(last.LR.l<(termination.threshold + Accept.threshold)/2){
          
          ylow.l = Accept.threshold
          yup.l = termination.threshold
          
        }else{
          
          ylow.l = Accept.threshold
          yup.l = Reject.threshold
        }
      }
      
      df.l = rbind.data.frame(data.frame('xval' = 1:nAnalyses.l,
                                         'yval' = Accept.threshold,
                                         'type' = 'A'),
                              data.frame('xval' = 1:nAnalyses.l,
                                         'yval' = Reject.threshold,
                                         'type' = 'R'),
                              data.frame('xval' = 1:nAnalyses.l,
                                         'yval' = LR.l[1:nAnalyses.l],
                                         'type' = 'LR'),
                              data.frame('xval' = nAnalyses.max,
                                         'yval' = termination.threshold,
                                         'type' = 'term.thresh'))
      
      df.l$type = factor(as.character(df.l$type),
                         levels = c('A', 'R', 'LR', 'term.thresh'))
      
      seqcompare.l = ggplot(data = df.l,
                            aes(x = xval, y = yval, group = type)) + 
        geom_point(aes(colour = type), size = 2) +
        geom_line(aes(colour = type), size = 1) +
        geom_segment(aes(x = nAnalyses.max, y = Accept.threshold,
                         xend = nAnalyses.max, yend = termination.threshold),
                     color="forestgreen", size = 1) +
        geom_segment(aes(x = nAnalyses.max, y = Reject.threshold,
                         xend = nAnalyses.max, yend = termination.threshold),
                     color="red2", size=1) +
        geom_point(aes(x = nAnalyses.max, y = termination.threshold),
                   colour = "black", size = 2) +
        xlim(c(1, nAnalyses.l)) +
        ylim(c(ylow.l, yup.l)) +
        scale_color_manual(labels = c('Acceptance threshold', 'Rejection threshold',
                                      'Likelihood ratio', 'Termination threshold'),
                           values = c('A' = 'forestgreen', 'R' = 'red2',
                                      'LR' = 'dodgerblue', 'term.thresh' = 'black'),
                           guide = guide_legend(nrow = 2)) +
        theme(plot.title = element_text(size=22, face="bold"),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              axis.text.x = element_text(color = "black", size = 22),
              axis.text.y = element_text(color = "black", size = 22),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1.2),
              panel.background = element_blank(), legend.title=element_blank(),
              legend.key.width = unit(1.4, "cm"),
              legend.spacing.x = unit(1, 'cm'), legend.text=element_text(size=22),
              legend.position = 'bottom') +
        labs(title = paste('Left-sided test at level ', Type1.target/2, sep = ''),
             y = 'Likelihood ratio',
             x = 'Steps in sequential analyses')
      
      suppressWarnings(print(annotate_figure(ggarrange(seqcompare.l, seqcompare.r, 
                                                       nrow = 1, ncol = 2,
                                                       legend = 'bottom', common.legend = T),
                                             top = text_grob(paste(testname, ': ', plot.subtitle, '\n'), face = "bold", 
                                                             size = 25, hjust = .5))))
    }
    
    return(list('n1' = n1, 'n2' = n2, 'decision' = decision, 
                'Accept.threshold' = Accept.threshold, 'Reject.threshold' = Reject.threshold,
                'LR' = list('right' = LR.r[1:nAnalyses.r], 'left' = LR.l[1:nAnalyses.l]),
                'theta.UMPBT' = theta.UMPBT))
    
  } # end both-sided twoZ
}


#### two-sample t test ####
implement.MSPRT_twoT = function(obs1, obs2, design.MSPRT.object, 
                                termination.threshold,
                                side = 'right', theta0 = 0, 
                                Type1.target =.005, Type2.target = .2,
                                N1.max, N2.max,
                                batch1.size, batch2.size,
                                verbose = T, plot.it = T){
  
  # side
  if(!missing(design.MSPRT.object)) side = design.MSPRT.object$side
  
  if(side!='both'){
    
    #################### two-sample t (right/left sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch1.size = design.MSPRT.object$batch1.size
      batch2.size = design.MSPRT.object$batch2.size
      N1.max = design.MSPRT.object$N1.max
      N2.max = design.MSPRT.object$N2.max
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      termination.threshold = design.MSPRT.object$termination.threshold
      nAnalyses.max = design.MSPRT.object$nAnalyses
      
      # msg
      if(verbose){
        
        if((batch1.size[1]>2)||any(batch1.size[-1]>1)||
           (batch2.size[1]>2)||any(batch2.size[-1]>1)){
          
          cat('\n')
          print("=========================================================================")
          print("Implementing the group sequential MSPRT for a two-sample t test:")
          print("=========================================================================")
          
        }else{
          
          cat('\n')
          print("=========================================================================")
          print("Implementing the sequential MSPRT for a two-sample t test:")
          print("=========================================================================")
        }
        
        print("Group 1:")
        print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
        print("Group 2:")
        print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
      batch1.size = c(0, cumsum(batch1.size))
      batch2.size = c(0, cumsum(batch2.size))
      
      nAnalyses = min(max(which(batch1.size<=length(obs1))),
                      max(which(batch2.size<=length(obs2)))) - 1
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(obs)) print("'obs' is ignored. Not required in two-sample tests.")
      
      ## ignoring batch.seq
      if(!missing(batch.size)) print("'batch.size' is ignored. Not required in two-sample tests.")
      
      ## ignoring N.max
      if(!missing(N.max)) print("'N.max' is ignored. Not required in two-sample tests.")
      
      ## checking if length(batch1.size) and length(batch2.size) are equal
      if((!missing(batch1.size)) && (!missing(batch2.size)) &&
         (length(batch1.size)!=length(batch2.size))) return("Lenghts of batch1.size and batch2.size should be same")
      
      ## batch sizes and N for group 1
      if(missing(batch1.size)){
        
        if(missing(N1.max)){
          
          return("Either 'batch1.size' or 'N1.max' needs to be specified")
          
        }else{batch1.size = c(2, rep(1, N1.max-2))}
        
      }else{
        
        if(batch1.size[1]<2){
          
          return("First batch size in Group 1 should be at least 2")
          
        }else{
          
          if(missing(N1.max)){
            
            N1.max = sum(batch1.size)
            
          }else{
            
            if(sum(batch1.size)!=N1.max) return("Sum of batch1.size should add up to N1.max")
          }
        }
      }
      
      ## batch sizes and N for group 2
      if(missing(batch2.size)){
        
        if(missing(N2.max)){
          
          return("Either 'batch2.size' or 'N2.max' needs to be specified")
          
        }else{batch2.size = c(2, rep(1, N2.max-2))}
        
      }else{
        
        if(batch2.size[1]<2){
          
          return("First batch size in Group 2 should be at least 2")
          
        }else{
          
          if(missing(N2.max)){
            
            N2.max = sum(batch2.size)
            
          }else{
            
            if(sum(batch2.size)!=N2.max) return("Sum of batch2.size should add up to N2.max")
          }
        }
      }
      
      nAnalyses.max = length(batch1.size)
      
      ## point H0
      if(missing(theta0)) theta0 = 0
      
      # msg
      if(verbose){
        
        if((batch1.size[1]>2)||any(batch1.size[-1]>1)||
           (batch2.size[1]>2)||any(batch2.size[-1]>1)){
          
          cat('\n')
          print("=========================================================================")
          print("Implementing the group sequential MSPRT for a two-sample t test:")
          print("=========================================================================")
          
        }else{
          
          cat('\n')
          print("=========================================================================")
          print("Implementing the sequential MSPRT for a two-sample t test:")
          print("=========================================================================")
        }
        
        print("Group 1:")
        print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
        print("Group 2:")
        print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
      batch1.size = c(0, cumsum(batch1.size))
      batch2.size = c(0, cumsum(batch2.size))
      
      nAnalyses = min(max(which(batch1.size<=length(obs1))),
                      max(which(batch2.size<=length(obs2)))) - 1
    }
    
    
    ###################### sequential comparison ######################
    # Wald's thresholds
    Accept.threshold = Type2.target/(1 - Type1.target)
    Reject.threshold = (1 - Type2.target)/Type1.target
    
    # cut-off (with sign) in fixed design one-sample t test
    signed_t.alpha = (2*(side=='right')-1)*
      qt(Type1.target, df = N1.max + N2.max -2, lower.tail = F)
    
    # required storages
    cumsum1_n = cumsum2_n = cumSS1_n = cumSS2_n = 0
    reached.decision = F
    rejectH0 = NA
    LR = theta.UMPBT = rep(NA, nAnalyses)
    
    for(n in 1:nAnalyses){
      
      if(!reached.decision){
        
        ## sum of observations until step n
        # Group 1
        cumsum1_n = cumsum1_n + sum(obs1[(batch1.size[n]+1):batch1.size[n+1]])
        
        # Group 2
        cumsum2_n = cumsum2_n + sum(obs2[(batch2.size[n]+1):batch2.size[n+1]])
        
        ## sum of squares of observations until step n
        # Group 1
        cumSS1_n = cumSS1_n + sum(obs1[(batch1.size[n]+1):batch1.size[n+1]]^2)
        
        # Group 2
        cumSS2_n = cumSS2_n + sum(obs2[(batch2.size[n]+1):batch2.size[n+1]]^2)
        
        ## xbar and (n-1)*(s^2) until step n
        xbar.diff_n = cumsum1_n/batch1.size[n+1] - cumsum2_n/batch2.size[n+1]
        divisor.pooled.sd_n.sq = 
          cumSS1_n - ((cumsum1_n)^2)/batch1.size[n+1] +
          cumSS2_n - ((cumsum2_n)^2)/batch2.size[n+1]
        
        # UMPBT alternative
        theta.UMPBT[n] = theta0 + signed_t.alpha*
          sqrt((divisor.pooled.sd_n.sq/(batch1.size[n+1] + batch2.size[n+1] -2))*
                 (1/N1.max + 1/N2.max))
        
        # likelihood ratio of observations until step n
        LR[n] = 
          ((1 + ((xbar.diff_n - theta0)^2)/
              (divisor.pooled.sd_n.sq*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
             (1 + ((xbar.diff_n - theta.UMPBT[n])^2)/
                (divisor.pooled.sd_n.sq*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
        
        # comparing with the thresholds
        AcceptedH0_n = LR[n]<=Accept.threshold
        RejectedH0_n = LR[n]>=Reject.threshold
        reached.decision = AcceptedH0_n||RejectedH0_n
        if(reached.decision){
          
          n1 = batch1.size[n+1]
          n2 = batch2.size[n+1]
          rejectH0 = RejectedH0_n
          if(rejectH0){
            
            decision = 'reject'
            
          }else{decision = 'accept'}
        }
      }
    }
    
    # inconclusive cases
    if(!reached.decision){
      
      if(nAnalyses==nAnalyses.max){
        
        n1 = N1.max
        n2 = N2.max
        rejectH0 = LR[nAnalyses]>=termination.threshold
        
        if(rejectH0){
          
          decision = 'reject'
          
        }else if(!rejectH0){decision = 'accept'}
        reached.decision = T
        
      }else{
        
        n1 = batch1.size[nAnalyses+1]
        n2 = batch2.size[nAnalyses+1]
        decision = 'continue'
      }
    }
    
    # msg
    if(verbose==T){
      
      if(decision=='continue'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Continue sampling')
        print(paste("Sample size used: Group 1 - ", n1, 
                    ", Group 2 - ", n2, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='reject'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Reject the null hypothesis')
        print(paste("Sample size used: Group 1 - ", n1, 
                    ", Group 2 - ", n2, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='accept'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Accept the null hypothesis')
        print(paste("Sample size used: Group 1 - ", n1, 
                    ", Group 2 - ", n2, sep = ''))
        print("=========================================================================")
        cat('\n')
      }
    }
    
    
    nAnalyses.test = max(which(!is.na(LR)))
    
    ## plotting
    if(plot.it==T){
      
      # title and decision
      if(side=='right'){
        
        testname = 'Right-sided two-sample t test'
        
      }else{testname = 'Left-sided two-sample t test'}
      
      if(decision=="accept"){
        
        plot.subtitle = 
          paste('Accept the null (n1 = ', n1, 
                ', n2 = ', n2, ')', sep = '')
        
        min.LR = min(LR, na.rm = T)
        if(min.LR<Accept.threshold){
          
          ylow = min.LR
          yup = max(LR, na.rm = T)
          
        }else if((min.LR>=Accept.threshold)&&
                 (min.LR<(termination.threshold + Accept.threshold)/2)){
          
          ylow = Accept.threshold
          yup = termination.threshold
          
        }else{
          
          ylow = Accept.threshold
          yup = Reject.threshold
        }
        
      }else if(decision=="reject"){
        
        plot.subtitle = 
          paste('Reject the null (n1 = ', n1, 
                ', n2 = ', n2, ')', sep = '')
        
        ylow = Accept.threshold
        yup = max(LR, na.rm = T)
        
      }else if(decision=="continue"){
        
        plot.subtitle = 
          paste('Continue sampling (n1 = ', n1, 
                ', n2 = ', n2, ')', sep = '')
        
        last.LR = LR[max(which(!is.na(LR)))]
        if(last.LR<(termination.threshold + Accept.threshold)/2){
          
          ylow = Accept.threshold
          yup = termination.threshold
          
        }else{
          
          ylow = Accept.threshold
          yup = Reject.threshold
        }
      }
      
      df = rbind.data.frame(data.frame('xval' = 1:nAnalyses.test,
                                       'yval' = Accept.threshold,
                                       'type' = 'A'),
                            data.frame('xval' = 1:nAnalyses.test,
                                       'yval' = Reject.threshold,
                                       'type' = 'R'),
                            data.frame('xval' = 1:nAnalyses.test,
                                       'yval' = LR[1:nAnalyses.test],
                                       'type' = 'LR'),
                            data.frame('xval' = nAnalyses.max,
                                       'yval' = termination.threshold,
                                       'type' = 'term.thresh'))
      
      df$type = factor(as.character(df$type),
                       levels = c('A', 'R', 'LR', 'term.thresh'))
      
      seqcompare = ggplot(data = df,
                          aes(x = xval, y = yval, group = type)) + 
        geom_point(aes(colour = type), size = 2) +
        geom_line(aes(colour = type), size = 1) +
        geom_segment(aes(x = nAnalyses.max, y = Accept.threshold,
                         xend = nAnalyses.max, yend = termination.threshold),
                     color="forestgreen", size = 1) +
        geom_segment(aes(x = nAnalyses.max, y = Reject.threshold,
                         xend = nAnalyses.max, yend = termination.threshold),
                     color="red2", size=1) +
        geom_point(aes(x = nAnalyses.max, y = termination.threshold),
                   colour = "black", size = 2) +
        xlim(c(1, nAnalyses.test)) +
        ylim(c(ylow, yup)) +
        scale_color_manual(labels = c('Acceptance threshold', 'Rejection threshold',
                                      'Bayes factor', 'Termination threshold'),
                           values = c('A' = 'forestgreen', 'R' = 'red2',
                                      'LR' = 'dodgerblue', 'term.thresh' = 'black'),
                           guide = guide_legend(nrow = 2, 
                                                override.aes = list(linetype=c(1,1,1,0), 
                                                                    shape=16))) +
        theme(plot.title = element_text(size=22, face="bold"),
              plot.subtitle = element_text(size=18, face="bold"),
              axis.title.x = element_text(size=15),
              axis.title.y = element_text(size=15),
              axis.text.x = element_text(color = "black", size = 15),
              axis.text.y = element_text(color = "black", size = 15),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1.2),
              panel.background = element_blank(), legend.title=element_blank(),
              legend.key.width = unit(1.4, "cm"),
              legend.spacing.x = unit(1, 'cm'), legend.text=element_text(size=20),
              legend.position = 'bottom') +
        labs(title = testname, subtitle = plot.subtitle,
             y = 'Bayes factor',
             x = 'Steps in sequential analyses')
      
      suppressWarnings(print(seqcompare))
    }
    
    return(list('n1' = n1, 'n2' = n2, 'decision' = decision, 
                'Accept.threshold' = Accept.threshold, 'Reject.threshold' = Reject.threshold,
                'LR' = LR[1:nAnalyses.test], 'theta.UMPBT' = theta.UMPBT[1:nAnalyses.test]))
    
    # end one-sided twoT
    
  }else{
    
    #################### two-sample t (both sided) ####################
    
    if(!missing(design.MSPRT.object)){
      
      batch1.size = design.MSPRT.object$batch1.size
      batch2.size = design.MSPRT.object$batch2.size
      N1.max = design.MSPRT.object$N1.max
      N2.max = design.MSPRT.object$N2.max
      Type1.target = design.MSPRT.object$Type1.target
      Type2.target = design.MSPRT.object$Type2.target
      theta0 = design.MSPRT.object$theta0
      termination.threshold = design.MSPRT.object$termination.threshold
      nAnalyses.max = design.MSPRT.object$nAnalyses
      
      # msg
      if(verbose){
        
        if((batch1.size[1]>2)||any(batch1.size[-1]>1)||
           (batch2.size[1]>2)||any(batch2.size[-1]>1)){
          
          cat('\n')
          print("=========================================================================")
          print("Implementing the group sequential MSPRT for a two-sample t test:")
          print("=========================================================================")
          
        }else{
          
          cat('\n')
          print("=========================================================================")
          print("Implementing the sequential MSPRT for a two-sample t test:")
          print("=========================================================================")
        }
        
        print("Group 1:")
        print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
        print("Group 2:")
        print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
      batch1.size = c(0, cumsum(batch1.size))
      batch2.size = c(0, cumsum(batch2.size))
      
      nAnalyses = min(max(which(batch1.size<=length(obs1))),
                      max(which(batch2.size<=length(obs2)))) - 1
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(obs)) print("'obs' is ignored. Not required in two-sample tests.")
      
      ## ignoring batch.seq
      if(!missing(batch.size)) print("'batch.size' is ignored. Not required in two-sample tests.")
      
      ## ignoring N.max
      if(!missing(N.max)) print("'N.max' is ignored. Not required in two-sample tests.")
      
      ## checking if length(batch1.size) and length(batch2.size) are equal
      if((!missing(batch1.size)) && (!missing(batch2.size)) &&
         (length(batch1.size)!=length(batch2.size))) return("Lenghts of batch1.size and batch2.size should be same")
      
      ## batch sizes and N for group 1
      if(missing(batch1.size)){
        
        if(missing(N1.max)){
          
          return("Either 'batch1.size' or 'N1.max' needs to be specified")
          
        }else{batch1.size = c(2, rep(1, N1.max-2))}
        
      }else{
        
        if(batch1.size[1]<2){
          
          return("First batch size in Group 1 should be at least 2")
          
        }else{
          
          if(missing(N1.max)){
            
            N1.max = sum(batch1.size)
            
          }else{
            
            if(sum(batch1.size)!=N1.max) return("Sum of batch1.size should add up to N1.max")
          }
        }
      }
      
      ## batch sizes and N for group 2
      if(missing(batch2.size)){
        
        if(missing(N2.max)){
          
          return("Either 'batch2.size' or 'N2.max' needs to be specified")
          
        }else{batch2.size = c(2, rep(1, N2.max-2))}
        
      }else{
        
        if(batch2.size[1]<2){
          
          return("First batch size in Group 2 should be at least 2")
          
        }else{
          
          if(missing(N2.max)){
            
            N2.max = sum(batch2.size)
            
          }else{
            
            if(sum(batch2.size)!=N2.max) return("Sum of batch2.size should add up to N2.max")
          }
        }
      }
      
      nAnalyses.max = length(batch1.size)
      
      ## point H0
      if(missing(theta0)) theta0 = 0
      
      # msg
      if(verbose){
        
        if((batch1.size[1]>2)||any(batch1.size[-1]>1)||
           (batch2.size[1]>2)||any(batch2.size[-1]>1)){
          
          cat('\n')
          print("=========================================================================")
          print("Implementing the group sequential MSPRT for a two-sample t test:")
          print("=========================================================================")
          
        }else{
          
          cat('\n')
          print("=========================================================================")
          print("Implementing the sequential MSPRT for a two-sample t test:")
          print("=========================================================================")
        }
        
        print("Group 1:")
        print(paste(" Maximum available sample sizes: ", N1.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch1.size, collapse = ', '), sep = ''))
        print("Group 2:")
        print(paste(" Maximum available sample sizes: ", N2.max, sep = ""))
        print(paste(' Batch sizes: ', paste(batch2.size, collapse = ', '), sep = ''))
        print(paste("Maximum number of sequential analyses: ", nAnalyses.max,
                    sep = ""))
        print(paste("Targeted Type I error probability: ", Type1.target, sep = ""))
        print(paste("Targeted Type II error probability: ", Type2.target, sep = ""))
        print(paste("Hypothesized value under H0: ", theta0, sep = ""))
        print(paste("Direction of the H1: ", side, sep = ""))
        print(paste("Termination threshold: ", termination.threshold, sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
      batch1.size = c(0, cumsum(batch1.size))
      batch2.size = c(0, cumsum(batch2.size))
      
      nAnalyses = min(max(which(batch1.size<=length(obs1))),
                      max(which(batch2.size<=length(obs2)))) - 1
    }
    
    
    ###################### sequential comparison ######################
    # Wald's thresholds
    Accept.threshold = Type2.target/(1 - Type1.target/2)
    Reject.threshold = (1 - Type2.target)/(Type1.target/2)
    
    # cut-off (with sign) in fixed design one-sample t test
    t.alpha = qt(Type1.target/2, df = N1.max + N2.max -2, lower.tail = F)
    
    # required storages
    cumsum1_n = cumsum2_n = cumSS1_n = cumSS2_n = 0
    reached.decision.r = reached.decision.l = reached.decision = F
    rejectH0 = NA
    LR.r = LR.l = theta.UMPBT.r = theta.UMPBT.l = rep(NA, nAnalyses)
    
    for(n in 1:nAnalyses){
      
      if(!reached.decision){
        
        ## sum of observations until step n
        # Group 1
        cumsum1_n = cumsum1_n + sum(obs1[(batch1.size[n]+1):batch1.size[n+1]])
        
        # Group 2
        cumsum2_n = cumsum2_n + sum(obs2[(batch2.size[n]+1):batch2.size[n+1]])
        
        ## sum of squares of observations until step n
        # Group 1
        cumSS1_n = cumSS1_n + sum(obs1[(batch1.size[n]+1):batch1.size[n+1]]^2)
        
        # Group 2
        cumSS2_n = cumSS2_n + sum(obs2[(batch2.size[n]+1):batch2.size[n+1]]^2)
        
        ## for right sided check
        if(!reached.decision.r){
          
          # xbar and (n-1)*(s^2) until step n
          xbar.diff_n.r = cumsum1_n/batch1.size[n+1] - cumsum2_n/batch2.size[n+1]
          divisor.pooled.sd_n.sq.r = cumSS1_n - ((cumsum1_n)^2)/batch1.size[n+1] +
            cumSS2_n - ((cumsum2_n)^2)/batch2.size[n+1]
          
          # UMPBT alternative
          theta.UMPBT.r[n] = theta0 + t.alpha*
            sqrt((divisor.pooled.sd_n.sq.r/(batch1.size[n+1] + batch2.size[n+1] -2))*
                   (1/N1.max + 1/N2.max))
          
          # likelihood ratio of observations until step n
          LR.r[n] = 
            ((1 + ((xbar.diff_n.r - theta0)^2)/
                (divisor.pooled.sd_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff_n.r - theta.UMPBT.r[n])^2)/
                  (divisor.pooled.sd_n.sq.r*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          # comparing with the thresholds
          AcceptedH0_n.r = LR.r[n]<=Accept.threshold
          RejectedH0_n.r = LR.r[n]>=Reject.threshold
          reached.decision.r = AcceptedH0_n.r||RejectedH0_n.r
        }
        
        ## for left sided check
        if(!reached.decision.l){
          
          # xbar and (n-1)*(s^2) until step n
          xbar.diff_n.l = cumsum1_n/batch1.size[n+1] - cumsum2_n/batch2.size[n+1]
          divisor.pooled.sd_n.sq.l = cumSS1_n - ((cumsum1_n)^2)/batch1.size[n+1] +
            cumSS2_n - ((cumsum2_n)^2)/batch2.size[n+1]
          
          # UMPBT alternative
          theta.UMPBT.l[n] = theta0 - t.alpha*
            sqrt((divisor.pooled.sd_n.sq.l/(batch1.size[n+1] + batch2.size[n+1] -2))*
                   (1/N1.max + 1/N2.max))
          
          # likelihood ratio of observations until step n
          LR.l[n] = 
            ((1 + ((xbar.diff_n.l - theta0)^2)/
                (divisor.pooled.sd_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1])))/
               (1 + ((xbar.diff_n.l - theta.UMPBT.l[n])^2)/
                  (divisor.pooled.sd_n.sq.l*(1/batch1.size[n+1] + 1/batch2.size[n+1]))))^((batch1.size[n+1] + batch2.size[n+1])/2)
          
          # comparing with the thresholds
          AcceptedH0_n.l = LR.l[n]<=Accept.threshold
          RejectedH0_n.l = LR.l[n]>=Reject.threshold
          reached.decision.l = AcceptedH0_n.l||RejectedH0_n.l
        }
        
        ## both-sided check
        if(AcceptedH0_n.r&&AcceptedH0_n.l){
          
          rejectH0 = F
          decision = 'accept'
          n1 = batch1.size[n+1]
          n2 = batch2.size[n+1]
          reached.decision = T
          
        }else if(RejectedH0_n.r||RejectedH0_n.l){
          
          rejectH0 = T
          decision = 'reject'
          n1 = batch1.size[n+1]
          n2 = batch2.size[n+1]
          reached.decision = T
        }
      }
    }
    
    # inconclusive cases
    if(!reached.decision){
      
      if(nAnalyses==nAnalyses.max){
        
        n1 = N1.max
        n2 = N2.max
        if(AcceptedH0_n.l&&(!reached.decision.r)){
          
          rejectH0 = LR.r[nAnalyses]>=termination.threshold
          
        }else if(AcceptedH0_n.r&&(!reached.decision.l)){
          
          rejectH0 = LR.l[nAnalyses]>=termination.threshold
          
        }else if((!reached.decision.r)&&(!reached.decision.l)){
          
          rejectH0 = max(LR.r[nAnalyses], LR.l[nAnalyses])>=termination.threshold
          
        }else{rejectH0 = F}
        
        if(rejectH0){
          
          decision = 'reject'
          
        }else if(!rejectH0){decision = 'accept'}
        reached.decision = T
        
      }else{
        
        n1 = batch1.size[nAnalyses+1]
        n2 = batch2.size[nAnalyses+1]
        decision = 'continue'
      }
    }
    
    # msg
    if(verbose==T){
      
      if(decision=='continue'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Continue sampling')
        print(paste("Sample size used: Group 1 - ", n1, 
                    ", Group 2 - ", n2, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='reject'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Reject the null hypothesis')
        print(paste("Sample size used: Group 1 - ", n1, 
                    ", Group 2 - ", n2, sep = ''))
        print("=========================================================================")
        cat('\n')
        
      }else if(decision=='accept'){
        
        cat('\n')
        print("=========================================================================")
        print("Sequential comparison summary:")
        print("=========================================================================")
        print('Decision: Accept the null hypothesis')
        print(paste("Sample size used: Group 1 - ", n1, 
                    ", Group 2 - ", n2, sep = ''))
        print("=========================================================================")
        cat('\n')
      }
    }
    
    
    nAnalyses.r = max(which(!is.na(LR.r)))
    nAnalyses.l = max(which(!is.na(LR.l)))
    
    ## plotting
    if(plot.it==T){
      
      # title and decision
      testname = 'Two-sided two-sample t test'
      if(decision=="accept"){
        
        plot.subtitle = 
          paste('Accept the null (n1 = ', n1, 
                ', n2 = ', n2, ')', sep = '')
        
      }else if(decision=="reject"){
        
        plot.subtitle = 
          paste('Reject the null (n1 = ', n1, 
                ', n2 = ', n2, ')', sep = '')
        
      }else if(decision=="continue"){
        
        plot.subtitle = 
          paste('Continue sampling (n1 = ', n1, 
                ', n2 = ', n2, ')', sep = '')
      }
      
      # right sided test
      if(AcceptedH0_n.r){
        
        min.LR.r = min(LR.r, na.rm = T)
        if(min.LR.r<Accept.threshold){
          
          ylow.r = min.LR.r
          yup.r = max(LR.r, na.rm = T)
          
        }else if((min.LR.r>=Accept.threshold)&&
                 (min.LR.r<(termination.threshold + Accept.threshold)/2)){
          
          ylow.r = Accept.threshold
          yup.r = termination.threshold
          
        }else{
          
          ylow.r = Accept.threshold
          yup.r = Reject.threshold
        }
        
      }else if(RejectedH0_n.r){
        
        ylow.r = Accept.threshold
        yup.r = max(LR.r, na.rm = T)
        
      }else if(!reached.decision.r){
        
        last.LR.r = LR.r[max(which(!is.na(LR.r)))]
        if(last.LR.r<(termination.threshold + Accept.threshold)/2){
          
          ylow.r = Accept.threshold
          yup.r = termination.threshold
          
        }else{
          
          ylow.r = Accept.threshold
          yup.r = Reject.threshold
        }
      }
      
      df.r = rbind.data.frame(data.frame('xval' = 1:nAnalyses.r,
                                         'yval' = Accept.threshold,
                                         'type' = 'A'),
                              data.frame('xval' = 1:nAnalyses.r,
                                         'yval' = Reject.threshold,
                                         'type' = 'R'),
                              data.frame('xval' = 1:nAnalyses.r,
                                         'yval' = LR.r[1:nAnalyses.r],
                                         'type' = 'LR'),
                              data.frame('xval' = nAnalyses.max,
                                         'yval' = termination.threshold,
                                         'type' = 'term.thresh'))
      
      df.r$type = factor(as.character(df.r$type),
                         levels = c('A', 'R', 'LR', 'term.thresh'))
      
      seqcompare.r = ggplot(data = df.r,
                            aes(x = xval, y = yval, group = type)) + 
        geom_point(aes(colour = type), size = 2) +
        geom_line(aes(colour = type), size = 1) +
        geom_segment(aes(x = nAnalyses.max, y = Accept.threshold,
                         xend = nAnalyses.max, yend = termination.threshold),
                     color="forestgreen", size = 1) +
        geom_segment(aes(x = nAnalyses.max, y = Reject.threshold,
                         xend = nAnalyses.max, yend = termination.threshold),
                     color="red2", size=1) +
        geom_point(aes(x = nAnalyses.max, y = termination.threshold),
                   colour = "black", size = 2) +
        xlim(c(1, nAnalyses.r)) +
        ylim(c(ylow.r, yup.r)) +
        scale_color_manual(labels = c('Acceptance threshold', 'Rejection threshold',
                                      'Bayes factor', 'Termination threshold'),
                           values = c('A' = 'forestgreen', 'R' = 'red2',
                                      'LR' = 'dodgerblue', 'term.thresh' = 'black'),
                           guide = guide_legend(nrow = 2,
                                                override.aes = list(linetype=c(1,1,1,0), 
                                                                    shape=16))) +
        theme(plot.title = element_text(size=22, face="bold"),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              axis.text.x = element_text(color = "black", size = 22),
              axis.text.y = element_text(color = "black", size = 22),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1.2),
              panel.background = element_blank(), legend.title=element_blank(),
              legend.key.width = unit(1.4, "cm"),
              legend.spacing.x = unit(1, 'cm'), legend.text=element_text(size=22),
              legend.position = 'bottom') +
        labs(title = paste('Right-sided test at level ', Type1.target/2, sep = ''),
             y = 'Bayes factor',
             x = 'Steps in sequential analyses')
      
      # left sided test
      if(AcceptedH0_n.l){
        
        min.LR.l = min(LR.l, na.rm = T)
        if(min.LR.l<Accept.threshold){
          
          ylow.l = min.LR.l
          yup.l = max(LR.l, na.rm = T)
          
        }else if((min.LR.l>=Accept.threshold)&&
                 (min.LR.l<(termination.threshold + Accept.threshold)/2)){
          
          ylow.l = Accept.threshold
          yup.l = termination.threshold
          
        }else{
          
          ylow.l = Accept.threshold
          yup.l = Reject.threshold
        }
        
      }else if(RejectedH0_n.l){
        
        ylow.l = Accept.threshold
        yup.l = max(LR.l, na.rm = T)
        
      }else if(!reached.decision.l){
        
        last.LR.l = LR.l[max(which(!is.na(LR.l)))]
        if(last.LR.l<(termination.threshold + Accept.threshold)/2){
          
          ylow.l = Accept.threshold
          yup.l = termination.threshold
          
        }else{
          
          ylow.l = Accept.threshold
          yup.l = Reject.threshold
        }
      }
      
      df.l = rbind.data.frame(data.frame('xval' = 1:nAnalyses.l,
                                         'yval' = Accept.threshold,
                                         'type' = 'A'),
                              data.frame('xval' = 1:nAnalyses.l,
                                         'yval' = Reject.threshold,
                                         'type' = 'R'),
                              data.frame('xval' = 1:nAnalyses.l,
                                         'yval' = LR.l[1:nAnalyses.l],
                                         'type' = 'LR'),
                              data.frame('xval' = nAnalyses.max,
                                         'yval' = termination.threshold,
                                         'type' = 'term.thresh'))
      
      df.l$type = factor(as.character(df.l$type),
                         levels = c('A', 'R', 'LR', 'term.thresh'))
      
      seqcompare.l = ggplot(data = df.l,
                            aes(x = xval, y = yval, group = type)) + 
        geom_point(aes(colour = type), size = 2) +
        geom_line(aes(colour = type), size = 1) +
        geom_segment(aes(x = nAnalyses.max, y = Accept.threshold,
                         xend = nAnalyses.max, yend = termination.threshold),
                     color="forestgreen", size = 1) +
        geom_segment(aes(x = nAnalyses.max, y = Reject.threshold,
                         xend = nAnalyses.max, yend = termination.threshold),
                     color="red2", size=1) +
        geom_point(aes(x = nAnalyses.max, y = termination.threshold),
                   colour = "black", size = 2) +
        xlim(c(1, nAnalyses.l)) +
        ylim(c(ylow.l, yup.l)) +
        scale_color_manual(labels = c('Acceptance threshold', 'Rejection threshold',
                                      'Bayes factor', 'Termination threshold'),
                           values = c('A' = 'forestgreen', 'R' = 'red2',
                                      'LR' = 'dodgerblue', 'term.thresh' = 'black'),
                           guide = guide_legend(nrow = 2)) +
        theme(plot.title = element_text(size=22, face="bold"),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              axis.text.x = element_text(color = "black", size = 22),
              axis.text.y = element_text(color = "black", size = 22),
              panel.border = element_rect(linetype = "solid", fill = NA, size = 1.2),
              panel.background = element_blank(), legend.title=element_blank(),
              legend.key.width = unit(1.4, "cm"),
              legend.spacing.x = unit(1, 'cm'), legend.text=element_text(size=22),
              legend.position = 'bottom') +
        labs(title = paste('Left-sided test at level ', Type1.target/2, sep = ''),
             y = 'Bayes factor',
             x = 'Steps in sequential analyses')
      
      suppressWarnings(print(annotate_figure(ggarrange(seqcompare.l, seqcompare.r, 
                                                       nrow = 1, ncol = 2,
                                                       legend = 'bottom', common.legend = T),
                                             top = text_grob(paste(testname, ': ', plot.subtitle, '\n'), face = "bold", 
                                                             size = 25, hjust = .5))))
    }
    
    return(list('n1' = n1, 'n2' = n2, 'decision' = decision, 
                'Accept.threshold' = Accept.threshold, 'Reject.threshold' = Reject.threshold,
                'LR' = list('right' = LR.r[1:nAnalyses.r], 'left' = LR.l[1:nAnalyses.l]),
                'theta.UMPBT' = list('right' = theta.UMPBT.r[1:nAnalyses.r],
                                     'left' = theta.UMPBT.l[1:nAnalyses.l])))
    
  } # end both-sided twoT
}


#### implementation of the MSPRT combined for all ####
implement.MSPRT = function(obs, obs1, obs2, design.MSPRT.object, 
                           termination.threshold, test.type, 
                           side = 'right', theta0, 
                           Type1.target =.005, Type2.target = .2,
                           N.max, N1.max, N2.max,
                           sigma = 1, sigma1 = 1, sigma2 = 1,
                           batch.size, batch1.size, batch2.size,
                           verbose = T, plot.it = T){
  
  if(!missing(design.MSPRT.object)) test.type = design.MSPRT.object$test.type
  
  if(test.type=='oneProp'){
    
    if(!missing(design.MSPRT.object)){
      
      return(implement.MSPRT_oneProp(obs = obs, design.MSPRT.object = design.MSPRT.object,
                                     verbose = verbose, plot.it = plot.it))
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(obs1)) print("'obs1' is ignored. Not required in one-sample tests.")
      if(!missing(obs2)) print("'obs2' is ignored. Not required in one-sample tests.")
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
      if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
      
      ## ignoring N1.max & N2.max
      if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
      if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
      
      ## batch sizes and N.max
      if(missing(batch.size)){
        
        if(missing(N.max)){
          
          return("Either 'batch.size' or 'N.max' needs to be specified")
          
        }else{batch.size = rep(1, N.max)}
        
      }else{
        
        if(missing(N.max)){
          
          N.max = sum(batch.size)
          
        }else{
          
          if(sum(batch.size)!=N.max) return("Sum of batch sizes should add up to N.max")
        }
      }
      
      ## point H0
      if(missing(theta0)) theta0 = 0.5
      
      return(implement.MSPRT_oneProp(obs = obs, 
                                     termination.threshold = termination.threshold,
                                     side = side, theta0 = theta0, 
                                     Type1.target = Type1.target,
                                     Type2.target = Type2.target,
                                     N.max = N.max, batch.size = batch.size,
                                     verbose = verbose, plot.it = plot.it))
    }
    
  }else if(test.type=='oneZ'){
    
    if(!missing(design.MSPRT.object)){
      
      return(implement.MSPRT_oneZ(obs = obs, design.MSPRT.object = design.MSPRT.object,
                                  verbose = verbose, plot.it = plot.it))
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(obs1)) print("'obs1' is ignored. Not required in one-sample tests.")
      if(!missing(obs2)) print("'obs2' is ignored. Not required in one-sample tests.")
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
      if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
      
      ## ignoring N1.max & N2.max
      if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
      if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
      
      ## batch sizes and N.max
      if(missing(batch.size)){
        
        if(missing(N.max)){
          
          return("Either 'batch.size' or 'N.max' needs to be specified")
          
        }else{batch.size = rep(1, N.max)}
        
      }else{
        
        if(missing(N.max)){
          
          N.max = sum(batch.size)
          
        }else{
          
          if(sum(batch.size)!=N.max) return("Sum of batch sizes should add up to N.max")
        }
      }
      
      ## point H0
      if(missing(theta0)) theta0 = 0
      
      return(implement.MSPRT_oneZ(obs = obs, 
                                  termination.threshold = termination.threshold,
                                  side = side, theta0 = theta0, sigma = sigma,
                                  Type1.target = Type1.target,
                                  Type2.target = Type2.target,
                                  N.max = N.max, batch.size = batch.size,
                                  verbose = verbose, plot.it = plot.it))
    }
    
  }else if(test.type=='oneT'){
    
    if(!missing(design.MSPRT.object)){
      
      return(implement.MSPRT_oneT(obs = obs, design.MSPRT.object = design.MSPRT.object,
                                  verbose = verbose, plot.it = plot.it))
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(obs1)) print("'obs1' is ignored. Not required in one-sample tests.")
      if(!missing(obs2)) print("'obs2' is ignored. Not required in one-sample tests.")
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(batch1.size)) print("'batch1.size' is ignored. Not required in one-sample tests.")
      if(!missing(batch2.size)) print("'batch2.size' is ignored. Not required in one-sample tests.")
      
      ## ignoring N1.max & N2.max
      if(!missing(N1.max)) print("'N1.max' is ignored. Not required in one-sample tests.")
      if(!missing(N2.max)) print("'N2.max' is ignored. Not required in one-sample tests.")
      
      ## batch sizes and N.max
      if(missing(batch.size)){
        
        if(missing(N.max)){
          
          return("Either 'batch.size' or 'N.max' needs to be specified")
          
        }else{batch.size = c(2, rep(1, N.max-2))}
        
      }else{
        
        if(batch.size[1]<2){
          
          return("First batch size should be at least 2")
          
        }else{
          
          if(missing(N.max)){
            
            N.max = sum(batch.size)
            
          }else{
            
            if(sum(batch.size)!=N.max) return("Sum of batch.size should add up to N.max")
          }
        }
      }
      
      ## point H0
      if(missing(theta0)) theta0 = 0
      
      return(implement.MSPRT_oneT(obs = obs, 
                                  termination.threshold = termination.threshold,
                                  side = side, theta0 = theta0,
                                  Type1.target = Type1.target,
                                  Type2.target = Type2.target,
                                  N.max = N.max, batch.size = batch.size,
                                  verbose = verbose, plot.it = plot.it))
    }
    
  }else if(test.type=='twoZ'){
    
    if(!missing(design.MSPRT.object)){
      
      return(implement.MSPRT_twoZ(obs1 = obs1, obs2 = obs2,
                                  design.MSPRT.object = design.MSPRT.object,
                                  verbose = verbose, plot.it = plot.it))
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(obs)) print("'obs' is ignored. Not required in two-sample tests.")
      
      ## ignoring batch.seq
      if(!missing(batch.size)) print("'batch.size' is ignored. Not required in two-sample tests.")
      
      ## ignoring N.max
      if(!missing(N.max)) print("'N.max' is ignored. Not required in two-sample tests.")
      
      ## checking if length(batch1.size) and length(batch2.size) are equal
      if((!missing(batch1.size)) && (!missing(batch2.size)) &&
         (length(batch1.size)!=length(batch2.size))) return("Lenghts of batch1.size and batch2.size should be same")
      
      ## batch sizes and N for group 1
      if(missing(batch1.size)){
        
        if(missing(N1.max)){
          
          return(print("Either 'batch1.size' or 'N1.max' needs to be specified"))
          
        }else{batch1.size = rep(1, N1.max)}
        
      }else{
        
        if(missing(N1.max)){
          
          N1.max = sum(batch1.size)
          
        }else{
          
          if(sum(batch1.size)!=N1.max) return(print("Sum of batch1.size should add up to N1.max"))
        }
      }
      
      ## batch sizes and N for group 2
      if(missing(batch2.size)){
        
        if(missing(N2.max)){
          
          return(print("Either 'batch2.size' or 'N2.max' needs to be specified"))
          
        }else{batch2.size = rep(1, N2.max)}
        
      }else{
        
        if(missing(N2.max)){
          
          N2.max = sum(batch2.size)
          
        }else{
          
          if(sum(batch2.size)!=N1.max) return(print("Sum of batch2.size should add up to N2.max"))
        }
      }
      
      ## point H0
      if(missing(theta0)) theta0 = 0
      
      return(implement.MSPRT_twoZ(obs1 = obs1, obs2 = obs2,
                                  termination.threshold = termination.threshold,
                                  side = side, theta0 = theta0, 
                                  Type1.target = Type1.target,
                                  Type2.target = Type2.target,
                                  N1.max = N1.max, N2.max = N2.max,
                                  sigma1 = sigma1, sigma2 = sigma2,
                                  batch1.size = batch1.size, batch2.size = batch2.size,
                                  verbose = verbose, plot.it = plot.it))
    }
    
  }else if(test.type=='twoT'){
    
    if(!missing(design.MSPRT.object)){
      
      return(implement.MSPRT_twoT(obs1 = obs1, obs2 = obs2,
                                  design.MSPRT.object = design.MSPRT.object,
                                  verbose = verbose, plot.it = plot.it))
      
    }else{
      
      ## ignoring batch1.seq & batch2.seq
      if(!missing(obs)) print("'obs' is ignored. Not required in two-sample tests.")
      
      ## ignoring batch.seq
      if(!missing(batch.size)) print("'batch.size' is ignored. Not required in two-sample tests.")
      
      ## ignoring N.max
      if(!missing(N.max)) print("'N.max' is ignored. Not required in two-sample tests.")
      
      ## checking if length(batch1.size) and length(batch2.size) are equal
      if((!missing(batch1.size)) && (!missing(batch2.size)) &&
         (length(batch1.size)!=length(batch2.size))) return("Lenghts of batch1.size and batch2.size should be same")
      
      ## batch sizes and N for group 1
      if(missing(batch1.size)){
        
        if(missing(N1.max)){
          
          return("Either 'batch1.size' or 'N1.max' needs to be specified")
          
        }else{batch1.size = c(2, rep(1, N1.max-2))}
        
      }else{
        
        if(batch1.size[1]<2){
          
          return("First batch size in Group 1 should be at least 2")
          
        }else{
          
          if(missing(N1.max)){
            
            N1.max = sum(batch1.size)
            
          }else{
            
            if(sum(batch1.size)!=N1.max) return("Sum of batch1.size should add up to N1.max")
          }
        }
      }
      
      ## batch sizes and N for group 2
      if(missing(batch2.size)){
        
        if(missing(N2.max)){
          
          return("Either 'batch2.size' or 'N2.max' needs to be specified")
          
        }else{batch2.size = c(2, rep(1, N2.max-2))}
        
      }else{
        
        if(batch2.size[1]<2){
          
          return("First batch size in Group 2 should be at least 2")
          
        }else{
          
          if(missing(N2.max)){
            
            N2.max = sum(batch2.size)
            
          }else{
            
            if(sum(batch2.size)!=N2.max) return("Sum of batch2.size should add up to N2.max")
          }
        }
      }
      
      ## point H0
      if(missing(theta0)) theta0 = 0
      
      return(implement.MSPRT_twoT(obs1 = obs1, obs2 = obs2,
                                  termination.threshold = termination.threshold,
                                  side = side, theta0 = theta0, 
                                  Type1.target = Type1.target,
                                  Type2.target = Type2.target,
                                  N1.max = N1.max, N2.max = N2.max,
                                  batch1.size = batch1.size, batch2.size = batch2.size,
                                  verbose = verbose, plot.it = plot.it))
    }
  }
}






