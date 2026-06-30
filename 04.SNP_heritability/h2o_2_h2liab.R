#### heritabity from obs to libility scale
h2o <- 0.176

pop_prev <- 0.05459509 # 53274/975802 CHIP

samp_prev <- 0.1105457

liab_T  <- 0.05 # median VAF>=0.05

h2_obs_2_liab <- function(h2o, pop_prev, 
                          samp_prev, liab_T){
  
  liab_thresh <- dnorm(qnorm(liab_T))
  
  h2_liab <- ( (h2o * pop_prev^2 * (1-pop_prev)^2 )/liab_thresh^2 ) / (samp_prev * (1-samp_prev) ) 
  
  return(h2_liab)
  
}
## 
h2_obs_2_liab(h2o = h2o, 
              pop_prev = pop_prev, 
              samp_prev = samp_prev, 
              liab_T =  liab_T)

