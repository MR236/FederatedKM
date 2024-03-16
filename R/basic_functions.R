#' Combines rmutil::int and stats::integrate to quickly and accurately integrate steep one-dimensional functions with vector-valued upper limit
#' @param f function
#' @param lower lower limit (one value)
#' @param upper upper limits (vector-valued)
#' @param disconts jump discontinuities (ignored if < lower or > pivot), increased accuracy if provided
#' @param pivot integration below pivot performed by Gauss-Konrad quadrature (accurate but expensive), above pivot by Romberg integration (fast but inaccurate for steep functions)
#' @return vector of definite integral values
#' @export
hybrid_integration <- function(f, lower, upper, disconts=NULL, pivot=lower) {
  if (pivot == lower) {return(int(f, lower, upper))}
  ord <- order(upper)
  upper <- upper[ord] # upper limits in ascending order
  disconts <- disconts[disconts<pivot & disconts>lower] # only relevant jump discontinuities
  partA <- c(lower, upper[upper<=pivot], pivot, disconts) # filter for upper lims before pivot and insert discontinuities
  ord_lower <- order(partA)
  partA <- partA[ord_lower] # ascending order
  integrals <- NULL
  integral <- 0
  for (i in seq(2, length(partA))) {
    integral <- integral + integrate(f, partA[i-1], partA[i])$value 
    integrals <- c(integrals, integral)
  } # compute integrals below pivot in sequential intervals
  integrals <- integrals[order(ord_lower[-c(1)])] # restore ordering
  low_part <- integrals[length(integrals)-length(disconts)] # integral from lower to pivot
  if (length(integrals) > 1) {
    integrals <- integrals[1:(length(integrals)-length(disconts) - 1)]
  } # remove the pivot and discontinuity value from the list of integrals
  else {integrals <- NULL} 
  partB <- upper[upper>pivot]
  if (length(partB) > 0) {
  PB_int <- rmutil::int(f, pivot, partB) + c(low_part)} # integrate above pivot with vectorized Romberg integration
  else {PB_int <- NULL}
  result <- c(integrals, PB_int) # rejoin lower and upper integrals
  return(result[order(ord)]) # return integrals in original ordering of upper limits
}


#' Computes times when survival function drops below given quantile values
#' @param t vector of time points
#' @param s_t corresponding vector of survival probabilites from KM estimate
#' @param qs quantiles
#' @return Vector of time points based on quantiles
#' @export
get_time_quantiles <- function(t, s_t, qs) {
  quantiles <- NULL
  for (i in qs) {
    quantiles <- c(quantiles, max(t[which(s_t >= i)]))
  }
  return(quantiles)
}

#' Expit function
#' @param eta
#' @return expit(eta)
#' @export
expit <- function(eta){return(exp(eta)/(1 + exp(eta)))}

#' Inverse of cloglog link function
#' @param eta
#' @return inverse cloglog(eta)
#' @export
cloglog_inv <- function(eta){return(1 - exp(-exp(eta)))}

#' Derivate of inverse of cloglog link function
#' @param eta
#' @return derivative of inverse cloglog(eta) wrt eta
#' @export
cloglog_inv_deriv <- function(eta){return(exp(eta-exp(eta)))}


#' Survival function at given times from I-spline parameters, cloglog link
#' @param t vector of times
#' @param knots knot locations
#' @param degree spline degrees of freedom
#' @param boundary_knots boundary knot locations
#' @param coefficients estimated spline coefficients
#' @return estimated survival probabilities at given times
#' @export
survival_function_cloglog.v <- function(t, knots, degree, boundary_knots, coefficients, restriction = 0) {
  final_knot <- boundary_knots[2]
  first_knot <- max(boundary_knots[1], restriction)
  basis_T1 <- c(1, iSpline(first_knot, knots=knots, degree=degree, Boundary.knots=boundary_knots))
  S_T1 <- cloglog_inv(coefficients %*% basis_T1)
  lambda_T1 <- -log(S_T1) / first_knot
  basis_TF <- c(1, iSpline(final_knot, knots=knots, degree=degree, Boundary.knots=boundary_knots))
  dbasis_TF <- c(0, mSpline(final_knot, knots=knots, degree=degree, Boundary.knots=boundary_knots, intercept=TRUE))
  S_TF <- cloglog_inv(coefficients %*% basis_TF)
  term1_TF <- cloglog_inv_deriv(coefficients %*% basis_TF)
  term2_TF <- coefficients %*% dbasis_TF
  haz_TF <- -term1_TF*term2_TF/S_TF
  haz_TF <- ifelse(haz_TF <= 0, 0, haz_TF)
  basis_all <-  cbind(rep(1, length(t)), iSpline(t, knots=knots, degree=degree, Boundary.knots = boundary_knots))
  est_1 <- exp(-t*c(lambda_T1))
  est_2 <- exp((final_knot-t)*c(haz_TF))*c(S_TF)
  est_3 <- cloglog_inv(coefficients %*% t(basis_all))
  final_est <- ifelse(t < first_knot, est_1, ifelse(t > final_knot, est_2, est_3))
  return(final_est)
}


#' Hazard function at given times from I-spline parameters, cloglog link
#' @param t vector of times
#' @param knots knot locations
#' @param degree spline degrees of freedom
#' @param boundary_knots boundary knot locations
#' @param coefficients estimated spline coefficients
#' @return estimated hazard at given times
#' @export
hazard_function_cloglog.v <- function(t, knots, degree, boundary_knots, coefficients, restriction=0) {
  final_knot <- boundary_knots[2]
  first_knot <- max(boundary_knots[1], restriction)
  basis_T1 <- c(1, iSpline(first_knot, knots=knots, degree=degree, Boundary.knots=boundary_knots))
  S_T1 <- cloglog_inv(coefficients %*% basis_T1)
  lambda_T1 <- -log(S_T1) / first_knot
  basis_TF <- c(1, iSpline(final_knot, knots=knots, degree=degree, Boundary.knots=boundary_knots))
  dbasis_TF <- c(0, mSpline(final_knot, knots=knots, degree=degree, Boundary.knots=boundary_knots, intercept=TRUE))
  S_TF <- cloglog_inv(coefficients %*% basis_TF)
  term1_TF <- cloglog_inv_deriv(coefficients %*% basis_TF)
  term2_TF <- coefficients %*% dbasis_TF
  haz_TF <- -term1_TF*term2_TF/S_TF
  haz_TF <- ifelse(haz_TF <= 0, 0, haz_TF)
  S_T_all <- survival_function_cloglog.v(t, knots, degree, boundary_knots, coefficients)
  dbasis <- cbind(rep(0, length(t)), mSpline(t, knots=knots, degree=degree, Boundary.knots = boundary_knots, intercept=TRUE, periodic=FALSE))
  basis <- cbind(rep(1, length(t)), iSpline(t, knots=knots, degree=degree, Boundary.knots = boundary_knots))
  term1 <- cloglog_inv_deriv(coefficients %*% t(basis))
  term2 <- coefficients %*% t(dbasis)
  haz <- -term1*term2/S_T_all
  haz <- ifelse(haz <= 0, 0, haz)
  final_haz <- ifelse(t < first_knot, lambda_T1, ifelse(t > final_knot, haz_TF, haz))
  return(final_haz)
}

#' Variance of survival function at given times from I-spline parameters, cloglog link
#' @param t_k vector of times
#' @param spline_params spline parameters, containing knots, degrees of freedom, boundary knots and coefficients for both at-risk and event functions
#' @param restriction pivot for hybrid integration (should be >0 if steep hazard observed near 0)
#' @return Variance of cloglog(S(t)) for each t in t_k
#' @export
variance_cloglog <- function(t_k, spline_params, restriction=0) {
  PHF.v <- function(t) {hazard_function_cloglog.v(t, spline_params$Event$Knots,
                                                  spline_params$Event$Degree, spline_params$Event$BKnots,
                                                  spline_params$Event$Coefficients)}
  
  PSF.v <- function(t) {survival_function_cloglog.v(t, spline_params$Event$Knots,
                                                    spline_params$Event$Degree, spline_params$Event$BKnots,
                                                    spline_params$Event$Coefficients)}
  PRF.v <- function(t) {survival_function_cloglog.v(t, spline_params$Risk$Knots,
                                                    spline_params$Risk$Degree, spline_params$Risk$BKnots,
                                                    spline_params$Risk$Coefficients)}
  term1 <- 1/(c(spline_params$N)*log(PSF.v(t_k))^2)
  integrand <- function(t) {
    return(PHF.v(t)/(PRF.v(t)*(1-PHF.v(t))))
  }
  term2 <- hybrid_integration(integrand, 0, t_k, disconts = NULL, pivot=restriction)
  return(term1*term2)
}


#' Numerator of the two-sample log-rank test (observed-expected) based on spline parameters 
#' @param spline_params_all spline parameters in the full population, containing knots, degrees of freedom, boundary knots and coefficients for both at-risk and event functions
#' @param spline_params_group1 spline parameters in one of the groups, containing knots, degrees of freedom, boundary knots and coefficients for both at-risk and event functions
#' @param restriction pivot for hybrid integration (should be >0 if steep hazard observed near 0)
#' @param limit integration limit for log-rank test to prevent error due to extrapolation, a good choice is the maximum value prior to which we have events in both groups
#' @return Numerator of log-rank test
#' @export
log_rank_numerator <- function(spline_params_all, spline_params_group1, restriction=0, limit=spline_params_all$Max*0.99) {
  PHF.v <- function(t) {hazard_function_cloglog.v(t, spline_params_all$Event$Knots,
                                                  spline_params_all$Event$Degree, spline_params_all$Event$BKnots,
                                                  spline_params_all$Event$Coefficients)}
  PRF.v <- function(t) {survival_function_cloglog.v(t, spline_params_all$Risk$Knots,
                                                    spline_params_all$Risk$Degree, spline_params_all$Risk$BKnots,
                                                    spline_params_all$Risk$Coefficients)}
  PSF1.v <- function(t) {survival_function_cloglog.v(t, spline_params_group1$Event$Knots,
                                                     spline_params_group1$Event$Degree, spline_params_group1$Event$BKnots,
                                                     spline_params_group1$Event$Coefficients)}
  PHF1.v <- function(t) {hazard_function_cloglog.v(t, spline_params_group1$Event$Knots,
                                                   spline_params_group1$Event$Degree, spline_params_group1$Event$BKnots,
                                                   spline_params_group1$Event$Coefficients)}
  PRF1.v <- function(t) {survival_function_cloglog.v(t, spline_params_group1$Risk$Knots,
                                                     spline_params_group1$Risk$Degree, spline_params_group1$Risk$BKnots,
                                                     spline_params_group1$Risk$Coefficients)}
  risk_integrand <- function(t) {
    return(PRF1.v(t)*PHF.v(t))
  }
  observed_integrand <- function(t){
    return(PRF1.v(t)*PHF1.v(t))
  }
  expected1 <- spline_params_group1$N*hybrid_integration(risk_integrand, 0, limit, disconts = NULL, pivot=restriction)
  observed1 <- spline_params_group1$N*hybrid_integration(observed_integrand, 0, limit, disconts = NULL, pivot=restriction)
  return(observed1 - expected1)
}



#' Log-rank test statistic (Chi-square) for two-sample logrank test based on spline parameters
#' @param spline_params_all spline parameters in the full population, containing knots, degrees of freedom, boundary knots and coefficients for both at-risk and event functions
#' @param spline_params_group1 spline parameters in one of the groups, containing knots, degrees of freedom, boundary knots and coefficients for both at-risk and event functions
#' @param spline_params_group2 spline parameters in the other group, containing knots, degrees of freedom, boundary knots and coefficients for both at-risk and event functions
#' @param restriction pivot for hybrid integration (should be >0 if steep hazard observed near 0)
#' @param limit integration limit for log-rank test to prevent error due to extrapolation, a good choice is the maximum value prior to which we have events in both groups
#' @return Log-rank test statistic (Chi-square with one DF under the null)
#' @export
log_rank_test <- function(spline_params_all, spline_params_group1, spline_params_group2, restriction=0, limit=spline_params_all$Max*0.99) {
  PHF.v <- function(t) {hazard_function_cloglog.v(t, spline_params_all$Event$Knots,
                                                  spline_params_all$Event$Degree, spline_params_all$Event$BKnots,
                                                  spline_params_all$Event$Coefficients)}
  PRF.v <- function(t) {survival_function_cloglog.v(t, spline_params_all$Risk$Knots,
                                                    spline_params_all$Risk$Degree, spline_params_all$Risk$BKnots,
                                                    spline_params_all$Risk$Coefficients)}
  PSF1.v <- function(t) {survival_function_cloglog.v(t, spline_params_group1$Event$Knots,
                                                     spline_params_group1$Event$Degree, spline_params_group1$Event$BKnots,
                                                     spline_params_group1$Event$Coefficients)}
  PHF1.v <- function(t) {hazard_function_cloglog.v(t, spline_params_group1$Event$Knots,
                                                   spline_params_group1$Event$Degree, spline_params_group1$Event$BKnots,
                                                   spline_params_group1$Event$Coefficients)}
  PRF1.v <- function(t) {survival_function_cloglog.v(t, spline_params_group1$Risk$Knots,
                                                     spline_params_group1$Risk$Degree, spline_params_group1$Risk$BKnots,
                                                     spline_params_group1$Risk$Coefficients)}
  PSF2.v <- function(t) {survival_function_cloglog.v(t, spline_params_group2$Event$Knots,
                                                     spline_params_group2$Event$Degree, spline_params_group2$Event$BKnots,
                                                     spline_params_group2$Event$Coefficients)}
  
  PHF2.v <- function(t) {hazard_function_cloglog.v(t, spline_params_group2$Event$Knots,
                                                   spline_params_group2$Event$Degree, spline_params_group2$Event$BKnots,
                                                   spline_params_group2$Event$Coefficients)}
  PRF2.v <- function(t) {survival_function_cloglog.v(t, spline_params_group2$Risk$Knots,
                                                     spline_params_group2$Risk$Degree, spline_params_group2$Risk$BKnots,
                                                     spline_params_group2$Risk$Coefficients)}
  risk_integrand <- function(t) {
    return(PRF1.v(t)*PHF.v(t))
  }
  observed_integrand <- function(t){
    return(PRF1.v(t)*PHF1.v(t))
  }
  expected1 <- spline_params_group1$N*hybrid_integration(risk_integrand, 0, limit, pivot=restriction)
  observed1 <- spline_params_group1$N*hybrid_integration(observed_integrand, 0, limit, pivot=restriction)
  risk_integrand <- function(t) {
    return(PRF2.v(t)*PHF.v(t))
  }
  observed_integrand <- function(t){
    return(PRF2.v(t)*PHF2.v(t))
  }
  expected2 <- spline_params_group2$N*hybrid_integration(risk_integrand, 0, limit, pivot=restriction)
  observed2 <- spline_params_group2$N*hybrid_integration(observed_integrand, 0, limit, pivot=restriction)
  LR_STAT <- (observed1-expected1)^2/expected1 + (observed2-expected2)^2/expected2
  return(list("X2"=LR_STAT, "O1" = observed1, "E1" = expected1, "O2"=observed2, "E2" = expected2))
}


#' Sum of squared KM influence function values for a set of observations (for variance estimation) for a given time
#' @param time vector of survival times
#' @param delta vector of censoring indicators
#' @param spline_params spline parameters, containing knots, degrees of freedom, boundary knots and coefficients for both at-risk and event functions
#' @param confint_time time we want to compute the influence function values for
#' @param restriction pivot for hybrid integration (should be >0 if steep hazard observed near 0)
#' @param weights vector of weights
#' @return sum of squared influence function values
#' @export
squared_sum_influence <- function(time, delta, spline_params, confint_time, restriction=0, weights=rep(1, length(time))) {
  PHF.v <- function(t) {hazard_function_cloglog.v(t, spline_params$Event$Knots,
                                                  spline_params$Event$Degree, spline_params$Event$BKnots,
                                                  spline_params$Event$Coefficients)}
  
  PSF.v <- function(t) {survival_function_cloglog.v(t, spline_params$Event$Knots,
                                                    spline_params$Event$Degree, spline_params$Event$BKnots,
                                                    spline_params$Event$Coefficients)}
  PRF.v <- function(t) {survival_function_cloglog.v(t, spline_params$Risk$Knots,
                                                    spline_params$Risk$Degree, spline_params$Risk$BKnots,
                                                    spline_params$Risk$Coefficients)}
  term1 <- delta*I(time <= confint_time) / PRF.v(time)
  time_trunc <- ifelse(time <= confint_time, time, confint_time)
  integrand <- function(u) return(PHF.v(u)/PRF.v(u))
  term2 <- hybrid_integration(integrand, 0, time_trunc, pivot=restriction)
  const <- weights/log(PSF.v(confint_time))
  influence_vec <- const*(term1-term2)
  return(sum(influence_vec^2))
}


#' Sum of squared logrank (observed-expected) influence function values for a set of observations (for variance estimation in weighted logrank test)
#' @param time vector of survival times
#' @param delta vector of censoring indicators
#' @param treatment binary vector for treatment assignment
#' @param spline_params spline parameters, containing knots, degrees of freedom, boundary knots and coefficients for both at-risk and event functions
#' @param spline_params_treated spline parameters for treatment group, containing knots, degrees of freedom, boundary knots and coefficients for both at-risk and event functions
#' @param restriction pivot for hybrid integration (should be >0 if steep hazard observed near 0)
#' @param weights vector of weights
#' @param p_treated proportion of population that is treated in the weighted sample
#' @param limit integration limit for log-rank test to prevent error due to extrapolation, a good choice is the maximum value prior to which we have events in both groups
#' @return sum of squared influence function values
#' @export
squared_sum_influence_logrank <- function(time, delta, treatment, spline_params, spline_params_treated, restriction=0, weights=rep(1, length(time)), p_treated=0.5, limit=spline_params$Max*0.99) {
  PHF.v <- function(t) {hazard_function_cloglog.v(t, spline_params$Event$Knots,
                                                  spline_params$Event$Degree, spline_params$Event$BKnots,
                                                  spline_params$Event$Coefficients)}
  PRF.v <- function(t) {survival_function_cloglog.v(t, spline_params$Risk$Knots,
                                                    spline_params$Risk$Degree, spline_params$Risk$BKnots,
                                                    spline_params$Risk$Coefficients)}
  PRFT.v <- function(t) {survival_function_cloglog.v(t, spline_params_treated$Risk$Knots,
                                                     spline_params_treated$Risk$Degree, spline_params_treated$Risk$BKnots,
                                                     spline_params_treated$Risk$Coefficients)}
  term1 <- I(time<limit)*delta*(treatment - p_treated*PRFT.v(time)/PRF.v(time))
  int_limit <- ifelse(time <= limit, time, limit)
  integrand <- function(u) {
    haz <- PHF.v(u)
    return(p_treated*haz*PRFT.v(u)/PRF.v(u))}
  term2 <- hybrid_integration(integrand, 0, int_limit, pivot=restriction)
  integrand <- function(u) return(PHF.v(u))
  term3 <- hybrid_integration(integrand, 0, int_limit, pivot=restriction)
  influence_vec <- weights*(term1+term2-term3*treatment)
  return(sum((influence_vec)^2))
}


#' Estimate risk and survival function splines at initial site, cloglog link
#' @param data dataframe containing censoring and survival time data
#' @param time name of survival time column (default is "time")
#' @param event name of censoring indicator column (default is "delta")
#' @param degree spline degrees of freedom
#' @param knots spline knot locations
#' @param max_time time to estimate spline to
#' @param bknots boundary knots (two values)
#' @param weights vector of weights
#' @param confint_times times to store squared sum of influence functions for variance estimation in weighted case
#' @param confint_restriction pivot for hybrid integration (should be >0 if steep hazard observed near 0) in weighted case where confidence intervals desired
#' @return spline parameters for risk and survival functions
#' @export
get_survival_splines_cloglog <- function(data, time="time", event="delta", degree=2, knots=5, max_time=max(data[,time]), bknots=c(0, max_time), weights=rep(1,nrow(data)),
                                         confint_times=NULL, confint_restriction=0) {
  
  n <- sum(weights)
  data[,time] <- ifelse(data[,time] > max_time, max_time, data[,time])
  Event_Curve <- survfit(Surv(data[,time], data[,event]) ~ 1, weights=weights)
  Event_Data <- data.frame(Time = Event_Curve$time, S_T = Event_Curve$surv)
  Event_Data <- filter(Event_Data, S_T != 1, Time > bknots[1], Time < bknots[2])
  max_ST <- max(Event_Data$S_T)
  min_ST <- min(Event_Data$S_T)
  Risk_Curve <- survfit(Surv(data[,time], rep(1,nrow(data))) ~ 1, weights=weights)
  Risk_Set_Data <- data.frame(Time = Risk_Curve$time, S_T = Risk_Curve$surv)
  Risk_Set_Data <- filter(Risk_Set_Data, S_T != 1, Time > bknots[1], Time < bknots[2])
  t_q <- get_time_quantiles(Event_Data$Time, Event_Data$S_T, seq(max_ST,min_ST,(-max_ST+min_ST)/(knots+1))[2:(knots+1)])
  times_event <- c(bknots[1], t_q, bknots[2])
  reg_points <- NULL
  for (i in seq(1, length(times_event)-1)){
    reg_points <- c(reg_points, seq(times_event[i], times_event[i+1], (times_event[i+1]-times_event[i])/10))
  }
  reg_data_event <- summary(Event_Curve, time = reg_points)
  Reg_Data <- data.frame(Time = reg_data_event$time, S_T = reg_data_event$surv)
  basis <- iSpline(Reg_Data$Time, knots=t_q, degree=degree, Boundary.knots = bknots)
  Spline_Event <- glm(S_T ~ basis, family=quasibinomial(link="cloglog"), data=Reg_Data)
  Event_Parameters <- list(t_q, degree, bknots, Spline_Event$coefficients)
  names(Event_Parameters) <- c("Knots", "Degree", "BKnots", "Coefficients")
  min_YT <- min(Risk_Set_Data$S_T)
  max_YT <- max(Risk_Set_Data$S_T)
  
  t_q_risk <- get_time_quantiles(Risk_Set_Data$Time, Risk_Set_Data$S_T, seq(max_YT,min_YT,(-max_YT+min_YT)/(knots+1))[2:(knots+1)])
  
  times_risk <- c(bknots[1], t_q_risk, bknots[2])
  reg_points <- NULL
  for (i in seq(1, length(times_risk)-1)){
    reg_points <- c(reg_points, seq(times_risk[i], times_risk[i+1], (times_risk[i+1]-times_risk[i])/10))
  }
  reg_data_risk <- summary(Risk_Curve, time = reg_points)
  Reg_Data <- data.frame(Time = reg_data_risk$time, Y_T = reg_data_risk$surv)
  basis <- iSpline(Reg_Data$Time, knots=t_q_risk, degree=degree, Boundary.knots = bknots)
  Spline_Risk <- glm(Y_T ~ basis, family=quasibinomial(link="cloglog"), data=Reg_Data)
  Risk_Parameters <- list(t_q_risk, degree, bknots, Spline_Risk$coefficients)
  names(Risk_Parameters) <- c("Knots", "Degree", "BKnots", "Coefficients")
  final_output <- list(Event_Parameters, Risk_Parameters, n, max_time)
  names(final_output) <- c("Event", "Risk", "N", "Max")
  influence_sums <- list()
  for (i in confint_times) {
    influence_sums[as.character(i)] <- squared_sum_influence(data[,time], data[,event], final_output, i, confint_restriction, weights)
  }
  final_output <- list(Event_Parameters, Risk_Parameters, n, max_time, influence_sums)
  names(final_output) <- c("Event", "Risk", "N", "Max", "Influence")
  return(final_output)
}


#' Update risk and survival function splines given survival observations, cloglog link
#' @param spline_params spline parameters, containing knots, degrees of freedom, boundary knots and coefficients for both at-risk and event functions
#' @param time vector of survival times
#' @param delta vector of censoring indicators
#' @param breakpoints number of breakpoints to use between each knot when adjusting spline
#' @param knots number of knots, can be larger than previous but should not increase by more than 1-2 at a given site
#' @param degree spline degrees of freedom
#' @param max_time time to estimate spline to
#' @param degree spline degrees of freedom
#' @param restriction pivot for hybrid integration (should be >0 if steep hazard observed near 0)
#' @return updated spline parameters for risk and survival functions
#' @export
update_survival_splines <- function(spline_params, time, delta,
                                    breakpoints=5, knots = length(spline_params$Event$Knots),
                                    degree = spline_params$Event$Degree, max_time=spline_params$Max, restriction=0) {
  time <- ifelse(time>max_time, max_time, time)
  n_to_add <- length(time)
  
  PSF.v <- function(t) {survival_function_cloglog.v(t, spline_params$Event$Knots,
                                                    spline_params$Event$Degree, spline_params$Event$BKnots,
                                                    spline_params$Event$Coefficients)}
  PRF.v <- function(t) {survival_function_cloglog.v(t, spline_params$Risk$Knots,
                                                    spline_params$Risk$Degree, spline_params$Risk$BKnots,
                                                    spline_params$Risk$Coefficients)}
  PHF.v <- function(t) {hazard_function_cloglog.v(t, spline_params$Event$Knots,
                                                  spline_params$Event$Degree, spline_params$Event$BKnots,
                                                  spline_params$Event$Coefficients)}
  PRHF.v <- function(t) {hazard_function_cloglog.v(t, spline_params$Risk$Knots,
                                                   spline_params$Risk$Degree, spline_params$Risk$BKnots,
                                                   spline_params$Risk$Coefficients)}
  # FIRST GET SPLINE FOR SURVIVAL FUNCTION
  
  t_q <- NULL
  all_knots <- c(spline_params$Event$BKnots[1], spline_params$Event$Knots, spline_params$Event$BKnots[2])
  for (i in seq(1,length(all_knots) - 1)){
    inc <-  (all_knots[i+1] - all_knots[i]) / (breakpoints + 1)
    t_q <- c(t_q, seq(all_knots[i] + inc, all_knots[i+1] - inc, inc))
  }
  current_values <- PSF.v(t_q)
  event_part <- NULL
  risk_t <- PRF.v(time)
  for (i in t_q) {
    delta_t <- delta*I(time <= i)
    event_part <- c(event_part, sum(delta_t/risk_t))
  }
  risk_set.v <- function(t) {
    return(colSums(outer(time, t, ">")))
  }
  integrand <- function(t) {
    hazard <- PHF.v(t)
    EY_T <- PRF.v(t)
    risk_t_new <- risk_set.v(t)
    return(hazard*risk_t_new/EY_T)
  }
  risk_part <- hybrid_integration(integrand, 0, t_q, disconts=c(unique(time)), pivot=restriction)
  divisor <- 0
  for (i in c(1:n_to_add)) {divisor <- divisor + 1/(i+spline_params$N)}
  divisor <- divisor / n_to_add
  new_values <- current_values - (current_values*divisor)*(event_part-risk_part)
  new_values <- ifelse(new_values > 1, 1, new_values)
  new_values <- ifelse(new_values < 0, 0, new_values)
  # need to transpose values when I don't have const hazard constraints on edges
  adjusted_data <- data.frame(Time = t_q, S_T = new_values)
  if (knots == length(spline_params$Event$Knots)) {t_knots <- spline_params$Event$Knots}
  else {
    full_prev <- c(spline_params$Event$BKnots[1], spline_params$Event$Knots, spline_params$Event$BKnots[2])
    inc <- (length(full_prev) - 1) / length(full_prev)
    placements <- seq(1+inc, length(full_prev)-inc, inc)
    t_knots <- NULL
    for (i in placements) {
      new_t <- full_prev[floor(i)]*(ceiling(i)-i) + full_prev[ceiling(i)]*(i-floor(i)) 
      t_knots <- c(t_knots, new_t)
    }
  }
  basis <- iSpline(adjusted_data$Time, knots=t_knots, degree=degree, Boundary.knots = spline_params$Event$BKnots)
  Spline_Event <- glm(S_T ~ basis, family=quasibinomial(link="cloglog"), data=adjusted_data)
  bknots <- spline_params$Event$BKnots
  Event_Parameters <- list(t_knots, degree, bknots, Spline_Event$coefficients)
  names(Event_Parameters) <- c("Knots", "Degree", "BKnots", "Coefficients")
  
  # THEN GET SPLINES FOR RISK FUNCTION
  t_q <- NULL
  all_knots <- c(spline_params$Risk$BKnots[1], spline_params$Risk$Knots, spline_params$Risk$BKnots[2])
  for (i in seq(1,length(all_knots) - 1)){
    inc <-  (all_knots[i+1] - all_knots[i]) / (breakpoints + 1)
    t_q <- c(t_q, seq(all_knots[i] + inc, all_knots[i+1] - inc, inc))
  }
  current_values <- PRF.v(t_q)
  event_part <- NULL
  risk_t <- PRF.v(time)
  for (i in t_q) {
    delta_t <- I(time <= i)
    event_part <- c(event_part, sum(delta_t/risk_t))
  }
  integrand <- function(t) {
    hazard <- PRHF.v(t)
    EY_T <- PRF.v(t)
    risk_t_new <- risk_set.v(t)
    return(hazard*risk_t_new/EY_T)
  }
  risk_part <- hybrid_integration(integrand, 0, t_q, disconts=c(unique(time)), pivot=restriction)
  new_values_risk <- current_values - (current_values*divisor)*(event_part-risk_part)
  new_values_risk <- ifelse(new_values_risk > 1, 1, new_values_risk)
  new_values_risk <- ifelse(new_values_risk < 0, 0, new_values_risk)
  adjusted_data <- data.frame(Time = t_q, Y_T = new_values_risk)
  
  if (knots == length(spline_params$Risk$Knots)) {t_knots <- spline_params$Risk$Knots}
  else {
    full_prev <- c(spline_params$Risk$BKnots[1], spline_params$Risk$Knots, spline_params$Risk$BKnots[2])
    inc <- (length(full_prev) - 1) / length(full_prev)
    placements <- seq(1+inc, length(full_prev)-inc, inc)
    t_knots <- NULL
    for (i in placements) {
      new_t <- full_prev[floor(i)]*(ceiling(i)-i) + full_prev[ceiling(i)]*(i-floor(i)) 
      t_knots <- c(t_knots, new_t)
    }
  }
  basis <- iSpline(adjusted_data$Time, knots=t_knots, degree=degree, Boundary.knots = spline_params$Risk$BKnots)
  adjusted_data$basis <- basis
  Spline_Risk <- glm(Y_T ~ basis, family=quasibinomial(link="cloglog"), data=adjusted_data)
  bknots <- spline_params$Risk$BKnots
  Risk_Parameters <- list(t_knots, degree, bknots, Spline_Risk$coefficients)
  names(Risk_Parameters) <- c("Knots", "Degree", "BKnots", "Coefficients")
  n <- length(time) + spline_params$N
  final_output <- list(Event_Parameters, Risk_Parameters, n, max_time)
  names(final_output) <- c("Event", "Risk", "N", "Max")
  return(final_output)
}

#' Update risk and survival function splines given survival observations with weights, cloglog link
#' @param spline_params spline parameters, containing knots, degrees of freedom, boundary knots and coefficients for both at-risk and event functions
#' @param time vector of survival times
#' @param delta vector of censoring indicators
#' @param weights vector of weights
#' @param breakpoints number of breakpoints to use between each knot when adjusting spline
#' @param knots number of knots, can be larger than previous but should not increase by more than 1-2 at a given site
#' @param degree spline degrees of freedom
#' @param max_time time to estimate spline to
#' @param degree spline degrees of freedom
#' @param restriction pivot for hybrid integration (should be >0 if steep hazard observed near 0)
#' @param confint_times times to store squared sum of influence functions for variance estimation in weighted case
#' @return updated spline parameters for risk and survival functions
#' @export
update_survival_splines_weighted <- function(spline_params, time, delta, weights=rep(1, length(time)),
                                             breakpoints=5, knots = length(spline_params$Event$Knots),
                                             degree = spline_params$Event$Degree, confint_times=names(spline_params$Influence), max_time=spline_params$Max, restriction=0) {
  time <- ifelse(time>max_time, max_time, time)
  n_to_add <- length(time)
  
  PSF.v <- function(t) {survival_function_cloglog.v(t, spline_params$Event$Knots,
                                                    spline_params$Event$Degree, spline_params$Event$BKnots,
                                                    spline_params$Event$Coefficients)}
  PRF.v <- function(t) {survival_function_cloglog.v(t, spline_params$Risk$Knots,
                                                    spline_params$Risk$Degree, spline_params$Risk$BKnots,
                                                    spline_params$Risk$Coefficients)}
  PHF.v <- function(t) {hazard_function_cloglog.v(t, spline_params$Event$Knots,
                                                  spline_params$Event$Degree, spline_params$Event$BKnots,
                                                  spline_params$Event$Coefficients)}
  PRHF.v <- function(t) {hazard_function_cloglog.v(t, spline_params$Risk$Knots,
                                                   spline_params$Risk$Degree, spline_params$Risk$BKnots,
                                                   spline_params$Risk$Coefficients)}
  # FIRST GET SPLINE FOR SURVIVAL FUNCTION
  t_q <- NULL
  all_knots <- c(spline_params$Event$BKnots[1], spline_params$Event$Knots, spline_params$Event$BKnots[2])
  for (i in seq(1,length(all_knots) - 1)){
    inc <-  (all_knots[i+1] - all_knots[i]) / (breakpoints + 1)
    t_q <- c(t_q, seq(all_knots[i] + inc, all_knots[i+1] - inc, inc))
  }
  current_values <- PSF.v(t_q)
  event_part <- NULL
  risk_t <- PRF.v(time)
  for (i in t_q) {
    delta_t <- weights*delta*I(time <= i)
    event_part <- c(event_part, sum(delta_t/risk_t))
  }
  risk_set.v <- function(t) {
    return(colSums(weights %*% outer(time, t, ">")))
  }
  integrand <- function(t) {
    hazard <- PHF.v(t)
    EY_T <- PRF.v(t)
    risk_t_new <- risk_set.v(t)
    return(hazard*risk_t_new/EY_T)
  }
  risk_part <- hybrid_integration(integrand, 0, t_q, disconts=c(unique(time)), pivot=restriction)
  divisor <- 0
  for (i in c(1:n_to_add)) {divisor <- divisor + 1/(sum(weights[1:i])+spline_params$N)}
  divisor <- divisor / n_to_add
  new_values <- current_values - (current_values*divisor)*(event_part-risk_part)
  new_values <- ifelse(new_values > 1, 1, new_values)
  new_values <- ifelse(new_values < 0, 0, new_values)
  # need to transpose values when I don't have const hazard constraints on edges
  adjusted_data <- data.frame(Time = t_q, S_T = new_values)
  if (knots == length(spline_params$Event$Knots)) {t_knots <- spline_params$Event$Knots}
  else {
    full_prev <- c(spline_params$Event$BKnots[1], spline_params$Event$Knots, spline_params$Event$BKnots[2])
    inc <- (length(full_prev) - 1) / length(full_prev)
    placements <- seq(1+inc, length(full_prev)-inc, inc)
    t_knots <- NULL
    for (i in placements) {
      new_t <- full_prev[floor(i)]*(ceiling(i)-i) + full_prev[ceiling(i)]*(i-floor(i)) 
      t_knots <- c(t_knots, new_t)
    }
  }
  basis <- iSpline(adjusted_data$Time, knots=t_knots, degree=degree, Boundary.knots = spline_params$Event$BKnots)
  Spline_Event <- glm(S_T ~ basis, family=quasibinomial(link="cloglog"), data=adjusted_data)
  bknots <- spline_params$Event$BKnots
  Event_Parameters <- list(t_knots, degree, bknots, Spline_Event$coefficients)
  names(Event_Parameters) <- c("Knots", "Degree", "BKnots", "Coefficients")
  
  # THEN GET SPLINES FOR RISK FUNCTION
  t_q <- NULL
  all_knots <- c(spline_params$Risk$BKnots[1], spline_params$Risk$Knots, spline_params$Risk$BKnots[2])
  for (i in seq(1,length(all_knots) - 1)){
    inc <-  (all_knots[i+1] - all_knots[i]) / (breakpoints + 1)
    t_q <- c(t_q, seq(all_knots[i] + inc, all_knots[i+1] - inc, inc))
  }
  current_values <- PRF.v(t_q)
  event_part <- NULL
  risk_t <- PRF.v(time)
  for (i in t_q) {
    delta_t <- weights*I(time <= i)
    event_part <- c(event_part, sum(delta_t/risk_t))
  }
  integrand <- function(t) {
    hazard <- PRHF.v(t)
    EY_T <- PRF.v(t)
    risk_t_new <- risk_set.v(t)
    return(hazard*risk_t_new/EY_T)
  }
  risk_part <- hybrid_integration(integrand, 0, t_q, disconts=c(unique(time)), pivot=restriction)
  new_values_risk <- current_values - (current_values*divisor)*(event_part-risk_part)
  new_values_risk <- ifelse(new_values_risk > 1, 1, new_values_risk)
  new_values_risk <- ifelse(new_values_risk < 0, 0, new_values_risk)
  adjusted_data <- data.frame(Time = t_q, Y_T = new_values_risk)
  if (knots == length(spline_params$Risk$Knots)) {t_knots <- spline_params$Risk$Knots}
  else {
    full_prev <- c(spline_params$Risk$BKnots[1], spline_params$Risk$Knots, spline_params$Risk$BKnots[2])
    inc <- (length(full_prev) - 1) / length(full_prev)
    placements <- seq(1+inc, length(full_prev)-inc, inc)
    t_knots <- NULL
    for (i in placements) {
      new_t <- full_prev[floor(i)]*(ceiling(i)-i) + full_prev[ceiling(i)]*(i-floor(i)) 
      t_knots <- c(t_knots, new_t)
    }
  }
  basis <- iSpline(adjusted_data$Time, knots=t_knots, degree=degree, Boundary.knots = spline_params$Risk$BKnots)
  adjusted_data$basis <- basis
  Spline_Risk <- glm(Y_T ~ basis, family=quasibinomial(link="cloglog"), data=adjusted_data)
  bknots <- spline_params$Risk$BKnots
  Risk_Parameters <- list(t_knots, degree, bknots, Spline_Risk$coefficients)
  names(Risk_Parameters) <- c("Knots", "Degree", "BKnots", "Coefficients")
  n <- sum(weights) + spline_params$N
  influence_sums <- spline_params$Influence
  confint_times <- as.numeric(confint_times)
  for (i in confint_times) {
    influence_sums[as.character(i)] <- as.numeric(influence_sums[as.character(i)]) + squared_sum_influence(time,delta, spline_params, i, restriction, weights)
  }
  final_output <- list(Event_Parameters, Risk_Parameters, n, max_time, influence_sums)
  names(final_output) <- c("Event", "Risk", "N", "Max", "Influence")
  return(final_output)
}

#' Update risk and survival function splines for a given site in chunks, cloglog link
#' @param spline_params spline parameters, containing knots, degrees of freedom, boundary knots and coefficients for both at-risk and event functions
#' @param time vector of survival times
#' @param delta vector of censoring indicators
#' @param chunksize number of observations to update at a time (would not recommend above 10)
#' @param breakpoints number of breakpoints to use between each knot when adjusting spline
#' @param knots_list number of knots for each chunk (should have length ceiling(length(time)/chunksize) and be non-decreasing)
#' @param degree spline degrees of freedom
#' @param max_time time to estimate spline to
#' @param loud whether to print out progress in number of observations passed
#' @param restriction pivot for hybrid integration (should be >0 if steep hazard observed near 0)
#' @return updated spline parameters for risk and survival functions
#' @export
update_survival_splines_chunk <- function(spline_params, time, delta, chunksize=5,
                                          breakpoints=5, knots_list = rep(length(spline_params$Event$Knots), ceiling(length(time)/chunksize)),
                                          degree = spline_params$Event$Degree, max_time=spline_params$Max, loud=TRUE, restriction=0) {
  updated_splines <- spline_params
  j <- 1
  for (i in seq(1, length(time), chunksize)) {
    if (loud) {
      print(i)
      print(i+chunksize-1)}
    if (i+chunksize <= length(time)) {
      chunk_time <- time[c(i:(i+chunksize-1))]
      chunk_delta <- delta[c(i:(i+chunksize-1))]
    }
    else{
      chunk_time <- time[c(i:length(time))]
      chunk_delta <- delta[c(i:length(time))]
    }
    updated_splines <- update_survival_splines(updated_splines, chunk_time, chunk_delta,
                                               breakpoints=breakpoints, knots = knots_list[j],
                                               degree = spline_params$Event$Degree, max_time=spline_params$Max, restriction=restriction)
    j <- j + 1
  }
  return(updated_splines)
}

#' Update risk and survival function splines for a given site in chunks with weighted data, cloglog link
#' @param spline_params spline parameters, containing knots, degrees of freedom, boundary knots and coefficients for both at-risk and event functions
#' @param time vector of survival times
#' @param delta vector of censoring indicators
#' @param chunksize number of observations to update at a time (would not recommend above 10)
#' @param weights vector of weights
#' @param breakpoints number of breakpoints to use between each knot when adjusting spline
#' @param knots_list number of knots for each chunk (should have length ceiling(length(time)/chunksize) and be non-decreasing)
#' @param degree spline degrees of freedom
#' @param confint_times times to store squared sum of influence functions for variance estimation in weighted case
#' @param max_time time to estimate spline to
#' @param loud whether to print out progress in number of observations passed
#' @param restriction pivot for hybrid integration (should be >0 if steep hazard observed near 0)
#' @return updated spline parameters for risk and survival functions
#' @export
update_survival_splines_chunk_weighted <- function(spline_params, time, delta, chunksize=5, weights=rep(1, length(time)),
                                                   breakpoints=5, knots_list = rep(length(spline_params$Event$Knots), ceiling(length(time)/chunksize)),
                                                   degree = spline_params$Event$Degree, confint_times=names(spline_params$Influence), max_time=spline_params$Max, loud=TRUE, restriction=0) {
  updated_splines <- spline_params
  j <- 1
  for (i in seq(1, length(time), chunksize)) {
    if (i+chunksize <= length(time)) {
      if (loud) {
        print(i)
        print(i+chunksize-1)}
      chunk_time <- time[c(i:(i+chunksize-1))]
      chunk_delta <- delta[c(i:(i+chunksize-1))]
      chunk_weights <- weights[c(i:(i+chunksize-1))]
    }
    else{
      if (loud) {
        print(i)
        print(length(time))}
      chunk_time <- time[c(i:length(time))]
      chunk_delta <- delta[c(i:length(time))]
      chunk_weights <- weights[c(i:length(time))]
    }
    updated_splines <- update_survival_splines_weighted(updated_splines, chunk_time, chunk_delta,
                                                        breakpoints=breakpoints, weights=chunk_weights, knots = knots_list[j],
                                                        degree = updated_splines$Event$Degree, confint_times = confint_times, max_time=updated_splines$Max, restriction=restriction)
    j <- j + 1
  }
  return(updated_splines)
}



#' Kaplan-Meier pseudo values derived from spline estimates and influence function
#' @param spline_params spline parameters, containing knots, degrees of freedom, boundary knots and coefficients for both at-risk and event functions
#' @param t_ps time of interest
#' @param time vector of survival times
#' @param delta vector of censoring indicators
#' @param restriction pivot for hybrid integration (should be >0 if steep hazard observed near 0)
#' @return pseudo values for each observation at time t_ps
#' @export
theoretical_KM_PS <- function(spline_params, t_ps, time, delta, restriction=0) {
  PSF.v <- function(t) {survival_function_cloglog.v(t, spline_params$Event$Knots,
                                                    spline_params$Event$Degree, spline_params$Event$BKnots,
                                                    spline_params$Event$Coefficients)}
  PRF.v <- function(t) {survival_function_cloglog.v(t, spline_params$Risk$Knots,
                                                    spline_params$Risk$Degree, spline_params$Risk$BKnots,
                                                    spline_params$Risk$Coefficients)}
  PHF.v <- function(t) {hazard_function_cloglog.v(t, spline_params$Event$Knots,
                                                  spline_params$Event$Degree, spline_params$Event$BKnots,
                                                  spline_params$Event$Coefficients)}
  S_T <- PSF.v(t_ps)
  EY_X <- PRF.v(time)
  term1 <- as.numeric(delta*I(time <= t_ps) / EY_X)
  upper_lims <- ifelse(t_ps < time, t_ps, time)
  lims_data <- data.frame(ID = seq(1, length(upper_lims)), upper_lims)
  lims_data <- lims_data %>% arrange(upper_lims)
  integrand <- function(t) {
    hazard <- PHF.v(t)
    EY_T <- PRF.v(t)
    return(hazard/EY_T)}
  term2 <- hybrid_integration(integrand, 0, lims_data$upper_lims, disconts=NULL, pivot=restriction)
  lims_data$term2 <- term2
  lims_data <- lims_data %>% arrange(ID)
  term2 <- lims_data$term2
  return(-c(S_T)*(term1 - term2) + c(S_T))
}
