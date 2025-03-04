#------------------------------------------------
#' Function that returns values from deterministic model based on user-specified series of EIR values
#'
#' \code{data_gen_piecewise} Generate simulated data based on a user-specified time series of EIR values
#'
#' @param EIR_vals Initial EIR value to kick of the random walk
#' @param EIR_times Minimum value over which the random walk will not exceed.
#' @param sample_size Average sample size for number of tested individuals in the simulated data set.
#' @param sample_sd Standard deviation of the normal distribution of
#'                   sample size for the number of tested individuals in the simulated data set.
#' @param out_step The number of days between each simulated data point.
#' @param duration The length of time in days the model will run.
#' @param prop_treated Proportion of the cases that are treated
#' @param init_age Age groups within the model
#' @param plot_instance If TRUE, plots the simulated prevalence.
#'
#' @export
data_gen_piecewise <- function(EIR_vals,
                               EIR_times,
                               sample_size=220,
                               sample_sd=40,
                               out_step=30,
                               duration= 5*365,
                               prop_treated = 0.4,
                               init_age = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80),
                               plot_instance = TRUE){
  model_file <- system.file("odin", "odin_model_stripped_matched.R", package = "mamasante")

  rA_preg <- 0.00512821
  rU_preg <- 0.00906627
  het_brackets <- 5
  comparison <- 'u5'
  init_EIR <- EIR_vals[1]

  if(length(EIR_vals)!=length(EIR_times)){
    errorCondition('Length of EIR_vals and EIR_times must match.')

  }
  mpl <- mamasante::model_param_list_create(init_EIR = init_EIR,
                                            init_ft = prop_treated,
                                            EIR_times=EIR_times,
                                            EIR_vals=EIR_vals,
                                            comparison=comparison,
                                            lag_rates=10
  )

  pars <- mamasante::equilibrium_init_create_stripped(age_vector = init_age,
                                                      init_EIR = init_EIR,
                                                      ft = prop_treated,
                                                      model_param_list = mpl,
                                                      het_brackets = het_brackets)

  ##The malaria model but only on human side (i.e. no mosquitoes to worry about)
  generator <- odin::odin(model_file)
  state_use <- pars[names(pars) %in% coef(generator)$name]

  # create model with initial values
  mod <- generator(user = state_use, use_dde = TRUE)
  tt <- seq(0, duration, out_step)

  # run the simulation to base the data
  start.time <- Sys.time()
  mod_run <- mod$run(tt, step_max_n = 1e7,
                     atol = 1e-5,
                     rtol = 1e-5)
  print(Sys.time()-start.time)

  # shape output
  out <- mod$transform_variables(mod_run)

  # plot data and generate data
  if(plot_instance){
    plot(out$t,out$prev,col="white")
    lines(out$t,out$prev,col="blue",lwd=4)
  }
  tested<-round(rnorm(length(out$prev_all),sample_size,sample_sd))
  positive<-rbinom(length(out$prev_all),tested,out$prev_all)
  month <- seq.Date(from = as.Date('2015-01-01'),by = 'month',length.out = length(tt))
  data_raw<-data.frame(t=out$t+30,
                       tested=tested,
                       positive=positive,
                       prev_true=out$prev_all,
                       EIR_true=out$EIR_out,
                       inc_true=out$incunder5)
  return(data_raw)
}
