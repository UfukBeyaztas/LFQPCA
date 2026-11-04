lfqpcr_simulate <- function(case = paste0("Case", 1:6), I_train = 60,
                            I_test = 20, J = 5, tau_F = 0.5,
                            tau_R = 0.9, replication = 1, seed = NULL,
                            settings = lfqpca_sim_settings()) {
  case <- match.arg(case)
  if (is.null(seed)) seed <- settings$seed_experiment2 + replication - 1
  lfqpcr_sim_data(case, seed, I_train = I_train, I_test = I_test, J = J,
                  tau_F = tau_F, tau_R = tau_R, settings = settings)
}