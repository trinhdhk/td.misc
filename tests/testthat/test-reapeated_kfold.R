full_data=list(
  N = nrow(iris),
  nX = 4,
  X = iris[,1:4],
  Y = as.integer(iris[,5])
)
cv_data = repeated_kfold(full_data, N_obs = nrow(iris), K = 2, N_rep = 2)
tmpfile = fs::file_temp(ext = '.stan')
stan_code = "
    data {
      int<lower=0> N;
      int<lower=0> nX;
      matrix[N, nX] X;
      int<lower=1, upper=3> Y[N];
      int<lower=0, upper=1> keptin[N];
    }

    parameters {
      matrix[nX, 3] a;
    }

    model {
      to_vector(a) ~ normal(0,5);
      for (n in 1:N) {
        if (keptin[n] == 1) {
          row_vector[3] C = X[n,:] * a;
          Y[n] ~ categorical_logit(C');
        }
      }
    }

    generated quantities {
      vector[N] log_lik;
      for (n in 1:N){
        row_vector[3] C = X[n,:] * a;
        log_lik[n] = categorical_logit_lpmf(Y[n] | C');
      }
    }
  "
writeLines(stan_code, tmpfile)
sampler = cmdstanr::cmdstan_model(tmpfile)
fits = stan_kfold(sampler = sampler, list_of_datas = cv_data$inputs,
                  cores = 15, chains = 2, seed = 122,
                  merge = T, adapt_delta=.65)
test_that("repeated kfold split works", {
  testthat::expect_identical(class(fits), 'array')
})
