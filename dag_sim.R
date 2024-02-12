library(MASS)
library(tidyverse)
library(latex2exp)

theme_set(theme_bw())

rmse <- function(yhat, y) {
  sqrt(mean((yhat-y)^2))
}

gen_data <- function(alpha, beta, gamma, delta, n=1e7) {
  
  s_B = sqrt(1 - beta^2)
  cov_B1_B2 = delta
  Sigma_B = matrix(c(s_B^2, cov_B1_B2, cov_B1_B2, s_B^2), nrow=2, ncol=2)
  U_B = mvrnorm(n, c(0,0), Sigma_B)
  
  s_A = sqrt(1 - alpha^2 - gamma^2 - 2*alpha*beta*gamma)

  d <- tibble(
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    delta = delta, 
    
    U_Z = rnorm(n),
    U_B0 = U_B[,1],
    U_B1 = U_B[,2],
    U_A0 = rnorm(n, sd = s_A),
    U_A1 = rnorm(n, sd = s_A),
    
    Z = U_Z,
    B0 = beta * Z + U_B0,
    B1 = beta * Z + U_B1,
    A0 = alpha * Z + gamma * B0 + U_A0,
    A1 = alpha * Z + gamma * B1 + U_A1
  )
  
  return(d)
}

gen_data(alpha = 0.1, beta = 0.2, gamma = 0.3, delta = 0.5) %>% 
  summarize(
    across(c(Z, B0, B1, A0, A1), .fns=list(s=sd)),
    cov_U_B = cov(U_B0, U_B1),
    cov_A0_Z = cov(A0, Z),
    cov_A1_Z = cov(A1, Z),
    cov_A_Z_ans = alpha + beta*gamma,
    cov_B1_Z = cov(B1, Z),
    cov_B1_Z_ans = beta,
    cov_A1_A0 = cov(A1, A0),
    cov_A1_A0_ans = alpha^2 + gamma^2*delta + 2*alpha*beta*gamma + beta^2*gamma^2,
    cov_B1_A0 = cov(B1, A0),
    cov_B1_A0_ans = gamma*delta + beta^2*gamma + alpha*beta
  ) %>%
  print()

perf <- function(beta, alpha, gamma, delta, n=1e4) {
  d <- gen_data(alpha, beta, gamma, delta, n)
  simple_m <- lm(A1 ~ A0, data=d)
  complex_m <- lm(A1 ~ A0 + Z, data=d)
  
  d <- d %>%
    mutate(
      simple_p = predict(simple_m),
      complex_p = predict(complex_m),
    )
  
  d %>%
    summarize(
      beta = beta[1],
      label = c('True label', 'True label', 'Proxy label', 'Proxy label'),
      Model = c('Simple', 'Complex', 'Simple', 'Complex'),
      RMSE = c(rmse(simple_p, B1), rmse(complex_p, B1), 
               rmse(simple_p, A1), rmse(complex_p, A1)),
  )
}

set.seed(1)
plot_data <- rep(seq(0, 0.6, .05), 10) %>%
  map_dfr(.f = ~ perf(alpha=0.4, beta=.x, gamma=0.4, delta=0.4, n=1e4)) %>%
  mutate(Model = fct_relevel(Model, 'Simple', 'Complex'))

ggplot(plot_data, aes(x=beta, y=RMSE)) + 
  geom_smooth(aes(color=Model), se=FALSE) + 
  scale_x_continuous(TeX('\\beta')) + 
  facet_grid(.~label)

ggsave('figs/dag_sim.pdf', height=4, width=9)

