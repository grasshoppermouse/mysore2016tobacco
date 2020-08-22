
library(mixtools)

mean1 <- m_mix$mu[1]
sd1 <- m_mix$sigma[1]
mean2 <- m_mix$mu[2]
sd2 <- m_mix$sigma[2]

thedata <- function(n, params){
  b0 <- params['b0']
  b1 <- params['b1']
  b2 <- params['b2']
  b3 <- params['b3']
  sigma <- params['sigma']
  d <- tibble(
    baseline_cotinine = 10^(c(rnorm(n/2, mean1, sd1), rnorm(n/2, mean2, sd2))),
    condition = rbinom(n, 1, p = 0.5),
    followup_cotinine = b0 + b1*baseline_cotinine + b2*condition + b3 * baseline_cotinine * condition + rnorm(n, 0, sigma)
  )
}

pvalue <- function(n, params){
  d <- thedata(n, params)
  m <- summary(lm(followup_cotinine ~ baseline_cotinine*condition, d))
  return(m$coefficients[4,4]) # Interaction term
}

nsim <- 2000
N <- seq(40, 100, 6)
names(N) <- N

m_pwr <- lm(post ~ pre * Presentation, data = df)
summary(m_pwr)

params1 <- c(b0=-11.8, b1=2.17, b2=53.9, b3=-1.53, sigma=169)
params2 <- c(b0=0, b1=1, b2=0, b3=-0.75, sigma=169)
params3 <- c(b0=0, b1=1, b2=0, b3=-0.5, sigma=169)

pwr1 <- map_dbl(N, ~sum(replicate(nsim, pvalue(., params1))<0.05)/nsim)
pwr1

pwr2 <- map_dbl(N, ~sum(replicate(nsim, pvalue(., params2))<0.05)/nsim)
pwr2

pwr3 <- map_dbl(N, ~sum(replicate(nsim, pvalue(., params3))<0.05)/nsim)
pwr3

df_pwr <- tibble(
  Power = unname(c(pwr1, pwr2, pwr3)),
  Sample_size = rep(N, 3),
  Parameters = c(rep('Estimated from data', length(N)), rep('Theoretically important', length(N)), rep('Theoretically small', length(N)))
)

power_curves <-
  ggplot(df_pwr, aes(Sample_size, Power, colour=Parameters)) + 
  geom_line() + 
  geom_vline(xintercept = nrow(samples)/2, linetype = 'dotted') +
  annotate(geom = 'text', label = 'Study sample size', x = 67, y = 0.25, hjust=0) +
  labs(x = '\nSample size', y = 'Power\n') +
  theme_bw(15)
power_curves
