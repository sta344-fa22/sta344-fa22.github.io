## Bayesian $R^2$

When we compute any statistic for our model we want to do so at each iteration so that we can obtain the posterior distribution of that particular statistic (e.g. the posterior distribution of $R^2$ in this case).

```{r}
df_R2 = df_pred %>%
  group_by(.iteration) %>%
  summarize(
    R2_classic = var(mu) / var(y),
    R2_bayes   = var(mu) / (var(mu) + var(resid))
  )

df_R2
```

## Uh oh ...

```{r}
summary(df_R2$R2_classic)
```

```{r}
summary(df_R2$R2_bayes)
```

## Visually

```{r}
df_R2 %>% tidyr::gather(method, R2, -.iteration) %>%
  ggplot(aes(x=R2, fill=method)) + 
  geom_density(alpha=0.5) +
  geom_vline(xintercept=modelr::rsquare(l, d), size=1)
```

## What if we collapsed first?

```{r echo=TRUE}
df_pred %>%
  group_by(i) %>%
  summarize(mu = mean(mu), y=mean(y), resid=mean(resid)) %>%
  summarize(
    R2_classic = var(mu) / var(y),
    R2_bayes   = var(mu) / (var(mu) + var(resid))
  )

modelr::rsquare(l, data=d)
```


## Some problems with $R^2$

Some new issues,

* $R^2$ doesn't really make sense in the Bayesian context

  * multiple possible definitions with different properties
  * fundamental equality doesn't hold anymore

. . .

Some old issues,

* $R^2$ always increases (or stays the same) when adding a predictor

* $R^2$ is highly susceptible to over fitting

* $R^2$ is sensitive to outliers

* $R^2$ depends heavily on current values of $Y$
  
  * $R^2$ can differ drastically for two equivalent models (i.e. nearly identical inferences about key parameters)


# Some Other Metrics

## Root Mean Square Error

The traditional definition of rmse  is as follows

$$ \text{RMSE} = \sqrt{ \frac{1}{n} \sum_{i=1}^n \left(Y_i - \hat{Y_i} \right)^2 } $$
  
  . . .

In the bayesian context, we have posterior samples from each parameter / prediction of interest so we can express this as

$$  \frac{1}{m} \sum_{s=1}^{m} \sqrt{  \frac{1}{n} \sum_{i=1}^n \left(Y_i - {\hat{Y}}^s_{i} \right)^2 } $$
  
  where $m$ is the number of iterations and $\hat{Y}^s_i$ is the prediction for $Y_i$ at iteration $s$.


## Continuous Rank Probability Score

Another approach is the continuous rank probability score which comes from the probabilistic forecasting literature, it compares the full posterior predictive distribution to the observation / truth.

$$ \text{CRPS} = \int_{-\infty}^\infty \left(F_{\hat{Y}}(z) - {1}_{z \geq Y}\right)^2 dz $$
  
  where $F_{\hat{Y}}$ is thes CDF of $\hat{Y}$ (the posterior predictive distribution for $Y$) and ${1}_{z \geq Y}$ is the indicator function which equals 1 when $z \geq Y$, the true/observed value of $Y$.

. . .

Since this calculates a score for a single probabilistic prediction we can naturally extend it to multiple predictions by calculating an average CRPS
$$
  \frac{1}{n} \sum_{i=1}^n \int_{-\infty}^\infty \left(F_{\hat{Y}_i}(z) - {1}_{z \geq Y_i}\right)^2 dz 
$$
  
  ## CDF vs Indicator
  
  ```{r echo=FALSE, warning=FALSE}
d_crps_ex = data_frame(
  value = rt(10000, df=2)
)

indicator = data.frame(value=seq(-12,12,len=1000)) %>% mutate(y = as.double(value >= 0))

ggplot(d_crps_ex, aes(x=value)) +
  geom_line(data=indicator, color="black", aes(y=y), size=1, alpha=0.5) +
  stat_ecdf(size=1, alpha=0.5, color="blue") +
  xlim(-10,10)
```


## Empirical CDF vs Indicator

```{r echo=FALSE, warning=FALSE}
d_crps_ex = data_frame(
  value = rt(10, df=2)
)

indicator = data.frame(value=seq(-12,12,len=1000)) %>% mutate(y = as.double(value >= 0))

ggplot(d_crps_ex, aes(x=value)) +
  geom_line(data=indicator, color="black", aes(y=y), size=1, alpha=0.5) +
  stat_ecdf(size=1, alpha=0.5, color="blue") +
  xlim(-5,5)
```


## Accuracy vs. Precision

```{r echo=FALSE, message=FALSE, warning=FALSE}
d_crps = data_frame(
  dist1 = rnorm(10000, sd=2)+0,
  dist2 = rnorm(10000, sd=2)+2,
  dist3 = rnorm(10000, sd=1)+0,
  dist4 = rnorm(10000, sd=1)+2
) %>% tidyr::gather(dist)

rmses = d_crps %>% group_by(dist) %>% summarise(rmse = (value-0)^2 %>% mean() %>% sqrt() %>% round(3))

rmse_lookup = rmses$rmse %>% setNames(rmses$dist)
rmse_labeler = function(variable, value)
  paste0(value, " (rmse = ", rmse_lookup[value],")")

g_up = ggplot(d_crps, aes(value, color=dist, fill=dist)) +
  geom_density(alpha=0.1) +
  facet_grid(~dist, labeller = rmse_labeler) + 
  geom_vline(xintercept=0)

crps = d_crps %>% group_by(dist) %>% summarise(crps = calc_crps(value, 0))

crps_lookup = crps$crps %>% setNames(crps$dist) %>% round(3)
crps_labeler = function(variable, value)
  paste0(value, " (crps = ", crps_lookup[value],")")

indicator = data.frame(value=seq(-10,10,len=1000)) %>% mutate(y = as.double(value >= 0))

g_low = ggplot(d_crps, aes(value, color=dist)) +
  geom_line(data=indicator, color="black", aes(y=y), size=1, alpha=0.5) +
  stat_ecdf(size=1, alpha=0.5) +
  facet_grid(~dist, labeller = crps_labeler)

gridExtra::grid.arrange(
  g_up,
  g_low
)

```


## Empirical Coverage

One final method, which assesses model calibration is to examine how well credible intervals, derived from the posterior predictive distributions of the $Y$s, capture the true/observed values.

```{r echo=FALSE}
set.seed(1111)
data.frame(
  x = rep(rnorm(20), 1000) + rnorm(20000),
  prediction = rep(1:20, 1000)
) %>%
  group_by(prediction) %>%
  ggplot(aes(y=prediction)) +
  tidybayes::stat_intervalh(aes(x = x), .width = c(.5, .90)) +
  geom_vline(xintercept = 0)

```

# Back to our example 

## RMSE

```{r}
rmse = df_pred %>%
  group_by(.iteration) %>%
  summarize(rmse = sqrt( sum( (y - y_pred)^2 ) / n())) %>%
  pull(rmse)

length(rmse)

head(rmse)

mean(rmse)

modelr::rmse(l, data = d)
```

## RMSE ($\mu$)

```{r}
rmse = df_pred %>%
  group_by(.iteration) %>%
  summarize(rmse = sqrt( sum( (y - mu)^2 ) / n())) %>%
  pull(rmse)

length(rmse)

head(rmse)

mean(rmse)

modelr::rmse(l, data = d)
```

## CRPS

```{r}
crps = df_pred %>% 
  group_by(i) %>% 
  summarise(crps = calc_crps(y_pred, y)) %>%
  pull(crps)

length(crps)

head(crps)

mean(crps)
```


## Empirical Coverage

\scriptoutput

```{r}
df_cover = df_pred %>% 
  group_by(x,y) %>% 
  tidybayes::mean_hdi(y_pred, .prob = c(0.5,0.9,0.95)) 

df_cover %>%
  mutate(contains = y >= .lower & y <= .upper) %>%
  group_by(prob=.width) %>%
  summarize(emp_cov = sum(contains)/n())
```

## Posterior predictive distribution ($y_{pred}$)

```{r echo=TRUE, fig.height=4}
df_pred %>% ungroup() %>%
  ggplot(aes(x=x)) + 
  tidybayes::stat_lineribbon(aes(y=y_pred), alpha=0.5) +
  geom_point(data=d, aes(y=y))
```


# Compared to what?

## Polynomial fit

```{r}
ggplot(d, aes(x=x,y=y)) + 
  geom_line() + 
  geom_smooth(
    method='lm', color="blue", se = FALSE, 
    formula = y~poly(x,5,simple=TRUE)
  )
```

## Model

```{r}
l_p = lm(y~poly(x,5,simple=TRUE), data=d)
summary(l_p)
```

## JAGS Model

```{r}
poly_model = 
  "model{
  # Likelihood
  for(i in 1:length(y)){
    y[i] ~ dnorm(mu[i], tau)
    y_pred[i] ~ dnorm(mu[i], tau)
    mu[i] = beta[1]        + beta[2]*x[i]   + beta[3]*x[i]^2 +
            beta[4]*x[i]^3 + beta[5]*x[i]^4 + beta[6]*x[i]^5
  }

  # Prior for beta
  for(j in 1:6){
    beta[j] ~ dnorm(0,1/1000)
  }

  # Prior for sigma / tau2
  tau ~ dgamma(1, 1)
  sigma2 = 1/tau
}"
```

```{r echo=FALSE, message=FALSE}
n_burn = 1000; n_iter = 5000

m = rjags::jags.model(
  textConnection(poly_model), data=d, 
  quiet=TRUE, n.chains = 1
) 
update(m, n.iter=n_burn, progress.bar="none")

df_poly = rjags::coda.samples(
  m, variable.names=c("beta","sigma2","mu","y_pred","y","x"), 
  n.iter=n_iter, progress.bar="none"
) %>%
  tidybayes::spread_draws(y_pred[i], y[i], x[i], mu[i]) %>%
  mutate(resid = y - mu)
```

## Posterior Predictive Distribution

```{r}
df_poly %>% ungroup() %>%
  ggplot(aes(x=x)) + 
  tidybayes::stat_lineribbon(aes(y=y_pred), alpha=0.5) +
  geom_point(data=d, aes(y=y))
```

## Comparing Results

```{r echo=FALSE}
rmse = rbind(
  mutate(df_pred, model="y~x"),
  mutate(df_poly, model="y~poly(x,5)")
) %>% 
  group_by(model,.iteration) %>%
  transmute(
    y_pred = (y-y_pred)^2  %>% mean() %>% sqrt(),
    mu     = (y-mu    )^2  %>% mean() %>% sqrt()
  ) %>% 
  tidyr::gather(parameter, rmse, -.iteration, -model) %>%
  group_by(parameter, model) %>%
  summarize(rmse = mean(rmse))
```

```{r echo=FALSE}
crps = rbind(
  mutate(df_pred, model="y~x"),
  mutate(df_poly, model="y~poly(x,5)")
) %>% 
  group_by(model,i) %>% 
  summarise(crps = calc_crps(y_pred, y)) %>%
  select(model,i,crps) %>%
  group_by(model) %>%
  summarize(crps = mean(crps)) %>%
  mutate(parameter = "y_pred")
```


```{r echo=FALSE}
empc = rbind(
  mutate(df_pred, model="y~x"),
  mutate(df_poly, model="y~poly(x,5)")
) %>% 
  group_by(x,y,model) %>% 
  tidybayes::mean_hdi(y_pred, .width = c(0.9)) %>%
  mutate(contains = y >= .lower & y <= .upper) %>%
  group_by(model) %>%
  summarize("emp_cov (90%)" = sum(contains)/n()) %>%
  mutate(parameter = "y_pred")
```

```{r echo=FALSE}
full_join(rmse, crps) %>% 
  full_join(empc) %>%
  tidyr::gather(metric, value, -parameter, -model) %>%
  tidyr::spread(model, value) %>%
  select(2,1,4,3) %>% 
  mutate(metric = factor(metric, levels = c("rmse", "crps", "emp_cov (90%)"))) %>%
  arrange(metric) %>%
  knitr::kable(digits = 3)
```