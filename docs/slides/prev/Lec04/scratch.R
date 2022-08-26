

# Gradient Boosting model

## Model

```{r}
y = anguilla$presence %>% as.integer()
x = model.matrix(presence~.-1, data=anguilla)

xg = xgboost::xgboost(data=x, label=y, nthead=4, nround=35, objective="binary:logistic")
```

## Residuals

```{r echo=FALSE}
d_xg = anguilla %>%
  mutate(p_hat = predict(xg, newdata=x)) %>%
  mutate(p_hat_bin = p_hat - (p_hat %% 0.025)) %>%
  mutate(
    resid = presence - p_hat,
    pearson = (presence - p_hat)/(p_hat*(1-p_hat)),
    deviance = sign(presence-p_hat) * sqrt(-2 * (presence * log(p_hat) + (1-presence) * log(1-p_hat)))
  )

resid_xg = d_xg %>%
  select(p_hat:deviance) %>%
  tidyr::gather(type, value, -p_hat, -p_hat_bin) %>%
  mutate(type = factor(type, levels=c("resid","pearson","deviance"))) 

ggplot(resid_xg, aes(x=p_hat, y=value, color=type)) +
  geom_point(alpha=0.1) +
  facet_wrap(~type, ncol=3, scale="free_y")

resid_xg %>%
  group_by(type, p_hat_bin) %>%
  summarize(mean = mean(value)) %>%
  ggplot(aes(x=p_hat_bin, y=mean, color=type)) +
  geom_point() +
  facet_wrap(~type, ncol=3, scale="free_y")
```


