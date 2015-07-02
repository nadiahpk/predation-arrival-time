# Brief Description

Octave code to find the fitness of migratory birds arriving one day earlier or one day later than the prevailing arrival-day strategy in a population subject to Holling Type II predation.

# Reference

Harts, A.M.F., Kristense, N.P, Kokko, H. Paper in preparation.

# Quickstart

In Octave,

```> params```

will load a parameter dictionary p into the workspace containing the default parameter values. 

The function ```calc_fit``` will find the invasion fitness of birds arriving one day earlier or later than the prevailing strategy.

```
> prevailing_arrival_day_strategy = 2
> [fit_e,fit_l]=calc_fit(prevailing_arrival_day_strategy,p)
fit_e =  1.6434
fit_l =  0.71954
```

Here we see that the population is invasible (fitness > 1) to mutant strategies one day earlier than the prevailing strategy.
```
> [fit_e,fit_l]=calc_fit(1,p)
fit_e = NaN
fit_l =  0.73186
```
When the population has an arrival-day strategy of day 1 (the earliest day) it can not be invaded by later-arriving strategies.



