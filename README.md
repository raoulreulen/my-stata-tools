# my-stata-tools

Stata commands by Raoul Reulen.

| Command | Purpose |
|---|---|
| `smraer` | SMRs and absolute excess risks with exact confidence intervals |
| `stexpectcr` | Expected cumulative incidence allowing for the competing risk of death (Ederer I, Ederer II, Hakulinen); companion to `stexpect` |

## Install

```stata
net install stexpectcr, from(https://raw.githubusercontent.com/raoulreulen/my-stata-tools/main) replace
net install smraer,     from(https://raw.githubusercontent.com/raoulreulen/my-stata-tools/main) replace
help stexpectcr
```

## Using stexpectcr

`stexpect` returns expected *survival* from one hazard. `stexpectcr` takes two reference rates, the
general-population incidence rate of the event and the all-cause mortality rate, and returns the expected
cumulative incidence F(t) = integral of S(u) * lambda(u) du, with S(u) = exp(-integral of (lambda + mu)).

You supply the reference rates yourself: they must be variables in your data, merged by sex, age band and
calendar period onto records that you have split with `stsplit`, so that both rates are constant within a record.
Rates are per person per unit of `_t`. No rate files are distributed with the package.

```stata
stset dox, fail(status==1) origin(dob) entry(survdate) scale(365.24) id(id)
stsplit ageband, at(0 1 5(5)85)
stsplit yeargrp, after(time=d(1/1/1900)) at(40(1)123)
replace yeargrp = yeargrp + 1900
* merge incidence rate  -> rate_cancer   (by sex ageband yeargrp)
* merge mortality rate  -> mortrate      (by sex ageband yeargrp)
stexpectcr Cumexp, ratevar(rate_cancer) mratevar(mortrate) method(2) at(1(1)75)
```

Cumulative incidence is returned as a proportion in `Cumexp`, with `Sfree` (expected event-free survival) and
`Dexp` (expected death without the event). `method(2)` allows late entry (attained-age time scale);
`method(1)` and `method(3)` require records contiguous from time 0.
Sources of general-population mortality rates include the Human Mortality Database and national statistics offices.
