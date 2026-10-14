{smcl}
{* ! version 1.0.0 RCR - STEXPECTCR}{...}
{* October 14, 2026}{...}
{cmd:help stexpectcr}{right:Version 1.0.0 (October 14, 2026)}
{hline}

{title:Title}

{phang}
{cmd:stexpectcr} - Expected cumulative incidence in a reference population allowing for the competing risk of death{p_end}

{title:Syntax}

{phang}
{cmd:stexpectcr} [{it:newvar}] [{cmd:if}] [{cmd:in}]{cmd:,} {opt ratevar(varname)} {opt mratevar(varname)}
[{opt at(numlist)} {opt np:oints(#)} {opt met:hod(#)} {opt by(varlist)} {opt out:put(filename[, replace])} {opt nolist}]{p_end}

{title:Description}

{pstd}
{cmd:stexpect} gives expected {it:survival} from one hazard. {cmd:stexpectcr} takes two reference
rates, the population incidence rate of the event ({opt ratevar()}) and the population mortality
rate ({opt mratevar()}), and returns the expected cumulative incidence
F(t) = {it:integral of} S(u) * lambda(u) du with S(u) = exp(-integral of (lambda + mu)),
computed exactly for piecewise-constant rates. Data must be {cmd:stset} with an {cmd:id()}, split
(e.g. {cmd:stsplit}) so that both rates are constant within a record, and rates must be per unit of {cmd:_t}.{p_end}

{title:Options}

{phang}{opt method(1)} Ederer I: average of individual F(t); needs {opt at()} and records through max(at), no late entry.{p_end}
{phang}{opt method(2)} Ederer II (conditional): person-time-weighted mean hazards over the observed overlap with each grid interval; allows late entry (attained-age time scale).{p_end}
{phang}{opt method(3)} Hakulinen (default): potential-follow-up weighted; no late entry.{p_end}
{phang}{opt at(numlist)} or {opt npoints(#)} define the time grid (mutually exclusive). Evaluate on a fine grid.{p_end}
{phang}{opt by(varlist)} separate curves by group.{p_end}
{phang}{opt output(filename[, replace])} saves results.{p_end}
{phang}{opt nolist} suppresses the listing.{p_end}

{title:Output variables}

{pstd}{cmd:t_exp}, {cmd:atrisk}, the cumulative incidence ({it:newvar}, default {cmd:Cumexp}),
{cmd:Sfree} (expected event-free survival) and {cmd:Dexp} (expected death without the event); the last three sum to 1.{p_end}

{title:Example}

{phang}{cmd:. stset dox, fail(status==1) origin(dob) entry(survdate) scale(365.24) id(indexno)}{p_end}
{phang}{cmd:. stsplit ageband, at(0 1 5(5)85)}{p_end}
{phang}{cmd:. * merge incidence rate (rate_cancer) and mortality rate (mortrate) by sex, ageband, year}{p_end}
{phang}{cmd:. stexpectcr Cumexp, ratevar(rate_cancer) mratevar(mortrate) method(2) at(1(1)75)}{p_end}

{title:Author}

{pstd}Raoul Reulen{p_end}
