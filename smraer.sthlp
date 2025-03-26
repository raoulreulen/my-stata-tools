{smcl}
{* ! version 1.1 RCR - SMRAER}{...}
{* March 26, 2025}{...}
{cmd:help smraer}{right:Version 1.1 (March 26, 2025)}
{hline}

{title:Title}

{phang}
{cmd:smraer} - Calculate Standardized Mortality Ratios (SMRs) and Absolute Excess Risks (AERs) with confidence intervals{p_end}

{title:Syntax}

{phang}
{cmd:smraer} {it:varlist} [{cmd:,} {opt dessmr(#)} {opt desaer(#)}]{p_end}

{phang}
where {it:varlist} must contain exactly three numeric variables in this order:{p_end}
{pmore}
1. Observed events ({it:obs}){break}
2. Expected events ({it:exp}){break}
3. Person-years ({it:pyrs}){p_end}

{title:Options}

{phang}
{opt dessmr(#)} specifies the number of decimal places for SMR results in the output string. Default is 1.{p_end}

{phang}
{opt desaer(#)} specifies the number of decimal places for AER results in the output string. Default is 1.{p_end}

{title:Description}

{phang}
{cmd:smraer} computes Standardized Mortality Ratios (SMRs) and Absolute Excess Risks (AERs) based on observed events, expected events, and person-years. It generates:{p_end}
{pmore}
- SMR with 95% exact confidence intervals (using inverse gamma distribution).{break}
- AER with 95% exact confidence intervals (using inverse gamma distribution for absolute differences).{break}
- Formatted string variables ({cmd:smrstr}, {cmd:aerstr}) displaying results as "value (lower, upper)".{break}
- Additional variables: {cmd:strobs} (observed events as string), {cmd:strexp} (expected events as string), and {cmd:obsexp} (observed/expected as string).{p_end}

{phang}
Missing or invalid inputs result in empty strings in the output.{p_end}

{title:Variables Generated}

{phang}
{cmd:smrstr}   String of SMR and its 95% CI, e.g., "1.23 (0.89, 1.67)".{p_end}
{phang}
{cmd:aerstr}   String of AER and its 95% CI, e.g., "0.45 (0.12, 0.78)".{p_end}
{phang}
{cmd:strobs}   Observed events as a string, e.g., "15".{p_end}
{phang}
{cmd:strexp}   Expected events as a string with 1 decimal, e.g., "12.3".{p_end}
{phang}
{cmd:obsexp}   Combined observed/expected, e.g., "15/12.3".{p_end}

{title:Examples}

{phang}
Suppose you have data with observed deaths, expected deaths, and person-years:{p_end}
{pmore}
{cmd:. input obs exp pyrs}{break}
    obs         exp         pyrs{break}
    10          8.5         1000{break}
    5           3.2         500{break}
    {cmd:end}{p_end}

{phang}
Run the command with default decimals:{p_end}
{pmore}
{cmd:. smraer obs exp pyrs}{break}
{cmd:. list smrstr aerstr obsexp}{break}
     +-----------------------------------------+{break}
     |    smrstr      |    aerstr    | obsexp |{break}
     +-----------------------------------------+{break}
     | 1.2 (0.6,2.2) | 0.0 (0.0,0.0)| 10/8.5 |{break}
     | 1.6 (0.5,3.6) | 0.0 (0.0,0.0)| 5/3.2  |{break}
     +-----------------------------------------+{p_end}

{phang}
With custom decimals:{p_end}
{pmore}
{cmd:. smraer obs exp pyrs, dessmr(2) desaer(3)}{break}
{cmd:. list smrstr aerstr}{break}
     +----------------------------+{break}
     |     smrstr    |   aerstr  |{break}
     +----------------------------+{break}
     | 1.18 (0.57,2.17)| 0.002 (0.000,0.004)|{break}
     | 1.56 (0.51,3.64)| 0.004 (0.000,0.008)|{break}
     +----------------------------+{p_end}

{title:Author}

{phang}
Raoul C Reulen, r.c.reulen@bham.ac.uk {p_end}

{title:Also See}

{phang}
Online: {cmd:net install smraer, from([your-URL])}{p_end}