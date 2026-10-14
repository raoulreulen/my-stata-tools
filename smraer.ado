*! version 1.7 RCR - SMRAER
cap program drop smraer
program define smraer, sortpreserve
syntax varlist(min=3 max=3 numeric) [, dessmr(integer 1) desaer(integer 1) ///
                                       LABel(varname) pyrsper(integer 1)]

	tokenize `varlist'

	tempvar obs exp pyrs
	quietly {
		gen `obs'  = `1'
		gen `exp'  = `2'
		gen `pyrs' = `3'

		*-----------------------------------------------------------------------
		* SMR (exact Poisson confidence intervals)
		*-----------------------------------------------------------------------
		gen smr   =  `obs' / `exp'
		gen smrll =  invgammap(`obs',       0.025) / `exp'
		gen smrul =  invgammap(`obs' + 1,   0.975) / `exp'

		*-----------------------------------------------------------------------
		* AER (exact Poisson confidence intervals derived from obs CI)
		*
		* CI for observed events: [invgammap(O, 0.025), invgammap(O+1, 0.975)]
		* AER CI follows by subtracting E, dividing by raw person-years, and
		* scaling by pyrsper. Correct whether AER is positive OR negative.
		*-----------------------------------------------------------------------
		gen aer   =  (`obs'                          - `exp') / `pyrs' * `pyrsper'
		gen aerll =  (invgammap(`obs',       0.025)  - `exp') / `pyrs' * `pyrsper'
		gen aerul =  (invgammap(`obs' + 1,   0.975)  - `exp') / `pyrs' * `pyrsper'

		*-----------------------------------------------------------------------
		* String presentation — strtrim() removes leading spaces for compact output
		*-----------------------------------------------------------------------
		gen str smrstr = strtrim(string(smr,   "%9.`dessmr'f")) + " (" ///
		               + strtrim(string(smrll, "%9.`dessmr'f")) + ", "  ///
		               + strtrim(string(smrul, "%9.`dessmr'f")) + ")"

		gen str aerstr = strtrim(string(aer,   "%9.`desaer'f")) + " (" ///
		               + strtrim(string(aerll, "%9.`desaer'f")) + ", "  ///
		               + strtrim(string(aerul, "%9.`desaer'f")) + ")"

		*-----------------------------------------------------------------------
		* Ancillary string variables
		*-----------------------------------------------------------------------
		gen strobs = string(`obs', "%9.0f")
		gen strexp = string(`exp', "%9.1f")
		gen obsexp = strtrim(strobs) + "/" + strtrim(strexp)

		* n (%) relative to first row (assumed to be the 'all' category)
		gen str strobsperc = strtrim(string(`obs')) + " (" ///
		                   + strtrim(string((`obs' / `obs'[1]) * 100, "%9.1f")) + "%)"

		* crude observed rate per pyrsper person-years
		gen str obsrate = strtrim(string(`obs' / `pyrs' * `pyrsper', "%9.1f"))

		replace obsexp     = "" if obsexp == "./."
		replace smrstr     = "" if smr    == .
		replace aerstr     = "" if aer    == .
		replace strobsperc = "" if `obs'  == .
		replace obsrate    = "" if `obs'  == .
	}

	*---------------------------------------------------------------------------
	* Format pyrsper for display (e.g. 100000 → "100,000")
	*---------------------------------------------------------------------------
	local pyrslab : display %12.0gc `pyrsper'
	local pyrslab = strtrim("`pyrslab'")

	*---------------------------------------------------------------------------
	* Print results table to Stata results window
	*
	* Column layout (adjust c_* locals if you need wider/narrower columns):
	*   c_lbl  : start of label column         (left-aligned,  16 chars)
	*   c_perc : start of n (%) column         (left-aligned,  12 chars)
	*   c_e    : start of _e column            (right-aligned,  7 chars)
	*   c_rate : start of Rate* column         (right-aligned,  7 chars)
	*   c_smr  : start of SMR (95% CI) column  (left-aligned,  20 chars)
	*   c_aer  : start of AER (95% CI) column  (left-aligned,  20 chars)
	*---------------------------------------------------------------------------
	local c_lbl   1
	local c_perc 18
	local c_e    31
	local c_rate 39
	local c_smr  47
	local c_aer  68
	local twd    87

	di as text "{hline `twd'}"
	if "`label'" != "" {
		di as text %-16s "Cause"                        ///
		           _col(`c_perc') %-12s "n (%)"         ///
		           _col(`c_e')     %7s "_e"              ///
		           _col(`c_rate')  %7s "Rate*"           ///
		           _col(`c_smr')  %-20s "SMR (95% CI)"  ///
		           _col(`c_aer')  %-20s "AER (95% CI)"
	}
	else {
		di as text _col(`c_perc') %-12s "n (%)"         ///
		           _col(`c_e')     %7s "_e"              ///
		           _col(`c_rate')  %7s "Rate*"           ///
		           _col(`c_smr')  %-20s "SMR (95% CI)"  ///
		           _col(`c_aer')  %-20s "AER (95% CI)"
	}
	di as text "{hline `twd'}"

	forvalues i = 1/`=_N' {
		local sperc = strobsperc[`i']
		local sexp  = strtrim(strexp[`i'])
		local srate = obsrate[`i']
		local ssmr  = smrstr[`i']
		local saer  = aerstr[`i']

		if "`label'" != "" {
			local lbl = `label'[`i']
			di as result %-16s "`lbl'"                      ///
			             _col(`c_perc') %-12s "`sperc'"      ///
			             _col(`c_e')     %7s "`sexp'"        ///
			             _col(`c_rate')  %7s "`srate'"       ///
			             _col(`c_smr')  %-20s "`ssmr'"       ///
			             _col(`c_aer')  %-20s "`saer'"
		}
		else {
			di as result _col(`c_perc') %-12s "`sperc'"      ///
			             _col(`c_e')     %7s "`sexp'"        ///
			             _col(`c_rate')  %7s "`srate'"       ///
			             _col(`c_smr')  %-20s "`ssmr'"       ///
			             _col(`c_aer')  %-20s "`saer'"
		}
	}
	di as text "{hline `twd'}"
	di as text " * rate per `pyrslab' person-years"

end
