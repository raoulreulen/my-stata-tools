*! stexpectcr 1.0.0  14OCT2026
*! Expected cumulative incidence of an event (e.g. cancer) in a reference
*! population, allowing for the competing risk of death.
*! Companion to stexpect 1.1.2 (Coviello, 2005), which gives expected
*! SURVIVAL from a single hazard.  Here two reference rates are used:
*!    ratevar()  = population incidence rate of the event   (lambda)
*!    mratevar() = population mortality rate                (mu)
*! and for each individual
*!    F_i(t) = int_0^t S_i(u) lambda_i(u) du,  S_i(u) = exp(-int_0^u (lambda_i+mu_i))
*! is computed exactly for piecewise-constant rates.
*!
*! Data requirements (same logic as stexpect):
*!   - stset with an id() variable; records split (stsplit) so that ratevar and
*!     mratevar are constant within each record; rates are per unit of _t.
*!   - methods 1 and 3: records within id contiguous from time 0 (no late entry, no gaps);
*!     method 2 (Ederer II) allows late entry, e.g. attained-age time scale.
*!   - method(1) Ederer I: every individual must have records through max(at()).
*!   - method(3) Hakulinen (default): individuals contribute while their
*!     POTENTIAL follow-up (end of last record) has not been reached.
*!
*! Ederer II (method(2)) uses person-time-weighted mean hazards over the observed
*! overlap of each record with each grid interval (allows late entry / attained
*! age time scale):  h_c = sum(lambda*ov)/sum(ov),  h_m likewise, a = h_c+h_m,
*!    S_E(t_k) = S_E(t_k-1) e^(-a dt) ;  F_E += S_E(t_k-1) h_c/a (1-e^(-a dt)).
*!
*! Hakulinen estimator on the time grid t_1<...<t_K, with R_k = {i: C_i >= t_k}:
*!    p_k = sum_R S_i(t_k) / sum_R S_i(t_k-1)
*!    q_k = sum_R [F_i(t_k)-F_i(t_k-1)] / sum_R S_i(t_k-1)
*!    S_E(t_k) = S_E(t_k-1) p_k ;  F_E(t_k) = F_E(t_k-1) + S_E(t_k-1) q_k
*! (the conditional-ratio structure of stexpect method 3, extended to the
*! cause-specific incidence).  Use a fine grid; it is evaluated on the grid only.
program define stexpectcr
    version 14
    st_is 2 analysis
    syntax [ newvarname ] [if] [in], RATEVAR(varname numeric) MRATEVAR(varname numeric) ///
        [ AT(numlist sort >0) NPoints(integer 0) METhod(integer 3) BY(varlist min=1 max=5) ///
          OUTput(string asis) noLIst ]

    marksample touse, novarlist
    markout `touse' `by' `ratevar' `mratevar', strok
    qui replace `touse' = 0 if _st==0 | _st>=.

    if `"`_dta[st_id]'"' == "" {
        di as err "stexpectcr requires that you have previously stset an id() variable"
        exit 198
    }
    cap confirm numeric variable `_dta[st_id]'
    if _rc {
        di as err "stexpectcr requires that the id() variable is numeric"
        exit 198
    }
    if `method'<1 | `method'>3 {
        di as err "method(1) Ederer I, method(2) Ederer II (conditional) or method(3) Hakulinen"
        exit 198
    }
    if "`at'"=="" & `npoints'<2 {
        di as err "specify at(numlist) or npoints(#) with # >= 2"
        exit 198
    }
    if "`at'"!="" & `npoints'>0 {
        di as err "at() and npoints() are mutually exclusive"
        exit 198
    }
    if `method'==1 & "`at'"=="" {
        di as err "Ederer I method requires an at(numlist) option"
        exit 198
    }
    if "`at'"!="" local at : list uniq at

    // output(filename[, replace]) -- as in stexpect/strate
    if `"`output'"' != "" {
        gettoken outnam output : output, parse(",")
        gettoken comma output : output, parse(",")
        if `"`comma'"' == "," {
            gettoken outrepl output : output, parse(" ,")
            gettoken comma output : output, parse(" ,")
            if `"`outrepl'"' != "replace" | `"`comma'"' != "" {
                di as err "option output() invalid"
                exit 198
            }
        }
        else if `"`comma'"' != "" {
            di as err "option output() invalid"
            exit 198
        }
        else confirm new file `"`outnam'.dta"'
    }

    local cif "`varlist'"
    if "`cif'"=="" local cif Cumexp
    local id `_dta[st_id]'

    preserve
    qui keep if `touse'
    qui count
    if r(N)==0 error 2000
    qui count if `ratevar'<0 | `mratevar'<0
    if r(N) {
        di as err "negative rates found in `ratevar' / `mratevar'"
        exit 459
    }

    qui keep `id' `ratevar' `mratevar' `by' _t0 _t
    tempvar bad a dt H S0 dF Fc F0 Cend
    if `method'!=2 {
        qui bysort `id' (_t): gen byte `bad' = cond(_n==1, _t0!=0, _t0!=_t[_n-1])
        qui count if `bad'
        if r(N) {
            di as err "records within id are not contiguous from time 0 " ///
                "(gaps or late entry, `r(N)' records); not supported"
            exit 459
        }
        qui {
            gen double `a'  = `ratevar' + `mratevar'
            gen double `dt' = _t - _t0
            bysort `id' (_t): gen double `H'  = sum(`a'*`dt')
            gen double `S0' = exp(-(`H' - `a'*`dt'))
            gen double `dF' = cond(`a'>0, `S0'*`ratevar'/`a'*(1-exp(-`a'*`dt')), 0)
            bysort `id' (_t): gen double `Fc' = sum(`dF')
            gen double `F0' = `Fc' - `dF'
            bysort `id' (_t): gen double `Cend' = _t[_N]
        }
    }

    // time grid
    qui su _t, meanonly
    local tmax = r(max)
    if "`at'"=="" {
        forvalues i = 1/`npoints' {
            local grid `grid' `=`tmax'*`i'/`npoints''
        }
    }
    else local grid `at'

    if `method'==1 {
        qui su `Cend', meanonly
        local tg : word count `grid'
        local gmax : word `tg' of `grid'
        if r(min) < `gmax' - 1e-9 {
            di as err "Ederer I: some individuals have follow-up records only to " r(min) ///
                " < max(at())=`gmax'. Extend the records (rates) to the common horizon."
            exit 459
        }
    }

    tempfile recs acc
    qui save `recs'

    local k 0
    local prev 0
    foreach t of local grid {
        local ++k
        if `method'==2 {
            // Ederer II: person-time-weighted mean hazards over the OBSERVED
            // overlap of each record with (prev, t]; late entry allowed.
            qui {
                use `recs', clear
                tempvar ov ec em u pos
                gen double `ov' = max(0, min(_t, `t') - max(_t0, `prev'))
                gen byte `pos' = `ov' > 1e-12
                replace `ov' = 0 if !`pos'
                gen double `ec' = `ratevar'*`ov'
                gen double `em' = `mratevar'*`ov'
                bysort `id' `pos': gen byte `u' = `pos' & _n==1
                collapse (sum) atrisk=`u' Ec=`ec' Em=`em' PY=`ov', by(`by')
                gen int __k = `k'
                gen double t_exp = `t'
                gen double __dt = `t' - `prev'
                if `k'>1 append using `acc'
                save `acc', replace
            }
            local prev `t'
            continue
        }
        qui {
            use `recs', clear
            tempvar inR cov cp Sc Fcu Sp Fp n
            gen byte `inR' = `Cend' >= `t' - 1e-9
            gen byte `cov' = (_t0 < `t' - 1e-9) & (_t >= `t' - 1e-9)
            gen double `Sc'  = cond(`cov' & `inR', `S0'*exp(-`a'*(`t'-_t0)), 0)
            gen double `Fcu' = cond(`cov' & `inR', `F0' + cond(`a'>0, ///
                `S0'*`ratevar'/`a'*(1-exp(-`a'*(`t'-_t0))), 0), 0)
            if `prev'==0 {
                gen double `Sp' = cond(`cov' & `inR', 1, 0)
                gen double `Fp' = 0
            }
            else {
                gen byte `cp' = (_t0 < `prev' - 1e-9) & (_t >= `prev' - 1e-9)
                gen double `Sp' = cond(`cp' & `inR', `S0'*exp(-`a'*(`prev'-_t0)), 0)
                gen double `Fp' = cond(`cp' & `inR', `F0' + cond(`a'>0, ///
                    `S0'*`ratevar'/`a'*(1-exp(-`a'*(`prev'-_t0))), 0), 0)
            }
            gen double `n' = `cov' & `inR'
            collapse (sum) atrisk=`n' Scur=`Sc' Sprev=`Sp' Fcur=`Fcu' Fprev=`Fp', by(`by')
            gen int __k = `k'
            gen double t_exp = `t'
            if `k'>1 append using `acc'
            save `acc', replace
        }
        local prev `t'
    }

    qui {
        use `acc', clear
        local gby `by'
        if "`by'"=="" {
            gen byte __g = 1
            local gby __g
        }
        sort `gby' __k
        if `method'==2 {
            gen double __hc = cond(PY>0, Ec/PY, 0)
            gen double __hm = cond(PY>0, Em/PY, 0)
            gen double __a  = __hc + __hm
            gen double __p  = exp(-__a*__dt)
            gen double __q  = cond(__a>0, __hc/__a*(1-__p), 0)
            bysort `gby' (__k): gen double Sfree = exp(sum(ln(max(__p,1e-300))))
            bysort `gby' (__k): gen double __Sp1 = cond(_n==1, 1, Sfree[_n-1])
            bysort `gby' (__k): gen double `cif' = sum(__Sp1*__q)
            drop __hc __hm __a __p __q __Sp1 Ec Em PY __dt
        }
        else if `method'==1 {
            gen double Sfree = Scur/atrisk
            gen double `cif' = Fcur/atrisk
        }
        else {
            gen double __p = cond(Sprev>0, Scur/Sprev, 1)
            gen double __q = cond(Sprev>0, (Fcur-Fprev)/Sprev, 0)
            bysort `gby' (__k): gen double __lnS = sum(ln(max(__p,1e-300)))
            gen double Sfree = exp(__lnS)
            bysort `gby' (__k): gen double __Sp1 = cond(_n==1, 1, Sfree[_n-1])
            bysort `gby' (__k): gen double `cif' = sum(__Sp1*__q)
            drop __p __q __lnS __Sp1
        }
        gen double Dexp = 1 - Sfree - `cif'
        cap drop Scur Sprev Fcur Fprev
        drop __k
        cap drop __g
        label var t_exp  "Time"
        label var atrisk "Number contributing"
        label var `cif'  "Expected cumulative incidence (competing risk of death)"
        label var Sfree  "Expected event-free survival"
        label var Dexp   "Expected cumulative mortality without event"
        order `by' t_exp atrisk `cif' Sfree Dexp
        sort `by' t_exp
    }

    if "`list'"=="" {
        format `cif' Sfree Dexp %6.4f
        cap assert t_exp==int(t_exp)
        if _rc format t_exp %9.4f
        di as txt _n "Expected cumulative incidence allowing for competing mortality" ///
            cond(`method'==1, " (Ederer I)", cond(`method'==2, " (Ederer II)", " (Hakulinen)"))
        if "`by'"!="" list `by' t_exp atrisk `cif' Sfree Dexp, sepby(`by') noobs
        else list t_exp atrisk `cif' Sfree Dexp, noobs
    }
    if "`outnam'" != "" {
        qui save `"`outnam'"', `outrepl'
    }
end
